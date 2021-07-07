#include <Maestro.H>
#include <Postprocess.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real QuadInterp(const Real x, const Real x0, const Real x1, const Real x2,
                const Real y0, const Real y1, const Real y2, bool limit);

// Given a multifab of data (phi), average down to a base state quantity, phibar.
// Assume we are spherical, and the averaging is done at constant radius.

void Postprocess::Average(const Vector<MultiFab>& phi, BaseState<Real>& phibar,
                          int comp) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocessing::Average()", PAverage);

    const auto nr_irreg = base_geom.nr_irreg;

    phibar.setVal(0.0);

    // spherical case with even base state spacing

    // For spherical, we construct a 1D array at each level, phisum, that has space
    // allocated for every possible radius that a cell-center at each level can
    // map into.  The radial locations have been precomputed and stored in radii.
    BaseState<Real> phisum_s(finest_level + 1, nr_irreg + 2);
    auto phisum = phisum_s.array();
    phisum_s.setVal(0.0);
    BaseState<Real> radii_s(finest_level + 1, nr_irreg + 3);
    auto radii = radii_s.array();
    BaseState<int> ncell_s(finest_level + 1, nr_irreg + 2);
    auto ncell = ncell_s.array();
    ncell_s.setVal(0);

    const auto& center_p = center;

    const int fine_lev = finest_level + 1;

    // radii contains every possible distance that a cell-center at the finest
    // level can map into
    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get the index space of the domain
        const auto dx = pgeom[lev].CellSizeArray();

        AMREX_PARALLEL_FOR_1D(nr_irreg + 1, r, {
            radii(lev, r + 1) = std::sqrt(0.75 + 2.0 * Real(r)) * dx[0];
        });
        Gpu::synchronize();

        radii(lev, nr_irreg + 2) = 1.e99;
        radii(lev, 0) = 0.0;
    }

    // loop is over the existing levels (up to finest_level)
    for (int lev = finest_level; lev >= 0; --lev) {
        // Get the grid size of the domain
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

        // get references to the MultiFabs at level lev
        const MultiFab& phi_mf = phi[lev];

        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = phi[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(phi_mf, fba, IntVect(2));

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
        for (MFIter mfi(phi_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const int> mask_arr = mask.array(mfi);
            const Array4<const Real> phi_arr = phi[lev].array(mfi, comp);

            bool use_mask = !(lev == fine_lev - 1);

            AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                // make sure the cell isn't covered by finer cells
                bool cell_valid = true;
                if (use_mask) {
                    if (mask_arr(i, j, k) == 1) cell_valid = false;
                }

                if (cell_valid) {
                    // compute distance to center
                    Real radius = sqrt(x * x + y * y + z * z);

                    // figure out which radii index this point maps into
                    int index = int(round(
                        ((radius / dx[0]) * (radius / dx[0]) - 0.75) / 2.0));

                    // due to roundoff error, need to ensure that we are in the proper radial bin
                    if (index < nr_irreg) {
                        if (fabs(radius - radii(lev, index + 1)) >
                            fabs(radius - radii(lev, index + 2))) {
                            index++;
                        }
                    }

                    amrex::HostDevice::Atomic::Add(&(phisum(lev, index + 1)),
                                                   phi_arr(i, j, k));
                    amrex::HostDevice::Atomic::Add(&(ncell(lev, index + 1)), 1);
                }
            });
        }
    }

    // reduction over boxes to get sum
    ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),
                                      (finest_level + 1) * (nr_irreg + 2));
    ParallelDescriptor::ReduceIntSum(ncell.dataPtr(),
                                     (finest_level + 1) * (nr_irreg + 2));

    // normalize phisum so it actually stores the average at a radius
    for (auto n = 0; n <= finest_level; ++n) {
        for (auto r = 0; r <= nr_irreg; ++r) {
            if (ncell(n, r + 1) != 0) {
                phisum(n, r + 1) /= Real(ncell(n, r + 1));
            }
        }
    }

    BaseState<int> which_lev_s(base_geom.nr_fine);
    auto which_lev = which_lev_s.array();
    BaseState<int> max_rcoord_s(fine_lev);
    auto max_rcoord = max_rcoord_s.array();

    // compute center point for the finest level
    phisum(finest_level, 0) = (11.0 / 8.0) * phisum(finest_level, 1) -
                              (3.0 / 8.0) * phisum(finest_level, 2);
    ncell(finest_level, 0) = 1;

    // choose which level to interpolate from
    const auto dr0 = base_geom.dr(0);
    const auto nrf = base_geom.nr_fine;

    AMREX_PARALLEL_FOR_1D(nrf, r, {
        Real radius = (Real(r) + 0.5) * dr0;
        // Vector<int> rcoord_p(fine_lev, 0);
        int rcoord_p[fine_lev];

        // initialize
        for (int& coord : rcoord_p) {
            coord = 0.0;
        }

        // for each level, find the closest coordinate
        for (auto n = 0; n < fine_lev; ++n) {
            for (auto j = rcoord_p[n]; j <= nr_irreg; ++j) {
                if (fabs(radius - radii(n, j + 1)) <
                    fabs(radius - radii(n, j + 2))) {
                    rcoord_p[n] = j;
                    break;
                }
            }
        }

        // make sure closest coordinate is in bounds
        for (auto n = 0; n < fine_lev - 1; ++n) {
            rcoord_p[n] = amrex::max(rcoord_p[n], 1);
        }
        for (auto n = 0; n < fine_lev; ++n) {
            rcoord_p[n] = amrex::min(rcoord_p[n], nr_irreg - 1);
        }

        // choose the level with the largest min over the ncell interpolation points
        which_lev(r) = 0;

        int min_all =
            amrex::min(ncell(0, rcoord_p[0]), ncell(0, rcoord_p[0] + 1),
                       ncell(0, rcoord_p[0] + 2));

        for (auto n = 1; n < fine_lev; ++n) {
            int min_lev =
                amrex::min(ncell(n, rcoord_p[n]), ncell(n, rcoord_p[n] + 1),
                           ncell(n, rcoord_p[n] + 2));

            if (min_lev > min_all) {
                min_all = min_lev;
                which_lev(r) = n;
            }
        }

        // if the min hit count at all levels is zero, we expand the search
        // to find the closest instance of where the hitcount becomes nonzero
        int j = 1;
        while (min_all == 0) {
            j++;
            for (auto n = 0; n < fine_lev; ++n) {
                int min_lev = amrex::max(
                    ncell(n, amrex::max(1, rcoord_p[n] - j) + 1),
                    ncell(n, amrex::min(rcoord_p[n] + j, nr_irreg - 1) + 1));
                if (min_lev != 0) {
                    which_lev(r) = n;
                    min_all = min_lev;
                    break;
                }
            }
        }
    });
    Gpu::synchronize();

    // squish the list at each level down to exclude points with no contribution
    for (auto n = 0; n <= finest_level; ++n) {
        int j = 0;
        for (auto r = 0; r <= nr_irreg; ++r) {
            while (ncell(n, j + 1) == 0) {
                j++;
                if (j > nr_irreg) {
                    break;
                }
            }
            if (j > nr_irreg) {
                for (auto i = r; i <= nr_irreg; ++i) {
                    phisum(n, i + 1) = 1.e99;
                }
                for (auto i = r; i <= nr_irreg + 1; ++i) {
                    radii(n, i + 1) = 1.e99;
                }
                max_rcoord(n) = r - 1;
                break;
            }
            phisum(n, r + 1) = phisum(n, j + 1);
            radii(n, r + 1) = radii(n, j + 1);
            ncell(n, r + 1) = ncell(n, j + 1);
            j++;
            if (j > nr_irreg) {
                max_rcoord(n) = r;
                break;
            }
        }
    }

    // compute phibar
    const Real drdxfac_loc = drdxfac;
    auto phibar_arr = phibar.array();

    AMREX_PARALLEL_FOR_1D(nrf, r, {
        Real radius = (Real(r) + 0.5) * dr0;
        int stencil_coord = 0;

        // find the closest coordinate
        for (auto j = stencil_coord; j <= max_rcoord(which_lev(r)); ++j) {
            if (fabs(radius - radii(which_lev(r), j + 1)) <
                fabs(radius - radii(which_lev(r), j + 2))) {
                stencil_coord = j;
                break;
            }
        }

        // make sure the interpolation points will be in bounds
        if (which_lev(r) != fine_lev - 1) {
            stencil_coord = amrex::max(stencil_coord, 1);
        }
        stencil_coord = amrex::min(stencil_coord, max_rcoord(which_lev(r)) - 1);

        bool limit = (r <= nrf - 1 - drdxfac_loc * pow(2.0, (fine_lev - 2)));

        phibar_arr(0, r) =
            QuadInterp(radius, radii(which_lev(r), stencil_coord),
                       radii(which_lev(r), stencil_coord + 1),
                       radii(which_lev(r), stencil_coord + 2),
                       phisum(which_lev(r), stencil_coord),
                       phisum(which_lev(r), stencil_coord + 1),
                       phisum(which_lev(r), stencil_coord + 2), limit);
    });
    Gpu::synchronize();
}

AMREX_GPU_DEVICE
Real QuadInterp(const Real x, const Real x0, const Real x1, const Real x2,
                const Real y0, const Real y1, const Real y2, bool limit) {
    Real y = y0 + (y1 - y0) / (x1 - x0) * (x - x0) +
             ((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0) *
                 (x - x0) * (x - x1);

    if (limit) {
        if (y > max(y0, max(y1, y2))) y = max(y0, max(y1, y2));
        if (y < min(y0, min(y1, y2))) y = min(y0, min(y1, y2));
    }

    return y;
}

// Given a multifab of data (phi), average down to quantity on a plane, phibar.
// Assume spherical case, and averaging is done in the psi direction,
// resulting in averages on an x-z plane cutting through the north and south poles.

void Postprocess::Average2d(const Vector<MultiFab>& phi,
                            Vector<MultiFab>& phibar, int barcomp,
                            const Vector<Box>& bardomain, int comp) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocessing::Average2d()", PAverage2d);

    // construct a 2D array at finest level for phi to be mapped onto
    Box domain = bardomain[0];
    const int xlen = domain.hiVect()[0] + 1;
    const int ylen = domain.hiVect()[2] + 1;
    // Print() << "HACK: " << xlen << ", " << ylen << std::endl;

    BaseState<Real> phisum_s(xlen, ylen);
    auto phisum = phisum_s.array();
    phisum_s.setVal(0.0);
    BaseState<int> ncell_s(xlen, ylen);
    auto ncell = ncell_s.array();
    ncell_s.setVal(0);

    const auto& center_p = center;

    const auto dxFine = pgeom[finest_level].CellSizeArray();

    // loop is over the existing levels (up to finest_level)
    for (int lev = finest_level; lev >= 0; --lev) {
        // Get the grid size of the domain
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

        // get references to the MultiFabs at level lev
        const MultiFab& phi_mf = phi[lev];

        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = phi[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(phi_mf, fba, IntVect(2));

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
        for (MFIter mfi(phi_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const int> mask_arr = mask.array(mfi);
            const Array4<const Real> phi_arr = phi[lev].array(mfi, comp);

            bool use_mask = !(lev == finest_level);

            AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                // make sure the cell isn't covered by finer cells
                bool cell_valid = true;
                if (use_mask) {
                    if (mask_arr(i, j, k) == 1) cell_valid = false;
                }

                if (cell_valid) {
                    // compute distance to center
                    Real radius = sqrt(x * x + y * y + z * z);

                    // compute distances in mapped coordinates
                    Real xp = sqrt(x * x + y * y);
                    Real yp = center_p[2] + z;

                    // figure out which (i,j) index this point maps into
                    //  TODO: could give smoother results if coarse grid data is
                    //  mapped into multiple fine grids
                    int index_i = int(round(xp / dxFine[0] - 0.5));
                    int index_j = int(round(yp / dxFine[2] - 0.5));

                    // TODO: currently does not take data in domain corners
                    if (index_i < xlen and index_j < ylen) {
                        amrex::HostDevice::Atomic::Add(
                            &(phisum(index_i, index_j)), phi_arr(i, j, k));
                        amrex::HostDevice::Atomic::Add(
                            &(ncell(index_i, index_j)), 1);
                    }
                }
            });
        }
    }

    // reduction over boxes to get sum
    ParallelDescriptor::ReduceRealSum(phisum.dataPtr(), xlen * ylen);
    ParallelDescriptor::ReduceIntSum(ncell.dataPtr(), xlen * ylen);

    // normalize phisum so it actually stores the averages
    for (auto i = 0; i < xlen; ++i) {
        for (auto k = 0; k < ylen; ++k) {
            if (ncell(i, k) != 0) {
                phisum(i, k) /= Real(ncell(i, k));
            }
        }
    }

    // Put average array into 2D MultiFab
    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
    for (MFIter mfi(phibar[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get the index space of the valid region
        const Box& tilebox = mfi.tilebox();
        const Array4<Real> phibar_arr = phibar[0].array(mfi, barcomp);

        AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
            phibar_arr(i, j, k) = 0.0;
            if (j == 0) phibar_arr(i, j, k) = phisum(i, k);
        });
    }
}
