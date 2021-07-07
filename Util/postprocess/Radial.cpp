#include <Maestro.H>
#include <Postprocess.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// ---------------------------------
// Write 1D radial diagnostics
// ---------------------------------
void Postprocess::WriteRadialFile(const BaseState<Real>& rho0_in,
                                  const BaseState<Real>& p0_in,
                                  const Vector<MultiFab>& u_in,
                                  const Vector<MultiFab>& w0_in) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::WriteRadialFile()", WriteRadialFile);

    // MakeVelrc
    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> circ_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rad_vel[lev].define(u_in[lev].boxArray(), u_in[lev].DistributionMap(),
                            1, 0);
        circ_vel[lev].define(u_in[lev].boxArray(), u_in[lev].DistributionMap(),
                             1, 0);
    }
    MakeVelrc(u_in, w0_in, rad_vel, circ_vel);

    // MakeConvectionVel
    BaseState<Real> convect_vel(base_geom.max_radial_level + 1,
                                base_geom.nr_fine);
    MakeConvectionVel(rad_vel, convect_vel);

    // MakeRotationRate
    Vector<MultiFab> omega(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        omega[lev].define(u_in[lev].boxArray(), u_in[lev].DistributionMap(), 1,
                          0);
    }
    MakeRotationRate(u_in, w0_in, omega);

    // MakeRadialRotationRatio
    BaseState<Real> ratio_omega(base_geom.max_radial_level + 1,
                                base_geom.nr_fine);
    MakeRadialRotationRatio(p0_in, omega, ratio_omega);

    // MakeLatitudinalShear
    BaseState<Real> latshear(base_geom.max_radial_level + 1, base_geom.nr_fine);
    MakeLatShear(omega, latshear, base_geom.r_cc_loc);
    //MakeLatShearAvg(omega, latshear);

    // MakeRadialNFreq
    BaseState<Real> Nfreq(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> s0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> grav0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    MakeRadialNFreq(p0_in, rho0_in, Nfreq, s0, grav0);

    std::string radialfilename = "radial_" + iFile;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    // write out the radial cell-centered diagnostics
    if (ParallelDescriptor::IOProcessor()) {
        for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
            std::ofstream RadialFile;
            RadialFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
                                          io_buffer.size());
            RadialFile.open(radialfilename.c_str(), std::ofstream::out |
                                                        std::ofstream::trunc |
                                                        std::ofstream::binary);
            if (!RadialFile.good()) {
                amrex::FileOpenFailed(radialfilename);
            }

            RadialFile.precision(12);

            RadialFile << "r_cc               "
                       << "rho0               "
                       << "p0                 "
                       << "convect_vel        "
                       << "omega_ratio        "
                       << "lat_shear          "
                       << "|N|\n";

            for (int i = 0; i < base_geom.nr(lev); ++i) {
                RadialFile << std::left << std::setw(18)
                           << base_geom.r_cc_loc(lev, i) << " " << std::setw(18)
                           << rho0_in.array()(lev, i) << " " << std::setw(18)
                           << p0_in.array()(lev, i) << " " << std::setw(18)
                           << convect_vel.array()(lev, i) << " "
                           << std::setw(18) << ratio_omega.array()(lev, i)
                           << " " << std::setw(18) << latshear.array()(lev, i)
                           << " " << std::setw(18) << Nfreq.array()(lev, i)
                           << "\n";
            }
        }
    }
}

// write additional diagnostics from model file
void Postprocess::WriteModelDiagFile(const std::string& plotfilename) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::WriteModelDiagFile()", WriteModelDiagFile);

    // Read from model file
    BaseState<Real> r_cc, rho0, p0;
    ReadModelFile(plotfilename, r_cc, rho0, p0);

    int npts = rho0.length();

    // make dimensionless entropy
    BaseState<Real> s0(1, npts);
    const auto rho0_arr = rho0.const_array();
    const auto p0_arr = p0.const_array();
    auto entropy = s0.array();
    Real gamma = 5.0 / 3.0;
    for (auto r = 0; r < npts; ++r) {
        entropy(0, r) = 1.0 / (gamma - 1.0) *
                        (log(p0_arr(0, r)) - gamma * log(rho0_arr(0, r)));
    }

    // make gravity
    BaseState<Real> grav0(1, npts);
    MakeRadialGrav(rho0, grav0, r_cc.array());

    // make B-V frequency
    BaseState<Real> Nfreq(1, npts);
    const auto r_cc_loc = r_cc.const_array();
    const auto grav = grav0.const_array();
    auto freq = Nfreq.array();
    for (auto r = 0; r < npts; ++r) {
        if (r == 0) {
            freq(0, r) = -(gamma - 1.0) / gamma * grav(0, r) *
                         (entropy(0, r + 1) - entropy(0, r)) /
                         (r_cc_loc(0, r + 1) - r_cc_loc(0, r));
        } else if (p0_arr(0, r) > 1.e12) {
            freq(0, r) = -(gamma - 1.0) / gamma * grav(0, r) *
                         (entropy(0, r + 1) - entropy(0, r - 1)) /
                         (r_cc_loc(0, r + 1) - r_cc_loc(0, r - 1));
        } else {
            freq(0, r) = 0.0;
        }
        freq(0, r) = std::sqrt(fabs(freq(0, r)));
    }

    std::string diagfilename = "diag_" + plotfilename;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    // write out the cell-centered diagnostics
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream DiagFile;
        DiagFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        DiagFile.open(diagfilename.c_str(), std::ofstream::out |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if (!DiagFile.good()) {
            amrex::FileOpenFailed(diagfilename);
        }

        DiagFile.precision(17);

        DiagFile << "r_cc  rho0  p0  s0  grav0  |N| \n";

        for (int i = 0; i < npts; ++i) {
            DiagFile << r_cc.array()(0, i) << " " << rho0.array()(0, i) << " "
                     << p0.array()(0, i) << " " << s0.array()(0, i) << " "
                     << grav0.array()(0, i) << " " << Nfreq.array()(0, i)
                     << "\n";
        }
    }
}

// read BaseCC file from plotfile directory
void Postprocess::ReadBaseCCFile(BaseState<Real>& r_s, BaseState<Real>& rho0_s,
                                 BaseState<Real>& p0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::ReadBaseCCFile()", ReadBaseCCFile);

    auto r = r_s.array();
    auto rho0 = rho0_s.array();
    auto p0 = p0_s.array();
    BaseState<Real> rhoh0_s(base_geom.max_radial_level + 1, base_geom.nr_fine);

    // read BaseCC
    std::string line, word;
    std::string File(iFile + "/BaseCC_0");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    // read in cell-centered base state
    std::getline(is, line);  // header line
    for (int i = 0; i < (base_geom.max_radial_level + 1) * base_geom.nr_fine;
         ++i) {
        std::getline(is, line);
        std::istringstream lis(line);
        lis >> word;
        r(i) = std::stod(word);
        lis >> word;
        rho0(i) = std::stod(word);
        lis >> word;
        rhoh0_s.array()(i) = std::stod(word);
        lis >> word;
        p0(i) = std::stod(word);
        lis >> word;
        // gamma1bar isn't needed for now
    }
}

// read input model file from disk
void Postprocess::ReadModelFile(const std::string& plotfilename,
                                BaseState<Real>& r_s, BaseState<Real>& rho0_s,
                                BaseState<Real>& p0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::ReadModelFile()", ReadModelFile);

    // read model file
    std::string line, word;
    std::string File(plotfilename);
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    // read in headers
    // get npts, num vars
    int npts, nvars, temp;
    for (int i = 0; i < 2; ++i) {
        std::getline(is, line);
        std::istringstream lis(line);
        while (lis) {
            lis >> word;
            if (std::stringstream(word) >> temp) {
                if (i == 0)
                    npts = temp;
                else if (i == 1)
                    nvars = temp;
            }
        }
    }
    // skip rest of headers
    for (int i = 0; i < nvars; ++i) std::getline(is, line);

    // read in data
    BaseState<Real> temp0_s(1, npts);
    r_s.resize(1, npts);
    rho0_s.resize(1, npts);
    p0_s.resize(1, npts);
    auto r = r_s.array();
    auto rho0 = rho0_s.array();
    auto p0 = p0_s.array();

    for (int i = 0; i < npts; ++i) {
        std::getline(is, line);
        std::istringstream lis(line);
        lis >> word;
        r(0, i) = std::stod(word);
        lis >> word;
        rho0(0, i) = std::stod(word);
        lis >> word;
        temp0_s.array()(0, i) = std::stod(word);
        lis >> word;
        p0(0, i) = std::stod(word);
        // species not needed
    }
}

void Postprocess::MakeRadialNFreq(const BaseState<Real>& p0_s,
                                  const BaseState<Real>& rho0_s,
                                  BaseState<Real>& freq0,
                                  BaseState<Real>& entropy_s,
                                  BaseState<Real>& grav_s) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeRadialNFreq()", MakeRadialNFreq);

    //BaseState<Real> entropy_s(1, base_geom.nr_fine);
    //BaseState<Real> grav_s(1, base_geom.nr_fine);

    // compute gravity
    // analogous to Maestro::MakeGravCell(grav_s,rho0_s);
    MakeRadialGrav(rho0_s, grav_s, base_geom.r_cc_loc);

    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto rho0 = rho0_s.const_array();
    const auto p0 = p0_s.const_array();
    const auto grav = grav_s.const_array();
    auto entropy = entropy_s.array();
    Real gamma = GetGamma(iFile);

    // dimensionless entropy = 1/(gam - 1)*( log(p) - gamma*log(rho) )
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        entropy(0, r) =
            1.0 / (gamma - 1.0) * (log(p0(0, r)) - gamma * log(rho0(0, r)));
    }

    // B-V freq = -(gam -1)/gam * g * grad(entropy)
    auto freq = freq0.array();
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        if (r == 0) {
            freq(0, r) = -(gamma - 1.0) / gamma * grav(0, r) *
                         (entropy(0, r + 1) - entropy(0, r)) /
                         (r_cc_loc(0, r + 1) - r_cc_loc(0, r));
        } else if (p0(0, r) > 1.e12) {
            freq(0, r) = -(gamma - 1.0) / gamma * grav(0, r) *
                         (entropy(0, r + 1) - entropy(0, r - 1)) /
                         (r_cc_loc(0, r + 1) - r_cc_loc(0, r - 1));
        } else {
            freq(0, r) = 0.0;
        }
        freq(0, r) = std::sqrt(fabs(freq(0, r)));
    }
}

void Postprocess::MakeRadialGrav(const BaseState<Real>& rho0_in,
                                 BaseState<Real>& grav,
                                 const BaseStateArray<Real>& r_cc_loc) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeRadialGrav()", MakeRadialGrav);

    auto grav_arr = grav.array();
    const auto rho0 = rho0_in.const_array();
    const auto nr_fine = r_cc_loc.length();

    // assume spherical
    BaseState<Real> m_state(1, nr_fine);
    auto m = m_state.array();

    m(0, 0) = 4.0 / 3.0 * M_PI * rho0(0, 0) * r_cc_loc(0, 0) * r_cc_loc(0, 0) *
              r_cc_loc(0, 0);
    grav_arr(0, 0) = -C::Gconst * m(0, 0) / (r_cc_loc(0, 0) * r_cc_loc(0, 0));

    for (auto r = 1; r < nr_fine; ++r) {
        // the mass is defined at the cell-centers, so to compute
        // the mass at the current center, we need to add the
        // contribution of the upper half of the zone below us and
        // the lower half of the current zone.

        // don't add any contributions from outside the star --
        // i.e.  rho < base_cutoff_density

        // assume constant radial spacing, dr
        Real r_edge_loc = (r_cc_loc(0, r - 1) + r_cc_loc(0, r)) / 2.0;

        Real term1 = 0.0;
        if (rho0(0, r - 1) > base_cutoff_density) {
            term1 = 4.0 / 3.0 * M_PI * rho0(0, r - 1) *
                    (r_edge_loc - r_cc_loc(0, r - 1)) *
                    (r_edge_loc * r_edge_loc + r_edge_loc * r_cc_loc(0, r - 1) +
                     r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
        }

        Real term2 = 0.0;
        if (rho0(0, r) > base_cutoff_density) {
            term2 = 4.0 / 3.0 * M_PI * rho0(0, r) *
                    (r_cc_loc(0, r) - r_edge_loc) *
                    (r_cc_loc(0, r) * r_cc_loc(0, r) +
                     r_cc_loc(0, r) * r_edge_loc + r_edge_loc * r_edge_loc);
        }

        m(0, r) = m(0, r - 1) + term1 + term2;

        grav_arr(0, r) =
            -C::Gconst * m(0, r) / (r_cc_loc(0, r) * r_cc_loc(0, r));
    }
}

void Postprocess::MakeConvectionVel(const Vector<MultiFab>& velr,
                                    BaseState<Real>& vel_conv) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeConvectionVel()", MakeConvectionVel);

    Vector<MultiFab> vel2(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel2[lev].define(velr[lev].boxArray(), velr[lev].DistributionMap(), 1,
                         0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(velr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> velr_arr = velr[lev].array(mfi);
            const Array4<Real> vel2_arr = vel2[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                vel2_arr(i, j, k) = 0.0;

                // square radial vel
                vel2_arr(i, j, k) = velr_arr(i, j, k) * velr_arr(i, j, k);
            });
        }
    }

    // average down and fill ghost cells
    // AverageDown(vel2, 0, 1);
    // FillPatch(t0, vel2, vel2, vel2, 0, 0, 1, 0, bcs_f);

    // radial average of square of radial vel
    Average(vel2, vel_conv, 0);

    // root-mean-squared radial velocity
    auto vel_conv_arr = vel_conv.array();

    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        vel_conv_arr(0, r) = std::sqrt(vel_conv_arr(0, r));
    }
}

void Postprocess::MakeRadialRotationRatio(const BaseState<Real>& s0_in,
                                          const Vector<MultiFab>& omega,
                                          BaseState<Real>& ratio_omega) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::RadialRotationRatio()",
                   MakeRadialRotationRatio);

    // radial average of rotation rate
    Average(omega, ratio_omega, 0);

    // find radius of surface
    int r_coord_surface = base_geom.nr_fine - 1;
    const auto s0 = s0_in.const_array();
    for (auto r = 0; r < base_geom.nr_fine - 1; ++r) {
        if (s0(0, r) == s0(0, r + 1)) {
            r_coord_surface = r;
            break;
        }
    }

    // divide by rotation rate at surface
    auto ratio_omega_arr = ratio_omega.array();

    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        if (r < r_coord_surface)
            ratio_omega_arr(0, r) /= ratio_omega_arr(0, r_coord_surface);
        else
            ratio_omega_arr(0, r) = 1.0;
    }
}

void Postprocess::MakeVelrc(const Vector<MultiFab>& vel,
                            const Vector<MultiFab>& w0rcart,
                            Vector<MultiFab>& rad_vel,
                            Vector<MultiFab>& circ_vel) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeVelrc()", MakeVelrc);

    const auto& center_p = center;

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
            const Array4<Real> circvel_arr = circ_vel[lev].array(mfi);
            const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                circvel_arr(i, j, k) = 0.0;
                radvel_arr(i, j, k) = 0.0;

                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
                Real inv_radius = 1.0 / sqrt(x * x + y * y + z * z);

                Vector<Real> normal(3);
                normal[0] = x * inv_radius;
                normal[1] = y * inv_radius;
                normal[2] = z * inv_radius;

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    radvel_arr(i, j, k) += vel_arr(i, j, k, n) * normal[n];
                }

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    Real circ_comp =
                        vel_arr(i, j, k, n) - radvel_arr(i, j, k) * normal[n];
                    circvel_arr(i, j, k) += circ_comp * circ_comp;
                }

                circvel_arr(i, j, k) = sqrt(circvel_arr(i, j, k));

                // add base state vel to get full radial velocity
                radvel_arr(i, j, k) += w0rcart_arr(i, j, k);
            });
        }
    }

    // average down and fill ghost cells
    // AverageDown(rad_vel, 0, 1);
    // FillPatch(t0, rad_vel, rad_vel, rad_vel, 0, 0, 1, 0, bcs_f);
    // AverageDown(circ_vel, 0, 1);
    // FillPatch(t0, circ_vel, circ_vel, circ_vel, 0, 0, 1, 0, bcs_f);
}

void Postprocess::MakeRotationRate(const Vector<MultiFab>& vel,
                                   const Vector<MultiFab>& w0cart,
                                   Vector<MultiFab>& omega) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeRotationRate()", MakeRotationRate);

    auto omega0 = GetRotationFreq(iFile);
    const auto& center_p = center;

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<const Real> w0cart_arr = w0cart[lev].array(mfi);
            const Array4<Real> omega_arr = omega[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                omega_arr(i, j, k) = 0;

                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
                Real inv_radius = 1.0 / sqrt(x * x + y * y + z * z);
                Real inv_rsin = 1.0 / sqrt(x * x + y * y);

                Real phi_x = -y * inv_radius * inv_rsin;
                Real phi_y = x * inv_radius * inv_rsin;

                omega_arr(i, j, k) =
                    phi_x * (vel_arr(i, j, k, 0) + w0cart_arr(i, j, k)) +
                    phi_y * (vel_arr(i, j, k, 1) + w0cart_arr(i, j, k));

                // add constant rotation to get full rotation
                omega_arr(i, j, k) += omega0;
            });
        }
    }

    // average down and fill ghost cells
    // AverageDown(omega, 0, 1);
    // FillPatch(t0, omega, omega, omega, 0, 0, 1, 0, bcs_f);
}

void Postprocess::MakeLatShear(const Vector<MultiFab>& omega_in,
                               BaseState<Real>& shear,
                               const BaseStateArray<Real>& r_cc_loc) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeLatShear()", MakeLatShear);

    // Vector<MultiFab> shear_cart(finest_level+1);
    // for (int lev = 0; lev <= finest_level; ++lev) {
    // 	shear_cart[lev].define(grid[lev], dmap[lev], 1, 0);
    // }

    auto shear_arr = shear.array();
    const auto nr_fine = r_cc_loc.length();
    const auto& center_p = center;

    const Real maxfac = 0.5;
    Real dx_fine;
    const auto probLo = pgeom[0].ProbLoArray();
    const auto probHi = pgeom[0].ProbHiArray();
    Real halfdom =
        0.5 * amrex::min(probHi[0] - probLo[0], probHi[1] - probLo[1]);
#if AMREX_SPACEDIM == 3
    halfdom = amrex::min(halfdom, 0.5 * (probHi[2] - probLo[2]));
#endif

    for (int r = 0; r < nr_fine; ++r) {
        Real rr = r_cc_loc(0, r);
        Real totkernel = 0.0;
        shear_arr(0, r) = 0.0;

        for (int lev = finest_level; lev >= 0; --lev) {
            // Get grid size of domain
            const auto dx = pgeom[lev].CellSizeArray();
            const auto prob_lo = pgeom[lev].ProbLoArray();

            // create mask assuming refinement ratio = 2
            int finelev = lev + 1;
            if (lev == finest_level) {
                finelev = finest_level;
                dx_fine = dx[0];
            }

            const BoxArray& fba = omega_in[finelev].boxArray();
            const iMultiFab& mask =
                makeFineMask(omega_in[lev], fba, IntVect(2));

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction)
#endif
            for (MFIter mfi(omega_in[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                const Array4<const int> mask_arr = mask.array(mfi);
                const Array4<const Real> omega_arr = omega_in[lev].array(mfi);

                bool use_mask = !(lev == finest_level);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
                    Real radius = std::sqrt(x * x + y * y + z * z);

                    // make sure the cell isn't covered by finer cells
                    bool cell_valid = true;
                    if (use_mask) {
                        if (mask_arr(i, j, k) == 1) {
                            cell_valid = false;
                        }
                    }

                    if (cell_valid) {
                        // Y_{2,0} spherical harmonic
                        Real Y20 = 0.25 * std::sqrt(5.0 / M_PI) *
                                   (2 * z * z - x * x - y * y) /
                                   (radius * radius);

                        // normalized Gaussian
                        Real width = amrex::min(maxfac, rr / halfdom) *
                                     dx_fine * dx_fine;
                        Real kernel =
                            std::exp(-(radius - rr) * (radius - rr) / width) /
                            std::sqrt(M_PI * width);

                        amrex::HostDevice::Atomic::Add(
                            &(shear_arr(0, r)),
                            omega_arr(i, j, k) * Y20 * kernel);
                        amrex::HostDevice::Atomic::Add(&totkernel, kernel);
                    }
                });
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(&(shear_arr(0, r)), 1);
        ParallelDescriptor::ReduceRealSum(&totkernel, 1);

        // normalize shear so it actually stores the average at radius r
        if (totkernel != 0.0) {
            shear_arr(0, r) /= totkernel;
        }
    }
}

void Postprocess::MakeLatShearAvg(const Vector<MultiFab>& omega_in,
                                  BaseState<Real>& shear) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::MakeLatShearAvg()", MakeLatShearAvg);

    Vector<MultiFab> omegaY(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        omegaY[lev].define(omega_in[lev].boxArray(),
                           omega_in[lev].DistributionMap(), 1, 0);
    }

    const auto& center_p = center;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get grid size of domain
        const auto dx = pgeom[lev].CellSizeArray();
        const auto prob_lo = pgeom[lev].ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(omega_in[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> omega_arr = omega_in[lev].array(mfi);
            const Array4<Real> omegaY_arr = omegaY[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
                Real radius = std::sqrt(x * x + y * y + z * z);

                // Y_{2,0} spherical harmonic
                Real Y20 = 0.25 * std::sqrt(5.0 / M_PI) *
                           (2 * z * z - x * x - y * y) / (radius * radius);

                // square radial vel
                omegaY_arr(i, j, k) = omega_arr(i, j, k) * Y20;
            });
        }
    }

    // average down and fill ghost cells
    // AverageDown(omegaY, 0, 1);
    // FillPatch(t0, omegaY, omegaY, omegaY, 0, 0, 1, 0, bcs_f);

    // radial average of omega*Y20
    Average(omegaY, shear, 0);
}
