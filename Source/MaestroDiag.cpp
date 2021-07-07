
#include <AMReX_buildInfo.H>
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// Choose precision of all output files
const int outfilePrecision = 10;
const int setwVal =
    outfilePrecision + 2 + 4 + 4;  // 0. + precision + 4 for exp + 4 for gap

// write diagnostics files to disk
// We hold many timesteps-worth of diagnostic information in a buffer
// and output to the files only when flush_diag() is called.  This
// gives better performance on large machines with slow filesystems.
void Maestro::DiagFile(const int step, const Real t_in,
                       const BaseState<Real>& rho0_in,
                       const BaseState<Real>& p0_in,
                       const Vector<MultiFab>& u_in,
                       const Vector<MultiFab>& s_in, int& index) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DiagFile()", DiagFile);

    // -- w0mac will contain an edge-centered w0 on a Cartesian grid,
    // -- for use in computing divergences.
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);
    // -- w0r_cart is w0 but onto a Cartesian grid in cell-centered as
    // -- a scalar.  Since w0 is the radial expansion velocity, w0r_cart
    // -- is the radial w0 in a zone
    Vector<MultiFab> w0r_cart(finest_level + 1);

    // rho_Hnuc and rho_Hext are used to determine energy generation
    Vector<MultiFab> stemp(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);
    Vector<MultiFab> rho_omegadot(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
#if NAUX_NET > 0
    Vector<MultiFab> rho_auxdot(finest_level + 1);
#endif
    Vector<MultiFab> sdc_source(finest_level + 1);

#if (AMREX_SPACEDIM == 3)
    if (spherical) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            AMREX_D_TERM(
                w0mac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                     dmap[lev], 1, 1);
                , w0mac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                       dmap[lev], 1, 1);
                , w0mac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                       dmap[lev], 1, 1););

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }

            w0r_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            w0r_cart[lev].setVal(0.);
        }

        // put w0 on Cartesian edges as a vector
        MakeW0mac(w0mac);

        // put w0 in Cartesian cell-centers as a scalar (the radial
        // expansion velocity)
        Put1dArrayOnCart(w0, w0r_cart, true, false, bcs_u, 0, 1);
    }
#endif

    // compute rho_Hext and rho_Hnuc
    for (int lev = 0; lev <= finest_level; ++lev) {
        stemp[lev].define(grids[lev], dmap[lev], Nscal, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
#if NAUX_NET > 0
        rho_auxdot[lev].define(grids[lev], dmap[lev], NumAux, 0);
#endif
        sdc_source[lev].define(grids[lev], dmap[lev], Nscal, 0);

        sdc_source[lev].setVal(0.0);
    }

#ifndef SDC
    if (dt < small_dt) {
        React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc,
#if NAUX_NET > 0
              rho_auxdot,
#endif
              p0_in, small_dt, t_in);
    } else {
        React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc,
#if NAUX_NET > 0
              rho_auxdot,
#endif
              p0_in, dt * 0.5, t_in);
    }
#else
    if (dt < small_dt) {
        ReactSDC(s_in, stemp, rho_Hext, p0_in, small_dt, t_in, sdc_source);
    } else {
        ReactSDC(s_in, stemp, rho_Hext, p0_in, dt * 0.5, t_in, sdc_source);
    }

    MakeReactionRates(rho_omegadot, rho_Hnuc, s_in);
#endif

    // initialize diagnosis variables
    // diag_temp.out
    Real T_max = 0.0, T_center = 0.0;
    int ncenter = 0;
    Vector<Real> coord_Tmax(AMREX_SPACEDIM, 0.0);
    Vector<Real> vel_Tmax(AMREX_SPACEDIM, 0.0);
    Real Rloc_Tmax = 0.0, vr_Tmax = 0.0;

    // diag_vel.out
    Real U_max = 0.0, Mach_max = 0.0;
    Real kin_ener = 0.0, int_ener = 0.0;
    Vector<Real> vel_center(AMREX_SPACEDIM, 0.0);

    // diag_enuc.out
    Real enuc_max = 0.0;
    Vector<Real> coord_enucmax(AMREX_SPACEDIM, 0.0);
    Vector<Real> vel_enucmax(AMREX_SPACEDIM, 0.0);
    Real Rloc_enucmax = 0.0, vr_enucmax = 0.0;
    Real nuc_ener = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // diagnosis variables at each level
        // diag_temp.out
        Real T_max_local = 0.0;
        Real T_center_level = 0.0;
        int ncenter_level = 0;
        Vector<Real> coord_Tmax_local(AMREX_SPACEDIM, 0.0);
        Vector<Real> vel_Tmax_local(AMREX_SPACEDIM, 0.0);

        // diag_vel.out
        Real U_max_level = 0.0;
        Real Mach_max_level = 0.0;
        Real kin_ener_level = 0.0;
        Real int_ener_level = 0.0;
        Vector<Real> vel_center_level(AMREX_SPACEDIM, 0.0);

        // diag_enuc.out
        Real enuc_max_local = 0.0;
        Vector<Real> coord_enucmax_local(AMREX_SPACEDIM, 0.0);
        Vector<Real> vel_enucmax_local(AMREX_SPACEDIM, 0.0);
        Real nuc_ener_level = 0.0;

        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // create mask assuming refinement ratio = 2
        int finelev = lev + 1;
        if (lev == finest_level) {
            finelev = finest_level;
        }

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in[lev], fba, IntVect(2));

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_ener_level) reduction(+:int_ener_level) reduction(+:nuc_ener_level) reduction(max:U_max_level) reduction(max:Mach_max_level)
#endif
        for (MFIter mfi(s_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const auto lo = amrex::lbound(tileBox);
            const auto hi = amrex::ubound(tileBox);

            const auto use_mask = !(lev == finest_level);

            const Array4<const Real> scal = s_in[lev].array(mfi);
            const Array4<const Real> rho_Hnuc_arr = rho_Hnuc[lev].array(mfi);
            const Array4<const Real> u = u_in[lev].array(mfi);
            const Array4<const int> mask_arr = mask.array(mfi);
            const auto w0_arr = w0.const_array();

            // weight is the factor by which the volume of a cell at the current level
            // relates to the volume of a cell at the coarsest level of refinement.
            const auto weight =
                AMREX_SPACEDIM == 2 ? 1.0 / pow(4.0, lev) : 1.0 / pow(8.0, lev);

#if (AMREX_SPACEDIM == 3)
            // for non-spherical we have to set these to be a valid array as
            // w0macx etc are not defined
            const Array4<const Real> w0macx =
                spherical ? w0mac[lev][0].array(mfi) : rho_Hnuc[lev].array(mfi);
            const Array4<const Real> w0macy =
                spherical ? w0mac[lev][1].array(mfi) : rho_Hnuc[lev].array(mfi);
            const Array4<const Real> w0macz =
                spherical ? w0mac[lev][2].array(mfi) : rho_Hnuc[lev].array(mfi);
            const Array4<const Real> normal_arr =
                spherical ? normal[lev].array(mfi) : rho_Hnuc[lev].array(mfi);
            const Array4<const Real> w0r =
                spherical ? w0r_cart[lev].array(mfi) : rho_Hnuc[lev].array(mfi);
#endif

            // The locations of the maxima here make trying to do this on the
            // GPU probably more trouble than it's worth.
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
                        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                        const Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];
#endif

                        // make sure the cell isn't covered by finer cells
                        bool cell_valid = true;
                        if (use_mask) {
                            if (mask_arr(i, j, k) == 1) cell_valid = false;
                        }

                        // For spherical, we only consider cells inside of where the
                        // sponging begins
                        if (cell_valid &&
                            (!spherical ||
                             scal(i, j, k, Rho) >=
                                 sponge_start_factor * sponge_center_density)) {
                            Real vel = 0.0;
                            if (spherical) {
#if (AMREX_SPACEDIM == 3)
                                // is it one of the 8 zones surrounding the center?
                                if (amrex::Math::abs(x - center[0]) < dx[0] &&
                                    amrex::Math::abs(y - center[1]) < dx[1] &&
                                    amrex::Math::abs(z - center[2]) < dx[2]) {
                                    ncenter_level++;

                                    T_center_level += scal(i, j, k, Temp);

                                    vel_center_level[0] +=
                                        u(i, j, k, 0) +
                                        0.5 * (w0macx(i, j, k) +
                                               w0macx(i + 1, j, k));
                                    vel_center_level[1] +=
                                        u(i, j, k, 1) +
                                        0.5 * (w0macy(i, j, k) +
                                               w0macy(i, j + 1, k));
                                    vel_center_level[2] +=
                                        u(i, j, k, 2) +
                                        0.5 * (w0macz(i, j, k) +
                                               w0macz(i, j, k + 1));
                                }

                                // velr is the projection of the velocity (including w0) onto
                                // the radial unit vector
                                // Real velr = u(i,j,k,0)*normal_arr(i,j,k,0) + \
                        //     u(i,j,k,1)*normal_arr(i,j,k,1) + \
                        //     u(i,j,k,2)*normal_arr(i,j,k,2) + w0r(i,j,k);

                                // vel is the magnitude of the velocity, including w0
                                vel = std::sqrt(
                                    (u(i, j, k, 0) +
                                     0.5 * (w0macx(i, j, k) +
                                            w0macx(i + 1, j, k))) *
                                        (u(i, j, k, 0) +
                                         0.5 * (w0macx(i, j, k) +
                                                w0macx(i + 1, j, k))) +
                                    (u(i, j, k, 1) +
                                     0.5 * (w0macy(i, j, k) +
                                            w0macy(i, j + 1, k))) *
                                        (u(i, j, k, 1) +
                                         0.5 * (w0macy(i, j, k) +
                                                w0macy(i, j + 1, k))) +
                                    (u(i, j, k, 2) +
                                     0.5 * (w0macz(i, j, k) +
                                            w0macz(i, j, k + 1))) *
                                        (u(i, j, k, 2) +
                                         0.5 * (w0macz(i, j, k) +
                                                w0macz(i, j, k + 1))));

                                // max T, location, and velocity at that location (including w0)
                                if (scal(i, j, k, Temp) > T_max_local) {
                                    T_max_local = scal(i, j, k, Temp);
                                    coord_Tmax_local[0] = x;
                                    coord_Tmax_local[1] = y;
                                    coord_Tmax_local[2] = z;
                                    vel_Tmax_local[0] =
                                        u(i, j, k, 0) +
                                        0.5 * (w0macx(i, j, k) +
                                               w0macx(i + 1, j, k));
                                    vel_Tmax_local[1] =
                                        u(i, j, k, 1) +
                                        0.5 * (w0macy(i, j, k) +
                                               w0macy(i, j + 1, k));
                                    vel_Tmax_local[2] =
                                        u(i, j, k, 2) +
                                        0.5 * (w0macz(i, j, k) +
                                               w0macz(i, j, k + 1));
                                }

                                // max enuc
                                if (rho_Hnuc_arr(i, j, k) / scal(i, j, k, Rho) >
                                    enuc_max_local) {
                                    enuc_max_local = rho_Hnuc_arr(i, j, k) /
                                                     scal(i, j, k, Rho);
                                    coord_enucmax_local[0] = x;
                                    coord_enucmax_local[1] = y;
                                    coord_enucmax_local[2] = z;
                                    vel_enucmax_local[0] =
                                        u(i, j, k, 0) +
                                        0.5 * (w0macx(i, j, k) +
                                               w0macx(i + 1, j, k));
                                    vel_enucmax_local[1] =
                                        u(i, j, k, 1) +
                                        0.5 * (w0macy(i, j, k) +
                                               w0macy(i, j + 1, k));
                                    vel_enucmax_local[2] =
                                        u(i, j, k, 2) +
                                        0.5 * (w0macz(i, j, k) +
                                               w0macz(i, j, k + 1));
                                }
#endif
                            } else {
                                // vel is the magnitude of the velocity, including w0
#if (AMREX_SPACEDIM == 2)
                                Real vert_vel =
                                    u(i, j, k, 1) +
                                    0.5 * (w0_arr(lev, j) + w0_arr(lev, j + 1));
                                vel = std::sqrt(u(i, j, k, 0) * u(i, j, k, 0) +
                                                vert_vel * vert_vel);
#else
                                Real vert_vel =
                                    u(i, j, k, 2) +
                                    0.5 * (w0_arr(lev, k) + w0_arr(lev, k + 1));
                                vel = std::sqrt(u(i, j, k, 0) * u(i, j, k, 0) +
                                                u(i, j, k, 1) * u(i, j, k, 1) +
                                                vert_vel * vert_vel);
#endif

                                // max T, location, and velocity at that location (including w0)
                                if (scal(i, j, k, Temp) > T_max_local) {
                                    T_max_local = scal(i, j, k, Temp);
                                    coord_Tmax_local[0] = x;
                                    coord_Tmax_local[1] = y;
#if (AMREX_SPACEDIM == 3)
                                    coord_Tmax_local[2] = z;
#endif
                                    vel_Tmax_local[0] = u(i, j, k, 0);
                                    vel_Tmax_local[1] = u(i, j, k, 1);
#if (AMREX_SPACEDIM == 2)
                                    vel_Tmax_local[1] +=
                                        0.5 *
                                        (w0_arr(lev, j) + w0_arr(lev, j + 1));
#else
                                    vel_Tmax_local[2] =
                                        u(i, j, k, 2) +
                                        0.5 * (w0_arr(lev, k) +
                                               w0_arr(lev, k + 1));
#endif
                                }

                                // max enuc
                                if (rho_Hnuc_arr(i, j, k) / scal(i, j, k, Rho) >
                                    enuc_max_local) {
                                    enuc_max_local = rho_Hnuc_arr(i, j, k) /
                                                     scal(i, j, k, Rho);
                                    coord_enucmax_local[0] = x;
                                    coord_enucmax_local[1] = y;
#if (AMREX_SPACEDIM == 3)
                                    coord_enucmax_local[2] = z;
#endif
                                    vel_enucmax_local[0] = u(i, j, k, 0);
                                    vel_enucmax_local[1] = u(i, j, k, 1);
#if (AMREX_SPACEDIM == 2)
                                    vel_enucmax_local[1] +=
                                        0.5 *
                                        (w0_arr(lev, j) + w0_arr(lev, j + 1));
#else
                                    vel_enucmax_local[2] =
                                        u(i, j, k, 2) +
                                        0.5 * (w0_arr(lev, k) +
                                               w0_arr(lev, k + 1));
#endif
                                }
                            }

                            eos_t eos_state;

                            // call the EOS to get the sound speed and internal energy
                            eos_state.T = scal(i, j, k, Temp);
                            eos_state.rho = scal(i, j, k, Rho);
                            for (auto comp = 0; comp < NumSpec; ++comp) {
                                eos_state.xn[comp] =
                                    scal(i, j, k, FirstSpec + comp) /
                                    eos_state.rho;
                            }
#if NAUX_NET > 0
                            for (auto comp = 0; comp < NumAux; ++comp) {
                                eos_state.aux[comp] =
                                    scal(i, j, k, FirstAux + comp) /
                                    eos_state.rho;
                            }
#endif

                            eos(eos_input_rt, eos_state);

                            // kinetic, internal, and nuclear energies
                            kin_ener_level +=
                                weight * scal(i, j, k, Rho) * vel * vel;
                            int_ener_level +=
                                weight * scal(i, j, k, Rho) * eos_state.e;
                            nuc_ener_level += weight * rho_Hnuc_arr(i, j, k);

                            // max vel and Mach number
                            U_max_level = amrex::max(U_max_level, vel);
                            Mach_max_level =
                                amrex::max(Mach_max_level, vel / eos_state.cs);
                        }
                    }
                }
            }
        }  // end MFIter

        // sum quantities over all processors
        ParallelDescriptor::ReduceRealSum(T_center_level);
        ParallelDescriptor::ReduceRealSum(vel_center_level.dataPtr(),
                                          AMREX_SPACEDIM);
        ParallelDescriptor::ReduceRealSum(kin_ener_level);
        ParallelDescriptor::ReduceRealSum(int_ener_level);
        ParallelDescriptor::ReduceRealSum(nuc_ener_level);
        ParallelDescriptor::ReduceIntSum(ncenter_level);

        // find the largest U and Mach number over all processors
        ParallelDescriptor::ReduceRealMax(U_max_level);
        ParallelDescriptor::ReduceRealMax(Mach_max_level);

        // for T_max, we want to know where the hot spot is, so we do a
        // gather on the temperature and find the index corresponding to
        // the maxiumum.  We then pack the coordinates and velocities
        // into a local array and gather that to the I/O processor and
        // pick the values corresponding to the maximum.
        int nprocs = ParallelDescriptor::NProcs();
        int ioproc = ParallelDescriptor::IOProcessorNumber();
        Vector<Real> T_max_data(nprocs);

        if (nprocs == 1) {
            T_max_data[0] = T_max_local;
        } else {
            ParallelDescriptor::Gather(&T_max_local, 1, &T_max_data[0], 1,
                                       ioproc);
        }

        // determine index of max global T
        int index_max = 0;
        Real T_max_level = 0.0;
        for (int ip = 0; ip < nprocs; ++ip) {
            if (T_max_data[ip] > T_max_level) {
                T_max_level = T_max_data[ip];
                index_max = ip;
            }
        }

        // T_max_coords will contain both the coordinate information and
        // the velocity information, so there are 2*dm values on each
        // proc
        Vector<Real> T_max_coords(2 * AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            T_max_coords[i] = coord_Tmax_local[i];
            T_max_coords[i + AMREX_SPACEDIM] = vel_Tmax_local[i];
        }

        Vector<Real> T_max_coords_level(2 * AMREX_SPACEDIM * nprocs);

        if (nprocs == 1) {
            for (int i = 0; i < 2 * AMREX_SPACEDIM; ++i) {
                T_max_coords_level[i] = T_max_coords[i];
            }
        } else {
            ParallelDescriptor::Gather(&T_max_coords[0], 2 * AMREX_SPACEDIM,
                                       &T_max_coords_level[0],
                                       2 * AMREX_SPACEDIM, ioproc);
        }

        // initialize global variables
        Vector<Real> coord_Tmax_level(AMREX_SPACEDIM);
        Vector<Real> vel_Tmax_level(AMREX_SPACEDIM);

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            coord_Tmax_level[i] =
                T_max_coords_level[2 * AMREX_SPACEDIM * index_max + i];
            vel_Tmax_level[i] =
                T_max_coords_level[2 * AMREX_SPACEDIM * index_max + i +
                                   AMREX_SPACEDIM];
        }

        // for enuc_max, we also want to know where the hot spot is, so
        // we do the same gather procedure as with the temperature
        Vector<Real> enuc_max_data(nprocs);

        if (nprocs == 1) {
            enuc_max_data[0] = enuc_max_local;
        } else {
            ParallelDescriptor::Gather(&enuc_max_local, 1, &enuc_max_data[0], 1,
                                       ioproc);
        }

        // determine index of max global enuc
        index_max = 0;
        Real enuc_max_level = 0.0;
        for (int ip = 0; ip < nprocs; ++ip) {
            if (enuc_max_data[ip] > enuc_max_level) {
                enuc_max_level = enuc_max_data[ip];
                index_max = ip;
            }
        }

        Vector<Real> enuc_max_coords(2 * AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            enuc_max_coords[i] = coord_enucmax_local[i];
            enuc_max_coords[i + AMREX_SPACEDIM] = vel_enucmax_local[i];
        }

        Vector<Real> enuc_max_coords_level(2 * AMREX_SPACEDIM * nprocs);
        if (nprocs == 1) {
            for (int i = 0; i < 2 * AMREX_SPACEDIM; ++i) {
                enuc_max_coords_level[i] = enuc_max_coords[i];
            }
        } else {
            ParallelDescriptor::Gather(&enuc_max_coords[0], 2 * AMREX_SPACEDIM,
                                       &enuc_max_coords_level[0],
                                       2 * AMREX_SPACEDIM, ioproc);
        }

        // initialize global variables
        Vector<Real> coord_enucmax_level(AMREX_SPACEDIM);
        Vector<Real> vel_enucmax_level(AMREX_SPACEDIM);

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            coord_enucmax_level[i] =
                enuc_max_coords_level[2 * AMREX_SPACEDIM * index_max + i];
            vel_enucmax_level[i] =
                enuc_max_coords_level[2 * AMREX_SPACEDIM * index_max + i +
                                      AMREX_SPACEDIM];
        }

        // reduce the current level's data with the global data
        if (ParallelDescriptor::IOProcessor()) {
            kin_ener += kin_ener_level;
            int_ener += int_ener_level;
            nuc_ener += nuc_ener_level;

            U_max = amrex::max(U_max, U_max_level);
            Mach_max = amrex::max(Mach_max, Mach_max_level);

            // if T_max_level is the new max, then copy the location as well
            if (T_max_level > T_max) {
                T_max = T_max_level;

                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    coord_Tmax[i] = coord_Tmax_level[i];
                    vel_Tmax[i] = vel_Tmax_level[i];
                }

#if (AMREX_SPACEDIM == 3)
                if (spherical) {
                    // compute the radius of the bubble from the center of the star
                    Rloc_Tmax = std::sqrt((coord_Tmax[0] - center[0]) *
                                              (coord_Tmax[0] - center[0]) +
                                          (coord_Tmax[1] - center[1]) *
                                              (coord_Tmax[1] - center[1]) +
                                          (coord_Tmax[2] - center[2]) *
                                              (coord_Tmax[2] - center[2]));

                    // use the coordinates of the hot spot and the velocity components
                    // to compute the radial velocity at the hotspot
                    vr_Tmax =
                        ((coord_Tmax[0] - center[0]) / Rloc_Tmax) *
                            vel_Tmax[0] +
                        ((coord_Tmax[1] - center[1]) / Rloc_Tmax) *
                            vel_Tmax[1] +
                        ((coord_Tmax[2] - center[2]) / Rloc_Tmax) * vel_Tmax[2];
                }
#endif
            }

            // if enuc_max_level is the new max, then copy the location as well
            if (enuc_max_level > enuc_max) {
                enuc_max = enuc_max_level;

                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    coord_enucmax[i] = coord_enucmax_level[i];
                    vel_enucmax[i] = vel_enucmax_level[i];
                }

#if (AMREX_SPACEDIM == 3)
                if (spherical) {
                    // compute the radius of the bubble from the center
                    Rloc_enucmax =
                        std::sqrt((coord_enucmax[0] - center[0]) *
                                      (coord_enucmax[0] - center[0]) +
                                  (coord_enucmax[1] - center[1]) *
                                      (coord_enucmax[1] - center[1]) +
                                  (coord_enucmax[2] - center[2]) *
                                      (coord_enucmax[2] - center[2]));

                    // use the coordinates of the hot spot and the velocity components
                    // to compute the radial velocity at the hotspot
                    vr_enucmax =
                        ((coord_enucmax[0] - center[0]) / Rloc_enucmax) *
                            vel_enucmax[0] +
                        ((coord_enucmax[1] - center[1]) / Rloc_enucmax) *
                            vel_enucmax[1] +
                        ((coord_enucmax[2] - center[2]) / Rloc_enucmax) *
                            vel_enucmax[2];
                }
#endif
            }

            T_center += T_center_level;
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                vel_center[i] += vel_center_level[i];
            }
            ncenter += ncenter_level;
        }
    }

    // compute the graviational potential energy too
    Real grav_ener = 0.0;
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    if (spherical) {
#if (AMREX_SPACEDIM == 3)

        auto rho0 = rho0_in.array();

        // m(r) will contain mass enclosed by the center
        BaseState<Real> m_s(base_geom.nr_fine);
        auto m = m_s.array();
        m(0) = 4.0 / 3.0 * M_PI * rho0(0, 0) * r_cc_loc(0, 0) * r_cc_loc(0, 0) *
               r_cc_loc(0, 0);

        // dU = - G M dM / r;  dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
        grav_ener = -4.0 * M_PI * Gconst * m(0) * r_cc_loc(0, 0) * rho0(0, 0) *
                    (r_edge_loc(0, 1) - r_edge_loc(0, 0));

        for (auto r = 1; r < base_geom.nr_fine; ++r) {
            // the mass is defined at the cell-centers, so to compute the
            // mass at the current center, we need to add the contribution
            // of the upper half of the zone below us and the lower half of
            // the current zone.

            // don't add any contributions from outside the star -- i.e.
            // rho < base_cutoff_density
            Real term1 = 0.0;
            if (rho0(0, r - 1) > base_cutoff_density) {
                term1 = 4.0 / 3.0 * M_PI * rho0(0, r - 1) *
                        (r_edge_loc(0, r) - r_cc_loc(0, r - 1)) *
                        (r_edge_loc(0, r) * r_edge_loc(0, r) +
                         r_edge_loc(0, r) * r_cc_loc(0, r - 1) +
                         r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
            }

            Real term2 = 0.0;
            if (rho0(0, r) > base_cutoff_density) {
                term2 = 4.0 / 3.0 * M_PI * rho0(0, r) *
                        (r_cc_loc(0, r) - r_edge_loc(0, r)) *
                        (r_cc_loc(0, r) * r_cc_loc(0, r) +
                         r_cc_loc(0, r) * r_edge_loc(0, r) +
                         r_edge_loc(0, r) * r_edge_loc(0, r));
            }

            m(r) = m(r - 1) + term1 + term2;

            // dU = - G M dM / r;
            // dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
            grav_ener -= 4.0 * M_PI * Gconst * m(r) * r_cc_loc(0, r) *
                         rho0(0, r) * (r_edge_loc(0, r + 1) - r_edge_loc(0, r));
        }
#endif
    } else {
        const auto rho0 = rho0_in.const_array();
        // diag_grav_energy(&grav_ener, rho0_in.dataPtr(), r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            Real dr_loc = r_edge_loc(0, r + 1) - r_edge_loc(0, r);
            grav_ener -= rho0(0, r) * r_cc_loc(0, r) * grav_const * dr_loc;
        }
    }

    // normalize
    if (ParallelDescriptor::IOProcessor()) {
        // the volume we normalize with is that of a single coarse-level
        // zone.  This is because the weight used in the loop over cells
        // was with reference to the coarse level
        const Real* dx = geom[0].CellSize();
        for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
            kin_ener *= dx[i];
            int_ener *= dx[i];
            nuc_ener *= dx[i];
        }

        if (spherical) {
            // for a full star ncenter should be 8 -- there are only 8 zones
            // that have a vertex at the center of the star.  For an octant,
            // ncenter should be 1
            if (!((ncenter == 8 && !octant) || (ncenter == 1 && octant))) {
                Abort("ERROR: ncenter invalid in Diag()");
            } else {
                T_center /= ncenter;
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    vel_center[i] /= ncenter;
                }
            }
        }
    }

    // write out diagnosis data if at initialization
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename1 = "diag_temp.out";
        std::ofstream diagfile1;
        const std::string& diagfilename2 = "diag_enuc.out";
        std::ofstream diagfile2;
        const std::string& diagfilename3 = "diag_vel.out";
        std::ofstream diagfile3;

        // num of variables in the outfile depends on geometry but not dimension
        const int ndiag1 = (spherical) ? 11 : 8;
        const int ndiag2 = (spherical) ? 11 : 9;
        const int ndiag3 = (spherical) ? 10 : 7;

        if (step == 0) {
            // create file after initialization
            diagfile1.open(diagfilename1, std::ofstream::out |
                                              std::ofstream::trunc |
                                              std::ofstream::binary);

            // diag_temp.out
            // write variable names
            diagfile1 << std::setw(setwVal) << std::left << "time";
            diagfile1 << std::setw(setwVal) << std::left << "max{T}";
            diagfile1 << std::setw(setwVal) << std::left << "x(max{T})";
            diagfile1 << std::setw(setwVal) << std::left << "y(max{T})";
            diagfile1 << std::setw(setwVal) << std::left << "z(max{T})";
            diagfile1 << std::setw(setwVal) << std::left << "vx(max{T})";
            diagfile1 << std::setw(setwVal) << std::left << "vy(max{T})";
            diagfile1 << std::setw(setwVal) << std::left << "vz(max{T})";
            if (spherical) {
                diagfile1 << std::setw(setwVal) << std::left << "R(max{T})";
                diagfile1 << std::setw(setwVal) << std::left << "vr(max{T})";
                diagfile1 << std::setw(setwVal) << std::left << "T_center"
                          << std::endl;
            } else {
                diagfile1 << std::endl;
            }

            // write data
            diagfile1.precision(outfilePrecision);
            diagfile1 << std::scientific;
            diagfile1 << std::setw(setwVal) << std::left << t_in;
            diagfile1 << std::setw(setwVal) << std::left << T_max;

            const Real coord_temp_max_y = coord_Tmax[1];
            const Real coord_temp_max_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : coord_Tmax[2];
            diagfile1 << std::setw(setwVal) << std::left << coord_Tmax[0];
            diagfile1 << std::setw(setwVal) << std::left << coord_temp_max_y;
            diagfile1 << std::setw(setwVal) << std::left << coord_temp_max_z;

            const Real vel_temp_max_y = vel_Tmax[1];
            const Real vel_temp_max_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : vel_Tmax[2];
            diagfile1 << std::setw(setwVal) << std::left << vel_Tmax[0];
            diagfile1 << std::setw(setwVal) << std::left << vel_temp_max_y;
            diagfile1 << std::setw(setwVal) << std::left << vel_temp_max_z;

            if (spherical) {
                diagfile1 << std::setw(setwVal) << std::left << Rloc_Tmax;
                diagfile1 << std::setw(setwVal) << std::left << vr_Tmax;
                diagfile1 << std::setw(setwVal) << std::left << T_center
                          << std::endl;
            } else {
                diagfile1 << std::endl;
            }

            // close files
            diagfile1.close();

            // diag_enuc.out
            diagfile2.open(diagfilename2, std::ofstream::out |
                                              std::ofstream::trunc |
                                              std::ofstream::binary);
            // write variable names
            diagfile2 << std::setw(setwVal) << std::left << "time";
            diagfile2 << std::setw(setwVal) << std::left << "max{enuc}";
            diagfile2 << std::setw(setwVal) << std::left << "x(max{enuc})";
            diagfile2 << std::setw(setwVal) << std::left << "y(max{enuc})";
            diagfile2 << std::setw(setwVal) << std::left << "z(max{enuc})";
            diagfile2 << std::setw(setwVal) << std::left << "vx(max{enuc})";
            diagfile2 << std::setw(setwVal) << std::left << "vy(max{enuc})";
            diagfile2 << std::setw(setwVal) << std::left << "vz(max{enuc})";
            if (spherical) {
                diagfile2 << std::setw(setwVal) << std::left << "R(max{enuc})";
                diagfile2 << std::setw(setwVal) << std::left << "vr(max{enuc})";
            }
            diagfile2 << std::setw(setwVal) << std::left
                      << "tot nuc ener(erg/s)" << std::endl;

            // write data
            diagfile2.precision(outfilePrecision);
            diagfile2 << std::scientific;
            diagfile2 << std::setw(setwVal) << std::left << t_in;
            diagfile2 << std::setw(setwVal) << std::left << enuc_max;

            const Real coord_enuc_y = coord_enucmax[1];
            const Real coord_enuc_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : coord_enucmax[2];
            diagfile2 << std::setw(setwVal) << std::left << coord_enucmax[0];
            diagfile2 << std::setw(setwVal) << std::left << coord_enuc_y;
            diagfile2 << std::setw(setwVal) << std::left << coord_enuc_z;

            const Real vel_enuc_y = vel_enucmax[1];
            const Real vel_enuc_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : vel_enucmax[2];
            diagfile2 << std::setw(setwVal) << std::left << vel_enucmax[0];
            diagfile2 << std::setw(setwVal) << std::left << vel_enuc_y;
            diagfile2 << std::setw(setwVal) << std::left << vel_enuc_z;

            if (spherical) {
                diagfile2 << std::setw(setwVal) << std::left << Rloc_enucmax;
                diagfile2 << std::setw(setwVal) << std::left << vr_enucmax;
            }
            diagfile2 << std::setw(setwVal) << std::left << nuc_ener
                      << std::endl;

            // close file
            diagfile2.close();

            // diag_vel.out
            diagfile3.open(diagfilename3, std::ofstream::out |
                                              std::ofstream::trunc |
                                              std::ofstream::binary);
            // write variable names
            diagfile3 << std::setw(setwVal) << std::left << "time";
            diagfile3 << std::setw(setwVal) << std::left << "max{U}";
            diagfile3 << std::setw(setwVal) << std::left << "max{Mach}";
            diagfile3 << std::setw(setwVal) << std::left << "tot kin energy";
            diagfile3 << std::setw(setwVal) << std::left << "tot grav energy";
            diagfile3 << std::setw(setwVal) << std::left << "tot int energy";
            if (spherical) {
                diagfile3 << std::setw(setwVal) << std::left << "velx_center";
                diagfile3 << std::setw(setwVal) << std::left << "vely_center";
                diagfile3 << std::setw(setwVal) << std::left << "velz_center";
            }
            diagfile3 << std::setw(setwVal) << std::left << "dt" << std::endl;

            // write data
            diagfile3.precision(outfilePrecision);
            diagfile3 << std::scientific;
            diagfile3 << std::setw(setwVal) << std::left << t_in;
            diagfile3 << std::setw(setwVal) << std::left << U_max;
            diagfile3 << std::setw(setwVal) << std::left << Mach_max;
            diagfile3 << std::setw(setwVal) << std::left << kin_ener;
            diagfile3 << std::setw(setwVal) << std::left << grav_ener;
            diagfile3 << std::setw(setwVal) << std::left << int_ener;
            if (spherical) {
                diagfile3 << std::setw(setwVal) << std::left << vel_center[0];
                diagfile3 << std::setw(setwVal) << std::left << vel_center[1];
                diagfile3 << std::setw(setwVal) << std::left << vel_center[2];
            }
            diagfile3 << std::setw(setwVal) << std::left << dt << std::endl;

            // close file
            diagfile3.close();

        } else {
            // store variable values in data array to be written later

            // temp
            diagfile1_data[index * ndiag1] = t_in;
            diagfile1_data[index * ndiag1 + 1] = T_max;
            const Real coord_temp_max_y = coord_Tmax[1];
            const Real coord_temp_max_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : coord_Tmax[2];
            diagfile1_data[index * ndiag1 + 2] = coord_Tmax[0];
            diagfile1_data[index * ndiag1 + 3] = coord_temp_max_y;
            diagfile1_data[index * ndiag1 + 4] = coord_temp_max_z;
            const Real vel_temp_max_y = vel_Tmax[1];
            const Real vel_temp_max_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : vel_Tmax[2];
            diagfile1_data[index * ndiag1 + 5] = vel_Tmax[0];
            diagfile1_data[index * ndiag1 + 6] = vel_temp_max_y;
            diagfile1_data[index * ndiag1 + 7] = vel_temp_max_z;
            if (spherical) {
                diagfile1_data[index * ndiag1 + 8] = Rloc_Tmax;
                diagfile1_data[index * ndiag1 + 9] = vr_Tmax;
                diagfile1_data[index * ndiag1 + 10] = T_center;
            }

            // enuc
            diagfile2_data[index * ndiag2] = t_in;
            diagfile2_data[index * ndiag2 + 1] = enuc_max;
            const Real coord_enuc_y = coord_enucmax[1];
            const Real coord_enuc_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : coord_enucmax[2];
            diagfile2_data[index * ndiag2 + 2] = coord_enucmax[0];
            diagfile2_data[index * ndiag2 + 3] = coord_enuc_y;
            diagfile2_data[index * ndiag2 + 4] = coord_enuc_z;
            const Real vel_enuc_y = vel_enucmax[1];
            const Real vel_enuc_z =
                (AMREX_SPACEDIM == 2) ? 0.0 : vel_enucmax[2];
            diagfile2_data[index * ndiag2 + 5] = vel_enucmax[0];
            diagfile2_data[index * ndiag2 + 6] = vel_enuc_y;
            diagfile2_data[index * ndiag2 + 7] = vel_enuc_z;
            diagfile2_data[index * ndiag2 + 8] = Rloc_enucmax;
            diagfile2_data[index * ndiag2 + 9] = vr_enucmax;
            diagfile2_data[index * ndiag2 + 10] = nuc_ener;

            // vel
            diagfile3_data[index * ndiag3] = t_in;
            diagfile3_data[index * ndiag3 + 1] = U_max;
            diagfile3_data[index * ndiag3 + 2] = Mach_max;
            diagfile3_data[index * ndiag3 + 3] = kin_ener;
            diagfile3_data[index * ndiag3 + 4] = grav_ener;
            diagfile3_data[index * ndiag3 + 5] = int_ener;
            if (spherical) {
                diagfile3_data[index * ndiag3 + 6] = vel_center[0];
                diagfile3_data[index * ndiag3 + 7] = vel_center[1];
                diagfile3_data[index * ndiag3 + 8] = vel_center[2];
            }
            const int idt = spherical ? 9 : 6;
            diagfile3_data[index * ndiag3 + idt] = dt;

            index += 1;
        }
    }  // } IOProcessor
}

// put together a vector of multifabs for writing
void Maestro::WriteDiagFile(int& index) {
    // num of variables in the outfile depends on geometry but not dimension
    const int ndiag1 = (spherical) ? 11 : 8;
    const int ndiag2 = (spherical) ? 11 : 9;
    const int ndiag3 = (spherical) ? 10 : 7;

    // timer for profiling
    BL_PROFILE_VAR("Maestro::WriteDiagFile()", WriteDiagFile);

    // write out diagnosis data
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename1 = "diag_temp.out";
        std::ofstream diagfile1(diagfilename1, std::ofstream::out |
                                                   std::ofstream::app |
                                                   std::ofstream::binary);
        // time
        // T_max
        // coord_Tmax (3)
        // vel_Tmax (3)
        // plus, if spherical:
        // -- Rloc_Tmax
        // -- vr_Tmax
        // -- T_center
        diagfile1.precision(outfilePrecision);
        diagfile1 << std::scientific;
        for (auto i = 0; i < index; ++i) {
            for (auto comp = 0; comp < ndiag1; ++comp) {
                diagfile1 << std::setw(setwVal) << std::left
                          << diagfile1_data[i * ndiag1 + comp];
            }
            diagfile1 << std::endl;
        }

        // close file
        diagfile1.close();

        const std::string& diagfilename2 = "diag_enuc.out";
        std::ofstream diagfile2(diagfilename2, std::ofstream::out |
                                                   std::ofstream::app |
                                                   std::ofstream::binary);
        // time
        // enuc_max
        // coord_enucmax (3)
        // vel_enucmax (3)
        // (if spherical):
        // -- Rloc_enucmax
        // -- vr_enucmax
        // nuc_ener
        diagfile2.precision(outfilePrecision);
        diagfile2 << std::scientific;
        for (auto i = 0; i < index; ++i) {
            for (auto comp = 0; comp < ndiag2; ++comp) {
                diagfile2 << std::setw(setwVal) << std::left
                          << diagfile2_data[i * ndiag2 + comp];
            }
            diagfile2 << std::endl;
        }

        // close file
        diagfile2.close();

        const std::string& diagfilename3 = "diag_vel.out";
        std::ofstream diagfile3(diagfilename3, std::ofstream::out |
                                                   std::ofstream::app |
                                                   std::ofstream::binary);
        // time
        // U_max
        // Mach_max
        // kin_ener
        // grav_ener
        // int_ener
        // vel_center (3) (only if spherical)
        // dt
        diagfile3.precision(outfilePrecision);
        diagfile3 << std::scientific;
        for (auto i = 0; i < index; ++i) {
            for (auto comp = 0; comp < ndiag3; ++comp) {
                diagfile3 << std::setw(setwVal) << std::left
                          << diagfile3_data[i * ndiag3 + comp];
            }
            diagfile3 << std::endl;
        }

        // close file
        diagfile3.close();

        // reset buffer array
        index = 0;
    }
}
