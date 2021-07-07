
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()", Evolve);

    Print() << "Calling Evolve()" << std::endl;

    // check to make sure spherical is only used for 3d
    if (spherical && AMREX_SPACEDIM != 3) {
        Abort("spherical = 1 and dm != 3");
    }

    // index for diag array buffer
    int diag_index = 0;

    for (istep = start_step; ((istep <= max_step || max_step < 0) &&
                              (t_old < stop_time || stop_time < 0.0));
         ++istep) {
        // check to see if we need to regrid, then regrid
        if (max_level > 0 && regrid_int > 0 && (istep - 1) % regrid_int == 0 &&
            istep != 1) {
            Regrid();
        }

        dtold = dt;

        // compute time step
        // if this is the first time step we already have a dt from either FirstDt()
        // or EstDt called during the divu_iters
        if (istep > 1) {
            EstDt();

            if (maestro_verbose > 0) {
                Print() << "Call to estdt at beginning of step " << istep
                        << " gives dt =" << dt << std::endl;
            }

            // fixme - add nuclear_dt_scalefac timestep limiter

            if (dt > max_dt_growth * dtold) {
                dt = max_dt_growth * dtold;
                if (maestro_verbose > 0) {
                    Print() << "dt_growth factor limits the new dt = " << dt
                            << std::endl;
                }
            }

            if (dt > max_dt) {
                if (maestro_verbose > 0) {
                    Print()
                        << "max_dt limits the new dt = " << max_dt << std::endl;
                }
                dt = max_dt;
            }

            if (fixed_dt != -1.) {
                dt = fixed_dt;
                if (maestro_verbose > 0) {
                    Print() << "Setting fixed dt = " << dt;
                }
            }

            if (stop_time >= 0. && t_old + dt > stop_time) {
                dt = amrex::min(dt, stop_time - t_old);
                Print() << "Stop time limits dt = " << dt << std::endl;
            }

            t_new = t_old + dt;
        }

        // wallclock time
        Real start_total = ParallelDescriptor::second();

        // advance the solution by dt
#ifdef SDC
        AdvanceTimeStepSDC(false);
#else
        if (use_exact_base_state || average_base_state) {
            // new temporal algorithm
            AdvanceTimeStepAverage(false);
        } else {
            // original temporal algorithm
            AdvanceTimeStep(false);
        }
#endif

        t_old = t_new;

        if ((sum_interval > 0 && istep % sum_interval == 0) ||
            (sum_per > 0 && std::fmod(t_new, sum_per) < dt) ||
            ((sum_interval > 0 || sum_per > 0) && t_old >= stop_time)) {
            Real diag_start_total = ParallelDescriptor::second();

            // save diag output into buffer
            DiagFile(istep, t_new, rho0_new, p0_new, unew, snew, diag_index);

            // wallclock time
            Real diag_end_total =
                ParallelDescriptor::second() - diag_start_total;
            ParallelDescriptor::ReduceRealMax(
                diag_end_total, ParallelDescriptor::IOProcessorNumber());

            Print() << "Diagnostic :" << diag_end_total << " seconds\n\n";
        }

        Real end_total = ParallelDescriptor::second() - start_total;
        ParallelDescriptor::ReduceRealMax(
            end_total, ParallelDescriptor::IOProcessorNumber());

        Print() << "Time to advance time step: " << end_total << '\n';

        if ((plot_int > 0 && istep % plot_int == 0) ||
            (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
            ((plot_int > 0 || plot_deltat > 0) &&
             (istep == max_step || t_old >= stop_time))) {
            // write a plotfile
            Print() << "\nWriting plotfile " << istep << std::endl;
            WritePlotFile(istep, t_new, dt, rho0_new, rhoh0_new, p0_new,
                          gamma1bar_new, unew, snew, S_cc_new);
        }

        if ((small_plot_int > 0 && istep % small_plot_int == 0) ||
            (small_plot_deltat > 0 &&
             std::fmod(t_new, small_plot_deltat) < dt) ||
            ((small_plot_int > 0 || small_plot_deltat > 0) &&
             (istep == max_step || t_old >= stop_time))) {
            // write a small plotfile
            Print() << "\nWriting small plotfile " << istep << std::endl;
            WriteSmallPlotFile(istep, t_new, dt, rho0_new, rhoh0_new, p0_new,
                               gamma1bar_new, unew, snew, S_cc_new);
        }

        if ((chk_int > 0 && istep % chk_int == 0) ||
            (chk_deltat > 0 && std::fmod(t_new, chk_deltat) < dt) ||
            ((chk_int > 0 || chk_deltat > 0) &&
             (istep == max_step || t_old >= stop_time))) {
            // write a checkpoint file
            Print() << "\nWriting checkpoint " << istep << std::endl;
            WriteCheckPoint(istep);
        }

        if ((diag_index == diag_buf_size || istep == max_step ||
             t_old >= stop_time) &&
            (sum_per > 0.0 || sum_interval > 0)) {
            // write out any buffered diagnostic information
            WriteDiagFile(diag_index);
        }

#ifdef DO_PROBLEM_POST_TIMESTEP

        // Provide a hook for the user to do things after all of
        // the normal updates have been applied. The user is
        // responsible for any actions after this point, like
        // doing a computeTemp call if they change the state data.

        ProblemPostTimestep();

#endif

        // move new state into old state by swapping pointers
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(sold[lev], snew[lev]);
            std::swap(uold[lev], unew[lev]);
            std::swap(S_cc_old[lev], S_cc_new[lev]);
        }

        rho0_old.swap(rho0_new);
        rhoh0_old.swap(rhoh0_new);
        p0_nm1.swap(p0_old);
        p0_old.swap(p0_new);

        beta0_old.swap(beta0_new);
        gamma1bar_old.swap(gamma1bar_new);
        grav_cell_old.swap(grav_cell_new);
    }
}
