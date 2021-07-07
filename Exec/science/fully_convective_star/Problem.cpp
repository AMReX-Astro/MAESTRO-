#include "Maestro.H"

using namespace amrex;

#ifdef DO_PROBLEM_POST_INIT

void Maestro::ProblemPostInit() {
    BL_PROFILE_VAR("Maestro::ProblemPostInit()", ProblemPostInit);

    // initialize radial diagnostic
    // 0: convective velocity
    // 1: meridional circulation (radial)
    postfile_data.resize(base_geom.max_radial_level + 1, base_geom.nr_fine,
                         NumPost);
    postfile_data.setVal(0.);

    // initialize multifab diagnostics

    // need single-level finest grid for whole domain
    const Box& domain0 = geom[0].Domain();
    auto dom0_lo = domain0.loVect();
    auto dom0_hi = domain0.hiVect();
    IntVect ncell(0, 1, 0);
    ncell[0] = (dom0_hi[0] + 1) * (finest_level + 1) / 2 - 1;
    ncell[2] = (dom0_hi[2] + 1) * (finest_level + 1) - 1;

    IntVect domLo(AMREX_D_DECL(dom0_lo[0], dom0_lo[1], dom0_lo[2]));
    IntVect domHi(AMREX_D_DECL(ncell[0], ncell[1], ncell[2]));
    Box domain(domLo, domHi);

    // we only need half x-z plane because theta is between (0,pi)
    auto prob_lo = geom[0].ProbLo();
    auto prob_hi = geom[0].ProbHi();
    const auto dx = geom[finest_level].CellSizeArray();

    Real probLo[] = {prob_lo[0] + (prob_hi[0] - prob_lo[0]) / 2.0, prob_lo[1],
                     prob_lo[2]};
    Real probHi[] = {prob_hi[0], prob_lo[1] + 2.0 * dx[1], prob_hi[2]};
    RealBox real_box(probLo, probHi);

    // get max_grid_size
    ParmParse pp("amr");
    int max_grid_size;
    pp.get("max_grid_size", max_grid_size);

    // make BoxArray and Geometry for single-level finest grid
    // previously declared in Maestro.H
    baFine.define(domain);
    baFine.maxSize(max_grid_size);
    dmFine.define(baFine);
    geomFine.define(domain, &real_box);

    // define multifab velocities on new single-level grid
    postfile_mf.resize(1);
    postfile_mf[0].define(baFine, dmFine, 2, 0);
    postfile_mf[0].setVal(0.);

    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> ang_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rad_vel[lev].define(grids[lev], dmap[lev], 1, 0);
        ang_vel[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    BaseState<Real> convect_vel(base_geom.max_radial_level + 1,
                                base_geom.nr_fine);
    Vector<MultiFab> circ_vel(1);
    circ_vel[0].define(baFine, dmFine, 2, 0);

    // initialize diagnostics
    if (restart_file.empty()) {
        // compute diagnostics for initial state
        // make radial and theta-component velocity
        MakeVelrth(uold, w0_cart, rad_vel, ang_vel);

        // make convection velocity
        MakeConvectionVel(rad_vel, convect_vel, postfile_data, 0, t_new);

        // make meridional circulation
        MakeMeridionalCirculation(rad_vel, ang_vel, circ_vel, postfile_mf,
                                  t_new);

    } else {
        // if restarting from checkpoint
        Print() << "Initializing postprocessing diagnostics from checkpoint "
                << restart_file << std::endl;

        std::string line, word;

        // read in cell-centered basestate
        std::string File(restart_file + "/PostCC");
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        for (int i = 0;
             i < (base_geom.max_radial_level + 1) * base_geom.nr_fine; ++i) {
            std::getline(is, line);
            std::istringstream lis(line);

            // NumPost components in postfile_data
            for (int ncomp = 0; ncomp < NumPost; ++ncomp) {
                lis >> word;
                postfile_data.array()(ncomp + 2 * i) = std::stod(word);
            }
        }

        // read in the multifab data
        VisMF::Read(postfile_mf[0],
                    amrex::MultiFabFileFullPrefix(0, restart_file,
                                                  "PostTimestep_", "diag"));

        // compute diagnostics output
        // set updated values to zero to get old values at t_old
        for (int lev = 0; lev <= finest_level; ++lev) {
            rad_vel[lev].setVal(0.);
            ang_vel[lev].setVal(0.);
        }

        // make convection velocity
        MakeConvectionVel(rad_vel, convect_vel, postfile_data, 0, t_old);

        // make meridional circulation
        MakeMeridionalCirculation(rad_vel, ang_vel, circ_vel, postfile_mf,
                                  t_old);
    }

    if (restart_file.empty()) {
        // wrote a plotfile
        if (plot_int > 0 || plot_deltat > 0) {
            std::string plotfilename = plot_base_name + "0000000";

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            // write out the cell-centered diagnostics
            if (ParallelDescriptor::IOProcessor()) {
                for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
                    std::ofstream PostFile;
                    PostFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
                                                io_buffer.size());
                    std::string PostFileName(plotfilename + "/PostTimestep_");
                    std::string levStr = std::to_string(lev);
                    PostFileName.append(levStr);
                    PostFile.open(PostFileName.c_str(),
                                  std::ofstream::out | std::ofstream::trunc |
                                      std::ofstream::binary);
                    if (!PostFile.good()) {
                        amrex::FileOpenFailed(PostFileName);
                    }

                    PostFile.precision(17);

                    PostFile << "r_cc  convection_vel   \n";

                    for (int i = 0; i < base_geom.nr(lev); ++i) {
                        PostFile << base_geom.r_cc_loc(lev, i) << " "
                                 << convect_vel.array()(lev, i) << "\n";
                    }
                }
            }

            // write out the planar-averaged diagnostics
            std::string PostFileName(plotfilename + "/PostTimestep");
            WriteSingleLevelPlotfile(PostFileName, circ_vel[0],
                                     {"vel_radial", "vel_theta"}, geomFine, 0,
                                     0);
        }

        // wrote a checkpoint file
        if (chk_int > 0 || chk_deltat > 0) {
            // checkpoint file name, e.g., chk00010
            const std::string& checkpointname =
                amrex::Concatenate(check_base_name, 0, 7);

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            // write out the cell-centered base state
            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream PostCCFile;
                PostCCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
                                              io_buffer.size());
                std::string PostCCFileName(checkpointname + "/PostCC");
                PostCCFile.open(PostCCFileName.c_str(),
                                std::ofstream::out | std::ofstream::trunc |
                                    std::ofstream::binary);
                if (!PostCCFile.good()) {
                    amrex::FileOpenFailed(PostCCFileName);
                }

                PostCCFile.precision(17);

                for (int i = 0;
                     i < (base_geom.max_radial_level + 1) * base_geom.nr_fine;
                     ++i) {
                    // NumPost components in postfile_data
                    for (int ncomp = 0; ncomp < NumPost; ++ncomp) {
                        PostCCFile << postfile_data.array()(ncomp + 2 * i)
                                   << " ";
                    }
                    PostCCFile << "\n";
                }
            }

            // write the MultiFab data
            amrex::UtilCreateCleanDirectory(checkpointname + "/PostTimestep_0",
                                            true);
            VisMF::Write(postfile_mf[0],
                         amrex::MultiFabFileFullPrefix(
                             0, checkpointname, "PostTimestep_", "diag"));
        }
    }
}

#endif

#ifdef DO_PROBLEM_POST_TIMESTEP

void Maestro::ProblemPostTimestep() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ProblemPostInit()", ProblemPostInit);

    // compute diagnostics after every time step
    // make radial and theta-component velocity
    Vector<MultiFab> rad_vel(finest_level + 1);
    Vector<MultiFab> ang_vel(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rad_vel[lev].define(grids[lev], dmap[lev], 1, 0);
        ang_vel[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    MakeVelrth(unew, w0_cart, rad_vel, ang_vel);

    // make convection velocity
    BaseState<Real> convect_vel(base_geom.max_radial_level + 1,
                                base_geom.nr_fine);
    MakeConvectionVel(rad_vel, convect_vel, postfile_data, 0, t_new);

    // make meridional circulation
    Vector<MultiFab> circ_vel(1);
    circ_vel[0].define(baFine, dmFine, 2, 0);
    MakeMeridionalCirculation(rad_vel, ang_vel, circ_vel, postfile_mf, t_new);

    // wrote a plotfile
    if ((plot_int > 0 && istep % plot_int == 0) ||
        (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
        ((plot_int > 0 || plot_deltat > 0) &&
         (istep == max_step || t_old >= stop_time))) {
        std::string plotfilename = plot_base_name;
        PlotFileName(istep, &plotfilename);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        // write out the cell-centered diagnostics
        if (ParallelDescriptor::IOProcessor()) {
            for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
                std::ofstream PostFile;
                PostFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
                                            io_buffer.size());
                std::string PostFileName(plotfilename + "/PostTimestep_");
                std::string levStr = std::to_string(lev);
                PostFileName.append(levStr);
                PostFile.open(PostFileName.c_str(), std::ofstream::out |
                                                        std::ofstream::trunc |
                                                        std::ofstream::binary);
                if (!PostFile.good()) {
                    amrex::FileOpenFailed(PostFileName);
                }

                PostFile.precision(17);

                PostFile << "r_cc  convection_vel   \n";

                for (int i = 0; i < base_geom.nr(lev); ++i) {
                    PostFile << base_geom.r_cc_loc(lev, i) << " "
                             << convect_vel.array()(lev, i) << "\n";
                }
            }
        }

        // write out the planar-averaged diagnostics
        std::string PostFileName(plotfilename + "/PostTimestep");
        WriteSingleLevelPlotfile(PostFileName, circ_vel[0],
                                 {"vel_radial", "vel_theta"}, geomFine, 0, 0);
    }

    // wrote a checkpoint file
    if ((chk_int > 0 && istep % chk_int == 0) ||
        (chk_deltat > 0 && std::fmod(t_new, chk_deltat) < dt) ||
        ((chk_int > 0 || chk_deltat > 0) &&
         (istep == max_step || t_old >= stop_time))) {
        // checkpoint file name, e.g., chk00010
        const std::string& checkpointname =
            amrex::Concatenate(check_base_name, istep, 7);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        // write out the cell-centered base state
        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream PostCCFile;
            PostCCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(),
                                          io_buffer.size());
            std::string PostCCFileName(checkpointname + "/PostCC");
            PostCCFile.open(PostCCFileName.c_str(), std::ofstream::out |
                                                        std::ofstream::trunc |
                                                        std::ofstream::binary);
            if (!PostCCFile.good()) {
                amrex::FileOpenFailed(PostCCFileName);
            }

            PostCCFile.precision(17);

            for (int i = 0;
                 i < (base_geom.max_radial_level + 1) * base_geom.nr_fine;
                 ++i) {
                // NumPost components in postfile_data
                for (int ncomp = 0; ncomp < NumPost; ++ncomp) {
                    PostCCFile << postfile_data.array()(ncomp + 2 * i) << " ";
                }
                PostCCFile << "\n";
            }
        }

        // write the MultiFab data
        amrex::UtilCreateCleanDirectory(checkpointname + "/PostTimestep_0",
                                        true);
        VisMF::Write(postfile_mf[0],
                     amrex::MultiFabFileFullPrefix(0, checkpointname,
                                                   "PostTimestep_", "diag"));
    }
}

#endif

void Maestro::MakeConvectionVel(const Vector<MultiFab>& velr,
                                BaseState<Real>& vel_conv,
                                BaseState<Real>& totsum, const int comp,
                                const Real t_interval) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeConvectionVel()", MakeConvectionVel);

    Vector<MultiFab> vel2(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel2[lev].define(grids[lev], dmap[lev], 1, 0);
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
    AverageDown(vel2, 0, 1);
    // FillPatch(t0, vel2, vel2, vel2, 0, 0, 1, 0, bcs_f);

    // radial average of square of radial vel
    Average(vel2, vel_conv, 0);

    // add vel^2 to total sum and compute convection vel
    auto vel_conv_arr = vel_conv.array();
    auto totsum_arr = totsum.array();

    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        // Zero the velocity if it is a negative value.
        // Note that the Average() subroutine can cause some interpolation
        // values near the center of star to be negative.
        // This should not affect values located further from the center.
        if (vel_conv_arr(0, r) < 0) {
            vel_conv_arr(0, r) = 0.0;
        }

        if (t_new - dt >= diag_skip_time) {
            totsum_arr(0, r, comp) += dt * vel_conv_arr(0, r);
        } else if (t_new >= diag_skip_time) {
            totsum_arr(0, r, comp) +=
                (t_new - diag_skip_time) * vel_conv_arr(0, r);
        } else {
            // do nothing
        }

        // root-mean-squared radial velocity
        vel_conv_arr(0, r) =
            std::sqrt(totsum_arr(0, r, comp) / (t_interval - diag_skip_time));
    }
}

void Maestro::MakeMeridionalCirculation(const Vector<MultiFab>& velr,
                                        const Vector<MultiFab>& velth,
                                        Vector<MultiFab>& vel_circ,
                                        Vector<MultiFab>& sum_circ,
                                        const Real t_interval) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeMeridionalCirculation()",
                   MakeMeridionalCirculation);

    vel_circ[0].setVal(0.);

    // average to 2D r-theta (x-z) plane
    Average2d(velr, vel_circ, 0, 0, {geomFine.Domain()});
    Average2d(velth, vel_circ, 1, 0, {geomFine.Domain()});

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(vel_circ[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get the index space of the valid region
        const Box& tileBox = mfi.tilebox();

        const Array4<Real> circvel_arr = vel_circ[0].array(mfi);
        const Array4<Real> circsum_arr = sum_circ[0].array(mfi);
        const auto dt_loc = dt;
        const auto tnew_loc = t_new;

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            for (int n = 0; n < 2; ++n) {
                if (tnew_loc - dt_loc >= diag_skip_time) {
                    circsum_arr(i, j, k, n) += dt_loc * circvel_arr(i, j, k, n);
                } else if (tnew_loc >= diag_skip_time) {
                    circsum_arr(i, j, k, n) +=
                        (tnew_loc - diag_skip_time) * circvel_arr(i, j, k, n);
                } else {
                    // do nothing
                }

                circvel_arr(i, j, k, n) =
                    circsum_arr(i, j, k, n) / (t_interval - diag_skip_time);
            }
        });
    }
}
