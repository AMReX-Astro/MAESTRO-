
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <Maestro_F.H>
#include <model_parser_F.H>
using namespace amrex;

// initialize AMR data
// perform initial projection
// perform divu iters
// perform initial (pressure) iterations
void Maestro::Init() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Init()", Init);

    Print() << "Calling Init()" << std::endl;

    if (restart_file.empty()) {
        start_step = 1;

        // fill in multifab and base state data
        InitData();

        if (plot_int > 0 || plot_deltat > 0) {
            // Need to fill normal vector to compute velrc in plotfile
            if (spherical) {
                MakeNormal();
            }

            Print() << "\nWriting plotfile " << plot_base_name
                    << "InitData after InitData" << std::endl;

            WritePlotFile(plotInitData, t_old, 0, rho0_old, rhoh0_old, p0_old,
                          gamma1bar_old, uold, sold, S_cc_old);
        } else if (small_plot_int > 0 || small_plot_deltat > 0) {
            // Need to fill normal vector to compute velrc in plotfile
            if (spherical) {
                MakeNormal();
            }

            Print() << "\nWriting small plotfile " << small_plot_base_name
                    << "InitData after InitData" << std::endl;
            WriteSmallPlotFile(plotInitData, t_old, 0, rho0_old, rhoh0_old,
                               p0_old, gamma1bar_old, uold, sold, S_cc_old);
        }

    } else {
        Print() << "Initializing from checkpoint " << restart_file << std::endl;

        input_model.ReadFile(model_file);

        // read in checkpoint file
        // this builds (defines) and fills the following MultiFabs:
        //
        // snew, unew, gpi, dSdt, S_cc_new
        //
        // and also fills in the 1D arrays:
        //
        // rho0_new, p0_new, gamma1bar_new, rhoh0_new, beta0_new, psi, tempbar, etarho_cc, tempbar_init
        ReadCheckPoint();

        // build (define) the following MultiFabs (that weren't read in from checkpoint):
        // snew, unew, S_cc_new, w0_cart, rhcc_for_nodalproj, normal, pi
        for (int lev = 0; lev <= finest_level; ++lev) {
            snew[lev].define(grids[lev], dmap[lev], Nscal, ng_s);
            unew[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
            S_cc_new[lev].define(grids[lev], dmap[lev], 1, 0);
            w0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 2);
            rhcc_for_nodalproj[lev].define(grids[lev], dmap[lev], 1, 1);
            if (spherical) {
                normal[lev].define(grids[lev], dmap[lev], 3, 1);
                cell_cc_to_r[lev].define(grids[lev], dmap[lev], 1, 0);
            }
            pi[lev].define(convert(grids[lev], nodal_flag), dmap[lev], 1,
                           0);  // nodal
#ifdef SDC
            intra[lev].define(grids[lev], dmap[lev], Nscal, 0);  // for sdc
            intra[lev].setVal(0.);
#endif
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            w0_cart[lev].setVal(0.);
            rhcc_for_nodalproj[lev].setVal(0.);
            pi[lev].setVal(0.);
            S_cc_new[lev].setVal(0.);
            unew[lev].setVal(0.);
            snew[lev].setVal(0.);
        }

        // put w0 on Cartesian cell-centers
        if (evolve_base_state) {
            Put1dArrayOnCart(w0, w0_cart, true, true, bcs_u, 0, 1);
        }

        if (!spherical) {
            // reset tagging array to include buffer zones
            TagArray();
        }

        // set finest_radial_level
        // compute numdisjointchunks, r_start_coord, r_end_coord
        BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                                   base_geom.nr_fine);
        base_geom.InitMultiLevel(finest_level, tag_array_b.array());
        base_geom.ComputeCutoffCoords(rho0_old.array());
    }

#if (AMREX_SPACEDIM == 3)
    if (spherical) {
        MakeNormal();
        MakeCCtoRadii();
    }
#endif

    if (do_sponge) {
        SpongeInit(rho0_old);
    }

    // make gravity
    MakeGravCell(grav_cell_old, rho0_old);

    if (restart_file.empty()) {
        // compute gamma1bar
        MakeGamma1bar(sold, gamma1bar_old, p0_old);

        // compute beta0
        MakeBeta0(beta0_old, rho0_old, p0_old, gamma1bar_old, grav_cell_old);

        // set beta0^{-1} = beta0_old
        beta0_nm1.copy(beta0_old);

        // initial projection
        if (do_initial_projection) {
            Print() << "Doing initial projection" << std::endl;
            InitProj();

            if (plot_int > 0 || plot_deltat > 0) {
                Print() << "\nWriting plotfile " << plot_base_name
                        << "after_InitProj after InitProj" << std::endl;

                WritePlotFile(plotInitProj, t_old, 0, rho0_old, rhoh0_old,
                              p0_old, gamma1bar_old, uold, sold, S_cc_old);

            } else if (small_plot_int > 0 || small_plot_deltat > 0) {
                Print() << "\nWriting small plotfile " << small_plot_base_name
                        << "after_InitProj after InitProj" << std::endl;

                WriteSmallPlotFile(plotInitProj, t_old, 0, rho0_old, rhoh0_old,
                                   p0_old, gamma1bar_old, uold, sold, S_cc_old);
            }
        }

        // compute initial time step
        FirstDt();

        // divu iters - also update dt at end of each divu_iter
        if (init_divu_iter > 0) {
            for (int i = 1; i <= init_divu_iter; ++i) {
                Print() << "Doing initial divu iteration #" << i << std::endl;
#ifdef SDC
                DivuIterSDC(i);
#else
                DivuIter(i);
#endif
            }

            if (plot_int > 0 || plot_deltat > 0) {
                Print() << "\nWriting plotfile " << plot_base_name
                        << "after_DivuIter after final DivuIter" << std::endl;
                WritePlotFile(plotDivuIter, t_old, dt, rho0_old, rhoh0_old,
                              p0_old, gamma1bar_old, uold, sold, S_cc_old);
            } else if (small_plot_int > 0 || small_plot_deltat > 0) {
                Print() << "\nWriting small plotfile " << small_plot_base_name
                        << "after_DivuIter after final DivuIter" << std::endl;
                WriteSmallPlotFile(plotDivuIter, t_old, dt, rho0_old, rhoh0_old,
                                   p0_old, gamma1bar_old, uold, sold, S_cc_old);
            }
        }

        if (stop_time >= 0. && t_old + dt > stop_time) {
            dt = amrex::min(dt, stop_time - t_old);
            Print() << "Stop time limits dt = " << dt << std::endl;
        }

        dtold = dt;
        t_new = t_old + dt;

        // copy S_cc_old into S_cc_new for the pressure iterations
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(S_cc_new[lev], S_cc_old[lev], 0, 0, 1, 0);
        }

        // initial (pressure) iters
        for (int i = 1; i <= init_iter; ++i) {
            Print() << "Doing initial pressure iteration #" << i << std::endl;
            InitIter();
        }

        if (plot_int > 0 || plot_deltat > 0) {
            Print() << "\nWriting plotfile 0 after all initialization"
                    << std::endl;
            WritePlotFile(0, t_old, dt, rho0_old, rhoh0_old, p0_old,
                          gamma1bar_old, uold, sold, S_cc_old);
        } else if (small_plot_int > 0 || small_plot_deltat > 0) {
            Print() << "\nWriting small plotfile 0 after all initialization"
                    << std::endl;
            WriteSmallPlotFile(0, t_old, dt, rho0_old, rhoh0_old, p0_old,
                               gamma1bar_old, uold, sold, S_cc_old);
        }

        if (chk_int > 0 || chk_deltat > 0) {
            Print() << "\nWriting checkpoint 0 after all initialization"
                    << std::endl;
            WriteCheckPoint(0);
        }

        if (sum_interval > 0 || sum_per > 0) {
            int index_dummy = 0;
            Print() << "\nWriting diagnosis file after all initialization"
                    << std::endl;
            DiagFile(0, t_old, rho0_old, p0_old, uold, sold, index_dummy);
        }
    }

#ifdef DO_PROBLEM_POST_INIT
    ProblemPostInit();
#endif
}

// fill in multifab and base state data
void Maestro::InitData() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    Print() << "Calling InitData()" << std::endl;

    Print() << "initdata model_File = " << model_file << std::endl;

    // read in model file and fill in s0_init and p0_init for all levels
    s0_init.setVal(0.0);
    p0_init.setVal(0.0);

    for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
        InitBaseState(rho0_old, rhoh0_old, p0_old, lev);
    }

    if (use_exact_base_state || !evolve_base_state) {
        psi.setVal(0.0);
    }

    // calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
    // that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
    InitFromScratch(t_old);

    if (!spherical) {
        // reset tagging array to include buffer zones
        TagArray();
    }

    // set finest_radial_level
    // compute numdisjointchunks, r_start_coord, r_end_coord
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());

    // average down data and fill ghost cells
    AverageDown(sold, 0, Nscal);
    FillPatch(t_old, sold, sold, sold, 0, 0, Nscal, 0, bcs_s);
    AverageDown(uold, 0, AMREX_SPACEDIM);
    FillPatch(t_old, uold, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

    if (fix_base_state) {
        // compute cutoff coordinates
        base_geom.ComputeCutoffCoords(rho0_old.array());
        MakeGravCell(grav_cell_old, rho0_old);
    } else {
        // first compute cutoff coordinates using initial density profile
        base_geom.ComputeCutoffCoords(rho0_old.array());

        if (do_smallscale) {
            // set rho0_old = rhoh0_old = 0.
            rho0_old.setVal(0.0);
            rhoh0_old.setVal(0.0);
        } else {
            // set rho0 to be the average
            Average(sold, rho0_old, Rho);
            base_geom.ComputeCutoffCoords(rho0_old.array());

            // compute gravity
            MakeGravCell(grav_cell_old, rho0_old);

            // compute p0 with HSE
            EnforceHSE(rho0_old, p0_old, grav_cell_old);

            // call eos with r,p as input to recompute T,h
            TfromRhoP(sold, p0_old, true);

            // set rhoh0 to be the average
            Average(sold, rhoh0_old, RhoH);
        }

        // set tempbar to be the average
        Average(sold, tempbar, Temp);
        tempbar_init.copy(tempbar);
    }

    // set p0^{-1} = p0_old
    p0_nm1.copy(p0_old);

    // initialize these since an initial plotfile needs valid data in here
    // gamma1bar_old.setVal(0.);
    // w0.setVal(0.);
}

// During initialization of a simulation, Maestro::InitData() calls
// AmrCore::InitFromScratch(), which calls
// a MakeNewGrids() function that repeatedly calls this function to build
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
                                      const DistributionMapping& dm) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNewLevelFromScratch()",
                   MakeNewLevelFromScratch);

    sold[lev].define(ba, dm, Nscal, ng_s);
    snew[lev].define(ba, dm, Nscal, ng_s);
    uold[lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
    unew[lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
    S_cc_old[lev].define(ba, dm, 1, 0);
    S_cc_new[lev].define(ba, dm, 1, 0);
    gpi[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define(ba, dm, 1, 0);
    w0_cart[lev].define(ba, dm, AMREX_SPACEDIM, 2);
    rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);

    pi[lev].define(convert(ba, nodal_flag), dm, 1, 0);  // nodal
    intra[lev].define(ba, dm, Nscal, 0);                // for sdc

    sold[lev].setVal(0.);
    snew[lev].setVal(0.);
    uold[lev].setVal(0.);
    unew[lev].setVal(0.);
    S_cc_old[lev].setVal(0.);
    S_cc_new[lev].setVal(0.);
    gpi[lev].setVal(0.);
    dSdt[lev].setVal(0.);
    w0_cart[lev].setVal(0.);
    rhcc_for_nodalproj[lev].setVal(0.);
    pi[lev].setVal(0.);
    intra[lev].setVal(0.);

    if (spherical) {
        normal[lev].define(ba, dm, 3, 1);
        cell_cc_to_r[lev].define(ba, dm, 1, 0);
    }

    if (!spherical) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Array4<Real> scal_arr = sold[lev].array(mfi);
            const Array4<Real> vel_arr = uold[lev].array(mfi);

            InitLevelData(lev, t_old, mfi, scal_arr, vel_arr);
        }
    } else {
#if (AMREX_SPACEDIM == 3)
        const auto dx_fine_vec = geom[max_level].CellSizeArray();
        const auto dx_lev = geom[lev].CellSizeArray();
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            InitBaseStateMapSphr(lev, mfi, dx_fine_vec, dx_lev);
        }

        InitLevelDataSphr(lev, t_old, sold[lev], uold[lev]);
#endif
    }

    if (lev > 0 && reflux_type == 2) {
        flux_reg_s[lev] = std::make_unique<FluxRegister>(
            ba, dm, refRatio(lev - 1), lev, Nscal);
    }
}

void Maestro::InitProj() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitProj()", InitProj);

    Vector<MultiFab> rho_omegadot(finest_level + 1);
    Vector<MultiFab> thermal(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);
    Vector<MultiFab> rhohalf(finest_level + 1);
    Vector<MultiFab> Tcoeff(finest_level + 1);
    Vector<MultiFab> hcoeff(finest_level + 1);
    Vector<MultiFab> Xkcoeff(finest_level + 1);
    Vector<MultiFab> pcoeff(finest_level + 1);
    Vector<MultiFab> delta_gamma1(finest_level + 1);
    Vector<MultiFab> delta_gamma1_term(finest_level + 1);

    BaseState<Real> Sbar(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> delta_gamma1_termbar(base_geom.max_radial_level + 1,
                                         base_geom.nr_fine);

    for (int lev = 0; lev <= finest_level; ++lev) {
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        thermal[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        rhohalf[lev].define(grids[lev], dmap[lev], 1, 1);
        Tcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        hcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        Xkcoeff[lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev], 1, 1);

        // we don't have a legit timestep yet, so we set rho_omegadot,
        // rho_Hnuc, and rho_Hext to 0
        rho_omegadot[lev].setVal(0.);
        rho_Hnuc[lev].setVal(0.);
        rho_Hext[lev].setVal(0.);
        delta_gamma1[lev].setVal(0.);
        delta_gamma1_term[lev].setVal(0.);

        // initial projection does not use density weighting
        rhohalf[lev].setVal(1.);
    }

    // compute thermal diffusion
    if (use_thermal_diffusion) {
        MakeThermalCoeffs(sold, Tcoeff, hcoeff, Xkcoeff, pcoeff);

        MakeExplicitThermal(thermal, sold, Tcoeff, hcoeff, Xkcoeff, pcoeff,
                            p0_old, temp_diffusion_formulation);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            thermal[lev].setVal(0.);
        }
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_old, delta_gamma1_term, delta_gamma1, sold, uold,
              rho_omegadot, rho_Hnuc, rho_Hext, thermal, p0_old, gamma1bar_old,
              delta_gamma1_termbar);

    // NOTE: not sure if valid for use_exact_base_state
    if (evolve_base_state && (!use_exact_base_state && !average_base_state)) {
        // average S into Sbar
        Average(S_cc_old, Sbar, 0);
    } else {
        Sbar.setVal(0.);
    }

    // make the nodal rhs for projection beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforNodalProj(rhcc_for_nodalproj, S_cc_old, Sbar, beta0_old,
                         delta_gamma1_term);

    // perform a nodal projection
#ifndef SDC
    NodalProj(initial_projection_comp, rhcc_for_nodalproj);
#else
    NodalProj(initial_projection_comp, rhcc_for_nodalproj, false);
#endif
}

void Maestro::DivuIter(int istep_divu_iter) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DivuIter()", DivuIter);

    Vector<MultiFab> stemp(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);
    Vector<MultiFab> rho_omegadot(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
    Vector<MultiFab> thermal(finest_level + 1);
    Vector<MultiFab> rhohalf(finest_level + 1);
    Vector<MultiFab> Tcoeff(finest_level + 1);
    Vector<MultiFab> hcoeff(finest_level + 1);
    Vector<MultiFab> Xkcoeff(finest_level + 1);
    Vector<MultiFab> pcoeff(finest_level + 1);
    Vector<MultiFab> delta_gamma1(finest_level + 1);
    Vector<MultiFab> delta_gamma1_term(finest_level + 1);

    BaseState<Real> Sbar(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> w0_force(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0_minus_peosbar(base_geom.max_radial_level + 1,
                                     base_geom.nr_fine);
    BaseState<Real> delta_gamma1_termbar(base_geom.max_radial_level + 1,
                                         base_geom.nr_fine);

    Sbar.setVal(0.);
    etarho_ec.setVal(0.0);
    w0_force.setVal(0.0);
    psi.setVal(0.0);
    etarho_cc.setVal(0.0);
    p0_minus_peosbar.setVal(0.);
    delta_gamma1_termbar.setVal(0.);

    for (int lev = 0; lev <= finest_level; ++lev) {
        stemp[lev].define(grids[lev], dmap[lev], Nscal, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
        thermal[lev].define(grids[lev], dmap[lev], 1, 0);
        rhohalf[lev].define(grids[lev], dmap[lev], 1, 1);
        Tcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        hcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        Xkcoeff[lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev], 1, 1);

        // divu_iters do not use density weighting
        rhohalf[lev].setVal(1.);
    }

    React(sold, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_old, 0.5 * dt,
          t_old);

    // WriteMF(sold,"a_sold_2levs");
    // Abort();

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(sold, Tcoeff, hcoeff, Xkcoeff, pcoeff);

        MakeExplicitThermal(thermal, sold, Tcoeff, hcoeff, Xkcoeff, pcoeff,
                            p0_old, temp_diffusion_formulation);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            thermal[lev].setVal(0.);
        }
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_old, delta_gamma1_term, delta_gamma1, sold, uold,
              rho_omegadot, rho_Hnuc, rho_Hext, thermal, p0_old, gamma1bar_old,
              delta_gamma1_termbar);

    // NOTE: not sure if valid for use_exact_base_state
    if (evolve_base_state) {
        if ((use_exact_base_state || average_base_state) &&
            use_delta_gamma1_term) {
            Sbar += delta_gamma1_termbar;
        } else {
            Average(S_cc_old, Sbar, 0);

            // compute Sbar = Sbar + delta_gamma1_termbar
            if (use_delta_gamma1_term) {
                Sbar += delta_gamma1_termbar;
            }

            const auto is_predictor = true;
            Makew0(w0, w0_force, Sbar, rho0_old, rho0_old, p0_old, p0_old,
                   gamma1bar_old, gamma1bar_old, p0_minus_peosbar, dt, dt,
                   is_predictor);

            // put w0 on Cartesian cell-centers
            Put1dArrayOnCart(w0, w0_cart, true, true, bcs_u, 0, 1);
        }
    }

    // make the nodal rhs for projection beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforNodalProj(rhcc_for_nodalproj, S_cc_old, Sbar, beta0_old,
                         delta_gamma1_term);

    // perform a nodal projection
    NodalProj(divu_iters_comp, rhcc_for_nodalproj, istep_divu_iter);

    Real dt_hold = dt;

    // compute new time step
    EstDt();

    if (maestro_verbose > 0) {
        Print() << "Call to estdt at end of istep_divu_iter = "
                << istep_divu_iter << " gives dt = " << dt << std::endl;
    }

    dt *= init_shrink;
    if (maestro_verbose > 0) {
        Print() << "Multiplying dt by init_shrink; dt = " << dt << std::endl;
    }

    if (dt > dt_hold) {
        if (maestro_verbose > 0) {
            Print() << "Ignoring this new dt since it's larger than the "
                       "previous dt = "
                    << dt_hold << std::endl;
        }
        dt = amrex::min(dt_hold, dt);
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }
}

// SDC
void Maestro::DivuIterSDC(int istep_divu_iter) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DivuIterSDC()", DivuIterSDC);

    Vector<MultiFab> stemp(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);
    Vector<MultiFab> rho_omegadot(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
    Vector<MultiFab> thermal(finest_level + 1);
    Vector<MultiFab> rhohalf(finest_level + 1);
    Vector<MultiFab> Tcoeff(finest_level + 1);
    Vector<MultiFab> hcoeff(finest_level + 1);
    Vector<MultiFab> Xkcoeff(finest_level + 1);
    Vector<MultiFab> pcoeff(finest_level + 1);
    Vector<MultiFab> delta_gamma1(finest_level + 1);
    Vector<MultiFab> delta_gamma1_term(finest_level + 1);
    Vector<MultiFab> sdc_source(finest_level + 1);

    BaseState<Real> Sbar(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> w0_force(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0_minus_pthermbar(base_geom.max_radial_level + 1,
                                       base_geom.nr_fine);
    BaseState<Real> delta_gamma1_termbar(base_geom.max_radial_level + 1,
                                         base_geom.nr_fine);

    Sbar.setVal(0.);
    etarho_ec.setVal(0.0);
    w0_force.setVal(0.0);
    psi.setVal(0.0);
    etarho_cc.setVal(0.0);
    p0_minus_pthermbar.setVal(0.);
    delta_gamma1_termbar.setVal(0.);

    for (int lev = 0; lev <= finest_level; ++lev) {
        stemp[lev].define(grids[lev], dmap[lev], Nscal, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
        thermal[lev].define(grids[lev], dmap[lev], 1, 0);
        rhohalf[lev].define(grids[lev], dmap[lev], 1, 1);
        Tcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        hcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        Xkcoeff[lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1[lev].define(grids[lev], dmap[lev], 1, 1);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev], 1, 1);
        sdc_source[lev].define(grids[lev], dmap[lev], Nscal, 0);

        // divu_iters do not use density weighting
        rhohalf[lev].setVal(1.);
        sdc_source[lev].setVal(0.);
    }

    ReactSDC(sold, stemp, rho_Hext, p0_old, 0.5 * dt, t_old, sdc_source);

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(sold, Tcoeff, hcoeff, Xkcoeff, pcoeff);

        MakeExplicitThermal(thermal, sold, Tcoeff, hcoeff, Xkcoeff, pcoeff,
                            p0_old, temp_diffusion_formulation);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            thermal[lev].setVal(0.);
        }
    }

    MakeReactionRates(rho_omegadot, rho_Hnuc, sold);

    // compute S at cell-centers
    Make_S_cc(S_cc_old, delta_gamma1_term, delta_gamma1, sold, uold,
              rho_omegadot, rho_Hnuc, rho_Hext, thermal, p0_old, gamma1bar_old,
              delta_gamma1_termbar);

    // NOTE: not sure if valid for use_exact_base_state
    if (evolve_base_state) {
        if ((use_exact_base_state || average_base_state) &&
            use_delta_gamma1_term) {
            Sbar += delta_gamma1_termbar;
        } else {
            Average(S_cc_old, Sbar, 0);

            // compute Sbar = Sbar + delta_gamma1_termbar
            if (use_delta_gamma1_term) {
                Sbar += delta_gamma1_termbar;
            }

            const auto is_predictor = true;
            Makew0(w0, w0_force, Sbar, rho0_old, rho0_old, p0_old, p0_old,
                   gamma1bar_old, gamma1bar_old, p0_minus_pthermbar, dt, dt,
                   is_predictor);
        }
    }

    // make the nodal rhs for projection beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforNodalProj(rhcc_for_nodalproj, S_cc_old, Sbar, beta0_old,
                         delta_gamma1_term);

    // perform a nodal projection
    NodalProj(divu_iters_comp, rhcc_for_nodalproj, istep_divu_iter, false);

    Real dt_hold = dt;

    // compute new time step
    EstDt();

    if (maestro_verbose > 0) {
        Print() << "Call to estdt at end of istep_divu_iter = "
                << istep_divu_iter << " gives dt = " << dt << std::endl;
    }

    dt *= init_shrink;
    if (maestro_verbose > 0) {
        Print() << "Multiplying dt by init_shrink; dt = " << dt << std::endl;
    }

    if (dt > dt_hold) {
        if (maestro_verbose > 0) {
            Print() << "Ignoring this new dt since it's larger than the "
                       "previous dt = "
                    << dt_hold << std::endl;
        }
        dt = amrex::min(dt_hold, dt);
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }
}

void Maestro::InitIter() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitIter()", InitIter);

    // wallclock time
    Real start_total = ParallelDescriptor::second();

    // advance the solution by dt
#ifndef SDC
    if (use_exact_base_state || average_base_state) {
        AdvanceTimeStepAverage(true);
    } else {
        AdvanceTimeStep(true);
    }
#else
    AdvanceTimeStepSDC(true);
#endif

    // wallclock time
    Real end_total = ParallelDescriptor::second() - start_total;
    ParallelDescriptor::ReduceRealMax(end_total,
                                      ParallelDescriptor::IOProcessorNumber());

    Print() << "Time to advance time step: " << end_total << '\n';

    // copy pi from snew to sold
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Copy(sold[lev], snew[lev], Pi, Pi, 1, ng_s);
    }
}
