
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void Maestro::AdvanceTimeStep(bool is_initIter) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStep()", AdvanceTimeStep);

    // timers
    Real advect_time = 0., advect_time_start;
    Real macproj_time = 0., macproj_time_start;
    Real ndproj_time = 0., ndproj_time_start;
    Real thermal_time = 0., thermal_time_start;
    Real react_time = 0., react_time_start;
    Real misc_time = 0., misc_time_start;
    Real base_time = 0., base_time_start;

    // HACK
    base_time_start = ParallelDescriptor::second();

    base_time += ParallelDescriptor::second() - base_time_start;
    ParallelDescriptor::ReduceRealMax(base_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&base_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    // features to be added later:
    // -ppm
    // -dpdt_factor

    // cell-centered MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab> rhohalf(finest_level + 1);
    Vector<MultiFab> macrhs(finest_level + 1);
    Vector<MultiFab> macphi(finest_level + 1);
    Vector<MultiFab> S_cc_nph(finest_level + 1);
    Vector<MultiFab> rho_omegadot(finest_level + 1);
#if NAUX_NET > 0
    Vector<MultiFab> rho_auxdot(finest_level + 1);
#endif
    Vector<MultiFab> thermal1(finest_level + 1);
    Vector<MultiFab> thermal2(finest_level + 1);
    Vector<MultiFab> rho_Hnuc(finest_level + 1);
    Vector<MultiFab> rho_Hext(finest_level + 1);
    Vector<MultiFab> s1(finest_level + 1);
    Vector<MultiFab> s2(finest_level + 1);
    Vector<MultiFab> s2star(finest_level + 1);
    Vector<MultiFab> delta_gamma1_term(finest_level + 1);
    Vector<MultiFab> delta_gamma1(finest_level + 1);
    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> delta_p_term(finest_level + 1);
    Vector<MultiFab> Tcoeff(finest_level + 1);
    Vector<MultiFab> hcoeff1(finest_level + 1);
    Vector<MultiFab> Xkcoeff1(finest_level + 1);
    Vector<MultiFab> pcoeff1(finest_level + 1);
    Vector<MultiFab> hcoeff2(finest_level + 1);
    Vector<MultiFab> Xkcoeff2(finest_level + 1);
    Vector<MultiFab> pcoeff2(finest_level + 1);
    Vector<MultiFab> scal_force(finest_level + 1);
    Vector<MultiFab> delta_chi(finest_level + 1);
    Vector<MultiFab> sponge(finest_level + 1);

    // face-centered in the dm-direction (planar only)
    Vector<MultiFab> etarhoflux(finest_level + 1);

    // face-centered
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > sedge(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > sflux(finest_level + 1);

    ////////////////////////
    // needed for spherical routines only

    // cell-centered
    Vector<MultiFab> w0_force_cart(finest_level + 1);

    // face-centered
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);

    // end spherical-only MultiFabs
    ////////////////////////

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    BaseState<Real> grav_cell_nph(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine);
    BaseState<Real> rho0_nph(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0_nph(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0_minus_peosbar(base_geom.max_radial_level + 1,
                                     base_geom.nr_fine);
    BaseState<Real> peosbar(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> w0_force(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> Sbar(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> beta0_nph(base_geom.max_radial_level + 1,
                              base_geom.nr_fine);
    BaseState<Real> gamma1bar_temp1(base_geom.max_radial_level + 1,
                                    base_geom.nr_fine);
    BaseState<Real> gamma1bar_temp2(base_geom.max_radial_level + 1,
                                    base_geom.nr_fine);
    BaseState<Real> delta_gamma1_termbar(base_geom.max_radial_level + 1,
                                         base_geom.nr_fine);

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    BaseState<Real> w0_old(base_geom.max_radial_level + 1,
                           base_geom.nr_fine + 1);
    BaseState<Real> rho0_predicted_edge(base_geom.max_radial_level + 1,
                                        base_geom.nr_fine + 1);

    bool is_predictor;

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl
            << std::endl;

    if (maestro_verbose > 0) {
        Print() << "Cell Count:" << std::endl;
        for (int lev = 0; lev <= finest_level; ++lev) {
            Print() << "Level " << lev << ", " << CountCells(lev) << " cells"
                    << std::endl;
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // cell-centered MultiFabs
        rhohalf[lev].define(grids[lev], dmap[lev], 1, 1);
        macrhs[lev].define(grids[lev], dmap[lev], 1, 0);
        macphi[lev].define(grids[lev], dmap[lev], 1, 1);
        S_cc_nph[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
#if NAUX_NET > 0
        rho_auxdot[lev].define(grids[lev], dmap[lev], NumAux, 0);
#endif
        thermal1[lev].define(grids[lev], dmap[lev], 1, 0);
        thermal2[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hnuc[lev].define(grids[lev], dmap[lev], 1, 0);
        rho_Hext[lev].define(grids[lev], dmap[lev], 1, 0);
        s1[lev].define(grids[lev], dmap[lev], Nscal, ng_s);
        s1[lev].setVal(0.);
        s2[lev].define(grids[lev], dmap[lev], Nscal, ng_s);
        s2star[lev].define(grids[lev], dmap[lev], Nscal, ng_s);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev], 1, 0);
        delta_gamma1[lev].define(grids[lev], dmap[lev], 1, 0);
        delta_p_term[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        Tcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
        hcoeff1[lev].define(grids[lev], dmap[lev], 1, 1);
        Xkcoeff1[lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff1[lev].define(grids[lev], dmap[lev], 1, 1);
        hcoeff2[lev].define(grids[lev], dmap[lev], 1, 1);
        Xkcoeff2[lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff2[lev].define(grids[lev], dmap[lev], 1, 1);
        if (ppm_trace_forces == 0) {
            scal_force[lev].define(grids[lev], dmap[lev], Nscal, 1);
        } else {
            // we need more ghostcells if we are tracing the forces
            scal_force[lev].define(grids[lev], dmap[lev], Nscal, ng_s);
        }
        delta_chi[lev].define(grids[lev], dmap[lev], 1, 0);
        sponge[lev].define(grids[lev], dmap[lev], 1, 0);

        // face-centered in the dm-direction (planar only)
        AMREX_D_TERM(etarhoflux[lev].define(convert(grids[lev], nodal_flag_x),
                                            dmap[lev], 1, 1);
                     , etarhoflux[lev].define(convert(grids[lev], nodal_flag_y),
                                              dmap[lev], 1, 1);
                     , etarhoflux[lev].define(convert(grids[lev], nodal_flag_z),
                                              dmap[lev], 1, 1););

        // face-centered arrays of MultiFabs
        AMREX_D_TERM(umac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                         dmap[lev], 1, 1);
                     , umac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                           dmap[lev], 1, 1);
                     , umac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                           dmap[lev], 1, 1););
        AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], Nscal, 0);
                     , sedge[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], Nscal, 0);
                     , sedge[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], Nscal, 0););
        AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], Nscal, 0);
                     , sflux[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], Nscal, 0);
                     , sflux[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], Nscal, 0););

        // initialize umac
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            umac[lev][d].setVal(0.);
            sedge[lev][d].setVal(0.);
            sflux[lev][d].setVal(0.);
        }

        w0_force_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

#if (AMREX_SPACEDIM == 3)
    for (int lev = 0; lev <= finest_level; ++lev) {
        w0mac[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                             1);
        w0mac[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                             1);
        w0mac[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                             1);
    }
#endif

    if (!evolve_base_state) {
        w0.setVal(0.0);
    }

    // make the sponge for all levels
    if (do_sponge) {
        SpongeInit(rho0_old);
        MakeSponge(sponge);
    }

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 1 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    react_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 1 : react state >>>" << std::endl;
    }

    React(sold, s1, rho_Hext, rho_omegadot, rho_Hnuc,
#if NAUX_NET > 0
          rho_auxdot,
#endif
          p0_old, 0.5 * dt, t_old);

    react_time += ParallelDescriptor::second() - react_time_start;
    ParallelDescriptor::ReduceRealMax(react_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&react_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- define average expansion at time n+1/2
    //////////////////////////////////////////////////////////////////////////////

    advect_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2 : make w0 >>>" << std::endl;
    }

    if (t_old == 0.) {
        // this is either a pressure iteration or the first time step
        // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev], 0.5, S_cc_old[lev], 0, 0.5,
                              S_cc_new[lev], 0, 0, 1, 0);
        }
    } else {
        // set S_cc_nph = S_cc_old + (dt/2) * dSdt
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev], 1.0, S_cc_old[lev], 0, 0.5 * dt,
                              dSdt[lev], 0, 0, 1, 0);
        }
    }
    // no ghost cells for S_cc_nph
    AverageDown(S_cc_nph, 0, 1);

    // compute p0_minus_peosbar = p0_old - peosbar_old (for making w0) and
    // compute delta_p_term = peos_old - p0_old (for RHS of projections)
    if (dpdt_factor > 0.0) {
        // peos_old (delta_p_term) now holds the thermodynamic p computed from sold(rho,h,X)
        PfromRhoH(sold, sold, delta_p_term);

        // compute peosbar = Avg(peos_old)
        Average(delta_p_term, peosbar, 0);

        // compute p0_minus_peosbar = p0_old - peosbar
        p0_minus_peosbar.copy(p0_old - peosbar);

        // put p0_old on cart
        Put1dArrayOnCart(p0_old, p0_cart, false, false, bcs_f, 0);

        // compute delta_p_term = peos_old - p0_old
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(delta_p_term[lev], p0_cart[lev], 0, 0, 1, 0);
        }
    } else {
        // these should have no effect if dpdt_factor <= 0
        p0_minus_peosbar.setVal(0.0);
        for (int lev = 0; lev <= finest_level; ++lev) {
            delta_p_term[lev].setVal(0.);
        }
    }

#if (AMREX_SPACEDIM == 3)
    // initialize MultiFabs and Vectors to ZERO
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            w0mac[lev][d].setVal(0.);
        }
    }
#endif
    for (int lev = 0; lev <= finest_level; ++lev) {
        w0_force_cart[lev].setVal(0.);
    }

    if (evolve_base_state) {
        // compute Sbar = average(S_cc_nph)
        Average(S_cc_nph, Sbar, 0);

        // save old-time value
        w0_old.copy(w0);

        base_time_start = ParallelDescriptor::second();

        ComputeCutoffCoords(rho0_old);
        base_geom.ComputeCutoffCoords(rho0_old.array());

        // compute w0, w0_force
        is_predictor = true;
        Makew0(w0_old, w0_force, Sbar, rho0_old, rho0_old, p0_old, p0_old,
               gamma1bar_old, gamma1bar_old, p0_minus_peosbar, dt, dtold,
               is_predictor);

        Put1dArrayOnCart(w0, w0_cart, true, true, bcs_u, 0, 1);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        // put w0 on Cartesian edges
#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            MakeW0mac(w0mac);
        }
#endif
        // put w0_force on Cartesian cells
        Put1dArrayOnCart(w0_force, w0_force_cart, false, true, bcs_f, 0);

    } else {
        // these should have no effect if evolve_base_state = false
        Sbar.setVal(0.);
        w0_force.setVal(0.);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 3 -- construct the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 3 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    is_predictor = true;
    AdvancePremac(umac, w0mac, w0_force_cart);

    for (int lev = 0; lev <= finest_level; ++lev) {
        delta_chi[lev].setVal(0.);
        macphi[lev].setVal(0.);
        delta_gamma1_term[lev].setVal(0.);
    }

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforMacProj(macrhs, rho0_old, S_cc_nph, Sbar, beta0_old,
                       delta_gamma1_term, gamma1bar_old, p0_old, delta_p_term,
                       delta_chi, is_predictor);

    advect_time += ParallelDescriptor::second() - advect_time_start;
    ParallelDescriptor::ReduceRealMax(advect_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&advect_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    macproj_time_start = ParallelDescriptor::second();

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac, macphi, macrhs, beta0_old, is_predictor);

    macproj_time += ParallelDescriptor::second() - macproj_time_start;
    ParallelDescriptor::ReduceRealMax(macproj_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&macproj_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4 -- advect the base state and full state through dt
    //////////////////////////////////////////////////////////////////////////////

    advect_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 4 : advect base >>>" << std::endl;
    }

    // advect the base state density
    if (evolve_base_state) {
        base_time_start = ParallelDescriptor::second();

        AdvectBaseDens(rho0_predicted_edge);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        ComputeCutoffCoords(rho0_new);
        base_geom.ComputeCutoffCoords(rho0_new.array());
    } else {
        rho0_new.copy(rho0_old);
    }

    // thermal is the forcing for rhoh or temperature
    if (use_thermal_diffusion) {
        MakeThermalCoeffs(s1, Tcoeff, hcoeff1, Xkcoeff1, pcoeff1);

        MakeExplicitThermal(thermal1, s1, Tcoeff, hcoeff1, Xkcoeff1, pcoeff1,
                            p0_old, temp_diffusion_formulation);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            thermal1[lev].setVal(0.);
        }
    }

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev = 0; lev <= finest_level; ++lev) {
        s2[lev].setVal(0.);
        MultiFab::Copy(s2[lev], s1[lev], Temp, Temp, 1, ng_s);
    }

    if (maestro_verbose >= 1) {
        Print() << "            :  density_advance >>>" << std::endl;
        Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // set etarhoflux to zero
    for (int lev = 0; lev <= finest_level; ++lev) {
        etarhoflux[lev].setVal(0.);
    }

    // need full UMAC velocities for DensityAdvance
    Addw0(umac, w0mac, 1.);

    // advect rhoX, rho, and tracers
    DensityAdvance(1, s1, s2, sedge, sflux, scal_force, etarhoflux, umac, w0mac,
                   rho0_predicted_edge);

    // subtract w0mac from umac
    Addw0(umac, w0mac, -1.);

    if (evolve_base_state) {
        if (use_etarho) {
            // compute the new etarho
            if (!spherical) {
                MakeEtarho(etarhoflux);
            } else {
                MakeEtarhoSphr(s1, s2, umac, w0mac);
            }

            // correct the base state density by "averaging"
            Average(s2, rho0_new, Rho);
            ComputeCutoffCoords(rho0_new);
            base_geom.ComputeCutoffCoords(rho0_new.array());
        }

        // update grav_cell_new
        base_time_start = ParallelDescriptor::second();

        MakeGravCell(grav_cell_new, rho0_new);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        // base state pressure update
        // set new p0 through HSE
        p0_new.copy(p0_old);

        base_time_start = ParallelDescriptor::second();

        EnforceHSE(rho0_new, p0_new, grav_cell_new);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        // make psi
        if (!spherical) {
            MakePsiPlanar();
        } else {
            // compute p0_nph
            p0_nph.copy(0.5 * (p0_old + p0_new));

            // compute gamma1bar^{(1)} and store it in gamma1bar_temp1
            MakeGamma1bar(s1, gamma1bar_temp1, p0_old);

            // compute gamma1bar^{(2),*} and store it in gamma1bar_temp2
            MakeGamma1bar(s2, gamma1bar_temp2, p0_new);

            // compute gamma1bar^{nph,*} and store it in gamma1bar_temp2
            gamma1bar_temp2.copy(0.5 * (gamma1bar_temp1 + gamma1bar_temp2));

            // make time-centered psi
            MakePsiSphr(gamma1bar_temp2, p0_nph, Sbar);
        }

        // base state enthalpy update
        // compute rhoh0_old by "averaging"
        Average(s1, rhoh0_old, RhoH);

        base_time_start = ParallelDescriptor::second();

        AdvectBaseEnthalpy(rho0_predicted_edge);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());
    } else {
        rhoh0_new.copy(rhoh0_old);
        grav_cell_new.copy(grav_cell_old);
        p0_new.copy(p0_old);
    }

    if (maestro_verbose >= 1) {
        Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    // need full UMAC velocities for EnthalpyAdvance
    Addw0(umac, w0mac, 1.);

    EnthalpyAdvance(1, s1, s2, sedge, sflux, scal_force, umac, w0mac, thermal1);

    // subtract w0mac from umac
    Addw0(umac, w0mac, -1.);

    advect_time += ParallelDescriptor::second() - advect_time_start;
    ParallelDescriptor::ReduceRealMax(advect_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&advect_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    thermal_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 4a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
        ThermalConduct(s1, s2, hcoeff1, Xkcoeff1, pcoeff1, hcoeff1, Xkcoeff1,
                       pcoeff1);
    }

    thermal_time += ParallelDescriptor::second() - thermal_time_start;
    ParallelDescriptor::ReduceRealMax(thermal_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&thermal_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Copy(s2[lev], s1[lev], Temp, Temp, 1, ng_s);
        MultiFab::Copy(s2[lev], s1[lev], Pi, Pi, 1, ng_s);
    }

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s2, p0_new, false);
    } else {
        TfromRhoH(s2, p0_new);
    }

    if (use_thermal_diffusion) {
        // make a copy of s2star since these are needed to compute
        // coefficients in the call to thermal_conduct_full_alg
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(s2star[lev], s2[lev], 0, 0, Nscal, ng_s);
        }
    }

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    react_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 5 : react state >>>" << std::endl;
    }

    React(s2, snew, rho_Hext, rho_omegadot, rho_Hnuc,
#if NAUX_NET > 0
          rho_auxdot,
#endif
          p0_new, 0.5 * dt, t_old + 0.5 * dt);

    react_time += ParallelDescriptor::second() - react_time_start;
    ParallelDescriptor::ReduceRealMax(react_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&react_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    if (evolve_base_state) {
        // compute beta0 and gamma1bar
        MakeGamma1bar(snew, gamma1bar_new, p0_new);

        base_time_start = ParallelDescriptor::second();

        MakeBeta0(beta0_new, rho0_new, p0_new, gamma1bar_new, grav_cell_new);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());
    } else {
        // Just pass beta0 and gamma1bar through if not evolving base state
        beta0_new.copy(beta0_old);
        gamma1bar_new.copy(gamma1bar_old);
    }

    beta0_nph.copy(0.5 * (beta0_old + beta0_new));

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 6 -- define a new average expansion rate at n+1/2
    //////////////////////////////////////////////////////////////////////////////

    advect_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 6 : make new S and new w0 >>>" << std::endl;
    }

    if (evolve_base_state) {
        // reset cutoff coordinates to old time value
        ComputeCutoffCoords(rho0_old);
        base_geom.ComputeCutoffCoords(rho0_old.array());
    }

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(snew, Tcoeff, hcoeff2, Xkcoeff2, pcoeff2);

        MakeExplicitThermal(thermal2, snew, Tcoeff, hcoeff2, Xkcoeff2, pcoeff2,
                            p0_new, temp_diffusion_formulation);
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            thermal2[lev].setVal(0.);
        }
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_new, delta_gamma1_term, delta_gamma1, snew, uold,
              rho_omegadot, rho_Hnuc, rho_Hext, thermal2, p0_old, gamma1bar_new,
              delta_gamma1_termbar);

    // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::LinComb(S_cc_nph[lev], 0.5, S_cc_old[lev], 0, 0.5,
                          S_cc_new[lev], 0, 0, 1, 0);
    }
    AverageDown(S_cc_nph, 0, 1);

    // compute p0_minus_peosbar = p0_new - peosbar_new (for making w0) and
    // set delta_p_term = peos_new - p0_new (for RHS of projection)
    if (dpdt_factor > 0.) {
        // peos now holds "peos_new", the thermodynamic p computed from snew(rho,h,X)
        PfromRhoH(snew, snew, delta_p_term);

        // compute peosbar = Avg(peos_new)
        Average(delta_p_term, peosbar, 0);

        // compute p0_minus_peosbar = p0_new - peosbar
        p0_minus_peosbar.copy(p0_new - peosbar);

        // put p0_new on cart
        Put1dArrayOnCart(p0_new, p0_cart, false, false, bcs_f, 0);

        // set delta_p_term = peos_new - p0_new
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(delta_p_term[lev], p0_cart[lev], 0, 0, 1, 0);
        }
    } else {
        // these should have no effect if dpdt_factor <= 0
        p0_minus_peosbar.setVal(0.);
        for (int lev = 0; lev <= finest_level; ++lev) {
            delta_p_term[lev].setVal(0.);
        }
    }

    if (evolve_base_state) {
        // compute Sbar = average(S_cc_nph)
        Average(S_cc_nph, Sbar, 0);

        // compute Sbar = Sbar + delta_gamma1_termbar
        if (use_delta_gamma1_term) {
            Sbar += delta_gamma1_termbar;
        }

        base_time_start = ParallelDescriptor::second();

        ComputeCutoffCoords(rho0_old);
        base_geom.ComputeCutoffCoords(rho0_old.array());

        // compute w0, w0_force
        is_predictor = false;
        Makew0(w0_old, w0_force, Sbar, rho0_old, rho0_new, p0_old, p0_new,
               gamma1bar_old, gamma1bar_new, p0_minus_peosbar, dt, dtold,
               is_predictor);

        Put1dArrayOnCart(w0, w0_cart, true, true, bcs_u, 0, 1);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            // put w0 on Cartesian edges
            MakeW0mac(w0mac);
        }
#endif
        // // put w0_force on Cartesian cells
        Put1dArrayOnCart(w0_force, w0_force_cart, false, true, bcs_f, 0);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 7 -- redo the construction of the advective velocity using the current w0
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 7 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    is_predictor = false;
    AdvancePremac(umac, w0mac, w0_force_cart);

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforMacProj(macrhs, rho0_new, S_cc_nph, Sbar, beta0_nph,
                       delta_gamma1_term, gamma1bar_new, p0_new, delta_p_term,
                       delta_chi, is_predictor);

    advect_time += ParallelDescriptor::second() - advect_time_start;
    ParallelDescriptor::ReduceRealMax(advect_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&advect_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    macproj_time_start = ParallelDescriptor::second();

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac, macphi, macrhs, beta0_nph, is_predictor);

    macproj_time += ParallelDescriptor::second() - macproj_time_start;
    ParallelDescriptor::ReduceRealMax(macproj_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&macproj_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8 -- advect the base state and full state through dt
    //////////////////////////////////////////////////////////////////////////////

    advect_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 8 : advect base >>>" << std::endl;
    }

    // advect the base state density
    if (evolve_base_state) {
        base_time_start = ParallelDescriptor::second();

        AdvectBaseDens(rho0_predicted_edge);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        ComputeCutoffCoords(rho0_new);
        base_geom.ComputeCutoffCoords(rho0_new.array());
    }

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Copy(s2[lev], s1[lev], Temp, Temp, 1, ng_s);
    }

    if (maestro_verbose >= 1) {
        Print() << "            :  density_advance >>>" << std::endl;
        Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // set etarhoflux to zero
    for (int lev = 0; lev <= finest_level; ++lev) {
        etarhoflux[lev].setVal(0.);
    }

    // need full UMAC velocities for DensityAdvance
    Addw0(umac, w0mac, 1.);

    // advect rhoX, rho, and tracers
    DensityAdvance(2, s1, s2, sedge, sflux, scal_force, etarhoflux, umac, w0mac,
                   rho0_predicted_edge);

    // subtract w0mac from umac
    Addw0(umac, w0mac, -1.);

    if (evolve_base_state) {
        if (use_etarho) {
            // compute the new etarho
            if (!spherical) {
                MakeEtarho(etarhoflux);
            } else {
                MakeEtarhoSphr(s1, s2, umac, w0mac);
            }

            // correct the base state density by "averaging"
            // call average(mla,s2,rho0_new,dx,rho_comp)
            Average(s2, rho0_new, Rho);
            ComputeCutoffCoords(rho0_new);
            base_geom.ComputeCutoffCoords(rho0_new.array());
        }

        // update grav_cell_new, rho0_nph, grav_cell_nph
        base_time_start = ParallelDescriptor::second();

        MakeGravCell(grav_cell_new, rho0_new);

        rho0_nph.copy(0.5 * (rho0_old + rho0_new));
        MakeGravCell(grav_cell_nph, rho0_nph);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        // base state pressure update
        // set new p0 through HSE
        p0_new.copy(p0_old);

        base_time_start = ParallelDescriptor::second();

        EnforceHSE(rho0_new, p0_new, grav_cell_new);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        p0_nph.copy(0.5 * (p0_old + p0_new));

        // make psi
        if (!spherical) {
            MakePsiPlanar();
        } else {
            // compute gamma1bar^{(2)} and store it in gamma1bar_temp2
            MakeGamma1bar(s2, gamma1bar_temp2, p0_new);

            base_time_start = ParallelDescriptor::second();

            // compute gamma1bar^{nph} and store it in gamma1bar_temp2
            gamma1bar_temp2.copy(0.5 * (gamma1bar_temp1 + gamma1bar_temp2));

            MakePsiSphr(gamma1bar_temp2, p0_nph, Sbar);

            base_time += ParallelDescriptor::second() - base_time_start;
            ParallelDescriptor::ReduceRealMax(
                base_time, ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&base_time, 1,
                                      ParallelDescriptor::IOProcessorNumber());
        }

        base_time_start = ParallelDescriptor::second();

        // base state enthalpy update
        AdvectBaseEnthalpy(rho0_predicted_edge);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());
    } else {
        rho0_nph.copy(rho0_old);
        grav_cell_nph.copy(grav_cell_old);
    }

    if (maestro_verbose >= 1) {
        Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    // need full UMAC velocities for EnthalpyAdvance
    Addw0(umac, w0mac, 1.);

    EnthalpyAdvance(2, s1, s2, sedge, sflux, scal_force, umac, w0mac, thermal1);

    // subtract w0mac from umac
    Addw0(umac, w0mac, -1.);

    advect_time += ParallelDescriptor::second() - advect_time_start;
    ParallelDescriptor::ReduceRealMax(advect_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&advect_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    thermal_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 8a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(s2star, Tcoeff, hcoeff2, Xkcoeff2, pcoeff2);

        ThermalConduct(s1, s2, hcoeff1, Xkcoeff1, pcoeff1, hcoeff2, Xkcoeff2,
                       pcoeff2);
    }

    thermal_time += ParallelDescriptor::second() - thermal_time_start;
    ParallelDescriptor::ReduceRealMax(thermal_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&thermal_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Copy(s2[lev], s1[lev], Temp, Temp, 1, ng_s);
        MultiFab::Copy(s2[lev], s1[lev], Pi, Pi, 1, ng_s);
    }

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s2, p0_new, false);
    } else {
        TfromRhoH(s2, p0_new);
    }

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 9 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    react_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 9 : react state >>>" << std::endl;
    }

    React(s2, snew, rho_Hext, rho_omegadot, rho_Hnuc,
#if NAUX_NET > 0
          rho_auxdot,
#endif
          p0_new, 0.5 * dt, t_old + 0.5 * dt);

    react_time += ParallelDescriptor::second() - react_time_start;
    ParallelDescriptor::ReduceRealMax(react_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&react_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    if (evolve_base_state) {
        // compute beta0 and gamma1bar
        MakeGamma1bar(snew, gamma1bar_new, p0_new);

        base_time_start = ParallelDescriptor::second();

        MakeBeta0(beta0_new, rho0_new, p0_new, gamma1bar_new, grav_cell_new);

        base_time += ParallelDescriptor::second() - base_time_start;
        ParallelDescriptor::ReduceRealMax(
            base_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&base_time, 1,
                                  ParallelDescriptor::IOProcessorNumber());
    }

    beta0_nph.copy(0.5 * (beta0_old + beta0_new));

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 10 -- compute S^{n+1} for the final projection
    //////////////////////////////////////////////////////////////////////////////

    ndproj_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 10: make new S >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(snew, Tcoeff, hcoeff2, Xkcoeff2, pcoeff2);

        MakeExplicitThermal(thermal2, snew, Tcoeff, hcoeff2, Xkcoeff2, pcoeff2,
                            p0_new, temp_diffusion_formulation);
    }

    Make_S_cc(S_cc_new, delta_gamma1_term, delta_gamma1, snew, uold,
              rho_omegadot, rho_Hnuc, rho_Hext, thermal2, p0_new, gamma1bar_new,
              delta_gamma1_termbar);

    if (evolve_base_state) {
        Average(S_cc_new, Sbar, 0);

        // compute Sbar = Sbar + delta_gamma1_termbar
        if (use_delta_gamma1_term) {
            Sbar += delta_gamma1_termbar;
        }
    }

    // define dSdt = (S_cc_new - S_cc_old) / dt
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::LinComb(dSdt[lev], -1. / dt, S_cc_old[lev], 0, 1. / dt,
                          S_cc_new[lev], 0, 0, 1, 0);
    }

    ndproj_time += ParallelDescriptor::second() - ndproj_time_start;
    ParallelDescriptor::ReduceRealMax(ndproj_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ndproj_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    //////////////////////////////////////////////////////////////////////////////
    // STEP 11 -- update the velocity
    //////////////////////////////////////////////////////////////////////////////

    advect_time_start = ParallelDescriptor::second();

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 11: update and project new velocity >>>"
                << std::endl;
    }

    // Define rho at half time using the new rho from Step 8
    for (int lev = 0; lev <= finest_level; ++lev) {
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        rhohalf[lev].setVal(0.);
    }
    FillPatch(0.5 * (t_old + t_new), rhohalf, sold, snew, Rho, 0, 1, Rho,
              bcs_s);

    VelocityAdvance(rhohalf, umac, w0mac, w0_force_cart, rho0_nph,
                    grav_cell_nph, sponge);

    if (evolve_base_state && is_initIter) {
        // throw away w0 by setting w0 = w0_old
        w0.copy(w0_old);
    }

    int proj_type;

    advect_time += ParallelDescriptor::second() - advect_time_start;
    ParallelDescriptor::ReduceRealMax(advect_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&advect_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    ndproj_time_start = ParallelDescriptor::second();

    // Project the new velocity field
    if (is_initIter) {
        proj_type = pressure_iters_comp;

        // rhcc_for_nodalproj needs to contain
        // (beta0^nph S^1 - beta0^n S^0 ) / dt

        Vector<MultiFab> rhcc_for_nodalproj_old(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            rhcc_for_nodalproj_old[lev].define(grids[lev], dmap[lev], 1, 1);
            MultiFab::Copy(rhcc_for_nodalproj_old[lev], rhcc_for_nodalproj[lev],
                           0, 0, 1, 1);
        }

        MakeRHCCforNodalProj(rhcc_for_nodalproj, S_cc_new, Sbar, beta0_nph,
                             delta_gamma1_term);

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(rhcc_for_nodalproj[lev],
                               rhcc_for_nodalproj_old[lev], 0, 0, 1, 1);
            rhcc_for_nodalproj[lev].mult(1. / dt, 0, 1, 1);
        }

    } else {
        proj_type = regular_timestep_comp;

        MakeRHCCforNodalProj(rhcc_for_nodalproj, S_cc_new, Sbar, beta0_nph,
                             delta_gamma1_term);

        // compute delta_p_term = peos_new - p0_new (for RHS of projection)
        if (dpdt_factor > 0.) {
            // peos now holds "peos_new", the thermodynamic p computed from snew(rho,h,X)
            PfromRhoH(snew, snew, delta_p_term);

            // put p0_new on cart
            Put1dArrayOnCart(p0_new, p0_cart, false, false, bcs_f, 0);

            // compute delta_p_term = peos_new - p0_new
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Subtract(delta_p_term[lev], p0_cart[lev], 0, 0, 1, 0);
            }

            CorrectRHCCforNodalProj(rhcc_for_nodalproj, rho0_new, beta0_nph,
                                    gamma1bar_new, p0_new, delta_p_term);
        }
    }

    // call nodal projection
    NodalProj(proj_type, rhcc_for_nodalproj);

    beta0_nm1.copy(0.5 * (beta0_old + beta0_new));

    ndproj_time += ParallelDescriptor::second() - ndproj_time_start;
    ParallelDescriptor::ReduceRealMax(ndproj_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ndproj_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    misc_time_start = ParallelDescriptor::second();

    if (!is_initIter) {
        if (!fix_base_state) {
            // compute tempbar by "averaging"
            Average(snew, tempbar, Temp);
        }
    }

    Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
            << " DT = " << dt << std::endl;

    misc_time += ParallelDescriptor::second() - misc_time_start;
    ParallelDescriptor::ReduceRealMax(misc_time,
                                      ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&misc_time, 1,
                              ParallelDescriptor::IOProcessorNumber());

    // print wallclock time
    if (maestro_verbose > 0) {
        Print() << "Timing summary:\n";
        Print() << "Advection  :" << advect_time << " seconds\n";
        Print() << "MAC Proj   :" << macproj_time << " seconds\n";
        Print() << "Nodal Proj :" << ndproj_time << " seconds\n";
        if (use_thermal_diffusion) {
            Print() << "Thermal    :" << thermal_time << " seconds\n";
        }
        Print() << "Reactions  :" << react_time << " seconds\n";
        Print() << "Misc       :" << misc_time << " seconds\n";
        Print() << "Base State :" << base_time << " seconds\n";
    }
}
