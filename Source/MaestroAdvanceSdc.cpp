
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStepSDC (bool is_initIter) {

    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStepSDC()",AdvanceTimeStepSDC);

    // cell-centered MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab>           shat(finest_level+1);
    Vector<MultiFab>        rhohalf(finest_level+1);
    Vector<MultiFab>         cphalf(finest_level+1);
    Vector<MultiFab>         xihalf(finest_level+1);
    Vector<MultiFab>         macrhs(finest_level+1);
    Vector<MultiFab>         macphi(finest_level+1);
    Vector<MultiFab>       S_cc_nph(finest_level+1);
    Vector<MultiFab>   rho_omegadot(finest_level+1);
    Vector<MultiFab>       diff_old(finest_level+1);
    Vector<MultiFab>       diff_new(finest_level+1);
    Vector<MultiFab>       diff_hat(finest_level+1);
    Vector<MultiFab> diff_hterm_new(finest_level+1);
    Vector<MultiFab> diff_hterm_hat(finest_level+1);
    Vector<MultiFab>       rho_Hnuc(finest_level+1);
    Vector<MultiFab>       rho_Hext(finest_level+1);
    Vector<MultiFab>     sdc_source(finest_level+1);
    Vector<MultiFab>           aofs(finest_level+1);
    Vector<MultiFab>    intra_rhoh0(finest_level+1);
    
    Vector<MultiFab> delta_gamma1_term(finest_level+1);
    Vector<MultiFab>      delta_gamma1(finest_level+1);
    Vector<MultiFab>          peos_old(finest_level+1);
    Vector<MultiFab>      peosbar_cart(finest_level+1);
    Vector<MultiFab>           p0_cart(finest_level+1);
    Vector<MultiFab>      delta_p_term(finest_level+1);
    
    Vector<MultiFab>      Tcoeff1(finest_level+1);
    Vector<MultiFab>      hcoeff1(finest_level+1);
    Vector<MultiFab>     Xkcoeff1(finest_level+1);
    Vector<MultiFab>      pcoeff1(finest_level+1);
    Vector<MultiFab>      Tcoeff2(finest_level+1);
    Vector<MultiFab>      hcoeff2(finest_level+1);
    Vector<MultiFab>     Xkcoeff2(finest_level+1);
    Vector<MultiFab>      pcoeff2(finest_level+1);
    Vector<MultiFab>   scal_force(finest_level+1);
    Vector<MultiFab>    delta_chi(finest_level+1);
    Vector<MultiFab>       sponge(finest_level+1);
    Vector<MultiFab>         w0cc(finest_level+1);

    // face-centered in the dm-direction (planar only)
    Vector<MultiFab> etarhoflux_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  umac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sedge(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sflux(finest_level+1);

    ////////////////////////
    // needed for spherical routines only

    // cell-centered
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac_dummy(finest_level+1);

    // end spherical-only MultiFabs
    ////////////////////////

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    RealVector     grav_cell_nph    ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     rho0_nph         ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     p0_nph           ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     p0_minus_peosbar ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     peosbar          ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     w0_force_dummy   ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     Sbar             ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    BaseState<Real> beta0_nph (base_geom.max_radial_level+1, base_geom.nr_fine);
    RealVector     gamma1bar_nph    ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector delta_gamma1_termbar ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector delta_chi_w0_dummy   ( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    RealVector     delta_rhoh0      ( (base_geom.max_radial_level+1)*base_geom.nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    RealVector   w0_old             ( (base_geom.max_radial_level+1)*(base_geom.nr_fine+1) );
    BaseState<Real> rho0_pred_edge_dummy (base_geom.max_radial_level+1, base_geom.nr_fine+1);

    // make sure C++ is as efficient as possible with memory usage
    grav_cell_nph.shrink_to_fit();
    rho0_nph.shrink_to_fit();
    p0_nph.shrink_to_fit();
    p0_minus_peosbar.shrink_to_fit();
    peosbar.shrink_to_fit();
    w0_force_dummy.shrink_to_fit();
    Sbar.shrink_to_fit();
    gamma1bar_nph.shrink_to_fit();
    delta_gamma1_termbar.shrink_to_fit();
    w0_old.shrink_to_fit();
    delta_chi_w0_dummy.shrink_to_fit();
    delta_rhoh0.shrink_to_fit();

    int is_predictor;

    bool split_projection = true;

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl << std::endl;

    if (maestro_verbose > 0) {
        Print() << "Cell Count:" << std::endl;
        for (int lev=0; lev<=finest_level; ++lev) {
            Print() << "Level " << lev << ", " << CountCells(lev) << " cells" << std::endl;
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // cell-centered MultiFabs
        shat        [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
        rhohalf     [lev].define(grids[lev], dmap[lev],       1,    1);
        cphalf      [lev].define(grids[lev], dmap[lev],       1,    1);
        xihalf      [lev].define(grids[lev], dmap[lev], NumSpec,    1);
        macrhs      [lev].define(grids[lev], dmap[lev],       1,    0);
        macphi      [lev].define(grids[lev], dmap[lev],       1,    1);
        S_cc_nph    [lev].define(grids[lev], dmap[lev],       1,    0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec,    0);
        diff_old    [lev].define(grids[lev], dmap[lev],       1,    0);
        diff_new    [lev].define(grids[lev], dmap[lev],       1,    0);
        diff_hat    [lev].define(grids[lev], dmap[lev],       1,    0);
        diff_hterm_new[lev].define(grids[lev], dmap[lev],     1,    0);
        diff_hterm_hat[lev].define(grids[lev], dmap[lev],     1,    0);
        rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1,    0);
        rho_Hext    [lev].define(grids[lev], dmap[lev],       1,    0);
        sdc_source  [lev].define(grids[lev], dmap[lev],   Nscal,    0);
        aofs        [lev].define(grids[lev], dmap[lev],   Nscal,    0);
        intra_rhoh0 [lev].define(grids[lev], dmap[lev],       1,    0);
        delta_gamma1_term[lev].define(grids[lev], dmap[lev],  1,    0);
        delta_gamma1[lev].define(grids[lev], dmap[lev],       1,    0);
        peos_old    [lev].define(grids[lev], dmap[lev],       1,    0);
        peosbar_cart[lev].define(grids[lev], dmap[lev],       1,    0);
        p0_cart     [lev].define(grids[lev], dmap[lev],       1,    0);
        delta_p_term[lev].define(grids[lev], dmap[lev],       1,    0);
        Tcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
        hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
        Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
        pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
        Tcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
        hcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
        Xkcoeff2    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
        pcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
        if (ppm_trace_forces == 0) {
            scal_force  [lev].define(grids[lev], dmap[lev],   Nscal,    1);
        } else {
            // we need more ghostcells if we are tracing the forces
            scal_force  [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
        }
        delta_chi   [lev].define(grids[lev], dmap[lev],       1,    0);
        sponge      [lev].define(grids[lev], dmap[lev],       1,    0);
        w0cc    [lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);

        // face-centered in the dm-direction (planar only)
        AMREX_D_TERM(etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

        // face-centered arrays of MultiFabs
        AMREX_D_TERM(umac [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1,     1); ,
                     umac [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1,     1); ,
                     umac [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1,     1); );
        AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
                     sedge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
                     sedge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );
        AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
                     sflux[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
                     sflux[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );

        // initialize umac
        for (int d=0; d < AMREX_SPACEDIM; ++d)
            umac[lev][d].setVal(0.);

        // initialize intra_rhoh0
        intra_rhoh0[lev].setVal(0.);

        rho_Hext[lev].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    for (int lev=0; lev<=finest_level; ++lev) {
        w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        w0mac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        w0mac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        w0mac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
    }
#endif
    
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }
    
    // set etarhoflux_dummy to zero
    for (int lev=0; lev<=finest_level; ++lev) {
        etarhoflux_dummy[lev].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    // initialize MultiFabs and Vectors to ZERO
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            w0mac[lev][d].setVal(0.);
            w0mac_dummy[lev][d].setVal(0.);
        }
    }
#endif

    // initialize to zero
    std::fill(Sbar.begin(), Sbar.end(), 0.);
    std::fill(w0.begin()  , w0.end()  , 0.);

    // set dummy variables to zero
    std::fill(w0_force_dummy.begin()      , w0_force_dummy.end()      , 0.);
    rho0_pred_edge_dummy.setVal(0.0);

    // make the sponge for all levels
    if (do_sponge) {
        SpongeInit(rho0_old);
        MakeSponge(sponge);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 1 -- Compute advection velocities
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 1 : Compute advection velocities >>>" << std::endl;
    }

    if (t_old == 0.) {
        // this is either a pressure iteration or the first time step
        // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
        }
    }
    else {
        // set S_cc_nph = S_cc_old + (dt/2) * dSdt
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev],1.0,S_cc_old[lev],0,0.5*dt,dSdt[lev],0,0,1,0);
        }
    }
    // no ghost cells for S_cc_nph
    AverageDown(S_cc_nph,0,1);

    // compute p0_minus_peosbar = p0_old - peosbar (for making w0) and
    // compute delta_p_term = peos_old - p0_old (for RHS of projections)
    if (dpdt_factor > 0.0) {
        // peos_old (delta_p_term) now holds the thermodynamic p computed from sold(rho,h,X)
        PfromRhoH(sold,sold,peos_old);

        // compute peosbar = Avg(peos_old)
        Average(peos_old,peosbar,0);

        // compute p0_minus_peosbar = p0_old - peosbar
        for (int i=0; i<p0_minus_peosbar.size(); ++i) {
            p0_minus_peosbar[i] = p0_old[i] - peosbar[i];
        }

        // compute peosbar_cart from peosbar
        Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);

        // compute delta_p_term = peos_old - peosbar_cart
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(delta_p_term[lev],1.0,peos_old[lev],0,-1.0,peosbar_cart[lev],0,0,1,0);
        }
    }
    else {
        // these should have no effect if dpdt_factor <= 0
        std::fill(p0_minus_peosbar.begin(), p0_minus_peosbar.end(), 0.);
        for (int lev=0; lev<=finest_level; ++lev) {
            delta_p_term[lev].setVal(0.);
        }
    }
    
    if (evolve_base_state) {
        if (split_projection) {

            // compute Sbar = average(S_cc_nph)
            Average(S_cc_nph,Sbar,0);

            // save old-time value
            w0_old = w0;

            // compute w0, w0_force, and delta_chi_w0
            is_predictor = 1;
            Makew0(w0_old, w0_force_dummy, Sbar, rho0_old, 
                   rho0_old, p0_old, p0_old, gamma1bar_old, 
                   gamma1bar_old, p0_minus_peosbar, 
                   delta_chi_w0_dummy, dt, dtold, is_predictor);

            if (spherical) {
                // put w0 on Cartesian edges
                MakeW0mac(w0mac);
            }
        }
    }

    
    // compute unprojected MAC velocities
    is_predictor = 1;
    AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);

    for (int lev=0; lev<=finest_level; ++lev) {
        delta_chi[lev].setVal(0.);
        macphi   [lev].setVal(0.);
        delta_gamma1_term[lev].setVal(0.);
    }

    if (evolve_base_state && !split_projection) {
        for (int i=0; i<Sbar.size(); ++i) {
                Sbar[i] = 1.0/(gamma1bar_old[i]*p0_old[i]) * (p0_old[i] - p0_nm1[i])/dtold;
        }
    }

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    MakeRHCCforMacProj(macrhs,rho0_old,S_cc_nph,Sbar,beta0_old,delta_gamma1_term,
                       gamma1bar_old,p0_old,delta_p_term,delta_chi,is_predictor);

    if (spherical == 1 && evolve_base_state && split_projection) {
        // subtract w0mac from umac
        Addw0(umac,w0mac,-1.);
    }

    // wallclock time
    Real start_total_macproj = ParallelDescriptor::second();

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac,macphi,macrhs,beta0_old,is_predictor);

    // wallclock time
    Real end_total_macproj = ParallelDescriptor::second() - start_total_macproj;
    ParallelDescriptor::ReduceRealMax(end_total_macproj,ParallelDescriptor::IOProcessorNumber());

    
    if (spherical == 1 && evolve_base_state && split_projection) {
        // add w0mac back to umac
        Addw0(umac,w0mac,1.);
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- Predictor
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2 : Predictor >>>" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2A -- compute advective flux divergences
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2A : compute advective flux div >>>" << std::endl;
    }

    // no need to advect the base state density
    rho0_new = rho0_old;

    // set diff to zero
    for (int lev=0; lev<=finest_level; ++lev) {
        diff_old[lev].setVal(0.);
        diff_new[lev].setVal(0.);
        diff_hat[lev].setVal(0.);
    }

    // diff is the forcing for rhoh or temperature
    if (use_thermal_diffusion) {
        MakeThermalCoeffs(sold,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1);

        MakeExplicitThermal(diff_old,sold,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1,p0_old,
                            temp_diffusion_formulation);
    }

    // copy sold into shat 
    // temperature will be overwritten later after enthalpy advance
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(shat[lev],sold[lev],0,0,Nscal,0);
    }

    if (maestro_verbose >= 1) {
        Print() << "            :  density_advance >>>" << std::endl;
    }

    // set sedge and sflux to zero
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
            sedge[lev][idim].setVal(0.);
            sflux[lev][idim].setVal(0.);
        }
    }

    // advect rhoX and rho
    DensityAdvanceSDC(1,sold,shat,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac,rho0_pred_edge_dummy);
    

    if (evolve_base_state) {
        // correct the base state density by "averaging"
        Average(shat, rho0_new, Rho);
        compute_cutoff_coords(rho0_new.dataPtr());
        ComputeCutoffCoords(rho0_new);
    }

    // update grav_cell_new
    if (evolve_base_state) {
        MakeGravCell(grav_cell_new, rho0_new);
    }
    else {
        grav_cell_new = grav_cell_old;
   }

    // base state pressure update
    if (evolve_base_state) {

        // set new p0 through HSE
        p0_new = p0_old;

        EnforceHSE(rho0_new, p0_new, grav_cell_new);

        // compute p0_nph
        for (int i=0; i<p0_nph.size(); ++i) {
            p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
        }

        // hold dp0/dt in psi for enthalpy advance
        auto psi_arr = psi.array();
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                psi_arr(l,r) = (p0_new[l+(base_geom.max_radial_level+1)*r] - p0_old[l+(base_geom.max_radial_level+1)*r])/dt;
            }
        }

    } else {
        p0_new = p0_old;
        p0_nph = p0_old;
    }

    // base state enthalpy update
    if (evolve_base_state) {
        // compute rhoh0_old by "averaging"
        Average(sold, rhoh0_old, RhoH);
    }
    else {
        rhoh0_new = rhoh0_old;
    }

    if (maestro_verbose >= 1) {
        Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    EnthalpyAdvanceSDC(1,sold,shat,sedge,sflux,scal_force,umac,w0mac,diff_old);

    // base state enthalpy update
    if (evolve_base_state) {
        // compute rhoh0_new by "averaging"
        Average(shat, rhoh0_new, RhoH);

        // store (rhoh0_hat - rhoh0_old)/dt in delta_rhoh0
        for (int i=0; i<rhoh0_new.size(); ++i) {
            delta_rhoh0[i] = (rhoh0_new[i] - rhoh0_old[i])/dt;
        }
    }

    // extract aofs = (shat - sold) / dt - intra
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::LinComb(aofs[lev],1.0/dt,shat[lev],0,-1.0/dt,sold[lev],0,0,Nscal,0);
        MultiFab::Subtract(aofs[lev],intra[lev],0,0,Nscal,0);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2B (optional) -- compute diffusive flux divergence
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2B : compute diffusive flux div >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
        // 1 = predictor, 2 = corrector
        ThermalConductSDC(1,sold,shat,snew,p0_old,p0_new,hcoeff1,Xkcoeff1,pcoeff1,
                          hcoeff2,Xkcoeff2,pcoeff2);

        // note p0_new => p0_hat if evolve_base_state = T
        MakeExplicitThermal(diff_hat,shat,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1,p0_new,
                            temp_diffusion_formulation);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2C -- advance thermodynamic variables
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2C : advance thermo variables >>>" << std::endl;
    }

    // build sdc_source
    for (int lev=0; lev<=finest_level; ++lev) {
        sdc_source[lev].setVal(0.);
        MultiFab::Add(sdc_source[lev],aofs[lev],FirstSpec,FirstSpec,NumSpec,0);
        MultiFab::LinComb(sdc_source[lev],0.5,diff_old[lev],0,0.5,diff_hat[lev],0,RhoH,1,0);
        MultiFab::Add(sdc_source[lev],aofs[lev],RhoH,RhoH,1,0);
    }
    
    // wallclock time
    Real start_total_react = ParallelDescriptor::second();
    
    ReactSDC(sold,snew,rho_Hext,p0_old,dt,t_old,sdc_source,0,1);

    // wallclock time
    Real end_total_react = ParallelDescriptor::second() - start_total_react;
    ParallelDescriptor::ReduceRealMax(end_total_react,ParallelDescriptor::IOProcessorNumber());

    // extract IR =  [ (snew - sold)/dt - sdc_source ]

    for (int lev=0; lev<=finest_level; ++lev) {
        intra[lev].setVal(0.);
        // species source term
        MultiFab::LinComb(intra[lev],1.0/dt,snew[lev],FirstSpec,-1.0/dt,sold[lev],FirstSpec,FirstSpec,NumSpec,0);
        MultiFab::Subtract(intra[lev],sdc_source[lev],FirstSpec,FirstSpec,NumSpec,0);
        // enthalpy source term
        MultiFab::LinComb(intra[lev],1.0/dt,snew[lev],RhoH,-1.0/dt,sold[lev],RhoH,RhoH,1,0);
        MultiFab::Subtract(intra[lev],sdc_source[lev],RhoH,RhoH,1,0);
    }
    
    // massage the rhoh intra term into the proper form, depending on
    // what we are predicting.  Note: we do this before we deal with
    // the species terms, since some enthalpy types need this default
    // species intra.

    // first create rhohalf -- a lot of forms need this. 
    FillPatch(0.5*(t_old+t_new), rhohalf, sold, snew, Rho, 0, 1, Rho, bcs_s);

    if (evolve_base_state) {
        // update base state density and pressure
        Average(snew, rho0_new, Rho);
        compute_cutoff_coords(rho0_new.dataPtr());
        ComputeCutoffCoords(rho0_new);
        
        if (use_etarho) {
            // compute the new etarho
            if (!spherical) {
                MakeEtarho(etarhoflux_dummy);
            } else {
                MakeEtarhoSphr(sold,snew,umac,w0mac_dummy);
            }
        }

        MakeGravCell(grav_cell_new, rho0_new);

        EnforceHSE(rho0_new, p0_new, grav_cell_new);

        // compute p0_nph
        for (int i=0; i<p0_nph.size(); ++i) {
            p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
        }

        // hold dp0/dt in psi for Make_S_cc
        auto psi_arr = psi.array();
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                psi_arr(l,r) = (p0_new[l+(base_geom.max_radial_level+1)*r] - p0_old[l+(base_geom.max_radial_level+1)*r])/dt;
            }
        }

        // update base state enthalpy
        Average(snew, rhoh0_new, RhoH);

        // compute intra_rhoh0 = (rhoh0_new - rhoh0_old)/dt 
        //                       - (rhoh0_hat - rhoh0_old)/dt
        for (int i=0; i<rhoh0_new.size(); ++i) {
            delta_rhoh0[i] = (rhoh0_new[i] - rhoh0_old[i])/dt - delta_rhoh0[i];
        }
        Put1dArrayOnCart(delta_rhoh0, intra_rhoh0, 0, 0, bcs_f, 0);
    }

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(snew,p0_new);
    }
    else {
        TfromRhoH(snew,p0_new);
    }
    

    if (enthalpy_pred_type == predict_rhohprime) {

        // intra is only different from predict_rhoh if rhoh0 is not constant
        if (evolve_base_state) {            
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Subtract(intra[lev],intra_rhoh0[lev],0,RhoH,1,0);
            }
        }

    } else if (enthalpy_pred_type == predict_h) {

        // we want this in terms of h, not (rho h)
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Divide(intra[lev],rhohalf[lev],0,RhoH,1,0);
        }

    } else if ((enthalpy_pred_type == predict_T_then_rhohprime) ||
               (enthalpy_pred_type == predict_T_then_h)) {

        // for predict_T_*, the intra force needs to be in the temp_comp
        // slot, since temperature is what is predicted.

        // first make the thermodynamic coefficients at the half-time
        MakeIntraCoeffs(sold,snew,cphalf,xihalf);

        // overwrite intra(temp_comp).  We want to create
        // I_T = (1 / (rho c_p)) [ (rhoh_new - rhoh_old)/dt - A_rhoh -
        //     sum_k xi_k ( (rhoX_new - rhoX_old)/dt - A_rhoX ) ]
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(intra[lev],intra[lev],RhoH,Temp,1,0);
            
            for (int comp=0; comp<NumSpec; ++comp) {
                // multiple xi by intra and store in xi
                MultiFab::Multiply(xihalf[lev],intra[lev],FirstSpec+comp,comp,1,0);

                // subtract from intra temp
                MultiFab::Subtract(intra[lev],xihalf[lev],comp,Temp,1,0);
            }

            MultiFab::Divide(intra[lev],rhohalf[lev],0,Temp,1,0);
            MultiFab::Divide(intra[lev],cphalf[lev],0,Temp,1,0);
        }

    }     

    // for some species_pred_types, we need to make intra in terms of
    // X, NOT rhoX
    if ( (species_pred_type == predict_rhoprime_and_X) ||
         (species_pred_type == predict_rho_and_X) ) {

        for (int lev=0; lev<=finest_level; ++lev) {
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                MultiFab::Divide(intra[lev],rhohalf[lev],0,comp,1,0);
            }
        }
    }

    // compute new-time coefficients and diffusion term
    if (use_thermal_diffusion) {
        MakeThermalCoeffs(snew,Tcoeff2,hcoeff2,Xkcoeff2,pcoeff2);

        MakeExplicitThermal(diff_new,snew,Tcoeff2,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
                            temp_diffusion_formulation);
    }

    if (evolve_base_state) {
        // compute beta0 and gamma1bar
        MakeGamma1bar(snew,gamma1bar_new,p0_new);

        MakeBeta0(beta0_new, rho0_new, p0_new, gamma1bar_new, 
                  grav_cell_new);
    }
    else {
        // Just pass beta0 and gamma1bar through if not evolving base state
        beta0_new.copy(beta0_old);
        gamma1bar_new = gamma1bar_old;
    }

    for(int i=0; i<gamma1bar_nph.size(); ++i) {
        gamma1bar_nph[i] = 0.5*(gamma1bar_old[i]+gamma1bar_new[i]);
    }
    beta0_nph.copy(0.5*(beta0_old + beta0_new));

    //////////////////////////////////////////////////////////////////////////////
    // Corrector loop
    //////////////////////////////////////////////////////////////////////////////

    for (int misdc=0; misdc<sdc_iters; ++misdc) {
        
        //////////////////////////////////////////////////////////////////////////////
        // STEP 3 -- Update advection velocities
        //////////////////////////////////////////////////////////////////////////////

        if (sdc_couple_mac_velocity) {
            if (maestro_verbose >= 1) {
                Print() << "<<< STEP 3 : Update advection velocities (MISDC iter = " 
                        << misdc << ") >>>" << std::endl;
            }

            MakeReactionRates(rho_omegadot,rho_Hnuc,snew); 
            
            // compute S at cell-centers
            Make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot,rho_Hnuc,
                      rho_Hext,diff_new,p0_new,gamma1bar_new,delta_gamma1_termbar);
        
            // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
            }
            AverageDown(S_cc_nph,0,1);

            // and delta_p_term = peos_new - p0_new (for RHS of projection)
            if (dpdt_factor > 0.) {
                // peos_new now holds the thermodynamic p computed from snew(rho,h,X)
                PfromRhoH(snew,snew,delta_p_term);

                // compute peos_nph = (1/2)*(peos_old+peos_new)
                for (int lev=0; lev<=finest_level; ++lev) {
                    MultiFab::Add(delta_p_term[lev],peos_old[lev],0,0,1,0);
                    delta_p_term[lev].mult(0.5);
                }
                
                // compute peosbar = Avg(peos_new)
                Average(delta_p_term,peosbar,0);

                // compute p0_nph
                for (int i=0; i<p0_nph.size(); ++i) {
                    p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
                }
                
                // compute p0_minus_peosbar = p0_new - peosbar
                for (int i=0; i<p0_minus_peosbar.size(); ++i) {
                    p0_minus_peosbar[i] = p0_nph[i] - peosbar[i];
                }

                // compute peosbar_cart from peosbar
                Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);
                
                // compute delta_p_term = peos_new - peosbar_cart
                for (int lev=0; lev<=finest_level; ++lev) {
                    MultiFab::Subtract(delta_p_term[lev],peosbar_cart[lev],0,0,1,0);
                }
            }
            else {
                // these should have no effect if dpdt_factor <= 0
                std::fill(p0_minus_peosbar.begin(), p0_minus_peosbar.end(), 0.);
                for (int lev=0; lev<=finest_level; ++lev) {
                    delta_p_term[lev].setVal(0.);
                }
            }

            if (evolve_base_state) {
                
                if (split_projection) {
                    
                    // compute Sbar = average(S_cc_nph)
                    Average(S_cc_nph,Sbar,0);
                    
                    // compute Sbar = Sbar + delta_gamma1_termbar
                    if (use_delta_gamma1_term) {
                        for(int i=0; i<Sbar.size(); ++i) {
                            Sbar[i] += delta_gamma1_termbar[i];
                        }
                    }
                    
                    // compute w0, w0_force, and delta_chi_w0
                    is_predictor = 0;
                    Makew0(w0_old, w0_force_dummy, Sbar, rho0_old, 
                           rho0_new, p0_old, p0_new, gamma1bar_old, 
                           gamma1bar_new, p0_minus_peosbar, 
                           delta_chi_w0_dummy, dt, dtold, is_predictor);
                                
                    if (spherical) {
                        // put w0 on Cartesian edges
                        MakeW0mac(w0mac);
                    }
                }
            }

            
    
            // compute unprojected MAC velocities
            is_predictor = 0;
            AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);

            if (evolve_base_state && !split_projection) {
                for (int i=0; i<Sbar.size(); ++i) {
                    Sbar[i] = 1.0/(gamma1bar_nph[i]*p0_nph[i]) * (p0_nph[i] - p0_old[i])/dt;
                }
            }

            // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
            MakeRHCCforMacProj(macrhs,rho0_new,S_cc_nph,Sbar,beta0_nph,delta_gamma1_term,
                               gamma1bar_new,p0_new,delta_p_term,delta_chi,is_predictor);

            if (spherical == 1 && evolve_base_state && split_projection) {
                // subtract w0mac from umac
                Addw0(umac,w0mac,-1.);
            }

            // wallclock time
            Real start_total_macproj_corrector = ParallelDescriptor::second();

            // MAC projection
            // includes spherical option in C++ function
            MacProj(umac,macphi,macrhs,beta0_nph,is_predictor);

            // wallclock time
            Real end_total_macproj_corrector = ParallelDescriptor::second() - start_total_macproj_corrector;
            ParallelDescriptor::ReduceRealMax(end_total_macproj_corrector,ParallelDescriptor::IOProcessorNumber());


            if (spherical == 1 && evolve_base_state && split_projection) {
                // add w0mac back to umac
                Addw0(umac,w0mac,1.);
            }

        } // end sdc_couple_mac_velocity
    
    //////////////////////////////////////////////////////////////////////////////
    // STEP 4A -- compute advective flux divergences
    //////////////////////////////////////////////////////////////////////////////

        if (maestro_verbose >= 1) {
            Print() << "<<< STEP 4 : Corrector loop (MISDC iter = " 
                    << misdc << ") >>>" << std::endl;
            Print() << "<<< STEP 4A : compute advective flux div (SDC iter = " 
                    << misdc << ") >>>" << std::endl;
        }
        
        // no need to advect the base state density
        rho0_new = rho0_old;
        
        if (maestro_verbose >= 1) {
            Print() << "            :  density_advance >>>" << std::endl;
        }

        // advect rhoX, rho, and tracers
        DensityAdvanceSDC(2,sold,shat,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac,rho0_pred_edge_dummy);
    

        if (evolve_base_state) {
            // correct the base state density by "averaging"
            Average(shat, rho0_new, Rho);
            compute_cutoff_coords(rho0_new.dataPtr());
            ComputeCutoffCoords(rho0_new);
        }

        // update grav_cell_new, rho0_nph, grav_cell_nph
        if (evolve_base_state) {
            MakeGravCell(grav_cell_new, rho0_new);
            
            for(int i=0; i<rho0_nph.size(); ++i) {
                rho0_nph[i] = 0.5*(rho0_old[i]+rho0_new[i]);
            }
            
            MakeGravCell(grav_cell_nph, rho0_nph);
        } else {
            rho0_nph = rho0_old;
            grav_cell_nph = grav_cell_old;
        }
        
        // base state pressure update
        if (evolve_base_state) {
            
            // set new p0 through HSE
            p0_new = p0_old;

            EnforceHSE(rho0_new, p0_new, grav_cell_new);
            
            for (int i=0; i<p0_nph.size(); ++i) {
                p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
            }
            
            // hold dp0/dt in psi for enthalpy advance
            auto psi_arr = psi.array();
            for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
                for (auto r = 0; r < base_geom.nr_fine; ++r) {
                    psi_arr(l,r) = (p0_new[l+(base_geom.max_radial_level+1)*r] - p0_old[l+(base_geom.max_radial_level+1)*r])/dt;
                }
            }
        }
        
        // enthalpy update
        if (maestro_verbose >= 1) {
            Print() << "            : enthalpy_advance >>>" << std::endl;
        }
        
        EnthalpyAdvanceSDC(2,sold,shat,sedge,sflux,scal_force,umac,w0mac,diff_old);

        // base state enthalpy update
        if (evolve_base_state) {
            Average(shat, rhoh0_new, RhoH);

            // store (rhoh0_hat - rhoh0_old)/dt in delta_rhoh0
            for (int i=0; i<rhoh0_new.size(); ++i) {
                delta_rhoh0[i] = (rhoh0_new[i] - rhoh0_old[i])/dt;
            }
        }


        // extract aofs = (shat - sold) / dt - intra
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(aofs[lev],1.0/dt,shat[lev],0,-1.0/dt,sold[lev],0,0,Nscal,0);
            MultiFab::Subtract(aofs[lev],intra[lev],0,0,Nscal,0);
        }
        
    //////////////////////////////////////////////////////////////////////////////
    // STEP 4B (optional) -- compute diffusive flux divergences
    //////////////////////////////////////////////////////////////////////////////

        if (maestro_verbose >= 1) {
            Print() << "<<< STEP 4B : compute diffusive flux div (SDC iter = " 
                    << misdc << ") >>>" << std::endl;
        }

        // set diff_hterm to zero
        for (int lev=0; lev<=finest_level; ++lev) {
            diff_hterm_new[lev].setVal(0.);
            diff_hterm_hat[lev].setVal(0.);
        }

        if (use_thermal_diffusion) {
            // 1 = predictor, 2 = corrector
            ThermalConductSDC(2,sold,shat,snew,p0_old,p0_new,hcoeff1,Xkcoeff1,pcoeff1,
                              hcoeff2,Xkcoeff2,pcoeff2);

            // compute diff_hat using shat, p0_new, and new coefficients from previous iteration
            // note p0_new = p0_hat if evolve_base_state = T
            MakeExplicitThermal(diff_hat,shat,Tcoeff2,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
                                temp_diffusion_formulation);

            // compute only the h term in diff_hat and diff_new
            MakeExplicitThermalHterm(diff_hterm_hat,shat,hcoeff2);
            MakeExplicitThermalHterm(diff_hterm_new,snew,hcoeff2);
        }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4C -- advance thermodynamic variables
    //////////////////////////////////////////////////////////////////////////////

        if (maestro_verbose >= 1) {
            Print() << "<<< STEP 4C: advance thermo variables (MISDC iter = " 
                    << misdc << ") >>>" << std::endl;
        }

        // build sdc_source
        for (int lev=0; lev<=finest_level; ++lev) {
            sdc_source[lev].setVal(0.);
            MultiFab::Copy(sdc_source[lev],diff_old[lev],0,RhoH,1,0);
            MultiFab::Add(sdc_source[lev],diff_new[lev],0,RhoH,1,0);
            MultiFab::Add(sdc_source[lev],diff_hterm_hat[lev],0,RhoH,1,0);
            MultiFab::Subtract(sdc_source[lev],diff_hterm_new[lev],0,RhoH,1,0);
            sdc_source[lev].mult(0.5,RhoH,1,0);
            MultiFab::Add(sdc_source[lev],aofs[lev],0,0,Nscal,0);
        }
    
        // wallclock time
        Real start_total_react_corrector = ParallelDescriptor::second();
    
        ReactSDC(sold,snew,rho_Hext,p0_new,dt,t_old,sdc_source,misdc,sdc_iters);

        // wallclock time
        Real end_total_react_corrector = ParallelDescriptor::second() - start_total_react_corrector;
        ParallelDescriptor::ReduceRealMax(end_total_react_corrector,ParallelDescriptor::IOProcessorNumber());

        // extract IR =  [ (snew - sold)/dt - sdc_source ]
        for (int lev=0; lev<=finest_level; ++lev) {
            intra[lev].setVal(0.);
            // species source term
            MultiFab::LinComb(intra[lev],1.0/dt,snew[lev],FirstSpec,-1.0/dt,sold[lev],FirstSpec,FirstSpec,NumSpec,0);
            MultiFab::Subtract(intra[lev],sdc_source[lev],FirstSpec,FirstSpec,NumSpec,0);
            // enthalpy source term
            MultiFab::LinComb(intra[lev],1.0/dt,snew[lev],RhoH,-1.0/dt,sold[lev],RhoH,RhoH,1,0);
            MultiFab::Subtract(intra[lev],sdc_source[lev],RhoH,RhoH,1,0);
        }
        
        // massage the rhoh intra term into the proper form, depending on
        // what we are predicting.  Note: we do this before we deal with
        // the species terms, since some enthalpy types need this default
        // species intra.

        // create rhohalf -- a lot of forms need this. 
        FillPatch(0.5*(t_old+t_new), rhohalf, sold, snew, Rho, 0, 1, Rho, bcs_s);

        if (evolve_base_state) {
            // update base state density and pressure
            Average(snew, rho0_new, Rho);
            compute_cutoff_coords(rho0_new.dataPtr());
            ComputeCutoffCoords(rho0_new);

            if (use_etarho) {
                // compute the new etarho
                if (!spherical) {
                    MakeEtarho(etarhoflux_dummy);
                } else {
                    MakeEtarhoSphr(sold, snew, umac, w0mac_dummy);
                }
            }

            MakeGravCell(grav_cell_new, rho0_new);
            
            EnforceHSE(rho0_new, p0_new, grav_cell_new);

            // compute p0_nph
            for (int i=0; i<p0_nph.size(); ++i) {
                p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
            }

            // hold dp0/dt in psi for Make_S_cc
            auto psi_arr = psi.array();
            for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
                for (auto r = 0; r < base_geom.nr_fine; ++r) {
                    psi_arr(l,r) = (p0_new[l+(base_geom.max_radial_level+1)*r] - p0_old[l+(base_geom.max_radial_level+1)*r])/dt;
                }
            }
            
            // also update base state enthalpy
            Average(snew, rhoh0_new, RhoH);
            
            // compute intra_rhoh0 = (rhoh0_new - rhoh0_old)/dt 
            //                       - (rhoh0_hat - rhoh0_old)/dt
            for (int i=0; i<rhoh0_new.size(); ++i) {
                delta_rhoh0[i] = (rhoh0_new[i] - rhoh0_old[i])/dt - delta_rhoh0[i];
            }
            Put1dArrayOnCart(delta_rhoh0, intra_rhoh0, 0, 0, bcs_f, 0);
        }

        // now update temperature
        if (use_tfromp) {
            TfromRhoP(snew,p0_new);
        }
        else {
            TfromRhoH(snew,p0_new);
        }
            
        
        if (is_initIter) {
            for (int lev=0; lev<=finest_level; ++lev) {
                intra[lev].setVal(0.);
                intra_rhoh0[lev].setVal(0.);
            }
        }

        if (enthalpy_pred_type == predict_rhohprime) {
            
            // intra is only different from predict_rhoh if rhoh0 is not constant
            if (evolve_base_state) {        
                for (int lev=0; lev<=finest_level; ++lev) {
                    MultiFab::Subtract(intra[lev],intra_rhoh0[lev],0,RhoH,1,0);
                }
            }

        } else if (enthalpy_pred_type == predict_h) {

            // we want this in terms of h, not (rho h)
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Divide(intra[lev],rhohalf[lev],0,RhoH,1,0);
            }

        } else if ((enthalpy_pred_type == predict_T_then_rhohprime) ||
                   (enthalpy_pred_type == predict_T_then_h)) {

        // for predict_T_*, the intra force needs to be in the temp_comp
        // slot, since temperature is what is predicted.

            // first make the thermodynamic coefficients at the half-time
            MakeIntraCoeffs(sold,snew,cphalf,xihalf);

            // overwrite intra(temp_comp).  We want to create
            // I_T = (1 / (rho c_p)) [ (rhoh_new - rhoh_old)/dt - A_rhoh -
            //     sum_k xi_k ( (rhoX_new - rhoX_old)/dt - A_rhoX ) ]
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Copy(intra[lev],intra[lev],RhoH,Temp,1,0);
                
                for (int comp=0; comp<NumSpec; ++comp) {
                    // multiple xi by intra and store in xi
                    MultiFab::Multiply(xihalf[lev],intra[lev],FirstSpec+comp,comp,1,0);
                    
                    // subtract from intra temp
                    MultiFab::Subtract(intra[lev],xihalf[lev],comp,Temp,1,0);
                }
                
                MultiFab::Divide(intra[lev],rhohalf[lev],0,Temp,1,0);
                MultiFab::Divide(intra[lev],cphalf[lev],0,Temp,1,0);
            }
            
        }     
        
        // for some species_pred_types, we need to make intra in terms of
        // X, NOT rhoX
        if ( (species_pred_type == predict_rhoprime_and_X) ||
             (species_pred_type == predict_rho_and_X) ) {

            for (int lev=0; lev<=finest_level; ++lev) {
                for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                    MultiFab::Divide(intra[lev],rhohalf[lev],0,comp,1,0);
                }
            }
        }

        // compute new-time coefficients and diffusion term
        if (use_thermal_diffusion) {
            MakeThermalCoeffs(snew,Tcoeff2,hcoeff2,Xkcoeff2,pcoeff2);

            MakeExplicitThermal(diff_new,snew,Tcoeff2,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
                                temp_diffusion_formulation);
        }

        if (evolve_base_state) {
            // compute beta0 and gamma1bar
            MakeGamma1bar(snew,gamma1bar_new,p0_new);

            MakeBeta0(beta0_new, rho0_new, p0_new, gamma1bar_new, 
                      grav_cell_new);
        }
        else {
            // Just pass beta0 and gamma1bar through if not evolving base state
            beta0_new = beta0_old;
            gamma1bar_new = gamma1bar_old;
      }

        for (int i=0; i<gamma1bar_nph.size(); ++i) {
            gamma1bar_nph[i] = 0.5*(gamma1bar_old[i]+gamma1bar_new[i]);
        }
        beta0_nph.copy(0.5*(beta0_old + beta0_new));
        
    } // end loop over misdc iterations
    
    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- Advance velocity and dynamic pressure
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 5 : Advance velocity and dynamic pressure >>>" << std::endl;
    }

    MakeReactionRates(rho_omegadot,rho_Hnuc,snew); 
    
    Make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot,rho_Hnuc,
              rho_Hext,diff_new,p0_new,gamma1bar_new,delta_gamma1_termbar);

    if (evolve_base_state) {

        if (split_projection) {

            // compute Sbar = average(S_cc_new)
            Average(S_cc_new,Sbar,0);

            // compute Sbar = Sbar + delta_gamma1_termbar
            if (use_delta_gamma1_term) {
                for(int i=0; i<Sbar.size(); ++i) {
                    Sbar[i] += delta_gamma1_termbar[i];
                }
            }

            // compute w0, w0_force, and delta_chi_w0
            is_predictor = 0;
            Makew0(w0_old, w0_force_dummy, Sbar, rho0_new, 
                   rho0_new, p0_new, p0_new, gamma1bar_new, 
                   gamma1bar_new, p0_minus_peosbar, 
                   delta_chi_w0_dummy, dt, dtold, is_predictor);

            if (spherical) {
                // put w0 on Cartesian cell-centers
                Put1dArrayOnCart(w0, w0cc, 1, 1, bcs_u, 0, 1);
            }
        }
    }

    
    // define dSdt = (S_cc_new - S_cc_old) / dt
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::LinComb(dSdt[lev],-1./dt,S_cc_old[lev],0,1./dt,S_cc_new[lev],0,0,1,0);
    }
    
    // Define rho at half time using the new rho from Step 4
    FillPatch(0.5*(t_old+t_new), rhohalf, sold, snew, Rho, 0, 1, Rho, bcs_s);

    VelocityAdvance(rhohalf,umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy,
                    rho0_nph,grav_cell_nph,sponge);

    if (evolve_base_state && is_initIter) {
        // throw away w0 by setting w0 = w0_old
        w0 = w0_old;
    }

    if (spherical == 1 && evolve_base_state && split_projection) {
        // subtract w0 from uold and unew for nodal projection
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(uold[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
            MultiFab::Subtract(unew[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
        }
    }

    if (evolve_base_state && !split_projection) {
        for (int i=0; i<Sbar.size(); ++i) {
            Sbar[i] = (p0_new[i] - p0_old[i])/(dt*gamma1bar_new[i]*p0_new[i]);
        }
    }

    int proj_type;

    // Project the new velocity field
    if (is_initIter) {

        proj_type = pressure_iters_comp;

        // rhcc_for_nodalproj needs to contain
        // (beta0^nph S^1 - beta0^n S^0 ) / dt

        Vector<MultiFab> rhcc_for_nodalproj_old(finest_level+1);
        for (int lev=0; lev<=finest_level; ++lev) {
            rhcc_for_nodalproj_old[lev].define(grids[lev], dmap[lev], 1, 1);
            MultiFab::Copy(rhcc_for_nodalproj_old[lev], rhcc_for_nodalproj[lev], 0, 0, 1, 1);
        }

        MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph,delta_gamma1_term);

        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Subtract(rhcc_for_nodalproj[lev], rhcc_for_nodalproj_old[lev], 0, 0, 1, 1);
            rhcc_for_nodalproj[lev].mult(1./dt,0,1,1);
        }

    }
    else {

        proj_type = regular_timestep_comp;

        MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph,delta_gamma1_term);

        // compute delta_p_term = peos_new - p0_new (for RHS of projection)
        if (dpdt_factor > 0.) {
            // peos_new now holds the thermodynamic p computed from snew(rho h X)
            PfromRhoH(snew,snew,delta_p_term);

            // compute peosbar = Avg(peos_new)
            Average(delta_p_term,peosbar,0);

            // no need to compute p0_minus_peosbar since make_w0 is not called after here

            // compute peosbar_cart from peosbar
            Put1dArrayOnCart(peosbar, peosbar_cart, 0, 0, bcs_f, 0);

            // compute delta_p_term = peos_new - peosbar_cart
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Subtract(delta_p_term[lev],peosbar_cart[lev],0,0,1,0);
            }

            CorrectRHCCforNodalProj(rhcc_for_nodalproj,rho0_new,beta0_nph,gamma1bar_new,
                                    p0_new,delta_p_term);
        }
    }

    // wallclock time
    const Real start_total_nodalproj = ParallelDescriptor::second();

    // call nodal projection
    NodalProj(proj_type,rhcc_for_nodalproj,0,false);

    // wallclock time
    Real end_total_nodalproj = ParallelDescriptor::second() - start_total_nodalproj;
    ParallelDescriptor::ReduceRealMax(end_total_nodalproj,ParallelDescriptor::IOProcessorNumber());

    if (spherical == 1 && evolve_base_state && split_projection) {
        // add w0 back to unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(unew[lev],w0cc[lev],0,0,AMREX_SPACEDIM,0);
        }
        AverageDown(unew,0,AMREX_SPACEDIM);
        FillPatch(t_new, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);
    }

    beta0_nm1.copy(0.5*(beta0_old + beta0_new));

    if (!is_initIter) {
        if (!fix_base_state) {
            // compute tempbar by "averaging"
            Average(snew,tempbar,Temp);
        }
    }

    Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
            << " DT = " << dt << std::endl;

    // print wallclock time
    if (maestro_verbose > 0) {
        Print() << "Time to solve mac proj   : " << end_total_macproj << '\n';
        Print() << "Time to solve nodal proj : " << end_total_nodalproj << '\n';
        Print() << "Time to solve reactions  : " << end_total_react << '\n';
    }

}
