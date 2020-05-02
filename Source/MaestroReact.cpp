
#include <Maestro.H>
#include <Maestro_F.H>
#include <Maestro_F.H>

using namespace amrex;

// compute heating term, rho_Hext, then
// react the state over dt_in and update rho_omegadot, rho_Hnuc
void
Maestro::React (const Vector<MultiFab>& s_in,
                Vector<MultiFab>& s_out,
                Vector<MultiFab>& rho_Hext,
                Vector<MultiFab>& rho_omegadot,
                Vector<MultiFab>& rho_Hnuc,
                const RealVector& p0,
                const Real dt_in,
                const Real time_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::React()",React);

    // external heating
    if (do_heating) {

        // computing heating term
        MakeHeating(rho_Hext,s_in);

        // if we aren't burning, then we should just copy the old state to the
        // new and only update the rhoh component with the heating term
        if (!do_burning) {
            for (int lev=0; lev<=finest_level; ++lev) {
                // copy s_in to s_out
                MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);

                // add in the heating term, s_out += dt_in * rho_Hext
                MultiFab::Saxpy(s_out[lev],dt_in,rho_Hext[lev],0,RhoH,1,0);
            }
        }
    }
    else {
        // not heating, so we zero rho_Hext
        for (int lev=0; lev<=finest_level; ++lev) {
            rho_Hext[lev].setVal(0.);
        }
    }

    // apply burning term
    if (do_burning) {
#ifndef SDC
        // do the burning, update rho_omegadot and rho_Hnuc
        // we pass in rho_Hext so that we can add it to rhoh in case we applied heating
        Burner(s_in,s_out,rho_Hext,rho_omegadot,rho_Hnuc,p0,dt_in,time_in);
#endif
        // pass temperature through for seeding the temperature update eos call
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],Temp,Temp,1,0);
        }
    }
    else {
        // not burning, so we zero rho_omegadot and rho_Hnuc
        for (int lev=0; lev<=finest_level; ++lev) {
            rho_omegadot[lev].setVal(0.);
            rho_Hnuc[lev].setVal(0.);
        }

    }

    // if we aren't doing any heating/burning, then just copy the old to the new
    if (!do_heating && !do_burning) {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);
        }
    }

    // average down and fill ghost cells
    AverageDown(s_out,0,Nscal);
    FillPatch(t_old,s_out,s_out,s_out,0,0,Nscal,0,bcs_s);

    // average down (no ghost cells)
    AverageDown(rho_Hext,0,1);
    AverageDown(rho_omegadot,0,NumSpec);
    AverageDown(rho_Hnuc,0,1);

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s_out,p0);
    }
    else {
        TfromRhoH(s_out,p0);
    }

}

// SDC subroutines
// compute heating term, rho_Hext, then
// react the state over dt_in and update s_out
void
Maestro::ReactSDC (const Vector<MultiFab>& s_in,
                   Vector<MultiFab>& s_out,
                   Vector<MultiFab>& rho_Hext,
                   const RealVector& p0,
                   const Real dt_in,
                   const Real time_in,
                   Vector<MultiFab>& source,
                   const int sdc_iter,
                   const int number_sdc_iterations)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ReactSDC()",ReactSDC);

    // external heating
    if (do_heating) {

        // computing heating term
        MakeHeating(rho_Hext,s_in);

        if (do_burning) {
            // add to source for enthalpy
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Add(source[lev],rho_Hext[lev],0,RhoH,1,0);
            }
        } else {
            // if we aren't burning, then we should just copy the old state to the
            // new and only update the rhoh component with the heating term
            // s_out = s_in + dt_in * rho_Hext
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);
                MultiFab::Saxpy(s_out[lev],dt_in,rho_Hext[lev],0,RhoH,1,0);
            }
        }
    }
    else {
        // not heating, so we zero rho_Hext
        for (int lev=0; lev<=finest_level; ++lev) {
            rho_Hext[lev].setVal(0.);
        }
    }

    // apply burning term
    if (do_burning) {
#ifdef SDC
        // do the burning, update s_out
        Burner(s_in,s_out,p0,dt_in,time_in,source,sdc_iter,number_sdc_iterations);
#endif
        // pass temperature through for seeding the temperature update eos call
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],Temp,Temp,1,0);
        }
    }

    // if we aren't doing any heating/burning, then just copy the old to the new
    if (!do_heating && !do_burning) {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);
        }
    }

    // average down and fill ghost cells
    AverageDown(s_out,0,Nscal);
    FillPatch(t_old,s_out,s_out,s_out,0,0,Nscal,0,bcs_s);

    // average down (no ghost cells)
    AverageDown(rho_Hext,0,1);

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s_out,p0);
    }
    else {
        TfromRhoH(s_out,p0);
    }

}


#ifndef SDC
void Maestro::Burner(const Vector<MultiFab>& s_in,
                     Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot,
                     Vector<MultiFab>& rho_Hnuc,
                     const RealVector& p0,
                     const Real dt_in,
                     const Real time_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Burner()",Burner);

    // Put tempbar_init on cart
    Vector<MultiFab> tempbar_init_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            tempbar_init_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            tempbar_init_cart[lev].setVal(0.);
        }

        if (drive_initial_convection == 1) {
            Put1dArrayOnCart(tempbar_init,tempbar_init_cart,0,0,bcs_f,0);
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab&         s_in_mf =         s_in[lev];
        MultiFab&        s_out_mf =        s_out[lev];
        const MultiFab&     rho_Hext_mf =     rho_Hext[lev];
        MultiFab& rho_omegadot_mf = rho_omegadot[lev];
        MultiFab&     rho_Hnuc_mf =     rho_Hnuc[lev];
        const MultiFab& tempbar_cart_mf = tempbar_init_cart[lev];

        // create mask assuming refinement ratio = 2
        int finelev = lev+1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in_mf, fba, IntVect(2));


        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s_in_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            int use_mask = !(lev==finest_level);

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 1) {
#pragma gpu box(tileBox)
                burner_loop_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                                 BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_Hext_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_omegadot_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                                 BL_TO_FORTRAN_ANYD(tempbar_cart_mf[mfi]), dt_in, time_in, 
                                 BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask);
            } else {
#pragma gpu box(tileBox)
                burner_loop(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                            lev,
                            BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(rho_Hext_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(rho_omegadot_mf[mfi]),
                            BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                            tempbar_init.dataPtr(), dt_in, time_in, 
                            BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask);
            }
        }
    }
}

#else
// SDC burner
void Maestro::Burner(const Vector<MultiFab>& s_in,
                     Vector<MultiFab>& s_out,
                     const RealVector& p0,
                     const Real dt_in,
                     const Real time_in,
                     const Vector<MultiFab>& source,
                     const int sdc_iter,
                     const int number_sdc_iterations)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::BurnerSDC()",BurnerSDC);

    // Put tempbar_init on cart
    Vector<MultiFab> p0_cart(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].setVal(0.);
        }

        if (drive_initial_convection == 1) {
            Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab&    s_in_mf =    s_in[lev];
              MultiFab&   s_out_mf =   s_out[lev];
        const MultiFab& p0_cart_mf = p0_cart[lev];
        const MultiFab&  source_mf =  source[lev];
        
        // create mask assuming refinement ratio = 2
        int finelev = lev+1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(s_in_mf, fba, IntVect(2));
        

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s_in_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            int use_mask = !(lev==finest_level);

            // call fortran subroutine
            
            if (spherical == 1) {
#pragma gpu box(tileBox)
                burner_loop_sphr(AMREX_INT_ANYD(tileBox.loVect()), 
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(source_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(p0_cart_mf[mfi]), dt_in, time_in,
                    BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask,
                    sdc_iter, number_sdc_iterations);
            } else {
#pragma gpu box(tileBox)
                burner_loop(AMREX_INT_ANYD(tileBox.loVect()), 
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    lev,
                    BL_TO_FORTRAN_ANYD(s_in_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(s_out_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(source_mf[mfi]), 
                    p0.dataPtr(), dt_in, time_in,
                    BL_TO_FORTRAN_ANYD(mask[mfi]), use_mask,
                    sdc_iter, number_sdc_iterations);
            }
        }
    }
}
#endif


// compute heating terms, rho_omegadot and rho_Hnuc
void
Maestro::MakeReactionRates (Vector<MultiFab>& rho_omegadot,
                            Vector<MultiFab>& rho_Hnuc,
                            const Vector<MultiFab>& scal)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeReactionRates()",MakeReactionRates);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
              MultiFab& rho_omegadot_mf = rho_omegadot[lev];
              MultiFab&     rho_Hnuc_mf = rho_Hnuc[lev];
        const MultiFab&         scal_mf =     scal[lev];

        if (!do_burning) {
            rho_omegadot[lev].setVal(0.);
            rho_Hnuc[lev].setVal(0.);
        } else {

            // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                
                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                instantaneous_reaction_rates(AMREX_INT_ANYD(tileBox.loVect()), 
                         AMREX_INT_ANYD(tileBox.hiVect()),
                                             BL_TO_FORTRAN_ANYD(rho_omegadot_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(rho_Hnuc_mf[mfi]),
                                             BL_TO_FORTRAN_ANYD(scal_mf[mfi]));
            }
        }
    }
    
}
