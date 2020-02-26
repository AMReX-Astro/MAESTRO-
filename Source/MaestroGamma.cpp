
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::MakeGamma1bar (const Vector<MultiFab>& scal,
                        RealVector& gamma1bar,
                        const RealVector& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGamma1bar()",MakeGamma1bar);

    Vector<MultiFab> gamma1(finest_level+1);
    Vector<MultiFab> p0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        gamma1[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1[lev].setVal(0.);

        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        p0_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(p0, p0_cart, 0, 0, bcs_f, 0);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& gamma1_mf = gamma1[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(gamma1[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> gamma1_arr = gamma1[lev].array(mfi);
            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                eos_t eos_state;

                eos_state.rho = scal_arr(i,j,k,Rho);

                if (use_pprime_in_tfromp) {
                    eos_state.p = p0_arr(i,j,k) + scal_arr(i,j,k,Pi);
                } else {
                    eos_state.p = p0_arr(i,j,k);
                }
                eos_state.T = scal_arr(i,j,k,Temp);
                for (auto n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = scal_arr(i,j,k,FirstSpec+n) / eos_state.rho;
                }

                // dens, pres, and xmass are inputs
                eos(eos_input_rp, eos_state);

                gamma1_arr(i,j,k) = eos_state.gam1;
            });
        }
    }

    // average fine data onto coarser cells
    AverageDown(gamma1, 0, 1);

    // call average to create gamma1bar
    Average(gamma1, gamma1bar, 0);
}
