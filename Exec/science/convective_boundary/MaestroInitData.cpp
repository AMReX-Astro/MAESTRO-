
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto s0_arr = s0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
#if NAUX_NET > 0
        for (auto comp = 0; comp < NumAux; ++comp) {
            scal(i, j, k, FirstAux + comp) = s0_arr(lev, r, FirstAux + comp);
        }
#endif
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr not implemented");
}