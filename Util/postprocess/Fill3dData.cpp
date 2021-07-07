/*
   A collection of routines for mapping 1D arrays onto multi-D cartesian MultiFabs
 */

#include <Postprocess.H>
#include "AMReX_MultiFabUtil.H"

using namespace amrex;

void Postprocess::PutDataOnFaces(
    const Vector<MultiFab>& s_cc,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& face,
    const bool harmonic_avg) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::PutDataOnFaces()", PutDataOnFaces);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // need one cell-centered MF for the MFIter
        const MultiFab& scc_mf = s_cc[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scc_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<const Real> scc = s_cc[lev].array(mfi);
            const Array4<Real> facex = face[lev][0].array(mfi);
            const Array4<Real> facey = face[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> facez = face[lev][2].array(mfi);
#endif

            if (harmonic_avg) {
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    Real denom = scc(i, j, k) + scc(i - 1, j, k);
                    Real prod = scc(i, j, k) * scc(i - 1, j, k);

                    if (denom != 0.0) {
                        facex(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facex(i, j, k) = 0.5 * denom;
                    }
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    Real denom = scc(i, j, k) + scc(i, j - 1, k);
                    Real prod = scc(i, j, k) * scc(i, j - 1, k);

                    if (denom != 0.0) {
                        facey(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facey(i, j, k) = 0.5 * denom;
                    }
                });
#if (AMREX_SPACEDIM == 3)
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    Real denom = scc(i, j, k) + scc(i, j, k - 1);
                    Real prod = scc(i, j, k) * scc(i, j, k - 1);

                    if (denom != 0.0) {
                        facez(i, j, k) = 2.0 * prod / denom;
                    } else {
                        facez(i, j, k) = 0.5 * denom;
                    }
                });
#endif
            } else {
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    facex(i, j, k) = 0.5 * (scc(i, j, k) + scc(i - 1, j, k));
                });
                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    facey(i, j, k) = 0.5 * (scc(i, j, k) + scc(i, j - 1, k));
                });
#if (AMREX_SPACEDIM == 3)
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    facez(i, j, k) = 0.5 * (scc(i, j, k) + scc(i, j, k - 1));
                });
#endif
            }
        }
    }

    // Make sure that the fine edges average down onto the coarse edges (edge_restriction)
    AverageDownFaces(face);
}

void Postprocess::AverageDownFaces(
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& edge) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::AverageDownFaces()", AverageDownFaces);

    for (int lev = finest_level - 1; lev >= 0; --lev) {
        Vector<const MultiFab*> edge_f(AMREX_SPACEDIM);
        Vector<MultiFab*> edge_c(AMREX_SPACEDIM);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            edge_f[dir] = &(edge[lev + 1][dir]);
            edge_c[dir] = &(edge[lev][dir]);
        }

        amrex::average_down_faces(edge_f, edge_c, IntVect(2));
    }
}
