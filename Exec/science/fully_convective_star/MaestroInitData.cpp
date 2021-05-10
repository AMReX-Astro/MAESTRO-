
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    Abort("Planar InitLevelData not implemented.");
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero
        ParallelFor(tileBox, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        vel_arr(i, j, k, n) = 0.0;
                    });

        ParallelFor(tileBox, Nscal,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        scal_arr(i, j, k, n) = 0.0;
                    });
    }

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    MultiFab p0_cart(scal.boxArray(), scal.DistributionMap(), 1, 0);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal
    MultiFab temp_mf(scal.boxArray(), scal.DistributionMap(), 1, 0);

    BaseState<Real> temp_vec(base_geom.max_radial_level + 1, base_geom.nr_fine);
    auto temp_arr = temp_vec.array();

    const auto s0_arr = s0_init.const_array();

    // initialize temperature
    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            temp_arr(l, r) = s0_arr(l, r, Temp);
        }
    }

    Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_f, 0);
    MultiFab::Copy(scal, temp_mf, 0, Temp, 1, 0);

    // initialize p0_cart
    Put1dArrayOnCart(lev, p0_init, p0_cart, 0, 0);

    // initialize species
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                temp_arr(l, r) = s0_arr(l, r, FirstSpec + comp);
            }
        }
        Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec + comp);
        MultiFab::Copy(scal, temp_mf, 0, FirstSpec + comp, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        const Array4<const Real> p0_arr = p0_cart.array(mfi);

        // initialize rho as sum of partial densities rho*X_i
        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i, j, k, Rho) += scal_arr(i, j, k, FirstSpec + comp);
            }

            // initialize (rho h) and T using the EOS
            eos_t eos_state;
            eos_state.T = scal_arr(i, j, k, Temp);
            eos_state.p = p0_arr(i, j, k);
            eos_state.rho = scal_arr(i, j, k, Rho);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] =
                    scal_arr(i, j, k, FirstSpec + comp) / eos_state.rho;
            }

            eos(eos_input_rp, eos_state);

            scal_arr(i, j, k, RhoH) = eos_state.rho * eos_state.h;
            scal_arr(i, j, k, Temp) = eos_state.T;
        });

        if (perturb_model) {
            // random numbers between -1 and 1
            GpuArray<Real, 27> alpha;
            GpuArray<Real, 27> beta;
            GpuArray<Real, 27> gamma;

            // random numbers between 0 and 2*pi
            GpuArray<Real, 27> phix;
            GpuArray<Real, 27> phiy;
            GpuArray<Real, 27> phiz;

            // compute the norm of k
            GpuArray<Real, 27> normk;

            for (auto k = 1; k <= 3; ++k) {
                for (auto j = 1; j <= 3; ++j) {
                    for (auto i = 1; i <= 3; ++i) {
                        const int n = i - 1 + 3 * (j - 1 + 3 * (k - 1));

                        alpha[n] = std::pow(0.5, i) * std::pow(0.7, j) *
                                   std::pow(0.3, k) * std::pow(-1.0, i);
                        beta[n] = std::pow(0.5, i) * std::pow(0.3, j) *
                                  std::pow(0.7, k) * std::pow(-1.0, j);
                        gamma[n] = std::pow(0.3, i) * std::pow(0.5, j) *
                                   std::pow(0.7, k) * std::pow(-1.0, k);

                        phix[n] = std::pow(0.3, i) * std::pow(0.7, j) *
                                  std::pow(0.5, k);
                        phix[n] *= 2.0 * M_PI;
                        phiy[n] = std::pow(0.7, i) * std::pow(0.3, j) *
                                  std::pow(0.5, k);
                        phiy[n] *= 2.0 * M_PI;
                        phiz[n] = std::pow(0.7, i) * std::pow(0.5, j) *
                                  std::pow(0.3, k);
                        phiz[n] *= 2.0 * M_PI;

                        normk[n] =
                            std::sqrt(Real(i * i) + Real(j * j) + Real(k * k));

                        // alpha[n] = 2.0 * amrex::Random() - 1.0;
                        // beta[n] = 2.0 * amrex::Random() - 1.0;
                        // gamma[n] = 2.0 * amrex::Random() - 1.0;

                        // phix[n] = 2.0 * M_PI * amrex::Random();
                        // phiy[n] = 2.0 * M_PI * amrex::Random();
                        // phiz[n] = 2.0 * M_PI * amrex::Random();
                    }
                }
            }

            const auto prob_lo = geom[lev].ProbLoArray();
            const auto prob_hi = geom[lev].ProbHiArray();
            const auto dx = geom[lev].CellSizeArray();

            const auto center_p = center;

            const auto velpert_scale_loc = velpert_scale;
            const auto velpert_amplitude_loc = velpert_amplitude;
            const auto velpert_steep_loc = velpert_steep;
            const auto velpert_r_outer = velpert_radius;
            const Real velpert_r_inner = 0.0;

            // now do the big loop over all the points in the domain
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const Real x =
                    prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                const Real y =
                    prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                const Real z =
                    prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                // set perturbational velocity to zero
                Real vpert[3];

                for (auto n = 0; n < 3; ++n) {
                    vpert[n] = 0.0;
                }

                // compute distance to the center of the star
                Real rloc = std::sqrt(x * x + y * y + z * z);

                // loop over the 27 combinations of Fourier components
                for (auto kk = 1; kk <= 3; ++kk) {
                    for (auto jj = 1; jj <= 3; ++jj) {
                        for (auto ii = 1; ii <= 3; ++ii) {
                            // array index
                            const int n = ii - 1 + 3 * (jj - 1 + 3 * (kk - 1));
                            // compute cosines and sines
                            Real cx = std::cos(2.0 * M_PI * Real(ii) *
                                                   (x + center_p[0]) /
                                                   velpert_scale_loc +
                                               phix[n]);
                            Real cy = std::cos(2.0 * M_PI * Real(jj) *
                                                   (y + center_p[1]) /
                                                   velpert_scale_loc +
                                               phiy[n]);
                            Real cz = std::cos(2.0 * M_PI * Real(kk) *
                                                   (z + center_p[2]) /
                                                   velpert_scale_loc +
                                               phiz[n]);

                            Real sx = std::sin(2.0 * M_PI * Real(ii) *
                                                   (x + center_p[0]) /
                                                   velpert_scale_loc +
                                               phix[n]);
                            Real sy = std::sin(2.0 * M_PI * Real(jj) *
                                                   (y + center_p[1]) /
                                                   velpert_scale_loc +
                                               phiy[n]);
                            Real sz = std::sin(2.0 * M_PI * Real(kk) *
                                                   (z + center_p[2]) /
                                                   velpert_scale_loc +
                                               phiz[n]);

                            // compute contribution from perturbation velocity from each mode
                            vpert[0] += (-gamma[n] * Real(jj) * cx * cz * sy +
                                         beta[n] * Real(kk) * cx * cy * sz) /
                                        normk[n];

                            vpert[1] += (gamma[n] * Real(ii) * cy * cz * sx -
                                         alpha[n] * Real(kk) * cx * cy * sz) /
                                        normk[n];

                            vpert[2] += (-beta[n] * Real(ii) * cy * cz * sx +
                                         alpha[n] * Real(jj) * cx * cz * sy) /
                                        normk[n];
                        }
                    }
                }

                // apply the cutoff function to the perturbational velocity
                // add perturbational velocity to background velocity
                for (auto n = 0; n < 3; ++n) {
                    vpert[n] *= velpert_amplitude_loc * 0.5 *
                                (1.0 + std::tanh((velpert_r_outer -
                                                  velpert_steep_loc - rloc) /
                                                 velpert_steep_loc)) *
                                0.5 *
                                (1.0 + std::tanh((rloc - velpert_r_inner -
                                                  velpert_steep_loc) /
                                                 velpert_steep_loc));

                    vel_arr(i, j, k, n) += vpert[n];
                }
            });
        }
    }
}
