#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeBeta0(BaseState<Real>& beta0_s, const BaseState<Real>& rho0_s,
                        const BaseState<Real>& p0_s,
                        const BaseState<Real>& gamma1bar_s,
                        const BaseState<Real>& grav_cell_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeBeta0()", MakeBeta0);

    const auto& dr = base_geom.dr;

    BaseState<Real> beta0_edge_state(base_geom.finest_radial_level + 1,
                                     base_geom.nr_fine + 1);
    auto beta0_edge = beta0_edge_state.array();

    beta0_s.setVal(0.0);

    auto beta0 = beta0_s.array();
    const auto rho0 = rho0_s.const_array();
    const auto p0 = p0_s.const_array();
    const auto gamma1bar = gamma1bar_s.const_array();
    const auto grav_cell = grav_cell_s.const_array();

    if (beta0_type == 1) {
        ///////////////////////////////////////////////////////////////////////
        // Compute beta0 on the edges and average to the center
        //
        // Multilevel Outline:
        //
        // First, compute beta0 on edges and centers at level 0 only
        // Obtain the starting value from rho0 at the bottom of the domain.
        // do n=1,base_geom.finest_radial_level
        //   Compute beta0 on edges and centers at level n
        //   Obtain the starting value of beta0_edge_lo from the coarser grid
        //   if n>0, compare the difference between beta0 at the top of level n to the
        //           corresponding point on level n-1
        //   do i=n-1,0,-1
        //     Offset the centered beta on level i above this point so the total integral
        //      is consistent
        //     Redo the anelastic cutoff part
        //   }
        // }
        // call restrict_base and fill_ghost_base
        //////////////////////////////////////////////////////////////////////

        for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
            for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
                // Compute beta0 on edges and centers at level n
                if (n == 0) {
                    beta0_edge(0, 0) = rho0(0, 0);
                } else {
                    // Obtain the starting value of beta0_edge_lo from the coarser grid
                    beta0_edge(n, base_geom.r_start_coord(n, j)) =
                        beta0_edge(n - 1, base_geom.r_start_coord(n, j) / 2);
                }

                // NOTE: the integral here prevents this from being done in parallel
                for (auto r = base_geom.r_start_coord(n, j);
                     r <= base_geom.r_end_coord(n, j); ++r) {
                    Real lambda = 0.0;
                    Real mu = 0.0;
                    Real nu = 0.0;

                    if (r < base_geom.anelastic_cutoff_density_coord(n)) {
                        Real drp = use_exact_base_state
                                       ? base_geom.r_edge_loc(n, r + 1) -
                                             base_geom.r_edge_loc(n, r)
                                       : dr(n);
                        Real drm = dr(n);
                        if (use_exact_base_state) {
                            drm = r > 0 ? base_geom.r_edge_loc(n, r) -
                                              base_geom.r_edge_loc(n, r - 1)
                                        : drp;
                        }

                        if (r == 0 || r == base_geom.nr(n) - 1) {
                            // lambda = 0.0;
                            // mu = 0.0;
                            // nu = 0.0;
                        } else {
                            Real drc = use_exact_base_state
                                           ? base_geom.r_cc_loc(n, r + 1) -
                                                 base_geom.r_cc_loc(n, r - 1)
                                           : dr(n);

                            // piecewise linear reconstruction of rho0,
                            // gamma1bar, and p0 -- see paper III, appendix C
                            Real del =
                                0.5 * (rho0(n, r + 1) - rho0(n, r - 1)) / drc;
                            Real dpls =
                                2.0 * (rho0(n, r + 1) - rho0(n, r)) / drp;
                            Real dmin =
                                2.0 * (rho0(n, r) - rho0(n, r - 1)) / drm;
                            Real slim = amrex::min(amrex::Math::abs(dpls),
                                                   amrex::Math::abs(dmin));
                            slim = slim == slim ? slim : 0.0;
                            slim = dpls * dmin > 0.0 ? slim : 0.0;
                            Real sflag = amrex::Math::copysign(1.0, del);
                            lambda =
                                sflag * amrex::min(slim, amrex::Math::abs(del));

                            del = 0.5 *
                                  (gamma1bar(n, r + 1) - gamma1bar(n, r - 1)) /
                                  drc;
                            dpls = 2.0 *
                                   (gamma1bar(n, r + 1) - gamma1bar(n, r)) /
                                   drp;
                            dmin = 2.0 *
                                   (gamma1bar(n, r) - gamma1bar(n, r - 1)) /
                                   drm;
                            slim = amrex::min(amrex::Math::abs(dpls),
                                              amrex::Math::abs(dmin));
                            slim = dpls * dmin > 0.0 ? slim : 0.0;
                            sflag = amrex::Math::copysign(1.0, del);
                            mu =
                                sflag * amrex::min(slim, amrex::Math::abs(del));

                            del = 0.5 * (p0(n, r + 1) - p0(n, r - 1)) / drc;
                            dpls = 2.0 * (p0(n, r + 1) - p0(n, r)) / drp;
                            dmin = 2.0 * (p0(n, r) - p0(n, r - 1)) / drm;
                            slim = amrex::min(amrex::Math::abs(dpls),
                                              amrex::Math::abs(dmin));
                            slim = dpls * dmin > 0.0 ? slim : 0.0;
                            sflag = amrex::Math::copysign(1.0, del);
                            nu =
                                sflag * amrex::min(slim, amrex::Math::abs(del));
                        }

                        if (use_exact_base_state) {
                            // edge-to-cell-center spacings
                            drp = 2.0 * (base_geom.r_edge_loc(n, r + 1) -
                                         base_geom.r_cc_loc(n, r));
                            drm = 2.0 * (base_geom.r_cc_loc(n, r) -
                                         base_geom.r_edge_loc(n, r));
                        }

                        Real integral = 0.0;

                        if (amrex::Math::abs(nu) == 0. ||
                            amrex::Math::abs(mu) == 0. ||
                            (nu * gamma1bar(n, r) - mu * p0(n, r)) == 0.0 ||
                            ((gamma1bar(n, r) + 0.5 * mu * drp) /
                             (gamma1bar(n, r) - 0.5 * mu * drm)) <= 0.0 ||
                            ((p0(n, r) + 0.5 * nu * drp) /
                             (p0(n, r) - 0.5 * nu * drm)) <= 0.0) {
                            // just do piecewise constant integration
                            integral = amrex::Math::abs(grav_cell(n, r)) *
                                       rho0(n, r) * 0.5 * (drp + drm) /
                                       (p0(n, r) * gamma1bar(n, r));

                        } else {
                            if (use_linear_grav_in_beta0 &&
                                !use_exact_base_state) {
                                // also do piecewise linear reconstruction of
                                // gravity -- not documented in publication yet.
                                Real del = 0.5 *
                                           (grav_cell(n, r + 1) -
                                            grav_cell(n, r - 1)) /
                                           dr(n);
                                Real dpls =
                                    2.0 *
                                    (grav_cell(n, r + 1) - grav_cell(n, r)) /
                                    dr(n);
                                Real dmin =
                                    2.0 *
                                    (grav_cell(n, r) - grav_cell(n, r - 1)) /
                                    dr(n);
                                Real slim = amrex::min(amrex::Math::abs(dpls),
                                                       amrex::Math::abs(dmin));
                                slim = dpls * dmin > 0.0 ? slim : 0.0;
                                Real sflag = amrex::Math::copysign(1.0, del);
                                Real kappa =
                                    sflag *
                                    amrex::min(slim, amrex::Math::abs(del));

                                Real denom =
                                    nu * gamma1bar(n, r) - mu * p0(n, r);
                                Real coeff1 =
                                    (lambda * gamma1bar(n, r) -
                                     mu * rho0(n, r)) *
                                    (kappa * gamma1bar(n, r) +
                                     mu * amrex::Math::abs(grav_cell(n, r))) /
                                    (mu * mu * denom);
                                Real coeff2 =
                                    (lambda * p0(n, r) - nu * rho0(n, r)) *
                                    (-kappa * p0(n, r) -
                                     nu * amrex::Math::abs(grav_cell(n, r))) /
                                    (nu * nu * denom);
                                Real coeff3 = kappa * lambda / (mu * nu);

                                integral =
                                    coeff1 * log((gamma1bar(n, r) +
                                                  0.5 * mu * dr(n)) /
                                                 (gamma1bar(n, r) -
                                                  0.5 * mu * dr(n))) +
                                    coeff2 *
                                        log((p0(n, r) + 0.5 * nu * dr(n)) /
                                            (p0(n, r) - 0.5 * nu * dr(n))) -
                                    coeff3 * dr(n);

                            } else {
                                // paper III, equation C2
                                Real denom =
                                    nu * gamma1bar(n, r) - mu * p0(n, r);
                                Real coeff1 =
                                    lambda * gamma1bar(n, r) / mu - rho0(n, r);
                                Real coeff2 =
                                    lambda * p0(n, r) / nu - rho0(n, r);

                                integral =
                                    (amrex::Math::abs(grav_cell(n, r)) /
                                     denom) *
                                    (coeff1 * log((gamma1bar(n, r) +
                                                   0.5 * mu * drp) /
                                                  (gamma1bar(n, r) -
                                                   0.5 * mu * drm)) -
                                     coeff2 * log((p0(n, r) + 0.5 * nu * drp) /
                                                  (p0(n, r) - 0.5 * nu * drm)));
                            }
                        }

                        beta0_edge(n, r + 1) =
                            beta0_edge(n, r) * exp(-integral);
                        beta0(n, r) =
                            0.5 * (beta0_edge(n, r) + beta0_edge(n, r + 1));

                    } else {  // r >= anelastic_cutoff_density

                        if (amrex::Math::abs(rho0(n, r - 1)) > rel_eps) {
                            beta0(n, r) =
                                beta0(n, r - 1) * (rho0(n, r) / rho0(n, r - 1));
                        } else {
                            beta0(n, r) = beta0(n, r - 1);
                        }
                        beta0_edge(n, r + 1) =
                            2.0 * beta0(n, r) - beta0_edge(n, r);
                    }
                }

                if (n > 0) {
                    // Compare the difference between beta0 at the top of level n to the
                    // corresponding point on level n-1
                    Real offset =
                        beta0_edge(n, base_geom.r_end_coord(n, j) + 1) -
                        beta0_edge(n - 1,
                                   (base_geom.r_end_coord(n, j) + 1) / 2);

                    for (int i = n - 1; i >= 0; --i) {
                        auto refrat = (int)amrex::Math::round(pow(2, n - i));

                        // Offset the centered beta on level i above this point so the total
                        // integral is consistent
                        for (int r = base_geom.r_end_coord(n, j) / refrat + 1;
                             r <= base_geom.nr(i); ++r) {
                            beta0(i, r) += offset;
                        }

                        // Redo the anelastic cutoff part
                        for (int r =
                                 base_geom.anelastic_cutoff_density_coord(i);
                             r <= base_geom.nr(i); ++r) {
                            if (rho0(i, r - 1) != 0.0) {
                                beta0(i, r) = beta0(i, r - 1) *
                                              (rho0(i, r) / rho0(i, r - 1));
                            }
                        }

                        // This next piece of code is needed for the case when the anelastic
                        // cutoff coordinate lives on level n.  We first average beta0 from
                        // level i+1 to level i in the region between the anelastic cutoff and
                        // the top of grid n.  Then recompute beta0 at level i above the top
                        // of grid n.
                        if (base_geom.r_end_coord(n, j) >=
                            base_geom.anelastic_cutoff_density_coord(n)) {
                            for (int r =
                                     base_geom.anelastic_cutoff_density_coord(
                                         i);
                                 r <=
                                 (base_geom.r_end_coord(n, j) + 1) / refrat - 1;
                                 ++r) {
                                beta0(i, r) = 0.5 * (beta0(i + 1, 2 * r) +
                                                     beta0(i + 1, 2 * r + 1));
                            }

                            for (int r =
                                     (base_geom.r_end_coord(n, j) + 1) / refrat;
                                 r <= base_geom.nr(i); ++r) {
                                if (rho0(i, r - 1) != 0.0) {
                                    beta0(i, r) = beta0(i, r - 1) *
                                                  (rho0(i, r) / rho0(i, r - 1));
                                }
                            }
                        }
                    }  // end loop over i=n-1,0,-1
                }      // } (n  >  0)
            }          // end loop over disjoint chunks
        }              // end loop over levels

        // 0.0 the beta0 where there is no corresponding full state array
        for (int n = 1; n <= base_geom.finest_radial_level; ++n) {
            for (int j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
                if (j == base_geom.numdisjointchunks(n)) {
                    for (int r = base_geom.r_end_coord(n, j) + 1;
                         r < base_geom.nr(n); ++r) {
                        beta0(n, r) = 0.0;
                    }
                } else {
                    for (int r = base_geom.r_end_coord(n, j) + 1;
                         r < base_geom.r_start_coord(n, j + 1); ++r) {
                        beta0(n, r) = 0.0;
                    }
                }
            }
        }
    } else if (beta0_type == 2) {
        // beta_0 = rho_0
        for (int n = 0; n <= base_geom.finest_radial_level; ++n) {
            for (int j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
                // for (int r = r_start_coord(n,j); r <= r_end_coord(n,j); ++r) {
                int lo = base_geom.r_start_coord(n, j);
                int hi = base_geom.r_end_coord(n, j);
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    beta0(n, r) = rho0(n, r);
                });
                Gpu::synchronize();
            }
        }
    } else if (beta0_type == 3) {
        // beta_0 = 1.0
        for (int n = 0; n <= base_geom.finest_radial_level; ++n) {
            for (int j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
                int lo = base_geom.r_start_coord(n, j);
                int hi = base_geom.r_end_coord(n, j);
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    beta0(n, r) = 1.0;
                });
                Gpu::synchronize();
            }
        }
    }

    RestrictBase(beta0, true);
    FillGhostBase(beta0, true);
}
