#include <Maestro.H>
#include <Maestro_F.H>
#if NAUX_NET > 0
#include <actual_network.H>
#endif

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0,
                            BaseState<Real>& p0, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (use_exact_base_state && !spherical) {
        Abort("Irregular base state not valid for planar");
    }

    if (!input_model.model_initialized) {
        // read model file
        input_model.ReadFile(model_file);
    }

    const int npts_model = input_model.npts_model;
    const Real TINY = 1.e-10;
    const int n = lev;
    const auto& dr = base_geom.dr;

    Real base_cutoff_density_loc = 1.e99;
    Real model_dr =
        (input_model.model_r[npts_model - 1] - input_model.model_r[0]) /
        Real(npts_model - 1);
    Real rmax = input_model.model_r[npts_model - 1];

    auto rhoh0_arr = rhoh0.array();
    auto rho0_arr = rho0.array();
    auto p0_arr = p0.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    // irregular base state dr
    BaseState<Real> model_dr_irreg(npts_model);
    auto model_dr_irreg_arr = model_dr_irreg.array();
    if (use_exact_base_state) {
        model_dr_irreg_arr(0) = input_model.model_r[0];
        for (auto i = 1; i < npts_model; ++i) {
            model_dr_irreg_arr(i) =
                input_model.model_r[i] - input_model.model_r[i - 1];
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        if (!spherical) {
            log_file.Log("model file mapping, level: ", n);
        } else {
            log_file.Log("model file mapping (spherical base state)");
        }

        log_file.Log("dr of MAESTRO base state =                            ",
                     dr(n));
        log_file.Log("dr of input file data =                               ",
                     model_dr);
        log_file.Log(" ");
        log_file.Log("maximum radius (cell-centered) of input model =       ",
                     rmax);

        Real mod_dr = 0.0;
        if (use_exact_base_state) {
            mod_dr = base_geom.r_cc_loc(lev, 0) < model_dr_irreg_arr(0)
                         ? remainder(model_dr_irreg_arr(0),
                                     base_geom.r_cc_loc(lev, 0))
                         : remainder(base_geom.r_cc_loc(lev, 0),
                                     model_dr_irreg_arr(0));
        } else {
            mod_dr = dr(n) < model_dr ? remainder(model_dr, dr(n))
                                      : remainder(dr(n), model_dr);
        }

        if (mod_dr > TINY) {
            log_file.Log(" ");
            log_file.Log(
                "WARNING: resolution of base state array is not an integer");
            log_file.Log(
                "         multiple of the initial model's resolution.     ");
            log_file.Log(
                "         make sure this is a desired property as this    ");
            log_file.Log(
                "         could lead to aliasing when performing the      ");
            log_file.Log(
                "         interpolation.                                  ");
            log_file.Log(" ");
            log_file.Log("modulus = ", mod_dr);
        }
    }

    Real starting_rad = spherical ? 0.0 : geom[0].ProbLo(AMREX_SPACEDIM - 1);

    Real rho_above_cutoff = s0_init_arr(n, 0, Rho);
    Real rhoh_above_cutoff = s0_init_arr(n, 0, RhoH);
    RealVector spec_above_cutoff(NumSpec);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        spec_above_cutoff[comp] = s0_init_arr(n, 0, FirstSpec + comp);
    }
#if NAUX_NET > 0
    RealVector aux_above_cutoff(NumAux);
    for (auto comp = 0; comp < NumAux; ++comp) {
        aux_above_cutoff[comp] = s0_init_arr(n, 0, FirstAux + comp);
    }
    // initialize bion in network
    network_init();
#endif
    Real temp_above_cutoff = s0_init_arr(n, 0, Temp);
    Real p_above_cutoff = p0_init_arr(n, 0);

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        Real rloc = use_exact_base_state
                        ? base_geom.r_cc_loc(n, r)
                        : starting_rad + (Real(r) + 0.5) * dr(n);

        // here we account for r > rmax of the model.hse array, assuming
        // that the state stays constant beyond rmax
        rloc = amrex::min(rloc, rmax);

        // also, if we've fallen below the cutoff density, just keep the
        // model constant
        if (rloc > base_cutoff_density_loc) {
            s0_init_arr(n, r, Rho) = rho_above_cutoff;
            s0_init_arr(n, r, RhoH) = rhoh_above_cutoff;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init_arr(n, r, FirstSpec + comp) = spec_above_cutoff[comp];
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                s0_init_arr(n, r, FirstAux + comp) = aux_above_cutoff[comp];
            }
#endif
            p0_init_arr(n, r) = p_above_cutoff;
            s0_init_arr(n, r, Temp) = temp_above_cutoff;

        } else {
            Real d_ambient =
                input_model.Interpolate(rloc, ModelParser::idens_model);
            Real t_ambient =
                input_model.Interpolate(rloc, ModelParser::itemp_model);
            Real p_ambient =
                input_model.Interpolate(rloc, ModelParser::ipres_model);

            RealVector xn_ambient(NumSpec);

            Real sumX = 0.0;

            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = amrex::max(
                    0.0, amrex::min(
                             1.0, input_model.Interpolate(
                                      rloc, ModelParser::ispec_model + comp)));
                sumX += xn_ambient[comp];
            }

            for (auto comp = 0; comp < NumSpec; ++comp) {
                if (sumX != 0.0) {
                    xn_ambient[comp] /= sumX;
                }
            }
#if NAUX_NET > 0
            RealVector aux_ambient(NumAux, 0.0);
#endif

            // initialize the aux variables
#ifdef NSE_THERMO
            for (auto comp = 0; comp < NumSpec; ++comp) {
                // set the aux quantities
                aux_ambient[iye] +=
                    xn_ambient[comp] * zion[comp] * aion_inv[comp];
                aux_ambient[iabar] += xn_ambient[comp] * aion_inv[comp];
                aux_ambient[ibea] +=
                    xn_ambient[comp] * aprox19::bion(comp) * aion_inv[comp];
            }

            aux_ambient[iabar] = 1.0_rt / aux_ambient[iabar];
#endif

            eos_t eos_state;

            // use the EOS to make the state consistent
            eos_state.T = t_ambient;
            eos_state.rho = d_ambient;
            eos_state.p = p_ambient;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = xn_ambient[comp];
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                eos_state.aux[comp] = aux_ambient[comp];
            }
#endif

            // (rho,T) --> p,h
            eos(eos_input_rt, eos_state);

            s0_init_arr(n, r, Rho) = d_ambient;
            s0_init_arr(n, r, RhoH) = d_ambient * eos_state.h;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init_arr(n, r, FirstSpec + comp) =
                    d_ambient * xn_ambient[comp];
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                s0_init_arr(n, r, FirstAux + comp) =
                    d_ambient * aux_ambient[comp];
            }
#endif
            p0_init_arr(n, r) = eos_state.p;  // p_ambient !
            s0_init_arr(n, r, Temp) = t_ambient;

            // keep track of the height where we drop below the cutoff density
            if (s0_init_arr(n, r, Rho) <= base_cutoff_density &&
                base_cutoff_density_loc == 1.e99) {
                Print() << ' ' << std::endl;
                Print() << "setting r_cutoff to " << r << std::endl;
                Print() << "radius at r_cutoff " << rloc << std::endl;

                base_cutoff_density_loc = rloc;

                rho_above_cutoff = s0_init_arr(n, r, Rho);
                rhoh_above_cutoff = s0_init_arr(n, r, RhoH);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    spec_above_cutoff[comp] =
                        s0_init_arr(n, r, FirstSpec + comp);
                }
#if NAUX_NET > 0
                for (auto comp = 0; comp < NumAux; ++comp) {
                    aux_above_cutoff[comp] = s0_init_arr(n, r, FirstAux + comp);
                }
#endif
                temp_above_cutoff = s0_init_arr(n, r, Temp);
                p_above_cutoff = p0_init_arr(n, r);
            }
        }
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        rho0_arr(lev, r) = s0_init_arr(lev, r, Rho);
        rhoh0_arr(lev, r) = s0_init_arr(lev, r, RhoH);
        tempbar_arr(lev, r) = s0_init_arr(lev, r, Temp);
        tempbar_init_arr(lev, r) = s0_init_arr(lev, r, Temp);
        p0_arr(lev, r) = p0_init_arr(lev, r);
    }

    // check whether we are in HSE

    Real mencl = 0.0;

    if (use_exact_base_state) {
        Real dr_irreg = base_geom.r_edge_loc(n, 1) -
                        base_geom.r_edge_loc(n, 0);  // edge-to-edge

        if (spherical || do_2d_planar_octant) {
            mencl = 4.0 / 3.0 * M_PI * dr_irreg * dr_irreg * dr_irreg *
                    s0_init_arr(n, 0, Rho);
        }
    } else {
        if (spherical || do_2d_planar_octant) {
            mencl = 4.0 / 3.0 * M_PI * dr(n) * dr(n) * dr(n) *
                    s0_init_arr(n, 0, Rho);
        }
    }

    Real max_hse_error = -1.e30;

    for (auto r = 1; r < base_geom.nr(n); ++r) {
        Real rloc = use_exact_base_state
                        ? base_geom.r_cc_loc(n, r)
                        : starting_rad + (Real(r) + 0.5) * dr(n);
        rloc = amrex::min(rloc, rmax);

        if (rloc < base_cutoff_density_loc) {
            Real r_r = starting_rad;
            r_r += use_exact_base_state ? base_geom.r_edge_loc(n, r + 1)
                                        : Real(r + 1) * dr(n);
            Real r_l = starting_rad;
            r_l += use_exact_base_state ? base_geom.r_edge_loc(n, r)
                                        : Real(r) * dr(n);

            Real dr_local = r_r - r_l;

            Real g = 0.0;

            if (spherical || do_2d_planar_octant) {
                g = -Gconst * mencl / (r_l * r_l);
                mencl += 4.0 / 3.0 * M_PI * dr_local *
                         (r_l * r_l + r_l * r_r + r_r * r_r) *
                         s0_init_arr(n, r, Rho);
            } else {
                if (!do_planar_invsq_grav) {
                    g = grav_const;
                } else {
                    g = -Gconst * planar_invsq_mass / (r_l * r_l);
                }
            }

            Real dpdr = 0.0;
            Real rhog = 0.0;
            if (use_exact_base_state) {
                Real dr_irreg =
                    base_geom.r_cc_loc(n, r) - base_geom.r_cc_loc(n, r - 1);
                dpdr = (p0_init_arr(n, r) - p0_init_arr(n, r - 1)) / dr_irreg;

                Real rfrac = (base_geom.r_edge_loc(n, r) -
                              base_geom.r_cc_loc(n, r - 1)) /
                             dr_irreg;
                rhog = ((1.0 - rfrac) * s0_init_arr(n, r, Rho) +
                        rfrac * s0_init_arr(n, r - 1, Rho)) *
                       g;
            } else {
                dpdr = (p0_init_arr(n, r) - p0_init_arr(n, r - 1)) / dr(n);
                rhog = 0.5 *
                       (s0_init_arr(n, r, Rho) + s0_init_arr(n, r - 1, Rho)) *
                       g;
            }

            if (print_init_hse_diag) {
                Print() << "r, dpdr, rhog, err: " << rloc << ", " << dpdr
                        << ", " << rhog << ", "
                        << amrex::Math::abs(dpdr - rhog) /
                               amrex::Math::abs(rhog)
                        << std::endl;
            }

            max_hse_error =
                amrex::max(max_hse_error, amrex::Math::abs(dpdr - rhog) /
                                              amrex::Math::abs(rhog));
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        log_file.Log(" ");
        log_file.Log("Maximum HSE Error = ", max_hse_error);
        log_file.Log(
            "   (after putting initial model into base state arrays, and");
        log_file.Log("    for density < base_cutoff_density)");
        log_file.Log(" ");
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
