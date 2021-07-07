
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeVelForce(
    Vector<MultiFab>& vel_force_cart,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& uedge_in,
    const Vector<MultiFab>& rho, const BaseState<Real>& rho0,
    const BaseState<Real>& grav_cell, const Vector<MultiFab>& w0_force_cart,
#ifdef ROTATION
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const bool is_final_update,
#endif
    int do_add_utilde_force) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVelForce()", MakeVelForce);

    Vector<MultiFab> gradw0_cart(finest_level + 1);
    Vector<MultiFab> grav_cart(finest_level + 1);
    Vector<MultiFab> rho0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gradw0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gradw0_cart[lev].setVal(0.);

        grav_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        grav_cart[lev].setVal(0.);

        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_cart[lev].setVal(0.);
    }

    BaseState<Real> gradw0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    gradw0.setVal(0.0);

    if (!(use_exact_base_state || average_base_state)) {
        auto w0_arr = w0.array();
        const Real dr0 = base_geom.dr_fine;
        auto gradw0_arr = gradw0.array();

        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE(int r) {
                gradw0_arr(l, r) = (w0_arr(l, r + 1) - w0_arr(l, r)) / dr0;
            });
            Gpu::synchronize();
        }
    }

    Put1dArrayOnCart(gradw0, gradw0_cart, false, false, bcs_u, 0, 1);
    Put1dArrayOnCart(rho0, rho0_cart, false, false, bcs_s, Rho);
    Put1dArrayOnCart(grav_cell, grav_cart, false, true, bcs_f, 0);

#ifdef ROTATION
    Put1dArrayOnCart(w0, w0_cart, 1, 1, bcs_f, 0);
#endif

    // Reset vel_force
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_force_cart[lev].setVal(0.);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get grid spacing
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(vel_force_cart[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get Array4 inputs
            const Array4<const Real> gpi_arr = gpi[lev].array(mfi);
            const Array4<const Real> rho_arr = rho[lev].array(mfi);
            const Array4<const Real> uedge = uedge_in[lev][0].array(mfi);
            const Array4<const Real> vedge = uedge_in[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wedge = uedge_in[lev][2].array(mfi);
#endif
            const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
            const Array4<const Real> gradw0_arr = gradw0_cart[lev].array(mfi);
            const Array4<const Real> w0_force = w0_force_cart[lev].array(mfi);
            const Array4<const Real> grav = grav_cart[lev].array(mfi);
            const Array4<const Real> rho0_arr = rho0_cart[lev].array(mfi);

            // output
            const Array4<Real> vel_force = vel_force_cart[lev].array(mfi);

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& domainBox = geom[lev].Domain();

            // y-direction in 2D, z-direction in 3D
            const int domhi = domainBox.hiVect()[AMREX_SPACEDIM - 1];

            // offload to GPU
            if (!spherical) {
#ifdef ROTATION
                Abort("MakeVelForce: rotation not implemented for planar");
#endif

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real rhopert = rho_arr(i, j, k) - rho0_arr(i, j, k);

                    // cutoff the buoyancy term if we are outside of the star
                    if (rho_arr(i, j, k) <
                        (buoyancy_cutoff_factor * base_cutoff_density)) {
                        rhopert = 0.0;
                    }

                    // note: if use_alt_energy_fix = T, then gphi is already
                    // weighted by beta0
                    for (int dim = 0; dim < AMREX_SPACEDIM - 1; ++dim) {
                        vel_force(i, j, k, dim) =
                            -gpi_arr(i, j, k, dim) / rho_arr(i, j, k);
                    }

                    vel_force(i, j, k, AMREX_SPACEDIM - 1) =
                        (rhopert * grav(i, j, k, AMREX_SPACEDIM - 1) -
                         gpi_arr(i, j, k, AMREX_SPACEDIM - 1)) /
                            rho_arr(i, j, k) -
                        w0_force(i, j, k, AMREX_SPACEDIM - 1);

                    if (do_add_utilde_force) {

#if (AMREX_SPACEDIM == 2)
                        if (j <= -1) {
                            // do not modify force since dw0/dr=0
                        } else if (j >= domhi) {
                            // do not modify force since dw0/dr=0
                        } else {
                            vel_force(i, j, k, 1) -=
                                (vedge(i, j + 1, k) + vedge(i, j, k)) *
                                (w0_arr(i, j + 1, k, 1) - w0_arr(i, j, k, 1)) /
                                (2.0 * dx[1]);
                        }
#else
                        if (k <= -1) {
                            // do not modify force since dw0/dr=0
                        } else if (k >= domhi) {
                            // do not modify force since dw0/dr=0
                        } else {
                            vel_force(i,j,k,2) -= (wedge(i,j,k+1)+wedge(i,j,k))*(w0_arr(i,j,k+1,2)-w0_arr(i,j,k,2)) / (2.0*dx[2]);
                        }
#endif
                    }
                });
            } else {  // spherical
#if (AMREX_SPACEDIM == 3)
                const Array4<const Real> normal_arr = normal[lev].array(mfi);

#ifdef ROTATION
                const Array4<const Real> w0macx_arr = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy_arr = w0mac[lev][1].array(mfi);
                const Array4<const Real> uold_arr = uold[lev].array(mfi);
                const auto omega_loc = omega;
#endif

                amrex::ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(
                                                int i, int j, int k) noexcept {
#ifdef ROTATION
                    const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
                    const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
                    const Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];

                    Real centrifugal_term[3];
                    if (use_centrifugal) {
                        centrifugal_term[0] = -omega_loc * omega_loc * x;
                        centrifugal_term[1] = -omega_loc * omega_loc * y;
                    } else {
                        centrifugal_term[0] = 0.0;
                        centrifugal_term[1] = 0.0;
                    }
                    centrifugal_term[2] = 0.0;
                    Real coriolis_term[3];

                    if (rho_arr(i, j, k) <
                        buoyancy_cutoff_factor * base_cutoff_density) {
                        for (auto d = 0; d < AMREX_SPACEDIM; ++d) {
                            centrifugal_term[d] = 0.0;
                        }
                    }

                    if (is_final_update) {
                        coriolis_term[0] =
                            -2.0 * omega_loc * 0.5 *
                            (vedge(i, j, k) + w0macy_arr(i, j, k) +
                             vedge(i, j + 1, k) + w0macy_arr(i, j + 1, k));
                        coriolis_term[1] =
                            2.0 * omega_loc * 0.5 *
                            (uedge(i, j, k) + w0macx_arr(i, j, k) +
                             uedge(i + 1, j, k) + w0macx_arr(i + 1, j, k));
                        coriolis_term[2] = 0.0;
                    } else {
                        coriolis_term[0] =
                            -2.0 * omega_loc *
                            (uold_arr(i, j, k, 1) + w0_arr(i, j, k, 1));
                        coriolis_term[1] =
                            2.0 * omega_loc *
                            (uold_arr(i, j, k, 0) + w0_arr(i, j, k, 0));
                        coriolis_term[2] = 0.0;
                    }
#endif
                    Real rhopert = rho_arr(i, j, k) - rho0_arr(i, j, k);

                    // cutoff the buoyancy term if we are outside of the star
                    if (rho_arr(i, j, k) <
                        (buoyancy_cutoff_factor * base_cutoff_density)) {
                        rhopert = 0.0;
                    }

                    // note: if use_alt_energy_fix = T, then gphi is already
                    // weighted by beta0
                    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                        vel_force(i, j, k, dim) =
                            (rhopert * grav(i, j, k, dim) -
                             gpi_arr(i, j, k, dim)) /
                                rho_arr(i, j, k) -
                            w0_force(i, j, k, dim);
                    }
#ifdef ROTATION
                    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                        vel_force(i, j, k, dim) +=
                            -coriolis_term[dim] - centrifugal_term[dim];
                    }
#endif

                    if (do_add_utilde_force == 1) {
                        Real Ut_dot_er =
                            0.5 * (uedge(i, j, k) + uedge(i + 1, j, k)) *
                                normal_arr(i, j, k, 0) +
                            0.5 * (vedge(i, j, k) + vedge(i, j + 1, k)) *
                                normal_arr(i, j, k, 1) +
                            0.5 * (wedge(i, j, k) + wedge(i, j, k + 1)) *
                                normal_arr(i, j, k, 2);

                        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                            vel_force(i, j, k, dim) -= Ut_dot_er *
                                                       gradw0_arr(i, j, k) *
                                                       normal_arr(i, j, k, dim);
                        }
                    }
                });
#endif
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(vel_force_cart, 0, AMREX_SPACEDIM);

    // note - we need to reconsider the bcs type here
    // it matches fortran MAESTRO but is that correct?
    FillPatch(t_old, vel_force_cart, vel_force_cart, vel_force_cart, 0, 0,
              AMREX_SPACEDIM, 0, bcs_u, 1);
}

void Maestro::ModifyScalForce(
    Vector<MultiFab>& scal_force, const Vector<MultiFab>& state,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac_in,
    const BaseState<Real>& s0_edge, const Vector<MultiFab>& s0_cart, int comp,
    const Vector<BCRec>& bcs, int fullform) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ModifyScalForce()", ModifyScalForce);

    Vector<MultiFab> s0_edge_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        s0_edge_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    if (!spherical) {
        Put1dArrayOnCart(s0_edge, s0_edge_cart, false, false, bcs_f, 0);
    }

    Vector<MultiFab> divu_cart(finest_level + 1);
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;

    if (spherical) {
        BaseState<Real> divw0(1, base_geom.nr_fine);
        divw0.setVal(0.0);
        auto divw0_arr = divw0.array();
        const auto w0_arr = w0.const_array();

        if (!use_exact_base_state) {
            for (int r = 0; r < base_geom.nr_fine - 1; ++r) {
                Real dr_loc = r_edge_loc(0, r + 1) - r_edge_loc(0, r);
                divw0_arr(0, r) =
                    (r_edge_loc(0, r + 1) * r_edge_loc(0, r + 1) *
                         w0_arr(0, r + 1) -
                     r_edge_loc(0, r) * r_edge_loc(0, r) * w0_arr(0, r)) /
                    (dr_loc * r_cc_loc(0, r) * r_cc_loc(0, r));
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            divu_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            divu_cart[lev].setVal(0.);
        }
        Put1dArrayOnCart(divw0, divu_cart, false, false, bcs_u, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get the index space and grid spacing of the domain
        const auto dx = geom[lev].CellSizeArray();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_force[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get Array4 inputs
            const Array4<const Real> scal = state[lev].array(mfi, comp);
            const Array4<const Real> umac = umac_in[lev][0].array(mfi);
            const Array4<const Real> vmac = umac_in[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wmac = umac_in[lev][2].array(mfi);
#endif
            const Array4<const Real> s0_arr = s0_cart[lev].array(mfi);
            const Array4<const Real> s0_edge_arr = s0_edge_cart[lev].array(mfi);
            const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);

            // output
            const Array4<Real> force = scal_force[lev].array(mfi, comp);

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // offload to GPU
            if (spherical) {
#if (AMREX_SPACEDIM == 3)
                // lo and hi of domain
                const Box& domainBox = geom[lev].Domain();
                const auto domlo = domainBox.loVect3d();
                const auto domhi = domainBox.hiVect3d();
                const Array4<const Real> divu_arr = divu_cart[lev].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // umac does not contain w0
                    Real divumac = (umac(i + 1, j, k) - umac(i, j, k)) / dx[0] +
                                   (vmac(i, j + 1, k) - vmac(i, j, k)) / dx[1] +
                                   (wmac(i, j, k + 1) - wmac(i, j, k)) / dx[2];

                    if (fullform == 1) {
                        force(i, j, k) -=
                            scal(i, j, k) * (divumac + divu_arr(i, j, k));
                    } else {
                        Real s0_xlo = 0.0;
                        Real s0_xhi = 0.0;
                        if (i < domhi[0]) {
                            s0_xhi =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i + 1, j, k));
                        } else {
                            s0_xhi = s0_arr(i, j, k);
                        }
                        if (i > domlo[0]) {
                            s0_xlo =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i - 1, j, k));
                        } else {
                            s0_xlo = s0_arr(i, j, k);
                        }

                        Real s0_ylo = 0.0;
                        Real s0_yhi = 0.0;
                        if (j < domhi[1]) {
                            s0_yhi =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i, j + 1, k));
                        } else {
                            s0_yhi = s0_arr(i, j, k);
                        }
                        if (j > domlo[1]) {
                            s0_ylo =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i, j - 1, k));
                        } else {
                            s0_ylo = s0_arr(i, j, k);
                        }

                        Real s0_zlo = 0.0;
                        Real s0_zhi = 0.0;
                        if (k < domhi[2]) {
                            s0_zhi =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i, j, k + 1));
                        } else {
                            s0_zhi = s0_arr(i, j, k);
                        }
                        if (k > domlo[2]) {
                            s0_zlo =
                                0.5 * (s0_arr(i, j, k) + s0_arr(i, j, k - 1));
                        } else {
                            s0_zlo = s0_arr(i, j, k);
                        }

                        Real divs0u = (umac(i + 1, j, k) * s0_xhi -
                                       umac(i, j, k) * s0_xlo) /
                                          dx[0] +
                                      (vmac(i, j + 1, k) * s0_yhi -
                                       vmac(i, j, k) * s0_ylo) /
                                          dx[1] +
                                      (wmac(i, j, k + 1) * s0_zhi -
                                       wmac(i, j, k) * s0_zlo) /
                                          dx[2];

                        force(i, j, k) +=
                            -divs0u - (scal(i, j, k) - s0_arr(i, j, k)) *
                                          (divumac + divu_arr(i, j, k));
                    }
                });
#else
                Abort("ModifyScalForce: Spherical is not valid for DIM < 3");
#endif
            } else {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // umac does not contain w0
#if (AMREX_SPACEDIM == 2)
                    Real divu = (umac(i + 1, j, k) - umac(i, j, k)) / dx[0] +
                                (vmac(i, j + 1, k) - vmac(i, j, k)) / dx[1];

                    // add w0 contribution
                    divu += (w0_arr(i, j + 1, k, AMREX_SPACEDIM - 1) -
                             w0_arr(i, j, k, AMREX_SPACEDIM - 1)) /
                            dx[AMREX_SPACEDIM - 1];
#elif (AMREX_SPACEDIM == 3)
                    Real divu = (umac(i+1,j,k) - umac(i,j,k)) / dx[0] 
                        +(vmac(i,j+1,k) - vmac(i,j,k)) / dx[1] 
                        +(wmac(i,j,k+1) - wmac(i,j,k)) / dx[2];

                    // add w0 contribution
                    divu += (w0_arr(i,j,k+1,AMREX_SPACEDIM-1)-w0_arr(i,j,k,AMREX_SPACEDIM-1))/dx[AMREX_SPACEDIM-1];
#endif
                    if (fullform == 1) {
                        force(i, j, k) -= scal(i, j, k) * divu;
                    } else {

#if (AMREX_SPACEDIM == 2)
                        Real divs0u =
                            s0_arr(i, j, k) *
                                (umac(i + 1, j, k) - umac(i, j, k)) / dx[0] +
                            (vmac(i, j + 1, k) * s0_edge_arr(i, j + 1, k) -
                             vmac(i, j, k) * s0_edge_arr(i, j, k)) /
                                dx[1];
#elif (AMREX_SPACEDIM == 3)
                        Real divs0u = s0_arr(i,j,k)*( (umac(i+1,j,k) - umac(i,j,k))/dx[0] 
                                                      +(vmac(i,j+1,k) - vmac(i,j,k))/dx[1] ) 
                            +(wmac(i,j,k+1) * s0_edge_arr(i,j,k+1) - 
                              wmac(i,j,k  ) * s0_edge_arr(i,j,k))/ dx[2];
#endif

                        force(i, j, k) +=
                            -(scal(i, j, k) - s0_arr(i, j, k)) * divu - divs0u;
                    }
                });
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(scal_force, comp, 1);

    // fill ghost cells
    FillPatch(t_old, scal_force, scal_force, scal_force, comp, comp, 1, 0,
              bcs_f);
}

void Maestro::MakeRhoHForce(
    Vector<MultiFab>& scal_force, const int is_prediction,
    const Vector<MultiFab>& thermal,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac_cart,
    const int add_thermal, const int& which_step)

{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoHForce()", MakeRhoHForce);

    // if we are doing the prediction, then it only makes sense to be in
    // this routine if the quantity we are predicting is rhoh', h, or rhoh
    if (is_prediction == 1 && !(enthalpy_pred_type == predict_rhohprime ||
                                enthalpy_pred_type == predict_h ||
                                enthalpy_pred_type == predict_rhoh)) {
        Abort(
            "ERROR: should only call mkrhohforce when predicting rhoh', h, or "
            "rhoh");
    }

    BaseState<Real> rho0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> grav(base_geom.max_radial_level + 1, base_geom.nr_fine);

    if (which_step == 1) {
        rho0.copy(rho0_old);
        p0.copy(p0_old);
    } else {
        rho0.copy(0.5 * (rho0_old + rho0_new));
        p0.copy(0.5 * (p0_old + p0_new));
    }

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> psi_cart(finest_level + 1);
    Vector<MultiFab> grav_cart(finest_level + 1);
    Vector<MultiFab> rho0_cart(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > p0mac(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        psi_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        grav_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        AMREX_D_TERM(p0mac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], 1, 1);
                     , p0mac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], 1, 1);
                     , p0mac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], 1, 1););
        psi_cart[lev].setVal(0.);
        grav_cart[lev].setVal(0.);
        rho0_cart[lev].setVal(0.);
        p0_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);
#if (AMREX_SPACEDIM == 3)
    if (spherical) {
        MakeS0mac(p0, p0mac);
    }
#endif
    Put1dArrayOnCart(psi, psi_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(rho0, rho0_cart, false, false, bcs_s, Rho);

    MakeGravCell(grav, rho0);

    Put1dArrayOnCart(grav, grav_cart, false, false, bcs_f, 0);

    // constants in Fortran
    const int enthalpy_pred_type_in = enthalpy_pred_type;
    const int predict_h_const = predict_h;
    const int predict_rhoh_const = predict_rhoh;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Get cutoff coord
        const auto base_cutoff_density_coord =
            base_geom.base_cutoff_density_coord(lev);

        // Get grid spacing
        const auto dx = geom[lev].CellSizeArray();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_force[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // inputs
            const Array4<const Real> umac = umac_cart[lev][0].array(mfi);
            const Array4<const Real> vmac = umac_cart[lev][1].array(mfi);
            const Array4<const Real> p0cart = p0_cart[lev].array(mfi);
            const Array4<const Real> p0macx = p0mac[lev][0].array(mfi);
            const Array4<const Real> p0macy = p0mac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wmac = umac_cart[lev][2].array(mfi);
            const Array4<const Real> p0macz = p0mac[lev][2].array(mfi);
#endif
            const Array4<const Real> thermal_arr = thermal[lev].array(mfi);
            const Array4<const Real> psicart = psi_cart[lev].array(mfi);
            const Array4<const Real> gravcart = grav_cart[lev].array(mfi);
            const Array4<const Real> rho0cart = rho0_cart[lev].array(mfi);

            // output
            const Array4<Real> rhoh_force = scal_force[lev].array(mfi, RhoH);

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& domainBox = geom[lev].Domain();

            const int domhi = domainBox.hiVect()[AMREX_SPACEDIM - 1];

            // if use_exact_base_state or average_base_state,
            // psi is set to dpdt in advance subroutine

            // For non-spherical, add wtilde d(p0)/dr
            // For spherical, we make u grad p = div (u p) - p div (u)
            if (!spherical) {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real gradp0 = 0.0;

#if (AMREX_SPACEDIM == 2)
                    if (j < base_cutoff_density_coord) {
                        gradp0 = rho0cart(i, j, k) * gravcart(i, j, k);
                    } else if (j == domhi) {
                        // NOTE: this should be zero since p0 is constant up here
                        gradp0 =
                            (p0cart(i, j, k) - p0cart(i, j - 1, k)) / dx[1];
                    } else {
                        // NOTE: this should be zero since p0 is constant up here
                        gradp0 =
                            (p0cart(i, j + 1, k) - p0cart(i, j, k)) / dx[1];
                    }

                    Real veladv = 0.5 * (vmac(i, j, k) + vmac(i, j + 1, k));
                    rhoh_force(i, j, k) = veladv * gradp0;
#else 
                    if (k < base_cutoff_density_coord) {
                        gradp0 = rho0cart(i,j,k) * gravcart(i,j,k);
                    } else if (k == domhi) {
                        // NOTE: this should be zero since p0 is constant up here
                        gradp0 = ( p0cart(i,j,k) - p0cart(i,j-1,k) ) / dx[2]; 
                    } else {
                        // NOTE: this should be zero since p0 is constant up here
                        gradp0 = ( p0cart(i,j+1,k) - p0cart(i,j,k) ) / dx[2];
                    }

                    Real veladv = 0.5*(wmac(i,j,k)+wmac(i,j,k+1));
                    rhoh_force(i,j,k) = veladv * gradp0;
#endif
                    // psi should always be in the force if we are doing the final update
                    // For prediction, it should not be in the force if we are predicting
                    // (rho h)', but should be there if we are predicting h or rhoh
                    //
                    // If use_exact_base_state or average_base_state is on, psi is instead dpdt term
                    if ((is_prediction == 1 &&
                         enthalpy_pred_type_in == predict_h_const) ||
                        (is_prediction == 1 &&
                         enthalpy_pred_type_in == predict_rhoh_const) ||
                        (!is_prediction)) {
                        rhoh_force(i, j, k) += psicart(i, j, k);
                    }

                    if (add_thermal == 1) {
                        rhoh_force(i, j, k) += thermal_arr(i, j, k);
                    }
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real divup = (umac(i + 1, j, k) * p0macx(i + 1, j, k) -
                                  umac(i, j, k) * p0macx(i, j, k)) /
                                     dx[0] +
                                 (vmac(i, j + 1, k) * p0macy(i, j + 1, k) -
                                  vmac(i, j, k) * p0macy(i, j, k)) /
                                     dx[1] +
                                 (wmac(i, j, k + 1) * p0macz(i, j, k + 1) -
                                  wmac(i, j, k) * p0macz(i, j, k)) /
                                     dx[2];

                    Real p0divu =
                        ((umac(i + 1, j, k) - umac(i, j, k)) / dx[0] +
                         (vmac(i, j + 1, k) - vmac(i, j, k)) / dx[1] +
                         (wmac(i, j, k + 1) - wmac(i, j, k)) / dx[2]) *
                        p0cart(i, j, k);

                    rhoh_force(i, j, k) = divup - p0divu;

                    if ((is_prediction == 1 &&
                         enthalpy_pred_type_in == predict_h_const) ||
                        (is_prediction == 1 &&
                         enthalpy_pred_type_in == predict_rhoh_const) ||
                        (is_prediction == 0)) {
                        rhoh_force(i, j, k) += psicart(i, j, k);
                    }

                    if (add_thermal == 1) {
                        rhoh_force(i, j, k) += thermal_arr(i, j, k);
                    }
                });
#endif
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(scal_force, RhoH, 1);
    FillPatch(t_old, scal_force, scal_force, scal_force, RhoH, RhoH, 1, 0,
              bcs_f);
}

void Maestro::MakeTempForce(
    Vector<MultiFab>& temp_force, const Vector<MultiFab>& scal,
    const Vector<MultiFab>& thermal,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeTempForce()", MakeTempForce);

    // if we are doing the prediction, then it only makes sense to be in
    // this routine if the quantity we are predicting is rhoh', h, or rhoh
    if (!(enthalpy_pred_type == predict_T_then_rhohprime ||
          enthalpy_pred_type == predict_T_then_h ||
          enthalpy_pred_type == predict_Tprime_then_h)) {
        Abort("ERROR: should only call mktempforce when predicting T' or T");
    }

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> psi_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        psi_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        p0_cart[lev].setVal(0.);
        psi_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(p0_old, p0_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(psi, psi_cart, false, false, bcs_f, 0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(temp_force[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Box& domainBox = geom[lev].Domain();
            const auto domhi = domainBox.hiVect3d();

            // Get grid spacing
            const auto dx = geom[lev].CellSizeArray();

            const Array4<Real> temp_force_arr = temp_force[lev].array(mfi);
            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<const Real> umac = umac_in[lev][0].array(mfi);
            const Array4<const Real> vmac = umac_in[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> wmac = umac_in[lev][2].array(mfi);
#endif
            const Array4<const Real> thermal_arr = thermal[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> psi_arr = psi_cart[lev].array(mfi);

            if (!spherical) {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                    if (j == 0) {
                        gradp0 =
                            (p0_arr(i, j + 1, k) - p0_arr(i, j, k)) / dx[1];
                    } else if (j == domhi[1]) {
                        gradp0 =
                            (p0_arr(i, j, k) - p0_arr(i, j - 1, k)) / dx[1];
                    } else {
                        gradp0 = 0.5 *
                                 (p0_arr(i, j + 1, k) - p0_arr(i, j - 1, k)) /
                                 dx[1];
                    }
#else 
                    if (k == 0) {
                        gradp0 = (p0_arr(i,j,k+1) - p0_arr(i,j,k)) / dx[2];
                    } else if (k == domhi[2]) {
                        gradp0 = (p0_arr(i,j,k) - p0_arr(i,j,k-1)) / dx[2];
                    } else {
                        gradp0 = 0.5 * (p0_arr(i,j,k+1) - p0_arr(i,j,k-1)) / dx[2];
                    }
#endif
                    eos_t eos_state;

                    eos_state.T = scal_arr(i, j, k, Temp);
                    eos_state.rho = scal_arr(i, j, k, Rho);
                    for (auto comp = 0; comp < NumSpec; ++comp) {
                        eos_state.xn[comp] =
                            scal_arr(i, j, k, FirstSpec + comp) /
                            scal_arr(i, j, k, Rho);
                    }
#if NAUX_NET > 0
                    for (auto comp = 0; comp < NumAux; ++comp) {
                        eos_state.aux[comp] =
                            scal_arr(i, j, k, FirstAux + comp) /
                            scal_arr(i, j, k, Rho);
                    }
#endif

                    // dens, temp, xmass inputs
                    eos(eos_input_rt, eos_state);

                    auto dhdp = 1.0 / scal_arr(i, j, k, Rho) +
                                (scal_arr(i, j, k, Rho) * eos_state.dedr -
                                 eos_state.p / scal_arr(i, j, k, Rho)) /
                                    (scal_arr(i, j, k, Rho) * eos_state.dpdr);

#if (AMREX_SPACEDIM == 2)
                    auto veladv = 0.5 * (vmac(i, j, k) + vmac(i, j + 1, k));
#else
                    auto veladv = 0.5*(wmac(i,j,k)+wmac(i,j,k+1));
#endif

                    temp_force_arr(i, j, k) =
                        thermal_arr(i, j, k) +
                        (1.0 - scal_arr(i, j, k, Rho) * dhdp) *
                            (veladv * gradp0 + psi_arr(i, j, k));
                    temp_force_arr(i, j, k) /=
                        (eos_state.cp * scal_arr(i, j, k, Rho));
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    eos_state.T = scal_arr(i, j, k, Temp);
                    eos_state.rho = scal_arr(i, j, k, Rho);
                    for (auto comp = 0; comp < NumSpec; ++comp) {
                        eos_state.xn[comp] =
                            scal_arr(i, j, k, FirstSpec + comp) /
                            scal_arr(i, j, k, Rho);
                    }

                    // dens, temp, xmass inputs
                    eos(eos_input_rt, eos_state);

                    auto dhdp = 1.0 / scal_arr(i, j, k, Rho) +
                                (scal_arr(i, j, k, Rho) * eos_state.dedr -
                                 eos_state.p / scal_arr(i, j, k, Rho)) /
                                    (scal_arr(i, j, k, Rho) * eos_state.dpdr);

                    auto p0_lox = 0.5 * (p0_arr(i, j, k) + p0_arr(i - 1, j, k));
                    auto p0_hix = 0.5 * (p0_arr(i, j, k) + p0_arr(i + 1, j, k));
                    auto p0_loy = 0.5 * (p0_arr(i, j, k) + p0_arr(i, j - 1, k));
                    auto p0_hiy = 0.5 * (p0_arr(i, j, k) + p0_arr(i, j + 1, k));
                    auto p0_loz = 0.5 * (p0_arr(i, j, k) + p0_arr(i, j, k - 1));
                    auto p0_hiz = 0.5 * (p0_arr(i, j, k) + p0_arr(i, j, k + 1));

                    auto divup =
                        (umac(i + 1, j, k) * p0_hix - umac(i, j, k) * p0_lox) /
                            dx[0] +
                        (vmac(i, j + 1, k) * p0_hiy - vmac(i, j, k) * p0_loy) /
                            dx[1] +
                        (wmac(i, j, k + 1) * p0_hiz - wmac(i, j, k) * p0_loz) /
                            dx[2];

                    auto p0divu =
                        ((umac(i + 1, j, k) - umac(i, j, k)) / dx[0] +
                         (vmac(i, j + 1, k) - vmac(i, j, k)) / dx[1] +
                         (wmac(i, j, k + 1) - wmac(i, j, k)) / dx[2]) *
                        p0_arr(i, j, k);

                    auto ugradp = divup - p0divu;

                    temp_force_arr(i, j, k) =
                        thermal_arr(i, j, k) +
                        (1.0 - scal_arr(i, j, k, Rho) * dhdp) *
                            (ugradp + psi_arr(i, j, k));
                    temp_force_arr(i, j, k) /=
                        (eos_state.cp * scal_arr(i, j, k, Rho));
                });
#endif
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(temp_force, Temp, 1);
    FillPatch(t_old, temp_force, temp_force, temp_force, Temp, Temp, 1, 0,
              bcs_f);
}
