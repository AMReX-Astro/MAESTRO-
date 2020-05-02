#ifndef _Maestro_F_H_
#define _Maestro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

    void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

    //////////////////////
    // in average.f90
    void average(const int* lev, const int* lo, const int* hi,
                 const amrex::Real* phi, const int* p_lo, const int* p_hi,
                 amrex::Real* phisum);

    void divide_phisum_by_ncell(amrex::Real* phisum, int* ncell, const int lev);

    void average_sphr_irreg(const int* lev, const int* lo, const int* hi,
                            const amrex::Real* phi, const int* p_lo, const int* p_hi,
                            amrex::Real* phisum, int* ncell,
                            const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);

    void divide_phisum_by_ncell_irreg(amrex::Real* phisum, int* ncell, const int lev);

    void average_sphr(amrex::Real* phisum, amrex::Real* phibar,
                      int* ncell, amrex::Real* radii, const int* finest_level);

    void sum_phi_3d_sphr(const int* lev, const int* lo, const int* hi,
                         const amrex::Real* phi, const int* p_lo, const int* p_hi,
                         amrex::Real* phisum,
                         const amrex::Real* radii, const int* finest_level,
                         const amrex::Real* dx, int* ncell,
                         const int* mask, const int* m_lo, const int* m_hi,
                         const int* use_mask);

    void compute_radii_sphr(const int* lev, amrex::Real* radii, const int* finest_level,
                            const amrex::Real* dx);
    //////////////////////

    //////////////////////
    // in base_state_geometry.f90

    void init_base_state_geometry(const int* max_radial_level_in,
                                  const int* nr_fine_in,
                                  const amrex::Real* dr_fine_in,
                                  amrex::Real* r_cc_loc,
                                  amrex::Real* r_edge_loc,
                                  const amrex::Real* dx_fine,
                                  int* nr_irreg_in);

    void init_base_state_map_sphr(const int* lo, const int* hi,
                  amrex::Real* cc_to_r, const int* c_lo, const int* c_hi,
                                  const amrex::Real* dx_fine,
                                  const amrex::Real* dx_lev);

    void compute_cutoff_coords(const amrex::Real* rho0);

    void init_multilevel(const int* tag_array,
                         const int* finest_radial_level_in);

    void destroy_base_state_geometry();
    
    void get_base_cutoff_density (const amrex::Real* base_cutoff_density_in);
    void get_base_cutoff_density_coord (const int lev, int* base_cutoff_density_coord_in);
    void get_buoyancy_cutoff_factor (const amrex::Real* buoyancy_cutoff_factor_in);
    void get_numdisjointchunks(int* nchunks);
    void get_r_start_coord(int* coord);
    void get_r_end_coord(int* coord);
    void get_finest_radial_level(int* lev);

    //////////////////////
    
    //////////////////////
    // in inlet_bc.F90

    void get_inlet_bcs (amrex::Real* params);

    //////////////////////

    //////////////////////
    // in burner.F90
    void burner_init();
    //////////////////////

    //////////////////////
    // in burner_loop.f90
#ifndef SDC
    void burner_loop(const int* lo, const int* hi,
                     const int lev,
                     const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                           amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                     const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi,
                           amrex::Real* rho_odot, const int* r_lo, const int* r_hi,
                           amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi,
                     const amrex::Real* tempbar_init_in,
                     const amrex::Real dt_in, const amrex::Real time_in, 
                     const int* mask, const int* m_lo, const int* m_hi,
                     const int use_mask);

    void burner_loop_sphr(const int* lo, const int* hi,
                          const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                                amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                          const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi,
                                amrex::Real* rho_odot, const int* r_lo, const int* r_hi,
                                amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi,
                          const amrex::Real* tempbar_init_cart, const int* t_lo, const int* t_hi,
                          const amrex::Real dt_in, const amrex::Real time_in, 
                          const int* mask, const int* m_lo, const int* m_hi,
                          const int use_mask);
    
#else
    void burner_loop(const int* lo, const int* hi,
                         const int lev,
                         const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                               amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                         const amrex::Real* source, const int* s_lo, const int* s_hi,
                         const amrex::Real* p0_in,
                         const amrex::Real dt_in,
                         const amrex::Real time_in,
                         const int* mask, const int* m_lo, const int* m_hi,
                         const int use_mask, const int sdc_iter, const int number_sdc_iterations);

    void burner_loop_sphr(const int* lo, const int* hi,
                              const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                                    amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                              const amrex::Real* source, const int* s_lo, const int* s_hi,
                              const amrex::Real* p0_cart, const int* t_lo, const int* t_hi,
                              const amrex::Real dt_in,
                              const amrex::Real time_in,
                              const int* mask, const int* m_lo, const int* m_hi,
                              const int use_mask, const int sdc_iter, const int number_sdc_iterations);
#endif
    
    void instantaneous_reaction_rates(const int* lo, const int* hi,
                                      amrex::Real* rho_omegadot, const int* o_lo, const int* o_hi,
                                      amrex::Real* rho_Hext, const int* h_lo, const int* h_hi,
                                      const amrex::Real* scal, const int* s_lo, const int* s_hi);
    //////////////////////

    //////////////////////
    // in compute_dt.F90
    void estdt(const int lev, amrex::Real* dt, amrex::Real* umax,
               const int* lo, const int* hi,
               const amrex::Real* dx,
               const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
               const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
               const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
               const amrex::Real* divu,  const int* d_lo, const int* d_hi,
               const amrex::Real* dSdt,  const int* t_lo, const int* t_hi,
               const amrex::Real* w0_cart, const int* w_lo, const int* w_hi,
               const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
               const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi);

    void estdt_sphr(amrex::Real* dt, amrex::Real* umax,
                    const int* lo, const int* hi, const amrex::Real* dx,
                    const amrex::Real* scal, const int* s_lo, const int* s_hi, const int* nc_s,
                    const amrex::Real* u, const int* u_lo, const int* u_hi, const int* nc_u,
                    const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
                    const amrex::Real* divu, const int* d_lo, const int* d_hi,
                    const amrex::Real* dSdt, const int* t_lo, const int* t_hi,
                    const amrex::Real* w0_cart, const int* w_lo, const int* w_hi,
                    const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                    const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
                    const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
                    const amrex::Real* gp0_cart, const int* g_lo, const int* g_hi);

    void firstdt(const int lev, amrex::Real* dt, amrex::Real* umax,
                 const int* lo, const int* hi,
                 const amrex::Real* dx,
                 const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
                 const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
                 const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
                 const amrex::Real* divu,  const int* d_lo, const int* d_hi,
                 const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
                 const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi);

    void firstdt_sphr(amrex::Real* dt, amrex::Real* umax,
                      const int* lo, const int* hi,
                      const amrex::Real* dx,
                      const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
                      const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
                      const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
                      const amrex::Real* divu,  const int* d_lo, const int* d_hi,
                      const amrex::Real* gp0_cart, const int* g_lo, const int* g_hi);

    void estdt_divu(amrex::Real* gp0,
                      const amrex::Real* p0, const amrex::Real* gamma1bar,
                      const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in enforce_HSE.f90
    void enforce_HSE(const amrex::Real* rho0,
                           amrex::Real* p0,
                     const amrex::Real* grav_cell,
                             const amrex::Real* r_cc_loc,
                     const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in fill_3d_data.F90
    void put_1d_array_on_cart(const int* lo, const int* hi, const int lev,
                              amrex::Real* s0_cart,
                              const int* s0_cart_lo, const int* s0_cart_hi, const int nc_s,
                              const amrex::Real* s0, const int is_input_edge_centered,
                              const int is_output_a_vector);

    void put_1d_array_on_cart_sphr(const int* lo, const int* hi,
                                   amrex::Real* s0_cart,
                                   const int* s0_cart_lo, const int* s0_cart_hi, const int nc_s,
                                   const amrex::Real* s0, const amrex::Real* dx,
                                   const int is_input_edge_centered,
                                   const int is_output_a_vector,
                                   const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
                                   const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);

    void addw0(const int* lo, const int* hi, const int lev, 
               amrex::Real* uedge, const int* u_lo, const int* u_hi,
               amrex::Real* w0, const amrex::Real mult);

    void addw0_sphr(const int* lo, const int* hi,
                    amrex::Real* umac, const int* u_lo, const int* u_hi,
                    const amrex::Real* w0mac, const int* x_lo, const int* x_hi,
                    const amrex::Real mult);

    void make_w0mac_sphr(const int* lo, const int* hi, const int idir,
                         const amrex::Real* w0,
                         amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                         amrex::Real* w0_cart, const int* w0_lo, const int* w0_hi, const int nc_w0,
                         amrex::Real* w0_nodal, const int* wn_lo, const int* wn_hi,
                         const amrex::Real* dx,
                         const amrex::Real* r_edge_loc);

    void make_w0mac_nodal(const int* lo, const int* hi,
                         const amrex::Real* w0,
                         amrex::Real* w0_nodal, const int* w0_lo, const int* w0_hi, 
                         const amrex::Real* dx);

    void make_s0mac_sphr(const int* lo, const int* hi, const int idir,
                         const amrex::Real* s0,
                         amrex::Real* s0mac, const int* x_lo, const int* x_hi,
                         amrex::Real* s0_cart, const int* s0_lo, const int* s0_hi,
                         const amrex::Real* dx,
                         const amrex::Real* r_cc_loc);

    void make_s0mac_sphr_irreg(const int* lo, const int* hi, const int idir,
                               const amrex::Real* s0,
                               amrex::Real* s0mac, const int* x_lo, const int* x_hi,
                               amrex::Real* s0_cart, const int* s0_lo, const int* s0_hi,
                               const amrex::Real* dx,
                               const amrex::Real* r_cc_loc);

    void make_normal(const int* lo, const int* hi,
                               amrex::Real* normal, const int* n_lo, const int* n_hi,
                     const amrex::Real* dx);

    void put_data_on_faces(const int* lo, const int* hi, const int idir,
                           const amrex::Real* scc, const int* cc_lo, const int* cc_hi,
                           amrex::Real* face, const int* x_lo, const int* x_hi,
                           const int harmonic_avg);
    //////////////////////

    //////////////////////
    // in initdata.f90
    void initdata(const int* lev,
                  const amrex::Real* time,
                  const int* lo, const int* hi,
                  const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
                  const amrex::Real* vel,  const int* vel_lo,  const int* vel_hi,  const int* nc_v,
                  const amrex::Real* s0_init,
                  const amrex::Real* p0_init,
                  const amrex::Real* dx);

    void initdata_sphr(const amrex::Real* time,
                       const int* lo, const int* hi,
                       const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
                       const amrex::Real* vel,  const int* vel_lo,  const int* vel_hi,  const int* nc_v,
                       const amrex::Real* s0_init,
                       const amrex::Real* p0_init,
                       const amrex::Real* dx,
                       const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
                       const amrex::Real* cell_cc_to_r, const int* ccr_lo, const int* ccr_hi);
    //////////////////////

    //////////////////////
    // in maestro_init.f90

    void maestro_network_init();
    void maestro_eos_init();
    void maestro_extern_init();
    void maestro_conductivity_init();
    void get_num_spec(int* nspec);
    void get_spec_names(int* spec_names, int* ispec, int* len);
    void get_spec_az(const int* ispec, amrex::Real* A, amrex::Real* Z);

    // set fortran-specific parameters in meth_params.F90
    void set_method_params (const amrex::Real* prob_lo_in,
                            const amrex::Real* prob_hi_in);

    void set_rel_eps (const amrex::Real* rel_eps_in);
    void get_rel_eps (const amrex::Real* rel_eps_in);

    void init_base_state (const amrex::Real* s0_init,
                          const amrex::Real* p0_init,
                          const amrex::Real* rho0,
                          const amrex::Real* rhoh0,
                          const amrex::Real* p0,
                          const amrex::Real* tempbar,
                          const amrex::Real* tempbar_init,
                          const int lev);

    void init_base_state_irreg (const amrex::Real* s0_init,
                                const amrex::Real* p0_init,
                                const amrex::Real* rho0,
                                const amrex::Real* rhoh0,
                                const amrex::Real* p0,
                                const amrex::Real* tempbar,
                                const amrex::Real* tempbar_init,
                                const amrex::Real* r_cc_loc,
                                const amrex::Real* r_edge_loc,
                                const int lev);
    //////////////////////

    //////////////////////
    // in make_intra_coeffs.F90
    void make_intra_coeffs(const int* lo,  const int* hi,
                           const amrex::Real* scalold, const int* s1_lo, const int* s1_hi,
                           const amrex::Real* scalnew, const int* s2_lo, const int* s2_hi,
                           amrex::Real* cp, const int* cp_lo, const int* cp_hi,
                           amrex::Real* xi, const int* xi_lo, const int* xi_hi);
    //////////////////////

    //////////////////////
    // in make_explicit_thermal.F90
    void make_thermal_coeffs(const int* lo,  const int* hi,
                             const amrex::Real* scal, const int* s_lo, const int* s_hi,
                             amrex::Real* Tcoeff, const int* t_lo, const int* t_hi,
                             amrex::Real* hcoeff, const int* h_lo, const int* h_hi,
                             amrex::Real* Xkcoeff, const int* xk_lo, const int* xk_hi,
                             amrex::Real* pcoeff, const int* p_lo, const int* p_hi);
    //////////////////////

    //////////////////////
    // in make_grav.f90
    void make_grav_cell(      amrex::Real* grav_cell,
                        const amrex::Real* rho0,
                        const amrex::Real* r_cc_loc,
                        const amrex::Real* r_edge_loc);

    void make_grav_edge(      amrex::Real* grav_edge,
                        const amrex::Real* rho0,
                        const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in make_heating.f90
    void make_heating(const int* lo, const int* hi,
                            amrex::Real* rho_Hext, const int* r_lo, const int* r_hi,
                      const amrex::Real* scal,     const int* s_lo, const int* s_hi, const int* nc_s,
                      const amrex::Real* dx, const amrex::Real* time);
    //////////////////////

    //////////////////////
    // in make_plot_variables.F90
    void make_magvel(const int* lo, const int* hi,
                     const amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const amrex::Real* w0, const int* w_lo, const int* w_hi, const int nc_w0,
                     amrex::Real* magvel, const int* m_lo, const int* m_hi);

    void make_magvel_sphr(const int* lo, const int* hi,
                          const amrex::Real* vel, const int* v_lo, const int* v_hi,
                          amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                          amrex::Real* w0macy, const int* y_lo, const int* y_hi,
                          amrex::Real* w0macz, const int* z_lo, const int* z_hi,
                          amrex::Real* magvel, const int* m_lo, const int* m_hi);

    void make_velrc(const int* lo, const int* hi,
                        const amrex::Real* vel, const int* v_lo, const int* v_hi,
                        const amrex::Real* w0rcart, const int* w_lo, const int* w_hi,
                        const amrex::Real* normal, const int* n_lo, const int* n_hi,
                        amrex::Real* rad_vel, const int* r_lo, const int* r_hi,
                        amrex::Real* circ_vel, const int* c_lo, const int* c_hi);

    void make_ad_excess(const int* lo, const int* hi,
                        const amrex::Real* state, const int* s_lo, const int* s_hi,
                        const int nc_s,
                        amrex::Real* ad_excess, const int* a_lo, const int* a_hi);

    void make_ad_excess_sphr(const int* lo, const int* hi,
                        const amrex::Real* state, const int* s_lo, const int* s_hi,
                        const int nc_s,
                        const amrex::Real* normal, const int* n_lo, const int* n_hi,
                        amrex::Real* ad_excess, const int* a_lo, const int* a_hi);

    void make_vorticity(const int* lo, const int* hi,
                        const int* domlo, const int* domhi,
                        const amrex::Real* vel, const int* v_lo, const int* v_hi,
                        const amrex::Real* dx,
                        amrex::Real* vorticity, const int* d_lo, const int* d_hi,
                        const int* bc);

    void make_deltagamma(const int* lo, const int* hi,
                         const amrex::Real* state, const int* s_lo, const int* s_hi,
                         const int nc_s,
                         const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
                         const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi,
                         amrex::Real* deltagamma, const int* d_lo, const int* d_hi);

    void make_entropy(const int* lo, const int* hi, const int lev,
                     const amrex::Real* state, const int* s_lo, const int* s_hi,
                     const int nc_s,
                     amrex::Real* entropy, const int* e_lo, const int* e_hi);

    void make_divw0(const int* lo, const int* hi, 
                    const amrex::Real* w0, const int* w_lo, const int* w_hi, const int nc_w0,
                    const amrex::Real* dx,
                    amrex::Real* divw0, const int* d_lo, const int* d_hi);

    void make_divw0_sphr(const int* lo, const int* hi,
                         const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                         const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
                         const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
                         const amrex::Real* dx,
                         amrex::Real* divw0, const int* d_lo, const int* d_hi);

     void make_pidivu(const int* lo, const int* hi,
                     const amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const amrex::Real* dx,
                     const amrex::Real* pi_cc, const int* p_lo, const int* p_hi, const int nc,
                     amrex::Real* pidivuw0, const int* d_lo, const int* d_hi);

    void make_abar(const int* lo, const int* hi,
                   const amrex::Real* state, const int* s_lo, const int* s_hi, const int nc_s,
                   amrex::Real* abar, const int* a_lo, const int* a_hi);
    //////////////////////

    //////////////////////
    // in meth_params.F90

    // read in shared C++/Fortran parameters from inputs file
    void read_method_params();

    // deallocate strings
    void finalize_meth_params();

    //////////////////////

    //////////////////////
    // in sponge.F90
    void init_sponge(const amrex::Real* rho0);

    void init_sponge_irreg(const amrex::Real* rho0,
                           const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc);

    void mk_sponge (const int* lo, const int* hi,
                    amrex::Real* sponge,  const int* s_lo, const int* s_hi,
                    const amrex::Real* dx, const amrex::Real dt);
    //////////////////////


#ifdef __cplusplus
}
#endif

#endif
