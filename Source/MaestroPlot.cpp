#include <Maestro.H>
#include <Maestro_F.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>
#include <iterator>     // std::istream_iterator

using namespace amrex;

// write a small plotfile to disk
void Maestro::WriteSmallPlotFile (const int step,
                                  const Real t_in,
                                  const Real dt_in,
                                  const RealVector& rho0_in,
                                  const RealVector& rhoh0_in,
                                  const RealVector& p0_in,
                                  const RealVector& gamma1bar_in,
                                  const Vector<MultiFab>& u_in,
                                  Vector<MultiFab>& s_in,
                                  const Vector<MultiFab>& S_cc_in)
{
    WritePlotFile (step, t_in, dt_in, rho0_in, rhoh0_in, p0_in,
                   gamma1bar_in, u_in, s_in, S_cc_in, true);
}

// write plotfile to disk
void
Maestro::WritePlotFile (const int step,
                        const Real t_in,
                        const Real dt_in,
                        const RealVector& rho0_in,
                        const RealVector& rhoh0_in,
                        const RealVector& p0_in,
                        const RealVector& gamma1bar_in,
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& s_in,
                        const Vector<MultiFab>& S_cc_in,
                        const bool is_small)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WritePlotFile()",WritePlotFile);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    std::string plotfilename;

    if (!is_small) {
        plotfilename = plot_base_name;
    } else {
        plotfilename = small_plot_base_name;
    }

    if (step == plotInitData) {
        if (plotfilename.back() == '_') {
            plotfilename += "InitData";
        } else {
            plotfilename += +"_InitData";
        }

    }
    else if (step == plotInitProj) {
        if (plotfilename.back() == '_') {
            plotfilename += "after_InitProj";
        } else {
            plotfilename += +"_after_InitProj";
        }
    }
    else if (step == plotDivuIter) {
        if (plotfilename.back() == '_') {
            plotfilename += "after_DivuIter";
        } else {
            plotfilename += +"_after_DivuIter";
        }
    }
    else {
        PlotFileName(step, &plotfilename);
    }

    // convert rho0 to multi-D MultiFab
    Vector<MultiFab> rho0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(rho0_in,rho0_cart,0,0);

    // convert rhoh0 to multi-D MultiFab
    Vector<MultiFab> rhoh0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rhoh0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(rhoh0_in,rhoh0_cart,0,0);

    // convert p0 to multi-D MultiFab
    Vector<MultiFab> p0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(p0_in,p0_cart,0,0);

    // convert gamma1bar to multi-D MultiFab
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    Put1dArrayOnCart(gamma1bar_in,gamma1bar_cart,0,0);

    int nPlot = 0;
    const auto& varnames = PlotFileVarNames(&nPlot);

    const auto& mf = PlotFileMF(nPlot,t_in,dt_in,rho0_cart,rhoh0_cart,p0_cart,
                                gamma1bar_cart,u_in,s_in,p0_in,gamma1bar_in,
                                S_cc_in);

    // WriteMultiLevelPlotfile expects an array of step numbers
    Vector<int> step_array;
    step_array.resize(maxLevel()+1, step);

    if (!is_small) {
        WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                Geom(), t_in, step_array, refRatio());
    } else {
        int nSmallPlot = 0;
        const auto& small_plot_varnames = SmallPlotFileVarNames(&nSmallPlot,
                                                                varnames);

        const auto& small_mf = SmallPlotFileMF(nPlot, nSmallPlot, mf, varnames,
                                               small_plot_varnames);

        WriteMultiLevelPlotfile(plotfilename, finest_level+1, small_mf,
                                small_plot_varnames, Geom(), t_in, step_array,
                                refRatio());

        for (int i = 0; i <= finest_level; ++i)
            delete small_mf[i];
    }

    WriteJobInfo(plotfilename);

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    // write out the cell-centered base state
    if (ParallelDescriptor::IOProcessor()) {

        for (int lev=0; lev<=base_geom.max_radial_level; ++lev) {

            std::ofstream BaseCCFile;
            BaseCCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            std::string BaseCCFileName(plotfilename + "/BaseCC_");
            std::string levStr = std::to_string(lev);
            BaseCCFileName.append(levStr);
            BaseCCFile.open(BaseCCFileName.c_str(), std::ofstream::out   |
                            std::ofstream::trunc |
                            std::ofstream::binary);
            if(!BaseCCFile.good()) {
                amrex::FileOpenFailed(BaseCCFileName);
            }

            BaseCCFile.precision(17);

            BaseCCFile << "r_cc  rho0  rhoh0  p0  gamma1bar \n";

            for (int i=0; i<base_geom.nr(lev); ++i) {
                BaseCCFile << base_geom.r_cc_loc(lev,i) << " "
                           << rho0_in[lev+(base_geom.max_radial_level+1)*i] << " "
                           << rhoh0_in[lev+(base_geom.max_radial_level+1)*i] << " "
                           << p0_in[lev+(base_geom.max_radial_level+1)*i] << " "
                           << gamma1bar_in[lev+(base_geom.max_radial_level+1)*i] << "\n";
            }
        }
    }

    // write out the face-centered base state
    if (ParallelDescriptor::IOProcessor()) {

        for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {

            std::ofstream BaseFCFile;
            BaseFCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            std::string BaseFCFileName(plotfilename + "/BaseFC_");
            std::string levStr = std::to_string(lev);
            BaseFCFileName.append(levStr);
            BaseFCFile.open(BaseFCFileName.c_str(), std::ofstream::out   |
                            std::ofstream::trunc |
                            std::ofstream::binary);
            if(!BaseFCFile.good()) {
                amrex::FileOpenFailed(BaseFCFileName);
            }

            BaseFCFile.precision(17);

            BaseFCFile << "r_edge  w0 \n";

            for (int i = 0; i <= base_geom.nr(lev); ++i) {
                BaseFCFile << base_geom.r_edge_loc(lev,i) << " "
                           << w0[lev+(base_geom.max_radial_level+1)*i] << "\n";
            }
        }
    }

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;

    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to write plotfile: " << end_total << '\n';
    }

    for (int i = 0; i <= finest_level; ++i) {
        delete mf[i];
    }

}


// get plotfile name
void
Maestro::PlotFileName (const int lev, std::string* plotfilename)
{
    *plotfilename = Concatenate(*plotfilename, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const int nPlot,
                     const Real t_in,
                     const Real dt_in,
                     const Vector<MultiFab>& rho0_cart,
                     const Vector<MultiFab>& rhoh0_cart,
                     const Vector<MultiFab>& p0_cart,
                     const Vector<MultiFab>& gamma1bar_cart,
                     const Vector<MultiFab>& u_in,
                     Vector<MultiFab>& s_in,
                     const RealVector& p0_in,
                     const RealVector& gamma1bar_in,
                     const Vector<MultiFab>& S_cc_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;

    // temporary MultiFab to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level+1);

    // temporary MultiFab for calculations
    Vector<MultiFab> tempmf(finest_level+1);
    Vector<MultiFab> tempmf_scalar1(finest_level+1);
    Vector<MultiFab> tempmf_scalar2(finest_level+1);
    RealVector tempbar_plot ((base_geom.max_radial_level+1)*base_geom.nr_fine);
    tempbar_plot.shrink_to_fit();
    std::fill(tempbar_plot.begin(), tempbar_plot.end(), 0.);

    int dest_comp = 0;

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot,0);
        tempmf[i].define(grids[i],dmap[i],AMREX_SPACEDIM,0);

        tempmf_scalar1[i].define(grids[i],dmap[i],1,0);
        tempmf_scalar2[i].define(grids[i],dmap[i],1,0);
    }

    // velocity
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(u_in[i],0,dest_comp,AMREX_SPACEDIM);
    }
    dest_comp += AMREX_SPACEDIM;

    // magvel
    MakeMagvel(u_in, tempmf);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    ++dest_comp;

    // momentum = magvel * rho
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        MultiFab::Multiply(*plot_mf_data[i], s_in[i], Rho, dest_comp, 1, 0);
    }
    ++dest_comp;

    // vorticity
    MakeVorticity(u_in, tempmf);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    ++dest_comp;

    // rho
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Rho,dest_comp,1);
    }
    ++dest_comp;

    // rhoh
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
    }
    ++dest_comp;

    // h
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
        MultiFab::Divide(*plot_mf_data[i], s_in[i], Rho, dest_comp, 1, 0);
    }
    ++dest_comp;

    // rhoX
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],FirstSpec,dest_comp,NumSpec);
    }
    dest_comp += NumSpec;

    if (plot_spec) {
        // X
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(s_in[i],FirstSpec,dest_comp,NumSpec);
            for (int comp=0; comp<NumSpec; ++comp) {
                MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
            }
        }
        dest_comp += NumSpec;

        // abar
        MakeAbar(s_in, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;

    }

    Vector<MultiFab> stemp             (finest_level+1);
    Vector<MultiFab> rho_Hext          (finest_level+1);
    Vector<MultiFab> rho_omegadot      (finest_level+1);
    Vector<MultiFab> rho_Hnuc          (finest_level+1);
    Vector<MultiFab> sdc_source        (finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        stemp             [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        rho_Hext          [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_omegadot      [lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc          [lev].define(grids[lev], dmap[lev],       1, 0);
        sdc_source        [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        
        sdc_source[lev].setVal(0.);
    }

#ifndef SDC
    if (dt_in < small_dt) {
        React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, small_dt, t_in);
    } else {
        React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, dt_in*0.5, t_in);
    }
#else   
    if (dt_in < small_dt) {
        ReactSDC(s_in, stemp, rho_Hext, p0_in, small_dt, t_in, sdc_source, 0, 1);
    } else {
        ReactSDC(s_in, stemp, rho_Hext, p0_in, dt_in*0.5, t_in, sdc_source, 0, 1);
    }
    
    MakeReactionRates(rho_omegadot,rho_Hnuc,s_in);
#endif

    if (plot_spec || plot_omegadot) {
        // omegadot
        if (plot_omegadot) {
            for (int i = 0; i <= finest_level; ++i) {
                plot_mf_data[i]->copy(rho_omegadot[i],0,dest_comp,NumSpec);
                for (int comp=0; comp<NumSpec; ++comp) {
                    MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
                }
            }
            dest_comp += NumSpec;
        }
    }

    if (plot_Hext) {
        // Hext
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(rho_Hext[i],0,dest_comp,1);
            MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp,1,0);
        }
        ++dest_comp;
    }

    if (plot_Hnuc) {
        // Hnuc
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(rho_Hnuc[i],0,dest_comp,1);
            MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp,1,0);
        }
        ++dest_comp;
    }

    if (plot_eta) {
        // eta_rho
        Put1dArrayOnCart(etarho_cc,tempmf,1,0,bcs_u,0,1);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;
    }

    // compute tfromp
    TfromRhoP(s_in,p0_in);
    // tfromp
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
    }
    ++dest_comp;

    // compute tfromh
    TfromRhoH(s_in,p0_in);
    for (int i = 0; i <= finest_level; ++i) {
        // tfromh
        plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
    }
    ++dest_comp;

    // deltap
    // compute & copy tfromp
    PfromRhoH(s_in,s_in,tempmf);
    for (int i = 0; i <= finest_level; ++i) {
        // tfromh
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        MultiFab::Subtract(*plot_mf_data[i],p0_cart[i],0,dest_comp,1,0);
    }
    ++dest_comp;

    // deltaT
    // compute & copy tfromp
    TfromRhoP(s_in,p0_in);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
    }
    // compute tfromh
    TfromRhoH(s_in,p0_in);
    // compute deltaT = (tfromp - tfromh) / tfromh
    for (int i = 0; i <= finest_level; ++i) {
        MultiFab::Subtract(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
        MultiFab::Divide(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
    }
    ++dest_comp;

    // restore tfromp if necessary
    if (use_tfromp) {
        TfromRhoP(s_in,p0_in);
    }

    // pi
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
    }
    ++dest_comp;

    // pioverp0
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
        MultiFab::Divide(*plot_mf_data[i], p0_cart[i], 0, dest_comp, 1, 0);
    }
    ++dest_comp;

    // p0pluspi
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
        MultiFab::Add(*plot_mf_data[i], p0_cart[i], 0, dest_comp, 1, 0);
    }
    ++dest_comp;

    if (plot_gpi) {
        // gpi
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(gpi[i],0,dest_comp,AMREX_SPACEDIM);
        }
        dest_comp += AMREX_SPACEDIM;
    }

    // rhopert
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],Rho,dest_comp,1);
        MultiFab::Subtract(*plot_mf_data[i],rho0_cart[i],0,dest_comp,1,0);
    }
    ++dest_comp;

    // rhohpert
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
        MultiFab::Subtract(*plot_mf_data[i],rhoh0_cart[i],0,dest_comp,1,0);
    }
    ++dest_comp;

    // tpert
    {
        Average(s_in, tempbar_plot, Temp);
        Put1dArrayOnCart(tempbar_plot,tempmf,0,0,bcs_f,0);

        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
            MultiFab::Subtract(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
        }
    }
    ++dest_comp;

    if (plot_base_state) {
        // rho0, rhoh0, h0 and p0
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy( rho0_cart[i],0,dest_comp,1);
            plot_mf_data[i]->copy(rhoh0_cart[i],0,dest_comp+1,1);
            plot_mf_data[i]->copy(rhoh0_cart[i],0,dest_comp+2,1);

            // we have to use protected_divide here to guard against division by zero
            // in the case that there are zeros rho0
            MultiFab& plot_mf_data_mf = *plot_mf_data[i];
            for ( MFIter mfi(plot_mf_data_mf); mfi.isValid(); ++mfi ) {
                plot_mf_data_mf[mfi].protected_divide<RunOn::Device>(plot_mf_data_mf[mfi], dest_comp, dest_comp+2);
            }

            plot_mf_data[i]->copy(p0_cart[i],0,dest_comp+3,1);
        }
        dest_comp += 4;
    }

    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
    Vector<MultiFab> w0r_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        if (spherical == 1) {
            // w0mac will contain an edge-centered w0 on a Cartesian grid,
            // for use in computing divergences.
            AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                         w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                         w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }
        }

        // w0r_cart is w0 but onto a Cartesian grid in cell-centered as
        // a scalar.  Since w0 is the radial expansion velocity, w0r_cart
        // is the radial w0 in a zone
        w0r_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        w0r_cart[lev].setVal(0.);
    }

    if (evolve_base_state == 1) {
        if (spherical == 1) {
            MakeW0mac(w0mac);
        }
        Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0);
    }

    // Mach number
    MachfromRhoH(s_in,u_in,p0_in,w0r_cart,tempmf);

    // MachNumber
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    ++dest_comp;

    // deltagamma
    MakeDeltaGamma(s_in, p0_in, p0_cart, gamma1bar_in, gamma1bar_cart, tempmf);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    ++dest_comp;

    // entropy
    MakeEntropy(s_in, tempmf);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    ++dest_comp;

    // entropypert = (entropy - entropybar) / entropybar
    {
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }

        Average(tempmf, tempbar_plot, 0);
        Put1dArrayOnCart(tempbar_plot,tempmf,0,0,bcs_f,0);

        for (int i = 0; i <= finest_level; ++i) {
            MultiFab::Subtract(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
            MultiFab::Divide(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
        }
    }
    ++dest_comp;

    if (plot_pidivu) {
        // pidivu
        MakePiDivu(u_in, s_in, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;
    }

    // processor number of each tile
    if (plot_processors) {
        for (int i = 0; i <= finest_level; ++i) {
            (*plot_mf_data[i]).setVal(ParallelDescriptor::MyProc());
        }
        ++dest_comp;
    }

    if (plot_ad_excess) {
        // ad_excess
        MakeAdExcess(s_in, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;
    }

    // S
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(S_cc_in[i],0,dest_comp,1);
    }
    ++dest_comp;

    // soundspeed
    if (plot_cs) {
        CsfromRhoH(s_in, p0_in, p0_cart, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;
    }

    // gravitational_acceleration
    if (plot_grav) {
        MakeGrav(rho0_new, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        ++dest_comp;
    }

    if (plot_base_state) {
        // w0
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(w0_cart[i],0,dest_comp,AMREX_SPACEDIM);
        }
        dest_comp += AMREX_SPACEDIM;

        // divw0
        MakeDivw0(w0mac, tempmf);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        }
        dest_comp++;
    }

    // thermal
    Vector<MultiFab> Tcoeff            (finest_level+1);
    Vector<MultiFab> hcoeff            (finest_level+1);
    Vector<MultiFab> Xkcoeff           (finest_level+1);
    Vector<MultiFab> pcoeff            (finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        Tcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
        hcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
        Xkcoeff           [lev].define(grids[lev], dmap[lev], NumSpec, 1);
        pcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
    }

    if (use_thermal_diffusion) {
        MakeThermalCoeffs(s_in,Tcoeff,hcoeff,Xkcoeff,pcoeff);
        MakeExplicitThermal(tempmf,s_in,Tcoeff,hcoeff,Xkcoeff,pcoeff,p0_in,0);
    } else {
        for (int lev=0; lev<=finest_level; ++lev) {
            Tcoeff[lev].setVal(0.);
            tempmf[lev].setVal(0.);
        }
    }
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
    }
    dest_comp++;

    // conductivity
    for (int i = 0; i <= finest_level; ++i) {
        tempmf[i].setVal(0.);
        plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
        MultiFab::Subtract(*plot_mf_data[i],Tcoeff[i],0,dest_comp,1,0);
    }
    dest_comp++;

    // radial and circular velocities
    if (spherical == 1) {
        MakeVelrc(u_in, w0r_cart, tempmf, tempmf_scalar1);
        for (int i = 0; i <= finest_level; ++i) {
            plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
            plot_mf_data[i]->copy(tempmf_scalar1[i],0,dest_comp+1,1);
        }
        dest_comp += 2;
    }

    if (do_sponge) {
        SpongeInit(rho0_old);
        MakeSponge(tempmf);

        if (plot_sponge_fdamp) {
            // compute f_damp assuming sponge=1/(1+dt*kappa*fdamp)
            // therefore fdamp = (1/sponge-1)/(dt*kappa)
            for (int i = 0; i <= finest_level; ++i) {
                // scalar1 = 1
                tempmf_scalar1[i].setVal(1.);
                // scalar2 = dt * kappa
                tempmf_scalar2[i].setVal(dt * sponge_kappa);
                // plot_mf = 1
                plot_mf_data[i]->copy(tempmf_scalar1[i],0,dest_comp,1);
                // plot_mf = 1/sponge
                MultiFab::Divide(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
                // plot_mf = 1/sponge - 1
                MultiFab::Subtract(*plot_mf_data[i],tempmf_scalar1[i],0,dest_comp,1,0);
                // plot_mf = (1/sponge-1)/(dt*kappa)
                MultiFab::Divide(*plot_mf_data[i],tempmf_scalar2[i],0,dest_comp,1,0);
            }
        } else {
            for (int i = 0; i <= finest_level; ++i) {
                plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
            }
        }
        dest_comp++;
    }

    // add plot_mf_data[i] to plot_mf
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf.push_back(plot_mf_data[i]);
    }

    return plot_mf;
}


// this takes the multifab of all variables and extracts those
// required for the small plot file
Vector<const MultiFab*>
Maestro::SmallPlotFileMF(const int nPlot, const int nSmallPlot,
                         Vector<const MultiFab*> mf,
                         const Vector<std::string> varnames,
                         const Vector<std::string> small_plot_varnames)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SmallPlotFileMF()",SmallPlotFileMF);

    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;

    // temporary MultiFabs to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level+1);

    int dest_comp = 0;

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] = new MultiFab(mf[i]->boxArray(),
                                       mf[i]->DistributionMap(),nSmallPlot,0);

    }

    for (auto it=small_plot_varnames.begin();
         it!=small_plot_varnames.end(); ++it) {

        for (auto n = 0; n < nPlot; n++) {
            if (*it == varnames[n]) {
                for (int i = 0; i <= finest_level; ++i)
                    plot_mf_data[i]->copy(*(mf[i]), n, dest_comp, 1);
                ++dest_comp;
                break;
            }
        }
    }

    // add plot_mf_data[i] to plot_mf
    for (int i = 0; i <= finest_level; ++i)
        plot_mf.push_back(plot_mf_data[i]);

    return plot_mf;
}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames (int * nPlot) const
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

    // velocities (AMREX_SPACEDIM)
    // magvel, momentum
    // rho, rhoh, h, rhoX, tfromp, tfromh, deltap, deltaT Pi (Nscal+4 -- the extra 4 are h, tfromh, deltap and deltaT)
    // rho' and rhoh' and t' (3)
    // pioverp0, p0pluspi (2)
    // MachNumber, deltagamma, entropy, entropypert, S
    // thermal, conductivity

    (*nPlot) = AMREX_SPACEDIM + Nscal + 19;

    if (plot_spec) (*nPlot) += NumSpec + 1; // X + 1 (abar)
    if (plot_spec || plot_omegadot) (*nPlot) += NumSpec; // omegadot

    if (plot_Hext) (*nPlot)++;
    if (plot_Hnuc) (*nPlot)++;
    if (plot_eta) (*nPlot)++;
    if (plot_gpi) (*nPlot) += AMREX_SPACEDIM;
    // rho0, rhoh0, h0, p0, w0, divw0 (5+AMREX_SPACEDIM)
    if (plot_base_state) (*nPlot) += AMREX_SPACEDIM + 5;
    if (plot_cs) (*nPlot)++;
    if (plot_grav) (*nPlot)++;
    if (plot_ad_excess) (*nPlot)++;
    if (plot_pidivu) (*nPlot)++;
    if (plot_processors) (*nPlot)++;
    if (spherical == 1) (*nPlot) += 2; // radial_velocity, circ_velocity
    if (do_sponge) (*nPlot)++;

    Vector<std::string> names(*nPlot);

    int cnt = 0;

    // add velocities
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        names[cnt++] = x;
    }

    names[cnt++] = "magvel";
    names[cnt++] = "momentum";

    names[cnt++] = "vort";

    // density and enthalpy
    names[cnt++] = "rho";
    names[cnt++] = "rhoh";
    names[cnt++] = "h";

    for (int i = 0; i < NumSpec; i++) {
        int len = 20;
        Vector<int> int_spec_names(len);
        //
        // This call return the actual length of each string in "len"
        //
        get_spec_names(int_spec_names.dataPtr(),&i,&len);
        char* spec_name = new char[len+1];
        for (int j = 0; j < len; j++)
            spec_name[j] = int_spec_names[j];
        spec_name[len] = '\0';
        std::string spec_string = "rhoX(";
        spec_string += spec_name;
        spec_string += ')';

        names[cnt++] = spec_string;
        delete [] spec_name;
    }

    if (plot_spec) {
        for (int i = 0; i < NumSpec; i++) {
            int len = 20;
            Vector<int> int_spec_names(len);
            //
            // This call return the actual length of each string in "len"
            //
            get_spec_names(int_spec_names.dataPtr(),&i,&len);
            char* spec_name = new char[len+1];
            for (int j = 0; j < len; j++) {
                spec_name[j] = int_spec_names[j];
            }
            spec_name[len] = '\0';
            std::string spec_string = "X(";
            spec_string += spec_name;
            spec_string += ')';

            names[cnt++] = spec_string;

            delete [] spec_name;
        }

        names[cnt++] = "abar";
    }

    if (plot_spec || plot_omegadot) {
        for (int i = 0; i < NumSpec; i++) {
            int len = 20;
            Vector<int> int_spec_names(len);
            //
            // This call return the actual length of each string in "len"
            //
            get_spec_names(int_spec_names.dataPtr(),&i,&len);
            char* spec_name = new char[len+1];
            for (int j = 0; j < len; j++) {
                spec_name[j] = int_spec_names[j];
            }
            spec_name[len] = '\0';
            std::string spec_string = "omegadot(";
            spec_string += spec_name;
            spec_string += ')';

            names[cnt++] = spec_string;

            delete [] spec_name;
        }
    }

    if (plot_Hext) names[cnt++] = "Hext";
    if (plot_Hnuc) names[cnt++] = "Hnuc";
    if (plot_eta) names[cnt++] = "eta_rho";

    names[cnt++] = "tfromp";
    names[cnt++] = "tfromh";
    names[cnt++] = "deltap";
    names[cnt++] = "deltaT";
    names[cnt++] = "Pi";
    names[cnt++] = "pioverp0";
    names[cnt++] = "p0pluspi";

    if (plot_gpi) {
        // add gpi
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            std::string x = "gpi";
            x += (120+i);
            names[cnt++] = x;
        }
    }

    names[cnt++] = "rhopert";
    names[cnt++] = "rhohpert";
    names[cnt++] = "tpert";
    if (plot_base_state) {
        names[cnt++] = "rho0";
        names[cnt++] = "rhoh0";
        names[cnt++] = "h0";
        names[cnt++] = "p0";
    }
    names[cnt++] = "MachNumber";
    names[cnt++] = "deltagamma";
    names[cnt++] = "entropy";
    names[cnt++] = "entropypert";
    if (plot_pidivu) names[cnt++] = "pi_divu";
    if (plot_processors) names[cnt++] = "processor_number";
    if (plot_ad_excess) names[cnt++] = "ad_excess";
    names[cnt++] = "S";

    if (plot_cs) names[cnt++] = "soundspeed";

    if (plot_grav) names[cnt++] = "maggrav";

    if (plot_base_state) {
        // w0 and divw0
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            std::string x = "w0";
            x += (120+i);
            names[cnt++] = x;
        }
        names[cnt++] = "divw0";
    }

    names[cnt++] = "thermal";
    names[cnt++] = "conductivity";

    if (spherical == 1) {
        names[cnt++] = "radial_velocity";
        names[cnt++] = "circ_velocity";
    }

    if (do_sponge) {
        if (plot_sponge_fdamp) {
            names[cnt++] = "sponge_fdamp";
        } else {
            names[cnt++] = "sponge";
        }
    }

    return names;

}

// set plotfile variable names
Vector<std::string>
Maestro::SmallPlotFileVarNames (int * nPlot, Vector<std::string> varnames) const
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SmallPlotFileVarNames()",SmallPlotFileVarNames);

    Vector<std::string> names(*nPlot);

    ParmParse pp("maestro");

    int nPltVars = pp.countval("small_plot_vars");

    if (nPltVars > 0) { // small_plot_vars defined in inputs file

        std::string nm;

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("small_plot_vars", nm, i);

            if (nm == "ALL")
                return varnames;
            else if (nm == "NONE") {
                names.clear();
                return names;
            } else {
                // test to see if it's a valid varname by iterating over
                // varnames
                auto found_name = false;
                for (auto it=varnames.begin(); it!=varnames.end(); ++it) {
                    if (nm == *it) {
                        names.push_back(nm);
                        found_name = true;
                        break;
                    }
                }

                if (!found_name)
                    Print() << "Small plot file variable " << nm << " is invalid\n";
            }
        }
    } else {
        // use default value of small_plot_vars which is a string that needs to be split
        std::stringstream sstream(small_plot_vars);
        std::string nm;

        while (sstream >> nm) {
            // test to see if it's a valid varname by iterating over
            // varnames
            auto found_name = false;
            for (auto it=varnames.begin(); it!=varnames.end(); ++it) {
                if (nm == *it) {
                    names.push_back(nm);
                    found_name = true;
                    break;
                }
            }

            if (!found_name)
                Print() << "Small plot file variable " << nm << " is invalid\n";
        }
    }

    names.shrink_to_fit();
    *nPlot = names.size();

    return names;

}

void
Maestro::WriteJobInfo (const std::string& dir) const
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WriteJobInfo()",WriteJobInfo);

    if (ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = dir;

        std::string PrettyLine = std::string(78, '=') + "\n";
        std::string OtherLine = std::string(78, '-') + "\n";
        std::string SkipSpace = std::string(8, ' ');

        FullPathJobInfoFile += "/job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " MAESTROeX Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "job name: " << job_name << "\n\n";
        jobInfoFile << "inputs file: " << inputs_name << "\n\n";

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";

        jobInfoFile << "tile size: ";
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            jobInfoFile << FabArrayBase::mfiter_tile_size[d] << " ";
        }
        jobInfoFile << "\n";
#endif
        jobInfoFile << "\n";
        jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
            getCPUTime()/3600.0;

        jobInfoFile << "\n\n";

        // plotfile information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Plotfile Information\n";
        jobInfoFile << PrettyLine;

        time_t now = time(0);

        // Convert now to tm struct for local timezone
        tm* localtm = localtime(&now);
        jobInfoFile << "output data / time: " << asctime(localtm);

        char currentDir[FILENAME_MAX];
        if (getcwd(currentDir, FILENAME_MAX)) {
            jobInfoFile << "output dir:         " << currentDir << "\n";
        }

        jobInfoFile << "\n\n";

#ifdef AMREX_USE_GPU
        // This output assumes for simplicity that every rank uses the
        // same type of GPU.

        jobInfoFile << PrettyLine;
        jobInfoFile << "GPU Information:       " << "\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "GPU model name: " << Gpu::Device::deviceName() << "\n";
        jobInfoFile << "Number of GPUs used: " << Gpu::Device::numDevicesUsed() << "\n";

        jobInfoFile << "\n\n";
#endif

        // build information

        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
        jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
        jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
        jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

        jobInfoFile << "\n";

        for (int n = 1; n <= buildInfoGetNumModules(); n++) {
            jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
        }

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        const char* githash3 = buildInfoGetGitHash(3);
        if (strlen(githash1) > 0) {
            jobInfoFile << "MAESTROeX git describe: " << githash1 << "\n";
        }
        if (strlen(githash2) > 0) {
            jobInfoFile << "AMReX git describe: " << githash2 << "\n";
        }
        if (strlen(githash3) > 0) {
            jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
        }

        const char* buildgithash = buildInfoGetBuildGitHash();
        const char* buildgitname = buildInfoGetBuildGitName();
        if (strlen(buildgithash) > 0) {
            jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
            {
                jobInfoFile << geom[i].Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << " Boundary conditions\n";
        Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
        ParmParse pp("maestro");
        pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
        pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


        // these names correspond to the integer flags setup in the
        // Castro_setup.cpp
        const char* names_bc[] =
            { "interior", "inflow", "outflow",
              "symmetry", "slipwall", "noslipwall" };


        jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
        jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
        if (BL_SPACEDIM >= 2) {
            jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
            jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
        }
        if (BL_SPACEDIM == 3) {
            jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
            jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
        }

        jobInfoFile << "\n\n";


        // species info
        Real Aion = 0.0;
        Real Zion = 0.0;

        int mlen = 20;

        jobInfoFile << PrettyLine;
        jobInfoFile << " Species Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile <<
            std::setw(6) << "index" << SkipSpace <<
            std::setw(mlen+1) << "name" << SkipSpace <<
            std::setw(7) << "A" << SkipSpace <<
            std::setw(7) << "Z" << "\n";
        jobInfoFile << OtherLine;

        for (int i = 0; i < NumSpec; i++)
        {

            int len = mlen;
            Vector<int> int_spec_names(len);
            //
            // This call return the actual length of each string in "len"
            //
            get_spec_names(int_spec_names.dataPtr(),&i,&len);
            char* spec_name = new char[len+1];
            for (int j = 0; j < len; j++)
                spec_name[j] = int_spec_names[j];
            spec_name[len] = '\0';

            // get A and Z
            get_spec_az(&i, &Aion, &Zion);

            jobInfoFile <<
                std::setw(6) << i << SkipSpace <<
                std::setw(mlen+1) << std::setfill(' ') << spec_name << SkipSpace <<
                std::setw(7) << Aion << SkipSpace <<
                std::setw(7) << Zion << "\n";
            delete [] spec_name;
        }
        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();

        // now the external parameters
        const int jobinfo_file_length = FullPathJobInfoFile.length();
        Vector<int> jobinfo_file_name(jobinfo_file_length);

        for (int i = 0; i < jobinfo_file_length; i++)
            jobinfo_file_name[i] = FullPathJobInfoFile[i];

        runtime_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);
    }
}

void
Maestro::WriteBuildInfo ()
{
    std::string PrettyLine = std::string(78, '=') + "\n";
    std::string OtherLine = std::string(78, '-') + "\n";
    std::string SkipSpace = std::string(8, ' ');

    // build information
    std::cout << PrettyLine;
    std::cout << " MAESTROeX Build Information\n";
    std::cout << PrettyLine;

    std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
    std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
    std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
    std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

    std::cout << "\n";

    std::cout << "COMP:          " << buildInfoGetComp() << "\n";
    std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

    std::cout << "\n";

    std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
    std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

    std::cout << "\n";

    std::cout << "Fortran comp:  " << buildInfoGetFName() << "\n";
    std::cout << "Fortran flags: " << buildInfoGetFFlags() << "\n";

    std::cout << "\n";

    std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
    std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

    std::cout << "\n";

    for (int n = 1; n <= buildInfoGetNumModules(); n++) {
        std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
    }

    const char* githash1 = buildInfoGetGitHash(1);
    const char* githash2 = buildInfoGetGitHash(2);
    const char* githash3 = buildInfoGetGitHash(3);
    if (strlen(githash1) > 0) {
        std::cout << "MAESTROeX git describe: " << githash1 << "\n";
    }
    if (strlen(githash2) > 0) {
        std::cout << "AMReX git describe: " << githash2 << "\n";
    }
    if (strlen(githash3) > 0) {
        std::cout << "Microphysics git describe: " << githash3 << "\n";
    }

    const char* buildgithash = buildInfoGetBuildGitHash();
    const char* buildgitname = buildInfoGetBuildGitName();
    if (strlen(buildgithash) > 0) {
        std::cout << buildgitname << " git describe: " << buildgithash << "\n";
    }

    std::cout << "\n\n";
}

void
Maestro::MakeMagvel (const Vector<MultiFab>& vel,
                     Vector<MultiFab>& magvel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeMagvel()", MakeMagvel);

#if (AMREX_SPACEDIM == 3)

    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        }
        MakeW0mac(w0mac);
    }
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        if (!spherical) {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                const Array4<const Real> vel_arr = vel[lev].array(mfi);
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<Real> magvel_arr = magvel[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    Real v_total = vel_arr(i,j,k,1) + 0.5 * (w0_arr(i,j,k,1) + w0_arr(i,j+1,k,1));
                    magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) + 
                        v_total*v_total);
#else
                    Real w_total = vel_arr(i,j,k,2) + 0.5 * (w0_arr(i,j,k,2) + w0_arr(i,j,k+1,2));
                    magvel_arr(i,j,k) = sqrt(vel_arr(i,j,k,0)*vel_arr(i,j,k,0) + 
                        vel_arr(i,j,k,1)*vel_arr(i,j,k,1) + 
                        w_total*w_total);
#endif
                });
            }
        } else {
#if (AMREX_SPACEDIM == 3)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                const Array4<const Real> vel_arr = vel[lev].array(mfi);
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
                const Array4<Real> magvel_arr = magvel[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real u_total = vel_arr(i,j,k,0) + 0.5 * (w0macx(i,j,k) + w0macx(i+1,j,k));
                    Real v_total = vel_arr(i,j,k,1) + 0.5 * (w0macy(i,j,k) + w0macy(i,j+1,k));
                    Real w_total = vel_arr(i,j,k,2) + 0.5 * (w0macz(i,j,k) + w0macz(i,j,k+1));
                    magvel_arr(i,j,k) = sqrt(u_total*u_total + 
                        v_total*v_total + w_total*w_total);
                });
            }
#endif
        }
    }

    // average down and fill ghost cells
    AverageDown(magvel, 0, 1);
    FillPatch(t_old, magvel, magvel, magvel, 0, 0, 1, 0, bcs_f);
}


void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

    for (int lev=0; lev<=finest_level; ++lev) {

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
            const Array4<Real> circvel_arr = circ_vel[lev].array(mfi);
            const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                circvel_arr(i,j,k) = 0.0;
                radvel_arr(i,j,k) = 0.0;

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    radvel_arr(i,j,k) += vel_arr(i,j,k,n) * normal_arr(i,j,k,n);
                }
                
                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    Real circ_comp = vel_arr(i,j,k,n) - radvel_arr(i,j,k) * normal_arr(i,j,k,n);
                    circvel_arr(i,j,k) += circ_comp * circ_comp;
                }

                circvel_arr(i,j,k) = sqrt(circvel_arr(i,j,k));

                // add base state vel to get full radial velocity
                radvel_arr(i,j,k) += w0rcart_arr(i,j,k);
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(rad_vel, 0, 1);
    FillPatch(t_old, rad_vel, rad_vel, rad_vel, 0, 0, 1, 0, bcs_f);
    AverageDown(circ_vel, 0, 1);
    FillPatch(t_old, circ_vel, circ_vel, circ_vel, 0, 0, 1, 0, bcs_f);
}


void
Maestro::MakeAdExcess (const Vector<MultiFab>& state,
                       Vector<MultiFab>& ad_excess)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeAdExcess()", MakeAdExcess);

    const auto base_cutoff_density_loc = base_cutoff_density;

    for (int lev=0; lev<=finest_level; ++lev) {

        // create MultiFabs to hold pressure and gradient
        MultiFab pres_mf(grids[lev], dmap[lev], 1, 0);
        MultiFab nabla_ad_mf(grids[lev], dmap[lev], 1, 0);

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<Real> ad_excess_arr = ad_excess[lev].array(mfi);
            const Array4<Real> pres = pres_mf.array(mfi);
            const Array4<Real> nabla_ad = nabla_ad_mf.array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<const Real> normal_arr = normal[lev].array(mfi);
#endif

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                eos_t eos_state;

                eos_state.rho   = state_arr(i,j,k,Rho);
                eos_state.T     = state_arr(i,j,k,Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
                }

                eos(eos_input_rt, eos_state);

                pres(i,j,k) = eos_state.p;
                // Print() << "pres = " << pres(i,j,k) << std::endl;

                Real chi_rho = eos_state.rho * eos_state.dpdr / eos_state.p;
                Real chi_t = eos_state.T * eos_state.dpdT / eos_state.p;
                nabla_ad(i,j,k) = (eos_state.gam1 - chi_rho) / (chi_t * eos_state.gam1);
            });

            const auto lo = tileBox.loVect3d();
            const auto hi = tileBox.hiVect3d();

            if (spherical == 0) {
                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real nabla = 0.0;

                    if (state_arr(i,j,k,Rho) > base_cutoff_density_loc) {
                        Real dtemp = 0.0;
                        Real dp = 0.0;
#if (AMREX_SPACEDIM == 2)
                        // forward difference
                        if (j == lo[1]) {
                            dtemp = state_arr(i,j+1,k,Temp) - state_arr(i,j,k,Temp);
                            dp = pres(i,j+1,k) - pres(i,j,k);
                        // backward difference
                        } else if (j == hi[1]) {
                            dtemp = state_arr(i,j,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp = pres(i,j,k) - pres(i,j-1,k);
                        // centered difference
                        } else {
                            dtemp = state_arr(i,j+1,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp = pres(i,j+1,k) - pres(i,j-1,k);
                        }
#else 
                        // forward difference
                        if (k == lo[2]) {
                            dtemp = state_arr(i,j,k+1,Temp) - state_arr(i,j,k,Temp);
                            dp = pres(i,j,k+1) - pres(i,j,k);
                        // backward difference
                        } else if (k == hi[2]) {
                            dtemp = state_arr(i,j,k,Temp) - state_arr(i,j,k-1,Temp);
                            dp = pres(i,j,k) - pres(i,j,k-1);
                        // centered difference
                        } else {
                            dtemp = state_arr(i,j,k+1,Temp) - state_arr(i,j,k-1,Temp);
                            dp = pres(i,j,k+1) - pres(i,j,k-1);
                        }
#endif
                        // prevent Inf
                        if (dp * state_arr(i,j,k,Temp) == 0.0) {
                            nabla = std::numeric_limits<Real>::min();
                        } else {
                            nabla = pres(i,j,k) * dtemp / (dp * state_arr(i,j,k,Temp));
                        }
                    }

                    ad_excess_arr(i,j,k) = nabla - nabla_ad(i,j,k);
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                RealVector dtemp_vec(AMREX_SPACEDIM, 0.0);
                RealVector dp_vec(AMREX_SPACEDIM, 0.0);

                Real * AMREX_RESTRICT dtemp = dtemp_vec.dataPtr();
                Real * AMREX_RESTRICT dp = dp_vec.dataPtr();

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    Real nabla = 0.0;

                    if (state_arr(i,j,k,Rho) > base_cutoff_density_loc) {
                        // compute gradient
                        // forward difference
                        if (i == lo[0]) {
                            dtemp[0] = state_arr(i+1,j,k,Temp) - state_arr(i,j,k,Temp);
                            dp[0] = pres(i+1,j,k) - pres(i,j,k);
                        // backward difference
                        } else if (i == hi[0]) {
                            dtemp[0] = state_arr(i,j,k,Temp) - state_arr(i-1,j,k,Temp);
                            dp[0] = pres(i,j,k) - pres(i-1,j,k);
                        // centered difference
                        } else {
                            dtemp[0] = state_arr(i+1,j,k,Temp) - state_arr(i-1,j,k,Temp);
                            dp[0] = pres(i+1,j,k) - pres(i-1,j,k);
                        }
                        // forward difference
                        if (j == lo[1]) {
                            dtemp[1] = state_arr(i,j+1,k,Temp) - state_arr(i,j,k,Temp);
                            dp[1] = pres(i,j+1,k) - pres(i,j,k);
                        // backward difference
                        } else if (j == hi[1]) {
                            dtemp[1] = state_arr(i,j,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp[1] = pres(i,j,k) - pres(i,j-1,k);
                        // centered difference
                        } else {
                            dtemp[1] = state_arr(i,j+1,k,Temp) - state_arr(i,j-1,k,Temp);
                            dp[1] = pres(i,j+1,k) - pres(i,j-1,k);
                        }
                        // forward difference
                        if (k == lo[2]) {
                            dtemp[2] = state_arr(i,j,k+1,Temp) - state_arr(i,j,k,Temp);
                            dp[2] = pres(i,j,k+1) - pres(i,j,k);
                        // backward difference
                        } else if (k == hi[2]) {
                            dtemp[2] = state_arr(i,j,k,Temp) - state_arr(i,j,k-1,Temp);
                            dp[2] = pres(i,j,k) - pres(i,j,k-1);
                        // centered difference
                        } else {
                            dtemp[2] = state_arr(i,j,k+1,Temp) - state_arr(i,j,k-1,Temp);
                            dp[2] = pres(i,j,k+1) - pres(i,j,k-1);
                        }

                        Real dp_dot = 0.0;
                        Real dtemp_dot = 0.0;
                        for (auto c = 0; c < AMREX_SPACEDIM; ++c) {
                            dp_dot += dp[c] * normal_arr(i,j,k,c);
                            dtemp_dot += dtemp[c] * normal_arr(i,j,k,c);
                        }

                        // prevent Inf
                        if (dp_dot * state_arr(i,j,k,Temp) == 0.0) {
                            nabla = std::numeric_limits<Real>::min();
                        } else {
                            nabla = pres(i,j,k) * dtemp_dot / (dp_dot * state_arr(i,j,k,Temp));
                        }
                    }

                    ad_excess_arr(i,j,k) = nabla - nabla_ad(i,j,k);
                });
#endif
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(ad_excess,0,1);
    FillPatch(t_old,ad_excess,ad_excess,ad_excess,0,0,1,0,bcs_f);
}


void
Maestro::MakeGrav (const RealVector& rho0,
                   Vector<MultiFab>& grav)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGrav()", MakeGrav);

    RealVector grav_cell( (base_geom.max_radial_level+1)*base_geom.nr_fine );
    grav_cell.shrink_to_fit();

    MakeGravCell(grav_cell, rho0);

    Put1dArrayOnCart(grav_cell,grav,0,0,bcs_f,0);

    // average down and fill ghost cells
    AverageDown(grav,0,1);
    FillPatch(t_old,grav,grav,grav,0,0,1,0,bcs_f);
}


void
Maestro::MakeVorticity (const Vector<MultiFab>& vel,
                        Vector<MultiFab>& vorticity)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVorticity()",MakeVorticity);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];

        const Real* dx = geom[lev].CellSize();
        const Box& domainBox = geom[lev].Domain();

        const Real hx = dx[0];
        const Real hy = dx[1];
#if (AMREX_SPACEDIM == 3)
        const Real hz = dx[2];
#endif
        const int ilo = domainBox.loVect()[0];
        const int ihi = domainBox.hiVect()[0];
        const int jlo = domainBox.loVect()[1];
        const int jhi = domainBox.hiVect()[1];
#if (AMREX_SPACEDIM == 3)
        const int klo = domainBox.loVect()[2];
        const int khi = domainBox.hiVect()[2];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            Array4<const Real> const u = vel[lev].array(mfi);
            Array4<Real> const vort = vorticity[lev].array(mfi);
            GpuArray<int,AMREX_SPACEDIM*2> physbc;
            for (int n = 0; n < AMREX_SPACEDIM*2; ++n) {
                physbc[n] = phys_bc[n];
            } 

#if (AMREX_SPACEDIM == 2)

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;

                if (i == ilo && 
                    (physbc[0] == Inflow || 
                     physbc[0] == SlipWall || 
                     physbc[0] == NoSlipWall)) 
                {
                    vx = (u(i+1,j,k,1) + 3.0*u(i,j,k,1) - 
                          4.0*u(i-1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;

                } else if (i == ihi+1 &&
                        (physbc[AMREX_SPACEDIM] == Inflow || 
                        physbc[AMREX_SPACEDIM] == SlipWall || 
                        physbc[AMREX_SPACEDIM] == NoSlipWall))
                {
                    vx = -(u(i-1,j,k,1) + 3.0*u(i,j,k,1) - 
                         4.0*u(i+1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;
                }

                if (j == jlo &&
                    (physbc[1] == Inflow || 
                     physbc[1] == SlipWall || 
                     physbc[1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,0)) / hx;
                    uy = (u(i,j+1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j-1,k,0)) / hy;

                } else if (j == jhi+1 && 
                           (physbc[AMREX_SPACEDIM+1] == Inflow || 
                            physbc[AMREX_SPACEDIM+1] == SlipWall || 
                            physbc[AMREX_SPACEDIM+1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,1)) / hx;
                    uy = -(u(i,j-1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j+1,k,0)) / hy;
                }

                vort(i,j,k) = vx - uy;
            });

#else 
            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;
                Real uz = 0.5*(u(i,j,k+1,0)-u(i,j,k-1,0))/hz;
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real vz = 0.5*(u(i,j,k+1,1)-u(i,j,k-1,1))/hz;
                Real wx = 0.5*(u(i+1,j,k,2)-u(i-1,j,k,2))/hx;
                Real wy = 0.5*(u(i,j+1,k,2)-u(i,j-1,k,2))/hy;

                bool fix_lo_x = (physbc[0] == Inflow || 
                                 physbc[0] == NoSlipWall);
                bool fix_hi_x = (physbc[AMREX_SPACEDIM] == Inflow || 
                                 physbc[AMREX_SPACEDIM] == NoSlipWall);

                bool fix_lo_y = (physbc[1] == Inflow || 
                                 physbc[1] == NoSlipWall);
                bool fix_hi_y = (physbc[AMREX_SPACEDIM+1] == Inflow ||
                                 physbc[AMREX_SPACEDIM+1] == NoSlipWall);

                bool fix_lo_z = (physbc[2] == Inflow || 
                                 physbc[2] == NoSlipWall);
                bool fix_hi_z = (physbc[AMREX_SPACEDIM+2] == Inflow ||
                                 physbc[AMREX_SPACEDIM+2] == NoSlipWall);

                // First do all the faces
                if (fix_lo_x && i == ilo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                } else if (fix_hi_x && i == ihi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                }

                if (fix_lo_y && j == jlo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                } else if (fix_hi_y && j == jhi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_z && k == klo) {
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_z && k == khi+1) {
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                // Next do all the edges
                if (fix_lo_x && fix_lo_y && i == ilo && j == jlo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_hi_x && fix_lo_y && i == ihi+1 && j == jlo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_hi_y && i == ilo && j == jhi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_lo_z && i == ilo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_z && i == ihi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_z && i == ilo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_z && i == ihi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_lo_z && j == jlo && k == klo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_lo_z && j == jhi+1 && k == klo) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_hi_z && j == jlo && k == khi+1) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_hi_z && j == jhi+1 && k == khi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }
                
                // Finally do all the corners
                if (fix_lo_x && fix_lo_y && fix_lo_z && 
                    i == ilo && j == jlo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_lo_z &&
                    i == ihi+1 && j == jlo && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_lo_z &&
                    i == ilo && j == jhi+1 && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_lo_z &&
                    i == ihi+1 && j == jhi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_lo_y && fix_hi_z &&
                    i == ilo && j == jlo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_hi_z &&
                    i == ihi+1 && j == jlo && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_hi_z &&
                    i == ilo && j == jhi+1 && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_hi_z &&
                    i == ihi+1 && j == jhi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                vort(i,j,k) = sqrt((wy-vz)*(wy-vz)+
                    (uz-wx)*(uz-wx)+(vx-uy)*(vx-uy));
            });
#endif
        }
    }

    // average down and fill ghost cells
    AverageDown(vorticity,0,1);
    FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeDeltaGamma (const Vector<MultiFab>& state,
                         const RealVector& p0,
                         const Vector<MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const Vector<MultiFab>& gamma1bar_cart,
                         Vector<MultiFab>& deltagamma)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDeltaGamma()", MakeDeltaGamma);

    const auto use_pprime_in_tfromp_loc = use_pprime_in_tfromp;

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr = gamma1bar_cart[lev].array(mfi);
            const Array4<Real> deltagamma_arr = deltagamma[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                eos_t eos_state;

                eos_state.rho   = state_arr(i,j,k,Rho);
                eos_state.T     = state_arr(i,j,k,Temp);
                if (use_pprime_in_tfromp_loc) {
                    eos_state.p     = p0_arr(i,j,k) + state_arr(i,j,k,Pi);
                } else {
                    eos_state.p     = p0_arr(i,j,k);
                }

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
                }

                eos(eos_input_rp, eos_state);

                deltagamma_arr(i,j,k) = eos_state.gam1 - gamma1bar_arr(i,j,k);
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(deltagamma, 0, 1);
    FillPatch(t_old, deltagamma, deltagamma, deltagamma, 0, 0, 1, 0, bcs_f);
}

void
Maestro::MakeEntropy (const Vector<MultiFab>& state,
                      Vector<MultiFab>& entropy)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEntropy()",MakeEntropy);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<Real> entropy_arr = entropy[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j ,k, {
                eos_t eos_state;

                eos_state.rho = state_arr(i,j,k,Rho);
                eos_state.T = state_arr(i,j,k,Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp) / state_arr(i,j,k,Rho);
                }

                eos(eos_input_rt, eos_state);

                entropy_arr(i,j,k) = eos_state.s;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(entropy,0,1);
    FillPatch(t_old,entropy,entropy,entropy,0,0,1,0,bcs_f);
}

void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDivw0()", MakeDivw0);

    for (int lev=0; lev<=finest_level; ++lev) {

        if (spherical == 0) {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    divw0_arr(i,j,k) = (w0_arr(i,j+1,k,1) - w0_arr(i,j,k,1)) / dx[1];
#else
                    divw0_arr(i,j,k) = (w0_arr(i,j,k+1,2) - w0_arr(i,j,k,2)) / dx[2];
#endif
                });
            }

        } else {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    divw0_arr(i,j,k) = (w0macx(i+1,j,k) - w0macx(i,j,k)) / dx[0] + 
                        (w0macy(i,j+1,k) - w0macy(i,j,k)) / dx[1] + 
                        (w0macz(i,j,k+1) - w0macz(i,j,k)) / dx[2];
                });
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(divw0,0,1);
    FillPatch(t_old,divw0,divw0,divw0,0,0,1,0,bcs_f);
}

void
Maestro::MakePiDivu (const Vector<MultiFab>& vel,
                     const Vector<MultiFab>& state,
                     Vector<MultiFab>& pidivu)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakePiDivu()",MakePiDivu);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(pidivu[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<const Real> pi_cc = state[lev].array(mfi, Pi);
            const Array4<Real> pidivu_arr = pidivu[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                pidivu_arr(i,j,k) = pi_cc(i,j,k) * 0.5 * (
                    (vel_arr(i+1,j,k,0) - vel_arr(i-1,j,k,0))/dx[0]
                  + (vel_arr(i,j+1,k,1) - vel_arr(i,j-1,k,1))/dx[1]
#if (AMREX_SPACEDIM == 3)
                  + (vel_arr(i,j,k+1,2) - vel_arr(i,j,k-1,2))/dx[2]
#endif
                    );
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(pidivu, 0, 1);
    FillPatch(t_old, pidivu, pidivu, pidivu, 0, 0, 1, 0, bcs_f);
}

void
Maestro::MakeAbar (const Vector<MultiFab>& state,
                   Vector<MultiFab>& abar)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeAbar()", MakeAbar);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(abar[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<Real> abar_arr = abar[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                Real denominator = 0.0;

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    denominator += state_arr(i,j,k,FirstSpec+comp) / aion[comp];
                }

                abar_arr(i,j,k) = state_arr(i,j,k,Rho) / denominator;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(abar, 0, 1);
    FillPatch(t_old, abar, abar, abar, 0, 0, 1, 0, bcs_f);
}
