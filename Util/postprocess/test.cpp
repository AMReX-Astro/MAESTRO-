#include <Maestro.H>
#include <Maestro_F.H>
#include <Postprocess.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// ---------------------------------
// Write 1D radial diagnostics
// ---------------------------------
void Postprocess::test() {
    // exact-solution problem to test subroutines that compute
    // diagnostics dependent on velocity only;

    // make BoxArray and Geometry
    int n_cell = 128;
    int max_grid_size = 32;
    BoxArray ba;
    Geometry tgeom;
    {
        IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
        IntVect dom_hi(AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1));
        Box domain(dom_lo, dom_hi);

        ba.define(domain);
        ba.maxSize(max_grid_size);

        // define physical box, [0, 200] in each direction
        RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                         {AMREX_D_DECL(200.0_rt, 200.0_rt, 200.0_rt)});

        tgeom.define(domain, &real_box);
    }
    DistributionMapping dm(ba);

    // define density and pressure
    MultiFab rho0(ba, dm, 1, 0);
    MultiFab p0(ba, dm, 1, 0);
    rho0.setVal(1.);
    p0.setVal(1.);

    // define species
    maestro_network_init();
    MultiFab rhoX(ba, dm, NumSpec, 0);
    rhoX.setVal(0.);

    // define velocities
    MultiFab u_mf(ba, dm, AMREX_SPACEDIM, 0);
    MultiFab w0_mf(ba, dm, AMREX_SPACEDIM, 0);
    w0_mf.setVal(0.);

    Vector<Real> center_p(3);
    center_p[0] = 100.0;
    center_p[1] = 100.0;
    center_p[2] = 100.0;

    const auto dx = tgeom.CellSizeArray();
    const auto prob_lo = tgeom.ProbLoArray();

    for (MFIter mfi(u_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Get the index space of the valid region
        const Box& tileBox = mfi.tilebox();

        const Array4<Real> rhoX_arr = rhoX.array(mfi);
        const Array4<Real> vel_arr = u_mf.array(mfi);

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            // first species: rhoX = 1
            rhoX_arr(i, j, k, 0) = 1.0;

            // initialize velocity
            Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
            Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
            Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
            Real radius = std::sqrt(x * x + y * y + z * z);

            Vector<Real> norm(3);
            norm[0] = x / radius;
            norm[1] = y / radius;
            norm[2] = z / radius;

            // v_r = 1/r       if r < 50
            // v_r = 50*(50/r) if r >= 50
            if (radius < 50.0) {
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    vel_arr(i, j, k, dim) = radius * norm[dim];
                }
            } else {
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    vel_arr(i, j, k, dim) = 50.0 * 50.0 / radius * norm[dim];
                }
            }
        });
    }

    // put all variables into a single multifab
    MultiFab s0(ba, dm, 2 + NumSpec + 2 * AMREX_SPACEDIM, 0);
    MultiFab::Copy(s0, rho0, 0, 0, 1, 0);
    MultiFab::Copy(s0, p0, 0, 1, 1, 0);
    MultiFab::Copy(s0, rhoX, 0, 2, NumSpec, 0);
    MultiFab::Copy(s0, u_mf, 0, 2 + NumSpec, AMREX_SPACEDIM, 0);
    MultiFab::Copy(s0, w0_mf, 0, 2 + NumSpec + AMREX_SPACEDIM, AMREX_SPACEDIM,
                   0);

    // variable names for plotfile
    Vector<std::string> varnames(2 + NumSpec + 2 * AMREX_SPACEDIM);
    int cnt = 0;
    varnames[cnt++] = "rho";
    varnames[cnt++] = "p0";

    for (int i = 0; i < NumSpec; i++) {
        int len = 20;
        Vector<int> int_spec_names(len);
        //
        // This call return the actual length of each string in "len"
        //
        get_spec_names(int_spec_names.dataPtr(), &i, &len);
        auto* spec_name = new char[len + 1];
        for (int j = 0; j < len; j++) {
            spec_name[j] = int_spec_names[j];
        }
        spec_name[len] = '\0';
        std::string spec_string = "rhoX(";
        spec_string += spec_name;
        spec_string += ')';

        varnames[cnt++] = spec_string;
        delete[] spec_name;
    }

    // add velocities
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120 + i);
        varnames[cnt++] = x;
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        std::string x = "w0";
        x += (120 + i);
        varnames[cnt++] = x;
    }

    // write plotfile
    Print() << "Writing exact solution plotfile" << std::endl;

    std::string basefilename = "test_plt";
    std::string testfilename = basefilename + "0000000";

    WriteSingleLevelPlotfile(testfilename, s0, varnames, tgeom, 0, 0);

    WriteTestJobInfo(testfilename, basefilename, max_grid_size);
}

void Postprocess::WriteTestJobInfo(const std::string& dir,
                                   const std::string& base,
                                   const int maxgridsize) {
    // timer for profiling
    BL_PROFILE_VAR("Postprocess::WritePartJobInfo()", WritePartJobInfo);

    if (ParallelDescriptor::IOProcessor()) {
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

        // // grid information
        // jobInfoFile << PrettyLine;
        // jobInfoFile << " Grid Information\n";
        // jobInfoFile << PrettyLine;

        // for (int i = 0; i <= finest_level; i++) {
        //     jobInfoFile << " level: " << i << "\n";
        //     jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
        //     jobInfoFile << "   maximum zones   = ";
        //     for (int n = 0; n < BL_SPACEDIM; n++) {
        //         jobInfoFile << geom_in[i].Domain().length(n) << " ";
        //     }
        //     jobInfoFile << "\n\n";
        // }

        // jobInfoFile << " Boundary conditions\n";
        // Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
        // ParmParse pp("maestro");
        // pp.getarr("lo_bc", lo_bc_out, 0, BL_SPACEDIM);
        // pp.getarr("hi_bc", hi_bc_out, 0, BL_SPACEDIM);

        // // these names correspond to the integer flags setup in the
        // // Castro_setup.cpp
        // const char* names_bc[] = {"interior", "inflow",   "outflow",
        //                           "symmetry", "slipwall", "noslipwall"};

        // jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
        // jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
        // if (BL_SPACEDIM >= 2) {
        //     jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
        //     jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
        // }
        // if (BL_SPACEDIM == 3) {
        //     jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
        //     jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
        // }

        // jobInfoFile << "\n\n";

        // // species info
        // Real Aion = 0.0;
        // Real Zion = 0.0;

        // int mlen = 20;

        // jobInfoFile << PrettyLine;
        // jobInfoFile << " Species Information\n";
        // jobInfoFile << PrettyLine;

        // jobInfoFile << std::setw(6) << "index" << SkipSpace
        //             << std::setw(mlen + 1) << "name" << SkipSpace
        //             << std::setw(7) << "A" << SkipSpace << std::setw(7) << "Z"
        //             << "\n";
        // jobInfoFile << OtherLine;

        // for (int i = 0; i < NumSpec; i++) {
        //     int len = mlen;
        //     Vector<int> int_spec_names(len);
        //     //
        //     // This call return the actual length of each string in "len"
        //     //
        //     get_spec_names(int_spec_names.dataPtr(), &i, &len);
        //     auto* spec_name = new char[len + 1];
        //     for (int j = 0; j < len; j++) {
        //         spec_name[j] = int_spec_names[j];
        //     }
        //     spec_name[len] = '\0';

        //     // get A and Z
        //     get_spec_az(&i, &Aion, &Zion);

        //     jobInfoFile << std::setw(6) << i << SkipSpace << std::setw(mlen + 1)
        //                 << std::setfill(' ') << spec_name << SkipSpace
        //                 << std::setw(7) << Aion << SkipSpace << std::setw(7)
        //                 << Zion << "\n";
        //     delete[] spec_name;
        // }
        // jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile << "eos_gamma = 1.666666667\n";
        jobInfoFile << "amr.max_grid_size = " << maxgridsize << "\n";
        jobInfoFile << "maestro.plot_base_name = " << base << "\n";
        jobInfoFile << "maestro.drdxfac = 5\n";
        jobInfoFile << "maestro.octant = false\n";
        jobInfoFile << "maestro.rotational_frequency = 1.e-6\n";

        jobInfoFile.close();
    }
}
