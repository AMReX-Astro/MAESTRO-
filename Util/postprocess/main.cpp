//
// Process a plotfile to produce diagnostics for rotational problems.
//
#include <fstream>
#include <iostream>

#include <Postprocess.H>
#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    Postprocess postproc;

    postproc.init();

    postproc.diag();

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;

    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total,
                                      ParallelDescriptor::IOProcessorNumber());
    Print() << "Total Time: " << end_total << '\n';

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
