#include <Maestro.H>
#include <Maestro_F.H>
#include <Postprocess.H>

using namespace amrex;

// helper IntVects used to define face/nodal MultiFabs
#if (AMREX_SPACEDIM == 2)
IntVect Postprocess::nodal_flag(1, 1);
IntVect Postprocess::nodal_flag_x(1, 0);
IntVect Postprocess::nodal_flag_y(0, 1);
#elif (AMREX_SPACEDIM == 3)
IntVect Postprocess::nodal_flag(1, 1, 1);
IntVect Postprocess::nodal_flag_x(1, 0, 0);
IntVect Postprocess::nodal_flag_y(0, 1, 0);
IntVect Postprocess::nodal_flag_z(0, 0, 1);
#endif

// constructor
Postprocess::Postprocess() = default;
// destructor
Postprocess::~Postprocess() = default;

// initialize data
void Postprocess::init() {
    // read in parameters from inputs file
    ParmParse pp;
    pp.query("infile", iFile);
    if (iFile.empty()) {
        Print() << "WARNING: Plotfile was not specified!\n";
        Print() << "Running exact solution test problem" << std::endl;
        test();
        iFile = "test_plt0000000";
    } else {
        // initialize species
        maestro_network_init();
    }

    // optional parameters
    pp.query("modelfile", imFile);
    pp.query("dt", deltat);
    pp.query("nfiles", numfiles);

    // read input grid
    amrex::PlotFileData pltfile(iFile);
    finest_level = pltfile.finestLevel();
    const auto probLo = pltfile.probLo();
    const auto probHi = pltfile.probHi();

    grid.resize(finest_level + 1);
    dmap.resize(finest_level + 1);
    pgeom.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        grid[lev] = pltfile.boxArray(lev);
        dmap[lev] = pltfile.DistributionMap(lev);

        Box domain = pltfile.probDomain(lev);
        RealBox real_box(probLo, probHi);
        pgeom[lev].define(domain, &real_box);
    }

    // setup geometry
    const auto domainBoxFine = pltfile.probDomain(finest_level);
    const auto dxFine = pltfile.cellSize(finest_level);

    // assume spherical
    int drdxfac = GetdrdxFac(iFile);
    bool octant = GetOctant(iFile);
    base_geom.max_radial_level = 0;
    base_geom.dr_fine = dxFine[0] / drdxfac;
    int domhi = domainBoxFine.bigEnd(0) + 1;

    // assume constant dr
    if (!octant) {
        base_geom.nr_irreg = int(
            round((3 * (domhi / 2 - 0.5) * (domhi / 2 - 0.5) - 0.75) / 2.0));
    } else {
        base_geom.nr_irreg =
            int(round((3 * (domhi - 0.5) * (domhi - 0.5) - 0.75) / 2.0));
    }

    double lenx, leny, lenz, max_dist;
    if (octant) {
        lenx = probHi[0] - probLo[0];
        leny = probHi[1] - probLo[1];
        lenz = probHi[2] - probLo[2];
    } else {
        lenx = 0.5 * (probHi[0] - probLo[0]);
        leny = 0.5 * (probHi[1] - probLo[1]);
        lenz = 0.5 * (probHi[2] - probLo[2]);
    }
    max_dist = sqrt(lenx * lenx + leny * leny + lenz * lenz);
    base_geom.nr_fine = int(max_dist / base_geom.dr_fine) + 1;

    base_geom.Init(base_geom.max_radial_level, base_geom.nr_fine,
                   base_geom.dr_fine, base_geom.nr_irreg, pgeom, finest_level,
                   center);
}

// write diagnostics
void Postprocess::diag() {
    amrex::PlotFileData pltfile(iFile);

    // read input MultiFabs
    Vector<MultiFab> rho_mf(finest_level + 1);
    Vector<MultiFab> p0_mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rho_mf[lev] = pltfile.get(lev, "rho");
        p0_mf[lev] = pltfile.get(lev, "p0");
    }

    Vector<MultiFab> u_mf(finest_level + 1);
    Vector<MultiFab> w0_mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        u_mf[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
        w0_mf[lev].define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }

    // velocities
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        std::string w = "w0";
        x += (120 + i);
        w += (120 + i);

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(u_mf[lev], pltfile.get(lev, x), 0, i, 1, 0);
            MultiFab::Copy(w0_mf[lev], pltfile.get(lev, w), 0, i, 1, 0);
        }
    }

    // Radial states
    BaseState<Real> rho0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> p0(base_geom.max_radial_level + 1, base_geom.nr_fine);

    Average(rho_mf, rho0, 0);
    Average(p0_mf, p0, 0);

    // Write radial output file
    Print() << "Writing radial diag file" << std::endl;
    WriteRadialFile(rho0, p0, u_mf, w0_mf);

    // Write 2D slice output file
    if (deltat.empty() || numfiles.empty()) {
        Print() << "Writing 2D slice file for single plot file" << std::endl;
        Write2dSliceFile(rho_mf, p0_mf, u_mf, w0_mf);
    } else {
        Print() << "Writing 2D slice file for " << stoi(numfiles) + 1
                << " plot files" << std::endl;
        Write2dSliceFile(rho_mf, p0_mf, u_mf, w0_mf, stoi(deltat),
                         stoi(numfiles));
    }

    // Write diag file for initial model if specified
    if (!imFile.empty()) {
        Print() << "Writing diag file for initial model" << std::endl;
        WriteModelDiagFile(imFile);
    }
}
