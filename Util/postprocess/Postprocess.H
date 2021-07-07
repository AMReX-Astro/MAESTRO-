#include <BaseState.H>
#include <BaseStateGeometry.H>

class Postprocess {
   public:
    /* public member functions */

    /// constructor
    Postprocess();

    /// destructor
    virtual ~Postprocess();

    /// initialize simulation by reading in plot file
    void init();

    /// write out diagnostics
    void diag();

    // ----------------------------
    // InputFile.cpp functions
    // ----------------------------

    // Read parameters from job info file
    std::string GetVarFromJobInfo(const std::string pltfile,
                                  const std::string& varname);

    int GetdrdxFac(const std::string pltfile);

    bool GetOctant(const std::string pltfile);

    amrex::Real GetGamma(const std::string pltfile);

    amrex::Real GetRotationFreq(const std::string pltfile);

    int GetTimeStep(const std::string plotfile, const std::string basefile);

    int GetMaxGridSize(const std::string pltfile);

    // ----------------------------
    // in Radial.cpp
    // ----------------------------

    // write radial plotfile to disk
    void WriteRadialFile(const BaseState<amrex::Real>& rho0_in,
                         const BaseState<amrex::Real>& p0_in,
                         const amrex::Vector<amrex::MultiFab>& u_in,
                         const amrex::Vector<amrex::MultiFab>& w0_in);

    // write additional diagnostics from model file
    void WriteModelDiagFile(const std::string& plotfilename);

    // read radial file from disk
    void ReadBaseCCFile(BaseState<amrex::Real>& r_in,
                        BaseState<amrex::Real>& rho0_in,
                        BaseState<amrex::Real>& p0_in);

    // read input model file from disk
    void ReadModelFile(const std::string& plotfilename,
                       BaseState<amrex::Real>& r_in,
                       BaseState<amrex::Real>& rho0_in,
                       BaseState<amrex::Real>& p0_in);

    // postprocessing subroutines
    void MakeRadialNFreq(const BaseState<amrex::Real>& p0_s,
                         const BaseState<amrex::Real>& rho0_s,
                         BaseState<amrex::Real>& freq0,
                         BaseState<amrex::Real>& entropy_s,
                         BaseState<amrex::Real>& grav_s);

    void MakeRadialGrav(const BaseState<amrex::Real>& rho0_in,
                        BaseState<amrex::Real>& grav,
                        const BaseStateArray<amrex::Real>& r_cc_loc);

    void MakeConvectionVel(const amrex::Vector<amrex::MultiFab>& velr,
                           BaseState<amrex::Real>& vel_conv);

    void MakeRadialRotationRatio(const BaseState<amrex::Real>& s0_in,
                                 const amrex::Vector<amrex::MultiFab>& omega,
                                 BaseState<amrex::Real>& ratio_omega);

    void MakeVelrc(const amrex::Vector<amrex::MultiFab>& vel,
                   const amrex::Vector<amrex::MultiFab>& w0rcart,
                   amrex::Vector<amrex::MultiFab>& rad_vel,
                   amrex::Vector<amrex::MultiFab>& circ_vel);

    void MakeRotationRate(const amrex::Vector<amrex::MultiFab>& vel,
                          const amrex::Vector<amrex::MultiFab>& w0cart,
                          amrex::Vector<amrex::MultiFab>& omega);

    void MakeLatShear(const amrex::Vector<amrex::MultiFab>& omega_in,
                      BaseState<amrex::Real>& shear,
                      const BaseStateArray<amrex::Real>& r_cc_loc);

    void MakeLatShearAvg(const amrex::Vector<amrex::MultiFab>& omega_in,
                         BaseState<amrex::Real>& shear);

    // ----------------------------
    // Slice.cpp functions
    // ----------------------------

    // write 2d slice plotfile to disk
    void Write2dSliceFile(const amrex::Vector<amrex::MultiFab>& rho_in,
                          const amrex::Vector<amrex::MultiFab>& p_in,
                          const amrex::Vector<amrex::MultiFab>& u_in,
                          const amrex::Vector<amrex::MultiFab>& w0_in,
                          const int deltat = 0, const int nfiles = 0);

    void PlotFileName(const int timestep, std::string* basefilename,
                      std::string* plotfilename);

    void MakeMeridionalCirculation(
        const amrex::Vector<amrex::MultiFab>& vel,
        const amrex::Vector<amrex::MultiFab>& w0rcart,
        amrex::Vector<amrex::MultiFab>& radii_vel,
        const amrex::Vector<amrex::Box>& domFine);

    void MakeBaroclinity(const amrex::Vector<amrex::MultiFab>& rho_s,
                         const amrex::Vector<amrex::MultiFab>& p_s,
                         amrex::Vector<amrex::MultiFab>& baroclinity,
                         const amrex::Vector<amrex::Box>& domFine);

    // ----------------------------
    // Average.cpp functions
    // ----------------------------

    // compute radially-averaged values
    void Average(const amrex::Vector<amrex::MultiFab>& phi,
                 BaseState<amrex::Real>& phibar, int comp);

    // compute averaged value onto r-theta (x-z) plane
    void Average2d(const amrex::Vector<amrex::MultiFab>& phi,
                   amrex::Vector<amrex::MultiFab>& phibar, int barcomp,
                   const amrex::Vector<amrex::Box>& bardomain, int comp);

    // ----------------------------
    // Fill3dData.cpp functions
    // ----------------------------

    void PutDataOnFaces(
        const amrex::Vector<amrex::MultiFab>& s_cc,
        amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM> >& face,
        const bool harmonic_avg);

    void AverageDownFaces(
        amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM> >& edge);

    // ----------------------------
    // test.cpp functions
    // ----------------------------

    // write plot file for exact solution test case
    void test();

    void WriteTestJobInfo(const std::string& dir, const std::string& base,
                          const int maxgridsize);

    /* private data members */

    // these are defined in `Postprocess.cpp`
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;

    // plotfile name
    std::string iFile;

    // input parameters
    std::string imFile;
    std::string deltat, numfiles;

    // geometry data
    int finest_level = 0;
    BaseStateGeometry base_geom;
    amrex::Vector<amrex::BoxArray> grid;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<amrex::Geometry> pgeom;
    amrex::GpuArray<amrex::Real, 3> center;
};
