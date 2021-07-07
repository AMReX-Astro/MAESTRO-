#include <Postprocess.H>
#include <regex>

using namespace amrex;

// ------------------------------------------
// Read input parameters from job info file
// ------------------------------------------

// Gets the variable ``varname`` from the ``job_info`` file and returns as a
// string
std::string Postprocess::GetVarFromJobInfo(const std::string pltfile,
                                           const std::string& varname) {
    std::string filename = pltfile + "/job_info";
    std::regex re("(?:[ \\t]*)" + varname + "\\s*=\\s*(.*)\\s*\\n");
    //std::regex re(varname + " = ([^ ]*)\\n");

    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
        std::stringstream buf;
        buf << jobfile.rdbuf();
        std::string file_contents = buf.str();

        if (std::regex_search(file_contents, m, re)) {
            return m[1];
        } else {
            Print() << "Unable to find " << varname << " in job_info file!"
                    << std::endl;
        }
    } else {
        Print() << "Could not open job_info file!" << std::endl;
    }

    return "";
}

// Get drdxfac from the job info file
int Postprocess::GetdrdxFac(const std::string pltfile) {
    auto drdxfac_str = GetVarFromJobInfo(pltfile, "maestro.drdxfac");
    // Print() << "drdxfac_str = " << drdxfac_str << std::endl;

    if (drdxfac_str == "") {
        return 5;
    }

    // retrieve first number
    std::istringstream iss{drdxfac_str};
    int drdxfac;
    std::string s;
    if (std::getline(iss, s, ' ')) drdxfac = stoi(s);

    return drdxfac;
}

// Get octant from the job info file
bool Postprocess::GetOctant(const std::string pltfile) {
    auto octant_str = GetVarFromJobInfo(pltfile, "maestro.octant");

    if (octant_str == "true")
        return true;
    else
        // Print() << "octant_str = " << octant_str << std::endl;
        return false;
}

// Get gamma from the job info file
Real Postprocess::GetGamma(const std::string pltfile) {
    auto gamma_str = GetVarFromJobInfo(pltfile, "eos_gamma");
    // Print() << "gamma_str = " << gamma_str << std::endl;

    // retrieve first number
    std::istringstream iss{gamma_str};
    Real gamma;
    std::string s;
    if (std::getline(iss, s, ' ')) gamma = stod(s);

    return gamma;
}

// Get rotation frequency from the job info file
Real Postprocess::GetRotationFreq(const std::string pltfile) {
    auto rotationfreq_str =
        GetVarFromJobInfo(pltfile, "maestro.rotational_frequency");
    // Print() << "rotationfreq_str = " << rotationfreq_str << std::endl;

    // retrieve first number
    std::istringstream iss{rotationfreq_str};
    Real freq0;
    std::string s;
    if (std::getline(iss, s, ' ')) freq0 = stod(s);

    return freq0;
}

// Get time step of plot file
int Postprocess::GetTimeStep(const std::string plotfilename,
                             const std::string basefilename) {
    int base = basefilename.length();
    std::string timestep_str = plotfilename.substr(base, 7);

    return stoi(timestep_str);
}

// Get max grid size from the job info file
int Postprocess::GetMaxGridSize(const std::string pltfile) {
    auto maxgridsize_str = GetVarFromJobInfo(pltfile, "amr.max_grid_size");
    // Print() << "maxgridsize_str = " << maxgridsize_str << std::endl;

    // retrieve first number
    std::istringstream iss{maxgridsize_str};
    int max_grid_size;
    std::string s;
    if (std::getline(iss, s, ' ')) max_grid_size = stoi(s);

    return max_grid_size;
}
