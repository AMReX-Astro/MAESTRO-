#include "Maestro.H"

#ifdef DO_PROBLEM_POST_TIMESTEP

using namespace amrex;

void Maestro::ProblemPostTimestep() {
    static const std::string script_name = "post_to_slack.sh";

    if ((plot_int > 0 && istep % plot_int == 0) ||
        (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
        (istep == max_step) || (t_old >= stop_time)) {
        // wrote a plotfile
        std::string plotfilename = plot_base_name;
        PlotFileName(istep, &plotfilename);
        // std::string command = "bash " + script_name + " " + std::to_string(t_new) + " " + std::to_string(istep) + " " + plotfilename + " &";

        // std::system(command.c_str());

        Print() << "Output file " << plotfilename << " at time "
                << std::to_string(t_new) << " and step "
                << std::to_string(istep) << std::endl;

    } else if ((small_plot_int > 0 && istep % small_plot_int == 0) ||
               (small_plot_deltat > 0 &&
                std::fmod(t_new, small_plot_deltat) < dt) ||
               (istep == max_step) || (t_old >= stop_time)) {
        // wrote a small plotfile
        std::string plotfilename = small_plot_base_name;
        PlotFileName(istep, &plotfilename);

        // std::string command = "bash " + script_name + " " + std::to_string(t_new) + " " + std::to_string(istep) + " " + plotfilename + " &";
        //
        // std::system(command.c_str());

        Print() << "Output file " << plotfilename << " at time "
                << std::to_string(t_new) << " and step "
                << std::to_string(istep) << std::endl;
    }
}

#endif
