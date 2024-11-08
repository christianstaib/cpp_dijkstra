#include "cli.hpp"

std::unique_ptr<CLI::App> cli::setup_app(int *step_size_hours, int *vis_step_size_hours, int *t_end, double *theta,
                                         std::string *bodies_file) {
  std::unique_ptr<CLI::App> app = std::make_unique<CLI::App>("Gravity Simulator");

  app->set_version_flag("--version", std::string(CLI11_VERSION));

  CLI::Option *opt0 = app->add_option("--file", *bodies_file, "File name");
  opt0->required();

  CLI::Option *opt1 = app->add_option("--dt", *step_size_hours, "Step size in hours")->capture_default_str();
  CLI::Option *opt4 =
      app->add_option("--vs", *vis_step_size_hours, "Visualization step size in hours")->capture_default_str();

  CLI::Option *opt2 = app->add_option("--t_end", *t_end, "Length of simulation in years")->capture_default_str();
  CLI::Option *opt3 = app->add_option("--theta", *theta, "Barnes-Hut theta")->capture_default_str();

  return app;
}

std::unique_ptr<indicators::ProgressBar> cli::setup_progressbar(size_t num_iterations) {
  using namespace indicators;
  show_console_cursor(true);

  std::unique_ptr<indicators::ProgressBar> bar = std::make_unique<indicators::ProgressBar>(
      option::BarWidth{50}, option::PrefixText{"Simulation"}, option::ShowElapsedTime{true},
      option::ShowRemainingTime{true}, indicators::option::MaxProgress{num_iterations});

  return bar;
}
