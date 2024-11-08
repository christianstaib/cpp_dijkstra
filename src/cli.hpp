#pragma once

#include <CLI/CLI.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>

#include "CLI/CLI.hpp"
namespace cli {

std::unique_ptr<CLI::App> setup_app(int *step_size_hours, int *vis_step_size_hours, int *t_end, double *theta,
                                    std::string *bodies_file);

std::unique_ptr<indicators::ProgressBar> setup_progressbar(size_t num_iterations);
}  // namespace cli
