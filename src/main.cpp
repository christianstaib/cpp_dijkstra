// Your First C++ Program

#include <omp.h>

#include <glm/ext/scalar_constants.hpp>
#include <memory>
#include <random>

#define GLM_ENABLE_EXPERIMENTAL

#include <CLI/CLI.hpp>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <glm/detail/qualifier.hpp>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/quaternion.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <unordered_map>
#include <vector>

#include "csv_writer.hpp"
#include "naive_calculations.hpp"
#include "octree.hpp"
#include "space.hpp"

void loging(std::vector<space::CelestialBody> &bodies, size_t &num_bodies, double *masses, glm::dvec3 *positions,
            glm::dvec3 *velocities, int &vis_step_size, double &time_step, std::ofstream &myfile, int &iteration,
            std::chrono::steady_clock::time_point *begin, indicators::ProgressBar *bar) {
  if (iteration % vis_step_size == 0) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double ms_per_it =
        ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - *begin).count() / iteration) / 1000.0;
    bar->set_option(indicators::option::PostfixText{std::to_string(ms_per_it) + "ms/it"});

    double ke = naive_calculations::get_kinetic_energy(num_bodies, masses, velocities);
    double pe = naive_calculations::get_potential_energy(num_bodies, masses, positions);
    csv_writer::write_data(myfile, positions, bodies);
  }
}

// x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
// Needs to be run in a parallel section
void update_positions(size_t num_bodies, glm::dvec3 *velocities, glm::dvec3 *forces, glm::dvec3 *positions,
                      double time_step, glm::dvec3 *min_edge, glm::dvec3 *max_edge) {
#pragma omp critical
  {
    *min_edge = glm::dvec3(std::numeric_limits<double>::max());
    *max_edge = glm::dvec3(std::numeric_limits<double>::min());
  }

  glm::dvec3 local_min_edge = glm::dvec3(std::numeric_limits<double>::max());
  glm::dvec3 local_max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp for simd schedule(guided)
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    positions[body_idx] += velocities[body_idx] * time_step + 0.5 * forces[body_idx] * time_step * time_step;

    local_min_edge = min(local_min_edge, positions[body_idx]);
    local_max_edge = max(local_max_edge, positions[body_idx]);
  }

  // better use custom reduction?
#pragma omp critical
  {
    *min_edge = min(*min_edge, local_min_edge);
    *max_edge = max(*max_edge, local_max_edge);
  }
#pragma omp barrier
}

void rebuild_tree(size_t num_bodies, glm::dvec3 *positions, double *masses, glm::dvec3 *min_edge, glm::dvec3 *max_edge,
                  octree::Octree &tree) {
  glm::dvec3 center = (*min_edge + *max_edge) * 0.5;
  glm::dvec3 diff = *max_edge - center;
  double size = std::max(std::max(diff.x, diff.y), diff.z);

  tree.clear(center, size);
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    tree.insert(positions[body_idx], masses[body_idx]);
  }
  tree.propagate();
}

std::unique_ptr<CLI::App> setup_app(int *step_size_hours, int *vis_step_size_hours, int *t_end, double *theta,
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

std::unique_ptr<indicators::ProgressBar> setup_progressbar(size_t num_iterations) {
  using namespace indicators;
  show_console_cursor(true);

  std::unique_ptr<indicators::ProgressBar> bar = std::make_unique<indicators::ProgressBar>(
      option::BarWidth{50}, option::PrefixText{"Simulation"}, option::ShowElapsedTime{true},
      option::ShowRemainingTime{true}, indicators::option::MaxProgress{num_iterations});

  return bar;
}

int main(int argc, char **argv) {
  int step_size_hours{1};
  int vis_step_size_hours{24};
  int t_end{1};
  double theta{1.0};
  std::string bodies_file;

  auto app = setup_app(&step_size_hours, &vis_step_size_hours, &t_end, &theta, &bodies_file);
  CLI11_PARSE(*app, argc, argv);

  std::cout << "Working on file: " << bodies_file << "\n";
  std::cout << "Step size in hours: " << step_size_hours << "\n";
  std::cout << "Vis step size in hours: " << vis_step_size_hours << "\n";
  std::cout << "Length of simulation in year: " << t_end << "\n";
  std::cout << "Barnes-Hut theta: " << theta << "\n";
  std::cout << "\n";

  //

  double squared_theta = theta * theta;
  double step_size_days = 1.0 / (24 * step_size_hours);
  size_t num_iterations = int((t_end * 365) / step_size_days);

  // Hide cursor
  std::unique_ptr<indicators::ProgressBar> bar = setup_progressbar(num_iterations);

  //

  std::vector<space::CelestialBody> bodies = space::read_bodies(bodies_file);
  // shfule to break order for better tree inserting performance
  auto rng = std::default_random_engine{};
  std::shuffle(std::begin(bodies), std::end(bodies), rng);

  // setup
  size_t num_bodies = bodies.size();
  double *masses = new double[num_bodies];
  glm::dvec3 *position = new glm::dvec3[num_bodies];
  glm::dvec3 *velocity = new glm::dvec3[num_bodies];
  glm::dvec3 *old_force = new glm::dvec3[num_bodies];

  glm::dvec3 min_edge(std::numeric_limits<double>::max());
  glm::dvec3 max_edge(std::numeric_limits<double>::min());

  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    masses[body_idx] = bodies[body_idx].mass;
    position[body_idx] = bodies[body_idx].pos;
    velocity[body_idx] = bodies[body_idx].vel;
    min_edge = min(min_edge, position[body_idx]);
    max_edge = max(max_edge, position[body_idx]);
  }

  octree::Octree test(glm::dvec3(0.0), 0.0);
  rebuild_tree(num_bodies, position, masses, &min_edge, &max_edge, test);

  std::ofstream myfile;
  myfile.open("data/data.txt");

#pragma omp parallel for schedule(static)
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    old_force[body_idx] = test.get_force(position[body_idx], theta);
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    bar->tick();
    loging(bodies, num_bodies, masses, position, velocity, vis_step_size_hours, step_size_days, myfile, iteration,
           &begin, bar.get());

#pragma omp parallel
    {
      update_positions(num_bodies, velocity, old_force, position, step_size_days, &min_edge, &max_edge);

#pragma omp single
      rebuild_tree(num_bodies, position, masses, &min_edge, &max_edge, test);

      // get_force performs a tree traversal with variable execution times,
      // therfore use schedule(guided) workload balancing
#pragma omp for simd schedule(guided)
      for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
        glm::dvec3 new_force = test.get_force(position[body_idx], squared_theta);
        // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
        velocity[body_idx] += 0.5 * (old_force[body_idx] + new_force) * step_size_days;
        old_force[body_idx] = new_force;
      }
    }
  }

  myfile.close();

  free(masses);
  free(position);
  free(velocity);
  free(old_force);

  return 0;
}
