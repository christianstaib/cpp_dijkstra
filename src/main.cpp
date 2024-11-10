// Your First C++ Program

#include <omp.h>

#include <glm/ext/scalar_constants.hpp>
#include <memory>
#include <random>

#include "cli.hpp"

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

void loging(std::vector<space::CelestialBody> &bodies, space::BodySystem &body_system, int &vis_step_size,
            double &time_step, std::ofstream &myfile, int &iteration, std::chrono::steady_clock::time_point *begin,
            indicators::ProgressBar *bar) {
  if (iteration % vis_step_size == 0) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double ms_per_it =
        ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - *begin).count() / iteration) / 1000.0;
    bar->set_option(indicators::option::PostfixText{std::to_string(ms_per_it) + "ms/it"});

    double ke = naive_calculations::get_kinetic_energy(body_system.num_bodies, body_system.mass, body_system.velocity);
    double pe =
        naive_calculations::get_potential_energy(body_system.num_bodies, body_system.mass, body_system.position);
    csv_writer::write_data(myfile, body_system.position, bodies);
  }
}

// x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
// Needs to be run in a parallel section
void update_positions(space::BodySystem &body_system, double time_step) {
  // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
#pragma omp parallel
  {
#pragma omp critical
    {
      body_system.min_edge = glm::dvec3(std::numeric_limits<double>::max());
      body_system.max_edge = glm::dvec3(std::numeric_limits<double>::min());
    }

    glm::dvec3 local_min_edge = glm::dvec3(std::numeric_limits<double>::max());
    glm::dvec3 local_max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp for simd schedule(guided)
    for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
      body_system.position[body_idx] +=
          body_system.velocity[body_idx] * time_step + 0.5 * body_system.acceleration[body_idx] * time_step * time_step;

      local_min_edge = min(local_min_edge, body_system.position[body_idx]);
      local_max_edge = max(local_max_edge, body_system.position[body_idx]);
    }

    // better use custom reduction?
#pragma omp critical
    {
      body_system.min_edge = min(body_system.min_edge, local_min_edge);
      body_system.max_edge = max(body_system.max_edge, local_max_edge);
    }
  }
}

void update_acceleration(space::BodySystem &body_system, octree::Octree &tree, double squared_theta) {
#pragma omp parallel for simd schedule(guided)
  for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
    body_system.acceleration_next_timestep[body_idx] =
        tree.get_acceleration(body_system.position[body_idx], squared_theta);
  }
}

// v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
void update_velocity(space::BodySystem &body_system, double step_size_days) {
#pragma omp parallel for simd
  for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
    body_system.velocity[body_idx] +=
        0.5 * (body_system.acceleration[body_idx] + body_system.acceleration_next_timestep[body_idx]) * step_size_days;
  }
}

void rebuild_tree(space::BodySystem &body_system, octree::Octree &tree) {
  glm::dvec3 center = (body_system.min_edge + body_system.max_edge) * 0.5;
  glm::dvec3 diff = body_system.max_edge - center;
  double size = std::max(std::max(diff.x, diff.y), diff.z);

  tree.clear(center, size);
  for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
    tree.insert(body_system.position[body_idx], body_system.mass[body_idx]);
  }
  tree.propagate();
}

int main(int argc, char **argv) {
  std::string bodies_file;
  int step_size_hours{1};
  int vis_step_size_hours{24};
  int t_end{1};
  double theta{1.0};

  auto app = cli::setup_app(&step_size_hours, &vis_step_size_hours, &t_end, &theta, &bodies_file);
  CLI11_PARSE(*app, argc, argv);

  //

  double squared_theta = theta * theta;
  double step_size_days = 1.0 / (24 * step_size_hours);
  size_t num_iterations = int((t_end * 365) / step_size_days);

  std::unique_ptr<indicators::ProgressBar> bar = cli::setup_progressbar(num_iterations);

  //

  std::vector<space::CelestialBody> bodies = space::read_bodies(bodies_file);
  // shfule to break order for better tree inserting performance
  auto rng = std::default_random_engine{};
  std::shuffle(std::begin(bodies), std::end(bodies), rng);

  // setup
  space::BodySystem body_system(bodies);

  std::ofstream myfile;
  myfile.open("data/data.txt");

  octree::Octree tree;
  rebuild_tree(body_system, tree);
  update_acceleration(body_system, tree, squared_theta);
  // naive_calculations::update_acceleration(body_system);
  std::swap(body_system.acceleration_next_timestep, body_system.acceleration);

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
    bar->tick();
    loging(bodies, body_system, vis_step_size_hours, step_size_days, myfile, iteration, &begin, bar.get());

    // TODO each MPI nodes updates its positions
    update_positions(body_system, step_size_days);
    // TODO each MPI nodes sends its positions to all other nodes via MPI_Allgather gg

    // // No need to send the tree, tree can be build on each node
    rebuild_tree(body_system, tree);
    update_acceleration(body_system, tree, squared_theta);
    // naive_calculations::update_acceleration(body_system);

    update_velocity(body_system, step_size_days);

    std::swap(body_system.acceleration_next_timestep, body_system.acceleration);
  }

  myfile.close();

  return 0;
}
