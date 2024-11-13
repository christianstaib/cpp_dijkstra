// Your First C++ Program

#include <mpi.h>
#include <omp.h>

#include <glm/ext/scalar_constants.hpp>
#include <memory>
#include <random>
#include <valarray>

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
            double &time_step, std::ofstream &myfile, int &iteration, size_t &num_iterations,
            std::chrono::steady_clock::time_point *begin, indicators::ProgressBar *bar) {
  if (iteration % vis_step_size == 0) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double ms_per_it =
        ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - *begin).count() / iteration) / 1000.0;
    bar->set_option(indicators::option::PostfixText{std::to_string(ms_per_it) + "ms/it"});

    printf("finished %f %%, %f ms per it\n", 100.0 * (float)iteration / (float)num_iterations, ms_per_it);

    double ke = naive_calculations::get_kinetic_energy(body_system.num_bodies, body_system.mass, body_system.velocity);
    double pe =
        naive_calculations::get_potential_energy(body_system.num_bodies, body_system.mass, body_system.position);
    csv_writer::write_data(myfile, body_system.position, bodies);
  }
}

void update_acceleration(space::BodySystem &body_system, octree::Octree &tree, double squared_theta) {
#pragma omp parallel for simd schedule(guided)
  for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
    body_system.acceleration_next_timestep[body_idx] =
        tree.get_acceleration(body_system.position[body_idx], squared_theta);
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

  //
  //
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);

  std::unique_ptr<indicators::ProgressBar> bar = cli::setup_progressbar(num_iterations);
  //

  std::vector<space::CelestialBody> bodies = space::read_bodies(bodies_file);
  // shfule to break order for better tree inserting performance
  auto rng = std::default_random_engine{};
  std::shuffle(std::begin(bodies), std::end(bodies), rng);

  printf("len of bodies is %zu\n", bodies.size());

  size_t chunk_size = bodies.size() / world_size;

  // setup

  std::vector<space::CelestialBody> local_bodies(bodies.begin() + world_rank * chunk_size,
                                                 bodies.begin() + (world_rank + 1) * chunk_size);
  space::BodySystem local_body_system(local_bodies);
  space::BodySystem global_body_system(bodies);
  glm::dvec3 *global_positions = new glm::dvec3[bodies.size()];

  std::ofstream myfile;
  myfile.open("data/data.txt");

  octree::Octree tree;
  rebuild_tree(global_body_system, tree);
  update_acceleration(local_body_system, tree, squared_theta);
  // naive_calculations::update_acceleration(body_system);
  std::swap(local_body_system.acceleration_next_timestep, local_body_system.acceleration);

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
    // TODO bar->tick();
    if (world_rank == 0) {
      loging(bodies, global_body_system, vis_step_size_hours, step_size_days, myfile, iteration, num_iterations, &begin,
             bar.get());
    }

    // TODO each MPI nodes updates its positions
    naive_calculations::update_positions(local_body_system, step_size_days);
    // TODO each MPI nodes sends its positions to all other nodes via MPI_Allgather gg

    MPI_Allgather(local_body_system.position, chunk_size * 3, MPI_DOUBLE, global_body_system.position, chunk_size * 3,
                  MPI_DOUBLE, MPI_COMM_WORLD);

    // // No need to send the tree, tree can be build on each node
    rebuild_tree(global_body_system, tree);
    update_acceleration(local_body_system, tree, squared_theta);
    // naive_calculations::update_acceleration(local_body_system, global_body_system);

    naive_calculations::update_velocity(local_body_system, step_size_days);

    std::swap(local_body_system.acceleration_next_timestep, local_body_system.acceleration);
  }

  myfile.close();

  return 0;
}
