// Your First C++ Program

#include <omp.h>

#include <glm/ext/scalar_constants.hpp>
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

#include "constants.hpp"
#include "naive_calculations.hpp"
#include "octree.hpp"
#include "space.hpp"

std::vector<space::CelestialBody> read_bodies(std::string path) {
  // Set up data structures needed for reading the asteroids.
  space::CelestialBody sun = space::CelestialBody::sun();

  std::unordered_map<std::string, space::CelestialBody> name_to_body;
  name_to_body.insert({"Sun", sun});

  std::vector<std::pair<glm::dvec3, space::CelestialBody>> pos_to_body;
  pos_to_body.push_back({sun.pos, sun});

  std::vector<space::CelestialBody> bodies;
  bodies.push_back(sun);

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line);  // skip header

    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_asteroid(line);
      space::CelestialBody body = row.to_body(name_to_body.size(), name_to_body);

      if (body.name == "Earth") {
        printf("Earth %f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
      }

      if (!body.name.empty()) {
        if (name_to_body.find(body.name) == name_to_body.end()) {
          name_to_body.insert({body.name, body});
        } else {
          printf(
              "Error: A bdoy with the name %s is already known (distance "
              "%f km)\n",
              body.name.c_str(),
              glm::distance(body.pos, name_to_body.at(body.name).pos) * constants::meters_per_astronomical_unit /
                  1000.0);
        }
      }

      for (const auto &entry : pos_to_body) {
        if (glm::distance(entry.second.pos, body.pos) <= 10000 * constants::astronomical_units_per_meter) {
          printf(
              "Error: A body at the position %f %f %f is already known (it "
              "is called %s)\n",
              body.pos.x, body.pos.y, body.pos.z, entry.second.name.c_str());
        }
      }
      pos_to_body.push_back({body.pos, body});

      bodies.push_back(body);
    }
    file.close();
  }

  return bodies;
}



void write_data(std::ofstream &myfile, glm::dvec3 *positions, std::vector<space::CelestialBody> const &bodies) {
  for (size_t i = 0; i < bodies.size(); ++i) {
    myfile << "(" << bodies[i].name.c_str() << "," << bodies[i].type.c_str() << "," << positions[i].x << ","
           << positions[i].y << ")";
    if (i != bodies.size() - 1) {
      myfile << ",";
    } else {
      myfile << "\n";
    }
  }
}

void loging(std::vector<space::CelestialBody> &bodies, size_t &num_bodies, double *masses, glm::dvec3 *positions,
            glm::dvec3 *velocities, int &vis_step_size, double &time_step, std::ofstream &myfile, int &iteration,
            std::chrono::steady_clock::time_point *begin, indicators::ProgressBar &bar) {
  if (iteration % vis_step_size == 0) {
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double ms_per_it =
        ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - *begin).count() / iteration) / 1000.0;
    bar.set_option(indicators::option::PostfixText{std::to_string(ms_per_it) + "ms/it"});

    double ke = naive_calculations::get_kinetic_energy(num_bodies, masses, velocities);
    double pe = naive_calculations::get_potential_energy(num_bodies, masses, positions);
    write_data(myfile, positions, bodies);
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

#pragma omp for simd schedule(static)
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
                  octree::Octree &test) {
  glm::dvec3 center = (*min_edge + *max_edge) * 0.5;
  glm::dvec3 diff = *max_edge - center;
  double size = std::max(std::max(diff.x, diff.y), diff.z);

  test.clear(center, size);
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    test.insert(positions[body_idx], masses[body_idx]);
  }
  test.propagate();
}

/// Naive approach to get gravitional force on all bodies in (Kg*AU)/d^2.
glm::dvec3 get_gravitational_force_ld(size_t body_idx_want_force, size_t num_bodies, glm::dvec3 *positions,
                                      double *masses) {
  long double x = 0;
  long double y = 0;
  long double z = 0;

  glm::dvec3 distance_vector;
  double squared_distance;

  for (size_t j = 0; j < num_bodies; ++j) {
    // Check body_idx_want_force == j can be skiped as distance_vector will be
    // zero in this case

    // Precompute distance vector
    distance_vector = positions[j] - positions[body_idx_want_force];
    squared_distance = glm::length2(distance_vector) + constants::squared_softening_factor;
    // x*sqrt(x) should be faster than pow(x, 3/2)
    x += (masses[j] * distance_vector.x) / (squared_distance * sqrt(squared_distance));
    y += (masses[j] * distance_vector.y) / (squared_distance * sqrt(squared_distance));
    z += (masses[j] * distance_vector.z) / (squared_distance * sqrt(squared_distance));
  }

  x *= constants::gravitational_constant_in_au3_per_kg_d2;
  y *= constants::gravitational_constant_in_au3_per_kg_d2;
  z *= constants::gravitational_constant_in_au3_per_kg_d2;

  glm::dvec3 force(x, y, z);

  // multipling once at the end is faster and also better for precision
  return force;
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
  int num_iterations = int((t_end * 365) / step_size_days);

  // Hide cursor
  using namespace indicators;
  show_console_cursor(true);

  indicators::ProgressBar bar{option::BarWidth{50}, option::PrefixText{"Simulation"}, option::ShowElapsedTime{true},
                              option::ShowRemainingTime{true}, indicators::option::MaxProgress{num_iterations}};

  //

  std::vector<space::CelestialBody> bodies = read_bodies(bodies_file);
  // shfule to break order for better tree inserting performance
  auto rng = std::default_random_engine{};
  std::shuffle(std::begin(bodies), std::end(bodies), rng);
  size_t num_bodies = bodies.size();
  printf("there are %zu bodies\n", num_bodies);

  // setup
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

  std::ofstream myfile;
  myfile.open("data/data.txt");

  octree::Octree test(glm::dvec3(0.0), 0.0);
  rebuild_tree(num_bodies, position, masses, &min_edge, &max_edge, test);
#pragma omp parallel for schedule(static)
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    old_force[body_idx] = test.get_force(position[body_idx], theta);
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  int n = 1000;
  for (int i = 0; i < n; ++i) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    octree::VecTreeNode root = octree::VecTreeNode::create_root(num_bodies, position, masses, test.nodes[0].cube);
    root.split_all();
    if (i == 0) {
      printf("root size is %f\n", test.nodes[0].cube.half_edge_length);
    }
    root.free_children();
  }
  // printf("root size is %f\n", root.cube.half_edge_length);
  // root.split();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double ms_per_it = ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000.0;
  printf("building tree took %f ms\n", ms_per_it / n);

  std::vector<double> build_tree_time;

  begin = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    bar.tick();
    loging(bodies, num_bodies, masses, position, velocity, vis_step_size_hours, step_size_days, myfile, iteration,
           &begin, bar);

#pragma omp parallel
    {
      update_positions(num_bodies, velocity, old_force, position, step_size_days, &min_edge, &max_edge);

#pragma omp single
      {
        std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
        rebuild_tree(num_bodies, position, masses, &min_edge, &max_edge, test);
        std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
        double xxx = ((double)std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count()) / 1000.0;
        build_tree_time.push_back(xxx);
      }

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

  double sum = accumulate(build_tree_time.begin(), build_tree_time.end(), 0.0);
  printf("build tree time is %f ms\n", sum / build_tree_time.size());

  myfile.close();

  free(masses);
  free(position);
  free(velocity);
  free(old_force);

  return 0;
}
