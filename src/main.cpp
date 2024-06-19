#include <nlohmann/json.hpp>
#include <stdio.h>

#include "test.hpp"

int main() {

#pragma omp parallel
  { printf("Hello World\n"); }
  return 0;
}
