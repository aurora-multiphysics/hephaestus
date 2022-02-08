#include <memory>
#include <iostream>
#include <fstream>
#include "joule_solver.hpp"
#include "mesh_extras.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

class BCMap
{
    public:
    BCMap(const std::string & bc_name,
          Array<int> boundary_ids,
          int num_boundaries);

    std::string name;
    Array<int> markers;
};