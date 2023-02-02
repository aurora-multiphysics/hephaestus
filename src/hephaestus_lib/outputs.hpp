#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "../common/pfem_extras.hpp"
#include "mesh_extras.hpp"

namespace hephaestus {

class Outputs {
public:
  Outputs();
  Outputs(std::map<std::string, mfem::DataCollection *> data_collections_);

  std::map<std::string, mfem::DataCollection *> data_collections;
  std::map<std::string, mfem::socketstream *> socks_;

  void RegisterOutputFields(
      mfem::DataCollection *dc_, mfem::ParMesh *pmesh_,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables) {
    dc_->SetMesh(pmesh_);
    for (auto var = _variables.begin(); var != _variables.end(); ++var) {
      dc_->RegisterField(var->first, var->second);
    }
  }

  void WriteConsoleSummary(int id, double t, int it) {
    // Write a summary of the timestep to console.
    if (id == 0) {
      std::cout << std::fixed;
      std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
                << std::setprecision(3) << t << std::endl;
    }
  }

  void WriteOutputFields(mfem::DataCollection *dc_, double t, int it) {
    if (dc_) {
      dc_->SetCycle(it);
      dc_->SetTime(t);
      dc_->Save();
    }
  }

  void
  InitializeGLVis(int id,
                  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables) {
    if (id == 0) {
      std::cout << "Opening GLVis sockets." << std::endl;
    }

    for (auto var = _variables.begin(); var != _variables.end(); ++var) {
      socks_[var->first] = new mfem::socketstream;
      socks_[var->first]->precision(8);
    }

    if (id == 0) {
      std::cout << "GLVis sockets open." << std::endl;
    }
  }

  void DisplayToGLVis(mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables) {
    char vishost[] = "localhost";
    int visport = 19916;

    int Wx = 0, Wy = 0;                 // window position
    int Ww = 350, Wh = 350;             // window size
    int offx = Ww + 10, offy = Wh + 45; // window offsets

    for (auto var = _variables.begin(); var != _variables.end(); ++var) {
      mfem::common::VisualizeField(*socks_[var->first], vishost, visport,
                                   *(var->second), (var->first).c_str(), Wx, Wy,
                                   Ww, Wh);
      Wx += offx;
    }
  }
};

} // namespace hephaestus
