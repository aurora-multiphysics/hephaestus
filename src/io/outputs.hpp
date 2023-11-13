#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "../common/pfem_extras.hpp"
#include "gridfunctions.hpp"
#include "mesh_extras.hpp"

namespace hephaestus {

class Outputs : public mfem::NamedFieldsMap<mfem::DataCollection> {
public:
  Outputs();
  Outputs(hephaestus::GridFunctions &gridfunctions);

  void SetGridFunctions(hephaestus::GridFunctions &gridfunctions) {
    _gridfunctions = &gridfunctions;
  }

  void EnableGLVis(bool &use_glvis) { _use_glvis = use_glvis; }

  void Reset() {
    // Reset cycle counter
    _cycle = 0;
    // Set up DataCollections to track fields of interest.
    RegisterOutputFields();
    // Write initial fields to disk
    WriteOutputFields(0.0);

    // Initialize GLVis _use_glvis and send the initial condition
    // by socket to a GLVis server.
    if (_use_glvis) {
      InitializeGLVis(_my_rank);
      DisplayToGLVis();
    }
  }

  void Write(double t = 1.0) {
    _cycle++;
    // Output timestep summary to console
    WriteConsoleSummary(_my_rank, t);
    // Make sure all ranks have sent their 'v' solution before initiating
    // another set of GLVis connections (one from each rank):
    // Send output fields to GLVis for visualisation
    MPI_Barrier(_my_comm);
    if (_use_glvis) {
      DisplayToGLVis();
    }
    // Save output fields at timestep to DataCollections
    WriteOutputFields(t);
  }

private:
  std::map<std::string, mfem::socketstream *> socks_;

  hephaestus::GridFunctions *_gridfunctions;
  int _cycle;
  bool _use_glvis;
  MPI_Comm _my_comm;
  int _n_ranks, _my_rank;

  void RegisterOutputFields() {
    for (auto output = begin(); output != end(); ++output) {
      auto const &dc_(output->second);
      mfem::ParMesh *pmesh_(
          _gridfunctions->begin()->second->ParFESpace()->GetParMesh());
      dc_->SetMesh(pmesh_);
      for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
           ++var) {
        dc_->RegisterField(var->first, var->second);
      }
    }
  }

  void WriteConsoleSummary(int _my_rank, double t) {
    // Write a summary of the timestep to console.
    if (_my_rank == 0) {
      std::cout << std::fixed;
      std::cout << "step " << std::setw(6) << _cycle
                << ",\tt = " << std::setw(6) << std::setprecision(3) << t
                << std::endl;
    }
  }

  void WriteOutputFields(double t) {
    // Write fields to disk
    for (auto output = begin(); output != end(); ++output) {
      auto const &dc_(output->second);
      dc_->SetCycle(_cycle);
      dc_->SetTime(t);
      dc_->Save();
    }
  }

  void InitializeGLVis(int _my_rank) {
    if (_my_rank == 0) {
      std::cout << "Opening GLVis sockets." << std::endl;
    }

    for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
         ++var) {
      socks_[var->first] = new mfem::socketstream;
      socks_[var->first]->precision(8);
    }

    if (_my_rank == 0) {
      std::cout << "GLVis sockets open." << std::endl;
    }
  }

  void DisplayToGLVis() {
    char vishost[] = "localhost";
    int visport = 19916;

    int Wx = 0, Wy = 0;                 // window position
    int Ww = 350, Wh = 350;             // window size
    int offx = Ww + 10, offy = Wh + 45; // window offsets

    for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
         ++var) {
      mfem::common::VisualizeField(*socks_[var->first], vishost, visport,
                                   *(var->second), (var->first).c_str(), Wx, Wy,
                                   Ww, Wh);
      Wx += offx;
    }
  }
};

} // namespace hephaestus
