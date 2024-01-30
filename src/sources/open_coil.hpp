#pragma once
#include "helmholtz_projector.hpp"
#include "scalar_potential_source.hpp"
#include "scaled_vector_gridfunction_aux.hpp"
#include "source_base.hpp"

namespace hephaestus
{

double highV(const mfem::Vector & x, double t);
double lowV(const mfem::Vector & x, double t);

double calcFlux(mfem::GridFunction * v_field, int face_attr);

double calcFlux(mfem::GridFunction * v_field, int face_attr, mfem::Coefficient & q);

void inheritBdrAttributes(const mfem::ParMesh * parent_mesh, mfem::ParSubMesh * child_mesh);

// Applies the HelmholtzProjector onto the J GridFunction to clean it of any
// divergences. This is for the simplest case with no BCs
void cleanDivergence(mfem::ParGridFunction & Vec_GF, hephaestus::InputParameters solve_pars);

// The more complicated case where BCs are needed
void cleanDivergence(const hephaestus::GridFunctions & gfs,
                     const hephaestus::BCMap & bcs,
                     const std::string vec_gf_name,
                     const std::string scalar_gf_name,
                     hephaestus::InputParameters solve_pars);

void attrToMarker(const mfem::Array<int> attr_list, mfem::Array<int> & marker_list, int max_attr);

class OpenCoilSolver : public hephaestus::Source
{

public:
  OpenCoilSolver(const hephaestus::InputParameters & params,
                 mfem::Array<int> coil_dom,
                 const std::pair<int, int> electrodes);

  ~OpenCoilSolver() override;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * lf) override;
  void SubtractSource(mfem::ParGridFunction * gf) override;

  // Initialises the child submesh.
  void InitChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void MakeFESpaces();

  // Creates the relevant GridFunctions for the child submesh.
  void MakeGridFunctions();

  // Sets up the boundary conditions to be used in the ScalarPotentialSource.
  // calculation.
  void SetBCs();

  // Solves for the divergence-free Hodge dual of the electric current based on
  // Dirichlet BCs.
  void SPSCurrent();

  // Creates a mass matrix with basis functions that will be used in the Apply()
  // method
  void BuildM1();

  // Sets the boundary attribute for the face to be used as reference in flux
  // calculation
  void SetRefFace(const int face);

private:
  // Parameters
  int _order_h1;
  int _order_hcurl;
  int _ref_face;
  bool _grad_phi_transfer;

  std::pair<int, int> _elec_attrs;
  mfem::Array<int> _coil_domains;
  mfem::Array<int> _coil_markers;
  hephaestus::InputParameters _solver_options;

  mfem::Coefficient * _sigma{nullptr};
  mfem::Coefficient * _itotal{nullptr};

  bool _owns_sigma{false};
  bool _owns_itotal{false};

  // Names
  std::string _grad_phi_name;
  std::string _v_gf_name;
  std::string _i_coef_name;
  std::string _cond_coef_name;
  std::string _source_efield_gf_name;
  std::string _source_jfield_gf_name;

  // Parent mesh, FE space, and current
  mfem::ParMesh * _mesh_parent{nullptr};

  mfem::ParGridFunction * _grad_phi_parent{nullptr};
  std::unique_ptr<mfem::ParGridFunction> _grad_phi_t_parent{nullptr};

  mfem::ParGridFunction * _v_parent{nullptr};
  std::unique_ptr<mfem::ParGridFunction> _vt_parent{nullptr};

  mfem::ParGridFunction * _source_electric_field{nullptr};
  mfem::ParGridFunction * _source_current_density{nullptr};

  // Child mesh and FE spaces
  std::unique_ptr<mfem::ParSubMesh> _mesh{nullptr};

  std::unique_ptr<mfem::ParFiniteElementSpace> _h1_fe_space{nullptr};
  std::unique_ptr<mfem::H1_FECollection> _h1_fe_space_fec{nullptr};

  std::unique_ptr<mfem::ParFiniteElementSpace> _h_curl_fe_space{nullptr};
  std::unique_ptr<mfem::ND_FECollection> _h_curl_fe_space_fec{nullptr};

  // Child GridFunctions
  std::unique_ptr<mfem::ParGridFunction> _grad_phi{nullptr};
  std::unique_ptr<mfem::ParGridFunction> _v{nullptr};

  // Child boundary condition objects
  mfem::FunctionCoefficient _high_src;
  mfem::FunctionCoefficient _low_src;

  mfem::Array<int> _high_terminal;
  mfem::Array<int> _low_terminal;

  // Mass Matrix
  std::unique_ptr<mfem::ParBilinearForm> _m1{nullptr};

  // BC Map
  hephaestus::BCMap _bc_maps;

  // Final LinearForm
  std::unique_ptr<mfem::ParLinearForm> _final_lf{nullptr};
};

} // namespace hephaestus