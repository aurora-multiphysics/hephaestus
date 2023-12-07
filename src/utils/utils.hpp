#pragma once
#include "../common/pfem_extras.hpp"
#include "coefficients.hpp"
#include "helmholtz_projector.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"

namespace hephaestus {

// Useful functions available to all classes

double calcFlux(mfem::GridFunction *v_field, int face_attr);

void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd,
                      mfem::Array<int> &arr);

void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr);

//template <typename T> void ifDelete(T *ptr);

void inheritBdrAttributes(const mfem::ParMesh *parent_mesh,
                          mfem::ParSubMesh *child_mesh);

void attrToMarker(const mfem::Array<int> attr_list,
                  mfem::Array<int> &marker_list, int max_attr);

void cleanDivergence(mfem::ParGridFunction &Vec_GF,
                     hephaestus::InputParameters solve_pars);

void cleanDivergence(hephaestus::GridFunctions &gfs, hephaestus::BCMap &bcs,
                     const std::string vec_gf_name,
                     const std::string scalar_gf_name,
                     hephaestus::InputParameters solve_pars);

} // namespace hephaestus