#pragma once
#include "essential_bcs.hpp"
#include "integrated_bcs.hpp"
#include "named_fields_map.hpp"
#include "robin_bcs.hpp"

namespace hephaestus
{

class BCMap : public hephaestus::NamedFieldsMap<hephaestus::BoundaryCondition>
{
public:
  mfem::Array<int> GetEssentialBdrMarkers(const std::string & name_, mfem::Mesh * mesh_);

  void ApplyEssentialBCs(const std::string & name_,
                         mfem::Array<int> & ess_tdof_list,
                         mfem::GridFunction & gridfunc,
                         mfem::Mesh * mesh_);

  void ApplyEssentialBCs(const std::string & name_,
                         mfem::Array<int> & ess_tdof_list,
                         mfem::ParComplexGridFunction & gridfunc,
                         mfem::Mesh * mesh_);

  void ApplyIntegratedBCs(const std::string & name_, mfem::LinearForm & lf, mfem::Mesh * mesh_);

  void ApplyIntegratedBCs(const std::string & name_,
                          mfem::ParComplexLinearForm & clf,
                          mfem::Mesh * mesh_);

  void ApplyIntegratedBCs(const std::string & name_,
                          mfem::ParSesquilinearForm & clf,
                          mfem::Mesh * mesh_);

  /// @brief Updates boundary conditions after a mesh refinement. This will
  /// iterate over all stored boundary conditions in the container and call
  /// their Update methods which will update the boundary attributes since each
  /// boundary element has been split into additional elements. The Apply BCs
  /// methods should then be called.
  void Update(mfem::Mesh & mesh);
};

} // namespace hephaestus
