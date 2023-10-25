#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

class HelmholtzProjector {
    public:

        HelmholtzProjector(const hephaestus::InputParameters &params);
        ~HelmholtzProjector();

        void Project(hephaestus::GridFunctions &gridfunctions,
                                 const hephaestus::FESpaces &fespaces,
                                 hephaestus::BCMap &bc_map);

        void setForms();
        void setGrad();
        void setBCs();
        void solveLinearSystem();

    private:

      std::string hcurl_fespace_name_;
      std::string h1_fespace_name_;
      std::string gf_grad_name_;
      std::string gf_name_;

      mfem::ParFiniteElementSpace *H1FESpace_;
      mfem::ParFiniteElementSpace *HCurlFESpace_;
      mfem::ParGridFunction *q_;
      mfem::ParGridFunction *g; // H(Curl) projection of user specified source
      mfem::ParGridFunction *div_free_src_gf_; // Divergence free projected source

      mfem::ParLinearForm *gDiv_;
      mfem::ParBilinearForm *a0_;
      mfem::ParMixedBilinearForm *weakDiv_;
      mfem::ParDiscreteLinearOperator *grad_;

      mfem::HypreParMatrix *A0_;
      mfem::Vector *X0_, *B0_;
      
      mfem::Array<int> *ess_bdr_tdofs_;
      hephaestus::BCMap *bc_map_;
};


} // namespace hephaestus
