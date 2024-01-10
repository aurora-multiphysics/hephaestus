#pragma once
#include "coefficient_aux.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus {

// The VectorGridFunctionDotProductCoefficient evaluates the dot product of two
// (possibly complex) gridfunctions
class VectorGridFunctionDotProductCoefficient : public mfem::Coefficient {
private:
  const mfem::ParGridFunction *_u_gf_re;
  const mfem::ParGridFunction *_u_gf_im;
  const mfem::ParGridFunction *_v_gf_re;
  const mfem::ParGridFunction *_v_gf_im;
  mfem::Coefficient &_coef;

public:
  VectorGridFunctionDotProductCoefficient(
      mfem::Coefficient &coef, const mfem::ParGridFunction *u_gf_re,
      const mfem::ParGridFunction *v_gf_re,
      const mfem::ParGridFunction *u_gf_im = nullptr,
      const mfem::ParGridFunction *v_gf_im = nullptr)
      : _u_gf_re(u_gf_re), _u_gf_im(u_gf_im), _v_gf_re(v_gf_re),
        _v_gf_im(v_gf_im), _coef(coef) {}
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);
  virtual ~VectorGridFunctionDotProductCoefficient() {}
};

// Auxsolver to project the dot product of two vector gridfunctions onto a third
// (scalar) GridFunction
class VectorGridFunctionDotProductAux : public CoefficientAux {
private:
  mfem::Coefficient *_scaling_coef;
  mfem::ParGridFunction *_u_gf_re;
  mfem::ParGridFunction *_u_gf_im;
  mfem::ParGridFunction *_v_gf_re;
  mfem::ParGridFunction *_v_gf_im;

  const std::string _u_gf_real_name;
  const std::string _v_gf_real_name;
  const std::string _u_gf_imag_name;
  const std::string _v_gf_imag_name;
  const std::string _scaling_coef_name;
  bool _complex_average;

public:
  VectorGridFunctionDotProductAux(const std::string &dot_product_gf_name,
                                  const std::string &dot_product_coef_name,
                                  const std::string &scaling_coef_name,
                                  const std::string &u_gf_real_name,
                                  const std::string &v_gf_real_name,
                                  const std::string &u_gf_imag_name = "",
                                  const std::string &v_gf_imag_name = "",
                                  const bool complex_average = false);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;
};

} // namespace hephaestus
