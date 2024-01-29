#pragma once
#include "coefficient_aux.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{

// The VectorGridFunctionDotProductCoefficient evaluates the dot product of two
// (possibly complex) gridfunctions
class VectorGridFunctionDotProductCoefficient : public mfem::Coefficient
{
private:
  const std::shared_ptr<mfem::ParGridFunction> _u_gf_re{nullptr};
  const std::shared_ptr<mfem::ParGridFunction> _u_gf_im{nullptr};
  const std::shared_ptr<mfem::ParGridFunction> _v_gf_re{nullptr};
  const std::shared_ptr<mfem::ParGridFunction> _v_gf_im{nullptr};
  mfem::Coefficient & _coef;

public:
  VectorGridFunctionDotProductCoefficient(
      mfem::Coefficient & coef,
      const std::shared_ptr<mfem::ParGridFunction> u_gf_re,
      const std::shared_ptr<mfem::ParGridFunction> v_gf_re,
      const std::shared_ptr<mfem::ParGridFunction> u_gf_im = nullptr,
      const std::shared_ptr<mfem::ParGridFunction> v_gf_im = nullptr)
    : _u_gf_re{std::move(u_gf_re)},
      _u_gf_im{std::move(u_gf_im)},
      _v_gf_re{std::move(v_gf_re)},
      _v_gf_im{std::move(v_gf_im)},
      _coef(coef)
  {
  }

  ~VectorGridFunctionDotProductCoefficient() override = default;

  double Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip) override;
};

// Auxsolver to project the dot product of two vector gridfunctions onto a third
// (scalar) GridFunction
class VectorGridFunctionDotProductAux : public CoefficientAux
{
private:
  std::shared_ptr<mfem::Coefficient> _scaling_coef{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _u_gf_re{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _u_gf_im{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _v_gf_re{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _v_gf_im{nullptr};

  const std::string _u_gf_real_name;
  const std::string _v_gf_real_name;
  const std::string _u_gf_imag_name;
  const std::string _v_gf_imag_name;
  const std::string _scaling_coef_name;
  bool _complex_average;

public:
  VectorGridFunctionDotProductAux(const std::string & dot_product_gf_name,
                                  const std::string & dot_product_coef_name,
                                  std::string scaling_coef_name,
                                  std::string u_gf_real_name,
                                  std::string v_gf_real_name,
                                  std::string u_gf_imag_name = "",
                                  std::string v_gf_imag_name = "",
                                  const bool complex_average = false);

  ~VectorGridFunctionDotProductAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;
};

} // namespace hephaestus
