#include "vector_gridfunction_dot_product_aux.hpp"

namespace hephaestus
{

double
VectorGridFunctionDotProductCoefficient::Eval(mfem::ElementTransformation & T,
                                              const mfem::IntegrationPoint & ip)
{
  double coefValue;
  coefValue = _coef.Eval(T, ip);

  mfem::Vector u_re;
  mfem::Vector u_im;
  mfem::Vector v_re;
  mfem::Vector v_im;

  _u_gf_re->GetVectorValue(T, ip, u_re);
  _v_gf_re->GetVectorValue(T, ip, v_re);
  if (_u_gf_im == nullptr || _v_gf_im == nullptr)
  {
    return coefValue * (u_re * v_re);
  }
  else
  {
    _u_gf_im->GetVectorValue(T, ip, u_im);
    _v_gf_im->GetVectorValue(T, ip, v_im);
    return 0.5 * coefValue * (u_re * v_re + u_im * v_im);
  }
}

VectorGridFunctionDotProductAux::VectorGridFunctionDotProductAux(
    const std::string & dot_product_gf_name,
    const std::string & dot_product_coef_name,
    const std::string & scaling_coef_name,
    const std::string & u_gf_real_name,
    const std::string & v_gf_real_name,
    const std::string & u_gf_imag_name,
    const std::string & v_gf_imag_name,
    const bool complex_average)
  : CoefficientAux(dot_product_gf_name, dot_product_coef_name),
    _u_gf_real_name(u_gf_real_name),
    _v_gf_real_name(v_gf_real_name),
    _u_gf_imag_name(u_gf_imag_name),
    _v_gf_imag_name(v_gf_imag_name),
    _scaling_coef_name(scaling_coef_name),
    _complex_average(complex_average),
    _scaling_coef(nullptr),
    _u_gf_re(nullptr),
    _u_gf_im(nullptr),
    _v_gf_re(nullptr),
    _v_gf_im(nullptr)
{
}

void
VectorGridFunctionDotProductAux::Init(const hephaestus::GridFunctions & gridfunctions,
                                      hephaestus::Coefficients & coefficients)
{

  _scaling_coef = coefficients.scalars.Get(_scaling_coef_name);
  if (_scaling_coef == NULL)
  {
    MFEM_ABORT("Conductivity coefficient not found for Joule heating");
  }

  if (_complex_average)
  {
    _u_gf_re = gridfunctions.Get(_u_gf_real_name);
    if (_u_gf_re == NULL)
    {
      MFEM_ABORT("GridFunction " << _u_gf_real_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
    _v_gf_re = gridfunctions.Get(_v_gf_real_name);
    if (_v_gf_re == NULL)
    {
      MFEM_ABORT("GridFunction " << _v_gf_real_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
    _u_gf_im = gridfunctions.Get(_u_gf_imag_name);
    if (_u_gf_im == NULL)
    {
      MFEM_ABORT("GridFunction " << _u_gf_imag_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
    _v_gf_im = gridfunctions.Get(_v_gf_imag_name);
    if (_v_gf_im == NULL)
    {
      MFEM_ABORT("GridFunction " << _v_gf_imag_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
  }
  else
  {
    _u_gf_re = gridfunctions.Get(_u_gf_real_name);
    if (_u_gf_re == NULL)
    {
      MFEM_ABORT("GridFunction " << _u_gf_real_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
    _v_gf_re = gridfunctions.Get(_v_gf_real_name);
    if (_v_gf_re == NULL)
    {
      MFEM_ABORT("GridFunction " << _v_gf_real_name
                                 << " not found when initializing VectorGridFunctionDotProductAux");
    }
  }

  coefficients.scalars.Register(_coef_name,
                                new VectorGridFunctionDotProductCoefficient(
                                    *_scaling_coef, _u_gf_re, _v_gf_re, _u_gf_im, _v_gf_im),
                                true);

  CoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus
