#include "coefficients.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class VectorPowerLawNLFIntegrator : public mfem::NonlinearFormIntegrator {
protected:
  mfem::Coefficient *Q;

private:
  mfem::DenseMatrix curlshape, curlshape_dFt, Jac;
  mfem::Vector J, vec, pointflux;
  double E0 = 0.0001;
  int n = 20;
  double Jc = 100000000;

public:
  VectorPowerLawNLFIntegrator(mfem::Coefficient &q) : Q(&q) {}

  virtual void AssembleElementGrad(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Ttr,
                                   const mfem::Vector &elfun,
                                   mfem::DenseMatrix &elmat) {

    int nd = el.GetDof();
    int dim = el.GetDim();

    mfem::Vector J(dim);
    mfem::DenseMatrix curlshape(nd, dim),
        curlshape_dFt(nd, dim); // both trial and test function in Nedelec
                                // space, represented with curlshape
    elmat.SetSize(nd * dim);
    elmat = 0.0;
    double w;

    const mfem::IntegrationRule *ir = IntRule;
    if (!ir) {
      ir = &(
          mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 3)); // <---
    }
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);
      Ttr.SetIntPoint(&ip);
      w = ip.weight / Ttr.Weight();
      w *= Q->Eval(Ttr, ip);           // multiply the PWconstantcoefficient
      el.CalcCurlShape(ip, curlshape); // curl operation on the shape function
      MultABt(curlshape, Ttr.Jacobian(), curlshape_dFt);
      // the curl operation of H(curl) space:  H(div)
      // u(x) = (J/w) * uh(xh)
      curlshape.MultTranspose(elfun, J); // compute the current density J

      double J_norm = J.Norml2();
      double J_de = E0 / Jc * (n - 1) * pow((J_norm / Jc), n - 2);
      // derivative factor (n-1)*E0/Jc*(CurlH.Norm/Jc)^(n-2)
      // the transpose may be needed AtA rather than AAt
      AddMult_a_AAt(w * J_de, curlshape_dFt, elmat); // (Curl u, curl v)*J_de*w
    }
  };

  virtual void AssembleElementVector(const mfem::FiniteElement &el,
                                     mfem::ElementTransformation &Ttr,
                                     const mfem::Vector &elfun,
                                     mfem::Vector &elvect) {
    int nd = el.GetDof(), dim = el.GetDim();

    mfem::DenseMatrix curlshape(nd, dim);
    // both trial and test function in Nedelec
    // space, represented with curlshape
    double w;
    J.SetSize(dim);
    Jac.SetSize(dim);
    pointflux.SetSize(dim);
    vec.SetSize(dim);
    elvect.SetSize(nd);
    elvect = 0.0;
    // cout << "elfun size " <<  elfun.Size() << endl;
    // cout << "Densemtrix row col " << nd <<" Elfun size " << elfun.Size() <<
    // endl;
    const mfem::IntegrationRule *ir = IntRule;
    if (!ir) {
      ir = &(
          mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + 3)); // <---
    }
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);
      Ttr.SetIntPoint(&ip);
      w = ip.weight / Ttr.Weight();
      w *= Q->Eval(Ttr, ip);           // multiply the PWconstantcoefficient
      el.CalcCurlShape(ip, curlshape); // curl operation on the shape function

      curlshape.MultTranspose(elfun, J); // compute the current density J
      Jac = Ttr.Jacobian(); // mapping Jacobian to the reference element

      curlshape.MultTranspose(elfun, vec); //
      Jac.MultTranspose(vec, pointflux);
      // double J_norm=  pow(J[0],2) + pow(J[1],2) ;
      double J_norm = J.Norml2();

      double J_f = E0 / Jc * pow((J_norm / Jc), n - 1);
      //  factor E0/Jc*(CurlH.Norm/Jc)^(n-1)
      // cout << "current level " <<  J_f << endl;
      pointflux *= w * J_f;
      Jac.Mult(pointflux, vec);
      curlshape.AddMult(vec, elvect); // (Curl u, curl v)*J_f*w
    }
  };
};

TEST(NonlinearIntegratorTest, CheckData) {
  mfem::Mesh mesh((std::string(DATA_DIR) + "coil.gen").c_str(), 1, 1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));

  mesh.EnsureNodes();
  int dim = mesh.Dimension();
  mfem::FiniteElementCollection *fec_ND;
  mfem::FiniteElementCollection *fec_RT;
  fec_ND = new mfem::ND_FECollection(1, pmesh->Dimension());
  fec_RT = new mfem::RT_FECollection(1, pmesh->Dimension());
  mfem::ParFiniteElementSpace HCurlFESpace(pmesh.get(), fec_ND);
  mfem::ParFiniteElementSpace HDivFESpace(pmesh.get(), fec_RT);
  mfem::ConstantCoefficient coeff(1.0);

  // Build the NonlinearForm
  const int vdim = 1;

  //* in weak form
  //* (ρ∇×H, ∇×v) + (μdH/dt, v) + (dBᵉ/dt, v) - <(ρ∇×H)×n, v>  = 0
  mfem::ParNonlinearForm k_test(&HCurlFESpace);
  k_test.AddDomainIntegrator(new VectorPowerLawNLFIntegrator(coeff));
  // k_test.Assemble();
  k_test.Setup();

  // Compare ceed with mfem.

  mfem::ParGridFunction x(&HCurlFESpace);
  mfem::ParLinearForm b(&HCurlFESpace);
  b = 0.0;
  mfem::Vector X(HCurlFESpace.GetTrueVSize()), B(HCurlFESpace.GetTrueVSize());
  b.ParallelAssemble(B);
  x.ParallelProject(X);

  k_test.Mult(B, X);
  x.Distribute(X);

  EXPECT_LT(x.Norml2(), 1.e-12);
}
