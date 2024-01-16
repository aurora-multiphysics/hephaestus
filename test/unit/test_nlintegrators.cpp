#include "coefficients.hpp"
#include <catch2/catch_test_macros.hpp>

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

/** NonlinearOperator operator of the form:
    k --> (M + dt*S)*k + H(x + dt*v + dt^2*k) + S*v,
    where M and S are given BilinearForms, H is a given NonlinearForm, v and x
    are given vectors, and dt is a scalar. */
class NonlinearOperator : public mfem::Operator {
private:
  mfem::ParBilinearForm *blf;
  mfem::ParNonlinearForm *nlf;
  mutable mfem::HypreParMatrix *Jacobian;
  double dt;
  const mfem::Vector *x0;
  mutable mfem::Vector x1;
  const mfem::Array<int> &ess_tdof_list;

public:
  NonlinearOperator(mfem::ParBilinearForm *blf_, mfem::ParNonlinearForm *nlf_,
                    const mfem::Array<int> &ess_tdof_list_)
      : Operator(blf_->ParFESpace()->TrueVSize()), blf(blf_), nlf(nlf_),
        Jacobian(NULL), dt(0.0), x0(NULL), x1(height),
        ess_tdof_list(ess_tdof_list_){};

  /// Set current dt, v, x values - needed to compute action and Jacobian.
  void SetParameters(double dt_, const mfem::Vector *x0_) {
    dt = dt_;
    x0 = x0_;
  };

  //* residual = (ρ(∇×H), ∇×v) + (μdH/dt, v) + (dBᵉ/dt, v) - <(ρ∇×H)×n, v>  = 0
  /// Compute y = H(x0 + dt* dx/dt) + M dx/dt
  virtual void Mult(const mfem::Vector &dx_dt, mfem::Vector &residual) const {
    add(*x0, dt, dx_dt, x1);
    nlf->Mult(x1, residual);
    blf->TrueAddMult(dx_dt, residual);
  };

  /// Compute J = dy/d(dx/dt)
  /// Compute J = M + dt grad_H(x0 + dt* dx/dt)
  virtual mfem::Operator &GetGradient(const mfem::Vector &dx_dt) const {
    add(*x0, dt, dx_dt, x1);
    delete Jacobian;
    mfem::SparseMatrix &localJ = blf->SpMat();
    localJ.Add(dt, nlf->GetLocalGradient(x1));
    Jacobian = blf->ParallelAssemble(&localJ);
    mfem::HypreParMatrix *Je = Jacobian->EliminateRowsCols(ess_tdof_list);
    delete Je;
    return *Jacobian;
  };

  virtual ~NonlinearOperator() { delete Jacobian; };
};

TEST_CASE("NonlinearIntegratorTest", "[CheckData]") {
  mfem::Mesh mesh((std::string(DATA_DIR) + "cylinder-hex-q2.gen").c_str(), 1,
                  1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));

  mesh.EnsureNodes();
  int dim = mesh.Dimension();
  mfem::FiniteElementCollection *fec_ND;
  mfem::FiniteElementCollection *fec_RT;
  fec_ND = new mfem::ND_FECollection(1, pmesh->Dimension());
  fec_RT = new mfem::RT_FECollection(1, pmesh->Dimension());
  mfem::ParFiniteElementSpace HCurlFESpace(pmesh.get(), fec_ND);
  mfem::ConstantCoefficient coeff(1.0);
  mfem::ConstantCoefficient mu(1.0);

  //* in weak form
  //* (ρ∇×H, ∇×v) + (μdH/dt, v) + (dBᵉ/dt, v) - <(ρ∇×H)×n, v>  = 0
  mfem::ParNonlinearForm nlf_test(&HCurlFESpace);
  mfem::ParBilinearForm blf_test(&HCurlFESpace);
  mfem::ParBilinearForm lf_test(&HCurlFESpace);
  nlf_test.AddDomainIntegrator(new VectorPowerLawNLFIntegrator(coeff));
  blf_test.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mu));

  //* Assemble and finalize
  nlf_test.Setup();
  blf_test.Assemble();
  blf_test.Finalize();

  mfem::ParGridFunction x(&HCurlFESpace);
  mfem::ParLinearForm b(&HCurlFESpace);
  b = 1.0;
  x = 0.0;
  mfem::Vector X(HCurlFESpace.GetTrueVSize()), B(HCurlFESpace.GetTrueVSize());
  b.ParallelAssemble(B);
  x.ParallelProject(X);

  mfem::Array<int> ess_tdof_list;
  NonlinearOperator nl_oper(&blf_test, &nlf_test, ess_tdof_list);
  nl_oper.SetParameters(0.1, &X);

  // Solver for the Jacobian solve in the Newton method
  mfem::Solver *jacobian_solver;
  // Set up the Jacobian solver
  mfem::HyprePCG j_pcg(HCurlFESpace.GetComm());
  mfem::HypreAMS ams(&HCurlFESpace);
  ams.SetPrintLevel(1);
  j_pcg.SetTol(1e-7);
  j_pcg.SetMaxIter(300);
  j_pcg.SetPrintLevel(1);
  j_pcg.SetPreconditioner(ams);
  jacobian_solver = &j_pcg;

  // Newton solver for the hyperelastic operator
  mfem::NewtonSolver newton_solver;
  newton_solver.iterative_mode = true;
  newton_solver.SetSolver(*jacobian_solver);
  newton_solver.SetOperator(nl_oper);
  newton_solver.SetPrintLevel(1); // print Newton iterations
  newton_solver.SetRelTol(0.0);
  newton_solver.SetAbsTol(1e-14);
  // newton_solver.SetAdaptiveLinRtol(2, 0.5, 0.9);
  newton_solver.SetMaxIter(10);

  int _my_rank;
  MPI_Comm_rank(HCurlFESpace.GetComm(), &_my_rank);
  std::cout << "Starting on rank:" << _my_rank << std::endl;

  // nlf_test.Mult(B, X);
  newton_solver.Mult(B, X);
  MFEM_VERIFY(newton_solver.GetConverged(), "Newton Solver did not converge.");
  std::cout << "Finished on rank:" << _my_rank << std::endl;
  x.Distribute(X);
}
