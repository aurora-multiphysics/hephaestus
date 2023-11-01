#pragma once
#include "mfem.hpp"

// MFEM integrator to calculate the Rhs term, when your test function is a vector value
//  represented by multiple copies of a scalar field i.e. H1^3
// (f, (∇ . v))



namespace mfem{
class DomainLFH1DivIntegrator : public mfem::LinearFormIntegrator
{
protected:
    mfem::Coefficient *Q;

private:
   mfem::Vector shape, divShape;
   mfem::DenseMatrix dshape, gshape;

public:
    /// Constructs the domain integrator (Q, grad v)
    DomainLFH1DivIntegrator (mfem::Coefficient &q)
        : Q(&q)
    {}

    ~DomainLFH1DivIntegrator () {};

    // virtual bool SupportsDevice() override { return true; }
    
    /** Given a particular Finite Element and a transformation (Tr)
         computes the element right hand side element vector, elvect. */
    virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                        mfem::ElementTransformation &Tr,
                                        mfem::Vector &elvect) override;

    // virtual void AssembleDeltaElementVect(const mfem::FiniteElement &fe,
    //                                      mfem::ElementTransformation &Trans,
    //                                      mfem::Vector &elvect) override;

    // using mfem::LinearFormIntegrator::AssembleRHSElementVect;

};
} //namespace mfem