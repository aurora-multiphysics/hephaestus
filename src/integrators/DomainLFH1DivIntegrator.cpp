#include "DomainLFH1DivIntegrator.hpp"


namespace mfem {
void DomainLFH1DivIntegrator ::AssembleRHSElementVect(
   const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
{
    const int dim = el.GetDim();
    const int dof = el.GetDof();
    double w, coeff;
    const int sdim = Tr.GetSpaceDim();
    
    dshape.SetSize(dof, dim);
    gshape.SetSize(dof, dim);
    divShape.SetSize(dof * dim);
    
    dshape = 0.0;
    gshape = 0.0;
    divShape = 0.0;

    elvect.SetSize(dof * dim);
    elvect = 0.0;

    const mfem::IntegrationRule *ir = IntRule;
    if (ir == NULL)
    {
        int intorder = 2 * el.GetOrder();
        ir = &mfem::IntRules.Get(el.GetGeomType(), intorder);
    }

    for (int q = 0; q < ir->GetNPoints(); q++)
    {
        const mfem::IntegrationPoint &ip = ir->IntPoint(q);

        Tr.SetIntPoint(&ip);
        el.CalcDShape(ip, dshape);
        
        w = ip.weight;
        
        coeff = Q->Eval(Tr, ip);  
        
        Mult(dshape, Tr.AdjugateJacobian(), gshape);
       
        gshape.GradToDiv(divShape);
        add(elvect, w * coeff, divShape, elvect);
    }
}
} //namespace mfem

// void DomainLFH1DivIntegrator ::AssembleDeltaElementVect(const mfem::FiniteElement &fe,
//                                         mfem::ElementTransformation &Trans,
//                                         mfem::Vector &elvect)
// {
//    MFEM_ABORT("Not implemented!");
// }