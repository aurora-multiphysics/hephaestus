#include "integrators.hpp"

void MixedWeakDivergenceIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &trial_fe,
                                const mfem::FiniteElement &test_fe,
                                mfem::ElementTransformation &Trans,
                                mfem::DenseMatrix &elmat)
{
    //  Get number of dofs/ basis functions
    dim = trial_fe.GetDim();
    int trial_dof = trial_fe.GetDof();
    int test_dof = test_fe.GetDof();
    double w, coeff, lambda, mu, alpha;
    Jadj.SetSize(dim);

    // vector to store values of basis functions
    trial_shape.SetSize(trial_dof);
    test_div_shape.SetSize(dim * test_dof);

    // Derivitive placeholder for getting div
    test_d_shape.SetSize(test_dof, dim);
    test_g_shape.SetSize(test_dof, dim);
    

    elmat.SetSize(dim*test_dof, trial_dof);
    elmat = 0.0;
    pelmat.SetSize(dim*test_dof, trial_dof);
    
    mfem::Geometry::Type geom = Trans.GetGeometryType();
    
    const mfem::IntegrationRule *ir = &GetIntRule(trial_fe, test_fe, Trans);
   
    for(int i = 0; i < ir->GetNPoints(); i++)
    {
        const mfem::IntegrationPoint &ip = ir->IntPoint(i);
        Trans.SetIntPoint(&ip);

        trial_fe.CalcShape(ip, trial_shape);
        test_fe.CalcDShape(ip, test_d_shape);

        w = ip.weight;

        if (a)
        {
            alpha = a -> Eval (Trans, ip);
        }

        alpha *= -1;
        
        Mult(test_d_shape, Trans.AdjugateJacobian(), test_g_shape);
        trial_shape *= w;
        trial_shape *= coeff;

        test_g_shape.GradToDiv(test_div_shape);
        AddMultVWt(test_div_shape, trial_shape, elmat);
    }
}

const mfem::IntegrationRule& MixedWeakDivergenceIntegrator::GetIntRule(const mfem::FiniteElement &trial_fe,
                                   const mfem::FiniteElement &test_fe,
                                   mfem::ElementTransformation &Trans)
{
    int order = trial_fe.GetOrder()
                + Trans.OrderGrad(&test_fe)
                + Trans.OrderJ();
    
    return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
} 

const mfem::IntegrationRule& MixedWeakDivergenceIntegrator::GetIntRule(const mfem::FiniteElement &trial_fe,
                                   const mfem::FiniteElement &test_fe,
                                   mfem::ElementTransformation &Trans)
{
    int order = trial_fe.GetOrder()
                + Trans.OrderGrad(&test_fe)
                + Trans.OrderJ();
    
    return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
}