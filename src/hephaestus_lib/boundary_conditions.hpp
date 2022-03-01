#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include "mesh_extras.hpp"


namespace hephaestus
{

class BoundaryCondition
{
    public:
    BoundaryCondition();
    BoundaryCondition(const std::string & boundary_name,
          mfem::Array<int> boundary_ids);
    mfem::Array<int> getMarkers(mfem::Mesh & mesh);

    std::string name;
    std::function<double(const mfem::Vector&, double)> scalar_func;
    std::function<double(const mfem::Vector&, double)> scalar_func_im;
    std::function<void(const mfem::Vector&, mfem::Vector&)> vector_func;
    std::function<void(const mfem::Vector&, mfem::Vector&)> vector_func_im;

    mfem::Array<int> bdr_attributes;
    mfem::Array<int> markers;

    virtual void applyBC(mfem::LinearForm& b){};
    virtual void applyBC(mfem::ComplexLinearForm& b)
    {
    };
    virtual void applyBC(mfem::ParComplexLinearForm& b)
    {
        std::cout <<"usshs type3" << std::endl;
    };
    virtual void testecho()
    {
        std::cout <<"usshs";
    };
};

class NeumannBC : public BoundaryCondition
{
    public:
    NeumannBC();
    NeumannBC(const std::string & boundary_name,
          mfem::Array<int> boundary_ids);


    mfem::LinearFormIntegrator * lfi_re;
    mfem::LinearFormIntegrator * lfi_im;

    virtual void applyBC(mfem::LinearForm& b) override
    {
        std::cout <<"lagagaga type1";

        b.AddBoundaryIntegrator(lfi_re, markers);
    };
    virtual void applyBC(mfem::ComplexLinearForm& b) override
    {
        std::cout <<"lagagaga type2";

        b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
    };
    virtual void applyBC(mfem::ParComplexLinearForm& b) override
    {
        std::cout <<"lagagaga type3";
        b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
    };
    virtual void testecho() override
    {
        std::cout <<"lagagaga";
    };
};

// class ComplexRobinBC : public BoundaryCondition
// {
//     public:
//     ComplexRobinBC(const std::string & boundary_name,
//           mfem::Array<int> boundary_ids);

//     BilinearFormIntegrator integrator_re; //Neumann term
//     BilinearFormIntegrator integrator_im; //Neumann term

//     Coefficient robin_coef_re;
//     Coefficient robin_coef_im;

//     VectorCoefficient ;

//     applyNeumannBC(ParSesquilinearForm );
// }

class BCMap
{
    public:
    BCMap();

    void setBC(std::string bc_name, BoundaryCondition bc);
    BoundaryCondition getBC(std::string bc_name);

    std::map<std::string, hephaestus::BoundaryCondition*> bc_map;
};

}