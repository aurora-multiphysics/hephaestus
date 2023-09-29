#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    // Refinement and order
    int par_ref_lvl = 1;
    int order = 2;

    // Total electrical current going around the coil. Must be nonzero, can be changed later. 
    double Jtotal = 10;

    // Attribute that defines the internal face over which we apply the potential difference
    int electrode_attr = 1;

    // Mesh file
    std::string mesh_filename = "torus.exo";


    mfem::OptionsParser args(argc, argv);
    args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                   "Directory storing input data for tests.");
    args.AddOption(&par_ref_lvl, "-ref", "--parallel-refinement",
                   "Parallel refinement level.");
    args.AddOption(&order, "-o", "--order",
                   "Base functions order");
    args.AddOption(&Jtotal, "-J", "--Jtotal",
                   "Total electrical current.");
    args.AddOption(&mesh_filename, "-f", "--mesh-filename",
                   "Mesh file name");
    args.AddOption(&electrode_attr, "-e", "--electrode",
                   "Boundary attribute of mesh face where potential difference is applied.");

    args.Parse();


    // Set Mesh
    mfem::Mesh mesh((std::string(DATA_DIR) + mesh_filename).c_str(), 1, 1);
    std::shared_ptr<mfem::ParMesh> pmesh = std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
    for (int l=0; l<par_ref_lvl; ++l) pmesh->UniformRefinement();
    
    // Defining the subdomains
    hephaestus::Subdomain coil1("coil1", 1);
    hephaestus::Subdomain coil2("coil2", 2);

    // This vector of subomains will form the coil that we pass to ClosedCoilSolver
    std::vector<hephaestus::Subdomain> coil_domains{coil1, coil2};
    
    // FES and GridFunctions
    mfem::ND_FECollection HCurl_Col(order, pmesh.get()->Dimension());
    mfem::ParFiniteElementSpace FES_HCurl(pmesh.get(), &HCurl_Col);
    mfem::ParGridFunction J(&FES_HCurl);

    // Input parameters and maps
    hephaestus::GridFunctions gfs;
    gfs.Register("J",&J,false);
    hephaestus::FESpaces fes;
    fes.Register("HCurl",&FES_HCurl,false);
    hephaestus::BCMap bcs;
    hephaestus::Coefficients coefs;

    hephaestus::InputParameters coilsolver_pars;
    coilsolver_pars.SetParam("HCurlFESpaceName", std::string("HCurl"));
    coilsolver_pars.SetParam("JGridFunctionName", std::string("J"));

    hephaestus::ClosedCoilSolver coil(coilsolver_pars, coil_domains, Jtotal, electrode_attr, order);
    coil.Init(gfs,fes,bcs,coefs);
    mfem::ParLinearForm dummy;
    coil.Apply(&dummy);

    mfem::VisItDataCollection* visit_DC = new mfem::VisItDataCollection("results", pmesh.get());
    visit_DC->RegisterField("J_new", &J);
    visit_DC->Save();

    MPI_Finalize();

}