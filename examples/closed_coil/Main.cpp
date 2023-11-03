#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

double current(const mfem::Vector &x, double t) { return 10.; }

int main(int argc, char *argv[]) {

  // Refinement and order
  int par_ref_lvl = 1;
  int order = 2;

  // Total electrical current going around the coil. Must be nonzero, can be
  // changed later.
  double Jtotal = 10;

  // Attribute that defines the internal face over which we apply the potential
  // difference
  int electrode_attr = 7;

  // Mesh file
  std::string mesh_filename = "team7.g";

  // Domain attributes of the coil to be solved
  std::string coil_attr = "3 4 5 6";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.AddOption(&par_ref_lvl, "-ref", "--parallel-refinement",
                 "Parallel refinement level.");
  args.AddOption(&order, "-o", "--order", "Base functions order");
  args.AddOption(&Jtotal, "-I", "--Itotal", "Total electrical current.");
  args.AddOption(&mesh_filename, "-f", "--mesh-filename", "Mesh file name");
  args.AddOption(
      &electrode_attr, "-e", "--electrode",
      "Boundary attribute of mesh face where potential difference is applied.");
  args.AddOption(
      &coil_attr, "-cd", "--coil-domains",
      "List of coil domain attributes separated by spaces, e.g. \'1 3 4\'");

  args.Parse();

  MPI_Init(&argc, &argv);

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + mesh_filename).c_str(), 1, 1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  // This vector of subdomains will form the coil that we pass to
  // ClosedCoilSolver
  mfem::Array<int> coil_domains;

  // Parsing the string of attributes
  std::stringstream ss(coil_attr);
  int att;
  while (ss >> att)
    coil_domains.Append(att);

  // FES and GridFunctions
  mfem::ND_FECollection HCurl_Col(order, pmesh.get()->Dimension());
  mfem::ParFiniteElementSpace FES_HCurl(pmesh.get(), &HCurl_Col);
  mfem::ParGridFunction J(&FES_HCurl);

  // Time-dependent current
  mfem::FunctionCoefficient I_val(current);

  // Input parameters and maps
  hephaestus::GridFunctions gfs;
  gfs.Register("J", &J, false);
  hephaestus::FESpaces fes;
  fes.Register("HCurl", &FES_HCurl, false);
  hephaestus::BCMap bcs;
  hephaestus::Coefficients coefs;
  coefs.scalars.Register("I", &I_val, false);

  hephaestus::InputParameters coilsolver_pars;
  coilsolver_pars.SetParam("HCurlFESpaceName", std::string("HCurl"));
  coilsolver_pars.SetParam("JGridFunctionName", std::string("J"));
  coilsolver_pars.SetParam("IFuncCoefName", std::string("I"));

  hephaestus::ClosedCoilSolver coil(coilsolver_pars, coil_domains,
                                    electrode_attr, order);
  coil.Init(gfs, fes, bcs, coefs);
  mfem::ParLinearForm ccs_rhs;
  coil.Apply(&ccs_rhs);

  mfem::VisItDataCollection *visit_DC =
      new mfem::VisItDataCollection("ClosedCoil_Results", pmesh.get());
  visit_DC->RegisterField("J", &J);
  visit_DC->Save();

  MPI_Finalize();
}