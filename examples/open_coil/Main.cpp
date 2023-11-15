#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

int main(int argc, char *argv[]) {

  // Refinement and order
  int par_ref_lvl = 1;
  int order = 1;

  // Total electrical current going around the coil. Must be nonzero, can be
  // changed later.
  double Itotal = 10;

  // Attribute that defines the internal faces over which we apply the potential
  // difference
  std::pair<int, int> elec_attrs;

  // Mesh file
  std::string mesh_filename = "coil.gen";

  // Domain attributes of the coil to be solved and boundary attributes of the electrodes
  std::string coil_attr = "1";
  std::string elec_bdr_attr = "1 2";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.AddOption(&par_ref_lvl, "-ref", "--parallel-refinement",
                 "Parallel refinement level.");
  args.AddOption(&order, "-o", "--order", "Base functions order");
  args.AddOption(&Itotal, "-I", "--Itotal", "Total electrical current.");
  args.AddOption(&mesh_filename, "-f", "--mesh-filename", "Mesh file name");
  args.AddOption(
      &coil_attr, "-cd", "--coil-domains",
      "List of coil domain attributes separated by spaces, e.g. \'1 3 4\'");
  args.AddOption(
      &elec_bdr_attr, "-e", "--electrode-attrs",
      "List of electrode attributes separated by spaces, e.g. \'1 2\'. Must be two values.");

  // ADD AN OPTION FOR DEFINING THE BOUNDARY ATTRIBUTES OF THE ELECTRODES ////
  args.Parse();

  MPI_Init(&argc, &argv);

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + mesh_filename).c_str(), 1, 1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  // This vector of subdomains will form the coil that we pass to
  // OpenCoilSolver
  mfem::Array<int> coil_domains;

  // This vector of attributes will form the electrodes that we pass to
  // OpenCoilSolver
  mfem::Array<int> elec_bdr_array;


  // Parsing the string of attributes
  std::stringstream ss(coil_attr);
  int att;
  while (ss >> att)
    coil_domains.Append(att);

  std::stringstream ss_bdr(elec_bdr_attr);
  while (ss_bdr >> att)
    elec_bdr_array.Append(att);

  if (elec_bdr_array.Size() != 2)
    mfem::mfem_error("Electrode boundary attribute list must contain two attributes.");

  elec_attrs.first = elec_bdr_array[0];
  elec_attrs.second = elec_bdr_array[1];

  // FES and GridFunctions
  mfem::ND_FECollection HCurl_Col(order, pmesh.get()->Dimension());
  mfem::ParFiniteElementSpace FES_HCurl(pmesh.get(), &HCurl_Col);
  mfem::ParGridFunction J(&FES_HCurl);

  // Time-dependent current
  mfem::ConstantCoefficient I_val(Itotal);

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
  coilsolver_pars.SetParam("SourceName", std::string("J"));
  coilsolver_pars.SetParam("PotentialName", std::string("V"));
  coilsolver_pars.SetParam("IFuncCoefName", std::string("I"));

  hephaestus::OpenCoilSolver coil(coilsolver_pars, coil_domains,
                                    elec_attrs);
  coil.Init(gfs, fes, bcs, coefs);
  mfem::ParLinearForm ocs_rhs(&FES_HCurl);
  coil.Apply(&ocs_rhs);

  mfem::VisItDataCollection *visit_DC =
      new mfem::VisItDataCollection("OpenCoil_Results", pmesh.get());
  visit_DC->RegisterField("J", &J);
  visit_DC->Save();

  MPI_Finalize();
}