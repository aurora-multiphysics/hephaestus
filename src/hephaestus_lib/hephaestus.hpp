#include "hephaestus_joule.hpp"
#include "hephaestus_hertz.hpp"

void run_hephaestus(int argc, char *argv[], hephaestus::Inputs inputs)
{
   if (inputs._formulation == "Joule")
   {
      joule_solve(argc, argv, inputs);
   }
   else if (inputs._formulation == "Hertz")
   {
      hertz_solve(argc, argv, inputs);
   }
   else
   {
      std::cout<<"Formulation name " << inputs._formulation << " not recognised. \n";
   }
}
