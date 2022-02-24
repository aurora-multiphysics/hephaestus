#pragma once
#include <memory>
#include <iostream>
#include <fstream>

namespace hephaestus
{

class Executioner
{
    public:
    Executioner(){}
    Executioner(const std::string & executioner_type,
             const double time_step,
             const double end_time);

    std::string type;
    double dt;
    double t_final;

};

} // namespace hephaestus