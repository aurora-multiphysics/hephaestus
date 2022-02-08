#include "boundary_conditions.hpp"



BCMap::BCMap(const std::string & bc_name,
    Array<int> bdr_attr,
    int num_boundaries)
    :name(bc_name),
    markers(num_boundaries)
{
    mfem::common::AttrToMarker(num_boundaries, bdr_attr, markers);
}
