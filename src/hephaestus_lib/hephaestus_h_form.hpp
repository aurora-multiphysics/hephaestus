// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
// mpirun -np 4 valgrind ./magnetodynamic_newton -m bulk1_fine.msh

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
//#ifndef MFEM_USE_PETSC
//#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
//#endif

int h_form_solve(int argc, char *argv[]);