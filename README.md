# Hephaestus

This project is intended to provide an interface to a selection of MFEM based electromagnetics
solvers. Hephaestus is still under active development and is being updated frequently.
# How to use this repository

1. Fork this repository using the `fork` button on the front page
2. Remove the fork relationship: `Settings > General > Advanced > Remove fork relationship`
3. Rename your repository to align with your project: `Settings > General > Advanced > Rename repository`
4. Clone the repository with `git clone --recursive <your-project-url-here>.git`
5. Run `bash rename.sh <your-project-name>`
6. (Optionally open the project in VS Code)
7. Start writing your software!

# Building software and tests
Hephaestus can be built with the following commands from the top level project directory:

    mkdir build
    cd build
    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ..
    make
    make test

(Or equivalently, you can use `cmake --build .` instead of `make`, and `ctest` instead of `make test`.)

