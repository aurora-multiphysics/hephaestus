# Build Hephaestus image
# Get base image
FROM alexanderianblair/hephaestus-deps:latest

# By default checkout master branch
ARG WORKDIR="opt"
ARG build_git_sha="master"

# Get Hephaestus
RUN cd /$WORKDIR && \
    git clone https://github.com/aurora-multiphysics/hephaestus.git && \
    cd hephaestus && \
    git checkout "$build_git_sha" && \
    git submodule update --init --recursive

# Build Hephaestus
RUN cd /$WORKDIR/hephaestus/ && \
    mkdir build && \
    cd build && \
    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/$WORKDIR/mfem/build -DMFEM_COMMON_INCLUDES=/$WORKDIR/mfem/miniapps/common  .. && \
    make -j1

# Test Hephaestus
RUN cd /$WORKDIR/hephaestus/build && \
    make test
