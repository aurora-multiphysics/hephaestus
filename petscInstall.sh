export WORKDIR=`pwd`
export PETSC_DIR=/opt/petsc/petsc-3.18.0
function _build_petsc() {
    cd $WORKDIR
    if [ -d "$WORKDIR/petsc" ] ; then
       return
    fi
    mkdir petsc
    cd petsc
    curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.18.0.tar.gz
    tar -xf petsc-3.18.0.tar.gz -C .
    cd petsc-3.18.0
    mkdir build
    ./configure \
        --prefix=$WORKDIR/petsc/petsc-3.18.0/build \
        --with-debugging=0 \
        --with-ssl=0 \
        --with-pic=1 \
        --with-openmp=1 \
        --with-mpi=1 \
        --with-shared-libraries=1 \
    --with-cxx-dialect=C++11 \
    --with-fortran-bindings=0 \
    --with-sowing=0 \
    --with-hypre-include=/usr/local/include \
    --with-hypre-lib=/usr/local/lib/libHYPRE-2.23.0.so \
    --with-superlu_dist-include=/usr/local/include \
    --with-superlu_dist-lib=/usr/local/lib/libsuperlu_dist.so \
    --download-metis=1 \
    --download-parmetis=1 \
    --download-fblaslapack=1 \
    --download-ptscotch=1 \
    --download-scalapack=1 \
    --download-mumps=0 \
    --download-slepc=1 \
    --with-64-bit-indices=1 \
#    --with-cc=/usr/local/Cluster-Apps/openmpi/gcc/9.3/4.0.4/bin/mpicc \
#    --with-cxx=/usr/local/Cluster-Apps/openmpi/gcc/9.3/4.0.4/bin/mpicxx \
#    --with-fc=/usr/local/Cluster-Apps/openmpi/gcc/9.3/4.0.4/bin/mpif90 \
#    --with-mpi-dir=/usr/local/software/spack/current/opt/spack/linux-rhel7-x86_64/gcc-7.2.0/openmpi-3.1.3-b5ihoskgy7ny3nty67rdeckedakesoqa \
    PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt
    echo "Running Make"
    make PETSC_DIR=$WORKDIR/petsc/petsc-3.18.0 PETSC_ARCH=arch-linux-c-opt    all
    make PETSC_DIR=$WORKDIR/petsc/petsc-3.18.0 PETSC_ARCH=arch-linux-c-opt    install
    make PETSC_DIR=$WORKDIR/petsc/petsc-3.18.0/build PETSC_ARCH= check
    cd ..
    cd ..
    export PETSC_DIR=$WORKDIR/petsc
}

_build_petsc

