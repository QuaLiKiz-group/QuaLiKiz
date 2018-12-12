# QuaLiKiz
QuaLiKiz is a quasilinear gyrokinetic code rapid enough to ease systematic interface with experiments. The quasilinear approximation is justified over a wide range of tokamak core parameters where nonlinear and quasilinear turbulent fluxes have been shown to agree. QuaLiKiz computes quasilinear gyrokinetic turbulent heat, particle and angular momentum fluxes. It is coupled in integrated modelling platforms such as CRONOS and JETTO. QuaLiKiz is available to all users. It allows for extensive stand-alone interpretative analysis and for first principle based integrated predictive modelling.

## Install
1. [Clone the repository from GitHub](https://help.github.com/articles/cloning-a-repository/).
    * If you want to install as submodule of QuaLiKiz (preferred)

            git clone git@github.com:QuaLiKiz-group/QuaLiKiz.git

      and then

          git submodule init
          git submodule update

2. Compile QuaLiKiz. To install QuaLiKiz you will need a fortran compiler and usually an openmpi library. The platform-specific files can be found in [src/make.inc/](src/make.inc/). It already contains files for most systems QuaLiKiz is run on. The [./install.sh](./install.sh) auto-detects which Makefile you need for your system, so:

``` bash
./install.sh
make
````

After making, you should have a binary called `QuaLiKiz` in your root directory.

## Usage
To run QuaLiKiz, one first has to create a folder `input` with the input binaries, as well as the folder `output`, `output/primitive` and `debug`. Then, one needs to generate input binaries. To make this easier, we have developed tools in [MATLAB](https://github.com/QuaLiKiz-group/QuaLiKiz-matlabtools) and [Python](https://github.com/QuaLiKiz-group/QuaLiKiz-pythontools) to set up, validate and plot QuaLiKiz runs. Please continue this guide on the respective GitHub page. There is also a [wiki](https://github.com/Karel-van-de-Plassche/QuaLiKiz/wiki) available with references and reading material.

## Disclaimer
QuaLiKiz is free and open-source software. If you have used QuaLiKiz in your own work, please cite our latest paper, [J. Citrin et al. PPCF 2017](http://iopscience.iop.org/article/10.1088/1361-6587/aa8aeb).
