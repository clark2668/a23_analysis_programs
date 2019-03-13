The contents of this folder deal with how to extract the Rayleigh parameters for noise in ARA

# Compiling

You need to have a functioning version of AraRoot intalled somewhere, and ROOT. Both need to be accesible from your shell environment as `ARA_UTIL_INSTALL_DIR` and `ROOTSYS`.

Compile like `make fit_rayleigh.cpp`

# Running

`save_fft.cpp` saves the FFTs to be run over later. `fit_rayleigh.cpp` performs the fitting on the outputs of the FFT saver.

I recommend computing the FFT and computing the rayleigh using a cluster, and the `.pbs` files will do that. Submit like `qsub fit_rayleigh.pbs`.