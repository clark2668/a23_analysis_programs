The contents of this folder deal with showing how to convolve a waveform with the diode response from ARA

# C++ Version

You need to have a functioning version of AraRoot intalled somewhere, and ROOT. Both need to be accesible from your shell environment as `ARA_UTIL_INSTALL_DIR` and `ROOTSYS`.

Compile like `make run_diode`

`run_diode.cpp` takes in data, finds a cal pulser, convolves with the diode response, and saves the output.

`diode_include.h` is a copy and clean-up of all the necessary functions from AraSim to make the diode convlution work.

The result will look like ![](waveform_and_diode.png = 200x200)

# Python Version

The python version is copied from Ben Hokanson-Fasig's version in PyRex. It shows the convlution with tunnel diode, but for a scipy generated unit_impulse instead of actual data.