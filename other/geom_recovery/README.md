The contents of this folder deal with the recovery of the A2/A3 calibrated geometry from araROOT revision 2908 to the most updated versions by Thomas.

This code basically supports the documentation contained at [DocDB 1791](http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1791).

# Compiling

You need to have a functioning version of araROOT intalled somewhere, and ROOT. Both need to be accesible from your shell environment as `ARA_UTIL_INSTALL_DIR` and `ROOTSYS`.

Compile like `make geomCheck_myl.cpp`