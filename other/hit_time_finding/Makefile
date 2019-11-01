LDFLAGS=-L${ARA_UTIL_INSTALL_DIR}/lib -L${shell root-config --libdir}
CXXFLAGS=-I${ARA_UTIL_INSTALL_DIR}/include -I${shell root-config --incdir}
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs) -lAraEvent -lRootFftwWrapper
