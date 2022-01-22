#-------------------------------------------------------------------------------
# compile the UMFPACK demos (for GNU make and original make)
#-------------------------------------------------------------------------------

# UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.
# All Rights Reserved.  See ../Doc/License for License.

CODE = Main

default: $(CODE)

include /mirror/home_nobck/sergiogus/SuiteSparse/UFconfig/UFconfig.mk

C = g++ $(CFLAGS) $(UMFPACK_CONFIG) -I/mirror/home_nobck/sergiogus/SuiteSparse/UMFPACK/Include -I/mirror/home_nobck/sergiogus/SuiteSparse/AMD/Include \
    -I/mirror/home_nobck/sergiogus/SuiteSparse/UFconfig \
    -I/mirror/home_nobck/sergiogus/libSLU

INC = /mirror/home_nobck/sergiogus/SuiteSparse/UMFPACK/Include/umfpack.h /mirror/home_nobck/sergiogus/SuiteSparse/AMD/Include/amd.h /mirror/home_nobck/sergiogus/SuiteSparse/UFconfig/UFconfig.h /mirror/home_nobck/sergiogus/libSLU/SLU.h

LIBS = -lblas $(LIB)

UMFPACK = /mirror/home_nobck/sergiogus/SuiteSparse/UMFPACK/Lib/libumfpack.a /mirror/home_nobck/sergiogus/SuiteSparse/AMD/Lib/libamd.a

#-------------------------------------------------------------------------------
# Create the demo programs, run them, and compare the output
#-------------------------------------------------------------------------------

/mirror/home_nobck/sergiogus/libSLU/SLU.o: /mirror/home_nobck/sergiogus/libSLU/SLU.cpp $(INC)
	$(C) -c -o /mirror/home_nobck/sergiogus/libSLU/SLU.o /mirror/home_nobck/sergiogus/libSLU/SLU.cpp

$(CODE): /mirror/home_nobck/sergiogus/libSLU/SLU.o $(CODE).cpp $(INC)
	$(C) -o $(CODE) $(CODE).cpp /home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosCPP/ReadTXT.cpp /home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosCPP/Mesh.cpp /home/sergiogus/Desktop/CTTC/NavierStokes/Codigos/ArchivosCPP/Solver.cpp /mirror/home_nobck/sergiogus/libSLU/SLU.o $(UMFPACK) $(LIBS)

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(CODE)
