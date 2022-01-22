#-------------------------------------------------------------------------------
# compile the UMFPACK demos (for GNU make and original make)
#-------------------------------------------------------------------------------

# UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.
# All Rights Reserved.  See ../Doc/License for License.

CODE = Main

CLASS1 = Memory
CLASS2 = ReadData
CLASS3 = ParPro
CLASS4 = Mesher
CLASS5 = Solver
CLASS6 = SLUSolver
CLASS7 = PostProcessing

DIRECTORIO_CLASES = /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/
DIRECTORIO_SOLVER = /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/

default: $(CODE)

include $(DIRECTORIO_SOLVER)/SuiteSparse/UFconfig/UFconfig.mk

C = mpic++ $(CFLAGS) $(UMFPACK_CONFIG) $(DIRECTORIO_SOLVER)/SuiteSparse/UMFPACK/Include  $(DIRECTORIO_SOLVER)/SuiteSparse/AMD/Include \
	$(DIRECTORIO_SOLVER)/SuiteSparse/UFconfig \
	$(DIRECTORIO_SOLVER)/libSLU


INC = $(DIRECTORIO_SOLVER)/SuiteSparse/UMFPACK/Include/umfpack.h  $(DIRECTORIO_SOLVER)/SuiteSparse/AMD/Include/amd.h  $(DIRECTORIO_SOLVER)/SuiteSparse/UFconfig/UFconfig.h  $(DIRECTORIO_SOLVER)/libSLU/SLU.h

LIBS = -lblas $(LIB)

UMFPACK = $(DIRECTORIO_SOLVER)/SuiteSparse/UMFPACK/Lib/libumfpack.a $(DIRECTORIO_SOLVER)/SuiteSparse/AMD/Lib/libamd.a
#-------------------------------------------------------------------------------
# Create the demo programs, run them, and compare the output
#-------------------------------------------------------------------------------

/home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/libSLU/SLU.o: /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/libSLU/SLU.cpp $(INC)
	$(C) -c -o /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/libSLU/SLU.o /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/libSLU/SLU.cpp

$(CODE): /home/sergiogus/Desktop/CFD/Incompressible3D/Codes/SuiteSparse/libSLU/SLU.o $(CODE).cpp $(INC)
	$(C) -o $(CODE) $(CODE).cpp 
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS1).cpp 
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS2).cpp  
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS3).cpp 
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS4).cpp 
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS5).cpp  
	$(DIRECTORIO_CLASES)/CppCodes/$(CLASS6).cpp  
	$(DIRECTORIO_SOLVER)/libSLU/SLU.o 
	$(UMFPACK) 
	$(LIBS)
#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(CODE)
