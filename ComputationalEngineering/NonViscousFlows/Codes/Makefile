#-------------------------------------------------------------------------------
# compile the UMFPACK demos (for GNU make and original make)
#-------------------------------------------------------------------------------

# UMFPACK Version 4.4, Copyright (c) 2005 by Timothy A. Davis.
# All Rights Reserved.  See ../Doc/License for License.

CODE = Main
DIRECTORIO = /home/sergiogus/Desktop/ComputationalEngineering/NonViscousFlows/Codes/
default: $(CODE)


C = mpic++

#-------------------------------------------------------------------------------
# Create the demo programs, run them, and compare the output
#-------------------------------------------------------------------------------

$(CODE): $(CODE).cpp $(INC)
	$(C) -o $(CODE) $(CODE).cpp $(DIRECTORIO)/CppCodes/Memory.cpp $(DIRECTORIO)/CppCodes/ReadData.cpp $(DIRECTORIO)/CppCodes/ParPro.cpp $(DIRECTORIO)/CppCodes/Mesher.cpp $(DIRECTORIO)/CppCodes/Solver.cpp  
#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(CODE)
