#-----------------------------------------------------------------------------------#
# Makefile para compilar todos los archivos de Navier Stokes en 3D incompressible con Gauss - Seidel en paralelo

CODE = Main

CLASS1 = Memory
CLASS2 = ReadData
CLASS3 = ParPro
CLASS4 = Mesher
CLASS5 = PostProcessing
CLASS6 = Solver

DIRECTORIO_CLASES = CppCodes/

default: $(CODE)

C = mpic++

$(CODE): $(CODE).cpp $(INC)
	$(C) -o $(CODE) $(CODE).cpp $(DIRECTORIO_CLASES)$(CLASS1).cpp $(DIRECTORIO_CLASES)$(CLASS2).cpp $(DIRECTORIO_CLASES)$(CLASS3).cpp $(DIRECTORIO_CLASES)$(CLASS4).cpp $(DIRECTORIO_CLASES)$(CLASS5).cpp $(DIRECTORIO_CLASES)$(CLASS6).cpp -O3
	
#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(CODE)
