# SuiteSparse library
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a 

# UMF PACK directory
U1 = -I SuiteSparse/UMFPACK/Include
U2 = -I SuiteSparse/SuiteSparse_config
U3 = -I SuiteSparse/AMD/Include

# Includes
U = -I includes

CC = cc
CFLAGS = -W -Wall
INCLUDES = $(U1) $(U2) $(U3) $(U)
LIBS = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack
SRC = $(wildcard ./src/*.c)
OBJ = $(SRC:.c=.o)

main : $(OBJ)
	@$(CC) $(LIBS) -o $@ $^ 

$(OBJ) : %.o: %.c
	@$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

clean :
	@rm -rf ./src/*.o
	@rm -rf ./graphics/*.png

mrproper : clean
	@rm main 