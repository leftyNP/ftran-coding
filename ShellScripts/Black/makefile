FC = gfortran	# compiler name
PN = expl

CFLAGS = -O2 		# compiler flags
# object files
OBJS = 	S1_Variables.F90 ../../S2_Subroutines.F90 ../../S3_Main.F90

all:
	$(FC) -o $(PN).out $(CFLAGS) $(OBJS)

.PHONY: clean
clean:	
	rm -rf *o *mod