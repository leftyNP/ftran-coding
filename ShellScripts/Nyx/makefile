FC=pgf90
OBJS = 	S1_Variables.F90 ../../S2_Subroutines.F90 ../../S3_Main.F90
CFLAGS= -fastsse -Munroll -Mconcur -O3 -tp x64,amd64e,barcelona-64,core2-64
all: expl.out
expl.out: $(OBJS)
	$(FC) $(OBJS) -o expl.out $(CFLAGS)

.PHONY: clean
clean:	
	rm -rf *o *mod