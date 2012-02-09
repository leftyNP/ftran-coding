! SOM3_Main.F90
!
! Written by Nick Patterson, January 2011
! 
! 3 files must be compiled for this program to run
!	SOM1_Variables.F90		: Loads all parameters, variables, initial/boundary conditions
!	SOM2_Subroutines.F90	: Contains all subroutines, functions, and modules called in program
!	SOM3_Main.F90			: Main program that runs the code
!
! Compile as follows:	gfortran SOM1_Variables.F90 SOM2_Subroutines.F90 SOM3_Main.F90
!
! CONVENTION: FORTRAN COMMANDS ARE ALL CAPS; variables I define may be lower case or MiXeD
! CONVENTION: Separate separate sections with line breaks, 
! CONVENTION: and try to bind sections together with !------! to fill an entire line
! 
! This code will solve the 2D diffusion equation, d(phi)/dt + DIV(Flux) = Source, Flux = - Kdiff*GRAD(phi)
! The discretization method for this code is the Support Operator Method (SOM)
! The numerical method is iterative, and since the matrix equation Z*phi=b has a symmetric positive defintied (SPD)
! matrix, the conjugate gradient method will be employed, perhaps with a multigrid preconditioner eventually.
!
! Each cell will have a diffusion contribution called a zone. The zone will track the current of the diffusion, 
! which can be calculated by the product of the area and flux of each face, A*F, and will equal some value b when 
! multiplied by phi, (AF)phi=b, and this is a 4x4 matrix multiplied by a 4x1 vector, one for each direction:
! left, right, top, and bottom. These are the 4 faces of each cell/zone.
! The matrix is given a fifth row and column, with the vector phi having a fifth entry as well, the cell-centered value
! Each entry in the fifth row/column will be the sum of the A*F contributions from that row/column
! The corner element (5,5) will be the sum of the other elements in the 5th row/column. Also, a term will be added
! which comes from the form of the diffusion equation.
! Boundary conditions are then applied to each zonal matrix. The system is now Z_zone*Phi=b_zone.
! The contributions from the fifth row can be substituted and removed from the system by solving for the cell-centered
! value. This makes the system fully described in the 4x4 section of the Z_zone system.
! The global system then assembles each row for the appropriate flux-face equation (such as F_left + F_right = 0)
! This system is then solved in a Z*Phi=b system where Phi is a vector with the number of elements equal to the 
! number of faces on the system.

PROGRAM SOM3_MAIN
! Equation is a*d(u)/dt + DEL(w)=f. We define w = -K*DEL(u) 
! This program solves the first order system of equations via the Support Operator Method (SOM)
    USE DefineParameters
    USE Types
    USE SOM_Utilities
	USE SOM_Subroutine_Module

	IMPLICIT NONE                     ! this prevents variables defined elsewhere to be used here
	!-------------------------------------------------------------------------------------------------------------------!
	REAL(PRCSN), DIMENSION(VecLength)		:: phi_vec, RHSvector, phi_vec_old	! Vector form of values with one row for each face
	REAL(PRCSN), DIMENSION(xmax,ymax)		:: phi_c, phi_old					! Cell-Centered values of phi (old and new)
	REAL(PRCSN), DIMENSION(6,xmax,ymax)		:: ZoneC							! Elements needed for Face_To_Cell
	REAL(PRCSN), DIMENSION(tmax,3)			:: error
	REAL(PRCSN)								:: time								! Time of physical system
	REAL(PRCSN)								:: x,y								! Spatial location of physical system
	INTEGER									:: n								! Time step (integer value)
	TYPE(Face_Centered)						:: Map								! Map of how each face corresponds to the global vector
	TYPE(Cell_Zone)			 				:: Zelement, Mflux					! Cell-wise matricies (4x4 for each cell)
	TYPE(Sparse_Matrix) 					:: Ztotal							! Coefficient matrix for global system (sparse)
	!TYPE(Global_Zone) 						:: Ztotal							! Coefficient matrix for global system (dense) ! This type is no longer used
    !-------------------------------------------------------------------------------------------------------------------!
	! Program Timer Variables
	REAL				    				:: etime, total   
	REAL, DIMENSION(2)	    				:: elapsed
	!-------------------------------------------------------------------------------------------------------------------!
	! System commands	
	CALL System("RM C/*.m")                    ! Deletes old files
	!-------------------------------------------------------------------------------------------------------------------!
	! Print to screen useful values for current run
	PRINT "(A55)", 'Boundary Condition: BC_alpha*u + BC_beta*(w,n) = BC_psi' ! Boundary conditions
	PRINT "(8X,A8,11X,A7,14X,A6)", 'BC_alpha','BC_beta','BC_psi'              ! Print the boundary condition values
	WRITE(*,"(3(4X,F12.6,4X))") (BC_alpha(i), BC_beta(i), BC_psi(i),i=1,4)
	! print total x/y/t points, values of dx/dy/dt, and convergence criteria dt/dx**2
	PRINT "(TR4,A12,TR7,A10,TR9,A10,TR9,/,TR4,3(I7,TR11))", "# time steps","# x points","# y points",tmax,xmax,ymax
	PRINT "(TR4,3(A2,TR17),A8,/,TR4,4(ES15.9,TR4))", "dt","dx","dy","dt/dx**2",dt,dx,dy,dt/(dx**2+dy**2)*2.D0
	CALL When
    !-------------------------------------------------------------------------------------------------------------------!
	! Initialize variables
	tnum = 0; n = 0										! Integers
	phi_vec = 0D0; phi_vec_old = 0D0; RHSvector = 0D0	! Vectors
	Map%xface = 0D0; Map%yface = 0D0					! Face-Centered
	Zelement%cz = 0D0; Mflux%cz = 0D0					! Cell-Zone
	ZoneC = 0D0; Ztotal%sgz = 0D0						! Others
	!-------------------------------------------------------------------------------------------------------------------!
	! Initialize phi
	DO i=1,xmax; DO j=1,ymax; 
		x = xmin + (DBLE(i)-0.5D0)*dx
		y = ymin + (DBLE(j)-0.5D0)*dy
		phi_c(i,j) = Phi_Value(x,y)
	ENDDO; ENDDO
	CALL Savefile(0,'a',phi_c)				! Save initial analytic conditions
	IF (MMS==1)		phi_c = 0D0 			! If using MMS, initial condition is 0
	CALL Savefile(0,'c',phi_c)				! Save initial conditions
	!-------------------------------------------------------------------------------------------------------------------!
	! Define variables needed to map from a given face to the correct location in the global system
	FORALL (i=1:VecLength) phi_vec_old(i) = i		! Ordered vector - this is recycled later
	CALL Vector_To_Face(phi_vec_old,Map)			! Map of how each face corresponds to the global vector
		!PRINT *, 'Layout of MapX (9x5)'
		!WRITE(*,'(5(I2,2X))') (FLOOR(Map%Xface(i,:)),i=1,xmax+1)
		!WRITE(*,*)
		!PRINT *, 'Layout of MapY (8x6)'
		!WRITE(*,'(6(I2,2X))') (FLOOR(Map%Yface(i,:)),i=1,xmax)
	!-------------------------------------------------------------------------------------------------------------------!
	CALL InverseFluxMatrix(Mflux)		! Calculate 4x4 system (per cell) of diffusion coefficients
	phi_old = phi_c
	!-------------------------------------------------------------------------------------------------------------------!		
	TimeLoop: DO tnum=1,tmax
		time = tmin + DBLE(tnum)*dt 	! this is equivalent to (time = time + dt)

		phi_old = phi_c
		phi_vec_old = phi_vec
		! NOTE: diffusion tensor values will eventually change with each time step and would need to be in this loop

		CALL Zonal_System(Mflux,phi_old,Zelement,RHSvector,ZoneC)		! Form each cell's zonal system - apply BC - remove cell-centered dependency

		CALL Assemble_Global_Matrix(Zelement,Ztotal,Map)				! Combine all zonal matricies into a global system
			IF (tnum==1) CALL SAVEGENERAL(8*VecLength,3,'C/A.m',Ztotal%sgz)	! Save coefficient matrix

		CALL Conjugate_Gradient(Ztotal,phi_vec,RHSvector,phi_vec_old)	! Solve the system of equations - use old solution as initial guess

		CALL Face_To_Cell(phi_vec,phi_c,ZoneC)							! Convert face-centered, vector-solution into a cell-centered, matrix form
			 !CALL Vector_To_Face(phi_vec,phi_f); CALL Face_To_Cell(phi_f,phi_c)
		!---------------------------------------------------------------------------------------------------------------!
		! Screen for NaNs
		IF (SUM(phi_c)+1D0 == SUM(phi_c)) THEN
			PRINT *,'NaNs detected - execution terminating'
			EXIT
		ENDIF
		!---------------------------------------------------------------------------------------------------------------!
    	! Save concentration variable
    	! 		If using MMS, one must run to steady state, and no other output files are needed.
    	! 		If not, one must run for a given allotment of time, and output files as often as desired    	
		SaveType: IF (MMS==0)	THEN		! Save phi_c to a file if n is a certain multiple of a number (see value in MOD)
			IF ((MOD(tnum,FLOOR(DBLE(tmax)/2D1))==0) .OR. (tnum==tmax))    THEN
			!IF ((MOD(tnum,400)==0) .OR. (tnum==tmax))    THEN
				PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",time,", which is ",(time-tmin)/(dt*tmax)*100,"% of the way done"
				n=n+1
				CALL Error_Calculation(phi_c,n,error)
				!CALL Savefile(tnum,'c',phi_c)
				CALL Savefile(n,'c',phi_c)
				CALL When
			ENDIF
		ELSEIF (MMS==1) THEN		! Steady-State end condition
			IF (MOD(tnum,2000)==0) PRINT *,tnum,SUM(ABS(phi_c - phi_old))*dx*dy
			IF ((tnum>1) .AND. (SUM(ABS(phi_c - phi_old))*dx*dy<SStol))	THEN
				PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",time,", which is ",(time-tmin)/(dt*tmax)*100,"% of the way done"
				CALL Savefile(1,'c',phi_c)	! Save file as c00001.m to make things easy/consistent
				CALL When
				PRINT *,tnum,SUM(ABS(phi_c-phi_old))*dx*dy,'Stead-State'
				EXIT
			ENDIF
		!---------------------------------------------------------------------------------------------------------------!
		ENDIF SaveType
	ENDDO TimeLoop
	!-------------------------------------------------------------------------------------------------------------------!
	IF (tnum >= tmax) CALL Savefile(2,'c',phi_c)	! Save very last time step if program has gone this far
	!-------------------------------------------------------------------------------------------------------------------!
	! Print the runtime of the program
	            ! Define(User time):  	time actually spent in your program;
	            ! Define(System time):	time spent in the operating system on the program's behalf
	total = etime(elapsed)
	PRINT *, 'Total time=',total,'User runtime =', elapsed(1), ' System runtime =', elapsed(2)
!-----------------------------------------------------------------------------------------------------------------------!
END PROGRAM SOM3_MAIN
