
PROGRAM MG_MAIN
!---------------------------------------------------------------------------------------------------
	USE DefineParameters
	USE Types
	USE MG_Utilities
	USE MG_Subroutines

	IMPLICIT NONE
	!----------------------------------------------------------------------------------------------!
	REAL(PRCSN),	DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1)	:: phi_new, phi_old	! Energy/Intensity/Concentration values
	REAL(PRCSN),	DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1)	:: phi_dlt0			! Change in concentration (phi_dlt = delta phi)
	REAL(PRCSN),	DIMENSION(0:GridSize1(1)+1,0:GridSize1(2)+1)	:: phi_dlt1			! 1st restricted grid size for phi_dlt
	REAL(PRCSN),	DIMENSION(0:GridSize2(1)+1,0:GridSize2(2)+1)	:: phi_dlt2			! 2nd restricted grid size for phi_dlt
	REAL(PRCSN),	DIMENSION(0:GridSize3(1)+1,0:GridSize3(2)+1)	:: phi_dlt3			! 3rd restricted grid size for phi_dlt
	REAL(PRCSN),	DIMENSION(GridSize0(1)+2,GridSize0(2)+2)		:: RHS, Res0		! RESidual and RHS of A*x=b for grid level 0
	REAL(PRCSN),	DIMENSION(GridSize1(1)+2,GridSize1(2)+2)		:: Res1				! RESidual for grid level 1
	REAL(PRCSN),	DIMENSION(GridSize2(1)+2,GridSize2(2)+2)		:: Res2				! RESidual for grid level 2
	REAL(PRCSN),	DIMENSION(GridSize3(1)+2,GridSize3(2)+2)		:: Res3				! RESidual for grid level 3
	REAL(PRCSN)														:: ResMag			! MAGnitude of the RESidual
	REAL(PRCSN)														:: time				! time of simulation
	INTEGER															:: n,k				! counting integer
	!----------------------------------------------------------------------------------------------!
	! variables to check elapsed time
	REAL				:: etime, total
	REAL, DIMENSION(2)	:: elapsed
	!----------------------------------------------------------------------------------------------!
	! System Commands
	CALL system("RM C/*.m")		! Delete old files
	!----------------------------------------------------------------------------------------------!
	! Print to screen useful values for current run
	PRINT "(A43,/,A25)", 'MultiGrid Solver for the Diffusion Equation','Written by Nick Patterson'
	PRINT "(A55)", 'Boundary Condition: BC_alpha*u + BC_beta*(w,n) = BC_psi' ! Boundary conditions
	PRINT "(8X,A8,11X,A7,14X,A6)", 'BC_alpha','BC_beta','BC_psi'              ! Print the boundary condition values
	WRITE(*,"(3(4X,F12.6,4X))") (BC_alpha(i), BC_beta(i), BC_psi(i),i=1,4)
	! print total x/y/t points, values of dx/dy/dt, and convergence criteria dt/dx**2
	PRINT "(TR4,A12,TR7,A10,TR9,A10,TR9,/,TR4,3(I7,TR11))", "# time steps","# x points","# y points",tmax,GridSize0(1),GridSize0(2)
	PRINT "(TR4,3(A2,TR17),/,TR4,3(ES15.9,TR4))", "dt","dx","dy",dt,GridSpace(1),GridSpace(2)
	CALL When
    !----------------------------------------------------------------------------------------------!
    ! Initialize variables to zero
    phi_new=0D0; phi_old=0D0
    Res0=0D0; Res1=0D0; Res2=0D0
    phi_dlt0=0D0; phi_dlt1=0D0; phi_dlt2=0D0
    n=0
    !----------------------------------------------------------------------------------------------!
	! Initialize nonzero variables
	CALL AnalyticFn(phi_old)				! Set initial analytic condition
	CALL Savefile(0,'a',phi_old)			! Save initial analytic condition
	IF (MMS/=1)		phi_new = phi_old		! If using MMS, initial condition is 0

	CALL Savefile(0,'c',phi_new)			! Save initial condition
	
	CALL DiffusionOperator(GridSize0,A0)	! Set coefficient matricies for all grid levels
	CALL DiffusionOperator(GridSize1,A1)
	CALL DiffusionOperator(GridSize2,A2)
	CALL DiffusionOperator(GridSize3,A3)

	CALL SaveGeneral((/(GridSize0(1)+2)*(GridSize0(2)+2),(GridSize0(1)+2)*(GridSize0(2)+2)/),'C/A0.m',A0)
	!CALL SaveGeneral((/(GridSize1(1)+2)*(GridSize1(2)+2),(GridSize1(1)+2)*(GridSize1(2)+2)/),'C/A1.m',A1)
	!CALL SaveGeneral((/(GridSize2(1)+2)*(GridSize2(2)+2),(GridSize2(1)+2)*(GridSize2(2)+2)/),'C/A2.m',A2)
	
	CALL SourceFn(Qsource)					! Calculate source function
	! Note: If A and/or Q is time dependent, simply move between TimeLoop and MG_Cycle
	IF 	(TypeMethod == 'SemiImplicit')	CALL DiffusionOperator(GridSize0,SI,1)
	!----------------------------------------------------------------------------------------------!
	! Iterate through time	
	TimeLoop: DO tnum = 1,tmax
		time = tmin + DBLE(tnum)*dt 	! this is equivalent to (time = time + dt)
		phi_old = phi_new				! Update phi_old 
		CALL RightHandSide(phi_old,RHS)	! Update right-hand-side of equation
		ResMag = MGtol*1D1				! Set ResMag > MGtol to enter MG_Cycle while loop
		!----------------------------------------------------------------------------------------------!
		! Go through multigrid V cycle until residual is smaller than the tolerance
		MG_Cycle: DO WHILE (ResMag > MGtol)
			!----------------------------------------------------------------------------------------------!
			! Pre Smooth
			CALL Smooth(GridSize0,Iter_Pre,phi_new,RHS)							! Calculate phi on level 0			
			Res0 = RHS-RESHAPE(MATMUL(A0,PACK(phi_new,.TRUE.)),GridSize0+2)		! Calculate current residual
			IF (Loudness/='Silent') PRINT '(TR4,I4,TR4,ES12.6,2X,A15)',tnum,SQRT(SUM(Res0**2)), 'Pre-Smooth'
			!----------------------------------------------------------------------------------------------!
			! Smooth
			! 	It is here where you can control the shape of the multigrid cycle. More complicated and automated stencils can be formed with ALLOCATABLE matricies
			phi_dlt1=0D0										! Initialize phi delta
			CALL Restriction(GridSize0+2,Res0,Res1)				! Restrict residual to level 1			
			CALL Smooth(GridSize1,Iter_Smooth,phi_dlt1,Res1)	! Calculate phi delta on level 1

			CALL Restriction(GridSize1+2,Res1,Res2)				! Restrict residual to level 2
			CALL Restriction(GridSize1+2,phi_dlt1,phi_dlt2)		! Restrict phi delta to level 2			
			CALL Smooth(GridSize2,Iter_Smooth,phi_dlt2,Res2)	! Calculate phi delta on level 2

			IF (ALL(GrdFctr>=4)) THEN
				CALL Restriction(GridSize2+2,phi_dlt2,phi_dlt3)	! Restrict phi delta to level 3
				CALL Restriction(GridSize2+2,Res2,Res3)			! Restrict residual to level 3
				CALL Smooth(GridSize3,Iter_Smooth,phi_dlt3,Res3)! Calculate phi delta on level 3

				CALL Prolongation(GridSize3+2,phi_dlt3,phi_dlt2)! Interpolate phi delta to level 2
				CALL Smooth(GridSize2,Iter_Smooth,phi_dlt2,Res2)! Calculate phi delta on level 2
			ENDIF

			CALL Prolongation(GridSize2+2,phi_dlt2,phi_dlt1)	! Interpolate phi delta to level 1
			CALL Smooth(GridSize1,Iter_Smooth,phi_dlt1,Res1)	! Calculate phi delta on level 1
			CALL Prolongation(GridSize1+2,phi_dlt1,phi_dlt0)	! Interpolate phi delta to level 0
			phi_new = phi_new + phi_dlt0						! Update phi_new
			!----------------------------------------------------------------------------------------------!			
			! Post Smooth
			IF (Loudness/='Silent') THEN
				Res0 = RHS-RESHAPE(MATMUL(A0,PACK(phi_new,.TRUE.)),GridSize0+2)
				PRINT '(TR4,I4,TR4,ES12.6,2X,A15)',tnum,SQRT(SUM(Res0**2)), 'Smooth'
			ENDIF
			CALL Smooth(GridSize0,Iter_Post,phi_new,RHS)			! Calculate phi_new on grid level 0
			!----------------------------------------------------------------------------------------------!
			! Calculate magnitude of residual
			Res0 = RHS-RESHAPE(MATMUL(A0,PACK(phi_new,.TRUE.)),GridSize0+2)
			ResMag = SQRT(SUM(Res0**2))
			IF (Loudness/='Silent') PRINT '(TR4,I4,TR4,ES12.6,2X,A15)',tnum,SQRT(SUM(Res0**2)), 'Post-Smooth'
			!----------------------------------------------------------------------------------------------!
			! Screen for NaNs
			IF (SUM(phi_new)+1D0 == SUM(phi_new)) THEN
				PRINT *,'NaNs detected - execution terminating'
				EXIT
			ENDIF
			!----------------------------------------------------------------------------------------------!
		ENDDO MG_Cycle
		!----------------------------------------------------------------------------------------------!

		!----------------------------------------------------------------------------------------------!
    	! Save concentration variable
    	! 		If using MMS, one must run to steady state, and no other output files are needed.
    	! 		If not, one must run for a given allotment of time, and output files as often as desired    	
		SaveType: IF (MMS==0)	THEN		! Save phi_c to a file if n is a certain multiple of a number (see value in MOD)
			IF ((MOD(tnum,FLOOR(DBLE(tmax)/1D1))==0) .OR. (tnum==tmax))    THEN	! gives a certain amount of total files
			!IF ((MOD(tnum,100)==0) .OR. (tnum==tmax))    THEN					! gives a file every certain amount of time steps
				PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",time,", which is ",(time-tmin)/(dt*tmax)*100,"% of the way done"
				n=n+1
				CALL Savefile(n,'c',phi_new)
				!CALL AnalyticFn(phi_old)
				!CALL Savefile(n,'a',phi_old)
				CALL When
			ENDIF
		ELSEIF (MMS==1) THEN		! Steady-State end condition
			IF (MOD(tnum,500)==0) PRINT *,tnum,SUM(ABS(phi_new - phi_old))*GridSpace(1)*GridSpace(2)

			IF ((tnum>1) .AND. (SUM(ABS(phi_new - phi_old))*GridSpace(1)*GridSpace(2)<SStol))	THEN
				PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",time,", which is ",(time-tmin)/(dt*tmax)*100,"% of the way done"

				CALL Savefile(1,'c',phi_new)	! Save file as c00001.m to make things easy/consistent
				CALL When
				PRINT *,tnum,SUM(ABS(phi_new-phi_old))*GridSpace(1)*GridSpace(2),'Stead-State'
				EXIT
			ENDIF
		ENDIF SaveType
		!---------------------------------------------------------------------------------------------------------------!
	ENDDO TimeLoop

	!-------------------------------------------------------------------------------------------------------------------!
	IF (MMS==1 .AND. tnum >= tmax) CALL Savefile(2,'c',phi_new)	! Save very last time step if program has gone this far	(useful if steady-state is not reached)
	!----------------------------------------------------------------------------------------------!
	! Print the runtime of the program
		! Define(User time):  	time actually spent in your program;
		! Define(System time):	time spent in the operating system on the program's behalf
	total = etime(elapsed)
	PRINT *, 'Total time=',total,'User runtime =', elapsed(1), ' System runtime =', elapsed(2)
!--------------------------------------------------------------------------------------------------!
END PROGRAM MG_MAIN
