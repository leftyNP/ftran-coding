PROGRAM ExplicitSolution
! Equation is 1/a*d(u)/dt + DEL(w)=f. We define w = -K*DEL(u) 
! In operator form, the equation is D*w + Omega*u = F, which in flux form is
! (G*OmI*D+I)*w = G*OmI*F, where w,F are vectors, OmI, I are diagonal matrices, and G,D are full matrices
    USE DefineParameters
    USE Types
	USE SubroutineModules
	
	! Declare variables
	IMPLICIT NONE							! this prevents variables defined elsewhere to be used here
	!-------------------------------------------------------------------------------------------------------------------!
	REAL(PRCSN), DIMENSION(xmax,ymax)	:: phi, phi_old, Qsource, RHS
    REAL(PRCSN)				            :: t
	TYPE(Tensor_C)						:: D_c
	TYPE(Tensor_F)						:: D_f
    !-------------------------------------------------------------------------------------------------------------------!
	! Program Timer Variables
	REAL								:: etime, total     
	REAL, DIMENSION(2)					:: elapsed
	!-------------------------------------------------------------------------------------------------------------------!
	! System commands	
	CALL System("RM C/*")					! Deletes old files (would overide files with same name anyways)
    !-------------------------------------------------------------------------------------------------------------------!
	! Print to screen useful values for current run
	PRINT "(A55)", 'Boundary Condition: BC_alpha*u + BC_beta*(w,n) = BC_psi' ! Boundary conditions
	PRINT "(4X,A8,7X,A7,10X,A6)", 'BC_alpha','BC_beta','BC_psi'              ! Print the boundary condition values
	WRITE(*,"(3(4X,F12.6,4X))") (BC_alpha(i), BC_beta(i), BC_psi(i),i=1,4)
	        ! print total x/y/t points, values of dx/dy/dt, and convergence criteria dt/dx**2
	PRINT "(TR4,A12,TR7,A10,TR9,A10,TR9,/,TR4,3(I7,TR11))", "# time steps","# x points","# y points",tmax,xmax,ymax
	PRINT "(TR4,3(A2,TR17),A5,/,TR4,4(ES15.9,TR4))", "dt","dx","dy","CFL #",dt,dx,dy,dt*2*(d1/dx**2+d2/dy**2)
    !-------------------------------------------------------------------------------------------------------------------!
    ! Initialize variables to zero
    phi=0D0; phi_old=0D0; Qsource=0D0; RHS=0D0
    D_c%xx=0D0; D_c%yy=0D0; D_c%xy=0D0
	D_f%xx_yy%Xface = 0.D0; D_f%xx_yy%Yface = 0.D0; D_f%xy_yx%Xface = 0.D0; D_f%xy_yx%Yface = 0.D0
    !-------------------------------------------------------------------------------------------------------------------!
    ! Specify constants used
	t     = tmin    ! Simulation time
    !-------------------------------------------------------------------------------------------------------------------!
    ! Define non-time dependent variables
	CALL Phi_Value(phi_old)
	CALL Savefile(0,'a',phi_old)			! Save initial analytic solution
    IF (MMS == 0) phi = phi_old				! Set phi to function value if not using MMS
	!CALL LoadFromFile('C/c00000.m',phi)	! Initialize phi from file
	CALL Savefile(0,'c',phi)				! Save initial conditions
	CALL DiffusionTensor(D_c)				! Load cell-centered diffusion coefficient
 	CALL Cell_To_Face(D_c,D_f)				! Calculate face-centered diffusion coefficients
    !-------------------------------------------------------------------------------------------------------------------!

    !-------------------------------------------------------------------------------------------------------------------!
    ! Begin time iteration
    DO tnum=1,tmax
        t = tmin + dt*DBLE(tnum)
        phi_old = phi
        !---------------------------------------------------------------------------------------------------------------!
        CALL SourceFn(Qsource,phi_old)
        CALL ExplicitDiffusionDerivative(phi,D_f,RHS)
        phi = phi_old + dt*( RHS + Qsource )*DiffVel
        !---------------------------------------------------------------------------------------------------------------!
        CALL BoundaryConditions(phi,D_f)
        !---------------------------------------------------------------------------------------------------------------!
        ! Save u to a file if n is a certain multiple of a number (see value in MOD)
!        IF ((MOD(tnum,1)==0) .OR. (tnum==tmax))    THEN
!        !IF ((MOD(tnum,FLOOR(DBLE(tmax)/2D1))==0) .OR. (tnum==tmax))    THEN
!	        PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",t,", which is ",(t-tmin)/(dt*tmax)*100,"% of the way done"
!	       	CALL When
!	       	CALL Phi_Value(phi_old)			! Recycle phi_old as analytic solution
!	       	CALL Savefile(tnum,'a',phi_old)	! Save analytic solution
!        	CALL Savefile(tnum,'c',phi)
!        	!CALL Savefile(FLOOR(DBLE(tnum)/(DBLE(tmax)/2D1)),'c',phi)
!        	!CALL Savefile(FLOOR(DBLE(tnum)/(DBLE(tmax)/2D1)),'a',phi_old)
!        	!CALL Savefile(FLOOR(DBLE(tnum)/50D0),'a',phi_old)
!        	!CALL Savefile(FLOOR(DBLE(tnum)/50D0),'c',phi)
!        ENDIF
		!-------------------------------------------------------------------------------------------------------------------!    	
		! Steady-State end condition
		IF (MOD(tnum,10000)==0) PRINT *,tnum,SUM(ABS(phi - phi_old))*dx*dy,SStol,dt*tnum+tmin
		IF ((tnum>1) .AND. (SUM(ABS(phi - phi_old))*dx*dy<SStol))	THEN
!			PRINT "(A16,F10.5,A11,F6.2,A17)", "The time is now ",t,", which is ",(t-tmin)/(dt*tmax)*100,"% of the way done"
		!	!!CALL Savefile(tnum,'c',phi_c)
			CALL Savefile(1,'c',phi)	! Save file as c00001.m to make things easy/consistent
			CALL When
			PRINT *,t,tnum,SUM(ABS(phi-phi_old))*dx*dy
			EXIT
		ENDIF
		!-------------------------------------------------------------------------------------------------------------------!
    ENDDO
	IF (tnum == tmax) CALL Savefile(2,'c',phi)	! if steady-state is not reached, save last value of phi
	!-------------------------------------------------------------------------------------------------------------------!
	! Print the runtime of the program
	            ! Define(User time):  	time actually spent in your program;
	            ! Define(System time):	time spent in the operating system on the program's behalf
	total = etime(elapsed)
	PRINT *, 'Total time=',total,'User runtime =', elapsed(1), ' System runtime =', elapsed(2)
!-----------------------------------------------------------------------------------------------------------------------!
END PROGRAM ExplicitSolution
