MODULE DefineParameters
!-----------------------------------------------------------------------------------------------------------------------!
! This module, stored in a separate file, initializes all variables in the simulation
!     Scalar constant in original equation:     InvDiffTime*d(u)/dt + DEL(w)=f    
!     Boundary condition constants:             BC_alpha*u+BC_beta*(W,normal)=BC_psi
!     grid spacing:                             dx, dy
!     time step spacing:                        dt
!     physical system dimensions:               xlength, ylength
!     total simulated time:                     tlength
!     total number of points used:              xmax, ymax
!     total number of time steps:               tmax
!     initial time (to avoid division by 0):    tmin
!-----------------------------------------------------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	! Real Precision 
	INTEGER, 					PARAMETER :: PRCSN			= 8					! Precision for real numbers - PRCSN is precision with no vowels
    !-------------------------------------------------------------------------------------------------------------------!
    ! Grid parameters
    REAL(PRCSN),				PARAMETER 	:: dx			= 1D-1/2D0**0		! spacing of mesh in x-direction
    REAL(PRCSN),				PARAMETER 	:: dy			= 1D-1/2D0**0		! spacing of mesh in y-direction
    REAL(PRCSN),				PARAMETER 	:: xlength		= 1D0				! total length of x-direction physical space being simulated
    REAL(PRCSN),				PARAMETER 	:: ylength		= 1D0				! total length of y-direction physical space being simulated
    REAL(PRCSN),				PARAMETER 	:: xmin			= 0D0				! minimum value of x-direction physical space
    REAL(PRCSN),				PARAMETER 	:: ymin			= 0D0				! minimum value of y-direction physical space
	INTEGER, 					PARAMETER 	:: xmax=CEILING(xlength/dx)+1		! total number of nodes in x-direction
	INTEGER, 					PARAMETER 	:: ymax=CEILING(ylength/dy)+1		! total number of nodes in y-direction
    !-------------------------------------------------------------------------------------------------------------------!
	! Other Parameters
	REAL(PRCSN),			 	PARAMETER 	:: d1 = 1.726D1, d2 = 6.209D0, d3 = 0.8723D0	! Diffusion coefficients 
	REAL(PRCSN),			 	PARAMETER	:: DiffVel	= 1D0					! Particle Velocity - Constant 'a' in: 1/a*du/dt = DEL[ Diff_Tensor * DEL[u] ] + Q
	REAL(PRCSN),				PARAMETER	:: PI = 3.14159265358979323846		! Transcendental Contant Pi
	REAL(PRCSN),				PARAMETER	:: SStol = 1D-10					! Tolerance for Steady-State
	INTEGER,					PARAMETER	:: MMS = 1							! Turns on Method of Manufactured Solution (non-zero source term, 0 initial condition)
	INTEGER									:: i,j,i0,j0,icell,jcell			! indexing parameters
	INTEGER,					PARAMETER	:: ProblemType = 7					! 1=exponential, 2=quadratic(x+y), 3=sinlog, 4=cubic, 5=quadratic(x*y), 6=fourth power, 7=time only
    !-------------------------------------------------------------------------------------------------------------------!
    ! Time Parameters
    !REAL(PRCSN),				PARAMETER	:: dt=.4D0/(2D0*(d1/dx**2+d2/dy**2)) ! time step size, in terms of diffusion coefficients and grid size  
    REAL(PRCSN),				PARAMETER	:: dt			= 1D-6
    REAL(PRCSN),				PARAMETER	:: tlength		= 1D2				! total length of physical time being simulated
    REAL(PRCSN),				PARAMETER	:: tmin			= 1D-1				! starting value for time (to avoid division by zero)
    INTEGER,			 		PARAMETER  	:: tmax=CEILING(tlength/dt)			! total number of time steps in simulation
	INTEGER									:: tnum								! Timestep number of system (time = DBLE(tnum)+tmin)
    !-------------------------------------------------------------------------------------------------------------------!
	! Mixed Cell Parameters
	INTEGER,					PARAMETER	:: MixedMethod = 0					! Mixed Cell index used: 1=a, 2=b, 3=g, 4=h, 5=r, 6=s, 0=off
	!! Linear Interface
	!REAL(PRCSN),				PARAMETER	:: m=1D0, b=0D0, alpha=0.2D0		! 2 corners hit 
	REAL(PRCSN),				PARAMETER	:: m=0.8D0, b=0D0, alpha=8D1		! a single corner
	!REAL(PRCSN),				PARAMETER	:: m=1.35D0, b=-0.2D0, alpha=8D1	! no corners hit
	!REAL(PRCSN),				PARAMETER	:: m=1D0/7D0, b=0.4D0, alpha=0.2D0	 
	!! Circular Interface
	!REAL(PRCSN),				PARAMETER	:: a=0.15D0, alpha=1.D0, beta=1.D-2
	!REAL(PRCSN)							:: r
    !-------------------------------------------------------------------------------------------------------------------!
    ! Boundary Conditions: alpha*u_concentration + beta*(w_flux,n_normal) = psi
    INTEGER, 					PARAMETER	:: BC_Fun		= 1	! Specifies if BC are a function of position (1=yes, 0=no)
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_alpha     = (/1D0, 1D0, 1D0, 1D0/)	! Dirichelet boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_beta      = (/0D0, 0D0, 0D0, 0D0/)	! Neumann boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_psi       = (/0D0, 0D0, 0D0, 0D0/)	! Value on the boundary
	!-------------------------------------------------------------------------------------------------------------------!

	CONTAINS
	!-------------------------------------------------------------------------------------------------------------------!

	SUBROUTINE Phi_Value(phi)
	!-------------------------------------------------------------------------------------------------------------------!
	! Specifies phi on a cell-centered mesh using an analytic solution
	! The only input is the current time step (n), which determines the time
		! The output is a file saved for each time this subroutine is called
			REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT) :: phi
			REAL(PRCSN) 									:: x,y,time
			!-----------------------------------------------------------------------------------------------------------!
			! Initialize variable to zero and calculate current time
			phi  = 0D0
			time = tmin + DBLE(tnum)*dt
			!-----------------------------------------------------------------------------------------------------------!
			DO i=1,xmax
				DO j=1,ymax
					x = xmin + DBLE(i-1)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
					y = ymin + DBLE(j-1)*dy					
					! Manufactured Solution:  phi = EXP(x*y)
						IF (ProblemType == 1) phi(i,j) = EXP(x*y)
					! Manufacutred Solution:  phi = d1*x**2 + d2*y**2
						IF (ProblemType == 2) phi(i,j) = d1*x**2 + d2*y**2
					! Manufacutred Solution:  phi = sin(pi*x)*log(y+1)
						IF (ProblemType == 3) phi(i,j) = SIN(PI*x)*LOG(y+1.D0)
					! Manufacutred Solution:  phi = (x**3 + y**3)/6
						IF (ProblemType == 4) phi(i,j) = (d2*x**3 + d1*y**3)/6D0
					! Manufacutred Solution:  phi = (x**2 * y**2)
						IF (ProblemType == 5) phi(i,j) = x**2 * y**2
					! Manufacutred Solution:  phi = (x**4 + y**4)/12
						IF (ProblemType == 6) phi(i,j) = (x**4 + y**4)/12D0
					! Manufacutred Solution:  phi = t (constant in space)
						IF (ProblemType == 7) phi(i,j) = time!**2/2D0
					! Manufactured Solution - Time Dependent
						!phi(i,j) = EXP(-x*y*time)
					! Mixed Cell Manufacutred solution: Linear Interface (4th order)
						!IF (y < m*x+b) THEN
						!	phi(i,j) = x**3*(x-(y-b)/m) + alpha
						!ELSE
						!	phi(i,j) = y**3*(y-(m*x+b)) + alpha
						!ENDIF
					! Mixed Cell Manufacutred solution: Linear Interface (Source fn is continuous)
						!IF (y < m*x+b) THEN
						!	phi(i,j) =    k2*x*(y-(m*x+b)) + alpha
						!ELSE
						!	phi(i,j) = -k1*m*y*(y-(m*x+b)) + alpha
						!ENDIF
					! Mixed Cell Manufacutred solution: Linear Interface (0 on all boundaries)
						!IF (y>m*x+b) THEN
						!	phi(i,j) = (y-m*x-b)   * (x-1D0) * x * y * (y-1D0) + alpha
						!ELSE
						!	phi(i,j) = (x-(y-b)/m) * (x-1D0) * x * y * (y-1D0) + alpha
						!ENDIF
					! Analytic solution - exponetial
						!!phi(i,j) = SQRT(tmin/t)*EXP(-x**2/(4.D0*K%xx(i,j)*t))*EXP(-y**2/(4.D0*K%yy(i,j)*t))
						!!phi(i,j) = SQRT(tmin/t)*EXP(-x**2/(4.D0*t))*EXP(-y**2/(4.D0*t))
						!phi(i,j) = (tmin/time)*EXP(-(x**2+y**2)/(4.D0*time))
						!!phi(i,j) = SQRT(tmin/t)*EXP(-x**2/(4.D0*t))
						!!phi(i,j) = SQRT(tmin/t)*EXP(-y**2/(4.D0*t))
					! Steady state solution for continuous coefficients
						!phi(i,j) = (x+2.D0*k1)/(1.D0+4.D0*k1)
					! Steady state solution for discontinuous coefficients
						!IF (x < jump) 	THEN
						!	phi(i,j) = (k2*x+2.D0*k1*k2)/(4.D0*k1*k2+5.D-1*(k1+k2))
						!ELSE
						!	phi(i,j) = (k1*x+2.D0*k1*k2+jump*(k2-k1))/(4.D0*k1*k2+5.D-1*(k1+k2))
						!	!phi(i,j) = (k1*x+2.D0*k1*k2+5.D-1*(k2-k1))/(4.D0*k1*k2+5.D-1*(k1+k2))
						!ENDIF
					! Steady state solution for discontinuous coefficients [2nd order]
						!IF (x < jump) 	THEN
						!	phi(i,j) = k2*x**2
						!ELSE
						!	phi(i,j) = k1*x**2+jump**2*(k2-k1)
						!ENDIF
					! Manufacutred Solution: phi = k1*x + k2*y
						!phi(i,j) = k1*x + k2*y
					! Mixed Cell Manufactured solution: Linear interface (sqrt)
						!IF (y < m*x+b)	THEN
						!	phi(i,j) = m*SQRT(x+1.D0) + alpha
						!ELSE
						!	phi(i,j) = (y-b+m)/SQRT(x+1.D0) + alpha
						!ENDIF
					! Mixed Cell Manufacutred solution: Circular interface
						!x = (DBLE(i)-0.5D0)*dx - 0.5D0
						!y = (DBLE(j)-0.5D0)*dy - 0.5D0
						!r = SQRT(x**2 + y**2)
						!IF (r < a) THEN
						!	phi(i,j) = COS( pi*r/(2.D0*a) ) + alpha
						!ELSE
						!	phi(i,j) = beta*(r/a - 1.D0)**2 + alpha
						!ENDIF
				ENDDO
			ENDDO
			!-----------------------------------------------------------------------------------------------------------!
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Phi_Value

	SUBROUTINE SourceFn(Qsource,phi_o)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine defines the source function and references the value of phi at the previous time-step
	! The source function is defined as Q in d/dt (phi) + DEL(Flux) = Q
	! This can depend on time and vary with space as well. Default is 0
	!-------------------------------------------------------------------------------------------------------------------!
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT( IN)	:: phi_o	! phi old - function value at previous time step
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT)	:: Qsource	! Source function
		REAL(PRCSN)										:: x,y,time

		Qsource = 0.D0
		x=0.D0; y=0.D0; time=0.D0
		
		IF (MMS == 1) THEN		! Control to turn on/off source loop
			time = tmin + DBLE(tnum)*dt	
			DO i=1,xmax
				DO j=1,ymax
					x = xmin + DBLE(i-1)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
					y = ymin + DBLE(j-1)*dy	
					! Manufactured Solutions - Exponential 
					IF (ProblemType == 1) Qsource(i,j) = -EXP(x*y)*( d1*y**2 + d2*x**2 + 2D0*d3*(1D0+x*y) )
					! Manufactured Solutions - Quadratic
					IF (ProblemType == 2) Qsource(i,j) = -2D0*(d1**2 + d2**2)
					! Manufactured Solutions - Sin-Log
					IF (ProblemType == 3) Qsource(i,j) = SIN(pi*x)*( d1*pi**2*LOG(y+1D0)+d2/(y+1D0)**2 )	&
						 				- COS(pi*x)*(2D0*pi*d3)/(y+1D0)
					! Manufactured Solutions - Cubic
					IF (ProblemType == 4) Qsource(i,j) = -d1*d2*(x+y)
					! Manufacutred Solution:  phi = (x**2 * y**2)
					IF (ProblemType == 5) Qsource(i,j) = -2D0*( (d1*y**2 + d2*x**2) + 2D0*(d3*x*y))
					! Manufacutred Solution:  phi = (x**4 + y**4)/12
					IF (ProblemType == 6) Qsource(i,j) = -(d1*x**2 + d2*y**2)
					! Manufacutred Solution:  phi = t (constant in space)
					IF (ProblemType == 7) Qsource(i,j) = 1D0!-time
					! Manufactured Solution - Time Dependent
						!Qsource(i,j) = -EXP(-x*y*time)*( d1*(y*time)**2 + d2*(x*time)**2 + x*y 	&
						!							   + 2D0*d3*time*(x*y*time-1D0) )
				ENDDO
			ENDDO
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SourceFn
!-----------------------------------------------------------------------------------------------------------------------!
END MODULE DefineParameters
