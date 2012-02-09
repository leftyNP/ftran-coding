MODULE DefineParameters
!-----------------------------------------------------------------------------------------------------------------------!
! This module initializes all variables in the simulation
!-----------------------------------------------------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	! Real Precision 
	INTEGER, 					PARAMETER	:: PRCSN	  = 8				! Precision for real numbers - PRCSN is precision with no vowels
	!-------------------------------------------------------------------------------------------------------------------!
	! Problem Parameters
	INTEGER, DIMENSION(2), 		PARAMETER	:: GrdFctr    	= (/3,3/)		! Integer factor to determine size of all grids - GrdFctr is GridFactor with no vowels
	CHARACTER (LEN=20),			PARAMETER	:: TypeMethod	='Implicit'		! String to determine mode for problem to run in {Implict, Explicit, SemiImplicit}
	INTEGER,					PARAMETER	:: MMS 			= 1				! Turns on Method of Manufactured Solution (non-zero source term, 0 initial condition)
	INTEGER, 					PARAMETER	:: ProblemType	= -3			! Select pre-defined problem/source
	CHARACTER (LEN=20),			PARAMETER	:: Loudness		='Silent'		! String to determine mode to print out statements
	CHARACTER (LEN=20),			PARAMETER	:: SolverType	='SGS'			! String to determine solver type, either SGS for symmetric Gauss-Seidel, or RBGS for Red-Black Gauss-Seidel
	!-------------------------------------------------------------------------------------------------------------------!
    ! Grid Parameters 	! Grid sizes must be of the form 2^n+1: {2,3,5,9,17,33,65,129,257,...}
	INTEGER,	 DIMENSION(2),	PARAMETER	:: GridSize0  =-1+(/2**GrdFctr(1), 2**GrdFctr(2)/)
	INTEGER,	 DIMENSION(2),	PARAMETER	:: GridSize1  =-1+(/2**(GrdFctr(1)-1), 2**(GrdFctr(2)-1)/)
	INTEGER, 	 DIMENSION(2),	PARAMETER	:: GridSize2  =-1+(/2**(GrdFctr(1)-2), 2**(GrdFctr(2)-2)/)
	INTEGER, 	 DIMENSION(2),	PARAMETER	:: GridSize3  =-1+(/2**(GrdFctr(1)-3), 2**(GrdFctr(2)-3)/)
	REAL(PRCSN), DIMENSION(2), 	PARAMETER 	:: Sys_Length = (/1D0, 1D0/)	! total length of physical space being simulated
    REAL(PRCSN),  				PARAMETER	:: xmin		  = 0D0				! miniumum value of x-direction physical space being simulated
	REAL(PRCSN),  				PARAMETER	:: ymin		  = 0D0				! miniumum value of y-direction physical space being simulated
	REAL(PRCSN), DIMENSION(2), 	PARAMETER 	:: GridSpace  = Sys_Length/(GridSize0+1)
    !-------------------------------------------------------------------------------------------------------------------!
    ! Time Parameters
	REAL(PRCSN),  				PARAMETER	:: dt = 1.5D-4					! time step size
	!REAL(PRCSN),  				PARAMETER	:: dt = 5D-1/(k1/dx**2 + k2/dy**2)/DiffVel
    REAL(PRCSN),  				PARAMETER	:: tlength      = 1.8D0			! total length of physical time being simulated
    REAL(PRCSN),  				PARAMETER	:: tmin			= 1D-2			! starting value for time (to avoid division by zero)
    INTEGER, 					PARAMETER	:: tmax=CEILING(tlength/dt)		! total number of time steps in simulation
	INTEGER									:: tnum							! Timestep number of system (time = DBLE(tnum)+tmin)
    !-------------------------------------------------------------------------------------------------------------------!
	! Other Parameters
	!REAL(PRCSN), 				PARAMETER 	:: k1 = 1D0, k2 = 1D0, k3 = 0D0
	REAL(PRCSN), 				PARAMETER 	:: k1 = 0.997D0, k2=0.574D0, k3=2.86D-2	! Diffusion coefficients 
	REAL(PRCSN), 				PARAMETER	:: DiffVel	= 1D0				! Particle Velocity - Constant 'a' in: 1/a*du/dt = DEL[ Diff_Tensor * DEL[u] ] + Q
	REAL(PRCSN),				PARAMETER	:: PI = 3.14159265358979323846	! Transcendental Contant Pi
	REAL(PRCSN),				PARAMETER	:: MGtol = 1D-10				! Tolerance for MultiGrid residual
	REAL(PRCSN),				PARAMETER	:: SStol = 1D-8					! Tolerance for Steady-State
	INTEGER									:: i,j,i0						! indexing parameters
	INTEGER, 					PARAMETER	:: Iter_Pre		= 2				! Number of pre-smoothing iterations
	INTEGER, 					PARAMETER	:: Iter_Smooth	= 1				! Number of rough-smooths for each call in the MG cycle
	INTEGER, 					PARAMETER	:: Iter_Post	= 1				! Number of post-smoothing iterations
    !-------------------------------------------------------------------------------------------------------------------!
    ! Boundary Conditions: alpha*u_concentration + beta*(w_flux,n_normal) = psi
    ! 	Note the 4 values correspond to the sides in this order: left, right, top, bottom
    ! 		Dirichlet condition: BC_alpha=1, BC_beta=0
    !		Neumann condition:	 BC_alpha=0, BC_beta=1
    ! 		Periodic:			 BC_alpha=BC_beta=0
    !	Note: if BC is a function of position, this can only be set with a Dirichlet condition and having BC_fun=1
    INTEGER,					PARAMETER	:: BC_fun		= 1							! Specifies if boundary value on a side is a function of position (BC_fun=1) or a single value (BC_fun=0)
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_alpha     = (/1D0, 1D0, 1D0, 1D0/)	! Dirichelet boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_beta      = (/0D0, 0D0, 0D0, 0D0/)	! Neumann boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_psi       = (/0D0, 0D0, 0D0, 0D0/)	! Value on the boundary
    !-------------------------------------------------------------------------------------------------------------------!
END MODULE DefineParameters
