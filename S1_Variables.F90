MODULE DefineParameters
!-----------------------------------------------------------------------------------------------------------------------!
! This module initializes all variables in the simulation
!     Scalar constant in original equation:     InvDiffTime*d(phi)/dt + DEL(F)=S
!     Boundary condition constants:             BC_alpha*u+BC_beta*(F,normal)=BC_psi
!     grid spacing:                             dx, dy
!     time step spacing:                        dt
!     physical system dimensions:               xlength, ylength
!     total simulated time:                     tlength
!     total number of points used:              xmax, ymax
!     total number of time steps:               tmax
!     initial time:								tmin ! This avoids division by zero in the analytic solution
!-----------------------------------------------------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	! Real Precision 
	INTEGER, 		PARAMETER	:: PRCSN		= 8				! Precision for real numbers - PRCSN is precision with no vowels
    !-------------------------------------------------------------------------------------------------------------------!
    ! Problem Control
	! ProblemType less than 1:	expressions with a full (spatial) taylor series
	! ProblemType less than 7: 	expressions for various (spatial) powers	[1 x**4 + y**4]  [2 x**3 + y**3]   [3 x**2 + y**2] [4 x**1.5D0 + y**1.5D0] [5 x + y] [6 x**0.5D0 + y**0.5D0]
	! ProblemType less than 14: expressions for various (temporal) powers	[7 time**3] [8 time**2] [9 time**1.5D0] [10 time] [11 time**0.5D0] [12 1/time**0.5D0] [13 1/time]
	! ProblemType in the 20s:	expressions for 2 solutions meeting along an interface
	! ProblemType in the 30s:	expressions from various papers' test problems
    INTEGER, 		PARAMETER	:: ProblemType	= 41
    INTEGER, 		PARAMETER	:: GridSize		= 0
    INTEGER, 		PARAMETER	:: TimeStep		= 0
	INTEGER,		PARAMETER	:: MMS 			= 1				! Turns on Method of Manufactured Solution (non-zero source term, 0 initial condition)
	INTEGER,		PARAMETER	:: MixedMethod  = 4				! Mixed Cell index used: 1=a, 2=b, 3=g, 4=h, 5=r, 6=s, 0=off
	!-------------------------------------------------------------------------------------------------------------------!
    ! Grid Parameters
	REAL(PRCSN), 	PARAMETER 	:: dx = 1D-1/2D0**GridSize		! x-grid spacing, complete with factor of 2 size control
	REAL(PRCSN), 	PARAMETER 	:: dy = 1D-1/2D0**(GridSize)	! y-grid spacing, complete with factor of 2 size control
    REAL(PRCSN),  	PARAMETER	:: xlength      = 1D0-dx		! total length of x-direction physical space being simulated
    REAL(PRCSN),  	PARAMETER	:: ylength      = 1D0-dy		! total length of y-direction physical space being simulated
    REAL(PRCSN),  	PARAMETER	:: xmin			= 0D0			! miniumum value of x-direction physical space being simulated
	REAL(PRCSN),  	PARAMETER	:: ymin			= 0D0			! miniumum value of y-direction physical space being simulated
	INTEGER, 		PARAMETER	:: xmax=CEILING(xlength/dx)+1	! total number of nodes in x-direction (not the maximum value of physical x, but maximum value of the indicies)
	INTEGER, 		PARAMETER	:: ymax=CEILING(ylength/dy)+1	! total number of nodes in y-direction
	INTEGER, 		PARAMETER	:: VecLength = (xmax+1)*ymax+xmax*(ymax+1)	! total number of elements in a vector
    !-------------------------------------------------------------------------------------------------------------------!
    ! Time Parameters
	REAL(PRCSN),  	PARAMETER	:: dt = 1D-4/4D0**TimeStep		! time step size
	!REAL(PRCSN),  	PARAMETER	:: dt = 5D-1/(k1/dx**2 + k2/dy**2)/DiffVel
	!REAL(PRCSN),  	PARAMETER	:: dt = 5D-2*(dx**2+dy**2)
    REAL(PRCSN),  	PARAMETER	:: tlength      = 1D1			! total length of physical time being simulated
    REAL(PRCSN),  	PARAMETER	:: tmin			= 1D-2			! starting value for time (to avoid division by zero)
    INTEGER, 		PARAMETER	:: tmax=CEILING(tlength/dt)		! total number of time steps in simulation
	INTEGER						:: tnum							! Timestep number of system (time = DBLE(tnum)+tmin)
    !-------------------------------------------------------------------------------------------------------------------!
	! Other Parameters
	!REAL(PRCSN), 	PARAMETER 	:: k1 = 1.726D1, k2 = 6.209D0, k3 = 0.8723D0	! Diffusion coefficients 
	REAL(PRCSN), 	PARAMETER 	:: k1 = 5D1, k2 = 1D0, k3 = 0D0
	!REAL(PRCSN), 	PARAMETER 	:: k1 = 0.997D0, k2=0.574D0, k3=2.86D-2	! Diffusion coefficients 
	REAL(PRCSN), 	PARAMETER	:: DiffVel	= 1D0				! Particle Velocity - Constant 'a' in: 1/a*du/dt = DEL[ Diff_Tensor * DEL[u] ] + Q
	REAL(PRCSN),	PARAMETER	:: PI = 3.14159265358979323846	! Transcendental Contant Pi
	REAL(PRCSN),	PARAMETER	:: CGtol = 1D-10				! Tolerance for Conjugate Gradient
	REAL(PRCSN),	PARAMETER	:: SStol = 1D-10				! Tolerance for Steady-State
	INTEGER						:: i,j,i0, icell, jcell			! indexing parameters
    !-------------------------------------------------------------------------------------------------------------------!
	! Mixed Cell Parameters
	REAL(PRCSN), 	PARAMETER	:: jump = .5D0
	!! Linear Interface
	REAL(PRCSN),	PARAMETER	:: m=1D0, b=0D0, alpha=0.2D0		! 2 corners hit 
	!REAL(PRCSN),	PARAMETER	:: m=0.8D0, b=0D0, alpha=8D1		! a single corner
	!REAL(PRCSN),	PARAMETER	:: m=1.35D0, b=-0.2D0, alpha=8D1	! no corners hit
	!REAL(PRCSN),	PARAMETER	:: m=1D0/7D0, b=0.4D0, alpha=0.2D0	 
	!REAL(PRCSN),	PARAMETER	:: m=1D-20, b=jump, alpha=1D0
	!! Circular Interface
	!REAL(PRCSN),	PARAMETER	:: a=0.15D0, alpha=1D0, beta=1D-2
	!REAL(PRCSN)					:: r
    !-------------------------------------------------------------------------------------------------------------------!
    ! Boundary Conditions: alpha*u_concentration + beta*(w_flux,n_normal) = psi
    INTEGER, 					PARAMETER	:: BC_Fun		= 1	! Specifies if BC are a function of position (1=yes, 0=no)
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_alpha     = (/1D0, 1D0, 1D0, 1D0/)	! Dirichelet boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_beta      = (/0D0, 0D0, 0D0, 0D0/)	! Neumann boundary weight
    REAL(PRCSN), DIMENSION(4), 	PARAMETER	:: BC_psi       = (/0D0, 0D0, 0D0, 0D0/)	! Value on the boundary
    !REAL(PRCSN), DIMENSION(4), PARAMETER   :: BC_alpha     = (/0D0, 0D0, 0D0, 0D0/)	
    !REAL(PRCSN), DIMENSION(4), PARAMETER   :: BC_beta      = (/1D0, 1D0, 1D0, 1D0/)
    !-------------------------------------------------------------------------------------------------------------------!

END MODULE DefineParameters
