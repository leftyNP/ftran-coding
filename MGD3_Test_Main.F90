MODULE MG_Unit_Tests
!-----------------------------------------------------------------------------------------------------------------------!
! This module contains a subroutine to call and test each subroutine used when running the actual program
! Normal Mode Call:	gfortran -O3 MG1* MG2* MG3*
! Debug Mode Call: 	gfortran -O3 MGD1* MG2* MGD3*
!-----------------------------------------------------------------------------------------------------------------------!
	USE DefineParameters
	USE Types
	USE MG_Utilities
	USE MG_Subroutines
	
	IMPLICIT NONE
	
    !-------------------------------------------------------------------------------------------------------------------!
    ! Declare any variables I need here
	!-------------------------------------------------------------------------------------------------------------------!
	
	CONTAINS

	SUBROUTINE	Random_Array(sz,matrix)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine creates a random matrix, based on the computer's clock time
	! The heart of this program is from a webpage about random numbers: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
		INTEGER, DIMENSION(2), 			INTENT( IN) :: sz
		REAL(PRCSN), DIMENSION(sz(1),sz(2)),	INTENT(OUT) :: matrix
		INTEGER :: i, n, clock
		INTEGER, DIMENSION(:), ALLOCATABLE :: seed
		
		CALL RANDOM_SEED(size = n)
		ALLOCATE(seed(n))
		
		CALL SYSTEM_CLOCK(COUNT=clock)
		
		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
		CALL RANDOM_SEED(PUT = seed)
		DEALLOCATE(seed)
		CALL Random_Number(matrix)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Random_Array

	SUBROUTINE Debug_SGS(sz)
	!-------------------------------------------------------------------------------------------------------------------!
	! Test the Symmetric Gauss-Seidel solver
		INTEGER,	DIMENSION(2),				INTENT( IN) :: sz
		REAL(PRCSN),DIMENSION(sz(1)*sz(2),sz(1)*sz(2))	:: rand, SPD	! rand=Random matrix, SPD = symmetric positive definite matrix
		REAL(PRCSN),DIMENSION(sz(1)*sz(2)) 				:: b, x, ex		! b=random RHS of matrix equation, ex = exact solution, x = solution from solver
		CHARACTER (LEN=30)									:: frmt			! format of write statement
	
		WRITE(frmt,'(A4,I4,A13)') '(4X,',sz(1)*sz(2),'(G12.6,2X))'  ! Write a format string: '(99(E16.9E2,2x))'

		CALL Random_Array((/sz(1)*sz(2),sz(1)*sz(2)/),rand)	! Create a random array 'rand'
		SPD = MATMUL(rand,TRANSPOSE(rand))	! Create a SPD array from rand
	
		CALL Random_Array((/sz(1)*sz(2),1/),ex)	! Create random vector to be identical solution
		
		b = MATMUL(SPD,ex)	! Calculate product of SPD*ex to get RHS of matrix equation: A*x=b
		
		WRITE(*,'(4X,80(A1))') ('=',i=1,80)

		WRITE(*,'(4X,A32,/4X,A37,/4X,A24,/4X,A106,/,A4)') 		&
		'Symmetric Gauss-Seidel Unit Test', 				&
		'Solves A*x=b with and without a guess',			&
		'Vectors will be of size ', sz(1)*sz(2), 			&
		' and matricies will be square of that same length',&
		'The matrix A is symmetric and positive-definite (SPD) and is generated randomly based on time of execution', &
		'A = '
		DO i=1,sz(1)*sz(2); WRITE(*,frmt) SPD(i,:); ENDDO
		WRITE(*,'(4X,A45)') 'The exact solution and RHS are given as [x b]'
		WRITE(*,'(4X,G12.6,2X,G12.6)') (x(i),b(i),i=1,sz(1)*sz(2))
	!	PRINT *,'SPD'
	!	WRITE(*,frmt) (SPD(i,:),i=1,sz(1)*sz(2))
	
	
	!print the following
	!input a matrix of size 'x' and 'y'
	!random SPD matrix is
	!random solution is 
	!rhs of matrix equation then is
	!output of gauss seidel is
	!error of compute solution with correct is
	!the GStol is 
	
	
	!	PRINT *,'SPD'
	!	WRITE(*,'(5(E16.9E2,2X))') (SPD(i,:),i=1,5)	
	!
	!	PRINT *,'b'
	!	WRITE(*,'(E16.9E2)') (b(i),i=1,5)	
	
		
		CALL SGS(sz,SPD,b,x)
	!	WRITE(*,'(E16.9E2)') (x(i),i=1,5)	
	!
	!	CALL SGS(sz,SPD,b,x,ex)
	!	PRINT *,x
	!
	!	PRINT *,'ex'
	!	WRITE(*,'(E16.9E2)') (ex(i),i=1,5)	
	!
	!	PRINT *,'SPD'
	!	WRITE(*,'(5(E16.9E2,2X))') (SPD(i,:),i=1,5)	
	!
	!	PRINT *,'b'
	!	WRITE(*,'(E16.9E2)') (b(i),i=1,5)	


	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Debug_SGS

!-----------------------------------------------------------------------------------------------------------------------!
END MODULE MG_Unit_Tests



PROGRAM MGD3_MAIN
!--------------------------------------------------------------------------------------------------!
	USE DefineParameters
	USE Types
	USE MG_Utilities
	USE MG_Subroutines
	USE MG_Unit_Tests

	IMPLICIT NONE
	!----------------------------------------------------------------------------------------------!
	REAL(PRCSN),	DIMENSION(GridSize0(1),GridSize0(2))	:: Enew, Eold	! Energy/Intensity/Concentration values
	REAL(PRCSN),	DIMENSION(GridSize0(1),GridSize0(2))	:: Res,  dE		! RESidual and change in concentration (dE)
	REAL(PRCSN),	DIMENSION(GridSize1(1),GridSize1(2))	:: Res1, dE1	! 1st restricted grid size for Res and dE
	REAL(PRCSN),	DIMENSION(GridSize2(1),GridSize2(2))	:: Res2, dE2	! 2nd restricted grid size for Res and dE
	REAL(PRCSN)												:: ResMag		! MAGnitude of the RESidual
	REAL(PRCSN)												:: time			! time of simulation
	INTEGER													:: n			! counting integer
	!----------------------------------------------------------------------------------------------!
	! variables to check elapsed time
	REAL				:: etime, total
	REAL, DIMENSION(2)	:: elapsed
	!----------------------------------------------------------------------------------------------!
	! Print to screen useful values for current run
	WRITE(*,'(4X,80(A1))') ('=',i=1,80)
	PRINT "(4X,A43,/,4X,A62,/,4X,A25)", 'MultiGrid Solver for the Diffusion Equation', 	&
		'DEBUG MODE - Runs Units Tests on Subroutines (MG2_Subroutines)' , 'Written by Nick Patterson'
	PRINT "(4X,A12,TR7,A10,TR9,A10,TR9,/,TR4,3(I7,TR11))", "# time steps","# x points","# y points",tmax,GridSize0(1),GridSize0(2)
	CALL When
    !----------------------------------------------------------------------------------------------!
    ! Initialize variables to zero
    
    !----------------------------------------------------------------------------------------------!
	CALL Debug_SGS((/2,3/))
	!----------------------------------------------------------------------------------------------!
	! Print the runtime of the program
	WRITE(*,'(4X,80(A1))') ('=',i=1,80)
	PRINT '(4X,A43,/,6X,A12,/,8X,A50,/,8X,A73)', 'Runtime Information for Program Performance','Definitions:', &
		'User time:  	time actually spent in your program',			&
		"System time:	time spent in the operating system on the program's behalf"
	total = ETIME(elapsed)
	PRINT *, '   ', 'Total time =',total,'User runtime =', elapsed(1), ' System runtime =', elapsed(2)
	WRITE(*,'(4X,80(A1))') ('=',i=1,80)
!--------------------------------------------------------------------------------------------------!
END PROGRAM MGD3_MAIN
