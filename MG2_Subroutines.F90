
MODULE Types
!-------------------------------------------------------------------------------------------------------------------!
! This module introduces the derived data types
    USE DefineParameters    ! Note: all variables are set in this module, which must be compiled first

    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!    
    ! Derived Data Types
    TYPE Tensor_C
    ! Cell-centered tensor matrix
    ! Note the built in assumption that yx = xy
        REAL(PRCSN)	:: xx, xy, yy
    END TYPE Tensor_C

    TYPE Tensor_F
    ! Face-centered tensor matrix, which is defined in terms of cell-centered values
    	REAL(PRCSN)	:: xx,xy,yx,yy
    END TYPE Tensor_F
    	
!    TYPE Sparse_Matrix
!    ! This is a sparse matrix format to represent the global matrix system
!    ! In 2D, it is a 7 banded matrix, meaning a vector of 8*VecLength will be sufficient
!	! The abbreviation, sgz, is in homage to the Global_Zone call of %gz.
!    	REAL(PRCSN), DIMENSION(8*VecLength,3) :: sgz
!    END TYPE Sparse_Matrix    
!-------------------------------------------------------------------------------------------------------------------!
END MODULE Types

MODULE MG_Utilities
!-----------------------------------------------------------------------------------------------------------------------!
! This module contains the subroutines called in the MultiGrid (MG) program that are not physics related
! These utility subroutines include how to name a file, save a file, or perform repetative simple mathematical functions
!-----------------------------------------------------------------------------------------------------------------------!
    USE DefineParameters
    USE Types

	IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	CHARACTER(LEN=12)   :: filename
	!-------------------------------------------------------------------------------------------------------------------!
	
	CONTAINS

	SUBROUTINE Namefile(k,LTTR,filename)
	!-------------------------------------------------------------------------------------------------------------------!
	! Names each file to be save. Files are named by the iteration count (k)
		INTEGER, 			INTENT( IN) :: k		! Time step index
		CHARACTER(LEN=1),	INTENT( IN) :: LTTR 	! LTTR = LeTTeR = first letter of filename
		CHARACTER(LEN=12),	INTENT(OUT) :: filename	! Name generated from time-step value

		IF (k<10)							THEN	! File number is 1 digit
			WRITE(filename,'(A3,A1,A4,I1,A2)'), "C//",LTTR,"0000",k,".m"
		ELSE IF (ABS(k-54.5)<45.5)			THEN	! File number is 2 digits
			WRITE(filename,'(A3,A1,A3,I2,A2)'), "C//",LTTR,"000", k,".m"
		ELSE IF (ABS(k-549.5)<450.5)		THEN	! File number is 3 digits
			WRITE(filename,'(A3,A1,A2,I3,A2)'), "C//",LTTR,"00",  k,".m"
		ELSE IF (ABS(k-5499.5)<4500.5)		THEN	! File number is 4 digits
			WRITE(filename,'(A3,A1,A1,I4,A2)'), "C//",LTTR,"0",   k,".m"
		ELSE IF (ABS(k-54999.5)<45000.5)	THEN	! File number is 5 digits
			WRITE(filename,'(A3,A1,I5,A2)'), 	"C//",LTTR,    	  k,".m"
		ELSE
			PRINT *, "Filename cannot be set: the iteration number is greater than 5 digits"
			STOP
		ENDIF 
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Namefile

	SUBROUTINE Savefile(k,LTTR,u)
	!-------------------------------------------------------------------------------------------------------------------!
	! Saves the current instance of u. File names are appended based on the time iteration (k)
	! The exponentials are restricted to 2 digits (hence the variable u_restricted). 
	! This restriction is needed since my graphic software (Octave) fails to properly read exponents with >2 digits
	USE Types
		INTEGER, 													INTENT( IN) :: k			! Time step index
		CHARACTER(LEN=1),											INTENT( IN) :: LTTR 		! LTTR = LeTTeR = first letter of filename
        REAL(PRCSN), DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1),	INTENT (IN) :: u			! input matrix
		REAL(PRCSN), DIMENSION(GridSize0(1)+2,GridSize0(2)+2)    				:: u_restricted	! altered matrix to be saved to file
		CHARACTER (LEN=20)                          							:: frmt			! format of write statement

        WRITE(frmt,'(A1,I4,A13)') '(',GridSize0(2)+2,'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))'
	    ! Meaning: GridSize0(2) gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that
	    
	    ! Leave input matrix unchanged while modifying values for saving
        u_restricted(:,:) = u(0:GridSize0(1)+1,0:GridSize0(2)+1)
        !-------------------------------------------------------------------------------------------------------!
		! Ensure exponent is only 2 digits (further restricted to exponent above or below 50)
        WHERE (ABS(u_restricted) < 1.0D-50)                     ! If number is ridiculously small (+/-)
            u_restricted = 0.0D0
        ELSEWHERE (u_restricted > 9.9D+50)                      ! If number is ridiculously positive
            u_restricted = 9.9D+50
        ELSEWHERE (u_restricted < -9.9D+50)                     ! If number is ridiculously negative
            u_restricted = -9.9D+50
        END WHERE
		!-------------------------------------------------------------------------------------------------------!		
		! Name a file, open it, write to it, and close it	    
		CALL Namefile(k,LTTR,filename)							! Call subroutine to generate sequential file name
		OPEN (unit=1, file=filename)							! Open said file
		WRITE(1, '(A7,D15.7,TR5,A7,I5/)') '% dt = ', dt, 'tnum = ',tnum   ! Set first line to give time-step
        WRITE(1,frmt) (u_restricted(i,:),i=1,GridSize0(1)+2)
		CLOSE (1)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Savefile

	SUBROUTINE SaveGeneral(size,filename,matrix)
	!-------------------------------------------------------------------------------------------------------------------!
	! Subroutine writes array to file, where user inputs the dimensions of the matrix, the filename, and the matrix.
	! The advantage to this is that this same subroutine can handle an array of any size.
	! The file type is expected to be '.m', which is an Octave extension, although it is equivalent to a '.txt'.  
	! The file will have the same layout as the actual array, in terms of rows and columns.
	! Function Arguments:
	!	xsize		number of rows of the matrix
	!	ysize		number of columns of the matrix
	!	k			amount of characters needed for filename
	!	filename	name of file to be saved. A path must be included if you want file to go elsewhere than the current directory
	!	matrix		the array to be saved	
		INTEGER, DIMENSION(2),						INTENT (IN) :: size
		CHARACTER (LEN=*),                  		INTENT (IN) :: filename
        REAL(PRCSN), DIMENSION(size(1),size(2)),	INTENT (IN) :: matrix
		REAL(PRCSN), DIMENSION(size(1),size(2))    			    :: m_restricted
		CHARACTER (LEN=20)                          		    :: frmt

        WRITE(frmt,'(A1,I6,A13)') '(',size(2),'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))
	    ! Meaning: ymax gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that
	    
        m_restricted = matrix
        !-------------------------------------------------------------------------------------------------------!
		! Ensure exponent is only 2 digits (further restricted to exponent above or below 50)
        WHERE (ABS(m_restricted) < 1.0D-50)                     ! If number is ridiculously small (+/-)
            m_restricted = 0.0D0
        ELSEWHERE (m_restricted > 9.9D+50)                      ! If number is ridiculously positive
            m_restricted = 9.9D+50
        ELSEWHERE (m_restricted < -9.9D+50)                     ! If number is ridiculously negative
            m_restricted = -9.9D+50
        END WHERE
		!-------------------------------------------------------------------------------------------------------!		
		OPEN (unit=1, file=filename)                            ! Open said file

        WRITE(1,frmt) (m_restricted(i0,:),i0=1,size(1))
       
		CLOSE (1)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SaveGeneral

    SUBROUTINE When
	!-------------------------------------------------------------------------------------------------------------------!
	! Prints current time and date to screen when called
        INTEGER, DIMENSION(8)   :: values
        CALL DATE_AND_TIME(VALUES=values)

        IF (values(5) < 12) PRINT '(TR4,A6,I2.2,A1,I2.2,A1,I4,A9,I2,A1,I2.2,A1,I2.2,A3)', 'Date: ',        &
                values(2),'/',values(3),'/',values(1),'  Time: ',values(5),':',values(6),':',values(7),' AM'

		IF (values(5) == 12) PRINT '(TR4,A6,I2.2,A1,I2.2,A1,I4,A9,I2,A1,I2.2,A1,I2.2,A3)', 'Date: ',        &
                values(2),'/',values(3),'/',values(1),'  Time: ',values(5),':',values(6),':',values(7),' PM'

        IF (values(5) >  12) PRINT '(TR4,A6,I2.2,A1,I2.2,A1,I4,A9,I2,A1,I2.2,A1,I2.2,A3)', 'Date: ',        &
                values(2),'/',values(3),'/',values(1),'  Time: ',values(5)-12,':',values(6),':',values(7),' PM'
	!-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE When	
    
	SUBROUTINE Restriction(sz,A,B)
	!----------------------------------------------------------------------------------------------!
	! Coarsens input matrix (A) into matrix (B)
	! A is of size N x M, and B is of size (N+1)/2 x (M+1)/2
	! This is a full weighting stencil, meaning all 9 points (i,j)+/-{0,1} contribute
	!----------------------------------------------------------------------------------------------!
		INTEGER,	 DIMENSION(2), 							INTENT( IN)	:: sz	! Size of input matrix
		REAL(PRCSN), DIMENSION(sz(1),sz(2)), 				INTENT( IN)	:: A	! Input fine matrix
		REAL(PRCSN), DIMENSION((sz(1)+1)/2,(sz(2)+1)/2),	INTENT(OUT)	:: B	! Output coarse matrix
		INTEGER															:: c,d	! Indexing integers
		!------------------------------------------------------------------------------------------!
		! Initialize variables to 0
		B=0D0; c=0; d=0
		!--------------------------------------------------------------------------------------!
		DO j=1,sz(2),2; DO i=1,sz(1),2
			c=(i+1)/2
			d=(j+1)/2
			!--------------------------------------------------------------------------------------!
			! All Interior Points
			IF (i/=1 .OR. i/=sz(1) .OR. j/=1 .OR. j/=sz(2)) 								&
				B(c,d)=((A(i-1,j-1)+A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))						&
						+2D0*(A(i,j-1)+A(i,j+1)+A(i+1,j)+A(i-1,j))+4D0*A(i,j))/16D0
			!--------------------------------------------------------------------------------------!
			! Boundary terms (non-corners)
			! Top Row
			IF (i==1 		.AND. (j/=1 .OR. j/=sz(2)))  									&
				B(c,d)=((0D0+0D0+A(i+1,j-1)+A(i+1,j+1)) 									&
					+2D0*(A(i,j-1)+A(i,j+1)+A(i+1,j)+0D0)+4D0*A(i,j))/12D0		
			! Bottom Row
			IF (i==sz(1)	.AND. (j/=1 .OR. j/=sz(2)))										&
				B(c,d)=((A(i-1,j-1)+A(i-1,j+1)+0D0+0D0)										&
					+2D0*(A(i,j-1)+A(i,j+1)+0D0+A(i-1,j))+4D0*A(i,j))/12D0
			! Left Column
			IF (j==1		.AND. (i/=1 .OR. i/=sz(1)))										&
				B(c,d)=((0D0+A(i-1,j+1)+0D0+A(i+1,j+1))										&
					+2D0*(0D0+A(i,j+1)+A(i+1,j)+A(i-1,j))+4D0*A(i,j))/12D0
			! Right Column
			IF (j==sz(2) 	.AND. (i/=1 .OR. i/=sz(1)))										&
				B(c,d)=((A(i-1,j-1)+0D0+A(i+1,j-1)+0D0)										&
					+2D0*(A(i,j-1)+0D0+A(i+1,j)+A(i-1,j))+4D0*A(i,j))/12D0
			!--------------------------------------------------------------------------------------!
			! The four corner points
			! Top Left
			IF (i==1       .AND. j==1      )	B(c,d)=((0D0+0D0+0D0+A(i+1,j+1)) 			&
										+2D0*(0D0+A(i,j+1)+A(i+1,j)+0D0) + 4D0*A(i,j))/9D0 
			! Top Right
			IF (i==1	   .AND. j==sz(2))	B(c,d)=((0D0+0D0+A(i+1,j-1)+0D0) 			&
										+2D0*(A(i,j-1)+0D0+A(i+1,j)+0D0) + 4D0*A(i,j))/9D0 
			! Bottom Left
			IF (i==sz(1) .AND. j==1      )	B(c,d)=((0D0+A(i-1,j+1)+0D0+0D0) 			&
										+2D0*(0D0+A(i,j+1)+0D0+A(i-1,j)) + 4D0*A(i,j))/9D0 
			! Bottom Right
			IF (i==sz(1) .AND. j==sz(2))	B(c,d)=((A(i-1,j-1)+0D0+0D0+0D0) 			&
										+2D0*(A(i,j-1)+0D0+0D0+A(i-1,j)) + 4D0*A(i,j))/9D0 
			!--------------------------------------------------------------------------------------!
		ENDDO; ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE Restriction
	
	SUBROUTINE Prolongation(sz,A,B)
	!----------------------------------------------------------------------------------------------!
	! Refines input matrix (A) into matrix (B) (interpolation)
	! A is of size N x M, and B is of size 2N-1 x 2M-1
	! This is a full weighting stencil, meaning all 9 points (i,j)+/-{0,1} contribute
	!----------------------------------------------------------------------------------------------!
		INTEGER, 	 DIMENSION(2), 						INTENT( IN)	:: sz	! Size of input matrix
		REAL(PRCSN), DIMENSION(sz(1),sz(2)), 			INTENT( IN)	:: A	! Input coarse matrix
		REAL(PRCSN), DIMENSION(2*sz(1)-1,2*sz(2)-1),	INTENT(OUT) :: B	! Output fine matrix
		INTEGER														:: c,d	! Indexing integers
		!------------------------------------------------------------------------------------------!
		! Initialize variables to 0
		B=0D0; c=0; d=0
		!------------------------------------------------------------------------------------------!
		DO j=1,sz(2); DO i=1,sz(1)
			c=2*i-1
			d=2*j-1			
											B(c  ,d  )= A(i,j)										! center
			IF (i/=1)						B(c-1,d  )=(A(i,j)+A(i-1,j  ))/2D0						! down
			IF (i/=sz(1))					B(c+1,d  )=(A(i,j)+A(i+1,j  ))/2D0						! up
			IF (j/=1)						B(c  ,d-1)=(A(i,j)+A(i  ,j-1))/2D0				  	  	! left
			IF (j/=sz(2))					B(c  ,d+1)=(A(i,j)+A(i  ,j+1))/2D0						! right
			IF (i/=1 	 .AND. j/=1)		B(c-1,d-1)=(A(i,j)+A(i-1,j-1)+A(i,j-1)+A(i-1,j))/4D0	! down left
			IF (i/=1 	 .AND. j/=sz(2))	B(c-1,d+1)=(A(i,j)+A(i-1,j+1)+A(i,j+1)+A(i-1,j))/4D0 	! down right
			IF (i/=sz(1) .AND. j/=1)		B(c+1,d-1)=(A(i,j)+A(i+1,j-1)+A(i,j-1)+A(i+1,j))/4D0 	! up left
			IF (i/=sz(1) .AND. j/=sz(2))	B(c+1,d+1)=(A(i,j)+A(i+1,j+1)+A(i,j+1)+A(i+1,j))/4D0 	! up right
		ENDDO; ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE Prolongation

	SUBROUTINE RBGS(sz,Niter,A,b,x,guess)
	!----------------------------------------------------------------------------------------------!
	! RBGS = Red-Black Gauss-Seidel
	! This subroutine is a solver for A*x=b using the Gauss-Seidel method
	! The Red-Black component simply means that first even entries of x are solved, and then odd ones are solved
	! The odd, which is second, uses the recently found even points for a more accurate update
	!
	! This subroutine only sweeps through the system Niter number of times (N iterations)
	! To function as a stand-alone solver, the DO loop with Niter must be replaced with a while loop, 
	! and the old value of x must be stored. 
	! The while loop stops only when the magnitude of x-xold is below a given tolerance
	!----------------------------------------------------------------------------------------------!
		INTEGER, 	 DIMENSION(2), 							INTENT( IN)	:: sz		! Size of system
		INTEGER,											INTENT( IN) :: Niter	! Number of smoothing iterations
		REAL(PRCSN), DIMENSION(sz(1)*sz(2),sz(1)*sz(2)),	INTENT( IN)	:: A		! matrix of system
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), 				INTENT( IN)	:: b		! RHS of system
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), 				INTENT(OUT)	:: x		! Solution of system (to be found)
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), OPTIONAL,		INTENT( IN)	:: guess 	! initial guess of solution
		REAL(PRCSN)														:: total	! Dummy variable to track sums
		INTEGER 														:: i0,j0,k	! Iteration integers
		!------------------------------------------------------------------------------------------!
		! Initialize variables to zero
		total = 0D0; i0=0; j0=0; x=0D0
		!------------------------------------------------------------------------------------------!
		! Initialize guess
		IF (present(guess)) THEN
			x = guess
		ELSE
			FORALL (i0=1:sz(1)*sz(2)) x(i0)=b(i0)/A(i0,i0) ! This would be the result of the first sweep if x=0D0, so this saves time
		ENDIF
		!------------------------------------------------------------------------------------------!
		! Begin Gauss-Seidel Iterations
		DO k=1,Niter	! Perform Niter number of sweeps
			!--------------------------------------------------------------------------------------!
			! Red Loop
			DO i0=1,sz(1)*sz(2)
				IF (MOD(i0,2)==0)	THEN ! EVEN
					total = 0D0
					DO j0=1,sz(1)*sz(2); total=total+A(i0,j0)*x(j0); ENDDO
					total = total - A(i0,i0)*x(i0)
					x(i0) = (b(i0)-total)/A(i0,i0)
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------------!
			! Black Loop
			DO i0=1,sz(1)*sz(2)
				IF (MOD(i0,2)==1)	THEN ! ODD
					total = 0D0
					DO j0=1,sz(1)*sz(2); total=total+A(i0,j0)*x(j0); ENDDO
					total = total - A(i0,i0)*x(i0)
					x(i0) = (b(i0)-total)/A(i0,i0)
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------------!
		ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE RBGS

	SUBROUTINE SGS(sz,Niter,A,b,x,guess)
	!----------------------------------------------------------------------------------------------!
	! Symmetric Gauss-Seidel (SGS) algorithm to iteratively solve the system Ax=b
	! SGS (distinct from GS) because the direction of iteration alternates forwards/backwards
	! Convergence is guaranteed if matrix A is SPD or diagonally dominant.
	!
	! Matrix A and vector b are inputs for Ax=b, and an optional vector guess can be used as an initial guess
	! x is current best guess for solution   	x	 --> x(k+1)
	! xold is previous best guess for solution  xold --> x(k)
	!
	! For multigrid purposes, the SGS solver is not run until residual is less than some tolerance,
	! but it is run niter times, each time being considered one sweep.
	! Typical values are niter=2 for presmooth, 1 for smooth, and 1 for post-smooth
	!----------------------------------------------------------------------------------------------!
		INTEGER, 	 DIMENSION(2), 							INTENT( IN)	:: sz		! Size of system
		INTEGER,											INTENT( IN) :: Niter	! Number of smoothing iterations
		REAL(PRCSN), DIMENSION(sz(1)*sz(2),sz(1)*sz(2)),	INTENT( IN)	:: A		! matrix of system
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), 				INTENT( IN)	:: b		! RHS of system
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), 				INTENT(OUT)	:: x		! Solution of system (to be found)
		REAL(PRCSN), DIMENSION(sz(1)*sz(2)), OPTIONAL,		INTENT( IN)	:: guess 	! initial guess of solution
		REAL(PRCSN), DIMENSION(sz(1)*sz(2))						 		:: xold !,Res	! Res=residual, xold is previous guess for x
		!REAL(PRCSN)														:: ResMag	! ResMag is magnitude of residual
		REAL(PRCSN)														:: total	! Variable used to track summations
		INTEGER															:: k !, counter
		!------------------------------------------------------------------------------------------!
		! Initialize variables to 0
		x=0D0; xold=0D0; total=0D0!; Res=0D0; ResMag=0D0; counter=0
		!------------------------------------------------------------------------------------------!
		! Set initial value for xold	
		IF (present(guess)) THEN
			xold = guess
		ELSE
			FORALL (i=1:sz(1)*sz(2)) xold(i)=b(i)/A(i,i)
		ENDIF
		!------------------------------------------------------------------------------------------!
		x=xold						! sets initial vector to xold
		!ResMag=1D1*GStol			! sets magnitude of residual larger than tolerance to enter while-loop
		!------------------------------------------------------------------------------------------!
		DO k=1,Niter
		!DO WHILE (ResMag>GStol)		! Iterate solver until residual magnitude (ResMag) is less than the Gauss-Seidel tolerance (GStol)
			xold=x					! Saves non-updated vector x to calculate residual for convergence test
			!------------------------------------------------------------------------------------------!
			! forward loop
			DO i=1,sz(1)*sz(2)
				total=0D0
				DO j=1,sz(1)*sz(2); total = total + A(i,j)*x(j); ENDDO
				total = total - A(i,i)*x(i)		! Remove j=i contribution
				x(i) = (b(i)-total)/A(i,i)
			ENDDO
			!------------------------------------------------------------------------------------------!
			! backward loop
			DO i=sz(1)*sz(2),1,-1
				total=0D0
				DO j=sz(1)*sz(2),1,-1; total = total + A(i,j)*x(j); ENDDO
				total = total - A(i,i)*x(i)		! Remove j=i contribution
				x(i) = (b(i)-total)/A(i,i)
			ENDDO
			!------------------------------------------------------------------------------------------!
			! Calculate Residual
!			Res = x - xold
!			ResMag = SQRT(SUM(Res**2))
			!------------------------------------------------------------------------------------------!
! 			! Counter status
!			counter = counter+1
!			IF (MOD(counter,5000)==0) WRITE(*,'(4X,A18,2X,A14,I3,A1,4X,A9,2X,E16.9)') &
!				'SGS Solver Status:','Iteration # = ',counter/1000,'K','ResMag = ',ResMag
		ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE SGS
!-----------------------------------------------------------------------------------------------------------------------!
END MODULE MG_Utilities

MODULE MG_Subroutines
!-----------------------------------------------------------------------------------------------------------------------!
! This module contains the subroutines called in the MG program that are used to solve the problem
!-----------------------------------------------------------------------------------------------------------------------!
	USE DefineParameters
	USE Types
	USE MG_Utilities

	IMPLICIT NONE
	!-------------------------------------------------------------------------------------------------------------------!
	! Shared Variables	
	REAL(PRCSN), DIMENSION((GridSize0(1)+2)*(GridSize0(2)+2),(GridSize0(1)+2)*(GridSize0(2)+2))	:: A0		! Coefficient matrix for grid 0
	REAL(PRCSN), DIMENSION((GridSize1(1)+2)*(GridSize1(2)+2),(GridSize1(1)+2)*(GridSize1(2)+2))	:: A1		! Coefficient matrix for grid 1
	REAL(PRCSN), DIMENSION((GridSize2(1)+2)*(GridSize2(2)+2),(GridSize2(1)+2)*(GridSize2(2)+2))	:: A2		! Coefficient matrix for grid 2
	REAL(PRCSN), DIMENSION((GridSize3(1)+2)*(GridSize3(2)+2),(GridSize3(1)+2)*(GridSize3(2)+2))	:: A3		! Coefficient matrix for grid 3
	REAL(PRCSN), DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1)									:: Qsource	! Source term for grid 0
	REAL(PRCSN), DIMENSION((GridSize0(1)+2)*(GridSize0(2)+2),(GridSize0(1)+2)*(GridSize0(2)+2))	:: SI		! Semi-Implicit matrix
	!-------------------------------------------------------------------------------------------------------------------!
	! Interfaces
	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE 	tensorC_to_asign_scalar, tensorC_to_asign_tensor, &
							tensorF_to_asign_scalar, tensorF_to_asign_tensor
    END INTERFACE
	!-------------------------------------------------------------------------------------------------------------------!
	
	CONTAINS

	!-------------------------------------------------------------------------------------------------------------------!
	! Interface subroutines	
	SUBROUTINE tensorC_to_asign_scalar(A,b)
		TYPE(Tensor_C), INTENT(INOUT)	:: A(:,:)
		REAL(PRCSN),	INTENT(IN)		:: b
		A%xx = b
		A%xy = b
		A%yy = b
	END SUBROUTINE tensorC_to_asign_scalar

	SUBROUTINE tensorC_to_asign_tensor(A,b)
		TYPE(Tensor_C), INTENT(INOUT)	:: A(:,:)
		TYPE(Tensor_C), INTENT(IN)		:: b(:,:)
		A%xx = b%xx
		A%xy = b%xy
		A%yy = b%yy
	END SUBROUTINE tensorC_to_asign_tensor

	SUBROUTINE tensorF_to_asign_scalar(A,b)
		TYPE(Tensor_F), INTENT(INOUT)	:: A(:,:)
		REAL(PRCSN),	INTENT(IN)		:: b
		A%xx = b
		A%xy = b
		A%yy = b
		A%yx = b
	END SUBROUTINE tensorF_to_asign_scalar

	SUBROUTINE tensorF_to_asign_tensor(A,b)
		TYPE(Tensor_F), INTENT(INOUT)	:: A(:,:)
		TYPE(Tensor_F), INTENT(IN)		:: b(:,:)
		A%xx = b%xx
		A%xy = b%xy
		A%yy = b%yy
		A%yx = b%yx
	END SUBROUTINE tensorF_to_asign_tensor
	!-------------------------------------------------------------------------------------------------------------------!

	SUBROUTINE Cell_To_Face(sz,Cell,Face)
	!----------------------------------------------------------------------------------------------!
	! This subroutine transforms coordinate systems from cell-centered to face-centered
	! The face-centered diffusion tensor values are needed for the discretization stencil, since
	! the derivatives move the data from cell to face
	! The face-centered variable, D_f, goes from 0:N+1,0:M+1 instead of the standard 1:N,1:M. 
	! A ghost zone surrounds the interior values (just repeating the nearest neighbor)
	!
	! The harmonic mean is used to change from cell- to face-centered for xx and yy values
	! However, for the xy and yx values, due to the fact that most are zeros, a straight sum is used
	! Not dividing by 2 allows a 0 and a non-zero cell-centered value to have the non-zero value on the face
	! Harmonic mean does not handle zeros well
	!----------------------------------------------------------------------------------------------!
		INTEGER,		DIMENSION(2),					INTENT( IN) :: sz		! Size of arrays
		TYPE(Tensor_C),	DIMENSION(1:sz(1),1:sz(2)),		INTENT( IN) :: Cell		! Input cell-centered diffusion tensor
		TYPE(Tensor_F),	DIMENSION(0:sz(1)+1,0:sz(2)+1),	INTENT(OUT) :: Face		! Output face-centered diffusion tensor
		TYPE(Tensor_C),	DIMENSION(0:sz(1)+1,0:sz(2)+1)				:: ghost	! Diffusion tensor with added ghost layers
		!------------------------------------------------------------------------------------------!
		! Initialize variables to zero
		Face=0D0	!Face%xx=0D0;  Face%xy=0D0;  Face%yy=0D0; Face%yx=0D0
		ghost=0D0	!ghost%xx=0D0; ghost%xy=0D0; ghost%yy=0D0
		!------------------------------------------------------------------------------------------!
		! Create cell-centered diffusion tensor with ghost layer
		ghost(1:sz(1),1:sz(2))	= Cell				! Interior (1:N, 1:M)
		ghost(0,1:sz(2)) 		= ghost(1,1:sz(1))			! row zero reflects row 1
		ghost(sz(1)+1,1:sz(2)) 	= ghost(sz(1),1:sz(2))		! row N+1 reflects for N
		ghost(0:sz(1)+1,0) 		= ghost(0:sz(1)+1,1)		! column 0 reflects column 1
		ghost(0:sz(1)+1,sz(2)+1)= ghost(0:sz(1)+1,sz(2))	! column M+1 reflects column M
		!------------------------------------------------------------------------------------------!
		! Define xx and yy as a harmonic mean (one knows that xx and yy are all non-zero and positive)
		Face(0:sz(1),0:sz(2)+1)%xx = 2D0/( 1D0/ghost(0:sz(1),0:sz(2)+1)%xx + 1D0/ghost(1:sz(1)+1,0:sz(2)+1)%xx)
		Face(0:sz(1)+1,0:sz(2))%yy = 2D0/( 1D0/ghost(0:sz(1)+1,0:sz(2))%yy + 1D0/ghost(0:sz(1)+1,1:sz(2)+1)%yy)
		!------------------------------------------------------------------------------------------!
		! Define xy and yx. Since these values are usually zero, the harmonic mean is not acceptable
		DO j=0,sz(2)+1; DO i=0,sz(1)+1
			! xy
			IF (i/=sz(1)+1)	THEN
				IF 		(ghost(i,j)%xy == 0D0 .AND. ghost(i+1,j)%xy == 0D0)	THEN	! Both zero
					Face(i,j)%xy = 0D0
				ELSEIF 	(ghost(i,j)%xy /= 0D0 .AND. ghost(i+1,j)%xy /= 0D0)	THEN	! Both non-zer0
					Face(i,j)%xy = 2D0/( 1D0/ghost(i,j)%xy + 1D0/ghost(i+1,j)%xy)
				ELSE																! Mixed
					Face(i,j)%xy = ghost(i,j)%xy + ghost(i+1,j)%xy
				ENDIF
			ENDIF
			
			! yx
			IF (j/=sz(2)+1)	THEN
				IF 		(ghost(i,j)%xy == 0D0 .AND. ghost(i,j+1)%xy == 0D0)	THEN	! Both zero
					Face(i,j)%yx = 0D0
				ELSEIF 	(ghost(i,j)%xy /= 0D0 .AND. ghost(i,j+1)%xy /= 0D0)	THEN	! Both non-zero
					Face(i,j)%yx = 2D0/( 1D0/ghost(i,j)%xy + 1D0/ghost(i,j+1)%xy)
				ELSE																! Mixed
					Face(i,j)%yx = ghost(i,j)%xy + ghost(i,j+1)%xy
				ENDIF
			ENDIF
		ENDDO; ENDDO
		!------------------------------------------------------------------------------------------!
		! Note that each type has an unused row or column, which is left as zeros
		! Row sz(1)+1 for %xx and %yy as well as column sz(2)+1 for %yy and %yx
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE Cell_To_Face	

	SUBROUTINE AnalyticFn(phi)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine specifices the intensity function E (or phi) on a cell-centered mesh
		REAL(PRCSN),DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1),	INTENT(OUT) :: phi	! Outputted solution, does not include ghost cells
		REAL(PRCSN)																:: time, x,y	! position and time variables
	
		phi = 0D0;
		time = tmin + DBLE(tnum)*dt
		
		DO j=0,GridSize0(2)+1; DO i=0,GridSize0(1)+1
			!x = xmin + (DBLE(i)-0.5D0)*GridSpace(1)	! note: (i-0.5)dx = (i-1)dx+dx/2
			!y = ymin + (DBLE(j)-0.5D0)*GridSpace(2)
			x = xmin + DBLE(i)*GridSpace(1)
			y = ymin + DBLE(j)*GridSpace(2)
			!---------------------------------------------------------------------------------------------------------------!
			! Full Taylor Series
			IF     (ProblemType == -3)	THEN	! Manufactured Solution - Full x series - SIN LOG
				phi(i,j) = SIN(PI*x)*LOG(y+1D0)
			ELSEIF (ProblemType == -2) 	THEN	! Manufactured Solution - Full x series - EXP (space and time)
				phi(i,j) = EXP(-x*y*time)
			ELSEIF (ProblemType == -1) 	THEN	! Manufactured Solution - Full x series - EXP
				phi(i,j) = EXP(x*y)
			!---------------------------------------------------------------------------------------------------------------!
			ELSEIF (ProblemType == 0) 	THEN
				phi(i,j) = 0D0
			!---------------------------------------------------------------------------------------------------------------!
			! Powers of x, y
			ELSEIF (ProblemType == 1)	THEN	! Manufactured Solution - Powers of x - x**4
				phi(i,j) = (x**4 + y**4)/12D0
			ELSEIF (ProblemType == 2)	THEN	! Manufactured Solution - Powers of x - x**3
				phi(i,j) = (x**3 + y**3)/6D0
			ELSEIF (ProblemType == 3)	THEN	! Manufactured Solution - Powers of x - x**2
				phi(i,j) = (x**2 + y**2)/2D0
			ELSEIF (ProblemType == 4)	THEN	! Manufactured Solution - Powers of x - x**1.5
				phi(i,j) = (x**1.5D0 + y**1.5D0)*(4D0/3D0)
			ELSEIF (ProblemType == 5)	THEN	! Manufactured Solution - Powers of x - x**1
				phi(i,j) = (x + y)
			ELSEIF (ProblemType == 6)	THEN	! Manufactured Solution - Powers of x - x**0.5
				phi(i,j) = (x**0.5D0 + y**0.5D0)*4D0
			!---------------------------------------------------------------------------------------------------------------!
			! Powers of time
			ELSEIF (ProblemType == 7)	THEN	! Manufactured Solution - Powers of t - t**3
				phi(i,j) = (time**3)/3D0
			ELSEIF (ProblemType == 8)	THEN	! Manufactured Solution - Powers of t - t**2
				phi(i,j) = (time**2)/2D0
			ELSEIF (ProblemType == 9)	THEN	! Manufactured Solution - Powers of t - t**1.5
				phi(i,j) = (time**1.5D0)*(2D0/3D0)
			ELSEIF (ProblemType == 10)	THEN	! Manufactured Solution - Powers of t - t**1
				phi(i,j) = (time)
			ELSEIF (ProblemType == 11)	THEN	! Manufactured Solution - Powers of t - t**0.5
				phi(i,j) = (time**0.5D0)*2D0
			ELSEIF (ProblemType == 12)	THEN	! Manufactured Solution - Powers of t - t**(-0.5)
				phi(i,j) = 1D0/(time**0.5D0)*2D0			
			ELSEIF (ProblemType == 13)	THEN	! Manufactured Solution - Powers of t - t**(-1)
				phi(i,j) = 1D0/time
			!---------------------------------------------------------------------------------------------------------------!
			ELSEIF (ProblemType == 30) THEN		! Analytic solution - exponetial
				phi(i,j) = (tmin/time)*EXP(-x**2/(4.D0*time))*EXP(-y**2/(4.D0*time))				
			ENDIF			
		ENDDO; ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE AnalyticFn

	SUBROUTINE SourceFn(Qsource)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine defines the source function and references the value of phi at the previous time-step
	! The source function is defined as Q in d/dt (phi) + DEL(Flux) = Q
	! This can depend on time and vary with space as well. Default is 0
	!-------------------------------------------------------------------------------------------------------------------!
		REAL(PRCSN),DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1), INTENT(OUT) 	:: Qsource	! Source term
		REAL(PRCSN)																:: x,y,time	! position and time variables
		!TYPE(Tensor_C)															:: D_c

		Qsource = 0.D0
		x=0.D0; y=0.D0; time=0.D0

		IF (MMS == 1) THEN		! Control to turn on/off source loop
			time = tmin + DBLE(tnum)*dt	
			DO j=0,GridSize0(2)+1; DO i=0,GridSize0(1)+1
				!x = xmin + (DBLE(i)-0.5D0)*GridSpace(1)	! note: (i-0.5)dx = (i-1)dx+dx/2
				!y = ymin + (DBLE(j)-0.5D0)*GridSpace(2)
				x = xmin + DBLE(i)*GridSpace(1)
				y = ymin + DBLE(j)*GridSpace(2)
				!---------------------------------------------------------------------------------------------------!
				! Full Taylor Series
				IF     (ProblemType == -3)	THEN	! Manufactured Solution - Full x series - SIN LOG
					Qsource(i,j) = SIN(PI*x)*( k1*PI**2*LOG(y+1D0)+k2/(y+1D0)**2 ) - COS(PI*x)*(2D0*PI*k3)/(y+1D0)
				ELSEIF (ProblemType == -2) 	THEN	! Manufactured Solution - Full x series - EXP (space and time)
					Qsource(i,j) = -EXP(-x*y*time)*( k1*(y*time)**2 + k2*(x*time)**2 + x*y 	&
												   + 2D0*k3*time*(x*y*time-1D0) )
				ELSEIF (ProblemType == -1) 	THEN	! Manufactured Solution - Full x series - EXP
					Qsource(i,j) = -EXP(x*y)*( k1*y**2 + k2*x**2 + 2D0*k3*(1D0+x*y) )
				!---------------------------------------------------------------------------------------------------!
				ELSEIF (ProblemType == 0) 	THEN
					Qsource(i,j) = 0D0
				!---------------------------------------------------------------------------------------------------!
				! Powers of x, y
				ELSEIF (ProblemType == 1)	THEN	! Manufactured Solution - Powers of x - x**4
					Qsource(i,j) = -(k1*x**2 + k2*y**2)
				ELSEIF (ProblemType == 2)	THEN	! Manufactured Solution - Powers of x - x**3
					Qsource(i,j) = -(k1*x + k2*y)
				ELSEIF (ProblemType == 3)	THEN	! Manufactured Solution - Powers of x - x**2
					Qsource(i,j) = -(k1 + k2)
				ELSEIF (ProblemType == 4)	THEN	! Manufactured Solution - Powers of x - x**1.5
					Qsource(i,j) = -(k1/SQRT(x) + k2/SQRT(y))
				ELSEIF (ProblemType == 5)	THEN	! Manufactured Solution - Powers of x - x**1
					Qsource(i,j) = 0D0
				ELSEIF (ProblemType == 6)	THEN	! Manufactured Solution - Powers of x - x**0.5
					Qsource(i,j) = (k1/x**1.5D0 + k2/y**1.5D0)
				!---------------------------------------------------------------------------------------------------!
				! Powers of time
				ELSEIF (ProblemType == 7)	THEN	! Manufactured Solution - Powers of t - t**3
					Qsource(i,j) = time**2
				ELSEIF (ProblemType == 8)	THEN	! Manufactured Solution - Powers of t - t**2
					Qsource(i,j) = time
				ELSEIF (ProblemType == 9)	THEN	! Manufactured Solution - Powers of t - t**1.5
					Qsource(i,j) = SQRT(time)
				ELSEIF (ProblemType == 10)	THEN	! Manufactured Solution - Powers of t - t**1
					Qsource(i,j) = 1D0
				ELSEIF (ProblemType == 11)	THEN	! Manufactured Solution - Powers of t - t**0.5
					Qsource(i,j) = 1D0/SQRT(time)
				ELSEIF (ProblemType == 12)	THEN	! Manufactured Solution - Powers of t - t**(-0.5)
					Qsource(i,j) = -1D0/time**1.5D0
				ELSEIF (ProblemType == 13)	THEN	! Manufactured Solution - Powers of t - t**(-1)
					Qsource(i,j) = -1D0/time**2
				!---------------------------------------------------------------------------------------------------!
				ENDIF
			ENDDO; ENDDO
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SourceFn
	
	SUBROUTINE DiffusionTensor(sz,CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine defines the diffusion tensor on the cell-centers to overlay with the initial concentration field
	! This version uses a linear interface (y=mx+b) to separate two materials, one below the line and the other above it.
	! NOTE: In the future, the Diffusion Tensor may be initialized from an outside source, say the material properties
	! 		are recovered from a volume-of-fluid, moment-of-fluid, or level set method.
	!-------------------------------------------------------------------------------------------------------------------!
		INTEGER, 		DIMENSION(2),			INTENT( IN)	:: sz				! Incoming matrix size (sz = size without vowels)
		TYPE(Tensor_C), DIMENSION(sz(1),sz(2)),	INTENT(OUT) :: CellCentered		! Outgoing diffusion tensor
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		CellCentered%xx = 0D0; CellCentered%yy = 0D0; CellCentered%xy = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Specify diffusion tensor as a function of position
		!IF (MixedMethod == 0)	THEN
			CellCentered%xx = k1
			CellCentered%yy = k2
			CellCentered%xy = k3			
		!ELSE !IF (MixedMethod /= 0) CALL MixDiffusionTensor(CellCentered)
		!	CALL MixChannels(CellCentered)				! Domain split into 4 channels with 8 lines
		!ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE DiffusionTensor

	SUBROUTINE DiffusionOperator(sz,L,state)
	!----------------------------------------------------------------------------------------------!
	! This outputs the coefficient matrix that discretizes the diffusion equation
	! The output is A = I-dt*L, where L = DIVERGENCE(D*GRADIENT(_))
	! When running semi-implicit, this becomes A=I-0.5*dt*L
	! This matrix is made in 1 step, so the identity matrix is imbedded naturally, and we multiply L by -dt
	! Boundary conditions are applied where they appear in the discretization of the diffuson operator
	! Ghost cell values are set when the RHS is defined
	!
	! There is a 3x3 stencil for 9 total points, centered at i,j and going i+/-{0,1} and j+/-{0,1}
	! The variable alpha(i,j,k) keeps track of all 9 coefficients for mesh point i,j.
	! The following layout correlates to the ordering of alpha:
	!	(i+1,j-1)	(i+1,j)		(i+1,j+1)	| alpha(i,j,6)	alpha(i,j,2)	alpha(i,j,7)
	!	(i,j-1)		(i,j)		(i,j+1)		| alpha(i,j,4)	alpha(i,j,1)	alpha(i,j,5)
	!	(i-1,j-1)	(i-1,j)		(i-1,j+1)	| alpha(i,j,8)	alpha(i,j,3)	alpha(i,j,9)
	!
	! The grid spacing (dx or dy) is calculated by Sys_Length/(sz+1), where both of those are 2 element vectors
	! So dx=Sys_Length(1)/(sz(1)+1) and dy=Sys_Length(2)/(sz(2)+1)
		INTEGER, 		DIMENSION(2),										INTENT( IN)	:: sz			! Incoming matrix size (sz = size without vowels)
		REAL(PRCSN), 	DIMENSION((sz(1)+2)*(sz(2)+2),(sz(1)+2)*(sz(2)+2)),	INTENT(OUT) :: L			! Output diffusion operator coefficient matrix 
		INTEGER, OPTIONAL,													INTENT( IN) :: state		! Optional output mode for Semi-Implicit RHS matrix
		TYPE(Tensor_C), DIMENSION(sz(1),sz(2))											:: D_cell		! Cell-centered diffusion tensor
		TYPE(Tensor_F), DIMENSION(0:sz(1)+1,0:sz(2)+1)									:: Diff			! Face-centered diffusion tensor
		REAL(PRCSN),	DIMENSION(9)													:: alpha		! Constant to keep track of each stencil's element's contribution
		REAL(PRCSN)																		:: dx, dy		! Grid spacing coefficients
		REAL(PRCSN)																		:: hx, hy, hxy	! Inverse grid spacing coefficients presistant throughout discretization
		REAL(PRCSN)																		:: Gamma		! Constant containing BC terms (like p,q,s)
		REAL(PRCSN)																		:: p, q			! Boundary condition variables: p*phi + q*d(phi)/d(n) = s
		
		! Initialize variables to 0
		L = 0D0; alpha = 0D0; Diff = 0D0
		hx = 0D0; hy = 0D0; hxy = 0D0
		Gamma=0D0; p=0D0; q=0D0

		! Define grid spacing parameters
		dx = Sys_Length(1)/DBLE(sz(1)+1)
		dy = Sys_Length(2)/DBLE(sz(2)+1)
		hx  = 1D0/dx**2		! 1/dx**2
		hy  = 1D0/dy**2		! 1/dy**2
		hxy = 0.25D0/dx/dy	! 1/(4*dx*dy)
	
		! Define diffusion tensor (cell-centered first, then transfer to face-centered)
		CALL DiffusionTensor(sz,D_cell)
		CALL Cell_To_Face(sz,D_cell,Diff)

		! Iterate through all points to create diffusion operator L
		DO j=0,sz(2)+1; DO i=0,sz(1)+1
			! Define global position (maps i,j components from a matrix to a vector)
			i0 = 1+i+j*(sz(1)+2)
			IF (i/=0 .AND. i/=sz(1)+1 .AND. j/=0 .AND. j/=sz(2)+1)	THEN	! Interior
				!----------------------------------------------------------------------------------------------!
				! Define alpha coefficients in terms of mesh spacing and diffusion tensor
				alpha(1) =-hx*(Diff(i,j)%xx + Diff(i-1,j)%xx) - hy*(Diff(i,j)%yy + Diff(i,j-1)%yy)	! i,j
				alpha(2) = hx*Diff(i,j)%xx 	 + hxy*(Diff(i+1,j-1)%yx-Diff(i+1,j)%yx)	! i+1,j
				alpha(3) = hx*Diff(i-1,j)%xx + hxy*(Diff(i-1,j)%yx-Diff(i-1,j-1)%yx)	! i-1,j
				alpha(4) = hy*Diff(i,j-1)%yy + hxy*(Diff(i,j-1)%xy-Diff(i-1,j-1)%xy)	! i,j-1
				alpha(5) = hy*Diff(i,j)%yy	 + hxy*(Diff(i-1,j+1)%xy-Diff(i,j+1)%xy)	! i,j+1
				alpha(6) =-hxy*(Diff(i,j-1)%xy 	 + Diff(i+1,j-1)%yx)					! i+1,j-1
				alpha(7) = hxy*(Diff(i,j+1)%xy	 + Diff(i+1,j)%yx)						! i+1,j+1
				alpha(8) = hxy*(Diff(i-1,j-1)%xy + Diff(i-1,j-1)%yx)					! i-1,j-1
				alpha(9) =-hxy*(Diff(i-1,j+1)%xy + Diff(i-1,j)%yx)						! i-1,j+1
				!----------------------------------------------------------------------------------------------!
				alpha = -dt*alpha ! This helps in forming the coefficient matrix A=I-dt*L 
				IF 	(TypeMethod == 'SemiImplicit')	alpha = 0.5D0*alpha	! (A=I-0.5D0*dt*L)
				IF (PRESENT(state)) alpha = -alpha	! A=I+0.5D0*dt*L for semi-implicit operator
				!----------------------------------------------------------------------------------------------!
				! Populate coefficient matrix L
				L(i0,i0) 			= alpha(1) + 1D0	! i,j	! This imbeds the indentity matrix
				L(i0,i0+1)			= alpha(2)			! i+1,j
				L(i0,i0-1)			= alpha(3)			! i-1,j
				L(i0,i0-sz(1)-2)	= alpha(4)			! i,j-1
				L(i0,i0+sz(1)+2)	= alpha(5)			! i,j+1
				L(i0,i0+1-sz(1)-2)	= alpha(6)			! i+1,j-1
				L(i0,i0+1+sz(1)+2)	= alpha(7)			! i+1,j+1
				L(i0,i0-1-sz(1)-2)	= alpha(8)			! i-1,j-1
				L(i0,i0-1+sz(1)+2)	= alpha(9)			! i-1,j+1
				!----------------------------------------------------------------------------------------------!
			ELSE	! Boundaries
				L(i0,i0) = 1D0

				! Left/Right sides
				IF ((BC_alpha(1)/=0 .OR. BC_beta(1)/=0) .AND. (BC_alpha(2)/=0 .OR. BC_beta(2)/=0)) THEN
					! Non-periodic (Robin)
					IF (i==0)		THEN	! Left side (i=0)
						p = BC_alpha(1); q = BC_beta(1)
						Gamma = q/dx
						IF (j/=0 .AND. j/=sz(2)+1)	THEN
							L(i0,i0+1) = -Gamma/(p+Gamma)	! phi(0,j)=(s+Gamma*phi(1,j))/(p+Gamma)
						ELSE
							L(i0,i0+1) = -Gamma/(p+Gamma)*5D-1
						ENDIF
					ENDIF
	
					IF (i==sz(1)+1)	THEN	! Right side (i=N+1)
						p = BC_alpha(2); q = BC_beta(2)
						Gamma = q/dx
						IF (j/=0 .AND. j/=sz(2)+1)	THEN
							L(i0,i0-1) = -Gamma/(p+Gamma)	! phi(N+1,j)=(s+Gamma*phi(N,j))/(p+Gamma)
						ELSE
							L(i0,i0-1) = -Gamma/(p+Gamma)*5D-1
						ENDIF
					ENDIF
				ELSE	
					! Periodic
					IF (j/=0 .AND. j/=sz(2)+1)	THEN
						IF (i==0) 		L(i0,i0+sz(1)) = -1D0		! 0,j 	=> N,j	! i0+N
						IF (i==sz(1)+1) L(i0,i0-sz(1)) = -1D0		! N+1,j	=> 1,j	! i0-N
					ENDIF
				ENDIF

				! Top/Bottom sides
				IF ((BC_alpha(3)/=0 .OR. BC_beta(3)/=0) .AND. (BC_alpha(4)/=0 .OR. BC_beta(4)/=0)) THEN
					! Non-periodic (Robin)
					IF (j==sz(2)+1)	THEN	! Top side (j=M+1)
						p = BC_alpha(3); q = BC_beta(3)
						Gamma = q/dy
						IF (i/=0 .AND. i/=sz(1)+1)	THEN
							L(i0,i0-sz(1)-2) = -Gamma/(p+Gamma)	! phi(i,M+1)=(s+Gamma*phi(i,M))/(p+Gamma)
						ELSE
							L(i0,i0-sz(1)-2) = -Gamma/(p+Gamma)*5D-1
						ENDIF
					ENDIF
	
					IF (j==0)		THEN	! Bottom side (j=0)
						p = BC_alpha(4); q = BC_beta(4)
						Gamma = q/dy
						IF (i/=0 .AND. i/=sz(1)+1)	THEN
							L(i0,i0+sz(1)+2) = -Gamma/(p+Gamma)	! phi(i,0)=(s+Gamma*phi(i,1))/(p+Gamma)
						ELSE
							L(i0,i0+sz(1)+2) = -Gamma/(p+Gamma)*5D-1
						ENDIF
					ENDIF
				ELSE	! Periodic
					IF (i/=0 .AND. i/=sz(1)+1)	THEN
						IF (j==0) 		L(i0,i0+sz(2)*(sz(1)+2)) = -1D0	! i,0 	=> i,M	! i0+M*(N+2)
						IF (j==sz(2)+1) L(i0,i0-sz(2)*(sz(1)+2)) = -1D0	! i,M+1	=> i,1	! i0-M*(N+2)
					ENDIF
				ENDIF
				
				! Handle corners when all sides are periodic
				IF (SUM(ABS(BC_alpha))+SUM(ABS(BC_beta))==0)	THEN
					IF (i==0 .AND. j==0)				L(i0,1+sz(1)+sz(2)*(sz(1)+2)) = -1D0	! 0,0 		=> N,M
					IF (i==0 .AND. j==sz(2)+1)			L(i0,1+sz(1)+1*(sz(1)+2)) = -1D0		! 0,M+1		=> N,1
					IF (i==sz(1)+1 .AND. j==0)			L(i0,1+1+sz(2)*(sz(1)+2)) = -1D0		! N+1,0		=> 1,M
					IF (i==sz(1)+1 .AND. j==sz(2)+1)	L(i0,1+1+1*(sz(1)+2)) = -1D0			! N+1,M+1	=> 1,1
				ENDIF
			ENDIF
			!----------------------------------------------------------------------------------------------!
		END DO; END DO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE DiffusionOperator

	SUBROUTINE RightHandSide(phi_old,RHS)
	!----------------------------------------------------------------------------------------------!
	! This calculates the RHS of the PDE discertization, and then applies boundary conditions
	! The source matrix (Q) and the semi-implicit matrix (SI) are pre-defined and shared in this module
	!----------------------------------------------------------------------------------------------!
		REAL(PRCSN),	DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1),	INTENT( IN) :: phi_old
		REAL(PRCSN),	DIMENSION(GridSize0(1)+2,GridSize0(2)+2),		INTENT(OUT)	:: RHS
		REAL(PRCSN),	DIMENSION((GridSize0(1)+2)*(GridSize0(2)+2))				:: b
		REAL(PRCSN),	DIMENSION((GridSize0(1)+2)*(GridSize0(2)+2))				:: DumVec
		REAL(PRCSN),	DIMENSION((GridSize0(1)+2),(GridSize0(2)+2))				:: DumMat
		REAL(PRCSN),	DIMENSION(0:GridSize0(1)+1,0:GridSize0(2)+1)				:: phi_analytic
		REAL(PRCSN)																	:: dx,dy,p,q,s,Gamma
		!------------------------------------------------------------------------------------------!
		! Initialize Variabls
		RHS=0D0; phi_analytic=0D0; DumVec=0D0; b=0D0
		dx=0D0; dy=0D0; p=0D0; q=0D0; s=0D0; Gamma=0D0
		IF (BC_fun==1) CALL AnalyticFn(phi_analytic)
		!------------------------------------------------------------------------------------------!
		! Calculate RHS
		IF 	(TypeMethod == 'Implicit') THEN
			b = PACK(phi_old+dt*Qsource,.TRUE.)
		ELSEIF 	(TypeMethod == 'SemiImplicit') THEN
			! Done in several steps which may not be necessary
			DumVec = PACK(phi_old,.TRUE.)
			b = MATMUL(SI,DumVec) + dt*PACK(Qsource,.TRUE.)
		ENDIF
		!------------------------------------------------------------------------------------------!
		! Apply BC to RHS
		dx = Sys_Length(1)/DBLE(GridSize0(1)+1)
		dy = Sys_Length(2)/DBLE(GridSize0(2)+1)
		DO j=0,GridSize0(2)+1; DO i=0,GridSize0(1)+1
			i0 = 1+i+j*(GridSize0(1)+2)
			IF (i==0)		THEN
				p = BC_alpha(1); q = BC_beta(1); s = BC_psi(1)
				IF (BC_fun==1) s = phi_analytic(i,j)
				Gamma = q/dx
				IF (p/=0 .OR. q/=0) THEN
					b(i0) = s/(p+Gamma)
					! If on corner points and y is not periodic
					IF ((BC_alpha(3)/=0 .OR. BC_beta(3)/=0) .AND. (BC_alpha(4)/=0 .OR. BC_beta(4)/=0))	THEN
						IF (j==0 .OR. j==GridSize0(2)+1) b(i0) = b(i0)*5D-1
					ENDIF
				ELSE
					b(i0) = 0D0
				ENDIF
			ENDIF
			
			IF (i==GridSize0(1)+1)	THEN
				p = BC_alpha(2); q = BC_beta(2); s = BC_psi(2)
				IF (BC_fun==1) s = phi_analytic(i,j)
				Gamma = q/dx
				IF (p/=0 .OR. q/=0) THEN
					b(i0) = s/(p+Gamma)
					! If on corner points and y is not periodic
					IF ((BC_alpha(3)/=0 .OR. BC_beta(3)/=0) .AND. (BC_alpha(4)/=0 .OR. BC_beta(4)/=0)) THEN 
						IF (j==0 .OR. j==GridSize0(2)+1) b(i0) = b(i0)*5D-1
					ENDIF
				ELSE
					b(i0) = 0D0
				ENDIF
			ENDIF
			
			IF (j==GridSize0(2)+1)	THEN
				p = BC_alpha(3); q = BC_beta(3); s = BC_psi(3)
				IF (BC_fun==1) s = phi_analytic(i,j)
				Gamma = q/dy
				IF (p/=0 .OR. q/=0)  THEN
					IF (i/=0 .AND. i/=GridSize0(1)+1) THEN
						b(i0) = s/(p+Gamma)
					ELSE
						IF ((BC_alpha(1)/=0 .OR. BC_beta(1)/=0) .AND. (BC_alpha(2)/=0 .OR. BC_beta(2)/=0)) THEN 
							! x not periodic
							b(i0) = b(i0) + s/(p+Gamma)*5D-1
						ELSE
							! x is periodic
							b(i0) = s/(p+Gamma)
						ENDIF
					ENDIF
				ELSE
					IF (i/=0 .AND. i/=GridSize0(1)+1) THEN
						b(i0)=0D0
					ELSE
						! Only set b(i0) to zero if x is also periodic
						IF (SUM(ABS(BC_alpha(1:2)))+SUM(ABS(BC_beta(1:2)))==0) b(i0)=0D0
					ENDIF
				ENDIF
			ENDIF
			
			IF (j==0)		THEN
				p = BC_alpha(4); q = BC_beta(4); s = BC_psi(4)
				IF (BC_fun==1) s = phi_analytic(i,j)
				Gamma = q/dy
				IF (p/=0 .OR. q/=0)  THEN
					IF (i/=0 .AND. i/=GridSize0(1)+1) THEN
						b(i0) = s/(p+Gamma)
					ELSE
						IF ((BC_alpha(1)/=0 .OR. BC_beta(1)/=0) .AND. (BC_alpha(2)/=0 .OR. BC_beta(2)/=0)) THEN 
							! x not periodic
							b(i0) = b(i0) + s/(p+Gamma)*5D-1
						ELSE
							! x is periodic
							b(i0) = s/(p+Gamma)
						ENDIF
					ENDIF
				ELSE
					IF (i/=0 .AND. i/=GridSize0(1)+1) THEN
						b(i0)=0D0
					ELSE
						! Only set b(i0) to zero if x is also periodic
						IF (SUM(ABS(BC_alpha(1:2)))+SUM(ABS(BC_beta(1:2)))==0) b(i0)=0D0
					ENDIF
				ENDIF
			ENDIF
		ENDDO; ENDDO
		!------------------------------------------------------------------------------------------!
		! Reform matrix into 
		RHS = RESHAPE(b,GridSize0+2)
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE RightHandSide
	
	SUBROUTINE Smooth(sz,Niter,sqx,sqb)
	!----------------------------------------------------------------------------------------------!
	! The subroutine sets up the matrix equation, A*x=b, and solves it to find x
	! By solving the system on different grids, this is effectively 'smoothing' the error
	! Subroutine inputs the current residual and the previously found dE (to be used as a guess)
	! The output is the newly found dE (so dE is an inout variable)
	!----------------------------------------------------------------------------------------------!
	! The matrix A is a globally shared variable already determined. Depending on the size of the system
	! (determined by variable sz), the appropriate coefficient matrix A will be used
	! The inout variable sqx is the square (i.e. matrix form) of x, which is used as an initial guess
	! for the solver, and then reused as the solution to the solver
	! sqb is the square b, or the RHS of the equation of interest.
	! Hence, I have A,b, and an initial guess for x, all I need to solve the system
	! The variable 'state' determines how many sweeps to preform, depending on grid level
	!----------------------------------------------------------------------------------------------!
	! For the fine grid (sz=GridSize0), A=A0, b=RHS, and x=phi_new
	! For any fine grid, A=A1 or A2 or A3, b=Res (i.e. residual) and x=delta_phi
	!----------------------------------------------------------------------------------------------!
		INTEGER,		DIMENSION(2),					INTENT( IN)		:: sz			! Input size of matricies
		INTEGER,										INTENT( IN) 	:: Niter		! Number of smoothing iterations
		REAL(PRCSN), 	DIMENSION(sz(1)+2,sz(2)+2),		INTENT( IN)		:: sqb			! Input RESidual. NOTE: in fine grid case, RES = phi_old
		REAL(PRCSN), 	DIMENSION(0:sz(1)+1,0:sz(2)+1),	INTENT(INOUT)	:: sqx			! Square x. Used as initial guess, then set as solution to Ax=b
		REAL(PRCSN), 	DIMENSION((sz(1)+2)*(sz(2)+2))					:: x, b, guess	! vectors
		REAL(PRCSN), 	DIMENSION(sz(1)+2,sz(2)+2)						:: DumMat		! dummy matrix
		!------------------------------------------------------------------------------------------!
		! Initialize variables to zero
		x=0D0; b=0D0; guess=0D0
		! Initialize vectors from input square (matrix) format
		guess = PACK(sqx,.TRUE.)
		b 	  = PACK(sqb,.TRUE.)
		!------------------------------------------------------------------------------------------!
		! Solve the matrix equation
		IF 		(ALL(sz==GridSize0)) 	THEN
			IF (SolverType=='SGS')  CALL SGS(sz+2,Niter,A0,b,x,guess)
			IF (SolverType=='RBGS') CALL RBGS(sz+2,Niter,A0,b,x,guess)
		ELSEIF 	(ALL(sz==GridSize1)) 	THEN
			IF (SolverType=='SGS')  CALL SGS(sz+2,Niter,A1,b,x,guess)
			IF (SolverType=='RBGS') CALL RBGS(sz+2,Niter,A1,b,x,guess)
		ELSEIF 	(ALL(sz==GridSize2)) 	THEN
			IF (SolverType=='SGS')  CALL SGS(sz+2,Niter,A2,b,x,guess)
			IF (SolverType=='RBGS') CALL RBGS(sz+2,Niter,A2,b,x,guess)
		ELSEIF 	(ALL(sz==GridSize3)) 	THEN
			IF (SolverType=='SGS')  CALL SGS(sz+2,Niter,A3,b,x,guess)
			IF (SolverType=='RBGS') CALL RBGS(sz+2,Niter,A3,b,x,guess)
		ENDIF
		!------------------------------------------------------------------------------------------!
		! Output solution (x) in matrix form
		DumMat = RESHAPE(x,sz+2)			! Put in standard matrix form
		sqx(0:sz(1)+1,0:sz(2)+1)=DumMat		! Put in matrix form with ghost cells
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE Smooth
!-----------------------------------------------------------------------------------------------------------------------!
END MODULE MG_Subroutines
