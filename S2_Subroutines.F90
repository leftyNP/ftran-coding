! SOM2_Subroutines.F90
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

MODULE Types
! This module introduces the derived data types
    USE DefineParameters    ! Note: all variables are set in this module, which must be compiled first

    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!    
    ! Derived Data Types
    TYPE Face_Centered
	! This derived data type represents a variable which is stored on the face-center of a cell
	! By keeping the x-direction and the y-direction together in this way, all the face-centered 
	! values are kept together as a single variable
	! An extension to 3D would simply add a zface to this type
        REAL(PRCSN), DIMENSION(xmax+1,ymax)   :: xface
        REAL(PRCSN), DIMENSION(xmax,ymax+1)   :: yface
    END TYPE Face_Centered
    
    TYPE Tensor_C
    ! Cell-centered tensor matrix
    ! Note the built in assumption that yx = xy
    ! If this deviates, the change would also necisitate proper propagation
        REAL(PRCSN), DIMENSION(xmax,ymax)  	:: xx, xy, yy
    END TYPE Tensor_C

    TYPE Tensor_F
    ! Face-centered tensor matrix. There are two options, iso = isotropic; aniso = anisotropic
    ! This leads to four terms:
    !	xx - K%xx_yy%xface
    !	yy - K%xx_yy%yface
    !	xy - K%xy_yx%xface
    !	yx - K%xy_yx%yface
    	TYPE(Face_Centered)	:: xx_yy, xy_yx
    END TYPE Tensor_F
    
    TYPE Cell_Zone
    ! Zonal matrix system for each cell of the matrix. Each zonal matrix is a 4x4 system,
    ! accouting for the left, right, top, bottom value of the concentration.
    ! At first, the cell-centered value is included, making this a 5x5 system.
    ! BC are applied to this, and then the cell-center value is solved for, allowing this
    ! system to then be reduced to a 4x4 matrix for each cell. 
    	REAL(PRCSN), DIMENSION(4,4,xmax,ymax)	:: cz
    END TYPE Cell_Zone
        
    TYPE Global_Zone
    ! Global matrix system that contains the contributions from each zone
    ! Each row represents one face of the system, where a face is a particular wall of a cell.
    	!REAL(PRCSN), DIMENSION((xmax+1)*ymax+xmax*(ymax+1),(xmax+1)*ymax+xmax*(ymax+1)) :: gz
    	REAL(PRCSN), DIMENSION(VecLength,VecLength) :: gz
    END TYPE Global_Zone
    	
    TYPE Sparse_Matrix
    ! This is a sparse matrix format to represent the global matrix system
    ! In 2D, it is a 7 banded matrix, meaning a vector of 8*VecLength will be sufficient
	! The abbreviation, sgz, is in homage to the Global_Zone call of %gz.
    	REAL(PRCSN), DIMENSION(8*VecLength,3) :: sgz
    END TYPE Sparse_Matrix    
    !-------------------------------------------------------------------------------------------------------------------!    
END MODULE Types

MODULE SOM_Utilities
!-----------------------------------------------------------------------------------------------------------------------!
! This module contains the subroutines called in the program SOM that are not physics related
! These utility subroutines include how to name a file, save a file, or perform repetative simple mathematical functions
!-----------------------------------------------------------------------------------------------------------------------!
    USE DefineParameters
    USE Types

	IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	CHARACTER(LEN=12)   :: filename
	!-------------------------------------------------------------------------------------------------------------------!
	! Interfaces must be specified for an overloaded function

	INTERFACE Invert4x4
	! To invert the array of 4x4 matricies, each cell's matrix must be inverted.
	! This is done for a (4,4) size, a (4,4,xmax,ymax) size, and the derived data type "Cell_Zone", (5,5,xmax,ymax)
		MODULE PROCEDURE Invert4x4_Small, Invert4x4Big_Exact_Size, Invert4x4Big_Cell_Zone
	END INTERFACE
	! Turns out I never call any of these because I don't invert a matrix any more

	INTERFACE Conjugate_Gradient
	! Conjugate Gradient Algorithm with need either a sparse of dense matrix input
		MODULE PROCEDURE Conjugate_Gradient_Full,Conjugate_Gradient_Sparse
	END INTERFACE
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
		INTEGER, 							INTENT( IN) :: k			! Time step index
		CHARACTER(LEN=1),					INTENT( IN) :: LTTR 		! LTTR = LeTTeR = first letter of filename
        REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT (IN) :: u			! input matrix
		REAL(PRCSN), DIMENSION(xmax,ymax)    			:: u_restricted	! altered matrix to be saved to file
		CHARACTER (LEN=20)                          	:: frmt			! format of write statement

        WRITE(frmt,'(A1,I4,A13)') '(',ymax,'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))'
	    ! Meaning: ymax gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that
	    
        u_restricted = u      ! Leave input matrix unchanged while modifying values for saving
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
	    
		CALL Namefile(k,LTTR,filename)							! Call subroutine to generate sequential file name

		OPEN (unit=1, file=filename)							! Open said file

		WRITE(1, '(A7,D15.7,TR5,A7,I5/)') '% dt = ', dt, 'tnum = ',tnum   ! Set first line to give time-step

        WRITE(1,frmt) (u_restricted(i,:),i=1,xmax)

		CLOSE (1)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Savefile

	SUBROUTINE SaveGeneral(xsize,ysize,filename,matrix)
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
		INTEGER, 						    INTENT (IN) :: xsize,ysize
		CHARACTER (LEN=*),                  INTENT (IN) :: filename
        REAL(PRCSN), DIMENSION(xsize,ysize),INTENT (IN) :: matrix
		REAL(PRCSN), DIMENSION(xsize,ysize)    		    :: u_restricted
		CHARACTER (LEN=20)                              :: frmt

        WRITE(frmt,'(A1,I6,A13)') '(',ysize,'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))
	    ! Meaning: ymax gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that
	    
        u_restricted = matrix
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
		OPEN (unit=1, file=filename)                            ! Open said file

        WRITE(1,frmt) (u_restricted(icell,:),icell=1,xsize)
       
		CLOSE (1)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SaveGeneral
	
	SUBROUTINE Invert4x4_Small(A,B)
	!-------------------------------------------------------------------------------------------------------------------!
	! Calculated the inverse of a 4x4 matrix exactly. Due to the small size of this matrix, the inverse can be done
	! by a finite number of algebra steps. An iterative solution is not necessary.
	! The steps used to calculate this inverse are from a Maple output
        REAL(PRCSN), DIMENSION(4,4), INTENT (IN) :: A
        REAL(PRCSN), DIMENSION(4,4), INTENT(OUT) :: B
		REAL(PRCSN)								 :: DetA 	! Determinant of A
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		B 	 = 0D0
		DetA = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate the determinant of A
		DetA =   ( A(1,1)*A(2,2)*A(3,3)*A(4,4) + A(1,1)*A(3,2)*A(2,4)*A(4,3) + A(1,1)*A(4,2)*A(2,3)*A(3,4)	&
				 + A(2,1)*A(1,2)*A(3,4)*A(4,3) + A(2,1)*A(3,2)*A(1,3)*A(4,4) + A(2,1)*A(4,2)*A(1,4)*A(3,3)	&
				 + A(3,1)*A(1,2)*A(2,3)*A(4,4) + A(3,1)*A(2,2)*A(1,4)*A(4,3) + A(3,1)*A(4,2)*A(1,3)*A(2,4)	&
				 + A(4,1)*A(1,2)*A(2,4)*A(3,3) + A(4,1)*A(2,2)*A(1,3)*A(3,4) + A(4,1)*A(3,2)*A(1,4)*A(2,3)	&
				 - A(1,1)*A(2,2)*A(3,4)*A(4,3) - A(1,1)*A(3,2)*A(2,3)*A(4,4) - A(1,1)*A(4,2)*A(2,4)*A(3,3)	&
				 - A(2,1)*A(1,2)*A(3,3)*A(4,4) - A(2,1)*A(3,2)*A(1,4)*A(4,3) - A(2,1)*A(4,2)*A(1,3)*A(3,4)	&
				 - A(3,1)*A(1,2)*A(2,4)*A(4,3) - A(3,1)*A(2,2)*A(1,3)*A(4,4) - A(3,1)*A(4,2)*A(1,4)*A(2,3)	&
				 - A(4,1)*A(1,2)*A(2,3)*A(3,4) - A(4,1)*A(2,2)*A(1,4)*A(3,3) - A(4,1)*A(3,2)*A(1,3)*A(2,4) )

		IF (ABS(DetA) < 1.D-10) PRINT '(/,TR4,A50)', 'Error: Cannot invert matrix since determinant is 0'
		
		! Calculate elements of matrix B
		B(1,1) =  A(2,2)*A(3,3)*A(4,4) + A(3,2)*A(2,4)*A(4,3) + A(4,2)*A(2,3)*A(3,4) 		&
			    - A(4,2)*A(2,4)*A(3,3) - A(2,2)*A(3,4)*A(4,3) - A(3,2)*A(2,3)*A(4,4)

		B(1,2) =  A(1,2)*A(3,4)*A(4,3) + A(3,2)*A(1,3)*A(4,4) + A(4,2)*A(1,4)*A(3,3)		&
 				- A(1,2)*A(3,3)*A(4,4) - A(3,2)*A(1,4)*A(4,3) - A(4,2)*A(1,3)*A(3,4)

		B(1,3) =  A(1,2)*A(2,3)*A(4,4) + A(2,2)*A(1,4)*A(4,3) + A(4,2)*A(1,3)*A(2,4)		&
				- A(4,2)*A(1,4)*A(2,3) - A(1,2)*A(2,4)*A(4,3) - A(2,2)*A(1,3)*A(4,4)

		B(1,4) =  A(1,2)*A(2,4)*A(3,3) + A(2,2)*A(1,3)*A(3,4) + A(3,2)*A(1,4)*A(2,3) 		&
				- A(1,2)*A(2,3)*A(3,4) - A(2,2)*A(1,4)*A(3,3) - A(3,2)*A(1,3)*A(2,4)

		B(2,1) =  A(2,1)*A(3,4)*A(4,3) + A(3,1)*A(2,3)*A(4,4) + A(4,1)*A(2,4)*A(3,3) 		&
				- A(2,1)*A(3,3)*A(4,4) - A(3,1)*A(2,4)*A(4,3) - A(4,1)*A(2,3)*A(3,4)

		B(2,2) =  A(1,1)*A(3,3)*A(4,4) + A(3,1)*A(1,4)*A(4,3) + A(4,1)*A(1,3)*A(3,4) 		&
				- A(4,1)*A(1,4)*A(3,3) - A(1,1)*A(3,4)*A(4,3) - A(3,1)*A(1,3)*A(4,4)

		B(2,3) =  A(1,1)*A(2,4)*A(4,3) + A(2,1)*A(1,3)*A(4,4) + A(4,1)*A(1,4)*A(2,3) 		&
				- A(1,1)*A(2,3)*A(4,4) - A(2,1)*A(1,4)*A(4,3) - A(4,1)*A(1,3)*A(2,4) 
		
		B(2,4) =  A(1,1)*A(2,3)*A(3,4) + A(2,1)*A(1,4)*A(3,3) + A(3,1)*A(1,3)*A(2,4)		&
				- A(3,1)*A(1,4)*A(2,3) - A(1,1)*A(2,4)*A(3,3) - A(2,1)*A(1,3)*A(3,4)

		B(3,1) =  A(2,1)*A(3,2)*A(4,4) + A(3,1)*A(2,4)*A(4,2) + A(4,1)*A(2,2)*A(3,4)		&
				- A(2,1)*A(3,4)*A(4,2) - A(3,1)*A(2,2)*A(4,4) - A(4,1)*A(2,4)*A(3,2) 
				
		B(3,2) =  A(1,1)*A(3,4)*A(4,2) + A(3,1)*A(1,2)*A(4,4) + A(4,1)*A(1,4)*A(3,2) 		&
				- A(1,1)*A(3,2)*A(4,4) - A(3,1)*A(1,4)*A(4,2) - A(4,1)*A(1,2)*A(3,4)
		
		B(3,3) =  A(1,1)*A(2,2)*A(4,4) + A(2,1)*A(1,4)*A(4,2) + A(4,1)*A(1,2)*A(2,4)		&
				- A(1,1)*A(2,4)*A(4,2) - A(2,1)*A(1,2)*A(4,4) - A(4,1)*A(1,4)*A(2,2) 
				
		B(3,4) =  A(1,1)*A(2,4)*A(3,2) + A(2,1)*A(1,2)*A(3,4) + A(3,1)*A(1,4)*A(2,2) 		&
				- A(1,1)*A(2,2)*A(3,4) - A(2,1)*A(1,4)*A(3,2) - A(3,1)*A(1,2)*A(2,4)
				
		B(4,1) =  A(2,1)*A(3,3)*A(4,2) + A(3,1)*A(2,2)*A(4,3) + A(4,1)*A(2,3)*A(3,2) 		&
				- A(2,1)*A(3,2)*A(4,3) - A(3,1)*A(2,3)*A(4,2) - A(4,1)*A(2,2)*A(3,3)

		B(4,2) =  A(1,1)*A(3,2)*A(4,3) + A(3,1)*A(1,3)*A(4,2) + A(4,1)*A(1,2)*A(3,3)		&
				- A(1,1)*A(3,3)*A(4,2) - A(3,1)*A(1,2)*A(4,3) - A(4,1)*A(1,3)*A(3,2) 
				
		B(4,3) =  A(1,1)*A(2,3)*A(4,2) + A(2,1)*A(1,2)*A(4,3) + A(4,1)*A(1,3)*A(2,2) 		&
				- A(1,1)*A(2,2)*A(4,3) - A(2,1)*A(1,3)*A(4,2) - A(4,1)*A(1,2)*A(2,3)
				
		B(4,4) =  A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(1,3)*A(3,2) + A(3,1)*A(1,2)*A(2,3)		&
		 		- A(1,1)*A(2,3)*A(3,2) - A(2,1)*A(1,2)*A(3,3) - A(3,1)*A(1,3)*A(2,2)		
		
		B = B/DetA
	!-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE Invert4x4_Small
    
	SUBROUTINE Invert4x4Big_Cell_Zone(A,B)
	!-------------------------------------------------------------------------------------------------------------------!
	! Calculated the inverse of a 4x4 matrix exactly. Due to the small size of this matrix, the inverse can be done
	! by a finite number of algebra steps. An iterative solution is not necessary.
	! The steps used to calculate this inverse are from a Maple output
	! This 4x4 matrix actually is a xmax by ymax array with a 4x4 matrix at each cell. The dimensions are (5,5,xmax,ymax)
	! The extra row/column is an artifact of the zonal matrix, which needs 5 rows/columns, but the dependency is removed
	! by this point in the program. Using the same variable and ignoring those values seemed more advantageous than making a 
	! new variable of size (4,4,xmax,ymax)
        TYPE(Cell_Zone), 		 INTENT (IN) :: A
        TYPE(Cell_Zone), 		 INTENT(OUT) :: B
		REAL(PRCSN), 	DIMENSION(xmax,ymax) :: DetA 	! Determinant of A
		REAL(PRCSN), 	DIMENSION(xmax,ymax) :: MarkZeros
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		DetA = 0D0
		B%cz = 0D0
		MarkZeros = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate the determinant of A
		DetA =   ( A%cz(1,1,:,:)*(A%cz(2,2,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:) + A%cz(3,2,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) 	&
								+ A%cz(4,2,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:) - A%cz(2,2,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) 	&
								- A%cz(3,2,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:) - A%cz(4,2,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:))	&

				 + A%cz(2,1,:,:)*(A%cz(1,2,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) + A%cz(3,2,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:) 	&
				 				+ A%cz(4,2,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:) - A%cz(1,2,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:) 	&
				 				- A%cz(3,2,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) - A%cz(4,2,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:))	&

				 + A%cz(3,1,:,:)*(A%cz(1,2,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:) + A%cz(2,2,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) 	&
				 				+ A%cz(4,2,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:) - A%cz(1,2,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) 	&
				 				- A%cz(2,2,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:) - A%cz(4,2,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:))	&

				 + A%cz(4,1,:,:)*(A%cz(1,2,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:) + A%cz(2,2,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:) 	&
				 				+ A%cz(3,2,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:) - A%cz(1,2,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:)		&
				 				- A%cz(2,2,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:) - A%cz(3,2,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:)))

		!---------------------------------------------------------------------------------------------------------------!
		! Calculate elements of matrix B
		B%cz(1,1,:,:) =   A%cz(2,2,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:) + A%cz(3,2,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) 	&
						+ A%cz(4,2,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:) - A%cz(4,2,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:) 	&
						- A%cz(2,2,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) - A%cz(3,2,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:)

		B%cz(1,2,:,:) =   A%cz(1,2,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) + A%cz(3,2,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:) 	&
						+ A%cz(4,2,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:)	- A%cz(1,2,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:) 	&
						- A%cz(3,2,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) - A%cz(4,2,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:)

		B%cz(1,3,:,:) =   A%cz(1,2,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:) + A%cz(2,2,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) 	&
						+ A%cz(4,2,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:) - A%cz(4,2,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:) 	&
						- A%cz(1,2,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) - A%cz(2,2,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:)

		B%cz(1,4,:,:) =   A%cz(1,2,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:) + A%cz(2,2,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:) 	&
						+ A%cz(3,2,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:) - A%cz(1,2,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:)	 	&
						- A%cz(2,2,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:) - A%cz(3,2,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:)

		B%cz(2,1,:,:) =   A%cz(2,1,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) + A%cz(3,1,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:) - A%cz(2,1,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:)		&
						- A%cz(3,1,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) - A%cz(4,1,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:)

		B%cz(2,2,:,:) =   A%cz(1,1,:,:)*A%cz(3,3,:,:)*A%cz(4,4,:,:) + A%cz(3,1,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:) - A%cz(4,1,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:)		&
						- A%cz(1,1,:,:)*A%cz(3,4,:,:)*A%cz(4,3,:,:) - A%cz(3,1,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:)

		B%cz(2,3,:,:) =   A%cz(1,1,:,:)*A%cz(2,4,:,:)*A%cz(4,3,:,:) + A%cz(2,1,:,:)*A%cz(1,3,:,:)*A%cz(4,4,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:) - A%cz(1,1,:,:)*A%cz(2,3,:,:)*A%cz(4,4,:,:) 	&
						- A%cz(2,1,:,:)*A%cz(1,4,:,:)*A%cz(4,3,:,:) - A%cz(4,1,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:) 
		
		B%cz(2,4,:,:) =   A%cz(1,1,:,:)*A%cz(2,3,:,:)*A%cz(3,4,:,:) + A%cz(2,1,:,:)*A%cz(1,4,:,:)*A%cz(3,3,:,:) 	&
						+ A%cz(3,1,:,:)*A%cz(1,3,:,:)*A%cz(2,4,:,:) - A%cz(3,1,:,:)*A%cz(1,4,:,:)*A%cz(2,3,:,:)		&
						- A%cz(1,1,:,:)*A%cz(2,4,:,:)*A%cz(3,3,:,:) - A%cz(2,1,:,:)*A%cz(1,3,:,:)*A%cz(3,4,:,:)

		B%cz(3,1,:,:) =   A%cz(2,1,:,:)*A%cz(3,2,:,:)*A%cz(4,4,:,:) + A%cz(3,1,:,:)*A%cz(2,4,:,:)*A%cz(4,2,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(2,2,:,:)*A%cz(3,4,:,:) - A%cz(2,1,:,:)*A%cz(3,4,:,:)*A%cz(4,2,:,:) 	&
						- A%cz(3,1,:,:)*A%cz(2,2,:,:)*A%cz(4,4,:,:) - A%cz(4,1,:,:)*A%cz(2,4,:,:)*A%cz(3,2,:,:) 
				
		B%cz(3,2,:,:) =   A%cz(1,1,:,:)*A%cz(3,4,:,:)*A%cz(4,2,:,:) + A%cz(3,1,:,:)*A%cz(1,2,:,:)*A%cz(4,4,:,:)		& 	
						+ A%cz(4,1,:,:)*A%cz(1,4,:,:)*A%cz(3,2,:,:) - A%cz(1,1,:,:)*A%cz(3,2,:,:)*A%cz(4,4,:,:) 	&
						- A%cz(3,1,:,:)*A%cz(1,4,:,:)*A%cz(4,2,:,:) - A%cz(4,1,:,:)*A%cz(1,2,:,:)*A%cz(3,4,:,:)
		
		B%cz(3,3,:,:) =   A%cz(1,1,:,:)*A%cz(2,2,:,:)*A%cz(4,4,:,:) + A%cz(2,1,:,:)*A%cz(1,4,:,:)*A%cz(4,2,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(1,2,:,:)*A%cz(2,4,:,:) - A%cz(1,1,:,:)*A%cz(2,4,:,:)*A%cz(4,2,:,:) 	&
						- A%cz(2,1,:,:)*A%cz(1,2,:,:)*A%cz(4,4,:,:) - A%cz(4,1,:,:)*A%cz(1,4,:,:)*A%cz(2,2,:,:) 
				
		B%cz(3,4,:,:) =   A%cz(1,1,:,:)*A%cz(2,4,:,:)*A%cz(3,2,:,:) + A%cz(2,1,:,:)*A%cz(1,2,:,:)*A%cz(3,4,:,:) 	&
						+ A%cz(3,1,:,:)*A%cz(1,4,:,:)*A%cz(2,2,:,:) - A%cz(1,1,:,:)*A%cz(2,2,:,:)*A%cz(3,4,:,:) 	&
						- A%cz(2,1,:,:)*A%cz(1,4,:,:)*A%cz(3,2,:,:) - A%cz(3,1,:,:)*A%cz(1,2,:,:)*A%cz(2,4,:,:)
				
		B%cz(4,1,:,:) =   A%cz(2,1,:,:)*A%cz(3,3,:,:)*A%cz(4,2,:,:) + A%cz(3,1,:,:)*A%cz(2,2,:,:)*A%cz(4,3,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(2,3,:,:)*A%cz(3,2,:,:) - A%cz(2,1,:,:)*A%cz(3,2,:,:)*A%cz(4,3,:,:) 	&
						- A%cz(3,1,:,:)*A%cz(2,3,:,:)*A%cz(4,2,:,:) - A%cz(4,1,:,:)*A%cz(2,2,:,:)*A%cz(3,3,:,:)

		B%cz(4,2,:,:) =   A%cz(1,1,:,:)*A%cz(3,2,:,:)*A%cz(4,3,:,:) + A%cz(3,1,:,:)*A%cz(1,3,:,:)*A%cz(4,2,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(1,2,:,:)*A%cz(3,3,:,:) - A%cz(1,1,:,:)*A%cz(3,3,:,:)*A%cz(4,2,:,:) 	&
						- A%cz(3,1,:,:)*A%cz(1,2,:,:)*A%cz(4,3,:,:) - A%cz(4,1,:,:)*A%cz(1,3,:,:)*A%cz(3,2,:,:) 
				
		B%cz(4,3,:,:) =   A%cz(1,1,:,:)*A%cz(2,3,:,:)*A%cz(4,2,:,:) + A%cz(2,1,:,:)*A%cz(1,2,:,:)*A%cz(4,3,:,:) 	&
						+ A%cz(4,1,:,:)*A%cz(1,3,:,:)*A%cz(2,2,:,:) - A%cz(1,1,:,:)*A%cz(2,2,:,:)*A%cz(4,3,:,:) 	&
						- A%cz(2,1,:,:)*A%cz(1,3,:,:)*A%cz(4,2,:,:) - A%cz(4,1,:,:)*A%cz(1,2,:,:)*A%cz(2,3,:,:)
				
		B%cz(4,4,:,:) =   A%cz(1,1,:,:)*A%cz(2,2,:,:)*A%cz(3,3,:,:) + A%cz(2,1,:,:)*A%cz(1,3,:,:)*A%cz(3,2,:,:) 	&
						+ A%cz(3,1,:,:)*A%cz(1,2,:,:)*A%cz(2,3,:,:) - A%cz(1,1,:,:)*A%cz(2,3,:,:)*A%cz(3,2,:,:) 	&
						- A%cz(2,1,:,:)*A%cz(1,2,:,:)*A%cz(3,3,:,:) - A%cz(3,1,:,:)*A%cz(1,3,:,:)*A%cz(2,2,:,:)		

		WHERE (ABS(DetA) > 1.D-10)
			! Must do each element separately (4x4 = 16 total) so that mask of DetA(i,j) fits the mask of B%cz(i,j,i0,j0)
			B%cz(1,1,:,:) = B%cz(1,1,:,:)/DetA
			B%cz(1,2,:,:) = B%cz(1,2,:,:)/DetA
			B%cz(1,3,:,:) = B%cz(1,3,:,:)/DetA
			B%cz(1,4,:,:) = B%cz(1,4,:,:)/DetA

			B%cz(2,1,:,:) = B%cz(2,1,:,:)/DetA
			B%cz(2,2,:,:) = B%cz(2,2,:,:)/DetA
			B%cz(2,3,:,:) = B%cz(2,3,:,:)/DetA
			B%cz(2,4,:,:) = B%cz(2,4,:,:)/DetA

			B%cz(3,1,:,:) = B%cz(3,1,:,:)/DetA
			B%cz(3,2,:,:) = B%cz(3,2,:,:)/DetA
			B%cz(3,3,:,:) = B%cz(3,3,:,:)/DetA
			B%cz(3,4,:,:) = B%cz(3,4,:,:)/DetA

			B%cz(4,1,:,:) = B%cz(4,1,:,:)/DetA
			B%cz(4,2,:,:) = B%cz(4,2,:,:)/DetA
			B%cz(4,3,:,:) = B%cz(4,3,:,:)/DetA
			B%cz(4,4,:,:) = B%cz(4,4,:,:)/DetA

			MarkZeros = 0D0
		ELSEWHERE
			MarkZeros = 1D0
		END WHERE
		
		IF (SUM(MarkZeros) > 0D0) THEN
			PRINT '(/,TR4,A10,I5,A31)', 'There are ',INT(SUM(MarkZeros)),' cells with a zero determinant:'
			IF (SUM(MarkZeros) < 30D0) THEN
				PRINT '(/,TR4,A64)','The following row, column pairs are cells with zero determinants'
				DO j=1,ymax
					DO i=i,xmax
						IF (MarKZeros(i,j) > 0D0) PRINT '(TR4,A1,I3,I3,A1)','(',i,j,')'
					ENDDO
				ENDDO
			ENDIF
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE Invert4x4Big_Cell_Zone
    
	SUBROUTINE Invert4x4Big_Exact_Size(B)
	!-------------------------------------------------------------------------------------------------------------------!
	! Calculated the inverse of a 4x4 matrix exactly. Due to the small size of this matrix, the inverse can be done
	! by a finite number of algebra steps. An iterative solution is not necessary.
	! The steps used to calculate this inverse are from a Maple output
	! This 4x4 matrix actually is a xmax by ymax array with a 4x4 matrix at each cell. The dimensions are (4,4,xmax,ymax)
	! That makes this the 'exact' size, because the zonal data structure is 'big', with a dimension of (5,5,xmax,ymax)
	! The 'small' version of this function is simply (4,4)
		REAL(PRCSN), DIMENSION(4,4,xmax,ymax), INTENT(INOUT) :: B 
		REAL(PRCSN), DIMENSION(4,4,xmax,ymax)				 :: A 
		REAL(PRCSN), DIMENSION(xmax,ymax) 				  	 :: DetA 	! Determinant of A
		REAL(PRCSN), DIMENSION(xmax,ymax) 				  	 :: MarkZeros
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		A 		  = 0D0
		DetA 	  = 0D0
		MarkZeros = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		A = B

		! Calculate the determinant of A
		DetA =   ( A(1,1,:,:)*(A(2,2,:,:)*A(3,3,:,:)*A(4,4,:,:) + A(3,2,:,:)*A(2,4,:,:)*A(4,3,:,:) 	&
							 + A(4,2,:,:)*A(2,3,:,:)*A(3,4,:,:) - A(2,2,:,:)*A(3,4,:,:)*A(4,3,:,:) 	&
							 - A(3,2,:,:)*A(2,3,:,:)*A(4,4,:,:) - A(4,2,:,:)*A(2,4,:,:)*A(3,3,:,:))	&

				 + A(2,1,:,:)*(A(1,2,:,:)*A(3,4,:,:)*A(4,3,:,:) + A(3,2,:,:)*A(1,3,:,:)*A(4,4,:,:) 	&
				 			 + A(4,2,:,:)*A(1,4,:,:)*A(3,3,:,:) - A(1,2,:,:)*A(3,3,:,:)*A(4,4,:,:) 	&
				 			 - A(3,2,:,:)*A(1,4,:,:)*A(4,3,:,:) - A(4,2,:,:)*A(1,3,:,:)*A(3,4,:,:))	&

				 + A(3,1,:,:)*(A(1,2,:,:)*A(2,3,:,:)*A(4,4,:,:) + A(2,2,:,:)*A(1,4,:,:)*A(4,3,:,:) 	&
							 + A(4,2,:,:)*A(1,3,:,:)*A(2,4,:,:) - A(1,2,:,:)*A(2,4,:,:)*A(4,3,:,:) 	&
			 				 - A(2,2,:,:)*A(1,3,:,:)*A(4,4,:,:) - A(4,2,:,:)*A(1,4,:,:)*A(2,3,:,:))	&

				 + A(4,1,:,:)*(A(1,2,:,:)*A(2,4,:,:)*A(3,3,:,:) + A(2,2,:,:)*A(1,3,:,:)*A(3,4,:,:) 	&
				 			 + A(3,2,:,:)*A(1,4,:,:)*A(2,3,:,:) - A(1,2,:,:)*A(2,3,:,:)*A(3,4,:,:)	&
				 			 - A(2,2,:,:)*A(1,4,:,:)*A(3,3,:,:) - A(3,2,:,:)*A(1,3,:,:)*A(2,4,:,:)))
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate elements of matrix B
		B(1,1,:,:) =   A(2,2,:,:)*A(3,3,:,:)*A(4,4,:,:) + A(3,2,:,:)*A(2,4,:,:)*A(4,3,:,:) 	&
					 + A(4,2,:,:)*A(2,3,:,:)*A(3,4,:,:) - A(4,2,:,:)*A(2,4,:,:)*A(3,3,:,:) 	&
					 - A(2,2,:,:)*A(3,4,:,:)*A(4,3,:,:) - A(3,2,:,:)*A(2,3,:,:)*A(4,4,:,:)

		B(1,2,:,:) =   A(1,2,:,:)*A(3,4,:,:)*A(4,3,:,:) + A(3,2,:,:)*A(1,3,:,:)*A(4,4,:,:) 	&
					 + A(4,2,:,:)*A(1,4,:,:)*A(3,3,:,:)	- A(1,2,:,:)*A(3,3,:,:)*A(4,4,:,:) 	&
					 - A(3,2,:,:)*A(1,4,:,:)*A(4,3,:,:) - A(4,2,:,:)*A(1,3,:,:)*A(3,4,:,:)

		B(1,3,:,:) =   A(1,2,:,:)*A(2,3,:,:)*A(4,4,:,:) + A(2,2,:,:)*A(1,4,:,:)*A(4,3,:,:) 	&
					 + A(4,2,:,:)*A(1,3,:,:)*A(2,4,:,:) - A(4,2,:,:)*A(1,4,:,:)*A(2,3,:,:) 	&
					 - A(1,2,:,:)*A(2,4,:,:)*A(4,3,:,:) - A(2,2,:,:)*A(1,3,:,:)*A(4,4,:,:)

		B(1,4,:,:) =   A(1,2,:,:)*A(2,4,:,:)*A(3,3,:,:) + A(2,2,:,:)*A(1,3,:,:)*A(3,4,:,:) 	&
					 + A(3,2,:,:)*A(1,4,:,:)*A(2,3,:,:) - A(1,2,:,:)*A(2,3,:,:)*A(3,4,:,:)	&
					 - A(2,2,:,:)*A(1,4,:,:)*A(3,3,:,:) - A(3,2,:,:)*A(1,3,:,:)*A(2,4,:,:)

		B(2,1,:,:) =   A(2,1,:,:)*A(3,4,:,:)*A(4,3,:,:) + A(3,1,:,:)*A(2,3,:,:)*A(4,4,:,:) 	&
					 + A(4,1,:,:)*A(2,4,:,:)*A(3,3,:,:) - A(2,1,:,:)*A(3,3,:,:)*A(4,4,:,:)	&
					 - A(3,1,:,:)*A(2,4,:,:)*A(4,3,:,:) - A(4,1,:,:)*A(2,3,:,:)*A(3,4,:,:)

		B(2,2,:,:) =   A(1,1,:,:)*A(3,3,:,:)*A(4,4,:,:) + A(3,1,:,:)*A(1,4,:,:)*A(4,3,:,:) 	&
					 + A(4,1,:,:)*A(1,3,:,:)*A(3,4,:,:) - A(4,1,:,:)*A(1,4,:,:)*A(3,3,:,:)	&
					 - A(1,1,:,:)*A(3,4,:,:)*A(4,3,:,:) - A(3,1,:,:)*A(1,3,:,:)*A(4,4,:,:)

		B(2,3,:,:) =   A(1,1,:,:)*A(2,4,:,:)*A(4,3,:,:) + A(2,1,:,:)*A(1,3,:,:)*A(4,4,:,:) 	&
					 + A(4,1,:,:)*A(1,4,:,:)*A(2,3,:,:) - A(1,1,:,:)*A(2,3,:,:)*A(4,4,:,:) 	&
					 - A(2,1,:,:)*A(1,4,:,:)*A(4,3,:,:) - A(4,1,:,:)*A(1,3,:,:)*A(2,4,:,:) 
		
		B(2,4,:,:) =   A(1,1,:,:)*A(2,3,:,:)*A(3,4,:,:) + A(2,1,:,:)*A(1,4,:,:)*A(3,3,:,:) 	&
					 + A(3,1,:,:)*A(1,3,:,:)*A(2,4,:,:) - A(3,1,:,:)*A(1,4,:,:)*A(2,3,:,:)	&
					 - A(1,1,:,:)*A(2,4,:,:)*A(3,3,:,:) - A(2,1,:,:)*A(1,3,:,:)*A(3,4,:,:)

		B(3,1,:,:) =   A(2,1,:,:)*A(3,2,:,:)*A(4,4,:,:) + A(3,1,:,:)*A(2,4,:,:)*A(4,2,:,:) 	&
					 + A(4,1,:,:)*A(2,2,:,:)*A(3,4,:,:) - A(2,1,:,:)*A(3,4,:,:)*A(4,2,:,:) 	&
					 - A(3,1,:,:)*A(2,2,:,:)*A(4,4,:,:) - A(4,1,:,:)*A(2,4,:,:)*A(3,2,:,:) 
				
		B(3,2,:,:) =   A(1,1,:,:)*A(3,4,:,:)*A(4,2,:,:) + A(3,1,:,:)*A(1,2,:,:)*A(4,4,:,:)	& 	
					 + A(4,1,:,:)*A(1,4,:,:)*A(3,2,:,:) - A(1,1,:,:)*A(3,2,:,:)*A(4,4,:,:) 	&
					 - A(3,1,:,:)*A(1,4,:,:)*A(4,2,:,:) - A(4,1,:,:)*A(1,2,:,:)*A(3,4,:,:)
		
		B(3,3,:,:) =   A(1,1,:,:)*A(2,2,:,:)*A(4,4,:,:) + A(2,1,:,:)*A(1,4,:,:)*A(4,2,:,:) 	&
					 + A(4,1,:,:)*A(1,2,:,:)*A(2,4,:,:) - A(1,1,:,:)*A(2,4,:,:)*A(4,2,:,:) 	&
					 - A(2,1,:,:)*A(1,2,:,:)*A(4,4,:,:) - A(4,1,:,:)*A(1,4,:,:)*A(2,2,:,:) 
				
		B(3,4,:,:) =   A(1,1,:,:)*A(2,4,:,:)*A(3,2,:,:) + A(2,1,:,:)*A(1,2,:,:)*A(3,4,:,:) 	&
					 + A(3,1,:,:)*A(1,4,:,:)*A(2,2,:,:) - A(1,1,:,:)*A(2,2,:,:)*A(3,4,:,:) 	&
					 - A(2,1,:,:)*A(1,4,:,:)*A(3,2,:,:) - A(3,1,:,:)*A(1,2,:,:)*A(2,4,:,:)
				
		B(4,1,:,:) =   A(2,1,:,:)*A(3,3,:,:)*A(4,2,:,:) + A(3,1,:,:)*A(2,2,:,:)*A(4,3,:,:) 	&
					 + A(4,1,:,:)*A(2,3,:,:)*A(3,2,:,:) - A(2,1,:,:)*A(3,2,:,:)*A(4,3,:,:) 	&
					 - A(3,1,:,:)*A(2,3,:,:)*A(4,2,:,:) - A(4,1,:,:)*A(2,2,:,:)*A(3,3,:,:)

		B(4,2,:,:) =   A(1,1,:,:)*A(3,2,:,:)*A(4,3,:,:) + A(3,1,:,:)*A(1,3,:,:)*A(4,2,:,:) 	&
					 + A(4,1,:,:)*A(1,2,:,:)*A(3,3,:,:) - A(1,1,:,:)*A(3,3,:,:)*A(4,2,:,:) 	&
					 - A(3,1,:,:)*A(1,2,:,:)*A(4,3,:,:) - A(4,1,:,:)*A(1,3,:,:)*A(3,2,:,:) 
				
		B(4,3,:,:) =   A(1,1,:,:)*A(2,3,:,:)*A(4,2,:,:) + A(2,1,:,:)*A(1,2,:,:)*A(4,3,:,:) 	&
					 + A(4,1,:,:)*A(1,3,:,:)*A(2,2,:,:) - A(1,1,:,:)*A(2,2,:,:)*A(4,3,:,:) 	&
					 - A(2,1,:,:)*A(1,3,:,:)*A(4,2,:,:) - A(4,1,:,:)*A(1,2,:,:)*A(2,3,:,:)
				
		B(4,4,:,:) =   A(1,1,:,:)*A(2,2,:,:)*A(3,3,:,:) + A(2,1,:,:)*A(1,3,:,:)*A(3,2,:,:) 	&
					 + A(3,1,:,:)*A(1,2,:,:)*A(2,3,:,:) - A(1,1,:,:)*A(2,3,:,:)*A(3,2,:,:) 	&
					 - A(2,1,:,:)*A(1,2,:,:)*A(3,3,:,:) - A(3,1,:,:)*A(1,3,:,:)*A(2,2,:,:)		

		WHERE (ABS(DetA) > 1.D-12)
			! Must do each element separately (4x4 = 16 total) so that mask of DetA(i,j) fits the mask of B(i,j,i0,j0)
			B(1,1,:,:) = B(1,1,:,:)/DetA
			B(1,2,:,:) = B(1,2,:,:)/DetA
			B(1,3,:,:) = B(1,3,:,:)/DetA
			B(1,4,:,:) = B(1,4,:,:)/DetA

			B(2,1,:,:) = B(2,1,:,:)/DetA
			B(2,2,:,:) = B(2,2,:,:)/DetA
			B(2,3,:,:) = B(2,3,:,:)/DetA
			B(2,4,:,:) = B(2,4,:,:)/DetA

			B(3,1,:,:) = B(3,1,:,:)/DetA
			B(3,2,:,:) = B(3,2,:,:)/DetA
			B(3,3,:,:) = B(3,3,:,:)/DetA
			B(3,4,:,:) = B(3,4,:,:)/DetA

			B(4,1,:,:) = B(4,1,:,:)/DetA
			B(4,2,:,:) = B(4,2,:,:)/DetA
			B(4,3,:,:) = B(4,3,:,:)/DetA
			B(4,4,:,:) = B(4,4,:,:)/DetA

			MarkZeros = 0D0
		ELSEWHERE
			MarkZeros = 1D0
		END WHERE

		!WRITE(*,'(4(4(F8.6,2X),/))') ((B(:,:,i,j),i=1,xmax),j=1,ymax)		

		IF (SUM(MarkZeros) > 0D0) THEN
			PRINT '(/,TR4,A10,I7,A30)', 'There are ',INT(SUM(MarkZeros)),'cells with a zero determinant:'
			IF (SUM(MarkZeros) < 30D0) THEN
				PRINT '(/,TR4,A64)','The following row, column pairs are cells with zero determinants'
				DO j=1,ymax
					DO i=i,xmax
						IF (MarKZeros(i,j) > 0D0) PRINT '(TR4,A1,I3,I3,A1)','(',i,j,')'
					ENDDO
				ENDDO
			ENDIF
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE Invert4x4Big_Exact_Size

	FUNCTION Sparse_Mat_Vec_Multiply(SpaZ,in_vec)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine computes the product of a matrix and a vector without specifically forming the matrix
	! The input matrix is in a sparse form [3 columns: row, column, value], 1 row for each entry
	! The output is a vector. This is needed for each step of the Krylov iteration in the ConjGrad solver
		REAL(PRCSN), DIMENSION(VecLength)	:: Sparse_Mat_Vec_Multiply		
		TYPE(Sparse_Matrix)					:: SpaZ		! input sparse matrix
		REAL(PRCSN), DIMENSION(VecLength) 	:: in_vec	! input vector to multiply by matrix
		INTEGER								:: k		! iteration integer
		!---------------------------------------------------------------------------------------------------------------!

		Sparse_Mat_Vec_Multiply = 0D0
		
		DO k=1,VecLength*8
			i = INT(SpaZ%sgz(k,1))
			j = INT(SpaZ%sgz(k,2))
			IF ((i/=0) .AND. (j/=0)) Sparse_Mat_Vec_Multiply(i) = Sparse_Mat_Vec_Multiply(i) + in_vec(j)*SpaZ%sgz(k,3)
		ENDDO

		!print *, 'magnitude of mat-vec multiplication:  ',DOT_PRODUCT(out_vec,out_vec)
	!-------------------------------------------------------------------------------------------------------------------!
	END FUNCTION Sparse_Mat_Vec_Multiply

!	SUBROUTINE Sparse_Mat_Vec_Multiply(SpaZ,in_vec,out_vec)
!	!-------------------------------------------------------------------------------------------------------------------!
!	! This subroutine computes the product of a matrix and a vector without specifically forming the matrix
!	! The input matrix is in a sparse form [3 columns: row, column, value], 1 row for each entry
!	! The output is a vector. This is needed for each step of the Krylov iteration in the ConjGrad solver
!		TYPE(Sparse_Matrix),					INTENT( IN) :: SpaZ		! input sparse matrix
!		REAL(PRCSN), DIMENSION(VecLength),		INTENT( IN) :: in_vec	! input vector to multiply by matrix
!		REAL(PRCSN), DIMENSION(VecLength),		INTENT(OUT) :: out_vec	! output result of matrix-vector product
!		INTEGER												:: k		! iteration integer
!		!---------------------------------------------------------------------------------------------------------------!
!
!		out_vec = 0D0
!		
!		DO k=1,VecLength*8
!			i = INT(SpaZ%sgz(k,1))
!			j = INT(SpaZ%sgz(k,2))
!			IF ((i/=0) .AND. (j/=0)) out_vec(i) = out_vec(i) + in_vec(j)*SpaZ%sgz(k,3)
!		ENDDO
!
!		!print *, 'magnitude of mat-vec multiplication:  ',DOT_PRODUCT(out_vec,out_vec)
!	!-------------------------------------------------------------------------------------------------------------------!
!	END SUBROUTINE Sparse_Mat_Vec_Multiply

    SUBROUTINE Conjugate_Gradient_Full(A_mat,x_vec,b_vec,guess)
	!-------------------------------------------------------------------------------------------------------------------!
	! Solves matrix equation, A*x=b, where A is a SPD matrix (Symmetric Positive-Definite)
	! 	A_mat is the matrix of size (VecLength x VecLength)
	!	There are several vectors involved, all of the same length
	!	b_vec is the RHS of the matrix equation, A*x=b
	!	x_vec is the solution to the matrix equation, which this subroutine calculates
	!	guess is an optional vector that can be used as a initial guess for x_vec
	! 	Ap simply stores the matrix-vector product of A_mat and p_vec (reduces computational work)
	! 	p_vec and r_vec are vectors needed for the CG algorithm
		TYPE(Global_Zone),							INTENT( IN) :: A_mat
		REAL(PRCSN), DIMENSION(VecLength), 			INTENT( IN) :: b_vec
		REAL(PRCSN), DIMENSION(VecLength), OPTIONAL,INTENT( IN) :: guess
		REAL(PRCSN), DIMENSION(VecLength), 			INTENT(OUT) :: x_vec
		REAL(PRCSN), DIMENSION(VecLength)						:: Ap, p_vec, r_vec
		REAL(PRCSN)												:: alpha, RmagO, RmagN
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		Ap = 0D0;	  p_vec = 0D0;	r_vec = 0D0
		alpha = 0D0; RmagO = 0D0; RmagN = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize - Define inital values of vectors x, r, and p (and magnitude of r)
		IF (PRESENT(guess)) THEN
			x_vec = guess				! Use guess for x_vec if it is supplied
			r_vec = b_vec - MATMUL(A_mat%gz,x_vec)
		ELSE
			x_vec = 0D0				! Use 0 for initial guess otherwise
			r_vec = b_vec				! Avoid one multiplication by the matrix
		ENDIF

		p_vec = r_vec
		RmagN = DOT_PRODUCT(r_vec,r_vec)

		!--------------------------------------------------------------------------------------------------------------!		
		! Conjugate Gradient Algorithm
		DO i=1, FLOOR(1.5*VecLength)				! Infinite precision would always converge at i=VecLength
			IF ((i==1) .AND. (RmagN < 1.D-10)) EXIT	! Check in case guess was really good
			RmagO = RmagN							! Set 'old' value to previous iterations 'new' value

			Ap = MATMUL(A_mat%gz,p_vec)				! Calculate this iteration's Matrix-Vector product
			
			alpha = RmagO / DOT_PRODUCT(p_vec,Ap)	! Update alpha
			
			x_vec = x_vec + alpha*p_vec				! Update x
			r_vec = r_vec - alpha*Ap				! Update r
			
			RmagN = DOT_PRODUCT(r_vec,r_vec)		! Calculate new magnitude of r
			
			p_vec = r_vec + RmagN/RmagO*p_vec		! Update p

			IF (RmagO < 1.D-10) EXIT				! Stopping criteria (NOTE: CG stopping criteria has many options)
		ENDDO		
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Conjugate_Gradient_Full
	
	SUBROUTINE Conjugate_Gradient_Sparse(A_mat,x_vec,b_vec,guess)
	!-------------------------------------------------------------------------------------------------------------------!
	! Solves matrix equation, A*x=b, where A is a SPD matrix (Symmetric Positive-Definite)
	! 	A_mat is the matrix of size (VecLength x VecLength)
	!	There are several vectors involved, all of the same length
	!	b_vec is the RHS of the matrix equation, A*x=b
	!	x_vec is the solution to the matrix equation, which this subroutine calculates
	!	guess is an optional vector that can be used as a initial guess for x_vec
	! 	Ap simply stores the matrix-vector product of A_mat and p_vec (reduces computational work)
	! 	p_vec and r_vec are vectors needed for the CG algorithm
		TYPE(Sparse_Matrix),						INTENT( IN) :: A_mat				! input matrix
		REAL(PRCSN), DIMENSION(VecLength), 			INTENT( IN) :: b_vec				! input RHS
		REAL(PRCSN), DIMENSION(VecLength), OPTIONAL,INTENT( IN) :: guess				! optional guess
		REAL(PRCSN), DIMENSION(VecLength), 			INTENT(OUT) :: x_vec				! output solution
		REAL(PRCSN), DIMENSION(VecLength)						:: Ap, p_vec, r_vec		! temporary vectors
		REAL(PRCSN)												:: alpha, RmagO, RmagN	! temporary scalars
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		Ap = 0D0;	  p_vec = 0D0;	r_vec = 0D0
		alpha = 0D0; RmagO = 0D0; RmagN = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize - Define inital values of vectors x, r, and p (and magnitude of r)
		IF (PRESENT(guess) .EQV. .FALSE.) THEN
			x_vec = 0D0				! Use 0 for initial guess otherwise
			r_vec = b_vec				! Avoid one multiplication by the matrix
		ELSE
			x_vec = guess				! Use guess for x_vec if it is supplied
			!CALl Sparse_Mat_Vec_Multiply(A_mat,guess,Ap)
			AP = Sparse_Mat_Vec_Multiply(A_mat,guess)
			r_vec = b_vec - Ap			
		ENDIF

		p_vec = r_vec
		RmagN = DOT_PRODUCT(r_vec,r_vec)

		!--------------------------------------------------------------------------------------------------------------!		
		! Conjugate Gradient Algorithm
		DO i0=1, FLOOR(1.5*VecLength)				! Infinite precision would always converge at i=VecLength

			IF ((i0>1) .AND. (RmagN < CGtol*1.D-5)) THEN ! Check in case guess was really good
				!PRINT *, 'Good Guess!'
				EXIT	
			ENDIF
			RmagO = RmagN							! Set 'old' value to previous iterations 'new' value

			!CALl Sparse_Mat_Vec_Multiply(A_mat,p_vec,Ap)
			AP = Sparse_Mat_Vec_Multiply(A_mat,p_vec)
			
			alpha = RmagO / DOT_PRODUCT(p_vec,Ap)	! Update alpha
			
			x_vec = x_vec + alpha*p_vec				! Update x
			r_vec = r_vec - alpha*Ap				! Update r
			
			RmagN = DOT_PRODUCT(r_vec,r_vec)		! Calculate new magnitude of r
			
			p_vec = r_vec + RmagN/RmagO*p_vec		! Update p

			!IF (RmagO < 1.D-10) EXIT				! Stopping criteria (NOTE: CG stopping criteria has many options)
			IF ((RmagN < CGtol) .AND. (RmagN/RmagO > 1.D-3)) EXIT
		ENDDO		
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Conjugate_Gradient_Sparse
	
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
!-----------------------------------------------------------------------------------------------------------------------!
END MODULE SOM_Utilities

MODULE SOM_Subroutine_Module
!-----------------------------------------------------------------------------------------------------------------------!
! This module contains all the subroutines called in the program SOM.
! This module also calls module Types to define several derived data types.
! CONVENTION: FORTRAN COMMANDS ARE ALL CAPS; variables I define may be lower case or MiXeD
! CONVENTION: Separate separate sections with line breaks, 
! CONVENTION: and try to bind sections together with !------! to fill an intire line
    USE DefineParameters
	USE SOM_Utilities
    USE Types
	IMPLICIT NONE
	
	!-------------------------------------------------------------------------------------------------------------------!
	! Interfaces must be specified for an overloaded function
	! To convert from cell-centered to face-centered, different procedures must be used whether the variable is a scalar
	! or a tensor. For simplicity, I will overload the function Cell_To_Face to handle either case, sicne I need both.
	INTERFACE Cell_To_Face
		MODULE PROCEDURE Cell_To_Face_Scalar, Cell_To_Face_Tensor
	END INTERFACE
! Turns out I never call Cell_To_Face_Scalar
	
	! To convert from Face-Centered to cell centered, one can simply take the average of the four faces,
	! or one can calculate phi_c from the four face values, as well as some info from the zonal matrix
	INTERFACE Face_To_Cell
		MODULE PROCEDURE Face_To_Cell_Average, Face_To_Cell_Zonal
	END INTERFACE
! Turns out I never call Face_To_Cell_Average

	INTERFACE Assemble_Global_Matrix
	! Assemble_Global_Matrix subroutine will output either a sparse or dense matrix
		MODULE PROCEDURE Assemble_Global_Matrix_Full, Assemble_Global_Matrix_Sparse
	END INTERFACE
	!-------------------------------------------------------------------------------------------------------------------!
	
	CONTAINS
	!-------------------------------------------------------------------------------------------------------------------!

	FUNCTION Phi_Value(x,y)
	!-------------------------------------------------------------------------------------------------------------------!
	! Specifies phi on a cell-centered mesh using an analytic solution
	! The only input is the current time step (n), which determines the time
	! The output is a file saved for each time this subroutine is called
		REAL(PRCSN) 				:: Phi_Value, x, y
		REAL(PRCSN)				 	:: time
		REAL(PRCSN) 				:: r, a, beta!, jump	! values that are usually not defined but cause compilation error
		TYPE(Tensor_C)				:: D_c
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero and calculate current time
		Phi_Value = 0.D0
		time = tmin + DBLE(tnum)*dt
		!---------------------------------------------------------------------------------------------------------------!
		! ProblemType less than 0:	expressions with a full (spatial) taylor series
		! ProblemType equal to  0:	zero solution
		! ProblemType less than 7: 	expressions for various (spatial) powers
		! ProblemType less than 14: expressions for various (temporal) powers
		! ProblemType in the 20s:	expressions for 2 solutions meeting along an interface
		! ProblemType in the 30s:	expressions from various papers' test problems
		!---------------------------------------------------------------------------------------------------------------!
		! Full Taylor Series
		IF     (ProblemType == -3)	THEN	! Manufactured Solution - Full x series - SIN LOG
			Phi_Value = SIN(PI*x)*LOG(y+1D0)
		ELSEIF (ProblemType == -2) 	THEN	! Manufactured Solution - Full x series - EXP (space and time)
			Phi_Value = EXP(-x*y*time)
		ELSEIF (ProblemType == -1) 	THEN	! Manufactured Solution - Full x series - EXP
			Phi_Value = EXP(x*y)
		!---------------------------------------------------------------------------------------------------------------!
		ELSEIF (ProblemType == 0) 	THEN
			Phi_Value = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Powers of x, y
		ELSEIF (ProblemType == 1)	THEN	! Manufactured Solution - Powers of x - x**4
			Phi_Value = (x**4 + y**4)/12D0
		ELSEIF (ProblemType == 2)	THEN	! Manufactured Solution - Powers of x - x**3
			Phi_Value = (x**3 + y**3)/6D0
		ELSEIF (ProblemType == 3)	THEN	! Manufactured Solution - Powers of x - x**2
			Phi_Value = (x**2 + y**2)/2D0
		ELSEIF (ProblemType == 4)	THEN	! Manufactured Solution - Powers of x - x**1.5
			Phi_Value = (x**1.5D0 + y**1.5D0)*(4D0/3D0)
		ELSEIF (ProblemType == 5)	THEN	! Manufactured Solution - Powers of x - x**1
			Phi_Value = (x + y)
		ELSEIF (ProblemType == 6)	THEN	! Manufactured Solution - Powers of x - x**0.5
			Phi_Value = (x**0.5D0 + y**0.5D0)*4D0
		!---------------------------------------------------------------------------------------------------------------!
		! Powers of time
		ELSEIF (ProblemType == 7)	THEN	! Manufactured Solution - Powers of t - t**3
			Phi_Value = (time**3)/3D0
		ELSEIF (ProblemType == 8)	THEN	! Manufactured Solution - Powers of t - t**2
			Phi_Value = (time**2)/2D0
		ELSEIF (ProblemType == 9)	THEN	! Manufactured Solution - Powers of t - t**1.5
			Phi_Value = (time**1.5D0)*(2D0/3D0)
		ELSEIF (ProblemType == 10)	THEN	! Manufactured Solution - Powers of t - t**1
			Phi_Value = (time)
		ELSEIF (ProblemType == 11)	THEN	! Manufactured Solution - Powers of t - t**0.5
			Phi_Value = (time**0.5D0)*2D0
		ELSEIF (ProblemType == 12)	THEN	! Manufactured Solution - Powers of t - t**(-0.5)
			Phi_Value = 1D0/(time**0.5D0)*2D0			
		ELSEIF (ProblemType == 13)	THEN	! Manufactured Solution - Powers of t - t**(-1)
			Phi_Value = 1D0/time
		!---------------------------------------------------------------------------------------------------------------!
		! Interface of 2 solutions
		ELSEIF (ProblemType == 20) 	THEN	! Linear Interface - 4th order
			IF (y < m*x+b) THEN
				Phi_Value = x**3*(x-(y-b)/m) + alpha
			ELSE
				Phi_Value = y**3*(y-(m*x+b)) + alpha
			ENDIF
		ELSEIF (ProblemType == 21) 	THEN	! Linear Interface - Source fn is continuous
			IF (y < m*x+b) THEN
				Phi_Value =    k2*x*(y-(m*x+b)) + alpha
			ELSE
				Phi_Value = -k1*m*y*(y-(m*x+b)) + alpha
			ENDIF
		ELSEIF (ProblemType == 22) 	THEN	! Linear Interface - 0 on all boundaries
			IF (y<m*x+b) THEN
				Phi_Value =   -(y-m*x-b)*(x-1D0)*x*y*(y-1D0) + alpha
			ELSE
				Phi_Value = -(x-(y-b)/m)*(x-1D0)*x*y*(y-1D0) + alpha
			ENDIF
		ELSEIF (ProblemType == 23) 	THEN	! Linear Interface - 0.5 order (sqrt)
			IF (y < m*x+b)	THEN
				Phi_Value = m*SQRT(x+1D0) + alpha
			ELSE
				Phi_Value = (y-b+m)/SQRT(x+1D0) + alpha
			ENDIF
		ELSEIF (ProblemType == 24) 	THEN	! Circular Interface
			r = SQRT(x**2 + y**2)
			IF (r < a) THEN
				Phi_Value = COS( PI*r/(2D0*a) ) + alpha
			ELSE
				Phi_Value = beta*(r/a - 1D0)**2 + alpha
			ENDIF
		ELSEIF (ProblemType == 25) 	THEN	! Linear Interface - Sin Log
			Phi_Value = SIN(PI*x)*LOG(y+1D0)
		!---------------------------------------------------------------------------------------------------------------!
		! Analytic and test-problem solutions
		ELSEIF (ProblemType == 30) 	THEN	! Analytic solution - exponetial
			!!Phi_Value = SQRT(tmin/time)*EXP(-x**2/(4D0*K%xx(i,j)*time))*EXP(-y**2/(4D0*K%yy(i,j)*time))
			!!Phi_Value = SQRT(tmin/time)*EXP(-x**2/(4D0*time))*EXP(-y**2/(4D0*time))
			!!Phi_Value = SQRT(tmin/time)*EXP(-x**2/(4D0*time))
			!!Phi_Value = SQRT(tmin/time)*EXP(-y**2/(4D0*time))
			!Phi_Value = (tmin/time)*EXP(-x**2/(4.D0*time))*EXP(-y**2/(4.D0*time))
			Phi_Value = (tmin/time)*EXP(-x**2/(4.D0*D_c%xx(i,j)*time))*EXP(-y**2/(4.D0*D_c%yy(i,j)*time))
		ELSEIF (ProblemType == 31) 	THEN	! Continuous coefficients
			Phi_Value = (x+2D0*k1)/(1D0+4D0*k1)
		ELSEIF (ProblemType == 32) 	THEN	! Discontinuous coefficients
			IF (x < jump) 	THEN		! jump \in [0,1]
				Phi_Value = (k2*x+2D0*k1*k2)/(4D0*k1*k2+5D-1*(k1+k2))
			ELSE
				Phi_Value = (k1*x+2D0*k1*k2+jump*(k2-k1))/(4D0*k1*k2+5D-1*(k1+k2))
			ENDIF
		ELSEIF (ProblemType == 33) 	THEN	! Discontinuous coefficients [2nd order]
			IF (x < jump) 	THEN
				Phi_Value = k2*x**2
			ELSE
				Phi_Value = k1*x**2+jump**2*(k2-k1)
			ENDIF
		!---------------------------------------------------------------------------------------------------------------!
		! Mixed-Cell Problem
		ELSEIF (ProblemType == 40) 	THEN	! Mixed-Cell - Slope 0 to 10
			Phi_Value = 5D0*(x+y)
		ELSEIF (ProblemType == 41) 	THEN	! Mixed-Cell - Corner is 10, rest is 0
			IF ((x>=5D-1) .AND. (y>=5D-1)) Phi_Value = 1D1
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END FUNCTION Phi_Value

	SUBROUTINE SourceFn(Qsource,phi_o)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine defines the source function and references the value of phi at the previous time-step
	! The source function is defined as Q in d/dt (phi) + DEL(Flux) = Q
	! This can depend on time and vary with space as well. Default is 0
	!-------------------------------------------------------------------------------------------------------------------!
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT( IN)	:: phi_o	! phi old - function value at previous time step
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT)	:: Qsource	! Source function
		REAL(PRCSN), DIMENSION(xmax,ymax)				:: phi
		REAL(PRCSN)										:: x,y,time
		TYPE(Tensor_C)									:: D_c

		Qsource = 0.D0
		phi = 0.D0; x=0.D0; y=0.D0; time=0.D0
		
		IF (MMS == 1) THEN		! Control to turn on/off source loop
			time = tmin + DBLE(tnum)*dt	
			DO i=1,xmax
				DO j=1,ymax
					x = xmin + (DBLE(i)-0.5D0)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
					y = ymin + (DBLE(j)-0.5D0)*dy
					!---------------------------------------------------------------------------------------------------!
					! ProblemType less than 0:	expressions with a full (spatial) taylor series
					! ProblemType equal to  0:	zero solution
					! ProblemType less than 7: 	expressions for various (spatial) powers
					! ProblemType less than 14: expressions for various (temporal) powers
					! ProblemType in the 20s:	expressions for 2 solutions meeting along an interface
					! ProblemType in the 30s:	expressions from various papers' test problems
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
					! Interface of 2 solutions
					ELSEIF (ProblemType == 20) 	THEN	! Linear Interface - 4th order
						IF (y < m*x+b) THEN
							Qsource(i,j) = -6D0*x/m*( k1*(2D0*m*x-y+b) - k3*x )
						ELSE
							Qsource(i,j) = 6D0*y*( k2*(-2D0*y+m*x+b) + k3*y*m )
						ENDIF
					ELSEIF (ProblemType == 21) 	THEN	! Linear Interface - Source fn is continuous
						Qsource(i,j) = 2D0*k1*k2*m
					ELSEIF (ProblemType == 22) 	THEN	! Linear Interface - 0 on all boundaries
						Qsource(i,j) = 2D0*(y**2*(1D0+3D0*m*x-m+b-y)+x**2*(1D0+m*x-m+b-3D0*y)	&
											 +3D0*x*y*(1D0-m)+y*(m-b)-x*(1D0+b))
						IF (y < m*x+b) THEN
							Qsource(i,j) = Qsource(i,j)*(-k1)
							!Qsource(i,j) = -2D0*k1*(y**2*(1D0+b-m+3D0*m*x-y)+x**2*(1D0+b-m+m*x-3D0*y)	&
							!			+ 3D0*x*y*(1D0-m)+m*y-x*(1D0+b))
						ELSE
							Qsource(i,j) = Qsource(i,j)*(k2/m)
							!Qsource(i,j) = -2D0*k2/m*(y**2*(1D0+b-m+3D0*m*x-y)+x**2*(1D0+b-m+m*x-3D0*y)	&
							!			+ 3D0*y*(1D0-m)+m*y-x*(1D0+b))
						ENDIF
					ELSEIF (ProblemType == 23) 	THEN	! Linear Interface - 0.5 order (sqrt)
						IF (y < m*x+b)	THEN
							Qsource(i,j) = k1*m/4D0*(x+1D0)**(-1.5D0)
						ELSE
							Qsource(i,j) = 3D0*k3*(y-b+m)/4D0*(x+1D0)**(-2.5D0)
						ENDIF
					ELSEIF (ProblemType == 24) 	THEN	! Circular Interface
						!r = SQRT(x**2 + y**2)
						!IF (r < a)	THEN
						!	Qsource(i,j) = PI/(4D0*a**2*r**3)*( 					&
						!	PI*r*COS(PI*r/(2D0*a))*(k1*x**2+k2*y**2+2D0*k3*x*y) 	&
						!	+2D0*a*SIN(PI*r/(2D0*a))*(k1*y**2+k2*x**2-2D0*k3*x*y))							
						!ELSE
						!	Qsource(i,j) = -2D0*beta/(a**2*r**3)*(					&
						!	k1*(r**3-a*y**2) + k2*(r**3-a*x**2) + k3*(2D0*a*x*y) )
						!ENDIF
					ELSEIF (ProblemType == 25) 	THEN	! Linear Interface - Sin Log
						IF (i+j==2) CALL MixDiffusionTensor(D_c)	! call only once
						Qsource(i,j) = SIN(PI*x)*( D_c%xx(i,j)*PI**2*LOG(y+1D0)+D_c%yy(i,j)/(y+1D0)**2 )	&
									- COS(PI*x)*(2D0*PI*D_c%xy(i,j))/(y+1D0)
					! If you can't call to find D_c, define Q as below
					!IF (y < m*x+b)	THEN
					!	Qsource(i,j) = SIN(PI*x)*k1*( PI**2*LOG(y+1D0)+1D0/(y+1D0)**2 )
					!ELSEIF (y > m*x+b)	THEN
					!	Qsource(i,j) = SIN(PI*x)*k2*( PI**2*LOG(y+1D0)+1D0/(y+1D0)**2 )
					!ELSE
					!	PRINT *,i,j,'hello hello can anyone hear me out there?'
					!	!SIN(PI*x)*( k1*PI**2*LOG(y+1D0)+k2/(y+1D0)**2 ) - COS(PI*x)*(2D0*PI*k3)/(y+1D0)
					!ENDIF
					!---------------------------------------------------------------------------------------------------!
					! Analytic and test-problem solutions
					ELSEIF (ProblemType == 30) 	THEN	! Analytic solution - exponetial
						Qsource(i,j) = 0D0	! analytic solution does not need a source to work
					ELSEIF (ProblemType == 31) 	THEN	! Continuous coefficients
						Qsource(i,j) = 0D0
					ELSEIF (ProblemType == 32) 	THEN	! Discontinuous coefficients
						Qsource(i,j) = 0D0
					ELSEIF (ProblemType == 33) 	THEN	! Discontinuous coefficients [2nd order]
						Qsource(i,j) = -2D0*k1*k2
					!---------------------------------------------------------------------------------------------------!
					ELSE	! Set to zero for all other problem types
						Qsource(i,j) = 0D0
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDDO
			ENDDO
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SourceFn

	SUBROUTINE Error_Calculation(phi_c,n,error)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine calculates error of the solution WRT to analytic solution
	! The Lp error is calculated for p={1,2,inf} for each time step
	! The final variable, saved after all time steps, gives the 3 error values for each step
	! This variable is saved to file on the last time step
	!
	! phi_c = phi calculated
	! phi_a = phi analytic
	! diff  = ABS(phi_c-phi_a)
		REAL(PRCSN), DIMENSION(xmax,ymax), 	INTENT( IN)		:: phi_c
		INTEGER,							INTENT( IN)		:: n
		REAL(PRCSN), DIMENSION(tmax,3), 	INTENT(INOUT)	:: error
		REAL(PRCSN), DIMENSION(xmax,ymax)				  	:: diff, phi_a
		REAL(PRCSN)											:: x,y
		CHARACTER(LEN=17)   								:: filename
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variables to zero
		IF ((tnum==1) .OR. (n==1)) error = 0D0
		diff = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate the error
		DO i=1,xmax; DO j=1,ymax; 
			x = xmin + (DBLE(i)-0.5D0)*dx
			y = ymin + (DBLE(j)-0.5D0)*dy
			phi_a(i,j) = Phi_Value(x,y)
		ENDDO; ENDDO
		!CALL Savefile(tnum,'a',phi_a)
		CALL Savefile(n,'a',phi_a)
		
		diff = ABS(phi_c - phi_a)
		error(n,1) = SUM(diff)*dx*dy			! L1 error	- SUM |x_i|
		error(n,2) = SQRT(SUM(diff*diff)*dx*dy)	! L2 error	- SQRT SUM |x_i|^2
		error(n,3) = MAXVAL(diff)				! Linfinite error = Lmax error = MAX |x_i|

		! Note: multiplication by volume of cell is a necessary step when going from continuum to discrete
		! Each cell represent an average over some volume, whereas a (continuum) point has no volume/dimension
		! The L1 and L2 errors need the volume term, otherwise having more points in the sum will unilaterally
		! create a higher error value. Dividing by size normalizes this effect.
		! However, the Lmax term does not need the volume term because a maximum point is just a point, and 
		! it is not meant to be divided by size, since it is a pointwise value.
		! error(n,3) = MAXVAL(diff)*dx*dy
		!---------------------------------------------------------------------------------------------------------------!
		! Save the variable
		IF (tnum==tmax) 	THEN
			IF ((xmax<100) .AND. (ymax<100)) THEN
				WRITE(filename,'(A6,I2,A1,I2,A2)'),'Error_',xmax,'x',ymax,'.m'
			ELSEIF ((xmax<1000) .AND. (ymax<1000)) THEN
				WRITE(filename,'(A6,I3,A1,I3,A2)'),'Error_',xmax,'x',ymax,'.m'
			ELSEIF ((xmax<10000) .AND. (ymax<10000)) THEN
				WRITE(filename,'(A6,I4,A1,I4,A2)'),'Error_',xmax,'x',ymax,'.m'
			ENDIF
			
			CALL SaveGeneral(n,3,filename,error(1:n,:))
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Error_Calculation

	SUBROUTINE DiffusionTensor(CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine defines the diffusion tensor on the cell-centers to overlay with the initial concentration field
	! This version uses a linear interface (y=mx+b) to separate two materials, one below the line and the other above it.
	! NOTE: In the future, the Diffusion Tensor may be initialized from an outside source, say the material properties
	! 		are recovered from a volume-of-fluid, moment-of-fluid, or level set method.
	!-------------------------------------------------------------------------------------------------------------------!
		TYPE(Tensor_C),	INTENT(OUT) :: CellCentered
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		CellCentered%xx = 0D0; CellCentered%yy = 0D0; CellCentered%xy = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Specify diffusion tensor as a function of position
		IF (MixedMethod == 0)	THEN
			CellCentered%xx = k1
			CellCentered%yy = k2
			CellCentered%xy = k3			
		ELSE !IF (MixedMethod /= 0) CALL MixDiffusionTensor(CellCentered)
			!CALL MixDiffusionTensor(CellCentered)		! Domain split into 2 regions with 1 line 
			CALL MixChannels(CellCentered)				! Domain split into 4 channels with 8 lines
		ENDIF
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE DiffusionTensor

	SUBROUTINE MixDiffusionTensor(CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine is called from DiffusionTensor if there is a linear interface.
	! This subroutine defines the diffusion tensor k1 above the line, k2 below the line, 
	! where k1 and k2 are diffusion coefficients defined with all the other initial variables.
	! If MixedMethod==5, all mixed cells are calculated via rotation matrix, creating non-zero xy values
	!-------------------------------------------------------------------------------------------------------------------!
		TYPE(Tensor_C), INTENT(OUT)	:: CellCentered
		REAL(PRCSN), PARAMETER		:: LineTolx = dx*.999D0, LineToly = dy*.999D0, CrossLine = 1D-8
		REAL(PRCSN)					:: x,y
		REAL(PRCSN)					:: arth, harm, theta
		REAL(PRCSN), DIMENSION(4)	:: MCxy				! Mixed-Cell x,y [xleft, xright, ytop, ybottom]
		REAL(PRCSN), DIMENSION(6,3)	:: MCcrit			! Mixed-Cell critical points - rows: [2 intersections, 4 corners] columns: x, y, material
		REAL(PRCSN), DIMENSION(2)	:: Vec1, Vec2		! Vector 1 and Vector 2
		REAL(PRCSN), DIMENSION(2)	:: VolFrac			! Volume Fraction
		REAL(PRCSN)					:: Length1, Length2	! Length 1 and Length 2
		INTEGER						:: count, iter		! Integers used in DO-Loops
		INTEGER, DIMENSION(2)		:: quad				! Integer used for identifying the corners of a material
		INTEGER, DIMENSION(0:2)		:: MatCount			! Counter for each material
		!---------------------------------------------------------------------------------------------------------------!
		DO icell=1,xmax
			DO jcell=1,ymax
				x = xmin + (DBLE(icell)-0.5D0)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
				y = ymin + (DBLE(jcell)-0.5D0)*dy
				!-------------------------------------------------------------------------------------------------------!
				! Linear interface problems
				IF (y < m*x+b) THEN
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ELSE!IF (y > m*x+b) THEN
					CellCentered%xx(icell,jcell) = k2
					CellCentered%yy(icell,jcell) = k2
				ENDIF
				!-------------------------------------------------------------------------------------------------------!
				! Circulare interface problems
				!r = SQRT(x**2 + y**2)
				!IF (r < a) THEN
				!	CellCentered%xx(icell,jcell) = k1
				!	CellCentered%yy(icell,jcell) = k1
				!ELSEIF (r > a) THEN
				!	CellCentered%xx(icell,jcell) = k2
				!	CellCentered%yy(icell,jcell) = k2
				!ELSE
				!	print *,icell,jcell,r
				!ENDIF
				!-------------------------------------------------------------------------------------------------------!
				! Determine which cells are mixed-cells
				MCxy 	  =	(/ DBLE(icell-1)*dx, DBLE(icell)*dx, DBLE(jcell)*dy, DBLE(jcell-1)*dy /) ! MCxy=[xl, xr, yt, yb]
	 MixedLoop: IF (((	   ABS(m*MCxy(1)+b-MCxy(4))<LineToly)	.AND. (ABS(MCxy(3)-(m*MCxy(1)+b))<LineToly)) 	&		! test: cross left face
					.OR. ((ABS(m*MCxy(2)+b-MCxy(4))<LineToly)	.AND. (ABS(MCxy(3)-(m*MCxy(2)+b))<LineToly)) 	&		! test: cross right face
					.OR. ((ABS((MCxy(3)-b)/m-MCxy(1))<LineTolx) .AND. (ABS(MCxy(2)-(MCxy(3)-b)/m)<LineTolx)) 	&		! test: cross top face
					.OR. ((ABS((MCxy(4)-b)/m-MCxy(1))<LineTolx) .AND. (ABS(MCxy(2)-(MCxy(4)-b)/m)<LineTolx))	&		! test: cross bottom face
					.OR. (y == m*x+b))	THEN	! test: line passes through center of cell (rectangular cell means crosses through 2 corners)
				!IF ((a-dx/2D0 < r) .AND. (a+dx/2D0 > r)) 	THEN		! inital attempt for circular interface
				!-------------------------------------------------------------------------------------------------------!
				! Calculate volume fractions of each component of a mixed-cell
					!---------------------------------------------------------------------------------------------------!
					! Determine the 2 intersection points
					count = 0
					MatCount = 0
					MCcrit = 0D0
					IF ((m*MCxy(1)+b-MCxy(4)>0D0)   .AND. (MCxy(3)-(m*MCxy(1)+b)>0D0))	THEN
						count = count + 1
						MCcrit(count,:) = (/ MCxy(1), m*MCxy(1)+b, 0D0 /)
					ENDIF

					IF ((m*MCxy(2)+b-MCxy(4)>0D0)   .AND. (MCxy(3)-(m*MCxy(2)+b)>0D0))	THEN
						count = count + 1
						MCcrit(count,:) = (/ MCxy(2), m*MCxy(2)+b, 0D0 /)
					ENDIF

					IF (((MCxy(3)-b)/m-MCxy(1)>0D0) .AND. (MCxy(2)-(MCxy(3)-b)/m>0D0))	THEN
						count = count + 1
						MCcrit(count,:) = (/ (MCxy(3)-b)/m, MCxy(3), 0D0 /)
					ENDIF

					IF (((MCxy(4)-b)/m-MCxy(1)>0D0) .AND. (MCxy(2)-(MCxy(4)-b)/m>0D0))	THEN
						count = count + 1
						MCcrit(count,:) = (/ (MCxy(4)-b)/m, MCxy(4), 0D0 /)
					ENDIF

					IF (count > 2) PRINT '(A23,I2,A6)', 'Interface crosses cell ', count, ' times'
					!---------------------------------------------------------------------------------------------------!
					! Determine the 4 corner values and which material they are in (clockwise from bottom left)
						! Note: the 3rd column says which material the point is: 0 is shared, 1 is under line, 2 is over line
					MCcrit(3,1:2) = (/ MCxy(1),MCxy(4) /)	! left bottom
					MCcrit(4,1:2) = (/ MCxy(1),MCxy(3) /)	! left top
					MCcrit(5,1:2) = (/ MCxy(2),MCxy(3) /)	! right top
					MCcrit(6,1:2) = (/ MCxy(2),MCxy(4) /)	! right bottom
					! Determine material property for each of the 4 corners
					DO iter = 3,6
						IF 		(MCcrit(iter,1)*m+b - MCcrit(iter,2) > CrossLine)	THEN
							MCcrit(iter,3) = 1D0
							MatCount(1) = MatCount(1) + 1
						ELSEIF	(MCcrit(iter,2) - (MCcrit(iter,1)*m+b) >CrossLine)	THEN
							MCcrit(iter,3) = 2D0
							MatCount(2) = MatCount(2) + 1
						ELSE
							MCcrit(iter,3) = 0D0
							MatCount(0) = MatCount(0) + 1
							IF (count<2) THEN
								IF 		(iter==3) THEN ! left-bottom
									count = count+1
									MCcrit(count,:) = (/ MCxy(1), MCxy(4), 0D0 /)	
								ELSEIF 	(iter==4) THEN ! left-top
									count = count+1
									MCcrit(count,:) = (/ MCxy(1), MCxy(3), 0D0 /)	
								ELSEIF 	(iter==5) THEN ! right-top
									count = count+1
									MCcrit(count,:) = (/ MCxy(2), MCxy(3), 0D0 /)	
								ELSEIF 	(iter==6) THEN ! right-bottom
									count = count+1
									MCcrit(count,:) = (/ MCxy(2), MCxy(4), 0D0 /)	
								ENDIF
							ENDIF
						ENDIF
					ENDDO
					!---------------------------------------------------------------------------------------------------!
					! Use the 6 critical points (MCcrit) to determine the volume fractions for that cell
						! Note: The values of MatCount determines the shape of at least one of the materials
						! Have at least 1 triangle: [0,3,1] [0,1,3] [1,1,2] [1,2,1]
						! Have exactly 2 triangles: [2,1,1]
						! Have 2 quadrilaterals:	[0,2,2]
						! Not possible:	[2,0,2],[2,2,0],[3,1,0],[3,0,1],[4,0,0],[0,4,0],[0,0,4]
						! This only works for a single, linear interface
					! Set all variables needed to calculate volume fractions to 0
					Vec1=0D0; Vec2=0D0; Length1=0D0; Length2=0D0; quad=0; VolFrac=0D0
	
					IF ((MatCount(1) == 2) .AND. (MatCount(2) == 2))	THEN
					! This case is when the line cuts the cell into 2 quadrilaterals
					! Note: Choose material 1 as the shape to be measured
						! Find the corners of material 1
						count = 0
						DO iter=3,6
							IF (MCcrit(iter,3) == 1D0)	THEN
								count = count+1
								quad(count) = iter
							ENDIF
						ENDDO
	
						! Calculate diagonal vectors (there are two options, assume one)
						Vec1 = MCcrit(quad(1),1:2)-MCcrit(1,1:2)
						Vec2 = MCcrit(quad(2),1:2)-MCcrit(2,1:2)
						! Switch pairing if you got sides of cell instead of diagonals
						IF (((Vec1(1)==0D0) .OR. (Vec1(2)==0D0)) .AND. ((Vec2(1)==0D0) .OR. (Vec2(2)==0D0)))	THEN
							Vec1 = (MCcrit(quad(2),1:2)-MCcrit(1,1:2))
							Vec2 = (MCcrit(quad(1),1:2)-MCcrit(2,1:2))
						ENDIF
						! Determine area of quadrilateral by taking half the ABS of cross product of the diagonals
						Length1 = ABS(Vec1(1)*Vec2(2)-Vec1(2)*Vec2(1))*0.5D0
						
						VolFrac = (/Length1, dx*dy-Length1/)/(dx*dy)
					ELSEIF ((MatCount(1) == 1) .OR. (MatCount(2) == 1))	THEN
					! All cases here have at least 1 triangle. [0,3,1],[0,1,3],[1,1,2],[1,2,1]
					! This case has exactly 2 triangles of equal size: [2,1,1]
					! We find the area of triangle since the triangle area is easier to calculate
						! Find the material that has only 1 corner
						IF (MatCount(0) < 2) THEN
							! This case is for a triangle and various non-trianglular shapes (quadrilaterals or pentagons)
							IF (MatCount(1) == 1) THEN; DO iter=3,6; IF (MCcrit(iter,3) == 1D0) EXIT; ENDDO; ENDIF
							IF (MatCount(2) == 1) THEN; DO iter=3,6; IF (MCcrit(iter,3) == 2D0) EXIT; ENDDO; ENDIF
						ELSE
							! This case has 2 corners intersected by the line, which must cut the cell into identical triangles
							! I will calculate volume fractions, although I know they should both be 0.50
							DO iter=3,6; IF (MCcrit(iter,3) == 1D0) EXIT; ENDDO
						ENDIF
	
						! Measure the sides of the right triangle (this will never be the hypotenuse)
							! Vector form					
							Vec1 = ABS(MCcrit(iter,1:2)-MCcrit(1,1:2))
							Vec2 = ABS(MCcrit(iter,1:2)-MCcrit(2,1:2))
	
							!Scalar length
							Length1 = SQRT(DOT_PRODUCT(Vec1,Vec1))
							Length2 = SQRT(DOT_PRODUCT(Vec2,Vec2))
										
						VolFrac = (/Length1*Length2, 2D0*dx*dy-Length1*Length2/)/(2D0*dx*dy)
						
						! If Volume calculated was for material 2, switch the order of the VolFrac vector
						IF (MatCount(2) == 1)	VolFrac = (/VolFrac(2), VolFrac(1)/)
					ENDIF
					
					IF (ABS(SUM(VolFrac)-1D0) > (LineTolx+LineToly)*5D-1) VolFrac = VolFrac/SUM(VolFrac)
					!---------------------------------------------------------------------------------------------------!
					! Formated outputs to view the Volume Fraction and MCcrit (mixed-cell critial points)
					!WRITE(*,'(2(I3,2X),2(F10.6,2X))') icell,jcell,VolFrac
					!WRITE(*,'(3(F8.6,2X))') (MCcrit(iter,:),iter=1,6)
					!---------------------------------------------------------------------------------------------------!
					! Treatments of Mixed-Cells
	MixedTreatment: IF 		(MixedMethod == 1)	THEN	! 1=a - Arithmetic mean
						CellCentered%xx(icell,jcell) = VolFrac(1)*k1 + VolFrac(2)*k2
						CellCentered%yy(icell,jcell) = VolFrac(1)*k1 + VolFrac(2)*k2
					ELSEIF 	(MixedMethod == 2)	THEN	! 2=b - Biggest value
						IF (VolFrac(1) > VolFrac(2))		THEN	! take bigger volume fraction
							CellCentered%xx(icell,jcell) = k1
							CellCentered%yy(icell,jcell) = k1
						ELSEIF (VolFrac(2) > VolFrac(1))	THEN
							CellCentered%xx(icell,jcell) = k2
							CellCentered%yy(icell,jcell) = k2
						ELSE	! if volume fractions are equal, take bigger diffusion coefficient
							CellCentered%xx(icell,jcell) = MAX( k1 , k2 )
							CellCentered%yy(icell,jcell) = MAX( k1 , k2 )						
						ENDIF
					ELSEIF 	(MixedMethod == 3)	THEN	! 3=g - Geometric mean
						CellCentered%xx(icell,jcell) = k1**VolFrac(1) * k2**VolFrac(2)
						CellCentered%yy(icell,jcell) = k1**VolFrac(1) * k2**VolFrac(2)
					ELSEIF 	(MixedMethod == 4)	THEN	! 4=h - Harmonic mean
						CellCentered%xx(icell,jcell) = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)
						CellCentered%yy(icell,jcell) = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)
					ELSEIF 	(MixedMethod == 5)	THEN	! 5=r - Rotated value
						arth = VolFrac(1)*k1 + VolFrac(2)*k2			!arth = ( k1 + k2 ) / 2D0
						harm = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)	!harm = 2D0/( 1D0/k1 + 1D0/k2 )
						theta = atan(m)									!theta=7D0*PI/4D0	!7pi/4=315˚=-45˚=-pi/4
						CellCentered%xx(icell,jcell) = harm*COS(theta)**2 + arth*SIN(theta)**2
						CellCentered%yy(icell,jcell) = harm*SIN(theta)**2 + arth*COS(theta)**2
						CellCentered%xy(icell,jcell) = 0.5D0*(harm - arth)*SIN(2D0*theta)
						!CellCentered%xx(icell,jcell) = arth*COS(theta)**2 + harm*SIN(theta)**2
						!CellCentered%yy(icell,jcell) = arth*SIN(theta)**2 + harm*COS(theta)**2
						!CellCentered%xy(icell,jcell) = 0.5D0*(arth - harm)*SIN(2D0*theta)
					ELSEIF 	(MixedMethod == 6)	THEN	! 6=s - Smallest Value
						IF (VolFrac(1) > VolFrac(2))		THEN	! take smaller volume fraction
							CellCentered%xx(icell,jcell) = k2
							CellCentered%yy(icell,jcell) = k2
						ELSEIF (VolFrac(2) > VolFrac(1))	THEN
							CellCentered%xx(icell,jcell) = k1
							CellCentered%yy(icell,jcell) = k1
						ELSE	! if volume fractions are equal, take smaller diffusion coefficient
							CellCentered%xx(icell,jcell) = MIN( k1 , k2 )
							CellCentered%yy(icell,jcell) = MIN( k1 , k2 )						
						ENDIF
					ENDIF MixedTreatment
				ENDIF MixedLoop
			!-----------------------------------------------------------------------------------------------------------!
			ENDDO
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE MixDiffusionTensor
	
	SUBROUTINE MixChannels(CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine is called from DiffusionTensor if there is a linear interface.
	! This subroutine defines the diffusion tensor k1 above the line, k2 below the line, 
	! where k1 and k2 are diffusion coefficients defined with all the other initial variables.
	! If MixedMethod==5, all mixed cells are calculated via rotation matrix, creating non-zero xy values
	!-------------------------------------------------------------------------------------------------------------------!
		TYPE(Tensor_C), INTENT(OUT)	:: CellCentered
		REAL(PRCSN), PARAMETER		:: LineTolx = dx*.999D0, LineToly = dy*.999D0, CrossLine = 1D-8
		REAL(PRCSN)					:: x,y
		REAL(PRCSN)					:: arth, harm, theta
		REAL(PRCSN), DIMENSION(4)	:: MCxy				! Mixed-Cell x,y [xleft, xright, ytop, ybottom]
		REAL(PRCSN), DIMENSION(4)	:: MCeff			! Mixed-Cell effective values of left,right,top,bottom 
		REAL(PRCSN), DIMENSION(6,3)	:: MCcrit			! Mixed-Cell critical points - rows: [2 intersections, 4 corners] columns: x, y, material
		REAL(PRCSN), DIMENSION(2)	:: Vec1, Vec2		! Vector 1 and Vector 2
		REAL(PRCSN), DIMENSION(2)	:: VolFrac			! Volume Fraction
		REAL(PRCSN), DIMENSION(8)	:: channel			! y-intercept for lines defining channels
		REAL(PRCSN)					:: Length1, Length2	! Length 1 and Length 2
		INTEGER						:: count, iter		! Integers used in DO-Loops
		INTEGER, DIMENSION(2)		:: quad				! Integer used for identifying the corners of a material
		INTEGER, DIMENSION(0:2)		:: MatCount			! Counter for each material
		!---------------------------------------------------------------------------------------------------------------!
		DO i=1,8; channel(i)=(2D0*DBLE(i)-9D0)*1D-1;END DO	! creates values: -.7, -.5, -.3, -.1, .1, .3, .5, .7
	
		VolFrac = 0.5D0	! all mixed cells are 50/50 mixes for this type of problem.
		CellCentered%xy = 0D0
	
		DO icell=1,xmax
			DO jcell=1,ymax
				x = xmin + (DBLE(icell)-0.5D0)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
				y = ymin + (DBLE(jcell)-0.5D0)*dy
				!-------------------------------------------------------------------------------------------------------!
				! Linear interface problems
				IF 	   (y > m*x+channel(8)) THEN
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ELSEIF (y > m*x+channel(7)) THEN
					CellCentered%xx(icell,jcell) = k2
					CellCentered%yy(icell,jcell) = k2
				ELSEIF (y > m*x+channel(6)) THEN
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ELSEIF (y > m*x+channel(5)) THEN
					CellCentered%xx(icell,jcell) = k2
					CellCentered%yy(icell,jcell) = k2
				ELSEIF (y > m*x+channel(4)) THEN
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ELSEIF (y > m*x+channel(3)) THEN
					CellCentered%xx(icell,jcell) = k2
					CellCentered%yy(icell,jcell) = k2
				ELSEIF (y > m*x+channel(2)) THEN
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ELSEIF (y > m*x+channel(1)) THEN
					CellCentered%xx(icell,jcell) = k2
					CellCentered%yy(icell,jcell) = k2
				ELSE
					CellCentered%xx(icell,jcell) = k1
					CellCentered%yy(icell,jcell) = k1
				ENDIF
				!-------------------------------------------------------------------------------------------------------!
				! Determine which cells are mixed-cells
				MCxy	=	(/ DBLE(icell-1)*dx, DBLE(icell)*dx, DBLE(jcell)*dy, DBLE(jcell-1)*dy /) ! MCxy=[xl, xr, yt, yb]
				DO i=1,8	! run through each of the 8 lines
					MCeff	=	(/ (y-channel(i))/m - dy/2D0/m,(y-channel(i))/m + dy/2D0/m, &
												 m*x+channel(i)+m*dx/2D0, m*x+channel(i)-m*dx/2D0/)
					! This calculates the 'effective' positions of the cell walls assuming the line crosses through the middle
					! These should match the MCxy, which are the actual positions of the cell walls.
					! This only works with positively sloped lines, and is not as robust as MixDiffusionTesnor,
					! which finds mixed cells for any (single) line
					IF ( (MCxy(1)-MCeff(1)<LineTolx) .AND. (MCxy(2)-MCeff(2)<LineTolx) .AND. &
						 (MCxy(3)-MCeff(3)<LineToly) .AND. (MCxy(4)-MCeff(4)<LineToly)) THEN
						!PRINT '((3(I2,2X),2(F8.4,2X)))',icell,jcell,i,x,y
						!-----------------------------------------------------------------------------------------------!
						! Treatments of Mixed-Cells
		MixedTreatment: IF 		(MixedMethod == 1)	THEN	! 1=a - Arithmetic mean
							CellCentered%xx(icell,jcell) = VolFrac(1)*k1 + VolFrac(2)*k2
							CellCentered%yy(icell,jcell) = VolFrac(1)*k1 + VolFrac(2)*k2
						ELSEIF 	(MixedMethod == 2)	THEN	! 2=b - Biggest value
							IF (VolFrac(1) > VolFrac(2))		THEN	! take bigger volume fraction
								CellCentered%xx(icell,jcell) = k1
								CellCentered%yy(icell,jcell) = k1
							ELSEIF (VolFrac(2) > VolFrac(1))	THEN
								CellCentered%xx(icell,jcell) = k2
								CellCentered%yy(icell,jcell) = k2
							ELSE	! if volume fractions are equal, take bigger diffusion coefficient
								CellCentered%xx(icell,jcell) = MAX( k1 , k2 )
								CellCentered%yy(icell,jcell) = MAX( k1 , k2 )						
							ENDIF
						ELSEIF 	(MixedMethod == 3)	THEN	! 3=g - Geometric mean
							CellCentered%xx(icell,jcell) = k1**VolFrac(1) * k2**VolFrac(2)
							CellCentered%yy(icell,jcell) = k1**VolFrac(1) * k2**VolFrac(2)
						ELSEIF 	(MixedMethod == 4)	THEN	! 4=h - Harmonic mean
							CellCentered%xx(icell,jcell) = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)
							CellCentered%yy(icell,jcell) = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)
						ELSEIF 	(MixedMethod == 5)	THEN	! 5=r - Rotated value
							arth = VolFrac(1)*k1 + VolFrac(2)*k2			!arth = ( k1 + k2 ) / 2D0
							harm = 1D0/( VolFrac(1)/k1 + VolFrac(2)/k2)	!harm = 2D0/( 1D0/k1 + 1D0/k2 )
							theta = atan(m)								!theta=7D0*PI/4D0	!7pi/4=315˚=-45˚=-pi/4
							CellCentered%xx(icell,jcell) = harm*COS(theta)**2 + arth*SIN(theta)**2
							CellCentered%yy(icell,jcell) = harm*SIN(theta)**2 + arth*COS(theta)**2
							CellCentered%xy(icell,jcell) = 0.5D0*(harm - arth)*SIN(2D0*theta)
							!CellCentered%xx(icell,jcell) = arth*COS(theta)**2 + harm*SIN(theta)**2
							!CellCentered%yy(icell,jcell) = arth*SIN(theta)**2 + harm*COS(theta)**2
							!CellCentered%xy(icell,jcell) = 0.5D0*(arth - harm)*SIN(2D0*theta)
						ELSEIF 	(MixedMethod == 6)	THEN	! 6=s - Smallest Value
							IF (VolFrac(1) > VolFrac(2))		THEN	! take smaller volume fraction
								CellCentered%xx(icell,jcell) = k2
								CellCentered%yy(icell,jcell) = k2
							ELSEIF (VolFrac(2) > VolFrac(1))	THEN
								CellCentered%xx(icell,jcell) = k1
								CellCentered%yy(icell,jcell) = k1
							ELSE	! if volume fractions are equal, take smaller diffusion coefficient
								CellCentered%xx(icell,jcell) = MIN( k1 , k2 )
								CellCentered%yy(icell,jcell) = MIN( k1 , k2 )						
							ENDIF
						ENDIF MixedTreatment
					ENDIF
				ENDDO
			!-----------------------------------------------------------------------------------------------------------!
			ENDDO
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE MixChannels

	SUBROUTINE InverseFluxMatrix(Kcell)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine takes the cell-centered diffusion tensor and outputs the 4x4 system for each cell
	! that specifies the diffusion coefficient values needed to define the flux
	! Flux = - D * GRAD(phi)
	! GRAD(phi) is a 4x1 vector, with elements Left Derivative, Right Derivative, Top Derivative, and Bottom Derivative
	! Flux is a 4x1 vector, with elements f_left, f_right, f_top, and f_bottom
	! D is the diffusion coefficient, often a scalar, but in this case a tensor, having xx, xy, yx, and yy values
	! D is the necessary values to write out this system and recover flux in the form
	! f_left = Kxx(left)*Diff_Phi(left) + 0.5*(Kxy(top)*Diff_Phi(top) + Kxy(bottom)*Diff_Phi(bottom) )
	!-------------------------------------------------------------------------------------------------------------------!
		TYPE(Tensor_C)					:: D_c				! Input cell-centered diffusion tensor
		TYPE(Cell_Zone),	INTENT(OUT) :: Kcell			! 4x4 system of diff tensor for each cell
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		Kcell%cz = 0D0
		D_c%xx = 0D0; D_c%yy = 0D0; D_c%xy = 0D0
		!---------------------------------------------------------------------------------------------------------------!		
		CALL DiffusionTensor(D_c)				! Specify cell-centered values of diffusion tensor
		
		CALL SaveGeneral(xmax,ymax,'C/Dxx.m',D_c%xx)
		CALL SaveGeneral(xmax,ymax,'C/Dyy.m',D_c%yy)
		CALL SaveGeneral(xmax,ymax,'C/Dxy.m',D_c%xy)

		Kcell%cz(1,1,:,:) = 	D_c%xx		! Left face Kxx values
		Kcell%cz(3,1,:,:) =.5D0*D_c%xy		! Left face Kxy values
		Kcell%cz(4,1,:,:) =.5D0*D_c%xy		! Left face Kxy values
	 
		Kcell%cz(2,2,:,:) = 	D_c%xx		! Right face Kxx values
		Kcell%cz(3,2,:,:) =.5D0*D_c%xy		! Right face Kxy values	
		Kcell%cz(4,2,:,:) =.5D0*D_c%xy		! Right face Kxy values	
	
		Kcell%cz(3,3,:,:) = 	D_c%yy		! Top face Kyy values
		Kcell%cz(1,3,:,:) =.5D0*D_c%xy		! Top face Kxy values
		Kcell%cz(2,3,:,:) =.5D0*D_c%xy		! Top face Kxy values
	
		Kcell%cz(4,4,:,:) = 	D_c%yy		! Bottom face Kyy values
		Kcell%cz(1,4,:,:) =.5D0*D_c%xy		! Bottom face Kxy values
		Kcell%cz(2,4,:,:) =.5D0*D_c%xy		! Bottom face Kxy values

		!Kcell%cz(1,1,:,:) = 	D_f%xx_yy%xface(1:xmax,:)	! Left face Kxx values
		!Kcell%cz(3,1,:,:) =.5D0*D_f%xy_yx%xface(1:xmax,:)	! Left face Kxy values
		!Kcell%cz(4,1,:,:) =.5D0*D_f%xy_yx%xface(1:xmax,:)	! Left face Kxy values
		!
		!Kcell%cz(2,2,:,:) = 	D_f%xx_yy%xface(2:xmax+1,:)	! Right face Kxx values
		!Kcell%cz(3,2,:,:) =.5D0*D_f%xy_yx%xface(2:xmax+1,:)	! Right face Kxy values	
		!Kcell%cz(4,2,:,:) =.5D0*D_f%xy_yx%xface(2:xmax+1,:)	! Right face Kxy values	
		!
		!Kcell%cz(3,3,:,:) = 	D_f%xx_yy%yface(:,2:ymax+1)	! Top face Kyy values
		!Kcell%cz(1,3,:,:) =.5D0*D_f%xy_yx%yface(:,2:ymax+1)	! Top face Kxy values
		!Kcell%cz(2,3,:,:) =.5D0*D_f%xy_yx%yface(:,2:ymax+1)	! Top face Kxy values
		!
		!Kcell%cz(4,4,:,:) = 	D_f%xx_yy%yface(:,1:ymax)	! Bottom face Kyy values
		!Kcell%cz(1,4,:,:) =.5D0*D_f%xy_yx%yface(:,1:ymax)	! Bottom face Kxy values
		!Kcell%cz(2,4,:,:) =.5D0*D_f%xy_yx%yface(:,1:ymax)	! Bottom face Kxy values

		Kcell%cz = -2D0/(dy*dx)*Kcell%cz		! Divide by volumetric weight

	! Note: theory shows 4/dx/dy, but direct discretization shows 2/dx/dy
	! 		results using 2 are better than 4, but i must determine why
	! The answer to this is that the factor of 4 is correct, but the xx and yy values appear twice.
	! If one assumes the left/right materials are the same, but the top/bottom differ, then the
	! approximation for the integral \int \vec H \cdot \vec J dV is 
	!		 dx*dy/4*( 2(h_r*j_r+h_l*j_l) + h_t*(j_t(1)+j_t(2)) + h_b*(j_b(1)+j_b(2)) )
	!		=dx*dy/2*( (h_r*j_r+h_l*j_l) + h_t*(j_t(1)+j_t(2))/2 + h_b*(j_b(1)+j_b(2))/2 )
	! Since I already used a factor of 2 in the division by 2 (0.5D0), I must only divide by 2, not 4.
	!-------------------------------------------------------------------------------------------------------------------!	
	END SUBROUTINE InverseFluxMatrix

	SUBROUTINE Cell_To_Face_Scalar(CellCentered,FaceCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input variable that is cell-centered, and output a face-centered variable.
	! However, there are two directions that can be a face: normal to the x-direction and normal to the y-direction
	! Both directions are contained in the derived data type "Face_Centered"
	! The additional two rows in xface or additional two columns in yface are for the boundary conditions,
	! therefore, they do not necessarily depend on CellCentered values
	!
	! This verion handles an input of a scalar cell-centered variable, namely phi, the concentration.
	! The method for calculating the face is taking the arithmetic average of the cell-centered values on either side of
	! the face. The values for the boundary-faces are calculated via a linear extrapolation, although these are not needed.
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT (IN) :: CellCentered
		TYPE(Face_Centered), 				INTENT(OUT) :: FaceCentered
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		FaceCentered%Xface = 0D0; FaceCentered%Yface = 0D0
		!---------------------------------------------------------------------------------------------------------------!		
		! Calculate all interior values via arithmetic mean
		! x-face
		FaceCentered%Xface(2:xmax,:) = ( CellCentered(1:xmax-1,:) + CellCentered(2:xmax,:) ) / 2D0
		! y-face
		FaceCentered%Yface(:,2:ymax) = ( CellCentered(:,1:ymax-1) + CellCentered(:,2:ymax) ) / 2D0
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate all boundary values via a linear extrapolation from the nearest 2 interior points
		! The math comes from y-y1 = m*(x-x1), knowing that x2-x1=1 unit, and x-x1=0.5 units. The rest is trivial.
		! x-face
		FaceCentered%Xface(1     ,:) = (3D0*CellCentered(1   ,:) - CellCentered(2     ,:) ) / 2D0
		FaceCentered%Xface(xmax+1,:) = (3D0*CellCentered(xmax,:) - CellCentered(xmax-1,:) ) / 2D0
		! y-face
		FaceCentered%Yface(:,1     ) = (3D0*CellCentered(:,1   ) - CellCentered(:,2     ) ) / 2D0
		FaceCentered%Yface(:,ymax+1) = (3D0*CellCentered(:,ymax) - CellCentered(:,ymax-1) ) / 2D0
		!---------------------------------------------------------------------------------------------------------------!
		! Note that there are other options for definining the boundary terms. When the CellCentered variable is the 
		! concentration function, these terms are unimportant since they will be calculated by the boundary conditions
		! and over-written. However, when CellCentered is the diffusion coefficient, this method can be useful. One might 
		! also rather choose to take the nearest value, or simply set as zero, or perhaps use a harmonic mean somehow.
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Cell_To_Face_Scalar

	SUBROUTINE Cell_To_Face_Tensor(CellCentered,FaceCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input variable that is cell-centered, and output a face-centered variable.
	! However, there are two directions that can be a face: normal to the x-direction and normal to the y-direction
	! Both directions are contained in the derived data type "Face_Centered"
	! The additional two rows in xface or additional two columns in yface are for the boundary conditions,
	! therefore, they do not necessarily depend on Cell-Centered values
	!
	! This version handles the a tensor variable for the input, namely the diffusion tensor, K.
	! The method for calculating the face value is to take the harmonic mean of the cell-centered values on either side
	! of the face, provided the terms are non-zero. Boundary faces are simply set equal to the nearest cell-centered value.
		TYPE(Tensor_C),	INTENT (IN) :: CellCentered
		TYPE(Tensor_F),	INTENT(OUT) :: FaceCentered
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		FaceCentered%xx_yy%Xface = 0D0	! xx
		FaceCentered%xx_yy%Yface = 0D0	! yy
		FaceCentered%xy_yx%Xface = 0D0	! xy
		FaceCentered%xy_yx%Yface = 0D0	! yx
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate all interior values via the harmonic mean
		! The factor of 2 in the numerator is an artifact from using a 50/50 volume fraction

		! This snippet below won't work because I can't do operations by entire rows if there is a single element that is 0
		! Therefore, I'll have to do this element by element
		!	WHERE CellCentered%xx /= 0D0
		!		FaceCentered%xx_yy%Xface(2:xmax,:) = 2D0 / ( 1D0/CellCentered(1:xmax-1,:) + 1D0/CellCentered(2:xmax,:) )
		!	ELSEWHERE
		!		FaceCentered%xx_yy%Xface(2:xmax,:) = ( CellCentered(1:xmax-1,:) + CellCentered(2:xmax,:) ) / 2D0
		!	END WHERE

		DO i=1,xmax
			DO j=1,ymax
				IF (i /= xmax)	THEN 	! all x terms
					!---------------------------------------------------------------------------------------------------!
					! Isotropic - x				
					IF ( (CellCentered%xx(i,j) /= 0D0) .AND. (CellCentered%xx(i+1,j) /= 0D0) ) THEN
						FaceCentered%xx_yy%Xface(i+1,j) = 2D0 / ( 1D0/CellCentered%xx(i,j) + 1D0/CellCentered%xx(i+1,j) )
					ELSE	! Avoids division by zero that skews the result
						PRINT *,'does this every happen?'
						FaceCentered%xx_yy%Xface(i+1,j) = ( CellCentered%xx(i,j) + CellCentered%xx(i+1,j) ) / 2D0
					ENDIF
					
					! Anisotropic - x
					IF ( (CellCentered%xy(i,j) /= 0D0) .AND. (CellCentered%xy(i+1,j) /= 0D0) ) THEN
						FaceCentered%xy_yx%Xface(i+1,j) = 2D0 / ( 1D0/CellCentered%xy(i,j) + 1D0/CellCentered%xy(i+1,j) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xy_yx%Xface(i+1,j) = ( CellCentered%xy(i,j) + CellCentered%xy(i+1,j) ) / 2D0
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDIF

				IF (j /= ymax) THEN		! all y terms
					!---------------------------------------------------------------------------------------------------!
					! Isotropic - y
					IF ( (CellCentered%yy(i,j) /= 0D0) .AND. (CellCentered%yy(i,j+1) /= 0D0) ) THEN
						FaceCentered%xx_yy%Yface(i,j+1) = 2D0 / ( 1D0/CellCentered%yy(i,j) + 1D0/CellCentered%yy(i,j+1) )
					ELSE	! Avoids division by zero that skews the result
						PRINT *,'does this every happen?'
						FaceCentered%xx_yy%Yface(i,j+1) = ( CellCentered%yy(i,j) + CellCentered%yy(i,j+1) ) / 2D0
					ENDIF
					
					! Anisotropic - y
					IF ( (CellCentered%xy(i,j) /= 0D0) .AND. (CellCentered%xy(i,j+1) /= 0D0) ) THEN
						FaceCentered%xy_yx%Yface(i,j+1) = 2D0 / ( 1D0/CellCentered%xy(i,j) + 1D0/CellCentered%xy(i,j+1) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xy_yx%Yface(i,j+1) = ( CellCentered%xy(i,j) + CellCentered%xy(i,j+1) ) / 2D0
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDIF
			ENDDO
		ENDDO

		! These are the harmonic mean done for the entire matrix at once. 
		! The downside is that any zeros will dominate the contribution for the 2 faces it's involved with
		! FaceCentered%Xface(2:xmax,:) = 2D0 / ( 1D0/CellCentered(1:xmax-1,:) + 1D0/CellCentered(2:xmax,:) )
		! FaceCentered%Yface(:,2:ymax) = 2D0 / ( 1D0/CellCentered(:,1:ymax-1) + 1D0/CellCentered(:,2:ymax) )
		!-----------------------------------------------------------------------------------------------------------!
		! Calculate all boundary values by equating them to their nearest cell-centered value
		! x-face
			! Isotropic
			FaceCentered%xx_yy%Xface(1     ,:) = CellCentered%xx(1   ,:)
			FaceCentered%xx_yy%Xface(xmax+1,:) = CellCentered%xx(xmax,:)
			! Anisotropic
			FaceCentered%xy_yx%Xface(1     ,:) = CellCentered%xy(1   ,:)
			FaceCentered%xy_yx%Xface(xmax+1,:) = CellCentered%xy(xmax,:)
		! y-face
			! Isotropic
			FaceCentered%xx_yy%Yface(:,1     ) = CellCentered%yy(:,1   )
			FaceCentered%xx_yy%Yface(:,ymax+1) = CellCentered%yy(:,ymax)
			! Anisotropic
			FaceCentered%xy_yx%Yface(:,1     ) = CellCentered%xy(:,1   )
			FaceCentered%xy_yx%Yface(:,ymax+1) = CellCentered%xy(:,ymax)
		!---------------------------------------------------------------------------------------------------------------!
		! Note that there are other options for definining the boundary terms. When the CellCentered variable is the 
		! concentration function, these terms are unimportant since they will be calculated by the boundary conditions
		! and over-written. However, when CellCentered is the diffusion coefficient, this method can be useful. One might 
		! also rather choose to take the nearest value, or simply set as zero, or perhaps use a harmonic mean somehow.
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Cell_To_Face_Tensor

	SUBROUTINE Face_To_Cell_Average(FaceCentered,CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input variable that is face centered into a cell-centered value.
		TYPE(Face_Centered),				INTENT( IN) :: FaceCentered	! Face-Centered input function
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT) :: CellCentered	! Cell-Centered output function
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		CellCentered = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		CellCentered = ( FaceCentered%Xface(1:xmax,:) + FaceCentered%Xface(2:xmax+1,:)		&
					   + FaceCentered%Yface(:,1:ymax) + FaceCentered%Yface(:,2:ymax+1) ) /4D0
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Face_To_Cell_Average

	SUBROUTINE Face_To_Cell_Zonal(FaceVector,CellCentered,ZoneCenter)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input variable that is face centered into a cell-centered value.
	! Note the face-centered value is in vector form initially, and must be converted to the standard layout
	! ZoneCenter's 6 elements are as followed: Z(5,1:5) (first 5), and RHSelement(5) (6th element)
		REAL(PRCSN), DIMENSION(VecLength),	INTENT( IN) :: FaceVector	! Face-Centered 1D array
		REAL(PRCSN), DIMENSION(6,xmax,ymax),INTENT( IN) :: ZoneCenter	! Elements from Zonal and RHS needed for calculation
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT) :: CellCentered	! Cell-Centered output function
		TYPE(Face_Centered)								:: FaceCentered	! Face-Centered 2D array
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		CellCentered = 0D0
		FaceCentered%xface = 0D0; FaceCentered%yface = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		CALL Vector_To_Face(FaceVector,FaceCentered)	! Convert face-centered value from vector to matrix
		!---------------------------------------------------------------------------------------------------------------!
		CellCentered = ( ZoneCenter(6,:,:) - &
						(ZoneCenter(1,:,:)*FaceCentered%Xface(1:xmax,:) + ZoneCenter(2,:,:)*FaceCentered%Xface(2:xmax+1,:) &
						+ZoneCenter(4,:,:)*FaceCentered%Yface(:,1:ymax) + ZoneCenter(3,:,:)*FaceCentered%Yface(:,2:ymax+1)))&
					   /ZoneCenter(5,:,:)
		! Note that ZoneCenter(4,:,:) goes with yface(:,1:ymax) because of the order (left,right,top,bottom)
		! In retrospect, it would have been slightly more intuitive if it went (left,right,bottom,top) 
		!---------------------------------------------------------------------------------------------------------------!
		!CALL SaveGeneral(Veclength,1,'phiV.m',FaceVector)
		!CALL SaveGeneral(xmax+1,ymax,'phix.m',FaceCentered%xface)
		!CALL SaveGeneral(xmax,ymax+1,'phiy.m',FaceCentered%yface)
		!CALL SaveGeneral(xmax,ymax,'ZC.m',ZoneCenter(6,:,:))
		
		!WRITE(*,'(6(F8.5,2X))') ((ZoneCenter(:,icell,jcell),icell=1,xmax),jcell=1,ymax)
		
		!print *,'ZoneCenter max,min',MAXVAL(ZoneCenter),MINVAL(ZoneCenter)
		!print *,'FaceVector max,min',MAXVAL(FaceVector),MINVAL(FaceVector)
		!print *,'FaceCentered%xface max,min',MAXVAL(FaceCentered%xface),MINVAL(FaceCentered%xface)
		!print *,'FaceCentered%yface max,min',MAXVAL(FaceCentered%yface),MINVAL(FaceCentered%yface)
		!print *,'CellCentered max,min',MAXVAL(CellCentered),MINVAL(CellCentered)


		! Note: A straight average (like below) would not account for material properties correctly
		! CellCentered = ( FaceCentered%Xface(1:xmax,:) + FaceCentered%Xface(2:xmax+1,:)		&
		!			   + FaceCentered%Yface(:,1:ymax) + FaceCentered%Yface(:,2:ymax+1) ) /4D0
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Face_To_Cell_Zonal
	
	SUBROUTINE Vector_To_Face(Vector,FaceCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine writes the global vector of concentration values, one for each face in the system, to a variable
	! which is face-centered. The face-centered variable, say FACE, has 2 components: FACE%xface and FACE%yface.
	! The order of Vector_To_Face is x-interior, y-interior, x-bounds, y-bounds. When these are mapped to the face-centered
	! variable, they are in the standard physical layout: an array with the boundaries on the edges, containing the 
	! interior values.
	! The subroutine works complementary with Face_To_Vector
		REAL(PRCSN), DIMENSION(VecLength), 	INTENT (IN) :: Vector
		TYPE(Face_Centered),				INTENT(OUT) :: FaceCentered
		INTEGER											:: istart, iend
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero		
		FaceCentered%xface = 0D0; FaceCentered%yface = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! x interior
		istart	= 1
		iend 	= ymax*(xmax-1)
		FaceCentered%xface(2:xmax,:) = RESHAPE(Vector(istart:iend),(/xmax-1,ymax/))
		
		! y interior	
		istart	= iend+1
		iend	= iend+xmax*(ymax-1)
		FaceCentered%yface(:,2:ymax) = RESHAPE(Vector(istart:iend),(/xmax,ymax-1/))
		
		! x boundary (first row)
		istart	= iend+1
		iend	= iend+ymax
		FaceCentered%xface(1,:)		 = Vector(istart:iend)	
		
		! x boundary (last row)
		istart	= iend+1
		iend	= iend+ymax
		FaceCentered%xface(xmax+1,:) = Vector(istart:iend)
		
		! y boundary (column 1)
		istart	= iend+1
		iend	= iend+xmax
		FaceCentered%yface(:,1)		 = Vector(istart:iend)
		
		! y boundary (last column)
		istart	= iend+1
		iend	= iend+xmax
		FaceCentered%yface(:,ymax+1) = Vector(istart:iend)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Vector_To_Face

	SUBROUTINE Face_To_Vector(FaceCentered,Vector)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input variable that is face centered into a single vector
	! This will convert the face-centered concentrations into a column for the global system
		TYPE(Face_Centered),				INTENT (IN) :: FaceCentered
		REAL(PRCSN), DIMENSION(VecLength), 	INTENT(OUT) :: Vector
		INTEGER											:: istart, iend
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		Vector = 0D0		
		!---------------------------------------------------------------------------------------------------------------!
		! x interior
		istart	= 1
		iend 	= ymax*(xmax-1)
		Vector(istart:iend) = PACK(FaceCentered%xface(2:xmax,:),.TRUE.)
		 
		! y interior
		istart  = iend+1
		iend 	= iend+xmax*(ymax-1)
		Vector(istart:iend) = PACK(FaceCentered%yface(:,2:ymax),.TRUE.)
		
		! x boundary (row 1)
		istart	= iend+1
		iend	= iend+ymax
		Vector(istart:iend) = FaceCentered%xface(1,:)
		
		! x boundary (last row)
		istart	= iend+1
		iend	= iend+ymax
		Vector(istart:iend) = FaceCentered%xface(xmax+1,:)
		
		! y boundary (column 1)
		istart	= iend+1
		iend	= iend+xmax
		Vector(istart:iend) = FaceCentered%yface(:,1)
		
		! y boundary (last column)
		istart	= iend+1
		iend	= iend+xmax
		Vector(istart:iend) = FaceCentered%yface(:,ymax+1)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Face_To_Vector

	SUBROUTINE Apply_Boundary_Conditions(Zonal,RHSelement)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine applies the boundary conditions (BC) to the matrix and RHS of the matrix equation: Ax=b
	! Each cell of the zonal matrix is a 5x5 system, and each cell of the RHS matrix is a 5x1 system
	! The elements {1,2,3,4} correspond to the sides of a cell {left,right,top,bottom}, with 5 for the cell-center
	! The left-face boundaries are when i=1, j=1,ymax. Right is when i=xmax, j=1,ymax.
	! The top-face boundaries are when j=ymax, i=1,xmax. Bottom is when j=1, i=1,xmax.
	! This subroutine has 3 DO-Loops: one steps through all 4 possible boundary faces
	! and the other two step through the (i,j) values, where the ranges for iteration are specified in the first loop
	!
	! A particular cell may be called more than once if it has more than one boundary 
	! (in 2D, the 4 corners will be called twice) (in 3D, the 4 corners will be called 3x, and edges will be called 2x)
	! The particular modification for the BC will depend on the BC terms defined in the S1_Variables.F90 file (makefile)
	! The general form of the boundary condition is as follows: BC_alpha*v*Phi + BC_beta*Area(k)*(Flux,n_normal) = BC_psi
	! Where the variables are defined as follows:
	! 		Phi			= concentration function
	! 		Flux		= concentration flux (Flux = - D_{diffusion tensor} GRADIENT(concentration)
	!		v			= volumetric velocity (A(k)*Kxx/d)
	! 		n_normal 	= unit normal to the face
	!		BC_alpha	= Dirichlet boundary weight
	!		BC_beta 	= Neumann boundary weight
	!		BC_psi		= Value of weighted sum on the boundary
	!
	! Two variables pass both in and out of this subroutine:
	! 	Zonal is a (5,5,xmax,ymax) dimensional array, with the last 2 indicies representing the cell of the system
	! 		Zonal is defined prior to this subroutine, and specific values are altered to enforce BC
	! 	RHSelement is a (5,xmax,ymax) dimensional array, with the last 2 indicies representing the cell of the system
	! 		RHSelement will be mostly zeros for the first 4 elements of each cell, with nonzeros related to BC
		REAL(PRCSN), DIMENSION(5,xmax,ymax),   INTENT(INOUT):: RHSelement	! RHS for each cell system (4 faces + center)
		REAL(PRCSN), DIMENSION(5,5,xmax,ymax), INTENT(INOUT):: Zonal		! 5x5 matrix for each cell (prior to removing center)
		REAL(PRCSN), DIMENSION(4)							:: FaceArea		! note: this implies all cells are assumed to be identical rectangles
		TYPE(Tensor_C)										:: D_c			! Cell Centered Diffusion Coefficient (for Robin BCs)
		INTEGER												:: BC_face		! Boundary face [LRTB]=[1234]
		INTEGER												:: istart, iend, jstart, jend
		REAL(PRCSN)											:: time, BC_v, dummy, dumX, dumY
		!-------------------------------------------------------------------------------------------------------------------!
		! Boundary Conditions: alpha*u_concentration + beta*(w_flux,n_normal) = psi
		REAL(PRCSN)											:: BC_Diric		! Dirichelet boundary weight
		REAL(PRCSN)											:: BC_Neuma		! Neumann boundary weight
		REAL(PRCSN)											:: BC_Value		! Value of weighted sum on the boundary
		!-------------------------------------------------------------------------------------------------------------------!
		FaceArea = 0D0; BC_v   = 0D0; dummy  = 0D0
		D_c%xx = 0D0; D_c%xy = 0D0; D_c%yy = 0D0
		time = tmin + DBLE(tnum)*dt
		!-------------------------------------------------------------------------------------------------------------------!
		! Define Face Area
		FaceArea = (/dy, dy, dx, dx/)	! note: this implies all cells are assumed to be identical rectangles
										! this variable would need to have dimension(4,xmax,ymax) to vary with space
										! it would need to be a full matrix if cells were not rectangular (i.e. more sides or non 90˚ corners)
		!-------------------------------------------------------------------------------------------------------------------!
		! Check BC for both being 0 and, if necessary [Robin], load diffusion tensor
		istart=0
		DO BC_face=1,4
			! Check that at least one of BC_alpha and BC_beta are non-zero (on each side)
			IF (ABS(BC_alpha(BC_face)) + ABS(BC_beta(BC_face)) == 0D0) PRINT '(/,TR4,A37)','Boundary conditions not set correctly'	
			
			! Call the cell-centered Diffusion Tensor if we have Robin BC - D_f is not needed on boundaries
			IF ((BC_alpha(BC_face) /= 0D0) .AND. (BC_beta(BC_face) /= 0D0))	THEN
				istart = istart+1	! This counter prevents multiple calls to D_c since only 1 is needed
				IF (istart==1) CALL DiffusionTensor(D_c)
			ENDIF
		ENDDO
		!-------------------------------------------------------------------------------------------------------------------!
		! Begin the 3 loop process of enforcing BC
		! Step through all four faces and choose the proper iteration indices for each face
		DO BC_face=1,4			! Note that for a given face, only 1 of (icell,jcell) will change
			IF (BC_face==1) 		THEN	! Left Face
				istart	= 1
				iend	= 1
				jstart	= 1
				jend	= ymax
			ELSEIF (BC_face==2)		THEN	! Right Face
				istart	= xmax
				iend	= xmax
				jstart	= 1
				jend	= ymax
			ELSEIF (BC_face==3) 	THEN	! Top Face
				istart	= 1
				iend	= xmax
				jstart	= ymax
				jend	= ymax
			ELSEIF (BC_face==4) 	THEN	! Bottom Face
				istart	= 1
				iend	= xmax
				jstart	= 1
				jend	= 1
			ENDIF
			!---------------------------------------------------------------------------------------------------------------!
			! Assign appropriate BC for a particular face
			BC_Diric = BC_alpha(BC_face)
			BC_Neuma = BC_beta(BC_face)
			BC_Value = BC_psi(BC_face)
			!---------------------------------------------------------------------------------------------------------------!
			DO jcell=jstart,jend
				DO icell=istart,iend
				!-----------------------------------------------------------------------------------------------------------!
				! Boundary values for when boundary is a function
					IF (BC_Fun == 1) 	THEN
						dumX = xmin + (DBLE(icell)-0.5D0)*dx
						dumY = ymin + (DBLE(jcell)-0.5D0)*dy
						IF (BC_face == 1)	dumX = xmin + 0D0
						IF (BC_Face == 2)	dumX = xmin + (xlength+dx)
						IF (BC_Face == 3)	dumY = ymin + (ylength+dy)
						IF (BC_Face == 4)	dumY = ymin + 0D0
						BC_Value = Phi_Value(dumX,dumY)
					ENDIF
				!-----------------------------------------------------------------------------------------------------------!
					! Dirichlet
					IF ((BC_Diric /= 0D0) .AND. (BC_Neuma == 0D0))		THEN
						DO i=1,5
							IF (i/=BC_face) 							THEN
								RHSelement(i,icell,jcell) = RHSelement(i,icell,jcell) - Zonal(i,BC_face,icell,jcell)*BC_Value/BC_Diric
								! Experimental line---------------------------------------------------------------------------------------------------------!
								! Zonal(5,i,icell,jcell) = Zonal(5,i,icell,jcell) - Zonal(BC_face,i,icell,jcell)
								! Experimental line---------------------------------------------------------------------------------------------------------!
								Zonal(i,BC_face,icell,jcell) = 0D0
								Zonal(BC_face,i,icell,jcell) = 0D0
							ELSE
								dummy = BC_Diric
								RHSelement(BC_face,icell,jcell) = BC_Value/dummy
								! Experimental line---------------------------------------------------------------------------------------------------------!
								! RHSelement(BC_face,icell,jcell) = BC_Value/dummy - Zonal(BC_face,BC_face,icell,jcell)
								! Experimental line---------------------------------------------------------------------------------------------------------!
								Zonal(BC_face,BC_face,icell,jcell) = 1D0
							ENDIF
						ENDDO
					!-------------------------------------------------------------------------------------------------------!
					! Neumann
					ELSEIF ((BC_Diric == 0D0) .AND. (BC_Neuma /= 0D0))	THEN
						dummy = BC_Neuma
						RHSelement(BC_face,icell,jcell) = RHSelement(BC_face,icell,jcell) - FaceArea(BC_face)*BC_Value/dummy
					!-------------------------------------------------------------------------------------------------------!
					! Robin
					ELSEIF ((BC_Diric /= 0D0) .AND. (BC_Neuma /= 0D0))	THEN
						dummy = BC_Neuma
						IF ((BC_face == 1) .OR. (BC_face == 2))			THEN
							BC_v = FaceArea(BC_face)*D_c%xx(icell,jcell)/dummy
						ELSE
							BC_v = FaceArea(BC_face)*D_c%yy(icell,jcell)/dummy
						ENDIF
						Zonal(BC_face,BC_face,icell,jcell) = Zonal(BC_face,BC_face,icell,jcell) - BC_Diric*BC_v
						RHSelement(BC_face,icell,jcell) = RHSelement(BC_face,icell,jcell) - BC_Value*BC_v
						! BC_Value = BC_Diric/BC_Neuma*BC_v*BC_phi+FaceArea(BC_face)*BC_flux
					!-------------------------------------------------------------------------------------------------------!
					! Periodic
					!ELSEIF ((BC_Diric == 0D0) .AND. (BC_Neuma == 0D0))	THEN
					! Where (e.g.) the flux from the far left face enters into the flux on the far right face
						!RHSelement(BC_face,icell,jcell) = RHSelement(BC_face,icell,jcell) - Zonal(BC_face,BC_face,icell,jcell)
					ENDIF
					!-------------------------------------------------------------------------------------------------------!
				ENDDO
			ENDDO
			!---------------------------------------------------------------------------------------------------------------!
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Apply_Boundary_Conditions

	SUBROUTINE Zonal_System(Mflux,phi_c_old,Zelement,RHSvector,ZoneCenter)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine calculates the 5x5 zonal matrix system for each cell, as well as the 
	! corresponding right-hand-side vector, which is 5x1. It is 'b' in the matrix equation Ax=b.
	! Each cell or element's 5x5 system is defined in terms of the diffusion coefficients and the areas of the cell's faces
	! The RHS is always 0 unless a boundary condition changes this
	! The boundary conditions are applied at the end of this subroutine via calls to Apply_Boundary_Conditions for each 
	! boundary face, of which there are 2*xmax + 2*ymax
	! Note: the dependency on 5th row/column value (phi cell-center) is removed by the end of the subroutine
		TYPE(Cell_zone),					INTENT (IN)	:: Mflux		! Coefficient matrix for each cell's flux
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT (IN) :: phi_c_old	! Old value of the concentration 
		TYPE(Cell_Zone),		  			INTENT(OUT)	:: Zelement		! 4x4 matrix for each cell, flux*area
		REAL(PRCSN), DIMENSION(VecLength),	INTENT(OUT)	:: RHSvector	! Vector form of RHS for each face
		REAL(PRCSN), DIMENSION(6,xmax,ymax),INTENT(OUT) :: ZoneCenter	! Elements from Zonal and RHS needed to for Face_To_Cell
		REAL(PRCSN), DIMENSION(5,xmax,ymax)				:: RHSelement	! RHS for each cell system (4 faces + center)
		TYPE(Face_Centered)								:: RHSface		! RHS for each face (no center value)
		REAL(PRCSN), DIMENSION(5,5,xmax,ymax) 			:: Zonal		! 5x5 matrix for each cell (prior to removing center)
		REAL(PRCSN), DIMENSION(4)			  			:: FaceArea		! note: this implies all cells are assumed to be identical rectangles
		REAL(PRCSN), DIMENSION(xmax,ymax)				:: Qsource		! Time/space dependent source function
		!---------------------------------------------------------------------------------------------------------------!
		! Set inital values to zero
		RHSelement 		= 0D0	! 5x1 in each cell
		RHSface%xface 	= 0D0	! one value on each face (x)
		RHSface%yface 	= 0D0	! one value on each face (y)
		RHSvector		= 0D0	! column form of each face		
		Zonal			= 0D0	! 5x5 in each cell
		Zelement%cz 	= 0D0	! 4x4 in each cell
		ZoneCenter		= 0D0	! 6x1 in each cell (used to calculate phi_c from phi_f later)
		Qsource			= 0D0	! cell-centered source function
		!---------------------------------------------------------------------------------------------------------------!
		FaceArea = (/dy, dy, dx, dx/)
		
		CALL SourceFn(Qsource,phi_c_old)

		RHSelement(5,:,:) =-dx*dy/(dt*DiffVel)*phi_c_old - Qsource*dx*dy		
		!RHSelement(5,:,:) =-dx*dy/(dt*DiffVel)*phi_c_old - Qsource	! Define diffusion equations RHS
		! Note that this is the only time the old value of phi is needed, and at this point operator splitting should be examined
		!---------------------------------------------------------------------------------------------------------------!
		! Formated print statement of Mflux to check values		
		! WRITE(*,'(4(4(F8.5,2x),/))') (((Mflux%cz(i,1:4,icell,jcell),i=1,4),icell=1,xmax),jcell=1,ymax)
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate the zonal system for all cells simultaneously
		DO i=1,4
			DO j=1,4
				! Inner 4x4
				Zonal(i,j,:,:) = Mflux%cz(i,j,:,:)*FaceArea(i)*FaceArea(j)
				! Fifth Column, Rows 1:4
				Zonal(i,5,:,:) =  Zonal(i,5,:,:) - Mflux%cz(i,j,:,:)*FaceArea(i)*FaceArea(j)
				
				! Alternative to Zonal(5,1:4,:,:) = Zonal(1:4,5,:,:)
				! Fifth Row, Columns 1:4
				Zonal(5,j,:,:) =  Zonal(5,j,:,:) - Mflux%cz(i,j,:,:)*FaceArea(i)*FaceArea(j)
			ENDDO
			! Fifth Column, Fifth Row
			Zonal(5,5,:,:) = Zonal(5,5,:,:) - Zonal(i,5,:,:)	! each row's contribution
			!Zonal(5,5,:,:) = Zonal(5,5,:,:) - .5D0*(Zonal(i,5,:,:) + Zonal(5,i,:,:))
		ENDDO
		! Fifth Column, Fifth Row
		!DO i=1,4; Zonal(5,5,:,:) = Zonal(5,5,:,:) - .5D0*(Zonal(i,5,:,:) + Zonal(5,i,:,:)); ENDDO
		!DO i=1,4; Zonal(5,5,:,:) = Zonal(5,5,:,:) - Zonal(i,5,:,:); ENDDO

		! Cell Centered Contribution - V/dt*phi_C
		Zonal(5,5,:,:) = Zonal(5,5,:,:) - dx*dy/(dt*DiffVel)	! this assumes the volume for every cell is the same
		!---------------------------------------------------------------------------------------------------------------!
		! Enforce Boundary Conditions
		CALL Apply_Boundary_Conditions(Zonal,RHSelement)
		!---------------------------------------------------------------------------------------------------------------!
		! Store Values Necessary for Face_To_Cell computation (i.e. to calculate phi_c from phi_f)
		ZoneCenter(1:5,:,:) = Zonal(5,:,:,:)		! Last column instead: ZoneCenter(1:5,:,:) = Zonal(:,5,:,:)
		ZoneCenter(6,:,:)	= RHSelement(5,:,:)
		!---------------------------------------------------------------------------------------------------------------!
		! Remove Cell-Centered Dependency
		DO i=1,4
			DO j=1,4
				Zelement%cz(i,j,:,:) = Zonal(i,j,:,:) - Zonal(i,5,:,:)*Zonal(5,j,:,:)/Zonal(5,5,:,:)
				!Zonal(i,j,:,:) = Zonal(i,j,:,:) - Zonal(i,5,:,:)*Zonal(5,j,:,:)/Zonal(5,5,:,:)					! modify Zonal instead of Zelement
				!Zelement%cz(i,j,:,:) = Zonal(i,j,:,:) - Zonal(i,5,:,:)*ZoneCenter(j,:,:)/ZoneCenter(5,:,:)		! use ZoneCenter values
			ENDDO
			RHSelement(i,:,:) = RHSelement(i,:,:) - Zonal(i,5,:,:)/Zonal(5,5,:,:)*RHSelement(5,:,:)
			!RHSelement(i,:,:) = RHSelement(i,:,:) - Zonal(i,5,:,:)/ZoneCenter(5,:,:)*ZoneCenter(6,:,:)			! use ZoneCenter values
		ENDDO
		! Zelement%cz = Zonal(1:4,1:4,:,:) ! Only needed if loop modifies Zonal instead of Zelement
		!---------------------------------------------------------------------------------------------------------------!
		! Some useful print statements for debugging:
			!WRITE(*,'(5(F10.5,2X))') (Zonal(i,:,4,4),i=1,5)
			!WRITE(*,'(5(F10.5,3X))') (RHSelement(i,4,4),i=1,5)
			!WRITE(*,'(/,4(F10.5,2X))') (Zelement%cz(i,:,4,4),i=1,4)
			!WRITE(*,'(5(F10.5,3X))') (RHSelement(i,4,4),i=1,5)
			!DO icell=1,xmax
			!	DO jcell=1,ymax
			!		WRITE(*,'(2(I3,2X))') icell,jcell
			!		WRITE(*,'(4(4(E16.9,2X),4X,E16.9,/))') (Zelement%cz(i,1:4,icell,jcell),RHSelement(i,icell,jcell),i=1,4)
			!		!WRITE(*,'(5(5(E16.9,2X),4X,E16.9,/))') (Zonal(i,1:5,icell,jcell),RHSelement(i,icell,jcell),i=1,5)
			!	ENDDO
			!ENDDO
		!---------------------------------------------------------------------------------------------------------------!
		! Convert the RHS into a face-centered vector by adding together shared faces
		! Boundary faces have only one term contributing, but the boundary value is already there
			! Note: first argument of RHSelement: 1=Left, 2=Right, 3=Top, 4=Bottom

		! x-interior
		RHSface%xface(2:xmax,:) =-RHSelement(1,2:xmax,:) - RHSelement(2,1:xmax-1,:)
		
		! x+boundaries
		RHSface%xface(1,:) 		=-RHSelement(1,1,:)
		RHSface%xface(xmax+1,:) =-RHSelement(2,xmax,:)

		IF ((BC_alpha(1)==0D0) .AND. (BC_beta(1)==0D0) .AND. (BC_alpha(2)==0D0) .AND. (BC_beta(2)==0D0)) THEN
			! x+boundaries
			RHSface%yface(:,ymax+1)	=-RHSelement(1,1,:) - RHSelement(2,xmax,:)
			RHSface%yface(:,1)		=-RHSelement(1,1,:) - RHSelement(2,xmax,:)
		ENDIF
		
		! y+interior
		RHSface%yface(:,2:ymax)	=-RHSelement(4,:,2:ymax) - RHSelement(3,:,1:ymax-1)
		
		! y+boundaries
		RHSface%yface(:,ymax+1)	=-RHSelement(3,:,ymax)
		RHSface%yface(:,1)		=-RHSelement(4,:,1)

		IF ((BC_alpha(3)==0D0) .AND. (BC_beta(3)==0D0) .AND. (BC_alpha(4)==0D0) .AND. (BC_beta(4)==0D0)) THEN
			! y+boundaries
			RHSface%yface(:,ymax+1)	=-RHSelement(3,:,ymax) - RHSelement(4,:,1)
			RHSface%yface(:,1)		=-RHSelement(4,:,1) - RHSelement(3,:,ymax)
		ENDIF
		!---------------------------------------------------------------------------------------------------------------!
		! Package RHS from face-centered to a single column form for output
		CALL Face_To_Vector(RHSface,RHSvector)
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Zonal_System
	
	SUBROUTINE Assemble_Global_Matrix_Full(Zelement,Ztotal,Map)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine maps assembles the global matrix from the cell-wise zonal matricies
	! On input, Zelement should have eliminated the dependency of phi_center
	! Calculation will be done one cell at a time, meaning all 4 faces done in parallel.
	! The variable Row_Global maps the current face being calculated to the correct row of the global system
	! The variables Col_Blank carry the 4 indices needed to correclty map the faces of the cell to the global system
	! NOTE: The columns for the bottom and left cells are the same, hence COL_LB
	! 		The 4 directions are represented by: Col_LB, Col_R, Col_T, representing Left/Bottom, Right, and Top
		TYPE(Face_Centered),	INTENT( IN) :: Map
		TYPE(Cell_Zone), 		INTENT( IN) :: Zelement
		TYPE(Global_Zone),		INTENT(OUT) :: Ztotal
		INTEGER								:: Row_Global
		INTEGER, DIMENSION(4)				:: Col_LB,Col_R,Col_T
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize Variables to zero
		Ztotal%gz = 0D0
		Row_Global= 0D0; Col_LB = 0D0; Col_R = 0D0; Col_T = 0D0
		!---------------------------------------------------------------------------------------------------------------!
		! Go to each cell and equate the appropriate faces
		! NOTE: face number {1,2,3,4} corresponds to {Left, Right, Top, Bottom}
		DO i=1,xmax
			DO j=1,ymax
			!-----------------------------------------------------------------------------------------------------------!
				! Columns for left and bottom cell
				Col_LB = (/Map%Xface(i,j),   Map%Xface(i+1,j),   Map%Yface(i,j+1),   Map%Yface(i,j)  /)
			!-----------------------------------------------------------------------------------------------------------!
				! Interior x: right face(2) of left cell(i,j) + left face(1) of right cell(i+1,j)
				IF (i<xmax)			THEN				
					! Current Row of Global Matrix for Left+Right=0
					Row_Global = Map%Xface(i+1,j)	
					! Columns for right cell (i+1,j)
					Col_R = (/Map%Xface(i+1,j), Map%Xface(i+2,j),   Map%Yface(i+1,j+1), Map%Yface(i+1,j)/)
					! Right face(2) of left cell(i,j)
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(2,1:4,i,j)
					! Left face(1) of right cell(i+1,j)
					Ztotal%gz(Row_Global,Col_R)  = Ztotal%gz(Row_Global,Col_R)  - Zelement%cz(1,1:4,i+1,j)
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Boundaries - X Face	
				IF (i==1) 			THEN
					! Current Row of Global Matrix for Left-Face (minimum row)
					Row_Global = Map%Xface(1,j)
					! Boundary x: left face(1) of first cell(1,j)							
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(1,1:4,1,j)
				ELSEIF (i==xmax)	THEN
					! Current Row of Global Matrix for Right-Face (maximum row)
					Row_Global = Map%Xface(xmax+1,j)
					! Boundary x: right face(2) of last cell(xmax,j)
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(2,1:4,xmax,j)
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Interior y: top face(3) of bottom cell(i,j) + bottom face(4) of top cell(i,j+1)
				IF (j<ymax)			THEN
					! Current Row of Global Matrix for Top+Bottom=0
					Row_Global = Map%Yface(i,j+1)
					! Columns for top cell (i,j+1)
					Col_T  = (/Map%Xface(i,j+1), Map%Xface(i+1,j+1), Map%Yface(i,j+2),   Map%Yface(i,j+1)/)
					! Top face(3) of bottom cell(i,j)
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(3,1:4,i,j)
					! Bottom face(4) top cell(i,j+1)
					Ztotal%gz(Row_Global,Col_T)  = Ztotal%gz(Row_Global,Col_T) - Zelement%cz(4,1:4,i,j+1)
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Boundaries - Y Face
				IF (j==1) 	THEN
					! Current Row of Global Matrix for Bottom-Face (minimum column)
					Row_Global = Map%Yface(i,1)
					! Boundary y: bottom face(4) of first cell(i,1)
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(4,1:4,i,1)
				ELSEIF (j==ymax) 		THEN
					! Current Row of Global Matrix for Top-Face (maximum column)
					Row_Global = Map%Yface(i,ymax+1)
					! Boundary y: top face(3) of last cell(i,ymax)
					Ztotal%gz(Row_Global,Col_LB) = Ztotal%gz(Row_Global,Col_LB) - Zelement%cz(3,1:4,i,ymax)
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
			ENDDO
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Assemble_Global_Matrix_Full

	SUBROUTINE Assemble_Global_Matrix_Sparse(Zelement,Ztotal,Map)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine maps assembles the global matrix from the cell-wise zonal matricies
	! On input, Zelement should have eliminated the dependency of phi_center
	! Calculation will be done one cell at a time, meaning all 4 faces done in parallel.
	! The variable Row_Global maps the current face being calculated to the correct row of the global system
	! The variables Col_Blank carry the 4 indices needed to correclty map the faces of the cell to the global system
	! NOTE: The columns for the bottom and left cells are the same, hence COL_LB
	! 		The 4 directions are represented by: Col_LB, Col_R, Col_T, representing Left/Bottom, Right, and Top
		TYPE(Face_Centered),	INTENT( IN) :: Map
		TYPE(Cell_Zone), 		INTENT( IN) :: Zelement
		TYPE(Sparse_Matrix),	INTENT(OUT) :: Ztotal
		INTEGER								:: Row_Global
		INTEGER, DIMENSION(4)				:: Col_LB,Col_R,Col_T
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize Variables to zero
		Ztotal%sgz = 0D0
		Row_Global= 0; Col_LB = 0; Col_R = 0; Col_T = 0
		i0 = 1		
		!---------------------------------------------------------------------------------------------------------------!
		! Go to each cell and equate the appropriate faces
		! NOTE: face number {1,2,3,4} corresponds to {Left, Right, Top, Bottom}
		DO i=1,xmax
			DO j=1,ymax
			!-----------------------------------------------------------------------------------------------------------!
				! Columns for left and bottom cell
				Col_LB = (/Map%Xface(i,j),   Map%Xface(i+1,j),   Map%Yface(i,j+1),   Map%Yface(i,j)  /)
			!-----------------------------------------------------------------------------------------------------------!
				! Interior x: right face(2) of left cell(i,j) + left face(1) of right cell(i+1,j)
				IF (i<xmax)			THEN				
					! Current Row of Global Matrix for Left+Right=0
					Row_Global = Map%Xface(i+1,j)	
					! Columns for right cell (i+1,j)
					Col_R = (/Map%Xface(i+1,j), Map%Xface(i+2,j),   Map%Yface(i+1,j+1), Map%Yface(i+1,j)/)
					! Right face(2) of left cell(i,j)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(2,1:4,i,j)
					i0=i0+4
					! Left face(1) of right cell(i+1,j)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_R
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(1,1:4,i+1,j)
					i0=i0+4
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Boundaries - X Face	
				IF (i==1) 			THEN
					! Current Row of Global Matrix for Left-Face (minimum row)
					Row_Global = Map%Xface(1,j)
					! Boundary x: left face(1) of first cell(1,j)							
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(1,1:4,1,j)
					i0=i0+4
				ELSEIF (i==xmax)	THEN
					! Current Row of Global Matrix for Right-Face (maximum row)
					Row_Global = Map%Xface(xmax+1,j)
					! Boundary x: right face(2) of last cell(xmax,j)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(2,1:4,xmax,j)
					i0=i0+4
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Interior y: top face(3) of bottom cell(i,j) + bottom face(4) of top cell(i,j+1)
				IF (j<ymax)			THEN
					! Current Row of Global Matrix for Top+Bottom=0
					Row_Global = Map%Yface(i,j+1)
					! Columns for top cell (i,j+1)
					Col_T  = (/Map%Xface(i,j+1), Map%Xface(i+1,j+1), Map%Yface(i,j+2),   Map%Yface(i,j+1)/)
					! Top face(3) of bottom cell(i,j)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(3,1:4,i,j)
					i0=i0+4
					! Bottom face(4) top cell(i,j+1)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_T
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(4,1:4,i,j+1)
					i0=i0+4
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
				! Boundaries - Y Face
				IF (j==1) 	THEN
					! Current Row of Global Matrix for Bottom-Face (minimum column)
					Row_Global = Map%Yface(i,1)
					! Boundary y: bottom face(4) of first cell(i,1)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(4,1:4,i,1)
					i0=i0+4

				ELSEIF (j==ymax) 		THEN
					! Current Row of Global Matrix for Top-Face (maximum column)
					Row_Global = Map%Yface(i,ymax+1)
					! Boundary y: top face(3) of last cell(i,ymax)
					Ztotal%sgz(i0:i0+3,1)=Row_Global
					Ztotal%sgz(i0:i0+3,2)=Col_LB
					Ztotal%sgz(i0:i0+3,3)= -Zelement%cz(3,1:4,i,ymax)
					i0=i0+4
				ENDIF
			!-----------------------------------------------------------------------------------------------------------!
			ENDDO
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Assemble_Global_Matrix_Sparse	
!-----------------------------------------------------------------------------------------------------------------------!
END MODULE SOM_Subroutine_Module
