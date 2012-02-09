MODULE TYPES
! This module introduces the variables defining the mesh sizes as well the derived data types.
! The grid parameters are set here so that the TYPES can be defined
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
    !	xx - K%xx_yy%Xface
    !	yy - K%xx_yy%Yface
    !	xy - K%xy_yx%Xface
    !	yx - K%xy_yx%Yface
    	TYPE(Face_Centered)	:: xx_yy, xy_yx
    END TYPE Tensor_F
    !-------------------------------------------------------------------------------------------------------------------!
END MODULE TYPES

MODULE SubroutineModules
! This module contains all the subroutines called in the program SOM.
! This module also calls module Types to define several derived data types.
! CONVENTION: FORTRAN COMMANDS ARE ALL CAPS; variables I define may be lower case or MiXeD
! CONVENTION: Separate separate sections with line breaks, 
! CONVENTION: and try to bind sections together with !------! to fill an intire line
    USE DefineParameters
    USE Types

	IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------!
	CHARACTER(LEN=12)   :: filename
    !-------------------------------------------------------------------------------------------------------------------!
    ! Variables are defined in DefineParameters.F90 which loads the module defined here
        ! Scalar constant in original equation:     1/DiffVel*d(u)/dt + DEL(w)=f    
        ! Boundary condition constants:             BC_alpha*u+BC_beta*(W,normal)=BC_psi
        ! grid spacing:                             dx, dy
        ! time step spacing:                        dt
        ! physical system dimensions:               xlength, ylength
        ! total simulated time:                     tlength
        ! total number of points used:              xmax, ymax, tmax
        ! initial time (to avoid division by 0):    tmin
    !-------------------------------------------------------------------------------------------------------------------!

	CONTAINS

    SUBROUTINE When
	!-------------------------------------------------------------------------------------------------------------------!
	! Prints current time and date to screen when called
        INTEGER, DIMENSION(8)   :: values
        CALL DATE_AND_TIME(VALUES=values)

        IF (values(5) <= 12) PRINT '(A6,I2.2,A1,I2.2,A1,I4,A9,I2,A1,I2.2,A1,I2.2,A3)', 'Date: ',        &
                values(2),'/',values(3),'/',values(1),'  Time: ',values(5),':',values(6),':',values(7),' AM'        


        IF (values(5) >  12) PRINT '(A6,I2.2,A1,I2.2,A1,I4,A9,I2,A1,I2.2,A1,I2.2,A3)', 'Date: ',        &
                values(2),'/',values(3),'/',values(1),'  Time: ',values(5)-12,':',values(6),':',values(7),' PM'
	!-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE When
    
	SUBROUTINE Namefile(k,prefix,filename)
	!-------------------------------------------------------------------------------------------------------------------!
	! Names each file to be save. Files are named by the iteration count (k)
		INTEGER, 			INTENT( IN) :: k
		CHARACTER(LEN=1),	INTENT( IN) :: prefix
		CHARACTER(LEN=12),	INTENT(OUT) :: filename

		IF (k<10)						THEN					! File number is 1 digit
			WRITE(filename,'(a3,a1,a4,i1,a2)'), "C//",prefix,"0000",k,".m"
		ELSE IF (abs(k-54.5)<45.5)		THEN					! File number is 2 digits
			WRITE(filename,'(a3,a1,a3,i2,a2)'), "C//",prefix,"000", k,".m"
		ELSE IF (abs(k-549.5)<450.5)	THEN					! File number is 3 digits
			WRITE(filename,'(a3,a1,a2,i3,a2)'), "C//",prefix,"00",  k,".m"
		ELSE IF (abs(k-5499.5)<4500.5)	THEN					! File number is 4 digits
			WRITE(filename,'(a3,a1,a1,i4,a2)'), "C//",prefix,"0",   k,".m"
		ELSE IF (abs(k-54999.5)<45000.5) THEN					! File number is 5 digits
			WRITE(filename,'(a3,a1,i5,a2)'),    "C//",prefix,       k,".m"
		ELSE
			PRINT *, "The number of points is greater than 5 digits"
			STOP
		ENDIF 
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Namefile
	
	SUBROUTINE Savefile(k,prefix,u)
	!-------------------------------------------------------------------------------------------------------------------!
	! Saves the current instance of u. File names are appended based on the time iteration (k)
	! The exponentials are restricted to 2 digits (hence the variable u_restricted). 
	! This restriction is needed since my graphic software (Octave) fails to properly read exponents with >2 digits
		INTEGER, 					    	INTENT( IN) :: k
		CHARACTER(LEN=1),					INTENT( IN) :: prefix
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT( IN) :: u
		REAL(PRCSN), DIMENSION(xmax,ymax)				:: u_restricted
		CHARACTER (LEN=20)                             	:: frmt

        WRITE(frmt,'(A1,I4,A13)') '(',ymax,'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))
	    ! Meaning: ymax gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that	    
        u_restricted = u						! Leave input matrix unchanged while modifying values for saving
        !---------------------------------------------------------------------------------------------------------------!
		! Ensure exponent is only 2 digits (further restricted to exponent above or below 50)
        WHERE (ABS(u_restricted) < 1.0D-50)                     ! If number is ridiculously small (+/-)
            u_restricted = 0D0
        ELSEWHERE (u_restricted > 9.9D+50)                      ! If number is ridiculously positive
            u_restricted = 9.9D+50
        ELSEWHERE (u_restricted < -9.9D+50)                     ! If number is ridiculously negative
            u_restricted = -9.9D+50
        END WHERE
		!---------------------------------------------------------------------------------------------------------------!
	    
		CALL Namefile(k,prefix,filename)                               	! Call subroutine to generate sequential file name
		OPEN (unit=1, file=filename)                            	! Open said file
		WRITE(1, '(A7,D15.7,TR5,A4,I5/)') '% dt = ', dt, 'k = ',k   ! Set first line to give time-step
        WRITE(1,frmt) (u_restricted(i,:),i=1,xmax)       			! Write variable to file
		CLOSE (1)													! Close said file
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Savefile

	SUBROUTINE SaveGeneral(xsize,ysize,filename,matrix)
	!-------------------------------------------------------------------------------------------------------------------!
	! Subroutine writes array to file, where user inputs the dimensions of the matrix, the filename, and the matrix.
	! The advantage to this is that this same subroutine can handle an array of any size.
	! The file type is expected to be '.m', which is an Octave extension, although it is equivalent to a '.txt'.  
	! The file will have the same layout as the actual array, in terms of rows and columns
		INTEGER,								INTENT (IN) :: xsize,ysize
		CHARACTER (LEN=*),                  	INTENT (IN) :: filename
        REAL(PRCSN), DIMENSION(xsize,ysize),	INTENT (IN) :: matrix
		REAL(PRCSN), DIMENSION(xsize,ysize)        		    :: u_restricted
		CHARACTER (LEN=20)									:: frmt

        WRITE(frmt,'(A1,I3,A13)') '(',ysize,'(E16.9E2,2x))'  ! Write a format string: '(99(E16.9E2,2x))
	    ! Meaning: ymax gives, say, 99 cases of a number (E16.9D2) followed by two spaces
	    ! After all columns of that row are done, the implied do loop moves onto the next row
	    ! E16.9 gives a single precision real number, exponential notation with 16 spaces and 9 dedicated 
	    ! to expressing decimal values. The E16.9E2 limits the exponent to having 2 digits or less (needed for Octave)
	    ! Note: double precision would be (D16.9E2,2x), but Octave doesn't seem to read that
	    
        u_restricted = matrix
        !---------------------------------------------------------------------------------------------------------------!
		! Ensure exponent is only 2 digits (further restricted to exponent above or below 50)
        WHERE (ABS(u_restricted) < 1.0D-50)                     ! If number is ridiculously small (+/-)
            u_restricted = 0.0D0
        ELSEWHERE (u_restricted > 9.9D+50)                      ! If number is ridiculously positive
            u_restricted = 9.9D+50
        ELSEWHERE (u_restricted < -9.9D+50)                     ! If number is ridiculously negative
            u_restricted = -9.9D+50
        END WHERE
		!---------------------------------------------------------------------------------------------------------------!
		OPEN (unit=1, file=filename)                            ! Open said file

        WRITE(1,frmt) (u_restricted(i,:),i=1,xsize)				! Write variable to file
       
		CLOSE (1)												! Close said file
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE SaveGeneral

    SUBROUTINE LoadFromFile(filename,variable)
    !-------------------------------------------------------------------------------------------------------------------!
    ! Loads file 'filename' and outputs results to variable 'variable'. 
    ! This does not load any file, only those saved with the format that I use in my SaveFile subroutine.
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(OUT)	:: variable
		CHARACTER (LEN=*),                  INTENT( IN)	:: filename
		CHARACTER (LEN=20)                      		:: frmt

        variable	= 0D0

		OPEN (unit=1, file=filename)
        WRITE(frmt,'(A3,I4,A13)') '(//',ymax,'(E16.9E2,2x))'	! Write a format string: '(####(E16.9E2,2x))
        READ(1,frmt),(variable(i,:),i=1,xmax)
		CLOSE (1)
		!variable = TRANSPOSE(variable) 	! I don't think a transpose is necessary
    !-------------------------------------------------------------------------------------------------------------------!
    END SUBROUTINE LoadFromFile

	SUBROUTINE DiffusionTensor(Kdiff)
	!-------------------------------------------------------------------------------------------------------------------!
	! Subroutine gives values to the 2D diffusion tensor, K. 
	! The assumption is that Kyx = Kxy in all cases, but that my not necessarily always be true.
    ! At some point, this will be time dependent, and will give specific values to specific cells
    ! based on the changing material locations
        Type(Tensor_C),  INTENT(OUT) :: Kdiff
        !---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		Kdiff%xx = 0D0; Kdiff%yy = 0D0; Kdiff%xy = 0D0		
		!---------------------------------------------------------------------------------------------------------------!
		! Specify diffusion tensor as a function of position
		Kdiff%xx = d1
		Kdiff%yy = d2
		Kdiff%xy = d3
		
		IF (MixedMethod /= 0) CALL MixDiffusionTensor(Kdiff)
    !-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE DiffusionTensor

	SUBROUTINE MixDiffusionTensor(CellCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine is called from DiffusionTensor if there is a linear interface.
	! This subroutine defines the diffusion tensor k1 above the line, k2 below the line, 
	! where k1 and k2 are diffusion coefficients defined with all the other initial variables.
	! If MixedMethod==5, all mixed cells are calculated via rotation matrix, creating non-zero xy values
	!-------------------------------------------------------------------------------------------------------------------!
		TYPE(Tensor_C)				:: CellCentered
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

		DO i=1,xmax
			DO j=1,ymax
				x = xmin + (DBLE(i)-0.5D0)*dx	! note: (i-0.5)dx = (i-1)dx+dx/2
				y = ymin + (DBLE(j)-0.5D0)*dy
				!-------------------------------------------------------------------------------------------------------!
				! Linear interface problems
				IF (y < m*x+b) THEN
					CellCentered%xx(i,j) = d1
					CellCentered%yy(i,j) = d1
				ELSE!IF (y > m*x+b) THEN
					CellCentered%xx(i,j) = d2
					CellCentered%yy(i,j) = d2
				ENDIF
				!-------------------------------------------------------------------------------------------------------!
				! Circulare interface problems
				!r = SQRT(x**2 + y**2)
				!IF (r < a) THEN
				!	CellCentered%xx(i,j) = d1
				!	CellCentered%yy(i,j) = d1
				!ELSEIF (r > a) THEN
				!	CellCentered%xx(i,j) = d2
				!	CellCentered%yy(i,j) = d2
				!ELSE
				!	print *,i,j,r
				!ENDIF
				!-------------------------------------------------------------------------------------------------------!
				IF 	(MixedMethod == 5)	THEN
					! Determine which cells are mixed-cells
					MCxy 	  =	(/ DBLE(i-1)*dx, DBLE(i)*dx, DBLE(j)*dy, DBLE(j-1)*dy /) ! MCxy=[xl, xr, yt, yb]
					IF ( 	 ((ABS(m*MCxy(1)+b-MCxy(4))<LineToly)	.AND. (ABS(MCxy(3)-(m*MCxy(1)+b))<LineToly)) 	&		! test: cross left face
						.OR. ((ABS(m*MCxy(2)+b-MCxy(4))<LineToly)	.AND. (ABS(MCxy(3)-(m*MCxy(2)+b))<LineToly)) 	&		! test: cross right face
						.OR. ((ABS((MCxy(3)-b)/m-MCxy(1))<LineTolx) .AND. (ABS(MCxy(2)-(MCxy(3)-b)/m)<LineTolx)) 	&		! test: cross top face
						.OR. ((ABS((MCxy(4)-b)/m-MCxy(1))<LineTolx) .AND. (ABS(MCxy(2)-(MCxy(4)-b)/m)<LineTolx))	&		! test: cross bottom face
						.OR. (y == m*x+b))	THEN	! test: line passes through center of cell (rectangular cell means crosses through 2 corners)
					!IF ((a-dx/2.D0 < r) .AND. (a+dx/2.D0 > r)) 	THEN				
					!IF ((m*x+b-dy/2.D0 < y) .AND. (y < m*x+b+dy/2.D0))	THEN	! We have a mixed cell
				!-------------------------------------------------------------------------------------------------------!
				! Calculate volume fractions of each component of a mixed-cell
						!---------------------------------------------------------------------------------------------------!
						! Determine the 2 intersection points
						count = 0
						MatCount = 0
						MCcrit = 0.D0
						IF ((m*MCxy(1)+b-MCxy(4)>0.D0)   .AND. (MCxy(3)-(m*MCxy(1)+b)>0.D0))	THEN
							count = count + 1
							MCcrit(count,:) = (/ MCxy(1), m*MCxy(1)+b, 0.D0 /)
						ENDIF
	
						IF ((m*MCxy(2)+b-MCxy(4)>0.D0)   .AND. (MCxy(3)-(m*MCxy(2)+b)>0.D0))	THEN
							count = count + 1
							MCcrit(count,:) = (/ MCxy(2), m*MCxy(2)+b, 0.D0 /)
						ENDIF
	
						IF (((MCxy(3)-b)/m-MCxy(1)>0.D0) .AND. (MCxy(2)-(MCxy(3)-b)/m>0.D0))	THEN
							count = count + 1
							MCcrit(count,:) = (/ (MCxy(3)-b)/m, MCxy(3), 0.D0 /)
						ENDIF
	
						IF (((MCxy(4)-b)/m-MCxy(1)>0.D0) .AND. (MCxy(2)-(MCxy(4)-b)/m>0.D0))	THEN
							count = count + 1
							MCcrit(count,:) = (/ (MCxy(4)-b)/m, MCxy(4), 0.D0 /)
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
								MCcrit(iter,3) = 1.D0
								MatCount(1) = MatCount(1) + 1
							ELSEIF	(MCcrit(iter,2) - (MCcrit(iter,1)*m+b) >CrossLine)	THEN
								MCcrit(iter,3) = 2.D0
								MatCount(2) = MatCount(2) + 1
							ELSE
								MCcrit(iter,3) = 0.D0
								MatCount(0) = MatCount(0) + 1
								IF (count<2) THEN
									IF 		(iter==3) THEN ! left-bottom
										count = count+1
										MCcrit(count,:) = (/ MCxy(1), MCxy(4), 0.D0 /)	
									ELSEIF 	(iter==4) THEN ! left-top
										count = count+1
										MCcrit(count,:) = (/ MCxy(1), MCxy(3), 0.D0 /)	
									ELSEIF 	(iter==5) THEN ! right-top
										count = count+1
										MCcrit(count,:) = (/ MCxy(2), MCxy(3), 0.D0 /)	
									ELSEIF 	(iter==6) THEN ! right-bottom
										count = count+1
										MCcrit(count,:) = (/ MCxy(2), MCxy(4), 0.D0 /)	
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
						Vec1=0D0; Vec2=0D0; Length1=0D0; Length2=0D0; quad=0; VolFrac=0.D0
		
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
							IF (((Vec1(1)==0.D0) .OR. (Vec1(2)==0.D0)) .AND. ((Vec2(1)==0.D0) .OR. (Vec2(2)==0.D0)))	THEN
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
								IF (MatCount(1) == 1) THEN; DO iter=3,6; IF (MCcrit(iter,3) == 1.D0) EXIT; ENDDO; ENDIF
								IF (MatCount(2) == 1) THEN; DO iter=3,6; IF (MCcrit(iter,3) == 2.D0) EXIT; ENDDO; ENDIF
							ELSE
								! This case has 2 corners intersected by the line, which must cut the cell into identical triangles
								! I will calculate volume fractions, although I know they should both be 0.50
								DO iter=3,6; IF (MCcrit(iter,3) == 1.D0) EXIT; ENDDO
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
						
						IF (ABS(SUM(VolFrac)-1.D0) > (LineTolx+LineToly)*5D-1) VolFrac = VolFrac/SUM(VolFrac)
					ENDIF
					!-------------------------------------------------------------------------------------------------------!
					! Formated outputs to view the Volume Fraction and MCcrit (mixed-cell critial points)
					!WRITE(*,'(2(I3,2X),2(F10.6,2X))') i,j,VolFrac
					!WRITE(*,'(3(F8.6,2X))') (MCcrit(iter,:),iter=1,6)
					!-------------------------------------------------------------------------------------------------------!
					! Treatments of Mixed-Cells
					IF 		(MixedMethod == 1)	THEN	! 1=a - Arithmetic mean
						CellCentered%xx(i,j) = VolFrac(1)*d1 + VolFrac(2)*d2
						CellCentered%yy(i,j) = VolFrac(1)*d1 + VolFrac(2)*d2
					ELSEIF 	(MixedMethod == 2)	THEN	! 2=b - Biggest value
						IF (VolFrac(1) > VolFrac(2))		THEN	! take bigger volume fraction
							CellCentered%xx(i,j) = d1
							CellCentered%yy(i,j) = d1
						ELSEIF (VolFrac(2) > VolFrac(1))	THEN
							CellCentered%xx(i,j) = d2
							CellCentered%yy(i,j) = d2
						ELSE	! if volume fractions are equal, take bigger diffusion coefficient
							CellCentered%xx(i,j) = MAX( d1 , d2 )
							CellCentered%yy(i,j) = MAX( d1 , d2 )						
						ENDIF
					ELSEIF 	(MixedMethod == 3)	THEN	! 3=g - Geometric mean
						! Not sure if it is better to take the volume-weighted mean
						! Might need to use SQRT( VolFrac(1)*d1 * VolFrac(2)*d2 )
						CellCentered%xx(i,j) = SQRT( d1 * d2 )
						CellCentered%yy(i,j) = SQRT( d1 * d2 )
					ELSEIF 	(MixedMethod == 4)	THEN	! 4=h - Harmonic mean
						CellCentered%xx(i,j) = 1.D0/( VolFrac(1)/d1 + VolFrac(2)/d2)
						CellCentered%yy(i,j) = 1.D0/( VolFrac(1)/d1 + VolFrac(2)/d2)
					ELSEIF 	(MixedMethod == 5)	THEN	! 5=r - Rotated value
						arth = VolFrac(1)*d1 + VolFrac(2)*d2			!arth = ( d1 + d2 ) / 2.D0
						harm = 1.D0/( VolFrac(1)/d1 + VolFrac(2)/d2)	!harm = 2.D0/( 1.D0/d1 + 1.D0/d2 )
						theta = atan(m)									!theta=7.D0*PI/4.D0	!7pi/4=315˚=-45˚=-pi/4
						CellCentered%xx(i,j) = arth*COS(theta)**2 + harm*SIN(theta)**2
						CellCentered%yy(i,j) = arth*SIN(theta)**2 + harm*COS(theta)**2
						CellCentered%xy(i,j) = 0.5D0*(arth - harm)*SIN(2.D0*theta)
						! Set values equal to some number to check mixed-cell placement
						!   CellCentered%xx(i,j) = 55
						!   CellCentered%yy(i,j) = 55
						!   CellCentered%xy(i,j) = 55
					ELSEIF 	(MixedMethod == 6)	THEN	! 6=s - Smallest Value
						IF (VolFrac(1) > VolFrac(2))		THEN	! take smaller volume fraction
							CellCentered%xx(i,j) = d2
							CellCentered%yy(i,j) = d2
						ELSEIF (VolFrac(2) > VolFrac(1))	THEN
							CellCentered%xx(i,j) = d1
							CellCentered%yy(i,j) = d1
						ELSE	! if volume fractions are equal, take smaller diffusion coefficient
							CellCentered%xx(i,j) = MIN( d1 , d2 )
							CellCentered%yy(i,j) = MIN( d1 , d2 )						
						ENDIF
					ENDIF
				ENDIF
				!-------------------------------------------------------------------------------------------------------!
			ENDDO
		ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE MixDiffusionTensor
		
	SUBROUTINE Cell_To_Face(CellCentered,FaceCentered)
	!-------------------------------------------------------------------------------------------------------------------!
	! This subroutine converts an input tensor variable that is cell-centered, and output a tensor face-centered variable.
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
		FaceCentered%xx_yy%Xface = 0.D0	! xx
		FaceCentered%xx_yy%Yface = 0.D0	! yy
		FaceCentered%xy_yx%Xface = 0.D0	! xy
		FaceCentered%xy_yx%Yface = 0.D0	! yx
		!---------------------------------------------------------------------------------------------------------------!
		! Calculate all interior values via the harmonic mean
		! The factor of 2 in the numerator is an artifact from using a 50/50 volume fraction, or more
		! simply, the fact that you weight each of 2 cells equally when calculating the face-centered value

		DO i=1,xmax
			DO j=1,ymax
				IF (i /= xmax)	THEN 	! all x terms
					!---------------------------------------------------------------------------------------------------!
					! Isotropic - x				
					IF ( (CellCentered%xx(i,j) /= 0.D0) .AND. (CellCentered%xx(i+1,j) /= 0.D0) ) THEN
						FaceCentered%xx_yy%Xface(i+1,j) = 2.D0 / ( 1.D0/CellCentered%xx(i,j) + 1.D0/CellCentered%xx(i+1,j) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xx_yy%Xface(i+1,j) = ( CellCentered%xx(i,j) + CellCentered%xx(i+1,j) ) / 2.D0
					ENDIF
					
					! Anisotropic - x
					IF ( (CellCentered%xy(i,j) /= 0.D0) .AND. (CellCentered%xy(i+1,j) /= 0.D0) ) THEN
						FaceCentered%xy_yx%Xface(i+1,j) = 2.D0 / ( 1.D0/CellCentered%xy(i,j) + 1.D0/CellCentered%xy(i+1,j) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xy_yx%Xface(i+1,j) = ( CellCentered%xy(i,j) + CellCentered%xy(i+1,j) ) / 2.D0
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDIF

				IF (j /= ymax) THEN		! all y terms
					!---------------------------------------------------------------------------------------------------!
					! Isotropic - y
					IF ( (CellCentered%yy(i,j) /= 0.D0) .AND. (CellCentered%yy(i,j+1) /= 0.D0) ) THEN
						FaceCentered%xx_yy%Yface(i,j+1) = 2.D0 / ( 1.D0/CellCentered%yy(i,j) + 1.D0/CellCentered%yy(i,j+1) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xx_yy%Yface(i,j+1) = ( CellCentered%yy(i,j) + CellCentered%yy(i,j+1) ) / 2.D0
					ENDIF
					
					! Anisotropic - y
					IF ( (CellCentered%xy(i,j) /= 0.D0) .AND. (CellCentered%xy(i,j+1) /= 0.D0) ) THEN
						FaceCentered%xy_yx%Yface(i,j+1) = 2.D0 / ( 1.D0/CellCentered%xy(i,j) + 1.D0/CellCentered%xy(i,j+1) )
					ELSE	! Avoids division by zero that skews the result
						FaceCentered%xy_yx%Yface(i,j+1) = ( CellCentered%xy(i,j) + CellCentered%xy(i,j+1) ) / 2.D0
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDIF
			ENDDO
		ENDDO

		! These are the harmonic mean done for the entire matrix at once. 
		! The downside is that any zeros will dominate the contribution for the 2 faces it's involved with
		! FaceCentered%Xface(2:xmax,:) = 2.D0 / ( 1.D0/CellCentered(1:xmax-1,:) + 1.D0/CellCentered(2:xmax,:) )
		! FaceCentered%Yface(:,2:ymax) = 2.D0 / ( 1.D0/CellCentered(:,1:ymax-1) + 1.D0/CellCentered(:,2:ymax) )
		!---------------------------------------------------------------------------------------------------------------!
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
		! Most likely, however, these 'boundary' values are not used since the boundary is treated differently.
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Cell_To_Face

	SUBROUTINE BoundaryConditions(phi,D_f)
	!-------------------------------------------------------------------------------------------------------------------!
    ! Subroutine calculates the boundary conditions for input scalar concentration
    ! As of now, only one type can be enforced on a side of the problem (left, right, top, bottom)
    ! BC: BC_alpha*phi(i,j) + d*BC_beta*( phi(i+1,j)-phi(i,j))/dx = BC_psi
    ! where d is the distance to the boundary. This is often set to the diffusion coefficient, 
    ! which has units of distance if there is a particle velocity in the differential equation
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT(INOUT)	:: phi		! Concentration function
        TYPE(Tensor_F), 					INTENT( IN)		:: D_f		! Face-Centered diffusion tensor (only needed for Robin BC)
		REAL(PRCSN), DIMENSION(xmax,ymax)					:: dumPhi
        REAL(PRCSN)											:: dummy
        INTEGER												:: BC_face,icell,jcell,istart,iend,jstart,jend
        REAL(PRCSN)											:: BC_Diric, BC_Neuma, BC_Value, BC_V
		!---------------------------------------------------------------------------------------------------------------!
		IF (BC_Fun == 1) CALL Phi_Value(dumPhi)		! Call function value if needed for boundary
		!---------------------------------------------------------------------------------------------------------------!
		! Iterate through all 4 faces, and then appropriate row/columns
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
			!-----------------------------------------------------------------------------------------------------------!
			! Assign appropriate BC for a particular face
			BC_Diric = BC_alpha(BC_face)
			BC_Neuma = BC_beta(BC_face)
			BC_Value = BC_psi(BC_face)
			!-----------------------------------------------------------------------------------------------------------!
			DO jcell=jstart,jend
				DO icell=istart,iend
				!-------------------------------------------------------------------------------------------------------!
					! Boundary values for when boundary is a function
					IF (BC_Fun == 1)	BC_Value = dumPhi(icell,jcell)
					!---------------------------------------------------------------------------------------------------!
					! Dirichlet
					IF ((BC_Diric /= 0.D0) .AND. (BC_Neuma == 0.D0))		THEN
						phi(icell,jcell) = BC_Value/BC_Diric
					!---------------------------------------------------------------------------------------------------!
					! Neumann
					ELSEIF ((BC_Diric == 0.D0) .AND. (BC_Neuma /= 0.D0))	THEN
						! 	The +/- dx/dy is determined by direction of the derivative
						IF (BC_face == 1) phi(icell,jcell) = phi(icell+1,jcell) - dx*BC_Value/BC_Neuma
						IF (BC_face == 2) phi(icell,jcell) = phi(icell-1,jcell) + dx*BC_Value/BC_Neuma
						IF (BC_face == 3) phi(icell,jcell) = phi(icell,jcell-1) + dy*BC_Value/BC_Neuma
						IF (BC_face == 4) phi(icell,jcell) = phi(icell,jcell+1) - dy*BC_Value/BC_Neuma
					!---------------------------------------------------------------------------------------------------!
					! Robin
					ELSE
						IF 		(BC_face == 1) 	THEN
							BC_V = -BC_Neuma*D_f%xx_yy%Xface(icell+1,jcell)/dx
							phi(icell,jcell) = (BC_Value+BC_V*phi(icell+1,jcell))/(BC_Diric+BC_V)
						ELSEIF (BC_face == 2) 	THEN
							BC_V =  BC_Neuma*D_f%xx_yy%Xface(icell,jcell)/dx
							phi(icell,jcell) = (BC_Value+BC_V*phi(icell-1,jcell))/(BC_Diric+BC_V)
						ELSEIF (BC_face == 3) 	THEN
							BC_V =  BC_Neuma*D_f%xx_yy%Yface(icell,jcell)/dy
							phi(icell,jcell) = (BC_Value+BC_V*phi(icell,jcell-1))/(BC_Diric+BC_V)
						ELSEIF (BC_face == 4) 	THEN
							BC_V = -BC_Neuma*D_f%xx_yy%Yface(icell,jcell+1)/dy
							phi(icell,jcell) = (BC_Value+BC_V*phi(icell,jcell+1))/(BC_Diric+BC_V)
						ENDIF
					ENDIF
					!---------------------------------------------------------------------------------------------------!
				ENDDO
			ENDDO
			!-----------------------------------------------------------------------------------------------------------!
		ENDDO
		!---------------------------------------------------------------------------------------------------------------!
		! Four corners 
		! 	Since BC are set by taking a full row or column, the corner values are not
		! 	well defined because there are (potentially) two values are trying to be assigned to them.
		! 	By taking the average of the nearest BC values, I hope to capture a more-true value.
!		phi(1,1)        = (phi(2,1)+phi(1,2))/2D0
!		phi(1,ymax)     = (phi(2,ymax)+phi(1,ymax-1))/2D0
!		phi(xmax,1)     = (phi(xmax-1,1)+phi(xmax,2))/2D0
!		phi(xmax,ymax)	= (phi(xmax,ymax-1)+phi(xmax-1,ymax))/2D0
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE BoundaryConditions
	
    SUBROUTINE ExplicitDiffusionDerivative(phi,D_f,RHS)
	!-------------------------------------------------------------------------------------------------------------------!
    ! Subroutine evaluates term DEL[ Kdiff DEL[u]], which has four components: d/dx d/dx, d/dy d/dy, d/dx d/dy, and d/dy d/dx
    ! DEL[ Kdiff DEL[u]] = DEL[(Kdiff%xx du/dx + Kdiff%xy du/dy)xhat+(Kdiff%xy du/dx + Kdiff%yy du/dy)yhat]
    !                    = d/dx[Kdiff%xx du/dx + Kdiff%xy du/dy] + d/dy[Kdiff%xy du/dx + Kdiff%yy du/dy]
    !                    = d/dx[Kdiff%xx du/dx] + d/dy[Kdiff%yy du/dy] + d/dx[Kdiff%xy du/dy] + d/dy[Kdiff%xy du/dx]
    !
    ! One major complication is that both the concentration values (u) and the diffusion tensor are located on the nodes
    ! However, after the first derivative, the values are now on half integer points.
    ! The result of this complication is that the diffusion tensors must be moved to the half integer points, also called the face-centers
    ! The subroutine FaceDiffusion accounts for any boundary effect complication in the constrution of K??_F_?
        REAL(PRCSN), DIMENSION(xmax,ymax),		INTENT( IN) :: phi
		TYPE(Tensor_F),							INTENT( IN) :: D_f
        REAL(PRCSN), DIMENSION(xmax,ymax),		INTENT(OUT) :: RHS
		!---------------------------------------------------------------------------------------------------------------!
		! Initialize variable to zero
		RHS = 0D0
        !---------------------------------------------------------------------------------------------------------------!
        DO j=2,ymax-1			! Start 1 column late and finish 1 column early
            DO i=2,xmax-1		! Start 1 row late and finish 1 row early
            ! Note that this hits every element of phi, but RHS is not defined on any boundary layer
            !-----------------------------------------------------------------------------------------------------------!
                ! d/dx (du_old/dx) -----------------  xx term ----------------------------------------------------------!
                RHS(i,j) = 	  RHS(i,j) +															&
                			( D_f%xx_yy%Xface(i+1,j  )*(phi(i+1,j  )-phi(i  ,j  )) 				&
                			- D_f%xx_yy%Xface(i  ,j  )*(phi(i  ,j  )-phi(i-1,j  )) )/(dx*dx)
                ! d/dy (du_old/dy) -----------------  yy term ----------------------------------------------------------!
                RHS(i,j) = 	  RHS(i,j) +															&
							( D_f%xx_yy%Yface(i  ,j+1)*(phi(i  ,j+1)-phi(i  ,j  )) 				&
							- D_f%xx_yy%Yface(i  ,j  )*(phi(i  ,j  )-phi(i  ,j-1)) )/(dy*dy)
                ! d/dx (du_old/dy) -----------------  xy term ----------------------------------------------------------!
                RHS(i,j) = 	RHS(i,j) +															&
							( D_f%xy_yx%Yface(i+1,j+1)*(phi(i+1,j+1)-phi(i+1,j  ))				&
							+ D_f%xy_yx%Yface(i+1,j  )*(phi(i+1,j  )-phi(i+1,j-1))				&
							- D_f%xy_yx%Yface(i-1,j+1)*(phi(i-1,j+1)-phi(i-1,j  ))				&
							- D_f%xy_yx%Yface(i-1,j  )*(phi(i-1,j  )-phi(i-1,j-1)) )/(4D0*dx*dy)
                ! d/dy (du_old/dx) -----------------  yx term ----------------------------------------------------------!
                RHS(i,j) = 	  RHS(i,j) +															&
							( D_f%xy_yx%Xface(i+1,j+1)*(phi(i+1,j+1)-phi(i  ,j+1))				&
							+ D_f%xy_yx%Xface(i  ,j+1)*(phi(i  ,j+1)-phi(i-1,j+1))				&
							- D_f%xy_yx%Xface(i+1,j-1)*(phi(i+1,j-1)-phi(i  ,j-1))				&
							- D_f%xy_yx%Xface(i  ,j-1)*(phi(i  ,j-1)-phi(i-1,j-1)) )/(4D0*dx*dy)
            !-----------------------------------------------------------------------------------------------------------!
            ENDDO
        ENDDO
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE ExplicitDiffusionDerivative

	SUBROUTINE Rotation(Kxx,Kyy,Kxy,theta)
	!-------------------------------------------------------------------------------------------------------------------!
	! Rotates the diffusion tensor, but only the elements passed to here
	! Take K(i,j) = | Kxx(i,j) & Kxy(i,j) |
	!               | Kyx(i,j) & Kyy(i,j) |
	! Multiply by R = | cos & -sin |
	!                 | sin &  cos |
	! Such that Keffective = R*K*R', where R' is the transpose
	! Since we're dealing with a single element of K, it's only a double precision number, not an array
        REAL(PRCSN), INTENT(INOUT)	:: Kxx, Kyy, Kxy
	    REAL(PRCSN), INTENT( IN)	:: theta            ! theta is in radians	    
        REAL(PRCSN)					:: RKxx, RKyy, RKxy
		
		! Rotated values ---  R*D*R' where R is rotation matrix, R' is transpose
        RKxx = Kxx*cos(theta)**2+Kyy*sin(theta)**2-Kxy*sin(2D0*theta)
		RKyy = Kxx*sin(theta)**2+Kyy*cos(theta)**2+Kxy*sin(2D0*theta)
		RKxy = Kxy*cos(2D0*theta)+0.5d0*(Kxx-Kyy)*sin(2D0*theta)
		!RKyx = RKxy ! Not necessary, but may be valuable at some point
        ! Now that rotation matrix is computed, pass the rotate values back out		
		Kxx = RKxx
		Kyy = RKyy
		Kxy = RKxy
	!-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE Rotation
	
	SUBROUTINE GD_FaceDiffusion(D_c,D_f,phi)
	!-------------------------------------------------------------------------------------------------------------------!
	! Gradient Dependent Face Diffusion subroutine
	! Subroutine outputs face-centered diffusion tensor from an input of cell-centered tensors.
	! This subroutine tries to take into acount concentration gradients as well as the value of the 
	! diffusion tensor between 2 points. Different behavior (arithmetic or harmonic mean) is executed
	! based on the concentration and diffusion-coefficient values.
	!
	! Since this takes the harmonic mean of the two cells neighboring a face, it depends on if you
	! considering an x-direction or y-direction face-centered value
	!
	! Complications arise when you take the harmonic mean of two numbers, when one or both is/are equal to zero
	! The inverse sum causes an infinite term in the denominator, which changes the harmonic mean to 0
	! This is avoided by simply setting the face value to the arithmetic mean
		REAL(PRCSN), DIMENSION(xmax,ymax),	INTENT (IN) :: phi
        TYPE(Tensor_C),                  	INTENT (IN) :: D_c
		TYPE(Tensor_F), 					INTENT(OUT) :: D_f
		INTEGER											:: MethodFlag
		
		!MethodFlag = 1	! Face-centered diffusion is defined from a single cell-centered value
		MethodFlag = 2	! Face-centered diffusion is defined as avg/harmonic depending on D_c and phi
		
		D_f%xx_yy%Xface=0D0; D_f%xx_yy%Yface=0D0; D_f%xy_yx%Xface=0D0; D_f%xy_yx%Yface=0D0

	    DO i=1,xmax-1
            DO j=1,ymax-1
	            IF (i>1)                        						THEN
					IF (phi(i,j) > phi(i-1,j)) 						THEN
						IF (MethodFlag == 2)							THEN
							IF (D_c%xx(i,j) >= D_c%xx(i-1,j))			THEN
								D_f%xx_yy%Xface(i,j) = 2D0/(1D0/D_c%xx(i-1,j) + 1D0/D_c%xx(i,j))							
							ELSE
								D_f%xx_yy%Xface(i,j) = (D_c%xx(i,j) + D_c%xx(i-1,j))/2D0
							ENDIF
	
							IF (D_c%xy(i,j) >= D_c%xy(i-1,j))			THEN
								D_f%xy_yx%Xface(i,j) = 2D0/(1D0/D_c%xy(i-1,j) + 1D0/D_c%xy(i,j))
							ELSE
								D_f%xy_yx%Xface(i,j) = (D_c%xy(i,j) + D_c%xy(i-1,j))/2D0
							ENDIF
						ELSEIF (MethodFlag == 1) 						THEN
							D_f%xx_yy%Xface(i,j) = D_c%xx(i-1,j)
							D_f%xy_yx%Xface(i,j) = D_c%xy(i-1,j)
						ENDIF
					ELSE	! phi(i,j) <= phi(i-1,j)
						IF (MethodFlag == 2)							THEN
							IF (D_c%xx(i,j) >= D_c%xx(i-1,j))			THEN
								D_f%xx_yy%Xface(i,j) = (D_c%xx(i,j) + D_c%xx(i-1,j))/2D0
							ELSE
								D_f%xx_yy%Xface(i,j) = 2D0/(1D0/D_c%xx(i-1,j) + 1D0/D_c%xx(i,j))							
							ENDIF					
	
							IF (D_c%xy(i,j) >= D_c%xy(i-1,j))			THEN
								D_f%xy_yx%Xface(i,j) = (D_c%xy(i,j) + D_c%xy(i-1,j))/2D0
							ELSE
								D_f%xy_yx%Xface(i,j) = 2D0/(1D0/D_c%xy(i-1,j) + 1D0/D_c%xy(i,j))
							ENDIF
						ELSEIF (MethodFlag == 1) 	THEN
							D_f%xx_yy%Xface(i,j) = D_c%xx(i,j)
							D_f%xy_yx%Xface(i,j) = D_c%xy(i,j)
						ENDIF
					ENDIF
                ENDIF

	            IF (j>1)                        						THEN
					IF (phi(i,j) > phi(i,j-1))						THEN
						IF (MethodFlag == 2)							THEN
							IF (D_c%yy(i,j) >= D_c%yy(i,j-1))			THEN
								D_f%xx_yy%Yface(i,j) = 2D0/(1D0/D_c%yy(i,j-1) + 1D0/D_c%yy(i,j))
							ELSE
								D_f%xx_yy%Yface(i,j) = (D_c%yy(i,j) + D_c%yy(i,j-1))/2D0
							ENDIF
	
							IF (D_c%xy(i,j) >= D_c%xy(i,j-1))			THEN
								D_f%xy_yx%Yface(i,j) = 2D0/(1D0/D_c%xy(i,j-1) + 1D0/D_c%xy(i,j))
							ELSE
								D_f%xy_yx%Yface(i,j) = (D_c%xy(i,j) + D_c%xy(i,j-1))/2D0
							ENDIF
						ELSEIF (MethodFlag == 1) 	THEN
							D_f%xx_yy%Yface(i,j) = D_c%yy(i,j-1)
							D_f%xy_yx%Yface(i,j) = D_c%xy(i,j-1)
						ENDIF
					ELSE	! (phi(i,j) <= phi(i,j-1))	
						IF (MethodFlag == 2)							THEN
							IF (D_c%yy(i,j) >= D_c%yy(i,j-1))			THEN
								D_f%xx_yy%Yface(i,j) = (D_c%yy(i,j) + D_c%yy(i,j-1))/2D0
							ELSE
								D_f%xx_yy%Yface(i,j) = 2D0/(1D0/D_c%yy(i,j-1) + 1D0/D_c%yy(i,j))
							ENDIF
	
							IF (D_c%xy(i,j) >= D_c%xy(i,j-1))			THEN
								D_f%xy_yx%Yface(i,j) = (D_c%xy(i,j) + D_c%xy(i,j-1))/2D0
							ELSE
								D_f%xy_yx%Yface(i,j) = 2D0/(1D0/D_c%xy(i,j-1) + 1D0/D_c%xy(i,j))
							ENDIF
						ELSEIF (MethodFlag == 1) 	THEN
							D_f%xx_yy%Yface(i,j) = D_c%yy(i,j)
							D_f%xy_yx%Yface(i,j) = D_c%xy(i,j)
						ENDIF
					ENDIF
                ENDIF
!IF (D_c%xy(i,j)/=0) PRINT '(2(I2,X),3(F8.4,X))', i,j,D_c%xy(i,j),D_f%xy_yx%Xface(i,j),D_f%xy_yx%Yface(i,j)
            ENDDO
        ENDDO
        
        !---------------------------------------------------------------------------------------------------------------!
        ! Define supplementary edge values of the x direction face-centered diffusion tensors
        ! This allows the explicit derivative to be more simplistic, since complications are accounted for here
        !---------------------------------------------------------------------------------------------------------------!
        ! Make the bottom interior row equal to just the cell-centered value (not the harmonic mean)
        ! This is because the finite difference is between the bottom interior row with the bottom boundary,
        ! and since there is no diffusion tensor on the boundary, the cell-centered value is the only diffusion term involved
        D_f%xx_yy%Xface(1,:) = D_c%xx(1,:)
        D_f%xy_yx%Xface(1,:) = D_c%xy(1,:)
        D_f%xy_yx%Yface(:,1) = D_c%xy(:,1)
        D_f%xx_yy%Yface(:,1) = D_c%yy(:,1)        
        ! Make any refrences to the top row (boundary) equal to the top interior cell-centered value (not the harmonic mean)
        ! This is because the finite difference is between the top boundary and the top interior row
        ! and since there is no diffusion tensor on the boundary, the cell-centered value is the only diffusion term involved
        D_f%xx_yy%Xface(xmax+1,:) =	D_c%xx(xmax,:)
        D_f%xy_yx%Xface(xmax+1,:) =	D_c%xy(xmax,:)
        D_f%xy_yx%Yface(:,ymax+1) =	D_c%xy(:,ymax)
        D_f%xx_yy%Yface(:,ymax+1) =	D_c%yy(:,ymax)        
        ! Make boundary columns duplicates of interior columns since diffusion tensor is not defined on the boundaries
        D_f%xy_yx%Xface(:,1)    = D_f%xy_yx%Xface(:,2)
        D_f%xx_yy%Xface(:,1)    = D_f%xx_yy%Xface(:,2)
        D_f%xy_yx%Xface(:,ymax) = D_f%xy_yx%Xface(:,ymax-1)
        D_f%xx_yy%Xface(:,ymax) = D_f%xx_yy%Xface(:,ymax-1)
        D_f%xy_yx%Yface(1,:) 	= D_f%xy_yx%Yface(2,:)
        D_f%xx_yy%Yface(1,:)    = D_f%xx_yy%Yface(2,:)
        D_f%xy_yx%Yface(xmax,:) = D_f%xy_yx%Yface(xmax-1,:)
        D_f%xx_yy%Yface(xmax,:) = D_f%xx_yy%Yface(xmax-1,:)        
        !---------------------------------------------------------------------------------------------------------------!
        
        !CALL SaveGeneral(xmax-1,ymax-1, 8,'C//Dxx.m',  D_c%xx)
        !CALL SaveGeneral(xmax-1,ymax-1, 8,'C//Dxy.m',  D_c%xy)
        !CALL SaveGeneral(xmax-1,ymax-1, 8,'C//Dyy.m',  D_c%yy)
        !CALL SaveGeneral(xmax  ,ymax+1,10,'C//DxxFx.m',D_f%xx_yy%Xface)
        !CALL SaveGeneral(xmax  ,ymax+1,10,'C//DxyFx.m',D_f%xy_yx%Xface)
        !CALL SaveGeneral(xmax+1,ymax  ,10,'C//DxyFy.m',D_f%xy_yx%Yface)
        !CALL SaveGeneral(xmax+1,ymax  ,10,'C//DyyFy.m',D_f%xx_yy%Yface)
    !-------------------------------------------------------------------------------------------------------------------!
	END SUBROUTINE GD_FaceDiffusion
END MODULE SubroutineModules
