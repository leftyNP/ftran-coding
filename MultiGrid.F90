!--------------------------------------------------------------------------------------------------!
MODULE COMMON
	IMPLICIT NONE
	!PRIVATE
	INTEGER 						:: i,j,iO
	INTEGER,				PUBLIC	:: xpoints, ypoints
	REAL*8, 				PUBLIC	:: hx, hy, dt,tol		! h is grid spacing, dt is time step
	REAL*8, DIMENSION(2), 	PUBLIC 	:: sys_length
	CHARACTER(LEN=6), 		PUBLIC	:: BCflag 
	CONTAINS
	!So I guess you are getting ~10x (keeping the accuracy the same) or more (if you let the accuracy go down a bit). Ê
	!There are also exponent-based methods that might give a better accuracy (but might be slower). Ê
	!Can you make a note, and we can discuss when you are back to this later.
	
	SUBROUTINE init(E)
	!----------------------------------------------------------------------------------------------!
	! Sets the initial conditions of the system (in this case, specifies energy everywhere)
	REAL*8, DIMENSION(xpoints,ypoints), INTENT(out)	:: E
	REAL*8 											:: x, y

	DO i=1,xpoints
		DO j=1,ypoints
			x = DBLE(i)/DBLE(xpoints)*sys_length(1)-sys_length(1)*.5d0
		 	y = DBLE(j)/DBLE(ypoints)*sys_length(2)-sys_length(2)*.5d0			
!			E(i,j) = exp( -(x-1.65d0)**2/4.d-1 )
			E(i,j) = exp( -(x**2+y**2)/4.d-1 )
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE init

	SUBROUTINE discretization(size,A)	
	!----------------------------------------------------------------------------------------------!
	! This forms the matrix A which comes from the discretization of the equation
	! F(E^{n+1}) = (I-dt*del^2) dE
	! Specify boundary conditions based on flag.  
	! flag may be the following:  BCflag ==
	!										NoFlux	- del*E = 0 on boundary
	!										Period	- periodic boundary conditions (like Pacman game)
	!										Dirich	- E = 0 on boundary
	INTEGER, DIMENSION(2), 								INTENT(in)	:: size
	REAL*8, DIMENSION(size(1)*size(2),size(1)*size(2)), INTENT(out)	:: A
	REAL*8, DIMENSION(size(1)*size(2),size(1)*size(2))		   		:: L,ID	! L-->Discretized Laplacian, ID-->identity matrix

	hx = sys_length(1)/size(1) ! grid spacing: 1/xpoints assumes a unit square grid
	hy = sys_length(2)/size(2)
	!	dt = 1e-5  	!2*(0.5*(hx+hy))**2 ! time step size. this should change depending on other variables

	! NOTE: is it necessary to specify 0's? or will Fortran automatically put non-specified entries as 0?
	L  = 0.d0				! sets all entries to 0 before specifying non-zero entries
	ID = 0.d0				! sets all entries to 0 before specifying non-zero entries

	DO i=1,size(1)
		DO j=1,size(2)
			iO = (j-1)*size(1)+i	! global position of each i,j component such that a matrix can be represented as a vector
							L(iO,iO		   ) =-2.d0*(1/hx**2+1/hy**2)	! center	! when hx=hy, this becomes -4/h**2
			IF (i/=1	  ) L(iO,iO-1      ) = 1.d0/hx**2				! down
			IF (i/=size(1)) L(iO,iO+1      ) = 1.d0/hx**2				! up
			IF (j/=1      ) L(iO,iO-size(1)) = 1.d0/hy**2				! left
			IF (j/=size(2)) L(iO,iO+size(1)) = 1.d0/hy**2				! right
						   ID(iO,iO		   ) = 1.d0						!  
		ENDDO
	ENDDO

	CALL BC(size,L)

	A = ID - 0.5d0*dt*L	! A is now the the matrix in r=Ax where r is the residual to be minimized, and x is the solution we seek
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE discretization

	SUBROUTINE BC(size,L)
	!----------------------------------------------------------------------------------------------!
	! Specify boundary conditions based on flag.  
	! flag may be the following:  BCflag ==
	!										NoFlux	- del*E = 0 on boundary
	!										Period	- periodic boundary conditions (like Pacman)
	!										Dirich	- E = 0 on boundary
	INTEGER,DIMENSION(2),								INTENT(in)		:: size
	REAL*8,DIMENSION(size(1)*size(2),size(1)*size(2)), 	INTENT(inout)	:: L

	hx = sys_length(1)/size(1) ! grid spacing: 1/xpoints assumes a unit square grid
	hy = sys_length(2)/size(2)
	L  = L 						! this is only necessary so L is assigned something in the case of Dirichlet BCs

	DO i=1,size(1)
		DO j=1,size(2)
			iO = (j-1)*size(1)+i
			IF (BCflag=='NoFlux')		THEN
			! Boundary Conditions: No Flux  dE/dx = 0 at i=1,size(1)   & dE/dy = 0 at j=1,size(2)
				IF (i==1	  )						L(iO,iO) =-1.d0/hx**2-2.d0/hy**2
				IF (i==size(1))						L(iO,iO) =-1.d0/hx**2-2.d0/hy**2
				IF (j==1	  )						L(iO,iO) =-2.d0/hx**2-1.d0/hy**2
				IF (j==size(2)) 					L(iO,iO) =-2.d0/hx**2-1.d0/hy**2
				IF (i==1       .AND. j==1      ) 	L(iO,iO) =-1.d0/hx**2-1.d0/hy**2
				IF (i==1	   .AND. j==size(2))	L(iO,iO) =-1.d0/hx**2-1.d0/hy**2
				IF (i==size(1) .AND. j==1      )	L(iO,iO) =-1.d0/hx**2-1.d0/hy**2
				IF (i==size(1) .AND. j==size(2))	L(iO,iO) =-1.d0/hx**2-1.d0/hy**2

			ELSEIF (BCflag=='Period') THEN
			! Boundary Conditions: Periodic		d^2 E(i,j)/dx^2 = (E(i+1,j)+E(i-1,j)-2E(i,j) / hx**2
			! For endpoints, the +/- term gets assigned to opposite boundary. Ex:  i==1: E(i-1,j)=E(0,j) ? no, set it to E(size(1),j)
				IF (i==1	  )	L(iO,       j   *size(1)  ) = 1.d0/hx**2		
				IF (i==size(1)) L(iO,(      j-1)*size(1)+1) = 1.d0/hx**2		
				IF (j==1	  )	L(iO,(size(2)-1)*size(1)+i) = 1.d0/hy**2		
				IF (j==size(2)) L(iO,   				 i) = 1.d0/hy**2		
				
			ELSEIF (BCflag=='Dirich') THEN
			! Boundary Conditions: Dirichlet	E(i,j) = 0 for any boundary:  i=:{1,size(1)}, j=:{1,size(2)}
			! Instead, set boundary to i==0 and i==size(1)+1 and j==0 and j==size(2)+1
			! This is identical to the current set up, except we do not compute the rows and columns of zeros bounding the non-zero matrix
				! No work is necessary	
			ELSE
				WRITE(*,*) 'Boundary Conditions not specified - Assuming Dirichlet'
			ENDIF
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE BC

	SUBROUTINE calc_res(length,E,Enew,res)
	!----------------------------------------------------------------------------------------------!
	! calculates the residual for the current iteration of E
	INTEGER, DIMENSION(2), 					INTENT(in)	:: length		! length of matrix side
	REAL*8, DIMENSION(length(1),length(2)), INTENT(in)	:: E			! Matrix found on the previous time step
	REAL*8, DIMENSION(length(1),length(2)), INTENT(in)	:: Enew			! Matrix being currently calculated
	REAL*8, DIMENSION(length(1),length(2)), INTENT(out)	:: res			! Residual in vector form 
	REAL*8, DIMENSION(length(1),length(2))				:: lapE			! Laplacian term of matrix inputed in laplace call
	REAL*8, DIMENSION(length(1),length(2))				:: square_res	! Residual in square matrix form


	CALL laplace(length,0.5d0*(Enew+E),lapE)		! when calling:
										!				input Enew 			for implicit
										!				input E 			for explicit
										!				input 0.5*(E+Enew) 	for midpoint implicit
	res = -((Enew-E) - dt*lapE)
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE calc_res

	SUBROUTINE laplace(size,E,lapE)
	!----------------------------------------------------------------------------------------------!
	! this calculates the laplacian of a given matrix
	! Specify boundary conditions for laplacian based on flag.  
	! flag may be the following:  BCflag ==
	!										NoFlux	- del*E = 0 on boundary
	!										Period	- periodic boundary conditions (like Pacman game)
	!										Dirich	- E = 0 on boundary
	INTEGER, DIMENSION(2), 				INTENT(in)	:: size 	! length of matrix (one side)
	REAL*8, DIMENSION(size(1),size(2)), INTENT(in)	:: E 		! input matrix
	REAL*8, DIMENSION(size(1),size(2)), INTENT(out)	:: lapE		! 2nd derivative of input matrix
	REAL*8											:: xs,ys 	! these are x-spacing, yspacing

	xs = 1.d0/hx**2
	ys = 1.d0/hy**2

	DO i=1,size(1)
		DO j=1,size(2)
		!------------------------------------------------------------------------------------------!
		! All Interior Points
			IF (i/=1 .OR. i/=size(1) .OR. j/=1 .OR. j/=size(2))	 THEN
				lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
			ENDIF
		!------------------------------------------------------------------------------------------!
			IF (BCflag=='NoFlux')		THEN
			! Boundary Conditions: No Flux  dE/dx = 0 at i=1,size(1)   & dE/dy = 0 at j=1,size(2)
			!--------------------------------------------------------------------------------------!
				! All Edge Points not on a corner  - Ex: (i==1) - since E(0,j)=E(1,j), one of the -4E(i,j) gets cancelled by the E(i-1,j) term				
				! Bottom Row ----------------------------------------------------------------------!
				IF (i==size(1) .AND. (j/=1 .OR. j/=size(2)))	&	
					lapE(i,j)=(  0.0d0 +E(i-1,j)-1.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Top Row -------------------------------------------------------------------------!
				IF (i==1 .AND. (j/=1 .OR. j/=size(2)))			&
					lapE(i,j)=(E(i+1,j)+  0.0d0 -1.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys			
				! Right Column --------------------------------------------------------------------!
				IF (j==size(2) .AND. (i/=1 .OR. i/=size(1)))	&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-1.d0*E(i,j))*ys
				! Left Column ---------------------------------------------------------------------!
				IF (j==1 .AND. (i/=1 .OR. i/=size(1)))			&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -1.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			! The four corner points  - similiar to above edge points, but here 2 factors of -4E(i,j) are cancelled 
				! Top Left ------------------------------------------------------------------------!
				IF (i==1 .AND. j==1)							&
					lapE(i,j)=(E(i+1,j)+  0.0d0 -1.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -1.d0*E(i,j))*ys
				! Top Right -----------------------------------------------------------------------!
				IF (i==1 .AND. j==size(2))						&
					lapE(i,j)=(E(i+1,j)+  0.0d0 -1.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-1.d0*E(i,j))*ys
				! Bottom Left ---------------------------------------------------------------------!
				IF (i==size(1) .AND. j==1) 						&
					lapE(i,j)=(  0.0d0 +E(i-1,j)-1.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -1.d0*E(i,j))*ys
				! Bottom Right --------------------------------------------------------------------!
				IF (i==size(1) .AND. j==size(2))				&
					lapE(i,j)=(  0.0d0 +E(i-1,j)-1.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-1.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			ELSEIF (BCflag=='Period') THEN
			! Boundary Conditions: Periodic		d^2 E(i,j)/dx^2 = (E(i+1,j)+E(i-1,j)-2E(i,j))/ hx**2
			! For endpoints, the +/- term gets assigned to opposite boundary. Ex:  i==1: E(i-1,j)=E(0,j) ? no, set it to E(size(1),j)				
			!--------------------------------------------------------------------------------------!
			! All Edge Points not on a corner  - E(0,j)=E(size(1),j), E(size(1)+1)=E(1,j), E(i,size(2)+1)=E(i,1), E(i,0)=E(i,size(2))
				! Bottom Row ----------------------------------------------------------------------!
				IF (i==size(1) .AND. (j/=1 .OR. j/=size(2)))	&
					lapE(i,j)=(E(1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Top Row -------------------------------------------------------------------------!
				IF (i==1 .AND. (j/=1.OR.j/=size(2)))			&
					lapE(i,j)=(E(i+1,j)+E(size(1),j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Right Column --------------------------------------------------------------------!
				IF (j==size(2) .AND. (i/=1 .OR. i/=size(1))) 	&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Left Column ---------------------------------------------------------------------!
				IF (j==1 .AND. (i/=1 .OR. i/=size(1)))			&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,size(2))-2.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			! The four corner points  - identical to above section
				! Top Left ------------------------------------------------------------------------!
				IF (i==1 .AND. j==1)							&
					lapE(i,j)=(E(i+1,j)+E(size(1),j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,size(2))-2.d0*E(i,j))*ys
				! Top Right -----------------------------------------------------------------------!
				IF (i==1 .AND. j==size(2))						&
					lapE(i,j)=(E(i+1,j)+E(size(1),j)-2.d0*E(i,j))*xs+(E(i,1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Bottom Left ---------------------------------------------------------------------!
				IF (i==size(1) .AND. j==1)						&
					lapE(i,j)=(E(1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,size(2))-2.d0*E(i,j))*ys
				! Bottom Right --------------------------------------------------------------------!
				IF (i==size(1) .AND. j==size(2))				&
					lapE(i,j)=(E(1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,1)+E(i,j-1)-2.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			ELSEIF (BCflag=='Dirich') THEN
			! Boundary Conditions: Dirichlet	E(i,j) = 0 for any boundary:  i=:{1,size(1)}, j=:{1,size(2)}
			! Instead, set boundary to i==0 and i==size(1)+1 and j==0 and j==size(2)+1
			! So, the spaces below with a zero represents a boundary point. 
			! If we had E(i,j)=f(i,j) on a boundary, we would have f(i,j) in place of the zeros
			!--------------------------------------------------------------------------------------!
			! All Edge Points not on a corner
				! Top Row -------------------------------------------------------------------------!
				IF (i==1 .AND. (j/=1 .OR. j/=size(2)))			&
					lapE(i,j)=(E(i+1,j)+  0.0d0 -2.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Bottom Row ----------------------------------------------------------------------!
				IF (i==size(1) .AND. (j/=1 .OR. j/=size(2))) 	&	
					lapE(i,j)=(  0.0d0 +E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+E(i,j-1)-2.d0*E(i,j))*ys
				! Left Column ---------------------------------------------------------------------!
				IF (j==1 .AND. (i/=1 .OR. i/=size(1)))			&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -2.d0*E(i,j))*ys
				! Right Column --------------------------------------------------------------------!
				IF (j==size(2) .AND. (i/=1 .OR. i/=size(1)))	&
					lapE(i,j)=(E(i+1,j)+E(i-1,j)-2.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-2.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			! The four corner points
				! Top Left ------------------------------------------------------------------------!
				IF (i==1 .AND. j==1)							&	
					lapE(i,j)=(E(i+1,j)+  0.0d0 -2.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -2.d0*E(i,j))*ys
				! Top Right -----------------------------------------------------------------------!
				IF (i==1 .AND. j==size(2))						&	
					lapE(i,j)=(E(i+1,j)+  0.0d0 -2.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-2.d0*E(i,j))*ys
				! Bottom Left ---------------------------------------------------------------------!
				IF (i==size(1) .AND. j==1) 						&
					lapE(i,j)=(  0.0d0 +E(i-1,j)-2.d0*E(i,j))*xs+(E(i,j+1)+  0.0d0 -2.d0*E(i,j))*ys
				! Bottom Right --------------------------------------------------------------------!
				IF (i==size(1) .AND. j==size(2)) 				&
					lapE(i,j)=(  0.0d0 +E(i-1,j)-2.d0*E(i,j))*xs+(  0.0d0 +E(i,j-1)-2.d0*E(i,j))*ys
			!--------------------------------------------------------------------------------------!
			ENDIF ! ends BCflag if statement
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE laplace
	
	SUBROUTINE Mat_Vec(size,A,v)
	!----------------------------------------------------------------------------------------------!
	! this converts a 2D array (A) to a 1D vector (v).  i.e. converts a MATrix to a VECtor
	INTEGER, DIMENSION(2), 				INTENT(in)	:: size
	REAL*8, DIMENSION(size(1),size(2)), INTENT(in) 	:: A
	REAL*8, DIMENSION(size(1)*size(2)), INTENT(out)	:: v
	
	DO j=1,size(2)
		DO i=1,size(1)
			iO    = (j-1)*size(1)+i	! global position
			v(iO) = A(i,j)
		ENDDO
	ENDDO	
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE Mat_Vec

	SUBROUTINE SGS(length,A,sqx,b,tol,guess)
	!----------------------------------------------------------------------------------------------!
	! Symmetric Gauss-Seidel (SGS) algorithm to iteratively solve the system Ax=b
	! SGS (distinct from GS) because the direction of iteration alternates forwards/backwards
	INTEGER, DIMENSION(2), 										INTENT(in)	:: length
	REAL*8, 													INTENT(in)	:: tol		! interation tolerance for error
	REAL*8, DIMENSION(length(1)*length(2),length(1)*length(2)), INTENT(in)	:: A		! matrix of system
	REAL*8, DIMENSION(length(1)*length(2)), 					INTENT(in)	:: b		! RHS of system
	REAL*8, DIMENSION(length(1)*length(2)), OPTIONAL, 			INTENT(in)	:: guess 	! initial guess of solution
	REAL*8, DIMENSION(length(1),length(2)), 					INTENT(out)	:: sqx		! reshaped solution of system (square x)
	REAL*8, DIMENSION(length(1)*length(2))						 			:: y,res	! res=residual,y is previous guess for x
	REAL*8, DIMENSION(length(1)*length(2))					 				:: x		! solution of system (to be found)
	REAL*8																	:: sum,mag	! mag is magnitude of residual
			
    IF (present(guess)) THEN
        y = guess
    ELSE
        FORALL (i=1:length(1)*length(2)) y(i)=b(i)/a(i,i)
    ENDIF
    
	! x is current best guess for solution   x --> x(k+1)
	! y is previous best guess for solution  y --> x(k)
    x=y						! sets initial vector to y
    mag=1.d0+tol				! sets magnitude of residual larger than tolerance to enter while-loop
    
    DO WHILE (mag>tol)
    	y=x					! saves non-updated vector x to calculate residual for convergence test
		!------------------------------------------------------------------------------------------!
    	! forward loop
    	DO i=1,length(1)*length(2)
    		sum=0.d0
     		DO j=1,i-1
     			sum = sum + A(i,j)*x(j)
     		ENDDO
			!--------------------------------------------------------------------------------------!
    		DO j=i+1,length(1)*length(2)
    			sum = sum + A(i,j)*x(j)
    		ENDDO
    		x(i)=(b(i)-sum)/A(i,i)
    	ENDDO
		!------------------------------------------------------------------------------------------!
		! backward loop
    	DO i=length(1)*length(2),1,-1
    		sum=0.d0
    		DO j=i-1,1,-1
    			sum = sum + A(i,j)*x(j)
    		ENDDO
			!--------------------------------------------------------------------------------------!    		
    		DO j=length(1)*length(2),i+1,-1
    			sum = sum + A(i,j)*x(j)
    		ENDDO
    		x(i)=(b(i)-sum)/A(i,i)
    	ENDDO
    	res=x-y
    	CALL norm(length,res,mag)	! calculated magnitude of residual (Euclidian norm = mag(res) )
    ENDDO
    sqx = reshape(x,(/length(1),length(2)/))
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE SGS

	subroutine smoother(size,res,dE)  ! note: this is smoother, NOT smother
	!----------------------------------------------------------------------------------------------!
	! smoothes acts as a solver for Ax=b, where A is called, b is an input, and x (dE) is solved for
	! this effectively smoothes the error, and may be used as a post/pre smoother, or a coarse grid solver
	!
	! Since 4 subroutines need to be called in a specific order, several times, this subroutine
	! simply makes all the calls in proper order
	INTEGER, DIMENSION(2), 				INTENT(in)		:: size
	REAL*8,DIMENSION(size(1),size(2)), 	INTENT(in)		:: res
	REAL*8,DIMENSION(size(1),size(2)),	INTENT(inout)	:: dE
	REAL*8,DIMENSION(size(1)*size(2))					:: b,guess
	REAL*8,DIMENSION(size(1)*size(2),size(1)*size(2))	:: A

	
	CALL Mat_Vec(size,res,b)				! convert residual to a vector, creating RHS
	CALL Mat_Vec(size,dE,guess)				! create a guess for SGS solver
	CALL discretization(size,A)				! create the matrix on LHS of Ax=b
	CALL SGS(size,A,dE,b,tol,guess)			! solve system Ax=b for x, where x is dE
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE smoother

	SUBROUTINE norm(length,vector,mag)
	!----------------------------------------------------------------------------------------------!
	! this calculates the L2 norm (Euclidian length) of a vector. this does not handle matricies
	INTEGER, DIMENSION(2), 					INTENT(in)	:: length
	REAL*8, DIMENSION(length(1)*length(2)), INTENT(in)	:: vector
	REAL*8, 								INTENT(out)	:: mag

	mag = 0.d0
	DO i=1,length(1)*length(2)
		mag = mag + vector(i)*vector(i)		! Multiplication is faster than flexible exponent
	ENDDO
	mag = sqrt(mag)
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE norm

	SUBROUTINE URfact(size,dE,E)
	!----------------------------------------------------------------------------------------------!
	! calculates under-relaxation factor "to robustly deal with convergence difficulties"
	INTEGER, DIMENSION(2), 				 INTENT(in)		:: size
	REAL*8,  DIMENSION(size(1),size(2)), INTENT(in)		:: dE
	REAL*8,  DIMENSION(size(1),size(2)), INTENT(in out)	:: E
	REAL*8,  DIMENSION(size(1)*size(2))					:: dvec,vec,rvec
	REAL*8												:: mag, xi

	CALL Mat_Vec(size,dE,dvec)
	CALL Mat_Vec(size,E,vec)
	rvec = dvec/vec
	CALL norm(size,rvec,mag)
	xi = min(1.d0,1.d0/mag)

	IF (xi>1.d-2) 	THEN
		E = E + xi*dE					! calculate under-relaxed Enew
	ELSE
		E = E+dE						! calculate Enew
	ENDIF
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE URfact

	SUBROUTINE magnitude(size,Eold,Enew,mag)
	!----------------------------------------------------------------------------------------------!
	! calculates the magnitude of the residual
	INTEGER, DIMENSION(2), 				INTENT(in)	:: size
	REAL*8, DIMENSION(size(1),size(2)), INTENT(in)	:: Eold, Enew
	REAL*8, 							INTENT(out)	:: mag
	REAL*8, DIMENSION(size(1),size(2))				:: res
	REAL*8, DIMENSION(size(1)*size(2))				:: b
		
	CALL calc_res(size,Eold,Enew,res)
	CALL Mat_Vec(size,res,b)
	CALL norm(size,b,mag)
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE magnitude

	SUBROUTINE restriction(size,A,B)
	!----------------------------------------------------------------------------------------------!
	! coarsens input matrix (A) into matrix (B)
	INTEGER, DIMENSION(2), 							INTENT(in)	:: size
	REAL*8, DIMENSION(size(1),size(2)), 			INTENT(in)	:: A
	REAL*8, DIMENSION((size(1)+1)/2,(size(2)+1)/2),	INTENT(out)	:: B
	INTEGER														:: c,d

	DO i=1,size(1),2
		DO j=1,size(2),2
			c=(i+1)/2
			d=(j+1)/2
		!------------------------------------------------------------------------------------------!
		! All Interior Points
			IF (i/=1 .OR. i/=size(1) .OR. j/=1 .OR. j/=size(2)) 								&
				B(c,d)=((A(i-1,j-1)+A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))							&
						+2.d0*(A(i,j-1)+A(i,j+1)+A(i+1,j)+A(i-1,j)) +4.d0*A(i,j))/16.d0
		!------------------------------------------------------------------------------------------!
		! All Edge Points not on a corner
			! Top Row -----------------------------------------------------------------------------!
			IF (i==1 .AND. (j/=1 .OR. j/=size(2)))  											&
				B(c,d)=((0.d0+0.d0+A(i+1,j-1)+A(i+1,j+1)) 										&
					+2.d0*(A(i,j-1)+A(i,j+1)+A(i+1,j)+0.d0)+4.d0*A(i,j))/12.d0		
			! Bottom Row --------------------------------------------------------------------------!
			IF (i==size(1) .AND. (j/=1 .OR. j/=size(2)))										&
				B(c,d)=((A(i-1,j-1)+A(i-1,j+1)+0.d0+0.d0)										&
					+2.d0*(A(i,j-1)+A(i,j+1)+0.d0+A(i-1,j))+4.d0*A(i,j))/12.d0
			! Left Column -------------------------------------------------------------------------!
			IF (j==1       .AND. (i/=1 .OR. i/=size(1)))										&
				B(c,d)=((0.d0+A(i-1,j+1)+0.d0+A(i+1,j+1))										&
					+2.d0*(0.d0+A(i,j+1)+A(i+1,j)+A(i-1,j))+4.d0*A(i,j))/12.d0
			! Right Column ------------------------------------------------------------------------!
			IF (j==size(2) .AND. (i/=1 .OR. i/=size(1)))										&
				B(c,d)=((A(i-1,j-1)+0.d0+A(i+1,j-1)+0.d0)										&
					+2.d0*(A(i,j-1)+0.d0+A(i+1,j)+A(i-1,j))+4.d0*A(i,j))/12.d0
		!------------------------------------------------------------------------------------------!
		! The four corner points
			! Top Left ----------------------------------------------------------------------------!
			IF (i==1       .AND. j==1      )	B(c,d)=((0.d0+0.d0+0.d0+A(i+1,j+1)) 			&
										+2.d0*(0.d0+A(i,j+1)+A(i+1,j)+0.d0) + 4.d0*A(i,j))/9.d0 
			! Top Right ---------------------------------------------------------------------------!
			IF (i==1	   .AND. j==size(2))	B(c,d)=((0.d0+0.d0+A(i+1,j-1)+0.d0) 			&
										+2.d0*(A(i,j-1)+0.d0+A(i+1,j)+0.d0) + 4.d0*A(i,j))/9.d0 
			! Bottom Left -------------------------------------------------------------------------!
			IF (i==size(1) .AND. j==1      )	B(c,d)=((0.d0+A(i-1,j+1)+0.d0+0.d0) 			&
										+2.d0*(0.d0+A(i,j+1)+0.d0+A(i-1,j)) + 4.d0*A(i,j))/9.d0 
			! Bottom Right ------------------------------------------------------------------------!
			IF (i==size(1) .AND. j==size(2))	B(c,d)=((A(i-1,j-1)+0.d0+0.d0+0.d0) 			&
										+2.d0*(A(i,j-1)+0.d0+0.d0+A(i-1,j)) + 4.d0*A(i,j))/9.d0 
		!------------------------------------------------------------------------------------------!
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE restriction
	
	SUBROUTINE prolongation(size,A,B)
	!----------------------------------------------------------------------------------------------!
	! refines input matrix (A) into matrix (B) (interpolation)
	INTEGER, DIMENSION(2), 						INTENT(in)	:: size
	REAL*8, DIMENSION(size(1),size(2)), 		INTENT(in)	:: A
	REAL*8, DIMENSION(2*size(1)-1,2*size(2)-1), INTENT(out)	:: B
	INTEGER													:: c,d

	DO i=1,size(1)
		DO j=1,size(2)
			c=2*i-1
			d=2*j-1			
									   			B(c  ,d  )= A(i,j)									  ! center
			IF (i/=1				       )	B(c-1,d  )=(A(i,j)+A(i-1,j  ))/2.d0				  	  ! down
			IF (i/=size(1)    			   )	B(c+1,d  )=(A(i,j)+A(i+1,j  ))/2.d0				  	  ! up
			IF (j/=1  				       )	B(c  ,d-1)=(A(i,j)+A(i  ,j-1))/2.d0				  	  ! left
			IF (j/=size(2)  			   )	B(c  ,d+1)=(A(i,j)+A(i  ,j+1))/2.d0				  	  ! right
			IF (i/=1       .AND. j/=1      )	B(c-1,d-1)=(A(i,j)+A(i-1,j-1)+A(i,j-1)+A(i-1,j))/4.d0 ! down left
			IF (i/=1       .AND. j/=size(2))	B(c-1,d+1)=(A(i,j)+A(i-1,j+1)+A(i,j+1)+A(i-1,j))/4.d0 ! down right
			IF (i/=size(1) .AND. j/=1      )	B(c+1,d-1)=(A(i,j)+A(i+1,j-1)+A(i,j-1)+A(i+1,j))/4.d0 ! up left
			IF (i/=size(1) .AND. j/=size(2))	B(c+1,d+1)=(A(i,j)+A(i+1,j+1)+A(i,j+1)+A(i+1,j))/4.d0 ! up right
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE prolongation
	
	SUBROUTINE savefile(k,size,B)
	!----------------------------------------------------------------------------------------------!
	! This subroutine writes the input matrix to a file for each instance of time, based on the 
	! time iteration number (k)
	INTEGER, 							INTENT(in)	:: k
	INTEGER,DIMENSION(2), 				INTENT(in)	:: size
	REAL*8, DIMENSION(size(1),size(2)), INTENT(in)	:: B
	REAL*8, DIMENSION(size(1),size(2))				:: A
	CHARACTER(LEN=12)								:: filename
	!----------------------------------------------------------------------------------------------!
	
	A = B	! allows input matrix to be modified to be altered 
	
	CALL namefile(k,filename)

	OPEN (unit=1, file=filename)
	DO i=1,size(1)
		DO j=1,size(2)
			!----------------------------------------------------------------!
			! Ensure exponent is only 2 digits (phi is assumed positive)
			IF (A(i,j) < 1.0d-99)	A(i,j) = 1.0d-99
			IF (A(i,j) > 9.9d+99)	A(i,j) = 9.9d+99
			!----------------------------------------------------------------!			
			IF (j==1) THEN
				IF (i==1) 	THEN
					WRITE(1, '(e15.7e2)',advance='no') A(i,j)
				ELSE
					WRITE(1,'(/e15.7e2)',advance='no') A(i,j)
				ENDIF
			ELSE
		 		WRITE(1,'(tr3,e15.7e2)', advance='no') A(i,j)
		 	ENDIF
		ENDDO
	ENDDO
	CLOSE (1)
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE savefile
	
	SUBROUTINE namefile(k,filename)
	!----------------------------------------------------------------------------------------------!
	! this names the file to be saved in subroutine savefile. 
	! the number of digits of the file number leads to different cases
	INTEGER, 			INTENT(in) 	:: k
	CHARACTER(LEN=12), 	INTENT(out)	:: filename
		!------------------------------------------------------------------------------------------!
		IF (k<10)							THEN	! File number is 1 digit
											WRITE(filename,'(a8,i1,a2)'), "R//r0000",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-54.5)<45.5)			THEN	! File number is 2 digits
											WRITE(filename,'(a7,i2,a2)'), "R//r000",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-549.5)<450.5)		THEN	! File number is 3 digits
											WRITE(filename,'(a6,i3,a2)'), "R//r00",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-5499.5)<4500.5)		THEN	! File number is 4 digits
											WRITE(filename,'(a5,i4,a2)'), "R//r0",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-54999.5)<45000.5)	THEN	! File number is 5 digits
											WRITE(filename,'(a4,i5,a2)'), "R//r",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE
			PRINT *, "The number of points is greater than 5 digits"
			STOP
		ENDIF
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE namefile

	SUBROUTINE analytical(k,t)
	!----------------------------------------------------------------------------------------------!
	REAL*8								:: x, y, tmin
	REAL*8, DIMENSION(xpoints,ypoints) 	:: phi
	REAL*8,  INTENT(in)					:: t
	INTEGER, INTENT(in)					:: k
	CHARACTER(LEN=12)					:: filename	
	INTEGER								:: n
 
	tmin = 1.d-1
	!----------------------------------------------------------------------------------------------!
	DO j=1,ypoints
		DO i=1,xpoints
			x = DBLE(i)/DBLE(xpoints)*sys_length(1)-sys_length(1)*.5d0	! DBLE converts an integer 
		 	y = DBLE(j)/DBLE(ypoints)*sys_length(2)-sys_length(2)*.5d0	! to double percision real
!			phi(i,j) = sqrt(tmin/t)*exp(-(x-1.65d0)**2/(4.d0*t))
			phi(i,j) = (tmin/t)*exp(-(x**2+y**2)/(4.d0*t))
		ENDDO
	ENDDO
	!----------------------------------------------------------------------------------------------!
	IF 		(mod(k,1)==0)			THEN	
	!----------------------------------------------------------------------------------------------!
		IF (k<10)							THEN	! File number is 1 digit
											WRITE(filename,'(a8,i1,a2)'), "A//a0000",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-54.5)<45.5)			THEN	! File number is 2 digits
											WRITE(filename,'(a7,i2,a2)'), "A//a000",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-549.5)<450.5)		THEN	! File number is 3 digits
											WRITE(filename,'(a6,i3,a2)'), "A//a00",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-5499.5)<4500.5)		THEN	! File number is 4 digits
											WRITE(filename,'(a5,i4,a2)'), "A//a0",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE IF (abs(k-54999.5)<45000.5)	THEN	! File number is 5 digits
											WRITE(filename,'(a4,i5,a2)'), "A//a",k,".m"
		!------------------------------------------------------------------------------------------!
		ELSE
			PRINT *, "The number of points is greater than 5 digits"
			STOP
		ENDIF
		!------------------------------------------------------------------------------------------!
		OPEN (unit=1, file=filename)
		DO i=1,ypoints
			DO j=1,xpoints
				!----------------------------------------------------------------!
				! Ensure exponent is only 2 digits (phi is assumed positive)
				IF (phi(i,j) < 1.0d-99)	phi(i,j) = 1.0d-99
				IF (phi(i,j) > 9.9d+99)	phi(i,j) = 9.9d+99
				!----------------------------------------------------------------!			
				IF (j==1) THEN
					IF (i==1) 	THEN
						WRITE(1, '(e15.7e2)',advance='no') phi(i,j)
					ELSE
						WRITE(1,'(/e15.7e2)',advance='no') phi(i,j)
					ENDIF
				ELSE
			 		WRITE(1,'(tr3,e15.7e2)', advance='no') phi(i,j)
		 		ENDIF
			ENDDO
		ENDDO
		CLOSE (1)
		!------------------------------------------------------------------------------------------!
	ENDIF
	!----------------------------------------------------------------------------------------------!
	END SUBROUTINE analytical
	
END MODULE COMMON
!--------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------
PROGRAM multi
	USE COMMON
	!----------------------------------------------------------------------------------------------!
	! declare variables
	IMPLICIT NONE		! this prevents variables defined elsewhere to be used here
	INTEGER, DIMENSION(2)					:: grid0,grid1,grid2	! 3 grid dimensions
	REAL*8, ALLOCATABLE, DIMENSION(:,:)		:: dE0,dE1,dE2			! 3 grids for change in energy
	REAL*8, ALLOCATABLE, DIMENSION(:,:)		:: res0,res1,res2			! 3 grids for change in residual
	REAL*8, ALLOCATABLE, DIMENSION(:,:)		:: Eold,Enew			! E stands for Energy
	REAL*8									:: tmax,t,tmin			! time variables
	REAL*8									:: mag,mag_old			! residual magnitues for convergence
	INTEGER									:: k,status,n			! status is 0 for success
	!----------------------------------------------------------------------------------------------!
	! variables to check elapsed time
	REAL				:: etime, total
	REAL, DIMENSION(2)	:: elapsed
	!----------------------------------------------------------------------------------------------!
	! Clear folders containing data
	CALL system("RM R/r*.m")		!Deletes old files (would overide files with same name anyways)
	CALL system("RM A/a*.m")		!Deletes old analytic files
	!----------------------------------------------------------------------------------------------!
	! Set Constants
	xpoints 		= 33				! must be of the form 2^n+1: {2,3,5,9,17,33,65,129,257,...}
	ypoints 		= xpoints			! x and y points are equal for a square grid
	sys_length(1) 	= 3.3d0		! length of physical system to be modeled.  Default is a unit square
	sys_length(2) 	= 3.3d0		! length of physical system to be modeled.  Default is a unit square
	hx 				= sys_length(1)/xpoints 	! this is the grid spacing: 1/xpoints assumes a unit square grid
	hy 				= sys_length(2)/ypoints
	tmax 			= 1.d0				! maximum time of simulation (seconds)
	tmin 			= 1.d-1
	t 				= tmin
	tol 			= 1.d-6					! tolerance of SGS solver
	k 				= 0
	BCflag 			= 'Period'				! Sets flag for boundary conditions  - Options: NoFlux, Period, Dirich
	dt				= 3.5d-2 !2*(0.5*(hx+hy))**2		! stability: D*dt/4/h^2 < 1  ----> dt<4*h^2 	since D=1, 1/h=xpoints
	! dt is the curxent time step.  now it is a constant, but a future subroutine can dynamically modify this on the fly
	!----------------------------------------------------------------------------------------------!
	! Set grid sizes
	grid0(1)=xpoints
	grid0(2)=ypoints
	grid1(1)=(1+grid0(1))/2
	grid1(2)=(1+grid0(2))/2
	grid2(1)=(1+grid1(1))/2
	grid2(2)=(1+grid1(2))/2
	!----------------------------------------------------------------------------------------------!
	! Allocate variables
	ALLOCATE(dE0(grid0(1),grid0(2)),stat=status)
	ALLOCATE(dE1(grid1(1),grid1(2)),stat=status)
	ALLOCATE(dE2(grid2(1),grid2(2)),stat=status)

	ALLOCATE(res0(grid0(1),grid0(2)),stat=status)
	ALLOCATE(res1(grid1(1),grid1(2)),stat=status)
	ALLOCATE(res2(grid2(1),grid2(2)),stat=status)

	ALLOCATE(Eold(grid0(1),grid0(2)),stat=status)
	ALLOCATE(Enew(grid0(1),grid0(2)),stat=status)
	!----------------------------------------------------------------------------------------------!
	! initialize variables
	CALL init(Eold)
	Enew = Eold
	CALL savefile(0,grid0,Enew)	
	CALL analytical(0,1.d-1)
	PRINT "(a16,f6.4,a11,f6.2,a17)", "The time is now ",t,", which is ",(t-tmin)/(dt*CEILING((tmax-tmin)/dt))*100,"% of the way done"
	CALL magnitude(grid0,Eold,Enew,mag)
	mag_old = mag
	!----------------------------------------------------------------------------------------------!
	! start time progression
	DO WHILE (t<tmax)
		! calculate new dt:	call t_update(dt)
		! if a new dt is found, then a new matrix A must be calculated: call discretization(A)
		t 	= t + dt
		k 	= k+1
		dE0 = 0.d0
		! CALL magnitude(grid0,Eold,Enew,mag)	! this output is overwritten on next line
		mag=1.d0	!mag_old = mag	
		! set mag_old to be high enough to enter while loop
		DO WHILE (mag>1.d-4) !(mag>1e-2*mag_old)
			!mag_old=mag
			!--------------------------------------------------------------------------------------!
			!presmoothing		! number of calls is equal to number of presmooth cycles
				CALL calc_res(grid0,Eold,Enew,res0)
				CALL smoother(grid0,res0,dE0)
				Enew = Enew + dE0
				!CALL URfact(grid0,dE0,Enew)		! update Enew with an under-relaxation factor	
				!----------------------------------------------------------------------------------!
				CALL calc_res(grid0,Eold,Enew,res0)
				CALL smoother(grid0,res0,dE0)
				Enew = Enew + dE0
				!CALL URfact(grid0,dE0,Enew)		! update Enew with an under-relaxation factor	
			!--------------------------------------------------------------------------------------!
			! coarse grid correction
!			DO n=1,15
				CALL calc_res(grid0,Eold,Enew,res0)
				CALL restriction(grid0,res0,res1)
				CALL restriction(grid0,dE0,dE1)
				!----------------------------------------------------------------------------------!
				CALL restriction(grid1,res1,res2)
				CALL restriction(grid1,dE1,dE2)	
				!----------------------------------------------------------------------------------!
				CALL smoother(grid2,res2,dE2)			
				!----------------------------------------------------------------------------------!			
				CALL prolongation(grid2,dE2,dE1)
				CALL prolongation(grid1,dE1,dE0)
				Enew = Enew + dE0
				!call URfact(grid0,dE0,Enew)
!			ENDDO
			!--------------------------------------------------------------------------------------!
			! postsmoothing
				CALL calc_res(grid0,Eold,Enew,res0)
				CALL smoother(grid0,res0,dE0)
				Enew = Enew + dE0
				!CALL URfact(grid0,dE0,Enew)
				!----------------------------------------------------------------------------------!
				CALL calc_res(grid0,Eold,Enew,res0)
				CALL smoother(grid0,res0,dE0)
				Enew = Enew + dE0
				!CALL URfact(grid0,dE0,Enew)			
			!--------------------------------------------------------------------------------------!
			! Calculate magnitude of residual for time step convergence
				CALL magnitude(grid0,Eold,Enew,mag)
		!------------------------------------------------------------------------------------------!
		ENDDO		! this ends the dE while loop: while (mag>1e-2*mag_old)
		!------------------------------------------------------------------------------------------!
		
		!------------------------------------------------------------------------------------------!
		! this section decides which files to be written	
		IF 		(mod(k,1)==0)			THEN
			PRINT "(a16,f6.4,a11,f6.2,a17)", "The time is now ",t,", which is ",(t-tmin)/(dt*CEILING((tmax-tmin)/dt))*100,"% of the way done"
			CALL savefile(k,grid0,Enew)			
		ELSE IF (t==tmax)				THEN
			PRINT "(a16,f6.4,a11,f6.2,a17)", "The time is now ",t,", which is ",(t-tmin)/(dt*CEILING((tmax-tmin)/dt))*100,"% of the way done"
			CALL savefile(k,grid0,Enew)	
		ENDIF
		!------------------------------------------------------------------------------------------!
		Eold = Enew 	! update the 'old' E as this time step is now complete
		CALL analytical(k,t)
	ENDDO				! this ends the time-stepping while loop: while (t<tmax)
	!----------------------------------------------------------------------------------------------!
	! Print the runtime of the program
	! Define(User time):  	time actually spent in your program;
	! Define(System time):	time spent in the operating system on the program's behalf
	total = etime(elapsed)
	PRINT *, 'Total time=',total,'User runtime =', elapsed(1), ' System runtime =', elapsed(2)
!--------------------------------------------------------------------------------------------------!
END PROGRAM multi