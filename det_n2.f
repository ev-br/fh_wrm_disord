!-----------------------------------------------------------
!  adding/dropping a row and a column via N**2 operations 
!   +  full N**3 inverse 
!   +  changing a single row
!  using BLAS / LAPACK
!-----------------------------------------------------------
	module det_n2
	implicit none; save

!------------------------------------
! Implementation issues:
!   1) inv_p1: need to set a vector = 0; 
!       Currnet implementation uses dscal; 
!       Will it be any better to use dcopy from a pre-set attay of zeros?
!   2) inv_p1 & inv_m1: outer product of two vectors;
!        Current implementation uses dger; Would dgemm() do better?
!   3) inv_m1: row & column to be dropped are swapped w/the last ones,
!         and the former are then being dropped. 
!         The SIGN of the determinant might get wrong, and it is NOT
!         being taken care of. 
!------------------------------------

	public :: full_inv, det_p1, inv_p1, det_m1, inv_m1, det_r, inv_r
	public :: det_p2

	private

	real*8, external :: ddot



!*****************************************************************
	contains



!-------------------------------
!--- Det: add a row & a column
!-------------------------------
	real*8 function  det_p1(pm,lda,a,u,v,z,s)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda)
      real*8  :: v(lda),u(lda),z(lda),s
!
!  Input: a(lda,lda), s, v(lda), u(lda)
!         all arrays are of the size pm<lda 
!
!  Output: z(lda)
!          det. ratio @ det_p1

!  z=a1 DOT u
	call dgemv('N',pm,pm,1.d0,a,lda,u,1,0.d0,z,1)


!  \lambda = v DOT z; \rho = s - \lambda; det = 1/\rho
	det_p1 = s - ddot(pm,v,1,z,1)

	end function det_p1


!-------------------------------
!--- Add row & column: update the inverse 
!-------------------------------
	subroutine  inv_p1(pm,lda,a,det,v,w,z)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda),det,v(lda),w(lda),z(lda) 

	real*8 :: rho
!
!  Input: a(lda,lda), v(lda)
!         z(lda), det --- as set by det_p1
!         
!          all arrays are of the size pm<lda 
!
!  Output: a(lda,lda) ==> (pm+1) * (pm+1) inverse
!          pm = pm( @ input) + 1

	rho = 1.d0/det	

!  w = transp(a1)*v
	call dgemv('T',pm,pm,1.d0,a,lda,v,1,0.d0,w,1)

! b^{-1} 
	call dger(pm,pm,rho,z,1,w,1,a,lda)

! last row
       call dscal(pm,0.d0,a(pm+1,1),lda) ! set =0 first
       call daxpy(pm, -rho, w,1,a(pm+1,1),lda)

! last column
       call dscal(pm,0.d0,a(1,pm+1),1) ! set =0 first
       call daxpy(pm, -rho, z,1,a(1,pm+1),1)

! (pm+1,pm+1)
	a(pm+1,pm+1)=rho

! update pm
	pm=pm+1

	end subroutine inv_p1



!-------------------------------
!--- Det: drop r-th row & c-th column
!-------------------------------
      real*8 function  det_m1(pm,lda,a,r,c)
      integer :: pm,lda       ! actual size & leading dimension of A
      real*8  :: a(lda,lda)    
	integer:: c,r
	integer :: s
!
!  Input:  a(lda,lda), pm
!
!  Output: det. ratio:  det BIG / det SMALL (=1/rho)
!

!!! ******** TODO: razobrat'sya s SIGN'om ********


c      det_m1 = 1.d0*(-1)**(c+r)/a(c,r)

	s=-1; if(c/=pm)s=-s; if(r/=pm)s=-s
	det_m1 = 1.d0*s/a(c,r)

      end function det_m1


!-------------------------------
!--- Inv: drop r-th row and c-th column
!-------------------------------
      subroutine  inv_m1(pm,lda,a,r,c)
      integer :: lda, pm       ! leading dimension and actual size of A
      real*8  :: a(lda,lda)
      integer :: c,r

	real*8 :: rh

!  Input: a(lda,lda) pm*pm
!         r, c --- row and  column to be dropped
!
!  Output: a(lda,lda) ==> (pm-1) * (pm-1)
!          pm = pm( @ input ) - 1
!
!   How:  swaps to-be-dropped and last row & cols, 
!         and then drops the former ones
!
!   The latter is done using: 
!       a(pm,pm)         is   \rho
!       a(1:pm,pm)       is  -\rho * z  -- last column
!       a(pm,1:pm)       is  -\rho * w  -- last row
!       a(1:pm-1,1:pm-1) is  a^{-1} + \rho z*w 

! swap c-th and last row (c-th column of A <==> c-th row of A^{-1})
	if(c/=pm) call dswap(pm,a(c,1),lda,a(pm,1),lda)
      
! swap r-th and last column
	if(r/=pm) call dswap(pm,a(1,r),1,a(1,pm),1)
	
!-------------------------------- drop the last row & column
	rh = -1.d0/a(pm,pm)

! stripe out the outer product z*w
	call dger(pm-1,pm-1,rh,a(1,pm),1,a(pm,1),lda,a,lda)

! update pm
	pm = pm - 1


      end subroutine inv_m1


!-------------------------------
!--- Det: change a single row
!-------------------------------
	real*8 function  det_r(pm,lda,a,r,v)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda), v(lda)
      integer :: r
!
!  Input: a(lda,lda) --- inverse matrix
!         v(lda)     --- a row to be added
!         r          --- a row number 
!           all arrays are of the size pm<lda 
!
!  Output:
!           det. ratio @ det_r

! \lambda =  ( v DOT A_r )
      det_r = 1.d0 + ddot(pm,v,1,a(1,r),1)

	end function det_r


!-------------------------------
!--- inv: change a single row
!-------------------------------
	subroutine  inv_r(pm,lda,a,r,det,v,w,z)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda), v(lda),w(lda),z(lda),det
      integer :: r

      real*8  :: rho
!
!  Input: a(lda,lda) --- inverse matrix
!         v(lda)     --- a row to be added
!         w(lda),z(lda)     --- working arrays
!         r          --- row number 
!         det        --- det ratio, as set by det_r
!           all arrays are of the size pm<lda 
!
!  Output:
!         a(lda,lda) contains an updated inverse matrix


      rho = -1.d0/det

! z_i = A_{i,r}
      call dcopy(pm,a(1,r),1,z,1)

! w = v DOT A 
      call dgemv('T',pm,pm,1.d0,a,lda,v,1,0.d0,w,1)


! A+ \rho* z DOT w^T
      call dger(pm,pm,rho,z,1,w,1,a,lda)

! one cannot get rid of z() array since if one otherwise plugs {a(1,r),1}
! directly into dger(...) instead of {z,1}, it all goes nuts.
! probably, dger() spoils the z array via blocking or the like.

	end subroutine inv_r



!-------------------------------
!--- inv & det of a matrix (honest N**3)
!-------------------------------
	real*8 function  full_inv(pm,lda,a)
	integer :: lda, pm  ! leading dimension and actual size of A
	real*8  :: a(lda,lda)

	integer, allocatable :: ipiv(:)
	real*8, allocatable :: work(:)
	integer :: info, i,lwork,icnt
!
!  Input: A(lda,lda), pm<lda 
!
!  Output: A contains the inverse
!          full_inv contains the determinant

	lwork=lda
	allocate(ipiv(1:pm), work(1:lwork))

	call dgetrf(pm,pm,a,lda,ipiv, info) 

        full_inv = 1d0
         icnt=0
         do i=1,pm
           full_inv = full_inv * a(i,i)
           if (ipiv(i).ne.i) then
             icnt = icnt+1
           endif
         enddo
         if (mod(icnt,2).eq.1) full_inv = -full_inv



	call dgetri(pm,a,lda,ipiv,work,lwork,info )
	if(info/=0)print*,'dgetri info = ',info

	deallocate( work,ipiv )

	end function full_inv



!-------------------------------
!--- Det: add a two rows & columns
!-------------------------------
	real*8 function  det_p2(pm,lda,a,u2,v2,s2,c2)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda)
      real*8  :: v2(2,lda),u2(lda,2),s2(2,2),c2(lda,2)

!
!  Input: a(lda,lda), v2(2,lda), u2(lda,2), s2(2,2)
!         all arrays are of the size pm<lda 
!
!  Output:    det. ratio @ det_p2
!

!  c2=a^{-1} DOT u2
	call dgemm('N','N',pm,2,pm,1.d0,a,lda,u2,lda,0.d0,c2,lda)

!  \Rho = s2 - v2 DOT c2  
	call dgemm('N','N',2,2,pm,-1.d0,v2,2,c2,lda,1.d0,s2,2)

! det = det \Rho [which is 2*2]
	det_p2 =  s2(1,1)*s2(2,2) - s2(1,2)*s2(2,1)

	end function det_p2





	end module det_n2
