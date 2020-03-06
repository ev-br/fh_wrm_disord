!
!   Analyze the statistics
!
	implicit none

!------------------- Phys. parameters --------------------

	real*8,parameter :: pi=3.1415926d0

      integer, parameter  :: d=3, dd=6            ! d is the global dimension, dd=d+d  
	integer, allocatable :: ddd(:)
      integer, allocatable :: N(:,:)           ! Array for Numbers of sites per dimensions
      integer, allocatable :: Nsite(:)            ! Total number of sites

	real*8,allocatable :: U(:)          ! MINUS interaction strength
	real*8,allocatable :: beta(:)       ! inverse temperature
	real*8,allocatable :: mu(:)         ! chem. potential

!--------------  Measurement related -------------------

	integer, parameter   :: b_n_max=100
	integer,allocatable  :: Z_b(:)           ! # block size
	integer,allocatable  :: b_n(:), i_b(:)    ! filled block nbr. & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block
! partition function
	real*8,allocatable :: Z(:)

! diagonal estimators:
	real*8,allocatable  :: PE(:), KE(:), ndens(:)     ! PE, KE and density estimators 
	real*8,allocatable  :: PE_stat(:,:), KE_stat(:,:), ndens_stat(:,:) ! files

! integrated ira-masha correlator:
	real*8,allocatable  :: im(:)               ! estimator
	real*8,allocatable  :: im_stat(:,:) ! file

! density-density
c	real*8, allocatable :: g_uu(:,:), g_ud(:,:)

	real*8, allocatable :: step(:)
 
!-----------  overall stat ---------------
	integer :: b_n_all, i_all
	real*8, allocatable :: PE_all(:), KE_all(:),im_all(:),ndens_all(:)

	real*8, allocatable :: ene_sc(:)  ! (PE + KE per fermion ) / E_F
	real*8, allocatable :: tef(:)     ! T/E_F


	real*8, allocatable :: g_uu_all(:), g_ud_all(:)

	real*8 :: step_all, Z_all

	real*8 :: av, err
	real*8 :: ef, def,kf

!-------------------- discard what -----------
	integer, allocatable :: istart(:,:)
	integer :: istart_num                ! how many tries
	integer :: klm, tmp
	logical, allocatable ::  empty(:)


!-------------------- working variables -----------------

	integer :: nfls     ! # of stat. files to process
	integer :: i,j,p

	integer :: who_max,Z_b_max

	real*8 :: eef, npart

	character*50 :: fname, fsuffix

!=========================================


!-------------- read the files ----------------

	print*,' Enter # of files: '; read*,nfls

	
	allocate( ddd(nfls),N(nfls,d),Nsite(nfls) )
	allocate( U(nfls),beta(nfls),mu(nfls) ) 
	allocate( Z_b(nfls),b_n(nfls),i_b(nfls),Z(nfls),step(nfls) )
	allocate( PE(nfls),KE(nfls),ndens(nfls),im(nfls) )
	allocate( PE_stat(nfls,b_n_max),KE_stat(nfls,b_n_max) )
	allocate( ndens_stat(nfls,b_n_max),im_stat(nfls,b_n_max) )
c	allocate( g_uu(nfls,0:100), g_ud(nfls,0:100) )


	do i=1,nfls; 
	   print*,' Enter filename #',i,' : '; read*,fname
	   call rd_stat(i,fname)
	   print*,i,' : ',trim(fname),' is read.'
	enddo

c	read*,fsuffix


!-------------- check files ---------------------
	do i=1,nfls
	  if( any ( mu/=mu(1) ) )then; print*,'mu...',mu(:)
	  endif
	  if( any ( beta/=beta(1) ) )then; print*,'beta...',beta(:)
	  endif
	  if( any ( U/=U(1) ) )then; print*,'U...',U(:)
	  endif
	  if( any ( N(:,1)/=N(1,1) ) )then; print*,'N(,1)...',N(:,1)
	  endif
	  if( any ( N(:,2)/=N(1,2) ) )then; print*,'N(,2)...',N(:,2)
	  endif
	  if( any ( N(:,3)/=N(1,3) ) )then; print*,'N(,3)...',N(:,3)
	  endif
	enddo


!--------  do density-density correlator ----------------
c	allocate( g_uu_all(0:100), g_ud_all(0:100) )
c
c	Z_all = sum(Z(:))
c
c	do i=0,N(1,1)/2; 
c	  g_uu_all(i)=sum ( g_uu(:,i) ) / nfls / Z_all
c	  g_ud_all(i)=sum ( g_ud(:,i) ) / nfls / Z_all
c	enddo
c
c	open(1,file='g_uu'//trim(fsuffix)//'.dat')
c	open(2,file='g_ud'//trim(fsuffix)//'.dat')
c	  do i=0,N(1,1)/2
c	       write(1,*)i,g_uu_all(i)
c	       write(2,*)i,g_ud_all(i)
c	  enddo
c	close(1); close(2)
c
c	deallocate(g_uu_all, g_ud_all)

!---------------------------------------------------------


!-------------- equate block sizes ----------------

	Z_b_max = maxval( Z_b )

c	print*,'before: ',Z_b, Z_b_max


	do i = 1,nfls
	
	  p = Z_b_max / Z_b(i); !print*,'i,p(i) =  ',i,p
	  Z_b(i) = Z_b(i)*p; b_n(i) = b_n(i)/p

	  do j = 1,b_n(i);  	 
	     PE_stat(i,j) = 1.d0*sum( PE_stat(i,p*(j-1)+1:p*j) ) / p
		 KE_stat(i,j) = 1.d0*sum( KE_stat(i,p*(j-1)+1:p*j) ) / p
	     im_stat(i,j) = 1.d0*sum( im_stat(i,p*(j-1)+1:p*j) ) / p
	     ndens_stat(i,j) = 1.d0*sum( ndens_stat(i,p*(j-1)+1:p*j) ) / p
	  enddo

	enddo

c	print*,'after: ',Z_b, Z_b_max

! ------------  where to start from --------------------
	istart_num=3
	allocate( istart(1:nfls, 1:istart_num) )
	
	do i=1,nfls
	  istart(i,:)=(/ 1,b_n(i)/4,b_n(i)/2 /)     ! everyth, drop 1st quarter, 1st half
	enddo


	  do klm=1,istart_num


! -------------  glue all the files -----------------
	print*; print*; 
	print*,'-------------------------------------------------' 
	print*


	allocate(empty(1:nfls)); empty=.false.

  	b_n_all = 0 !sum(b_n)
	do i=1,nfls
	   tmp = b_n(i) - istart(i,klm) + 1
	   b_n_all = b_n_all + tmp
	   if(tmp==0)then; print*,'No stat. left in #', i
	                   empty(i)=.true. 
	   endif
	enddo
	if(b_n_all==0)then; print*,'No stat. left whatsoever...'
	                    stop
	endif


	allocate( PE_all(b_n_all), KE_all(b_n_all), im_all(b_n_all) )
	allocate( ndens_all(b_n_all), ene_sc(b_n_all),tef(b_n_all) )

	i_all=0; Z_all = 0

	do i=1,nfls;  
	
	   if(empty(i))cycle
	  
	  do j=istart(i,klm),b_n(i)
	     i_all = i_all + 1

	     ndens_all(i_all) = ndens_stat(i,j)

	     eef = (3.d0*pi*pi*ndens_all(i_all))**(2.d0/3.d0)
	     npart = ndens_all(i_all)*Nsite(1)

	     PE_all(i_all) = -PE_stat(i,j) / npart / eef

	     KE_all(i_all) = ( 6.d0 - 2.d0*KE_stat(i,j)/npart ) / eef

	     ene_sc(i_all) = KE_all(i_all) + PE_all(i_all)

	     tef(i_all) = 1.d0 / beta(1) / eef

		 
	  enddo

	enddo





	Z_all = Z_b(1) * i_all

	step_all = sum(step); 



! ---------------  prntout ------------------
	
	print*,d, N(1,:)
	print*,'beta = ', beta(1), ' U = ', U(1), ' mu = ', mu(1)
	print*

	print*,'MC step (mln) = ', step_all/1.d6
	print*,' Z  = ',Z_all/1.d6


!--- pot. energy -----------------------------------
	print*; print*,'doing PE: '
	call mrg(PE_all(1:b_n_all),b_n_all,1)


      
!--- kin. energy -----------------------------------
	print*; print*,'doing KE: '
	call mrg(KE_all(1:b_n_all),b_n_all,1)



!---------  scaled energy ------------------------
	print*; print*,'doing scaled energy: '
	call mrg(ene_sc(1:b_n_all),b_n_all,1)
	

!--- ndens
	print*; print*,'doing denity: '
	call mrg(ndens_all(1:b_n_all),b_n_all,1)


!--- T / E_F
	print*; print*,'doing T / E_F: '
	call mrg(tef(1:b_n_all),b_n_all,1)



!--- condensate
	print*; print*,'doing g_im(w=0,k=0): '
	call mrg(im_all(1:b_n_all),b_n_all,1)



	deallocate(PE_all,KE_all,ndens_all,im_all, empty, ene_sc,tef)

		enddo   ! klm



	contains  !========================================

!---------------------
!--- Read statistics
!---------------------
	subroutine rd_stat(i,fname)
	integer :: i
	character*50 :: fname
	character*6  :: cvoid


	open(1,file=trim(fname))
	  read(1,*)ddd(i),N(i,:)
	  read(1,*)beta(i),U(i),mu(i)
	  read(1,*)step(i), Z(i)
	  read(1,*)Z_b(i), b_n(i), i_b(i)
	  read(1,*)PE(i)
	  read(1,*)PE_stat(i,1:b_n(i));   
	  read(1,*)KE(i)
	  read(1,*)KE_stat(i,1:b_n(i))
	  read(1,*)ndens(i)
	  read(1,*)ndens_stat(i,1:b_n(i))
	  read(1,*)im(i)
	  read(1,*)im_stat(i,1:b_n(i))
c	  read(1,*)cvoid
c	  read(1,*)g_uu(i,0:N(i,1)/2)
c	  read(1,*)g_ud(i,0:N(i,1)/2)
	close(1)

	Nsite(i)=N(i,1)*N(i,2)*N(i,3)

	PE_stat(i,1:b_n(i)) = PE_stat(i,1:b_n(i)) / Z_b(i)
	KE_stat(i,1:b_n(i)) = KE_stat(i,1:b_n(i)) / Z_b(i)
	ndens_stat(i,1:b_n(i)) = ndens_stat(i,1:b_n(i)) / Z_b(i)
	im_stat(i,1:b_n(i)) = im_stat(i,1:b_n(i)) / Z_b(i)

	end subroutine rd_stat


!-------------------------------
!--- Analyze block statistics
!-------------------------------
	subroutine bSTAT(arr,n,Zb,av,err)
	integer              :: n, Zb
	real*8, dimension(1:n) :: arr
	real*8               :: av, err

	real*8 :: av2

	av  = sum( arr(1:n) )/Zb/n
	av2 = sum( arr(1:n)**2 )/Zb/Zb/n

				!av2 = av2 + (arr(j)/Zb)**2

	err = sqrt( av2 - av*av ) / sqrt(1.d0*n)


	end subroutine bSTAT


!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
	subroutine mrg(arr,n,Zb)
	integer, intent(in)              :: n, Zb
	real*8, dimension(1:n), intent(in) :: arr

	real*8  :: av, err, arr1(1:n)
	integer :: i, n1, zb1


	zb1 = zb; 	arr1(1:n) = arr(1:n); n1=n

	print*,'-----------'
	
	do;
	
! emit
	call bSTAT(arr1,n1,zb1,av,err)
      print 777, av, err,n1
 777  format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
	if(n1<3)exit

! merge
	n1=INT(n1/2); zb1=zb1*2
	do i=1,n1
	    arr1(i) =  arr1(2*i-1) + arr1(2*i)
	enddo

	enddo

	print*,'------------'; print*; 


	end subroutine mrg




	end
