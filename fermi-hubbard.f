!
! Determinant diagrammatic MC code for the Fermi-Hubbard model
! at finite-temperature.
!
!
!  for details see
! 
!  http://montecarlo.csi.cuny.edu/umass/fermion.html
!
! **********************************
!    This is a non-MPI version 
! **********************************

!---------------------
!--- Variables
!---------------------
	module vrbls
	implicit none

	real*8,  parameter :: nul=0.0d0, un1=1.0d0, un2 = 2.0d0
	real*8   :: pi
	integer, parameter :: i_hu = huge(1)

	integer, parameter :: nkink_m = 512  ! max # of kinks on one site
	integer, parameter :: nmnm_max = 8192 ! max total # of kinks

!------------------------------------------------------
!     Phys. parameters
!------------------------------------------------------

      integer, parameter :: d=3, dd=6       ! d is the global dimension; 
	                                      ! if you need d/=3, change it in here and in the parameter file  
										  ! It is essential for performance reasons to have the 
										  ! dimension as a static variable (othewise some functions 
										  ! can't be inlined etc)
      integer :: N(1:d)     ! Array for Numbers of sites per dimensions
	integer :: N2(1:d)
      integer :: Nsite      ! Total number of sites

	real*8 :: U, U_in    ! MINUS interaction strength, initial for thermalization
	real*8 :: beta       ! inverse temperature
	real*8 :: mu         ! chem. potential
	real*8 :: disorder_amp   ! disorder amplitude
	real*8, allocatable  :: musite(:)    ! site-dependent chem.potential
	integer              :: replica_id  ! replica ID (read from the disorder file) 

      real*8 :: eta        ! GF- vs Z-sector weight

	real*8,parameter :: gf_eta=1.038d0   ! 1 + GF anomalous dimension, U(1) universality class

!--- Green function array
!d	integer, parameter :: mt_max=500    ! Max. number of the mesh points for tau
!d	integer, parameter :: ntab_max=20 ! Max. number of of sites per dimension for tabulation
	integer :: Ntab                   ! Actual Number of sites per dimension for tabulation
	integer :: mtau

	real*8, allocatable :: GR_DAT(:,:,:)	! Green function array
	real*8, allocatable :: GRD_DAT(:,:,:)     ! Green function derivative
 
	real*8 :: bmt, bmt1     ! a shorthand for beta/mtau, and its inverse

	real*8, allocatable :: g0000(:)       ! green function w/ all arg = 0, non-interacting density

!------------------------------------------------------
!     Sites and Associations
!------------------------------------------------------

      integer, allocatable :: ass(:,:)
      integer, allocatable :: back(:)
      integer, allocatable :: x(:,:)  ! coordinates x(1:d,site)


!     The array ass(...) specifies the nearest-neighbor sites in both positive and
!     negative i-direction. For example, in a 2-D lattice, like the one below,
!     the associations to the io site are:  
!     (here +(-)1 = positive (negative) x-direction, 
!     +(-)2 = positive (negative) y-direction)
      
!     ass(+1,i0)=i1;    ass(-1,i0)=i3;    ass(+2,i0)=i4;    ass(-2,i0)=i2
      
!                            ...    
!                             |  
!                     ... -- i2 -- ...
!                       |     |     |
!                ...-- i3 -- i0 -- i1 --.... 
!                       |     |     |
!                     ... -- i4 -- ... 
!                             |  
!                            ...
!           
!     * The array back(..) specifies the 'opposite direction', e.g. in 3D
!           back(1) = 4, back(4)=1 (positive x is 1, negative x is 4), and so forth
! 


!------------------------------------------------------
!     Update related
!------------------------------------------------------

	real*8 :: bun           !  = beta*U*Nsite 
	real*8 :: udt           !  = U*ddt_same

	integer :: n_ca, n_le             ! # of sites within the cre/ann, leap cube
	integer, allocatable :: uzli(:)   ! site list for cube()
	logical, allocatable :: yes(:)    ! flags for cube()

	integer :: rx_ca, rx_le            ! x- and \tau- radii for cre/ann, leaps
	real*8  :: rt_ca, rt_le

	integer, allocatable :: uzli_cube_ca(:,:), uzli_cube_le(:,:) ! site lists 

! mega-leaps
	integer :: nt_same             ! # of \tau points for leaps
	real*8  :: ddt_same            ! corresp. interval length
	integer :: n_ev_cut            ! # of events after the cutoff
	integer, allocatable :: dx_le(:,:)
	real*8,  allocatable :: dt_le(:), w_le(:)
	real*8 :: w_norm            ! normalization factor for w_le(:)

! nearest kink:
	integer, allocatable :: name_best(:)
	real*8, allocatable :: dist_best(:), distance(:) 

!------------------------------------------------------
!     Configuration
!------------------------------------------------------

! The data structure is purposedly redundant, for in all cases where 
! there is a memory/speed tradeoff, it is solved in favor of speed.
!
! The interaction vertices (a.k.a 'kinks')  live on the sites of a d-dimensional lattice
! and on the [0,\beta] interval of the \tau-continuum. The kinks are labelled by
! the 'names', e.g. unique numbers. The configuration consists of:
!
! 1. Overall list of kinks' names, managed by the name manager, see GetName() and DropName() subroutines
!    Given the kink name, one has access to the kink's site, temporal position (a.k.a 'time', also 'tau'),
!    and the position in TheMatrix, see below. 
!
! 2. Each site hosts a list of names of kinks which belong to this site. The order of kinks in these
!    lists is arbitrary.
!
! 3. The configuration weight matrix is not stored itself; one always deals with its inverse, 
!    a.k.a. TheMatrix :). 
!
! 4. ira/masha are not kinks. They are stored separately, see below.
!
	integer,   allocatable :: kname(:,:)   ! list of kinks's names @ a given site: kname(j,site)
	integer*2, allocatable :: nkink(:)     ! # of kinks @ site: nkink(site)

	integer, allocatable   :: ksite(:)      ! ksite(name) => site of a kink 'name'
	real*8, allocatable    :: ktau(:)       ! ktau(name) => tau of a kink 'name'

	integer, allocatable   :: row(:), clmn(:) ! row(name), clmn(name) => position of a 'name' in TheMatrix
	integer, allocatable   :: nm_row(:),nm_clmn(:)  ! nm_row(row) => name of the kink associated with the row

! NameManager data
	integer                        :: nmnm         ! total number of kinks
      INTEGER, DIMENSION(-2:nmnm_max) :: namelist, numnam

      logical :: present             ! =.true. is ira & masha are present in the configuration
      integer, parameter :: ira=-2, masha=-1

! TheMatrix
	integer             :: pm,lda      ! actual size & leading dimension
	real*8, allocatable :: matr(:,:)   ! TheMatrix, inverse
	real*8, allocatable :: m_u(:), m_v(:), m_w(:),m_z(:)       ! working arrays for rank-1 updates
	real*8              :: m_s

	real*8, allocatable :: m_u2(:,:),m_v2(:,:),m_s2(:,:),m_c2(:,:)  ! for rank-2 updates

!-------------------------------------------------------
!     Measurement related
!-------------------------------------------------------

!
!   Errorbar analysis is performed via blocking method. The maximum number of blocks is
!   fixed, once this number is reached, the blocks are collated so that one has twice
!   as few blocks of twice the size. 
!
	integer, parameter :: b_n_max=1000   ! max # of blocks
	integer  :: Z_b                     ! block size

	integer :: b_n, i_b    ! # of filled blocks  & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block
! Partition function
	real*8 :: Z            

! Diagonal estimators:
	real*8  :: PE, KE, ndens     ! PE, KE and density estimators 
	real*8  :: PE_stat(b_n_max), KE_stat(b_n_max), ndens_stat(b_n_max) ! files

	real*8  :: g_uu(0:100),g_ud(0:100)  ! Density-densiy correltaors, a.k.a 'dance-dance'


! Integrated ira-masha correlator:
	real*8  :: im               ! estimator
	real*8  :: im_stat(b_n_max) ! file
 
! Service histograms:
	integer, parameter :: im_t_max=15
	real*8  :: im_t(-im_t_max:im_t_max)          ! t(ira-masha) histogram

	integer, parameter :: nmnm_sup = 7500        ! kink # distribution
	real*8, dimension(0:nmnm_sup) :: nmnm_distr

	integer, parameter :: det_distr_max=25       ! det. distribution
	real*8 :: det_distr(-det_distr_max:det_distr_max)


!------------------------------------------------------
!      Flow control and output
!------------------------------------------------------

!
!  Since configuration weight determinant ratios (a.k.a.'det') are
!  recalculated incrementally, the roundoff error buildup is possible
!  Also, if the current configuration is close to the nodal surface,
!  the det-s might get extravagantly large/small, which also might lead to
!  error acumulation. In order to deal with these, determinants are recalculated 
!  from scratch (which costs N**3 operations!) (i) every so often (once per N**4 
!  proves to be OK), and (ii) every time unusually large det pops out (the latter
!  signals that the pre-update configuration was close to a nodal surface)
!
!
	real*8 :: tolerance    ! determinant recalculation tolerance level
	real*8 :: recalc_rate  ! recalc. rate counter

! Update manager probabilities
	real*8, dimension(1:10) :: prob

      real*8 :: step, step_p, step_w, step_m, step_t, step_r
!     step    counts MC steps
!     step_p  the number of steps between printouts
!     step_w  the number of steps between writing to disk
!     step_m  the number of steps between measurements
!     step_t  the number of steps for thermolization
!     step_r  the number of steps for det. recalculation

!  Counters for printing, writing, measuring, thermolization, det recalculation
      real*8 :: i_p, i_w, i_m, i_t, i_r

      real*8  :: n_sw   ! number of sweeps for thermalization, 1 sweep = \beta*U*Nsite updates

! Address/accept counters	for various updates
	real*8 :: c_a_v, a_a_v, c_d_v, a_d_v
	real*8 :: c_ann, a_ann, c_cre, a_cre
	real*8 :: c_l_a, a_l_a, c_l_d, a_l_d
	real*8 :: c_r, a_r

	real*8 :: i_backsc     ! 'backscattering' counter: leap_add immediately followed by a leap_drop

! Kink number output
      integer :: nm_min=0, nm_max=i_hu
      real*8  :: nm_av=0

! 'status' variable holds the name of the update; proves useful for debugging
	character*6 :: status, prevstatus

! Output channel #:
	integer, parameter :: OUT_UNIT=22        ! For MPI code, writing to the stdout is useless, 
											 ! thus all the output goes to file(s)

! rndm() seed, see the 'mumbers' module
	integer :: r_ij, r_kl
	integer, allocatable :: ar_ij(:), ar_kl(:)


! Walltime constraints
	character*8  :: date               ! "how much on your watch?" -- "Six o'cloch"
	character*10 :: time
	character*5  :: zone

	integer :: tvalues(8)
      
	real*8 :: time_prev, time_curr, time_elapsed, time_limit
	integer :: hr_curr, hr_prev


! MPIzation related, see below:
	character*50 :: myparfname, fnamesuffix, outputfname, rndmstr

	end module vrbls




!===============================================
!=== Main block :)
!===============================================
	program MAIN
	use vrbls; use det_n2
	use mumbers; use tree
	use bstats
	implicit none 

	logical :: acpt                   ! accepted update flag
	integer :: ddd, i, ans1, ans2, ans3
	real*8  :: r, dp, det , dt
	character*6 :: cvoid, version_tag


! ------- pi pi pi pi ------------
	pi = 4.d0*atan(1.d0)

!--------- Set the clock --------
	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	time_prev=0.d0; time_elapsed = 0.d0; hr_prev = 0


!--------- Get the par-file name from the command line --------------
        call  GETARG(1 , myparfname)

        ! determine the suffix for output files
        fnamesuffix = myparfname(4:)

        outputfname = 'out'//trim(fnamesuffix)   ! send output here


!------- Reading parameter file ------------------

        open(OUT_UNIT,file=outputfname,position='append')      ! to be open until the MC starts
        write(OUT_UNIT,*) ' reading ', myparfname,' ...'

        open( 1, file=trim(myparfname) )
      
      read(1, *)version_tag
      if(version_tag /= "__v2__")then
	 print*,'parameters: expected __v2__, got ', version_tag, '. STOP.'
	 call mystop
      endif
      read(1,*) ddd        ! Global dimension
	if(ddd/=d)then; print*,'d/=3'; call mystop; endif

      Nsite=1
      do i=1,d 
         read(1,*) N(i)  ! Number of sites per given dimension
         Nsite=Nsite*N(i)
	   N2(i)=N(i)/2
      end do
      
	read(1,*) ans1  ! 0 if new configuration, 1 if old one
	read(1,*) ans2  ! 0 if new statistics,    1 if old one
	read(1,*) ans3  ! 0 if new rndm() seed, 1 if read one
	
	read(1,*) mu            ! Chemical potential
	read(1,*) disorder_amp  ! Disorder amplitudean
	read(1,*) U,U_in        ! - U, interaction, & initial for thermalization
	read(1,*) beta          ! Inverse temperature
	    bun = beta*U*Nsite  ! Shorthand
	    step_r = bun*bun    ! step for determinant recalculation
		
	read(1,*) eta          ! G- vs. Z-sector weight
	read(1,*) mtau         ! # of points in \tau for tabulation of GFs
	   bmt = beta/mtau; bmt1 = 1.d0/bmt                   ! shorthands
	read(1,*) tolerance    ! det recalculation tolerance
	read(1,*) rx_ca, rt_ca ! cre/ann x- and \tau- radii
	    if(2*rt_ca>beta)then; print*,'rt_ca > beta /2'; call mystop; 
		endif
	read(1,*) rx_le, rt_le ! cre/ann x- and \tau- radii
	    if(2*rt_le>beta)then; print*,'rt_le > beta /2'; call mystop; 
		endif 
	read(1,*) nt_same              ! number of points for tabulation of weights for leaps
	    ddt_same = rt_le/nt_same   ! shorthands
	    udt=U*ddt_same
      read(1,*) n_sw                 ! number of sweeps for thermolization
      read(1,*) step_p               ! step for printing
      read(1,*) step_w               ! step for writing to disk
      read(1,*) step_m               ! step for measuring
	read(1,*) time_limit           ! time limit [hrs]
 	    time_limit = time_limit*3600   ! seconds      
	read(1,*) cvoid                ! horizontal separator :)
! Update manager probabilities: 
      read(1,*) prob(1); prob(2)=prob(1)+prob(1)                ! add/drop
      read(1,*) dp; prob(3)=prob(2)+dp;  prob(4)=prob(3)+dp     ! cre/ann
	read(1,*) dp; prob(5)=prob(4)+dp;  prob(6)=prob(5)+dp     ! leap_add/drop
	read(1,*) dp; prob(7)=prob(6)+dp                          ! hop
	if(abs(prob(7)-1.d0)>1.0d-10)then
	   print*,'Update manager probs. normalization problem...'
	   call mystop
	endif

!--- Random number generator:
!
!  Number of seeds provided MUST be equal to the groupsize (one seed per process)
!  Otherwise the program crashes with mysterious error message
!
!
       if(ans3==0)then     ! init rndm() state

          read(1,*)cvoid
          read(1,*) r_ij,r_kl
          call init_rndm(r_ij,r_kl)

	else                ! read rndm() state

	  rndmstr='rndm'//trim(fnamesuffix)//'.dat'
	  call read_rndm(rndmstr)

	endif



      close(1)

	write(OUT_UNIT,*)'....done'



!---------- DONE reading, allocating memory:

!--- Lattice
	allocate (ass(dd,Nsite),back(dd),x(1:d,1:Nsite))
      CALL ASSA

!--- Green functions
	allocate( GR_DAT(0:mtau+1,1:Nsite,1:Nsite), GRD_DAT(0:mtau+1,1:Nsite,1:Nsite) )

!--- site-dependent chemical potential
	allocate( musite(1:Nsite) )
	allocate( g0000(1:Nsite) )

!--- Configuration
	allocate( nkink(1:Nsite), kname(nkink_m,1:Nsite) )  ! site data
	allocate( ksite(-2:nmnm_max),ktau(-2:nmnm_max) )    ! global lists

	allocate( row(-2:nmnm_max),clmn(-2:nmnm_max) )      ! c/r <-> name links
	
	lda=128; allocate(nm_row(lda),nm_clmn(lda))         ! Allocate 128*128 matrix first. 
														! more memory will be allocated later if necessary	
	allocate(matr(1:lda,1:lda),m_u(1:lda),m_v(1:lda))   ! TheMatrix
	allocate(m_w(1:lda),m_z(1:lda))

	allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )
 

! Name Manager is initialized in either init_cnf or rd_cnf

! Site lists that contain sites within cubes with given edge lengths around given site
! it proves convenient to have TWO sets of lists: one for leap-type updates, one for cre/ann updates 
	allocate( uzli(1:Nsite),yes(1:Nsite) )
	n_ca=cube(1,rx_ca,uzli); allocate(uzli_cube_ca(1:n_ca,1:Nsite))
	 call tab_cube(rx_ca,n_ca,uzli_cube_ca)
	n_le=cube(1,rx_le,uzli); allocate(uzli_cube_le(1:n_le,1:Nsite))
	 call tab_cube(rx_le,n_le,uzli_cube_le)
	deallocate(uzli,yes)

! Working arrays for the closest neighbor selection in leap_drop.
	allocate( name_best(1:n_le),dist_best(1:n_le) )
	allocate( distance(1:nkink_m) )

!------------ DONE with memory allocation, proceed to tabulations

!--- Site-dependent chemical potential
	call RD_DISORD         ! read the replica


!--- Tabulate Green functions
	write(OUT_UNIT,*)'tabulating GFs...'  
	Ntab=maxval(N(:)/2+1)                          ! the longest jump w/PBC is N/2
	call TABULATE 
	write(OUT_UNIT,*)'...done'

	
!--- Tabulate the weights for leap_add/drop updates
	call tab_leaps

!--- Structure output
      nm_av=0; nm_min=i_hu; nm_max=0
      nmnm_distr=0; det_distr=0.d0; 


! 'how much on your watch?' a.k.a initialize the clock
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! take care of the runs that span across midnight 


	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' init. done: ',dt,' sec'	
	
!--- close output
	close(OUT_UNIT)


!----------- DONE with initializations



!=========== Read / init conf & thermalize =========

      if (ans1 == 1) then
	   call RD_CNF         ! read the condiguration
	   call therm2         ! thermalize if requested (via worm scheme)
      else
         call INIT_CNF       ! init the configuration 
	   call therm1         ! thermalize (via diagonal scheme)
      end if

! nullify counters
	i_p=0.d0; i_w=0.d0; i_m=0.d0; step=0.d0; 
      c_a_v=0.d0; a_a_v=0.d0; c_d_v=0.d0; a_d_v=0.d0
      c_cre=0.d0; a_cre=0.d0; c_ann=0.d0; a_ann=0.d0
	c_l_a=0.d0; a_l_a=0.d0; c_l_d=0.d0; a_l_d=0.d0
	c_r=0.d0; a_r=0.d0
	recalc_rate=0.d0

! Construct TheMatrix
	det=recalc_matrix(pm,lda,matr)


!============== Read / init  statistics ==============
      if (ans2 == 1) then
         call RD_STAT
      else
         call INIT_STAT
      end if


!==================== main MC loop  ===============
	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'  '
	write(OUT_UNIT,*)' Start MC loop'
	close(OUT_UNIT)

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_m=i_m+un1; i_r=i_r+1.d0


!------------------------ update --------------------------
	   r = rndm()
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call create(acpt,det)
	   else if(r<prob(4))then;    call annihilate(acpt,det)
	   else if(r<prob(5))then;    call leap_add(acpt,det)
	   else if(r<prob(6))then;    call leap_drop(acpt,det)
	   else if(r<=prob(7))then;   call hop(acpt,det)
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

! det. distribution
            i = floor(log10(abs(det)))
	      if(abs(i)<det_distr_max) det_distr(i) = det_distr(i)+1.d0

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-----------------------------------------------------


         nm_av=nm_av+(un1*nmnm)/step_p        
         IF(nmnm>nm_max)nm_max=nmnm
         IF(nmnm<nm_min)nm_min=nmnm

	   if(nmnm<nmnm_sup) nmnm_distr(nmnm) = nmnm_distr(nmnm) + 1.d0


         if (i_m == step_m)  then; i_m=nul; call MEASURE; end if
         if (i_p == step_p)  then; i_p=nul; call PRNT;
	       end if  
         if (i_w == step_w)  then;  i_w=nul; call WRT;     
	   end if 

	enddo


	contains  

!*************************  subroutines below ONLY ************************

!-----------------
!--- Add a kink
!-----------------
	subroutine add(acpt,det)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio
	integer :: name, site,j,nk, xnew(d), vova
	real*8  :: tnew, ti, tm, tv

	c_a_v = c_a_v + un1; acpt=.false.
	prevstatus=status; status='_add__' 

	if(pm==lda)then; call resize_matrix(-1); endif

!------- check
	if(nmnm>=nmnm_max)then;
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif
!---------- end check

	site=Nsite*rndm()+1.d0
	if(site>Nsite)site=Nsite    ! select a site to play on
	xnew(:) = x(:,site) 
	tnew = beta*rndm()                    ! select tau		

!------------- determinant
	if(pm==0)then; det=g0000(site)  !ratio = g0000**2
	else

! ira-masha 
      if(present)then
          tm=ktau(masha)
          m_u(row(masha)) = GREENFUN(site, tnew, ksite(masha), tm)
          ti=ktau(ira)
          m_v(clmn(ira)) = GREENFUN(ksite(ira), ti, site, tnew)
      endif

! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j); 
	   tv=ktau(vova)
	   m_v(clmn(vova)) = GREENFUN(ksite(vova), tv, site, tnew)
	   m_u(row(vova))  = GREENFUN(site, tnew, ksite(vova), tv)
	enddo
	m_s = g0000(site)
	det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*bun/(nmnm+1)


! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 

!	call CheckGetName;    ! this slows things down a lot 

	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; ktau(name)=tnew                ! global list

	if(pm==0)then                                    ! update TheMatrix
	   pm=1; matr(pm,pm)=1.d0/g0000(site)
	else
	   call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	a_a_v = a_a_v + un1

	end subroutine add



!-----------------
!--- Drop a kink
!-----------------
	subroutine drop(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio 
	integer :: site,j,nk, name, num,r,c, vova

	c_d_v = c_d_v + un1;   acpt=.false.
	prevstatus=status; status='_drop_'

	if(nmnm==0)return    ! nothing to drop yet :)
 
	num=nmnm*rndm()+1.d0                           ! play a kink to be dropped
	if(num>nmnm)num=nmnm; name=namelist(num); 
	site=ksite(name); nk = nkink(site)
	
!---------- determinant
	if(pm==1)then; det=1.d0/matr(pm,pm)
	else
	    r = row(name); c=clmn(name)
	    det = det_m1(pm,lda,matr,r,c)
	endif
!----------------------

	ratio = det*det*bun/nmnm
	ratio = un1/ratio

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
	if(pm==1)then; pm=0; nkink(site)=0
	else

	   do j=1,nk; if(kname(j,site)==name)exit   ! find name
	   enddo

	   kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site

	   vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
	   vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	   call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix

	endif

!	call CheckDropName(name);
	call DropName(name)

	a_d_v = a_d_v + un1

	end subroutine drop


!------------------
!--- Create ira & masha
!------------------
      subroutine create(acpt,det)
	logical :: acpt
	real*8  :: det

      real*8 :: ti,tm,tv,ratio, gim
      integer :: si,sm, sv, xi(d),xm(d) ,j, vova, jm

	prevstatus=status; status='create'; acpt=.false.

      if(present)return         ! one masha at a time, please

	if(pm==lda)then; call resize_matrix(-1); endif       ! check whether there is enough memory allocated

	c_cre=c_cre+1.d0

! seed ira uniformly over the volume & beta
	si=Nsite*rndm()+1.d0; if(si>Nsite)si=Nsite; 
	xi=x(:,si); ti=beta*rndm()

! seed masha within the cube around ira
	jm=n_ca*rndm()+1.d0
	if(jm>n_ca)jm=n_ca
	sm=uzli_cube_ca(jm,si); xm=x(:,sm)
	
	tm=ti + (2.d0*rndm()-1.d0)*rt_ca
	if(tm<0.d0)tm=tm+beta; if(tm>beta)tm=tm-beta   ! \beta-periodic
	
!----------- determinant
      gim=GREENFUN(si, ti, sm, tm)

      if(nmnm==0)then; det=gim
      else

          do j=1,nmnm; vova=namelist(j);
             sv=ksite(vova); tv=ktau(vova)
             m_u(row(vova))  = GREENFUN(si, ti, sv, tv)   	
             m_v(clmn(vova)) = GREENFUN(sv, tv, sm, tm)
          enddo
          m_s = gim
          det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

      endif  ! nmnm==0
!---------------------------

      ratio = det*det*eta             ! Note that \eta has been implicitely rescaled to absorb
	                                ! the macroscopic factos, c.f. measure(), where the reascaling 
	                                ! is undone for the statistics

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif
      
	acpt=.true.

! accepted.. update configuration

	present=.true.
	ksite(ira)=si; ktau(ira)=ti
	ksite(masha)=sm; ktau(masha)=tm

	if(pm==0)then
	   pm=1; matr(pm,pm)=1.d0/gim
      else
         call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)         ! TheMatrix
      endif

      clmn(ira) = pm; nm_clmn(pm) = ira
      row(masha)  = pm; nm_row(pm)  = masha

      a_cre=a_cre + 1.d0


	end subroutine create



!------------------
!--- Annihilate ira & masha
!------------------
      subroutine annihilate(acpt,det)
	logical :: acpt
	real*8  :: det

!
! This update is complementary to create() above
!

      real*8 :: ratio, dt
      integer :: vova, r,c,j, sm
	logical :: within

      prevstatus=status; status='annih_'; acpt=.false.

	if(.not.present)return   !  nothing to annihilate

!--------------- detailed balance: only erase if ira & masha are close enough 
	within=.false.
	sm=ksite(masha)
	do j=1,n_ca
	   if(sm==uzli_cube_ca(j,ksite(ira)))then; within=.true.; exit
	   endif
	enddo
	if(.not.within)return

	dt = abs(ktau(ira)-ktau(masha))
	if(min(dt,beta-dt)>rt_ca)return     ! min() here is the short way of handling PBC

!---------------------------------

	c_ann=c_ann+1.d0

!---------- determinant
      if(pm==1)then; det=1.d0/matr(pm,pm)
      else
           r = row(masha); c=clmn(ira)
           det = det_m1(pm,lda,matr,r,c)
      endif
!----------------------

	ratio=det*det*eta               ! Note the rescaling of \eta, see create() above
	ratio=1.d0/ratio

! Metropolis
	if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! accepted.. update configuration

      present=.false.

      if(pm==1)then; pm=0
      else
         vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
         vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
         call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix
      endif

      a_ann=a_ann+1.d0


	end subroutine annihilate



!----------------------------
!--- Leap masha & add a kink
!----------------------------
	subroutine leap_add(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is the heart of the method: worm-type update that creates an extra kink
!
	real*8  :: ratio, tnew, tm, tv, best, sign
	integer :: site,sm, sv, i,j,nk, xnew(d), xm(d), xv(d)
	integer :: vova, name, site1,ev, nk1

	prevstatus=status; status='leap_a'; acpt=.false.

	if(.not.present)return

	if(pm==lda)then; call resize_matrix(-1); endif

!------- check --------------------------------------
	if(nmnm>=nmnm_max)then;                 ! too many kinks
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'leap_add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif

	sm=ksite(masha)                         ! too many kinks on the particular site
	if(nkink(sm)==nkink_m)then
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'leap_add: nk(masha)= ',nkink(sm), step
	   close(OUT_UNIT)
	   call mystop
	endif
!---------- end check ----------------------------------

	c_l_a = c_l_a + 1.d0

	sm=ksite(masha); xm=x(:,sm); tm=ktau(masha)

!----------------------------- select where to leap
	ev=EVENT(rndm())                                ! EVENT() selects where to add new masha 
													! according to the weights w_le(...);
													! the function itself is defined in the 
													! 'tree' module, file event_mod.f
! PBC: coords
	xnew=xm+dx_le(:,ev)                             ! select coordinates according to the event number
	do i=1,d
	   if( xnew(i) < 1 )    xnew(i)=xnew(i)+N(i)
	   if( xnew(i) > N(i) ) xnew(i)=xnew(i)-N(i)
	enddo
	site = xnew(1) + N(1)*(xnew(2)-1)+N(1)*N(2)*(xnew(3)-1)

! PBC: \tau-s
	sign=1.d0; !if(rndm()>0.5d0)sign=-1.d0
	tnew=tm + sign* (dt_le(ev) + ddt_same*(rndm()-0.5d0) )        ! unidirectional leaps 
	if(tnew<0)tnew=tnew+beta; if(tnew>beta)tnew=tnew-beta         
!------------------------------------------------

!------------------------ check detailed balance
!
! The search for the closest kink is made in two steps: first
!  the best candidate is looked for on each site within the cube,
!  and then, the very best kink is selected among the site-best ones.
!  This way the sites may be processed in parallel.


! pick the best candidate on each site
	do i=1,n_le; site1=uzli_cube_le(i,site); xv = x(:,site1)

	   dist_best(i)=1.d200
	   
	   nk1 = nkink(site1)
	   if(nk1>0)then

	     do j=1,nk1; tv=ktau(kname(j,site1))          ! get the list of candidates on the site
	         distance(j) = GOK(xnew,tnew,xv,tv)       ! GOK() stands for goodness-of-kink
	     enddo
	     dist_best(i)= minval(distance(1:nk1))        ! best candidate on this site

	   endif   ! nk1>0

	enddo

! 2nd round: pick the very best of the best
	best=minval(dist_best) 

	if(best<GOK(xnew,tnew,xm,tm))return               ! detailed balance cannot be satisfied
													  ! the update can't be completed
!-------------------------------------------
 


!------------- determinant
	do j=1,pm; vova=nm_row(j)               ! fill a column
	   sv=ksite(vova); tv=ktau(vova)
	   m_u(j)  = GREENFUN(sm, tm, sv, tv)
	enddo

	do j=1,pm; vova=nm_clmn(j)              ! fill a row
	   sv = ksite(vova); tv=ktau(vova)
	   m_v(j) = GREENFUN(sv, tv, site, tnew)
	enddo

	m_s = GREENFUN(sm, tm, site, tnew)          ! a corner & det
	det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself
!---------------------------

	ratio = det*det*udt/w_le(ev) 

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 

! 	call CheckGetName; 
	call GetName(name)
	nk=nkink(sm); 
	nkink(sm) = nk+1; kname(nk+1,sm) = name           ! on-site info
	ksite(name)=sm;   ktau(name)=ktau(masha)          ! global list

	call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)          ! TheMatrix

	row(name)  = row(masha); nm_row(row(masha))  = name;   ! TheMatrix <-> names links
	clmn(name) = pm; nm_clmn(pm) = name

! update masha
	row(masha)=pm; nm_row(pm)=masha
	ksite(masha)=site; ktau(masha)=tnew

	a_l_a = a_l_a + 1.d0

	end subroutine leap_add



!-----------------
!--- Leap masha & drop a kink
!-----------------
	subroutine leap_drop(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to leap_add() above
!
	real*8 :: ratio, tau, tm,t1,tv,www, dovesok
	integer :: site,i,j,nk, name, r,c, vova, sm
	integer ::  site1, nk1, best_i(1), dx(d), xm(d), xv(d)

	prevstatus=status; status='leap_d'; acpt=.false.; 

	if(.not.present)return   
	if(nmnm==0)return          ! nothing to drop yet :)

	c_l_d = c_l_d + 1.d0

	sm = ksite(masha); xm=x(:,sm); tm=ktau(masha)

!------------------- find the closest kink [which will be deleted]
!
! The selection proceeds in two rounds, see leap_add() above 
!

! pick the best candidate on each site
	do i=1,n_le; site1=uzli_cube_le(i,sm); xv=x(:,site1)

	   dist_best(i)=1.d200; name_best(i) = -101     ! -101 stands for 'no such kink'
	   
	   nk1 = nkink(site1)
	   if(nk1>0)then

	     do j=1,nk1; tv = ktau(kname(j,site1))   ! get the list of candidates on the site
	         distance(j) = GOK(xm,tm,xv,tv)
	     enddo
	     
		best_i= minloc(distance(1:nk1))        ! best candidate on this site
	    dist_best(i) = distance(best_i(1))
	    name_best(i) = kname(best_i(1),site1)

	   endif   ! nk1>0

	enddo

! 2nd round: pick the very best of the best
	best_i=minloc(dist_best)
	name = name_best(best_i(1)); if(name==-101)return       ! can't do anything if there is nobody around

	site = ksite(name); tau = ktau(name) 
!------------------------------------
 

!-------- detailed balance requires using the proper weight, see leap_add() and tab_leaps() for details
	t1 = abs(tm-tau); t1 = min(t1,beta-t1)
	j = floor(t1/ddt_same)
	t1 = -ddt_same*(j+0.5d0)
 
 	dx = x(:,site)-x(:,sm)
	do i=1,d;  if(dx(i)>N2(i))dx(i)=dx(i)-N(i)   ! PBC
	enddo
	dovesok = ( dx(1)+dx(2)*N(1)+dx(3)*N(1)*N(2) ) * 1.d-14      ! see tab_leaps()

	www = GREENFUN(site, 0.0d0, sm, t1)**2 + dovesok
	www = www / w_norm
!----------------------------------------
 

!---------- determinant
	r = row(masha); c=clmn(name)
	det = det_m1(pm,lda,matr,r,c)
!------------------------

	ratio =  det*det*udt/www 
	ratio = 1.d0/ratio


!------------ Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

!---------- update configuration 

! on-site info
	nk=nkink(site)
	do j=1,nk; if(kname(j,site)==name)exit   ! find name
	enddo

	if(nmnm==1)then; nkink(site)=0            ! erase the kink from the site list
	else;  kname(j,site) = kname(nk,site); nkink(site) = nk-1 
	endif

! TheMatrix
	vova = nm_row(pm); row(vova) = r; nm_row(r) = vova    
	vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	call inv_m1(pm,lda,matr,r,c)                    

! masha
	ksite(masha)=site; ktau(masha)=tau
	row(masha)=row(name); nm_row(row(masha))=masha

! erase the kink from the global list
!	call CheckDropName(name); 
	call DropName(name)

	a_l_d = a_l_d + 1.d0

	if(prevstatus=='leap_a')i_backsc = i_backsc + 1.d0      ! this is a 'time wasted' counter: since dropping
															! is deterministic, the pair of add() 
															! immediately followed by drop() does
															! not change a configuration.
	end subroutine leap_drop


!------------------
!--- Hop masha
!------------------
	subroutine hop(acpt,det)
	logical :: acpt
	real*8  :: det

	real*8 :: ratio, tm, tv, tnew  
	integer :: i,j,sm, sv, sm_new,xm(d),xnew(d), r, vova

	prevstatus=status; status='__hop_'; acpt=.false.

	if(.not.present)return

	c_r = c_r + 1.d0

! propose the hop; uniform in <nn> & same \tau
	sm=ksite(masha); tm=ktau(masha); xm=x(:,sm)

	i=dd*rndm()+1; if(i>dd)i=dd
	sm_new = ass(i,sm); xnew=x(:,sm_new)

	tnew = tm + (beta/5.d0)*(2.d0*rndm()-1.d0)
	if(tnew<0)tnew=beta+tnew; if(tnew>beta)tnew=tnew-beta 

!--------- the row & det
	do j=1,pm; vova=nm_clmn(j)
	   sv = ksite(vova); tv=ktau(vova)
	   m_v(j) = GREENFUN(sv,tv,sm_new,tnew) - GREENFUN(sv,tv,sm,tm)
	enddo

	r=row(masha); det = det_r(pm,lda,matr,r,m_v)
!---------------------------

	ratio =  det*det

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 

	call inv_r(pm,lda,matr,r,det,m_v,m_w,m_z)   ! TheMatrix

	ksite(masha)=sm_new; ktau(masha)=tnew       ! masha

	a_r = a_r + 1.d0

	end subroutine hop



!------------------
!--- Measurements
!------------------
	subroutine measure

	integer :: j
	real*8  :: eta1

	if(present)then       !   GF sector

	    j = floor((ktau(masha)-ktau(ira))*im_t_max/beta)    
	    im_t(j) = im_t(j)+1.d0                       ! t(ira-masha) histogram

		im = im + 1.d0                               ! integrated correlator

	else                  !   Z sector

	i_b = i_b + 1 ; 	Z = Z + un1
	PE = PE + 1.d0*nmnm/beta               ! PE
	KE = KE + diag_KE()                    ! KE
	ndens = ndens + diag_dens()            ! number density 
!	call dance_dance                       ! dance-dance correlators

	endif  ! present


!-------------- is a current block full?
	if( i_b == Z_b )then    ! wrap it up
	    b_n = b_n + 1; i_b=0
		PE_stat(b_n) = PE; PE = 0.d0
		KE_stat(b_n) = KE; KE = 0.d0
		ndens_stat(b_n) = ndens; ndens = 0.d0


! normalize [see create()] & rescale via *L^{gf_eta}/(beta*Nsite)**2
		eta1 = n_ca*2.d0*rt_ca/eta
	    eta1 = eta1 * N(1)**(gf_eta)/beta/Nsite
		im_stat(b_n) = im*eta1; im = 0.d0      ! normalize via eta1*Z

		if(b_n == b_n_max)then                  
	        call collate( im_stat, b_n)
		    call collate( PE_stat, b_n )
			call collate( KE_stat, b_n )
			call collate( ndens_stat, b_n )
			b_n = b_n/2; Z_b = Z_b * 2.d0
		endif


	endif



	end subroutine measure


!--------------------------------
!--- measure density: insert an extra kink @ random place
!--------------------------------
	real*8 function diag_dens()

	real*8 :: det,ti, tn
	integer :: site, si, xn(d), j, vova 

! select where to insert a kink
	site=Nsite*rndm()+1.d0; if(site>Nsite)site=Nsite
	xn(:)=x(:,site); tn=beta*rndm()
	
!------------- determinant
	if(nmnm==0)then; det=g0000(site)
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     si = ksite(vova); ti=ktau(vova)
	     m_u(j)  = GREENFUN(site, tn, si, ti)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     si = ksite(vova); ti=ktau(vova)
	     m_v(j) = GREENFUN(si, ti, site, tn)
	  enddo

	  m_s = g0000(site)

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! nmnm==0

! estimator
	diag_dens = det*2.d0   ! a factor of 2 reflects the spin summation

	end function diag_dens




!-----------------------------------------------
!--- measure kinetic energy from the Z-sector
!-----------------------------------------------
	real*8 function diag_KE()

	integer :: site1, site2, sv, x1(d),x2(d), vova, j, i
	real*8  :: t, det,tv
!
!  This estmator is for the tight-binding dispersion ONLY (n.n. sites)
!
	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite; 

	i=dd*rndm()+1; if(i>dd)i=dd; site2 = ass(i,site1)
		
	x1=x(:,site1); x2=x(:,site2)
	t=beta*rndm()
	

!------------- determinant
	m_s = GREENFUN(site2, t, site1, t)        ! a corner 

	if(pm==0) then; det = m_s        
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     sv = ksite(vova); tv=ktau(vova)
	     m_u(j)  = GREENFUN(site2, t, sv, tv)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     sv = ksite(vova); tv=ktau(vova)
	     m_v(j) = GREENFUN(sv, tv, site1, t)
	  enddo

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif

! return value
	diag_KE = dd*Nsite*det  ! d*Nsite is the # of bonds,       
							! a factor of 2 is due to the spin summation


	end function diag_KE



!--------------------------------
!--- measure density-density correlators : FIXME g0000
!--------------------------------
	subroutine dance_dance()
!
! This is a very time-consuming version since it calculates 
! the contributions to all L/2 distances
!
	real*8 :: det1,det2, t1, tv 
	integer :: site1,site2, sv, x1(d),x2(d),j,vova, dir,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	print*, "DANCE DANCE g0000 FIXME"
	stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! play two extra half-kinks
	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite
	x1(:)=x(:,site1); t1=beta*rndm()

	          do i=0,N2(1)

	site2=site1;  dir=1 
	     do j=0,i-1; site2= ass(dir,site2); enddo; 
	x2=x(:,site2)



!--- determinantz
!!	m_s = g0000

!!	m_s2(1,1)=g0000; m_s2(1,2)=GREENFUN(site2, t1, site1, t1)
!!	m_s2(2,2)=g0000; m_s2(2,1)=GREENFUN(site1, t1, site2, t1)


	if(pm==0)then; det1=g0000(site1)**2; det2 = det1 - m_s2(1,2)*m_s2(2,1)
	else

! 1st for up-down
	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     sv = ksite(vova); tv=ktau(vova)
	     m_u(j)  = GREENFUN(site1, t1, sv, tv)
	  enddo
	  m_u2(1:pm,1)=m_u(1:pm)

	  do j=1,pm; vova=nm_clmn(j)              ! fill a row
	     sv=ksite(vova); tv=ktau(vova)
	     m_v(j) = GREENFUN(sv, tv, site1, t1)
	  enddo
	  m_v2(1,1:pm)=m_v(1:pm)

	  det1 = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)     

! 2nd for up-down
	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     sv = ksite(vova); tv=ktau(vova)
	     m_u(j)  = GREENFUN(site2, t1, sv, tv)
	  enddo
	  m_u2(1:pm,2)=m_u(1:pm)

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     sv = ksite(vova); tv=ktau(vova)
	     m_v(j) = GREENFUN(sv, tv, site2, t1)
	  enddo
	  m_v2(2,1:pm)=m_v(1:pm)
	  
	  det2 = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)     


! up-down
	  det1 = det1*det2


! up-up
	  det2 = det_p2(pm,lda,matr,m_u2,m_v2,m_s2,m_c2)

	endif


	g_ud(i) = g_ud(i) + det1 
	g_uu(i) = g_uu(i) + det2 

	                    enddo   ! i



	end subroutine dance_dance


!================  DONE with measurements; various service functions below 

!------------------------------------
!--- Print out and check the runtime
!------------------------------------
	subroutine prnt
	implicit none
	integer :: i,j
	real*8 :: xxx, yyy, dt 

	real*8 :: PE_av, PE_err, KE_av, KE_err, im_av, im_err
	real*8 :: dens_av, dens_err

	logical :: lastone, dens_conv, im_conv


! 'how much on your watch?'
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! across the midnight

	time_elapsed = time_elapsed + dt

	lastone=.false.
	if( time_elapsed > time_limit )  lastone=.true.   ! time to wrap up
!-------------------------------------------

	open(1,file=outputfname,position='append')

	write(1,*)''
      write(1,*)'-------------------------', time_elapsed/3600,' hrs'
	
	if(i_t<step_t)then
	  write(1,*)' thermalization: ',1.d2*i_t/step_t,' % done'
	  write(1,*)' current U = ', bun/beta/Nsite
	  write(1,*)'  '
	endif

! is there enough statistics?
	if(b_n>3)then

      do i=1,d
      write(1,fmt=700)i,N(i)
	end do
 700  format (8x,'N(',i1,') =', i4)

      write(1,*)' '
	write(1,*)' 1/T  = ', beta, '  -U  = ', bun/beta/Nsite
	write(1,*)' \mu  = ', mu, ' av nnn = ', sum(g0000)/Nsite
	write(1, *)' amp = ', disorder_amp, '  replica id = ', replica_id
	write(1,*)'  '

      write(1,*)' MC step (mln) =', step/1.d6,'  Z(mln) = ',Z/1.d6
	write(1,*)' Z_b(mln) = ',Z_b/1d6,' b_n = ', b_n
      write(1,*)' '


	write(1,fmt=771)nmnm,nm_max,nm_min,nm_av
 771  format(' nmnm-> ',I5,' now [ max = ',I5,' min = ',I5,  ' average = ',G11.5,' ]')


!--- pot. energy -----------------------------------
	call bSTAT(PE_stat(1:b_n),b_n, 1.d0*Z_b,PE_av,PE_err)

	write(1,*)'  '
      write(1,fmt=701) PE_av, PE_err
 701  format(8x,'PE =',g12.5,4x,' +/- ',g10.3)

      
!--- kin. energy -----------------------------------
	call bSTAT(KE_stat(1:b_n),b_n, 1.d0*Z_b,KE_av,KE_err)

	write(1,*)'  '
      write(1,fmt=711) KE_av, KE_err
 711  format(8x,'KE =',g12.5,4x,' +/- ',g10.3)
      

!--- number density ---------------------------------
	call mrg_conv(ndens_stat(1:b_n), b_n, Z_b, dens_av, dens_err, dens_conv)
	write(1,*)'  '
    write(1, fmt=784) 'dens ', dens_av, dens_err, bool_to_str(dens_conv)

    call mrg(ndens_stat, b_n, Z_b, OU=1)

!--- integrated correlator ---------------------------
	call mrg_conv(im_stat(1:b_n), b_n, Z_b, im_av, im_err, im_conv)
	write(1,*)'  '
	write(1, fmt=784) 'g_im(w=0, k=0)', im_av, im_err, bool_to_str(im_conv)

    call mrg(im_stat, b_n, Z_b, OU=1)


 784    format(8x, A, ' = ', g12.5,2x, ' +/- ', g12.5,8x, '   ', A)


! write the "results" line
    open(2, file='res'//trim(fnamesuffix)//'.dat')
    write(2, *) replica_id, dens_av, dens_err, bool_to_char(dens_conv), im_av, im_err, bool_to_char(im_conv), Z/1d6
    close(2)

!-----------------------------------------------------
	write(1,*)
	write(1,*)' address/accept:'
	write(1,*)' add/drop ',a_a_v/(c_a_v+.001),' / ',a_d_v/(c_d_v+.001)
	write(1,*)' cre/ann  ',a_cre/(c_cre+.001),' / ',a_ann/(c_ann+.001)
	write(1,*)' leap a/d ',a_l_a/(c_l_a+.001),' / ',a_l_d/(c_l_d+.001)
	write(1,*)' hop        ',a_r/(c_r+.001)
	write(1,*)' recalc / backsc ', recalc_rate/(step + 0.0001),' / ', i_backsc/(a_l_a + .0001)

!---- nullify counters
      nm_max=0; nm_min=i_hu; nm_av=nul

!=============================  writeouts: various service distributions


! write t(ira-masha) distribution
	if(im_t(0)>0)then;	yyy = PE_av/U/Nsite/im_t(0)
	else; yyy=1.d0
	endif

	open(2,file='ct3d'//trim(fnamesuffix)//'.dat')
	do j=-im_t_max,im_t_max-1
	  write(2,*)beta*(j+0.5d0)/im_t_max, im_t(j)*yyy
	enddo
  	write(2,*)beta, im_t(im_t_max)*yyy  ! the last point
	close(2)


! write kink number distrr
      xxx=sum(nmnm_distr(1:nmnm_sup))
      if(xxx==0.d0)xxx=1.d0

      open(2,file='nmnm'//trim(fnamesuffix)//'.dat')
	do j=0,nmnm_sup; write(2,*)j,nmnm_distr(j)/xxx
      enddo
	close(2)

! write det distr
	xxx=sum(det_distr(-det_distr_max:det_distr_max))
	if(xxx==0.d0)xxx=1.d0

	open(2,file='det'//trim(fnamesuffix)//'.dat')
	do j=-det_distr_max,det_distr_max
	   write(2,*)j, det_distr(j)/xxx
	enddo
	close(2)


!---------  uncomment this if you want to see the configuration -----------

! write kink visualization
!	open(2,file=trim(workdirpath)//'viz'//trim(fnamesuffix)//'.dat')
!	do i=1,nmnm; name=namelist(i)
!		write(2,*)x(:,ksite(name)),ktau(name)
!	enddo
!	close(2)
!
!	if(present)then
!
!	open(2,file=trim(workdirpath)//'viz_i'//trim(fnamesuffix)//'.dat')
!		write(2,*)x(:,ksite(ira)), ktau(ira)
!	close(2)
!
!	open(2,file=trim(workdirpath)//'viz_m'//trim(fnamesuffix)//'.dat')
!	    write(2,*)x(:,ksite(masha)),ktau(masha)
!	close(2)
!
!	endif
!---------------------------------------------------------------------------

	endif   ! b_n>3

	close(1)


! time to wrap up? --- write everything, release allocated memory and wrap up.
	if(lastone)then; 
	    open(1,file=outputfname,position='append')
		write(1,*)'Time to wrap up.'; 
		close(1)
	    call wrt
		call mystop
	endif


	end subroutine prnt
	


!-------------------
!--- Write data
!-------------------
	subroutine wrt
	integer :: j, name

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'write .....'

! configuration
	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    write(1,*)beta
          write(1,*)present
           if(present)then
              write(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
              write(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
           endif
	  write(1,*)nmnm,pm
	  write(1,*)nm_row(1:pm)
	  write(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	     name=namelist(j)
	     write(1,*)name, ksite(name),ktau(name), row(name), clmn(name)
	  enddo
	  write(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! rndm state
	rndmstr = 'rndm'//trim(fnamesuffix)//'.dat'
	call save_rndm(rndmstr)
!
! There is a switch in the parfile as to read 
! the generator state or to seed it anew.
!

! statistics
	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  write(1,*)d,N
	  write(1,*)beta, U, mu
	  write(1,*)step, Z
	  write(1,*)Z_b, b_n, i_b
	  write(1,*)PE
	  write(1,*)PE_stat(1:b_n)
	  write(1,*)KE
	  write(1,*)KE_stat(1:b_n)
	  write(1,*)ndens
	  write(1,*)ndens_stat(1:b_n)
	  write(1,*)im
	  write(1,*)im_stat(1:b_n)
	  write(1,*)'------'
	  write(1,*)g_uu(0:N2(1))
	  write(1,*)g_ud(0:N2(1))
	close(1)


	write(OUT_UNIT,*)'writing done!'
	close(OUT_UNIT)

	end subroutine wrt


!---------------------
!--- Read statistics
!---------------------
	subroutine rd_stat
	integer :: ddd
	real*8  :: dummy
	character*6 :: cvoid


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading stat .....'

	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  read(1,*)ddd,N
	      if(ddd/=d) call mystop
	  read(1,*)dummy, dummy, dummy
	  read(1,*)step, Z
	  read(1,*)Z_b, b_n, i_b
	  read(1,*)PE
	  read(1,*)PE_stat(1:b_n)
	  read(1,*)KE
	  read(1,*)KE_stat(1:b_n)
	  read(1,*)ndens
	  read(1,*)ndens_stat(1:b_n)
	  read(1,*)im
	  read(1,*)im_stat(1:b_n)
	  read(1,*)cvoid
	  read(1,*)g_uu(0:N2(1))
	  read(1,*)g_ud(0:N2(1))
	close(1)

	det_distr=0.d0; nmnm_distr=0.d0

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_stat


!---------------------
!--- Init statistics
!---------------------
	subroutine init_stat

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init stat .....'

	Z_b = ceiling(beta*U*Nsite) + 1    ! The block size for statistics 
	Z_b = Z_b / 8                      ! [something to start with]
	i_b = 0; b_n = 0

	Z=0.d0
	im = 0.d0; im_stat=0.d0
	PE = 0.d0; PE_stat=0.d0
	KE = 0.d0; KE_stat=0.d0
	ndens = 0.d0; ndens_stat=0.d0
	g_uu=0.d0; g_ud=0.d0
!---------------------------------

	im_t=0.d0; det_distr=0.d0; nmnm_distr=0.d0; 

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)


	end subroutine init_stat


!--------------------------
!--- Read configuration
!--------------------------
	subroutine rd_cnf
	integer :: j,name, site
	real*8  :: bbb, f

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading conf.....'
	
! configuration
	namelist(ira)=ira; namelist(masha)=masha
      do j=1,nmnm_max; namelist(j)=j; numnam(j)=j
      enddo; nmnm=0; pm=0

	nkink=0; kname=0

	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    read(1,*)bbb;        f = beta / bbb    
					! If the configuration if for some \beta_1, not equal to \beta
					! then all the times are just 'stretched' by a factor of f.
					! This feature is useful for re-thermalization of configurations
					! with close \beta-s
          read(1,*)present
           if(present)then
              read(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
              read(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
           endif
	  read(1,*)nmnm,pm
	     if(nmnm>nmnm_max)then
	        print*,'rd_cnf: nmnm= ',nmnm,' while nmnm_max= ',nmnm_max
			call mystop 
	     endif
!================== allocate enough memory to fit TheMatrix
	    deallocate(matr,m_u,m_v,m_w,m_z,nm_row,nm_clmn)    ! deallocate first
	    deallocate(m_u2,m_v2,m_s2,m_c2)
		lda=(int(pm/64)+1)*64+128
		allocate(nm_row(lda),nm_clmn(lda))
		allocate(matr(lda,lda))
		allocate(m_u(lda),m_v(lda),m_w(lda),m_z(lda))
	    allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )
!=========================================
	  read(1,*)nm_row(1:pm)
	  read(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	      read(1,*)name, site,ktau(name), row(name), clmn(name)
	      namelist(j)=name; numnam(name)=j
	      ktau(name) = ktau(name) * f          ! rescale \tau
	         ksite(name)=site
	         nkink(site)=nkink(site)+1
	         kname(nkink(site),site)=name
	  enddo
	  read(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! calculate TheMatrix
	bbb = recalc_matrix(pm,lda,matr)
	

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_cnf


!--------------------------
!--- Init configuration
!--------------------------
	subroutine init_cnf
	integer :: j

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init conf.....'

! configuration
	nkink=0; kname=0
	ksite=0; ktau=0.d0
	row=0; clmn=0; pm=0
	nm_row=0; nm_clmn=0

! name manager
	namelist(ira)=ira; namelist(masha)=masha	
      do j=1,nmnm_max
          namelist(j)=j; numnam(j)=j
      enddo; nmnm=0; 

! no ira, masha @ the beginning
      present=.false.

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine init_cnf



!----------------------------------------------
!---  Init the site-dependent chemical potential
!----------------------------------------------
!!	subroutine init_disord
!!	implicit none
!!	integer :: site
!!	
!!      DO NOT USE: ONLY USE PRE-GENERATED REPLICAS
!!
!!	do site=1, Nsite
!!	   musite(site) = mu + (2.d0*rndm() - 1.d0)*disorder_amp
!!	enddo
!!	
!!	end subroutine init_disord



!----------------------------------------------
!---  Read in the site-dependent chemical potential
!----------------------------------------------
	subroutine rd_disord
	implicit none
	integer :: L

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading the replica.....'
	
	open(1,file='disord'//trim(fnamesuffix)//'.dat')
	read(1, *) replica_id     ! np.random.seed the replica is generated with

	read(1, *) L
		if(L /= N(1))then
			print*, "Disorder replica: L = ", L, ' /= N(1) =', N(1);
			call mystop
		endif
	read(1, *) musite
	close(1)

!	print*, '-----------'
!	print*, musite
!	pause

	musite = musite*disorder_amp + mu

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)
	
	end subroutine rd_disord



!----------------------------------------------
!---  Green Function (former selector + spline interp.)
!----------------------------------------------
      real*8 function GREENFUN(site1,tau1,site2,tau2)
      implicit none
      integer :: site1,site2, sgn
      double precision :: tau, tau1, tau2, dt, gre

      integer :: nta  !, ntb
      double precision :: tta,ttb,ga,gb,c, gra,grb   !,p


! prepare \tau
      tau=tau1-tau2
      dt=tau; sgn=1

	if(tau < 1.d-14)then; dt=beta+tau; sgn=-1; endif
! Explanation: G(t=0) must be understood as G(t-> -0) = -G(t=\beta)
! A long way to accomplish this is below, commented out. A short way is above :).
!
!     if (abs(tau) < 1.d-14) then; dt = beta; sgn=-1
!     else if (tau > 0.d0) then; dt=tau; sgn=1
!	else; dt=beta+tau; sgn=-1
!	end if
!----------------------------------------



! prepare coords
!	j=1; nx = abs(x1(j)-x2(j)); nx = min(nx,N(j)-nx)
!	j=2; ny = abs(x1(j)-x2(j)); ny = min(ny,N(j)-ny)
!	j=3; nz = abs(x1(j)-x2(j)); nz = min(nz,N(j)-nz)

!----------------------------------- spline

	nta=dt*bmt1 !*p

      tta=dt-nta*bmt 
	ttb=tta - bmt     !dt-ntb*(beta/mtau) 
!cccccccccccccccccccccccccccccccccccccc
      
	ga=GR_DAT(nta,site1,site2)
	gb=GR_DAT(nta+1,site1,site2)

	gra=GRD_DAT(nta,site1,site2)
	grb=GRD_DAT(nta+1,site1,site2)

      c=(ga-gb)*bmt1

      gre=(c+gra)*ttb + (c+grb)*tta
      gre=gre*tta*ttb*bmt1 + gb*tta-ga*ttb
      gre=gre*bmt1


	GREENFUN = gre*sgn

      end function GREENFUN


!-------------------------------------
!     Tabulates Green function and its time derivate at positive tau.
!-------------------------------------
      subroutine TABULATE
      implicit none
	real*8, allocatable :: ham(:,:)
	integer :: site,site1,j
	real*8 :: factor, ww,ttt,term, gamma, expet(0:mtau)
	integer :: nt

	! lapack stuff
	character*1 :: jobz,uplo
	integer     :: ldh, lwork,info
	real*8, allocatable  :: work(:), eps(:)

	print*,'start with TABULATE'

! check the array sizes
!	if(Nsite>ntab_max)then
!	    print*,' ******************************'
!	    print*,' TABULATE: need to have ntab_max >= Nsite '
!	    print*,' currently Nsite = ',Nsite, 'and ntab_max = ', ntab_max
!	    call mystop
!	endif

!        do site=1,Nsite
!                print*, site, musite(site)
!        enddo
!        pause


! build the hamiltonian
	allocate(ham(1:Nsite,1:Nsite)) ; 
	ham=0.d0
	do site=1,Nsite
	   ham(site,site)=ham(site,site)-musite(site)
	   do j=1,dd; site1=ass(j,site); 
	      ham(site,site1)=ham(site,site1)-1.d0  
	   enddo
	enddo;	

!  compute eigenvalues; for LAPACK parameters and arguments, see
!  http://www.netlib.org/lapack/double/dsyev.f
!  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

	jobz='V'  ! compute eigenvectorz
	uplo='U'  ! use upper diag of ham(:,:) --- doesn't matter really
	ldh=Nsite
	lwork=12*Nsite
	allocate( work(lwork), eps(Nsite) )

! query the optimal workspace size
	call dsyev(jobz,uplo,Nsite,ham,ldh,eps,work,-1,info)
	lwork=work(1)
	deallocate(work); allocate(work(lwork))

! diagonalize
	call dsyev(jobz,uplo,Nsite,ham,ldh,eps,work,lwork,info)

	if(info/=0)then; print*,'*** dsyev returns info = ', info
			 print*,'*** check the TABULATE routine'
			 call mystop
	endif


!	print*,'optimal lwork =', work(1),lwork
!
!	print*,'eigenvalues:'
!	do j=1,Nsite
!		print*,eps(j)
!	enddo
!
!	print*; print*,'----------- eigenvector # '
!	do j=1,Nsite;
!		!write(1,*)j,ham(j,6)
!		print*,j,ham(j,1)
!	enddo
!	print*; print*; print*

!------------- have the spectrum, proceed to the GFs
	GR_DAT=0.d0; GRD_DAT=0.d0	

	do j=1,Nsite

	  gamma=-eps(j)*beta
	  gamma=exp(gamma)+1.d0
	  gamma=-1.d0/gamma

	  ww = exp(-eps(j)*bmt) ! bmt=beta/mtau
	  do nt=0,mtau; expet(nt)=ww**nt
	  enddo

          do site=1,Nsite ; do site1=1,Nsite
	    factor = ham(site,j)*ham(site1,j)
            do nt=0,mtau
               term = factor*expet(nt)*gamma    !/Nsite
               GR_DAT(nt,site,site1) = GR_DAT(nt,site,site1) + term
               GRD_DAT(nt,site,site1) = GRD_DAT(nt,site,site1) -eps(j)*term
            enddo

          enddo; enddo ! site, site1

	enddo   ! j: eigenvalues
!------------------------


! fill in fictitious nt=mtau+1, see GREENFUN for explanation
	GR_DAT(mtau+1,:,:)=0.d0; GRD_DAT(mtau+1,:,:)=0.d0

! fill g0000
	ttt=0.d0
	do site=1, Nsite
		g0000(site) = GREENFUN(site,ttt,site,ttt)
	enddo
	
!	do site=1, Nsite
!	print*, site, GREENFUN(site, ttt, site, ttt)
!	enddo
	
!	print*, "mean value = ", sum(g0000) / Nsite
!	pause
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!	site0=1; site1=site0
!	do j=1,N(1)
!	  print*,site0,site1,GR_DAT(100,site0,site1),GRD_DAT(100,site0,site1)
!	  site1=ass(1,site1); site1=ass(2,site1)
!	enddo

	print*,'done with TABULATE'

      end subroutine TABULATE



!-------------------------------
!--- Increase LDA (matrix storage): 
!      if lnew == -1 then just add 512
!-------------------------------
	subroutine resize_matrix(lnew)
	integer :: lnew

	integer :: lda_new 
	real*8, allocatable :: matr_new(:,:)
	real*8, allocatable :: nm_r_new(:), nm_c_new(:)

	if(lnew==-1)then; lda_new=lda+512
	else; lda_new = lnew
	endif

	allocate(matr_new(lda_new,lda_new))
	allocate(nm_r_new(lda_new),nm_c_new(lda_new))

! save TheMatrix as is
	matr_new(1:pm,1:pm)=matr(1:pm,1:pm)
	nm_r_new(1:pm)=nm_row(1:pm); nm_c_new(1:pm)=nm_clmn(1:pm)

! Resize
	deallocate(matr); deallocate(m_u,m_v,m_w,m_z,m_u2,m_v2,m_s2,m_c2)
	deallocate(nm_row,nm_clmn)

	lda=lda_new

	allocate(matr(lda,lda),m_u(lda),m_v(lda))
	allocate(m_w(lda),m_z(lda))

	allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )

	allocate(nm_row(lda),nm_clmn(lda))


! Restore
	matr(1:pm,1:pm)=matr_new(1:pm,1:pm)
	nm_row(1:pm)=nm_r_new(1:pm); nm_clmn(1:pm)=nm_c_new(1:pm)

! Cleanup
	deallocate(matr_new,nm_r_new,nm_c_new)


	end subroutine resize_matrix


!-------------------------------
!--- Recalculate the inverse from scratch (N**3 operations)
!-------------------------------
	real*8 function recalc_matrix(pm,lda,a)
	integer :: pm,lda     ! actual size & leading dimension of A
      real*8  :: a(lda,lda)

	integer :: i,j, vova, lesha, xi(d), xj(d), si, sj
	real*8  :: ti,tj, det_big

	if(pm<2)return   ! no use to recalculate

	recalc_rate = recalc_rate + 1.d0

! build the matrix
	do j=1,pm; do i=1,pm
	  vova=nm_row(i); lesha=nm_clmn(j)
	  si = ksite(vova); sj = ksite(lesha)
	  xi=x(:,ksite(vova)); xj=x(:,ksite(lesha))
	  ti=ktau(vova); tj=ktau(lesha)
	  a(i,j)=GREENFUN(sj, tj, si, ti)
	enddo; enddo

! invert
	det_big = full_inv(pm,lda,a)
	
	recalc_matrix = det_big       ! return the absolute value of the determinant

	end function recalc_matrix


!--------------------------------------------
!--- Arranges associations between sites
!--------------------------------------------
      subroutine ASSA  
 
      integer :: site, site1, i, i1, i2 !, i3 
      integer :: ic(d), NN(d+1) 
!     ic(i) is the i-coordinate of a site, 
!     the relation between site indexes and coordinates reads:
!     site = 1 + (ic(1)-1) + N(1)*(ic(2)-1) +...
!               + N(1)*...*N(d-1)*(ic(d)-1)

      NN(1)=1; do i=2,d+1; NN(i)=NN(i-1)*N(i-1); enddo
       
      do i=1,d; back(i)=i+d; back(i+d)=i; enddo 
      
      ic=1; ic(1)=0
      DO site=1, Nsite

!------------ Coordinates for site ----------
         i1=1 
         DO
         if (ic(i1) < N(i1)) then
            ic(i1)=ic(i1)+1
            DO i2=1,i1-1; ic(i2)=1; END DO
            EXIT
         else; i1=i1+1;  end if 

         END DO

!-------------------------------------------------------


         DO i=1, d
            if (ic(i) < N(i)) then
               site1=site+NN(i)
            else
               site1=site+NN(i)-NN(i+1)
            end if

            ass(i,site)=site1
            ass(back(i),site1)=site
	      x(:,site)=ic(:)
                 
         END DO
      END DO
      
      end subroutine ASSA


!--------------------------------------------
!--- Returns the names of the sites around the given one
!--------------------------------------------
      integer function cube(site0,half_edge,uzli)
	integer :: site0           ! center of the cube
	integer :: half_edge       ! half edge
	integer :: uzli(1:Nsite)   ! [output] names of the sites
                                 ! number of sites is returned via cube()
	integer :: i1,i2,i3, site1, site2, site3, ncurr, edge

! fill yes(:) array
	yes=.false.

! go to the far negative corner
	site1=site0
	do i1=1,half_edge; site1=ass(4,site1); enddo
	do i2=1,half_edge; site1=ass(5,site1); enddo
	do i3=1,half_edge; site1=ass(6,site1); enddo
	
	edge=2*half_edge

! now traverse the cube
	ncurr=1; site2=site1; site3=site1

	do i3=0,edge
	  do i2=0,edge
	     do i1=0,edge; !print*,ncurr
	        if(.not.yes(site3))then
			   uzli(ncurr)=site3; ncurr=ncurr+1; yes(site3)=.true.
			endif
	        site3=ass(1,site3)
	     enddo
	     site3=ass(2,site2); site2=site3
	  enddo
	  site2=ass(3,site1); site1=site2; site3=site2
	enddo

	cube=ncurr-1

	end function cube


!----------------------
!---- Tabulate site lists
!----------------------
	subroutine tab_cube(rx,n_in,uzli_cube)
	integer :: rx, n_in
	integer, dimension(1:n_in,1:Nsite) :: uzli_cube

	integer :: n_uz, site

	do site=1,Nsite

	n_uz = cube(site,rx,uzli)

	 if(n_uz/=n_in)then                  ! these must be equal unless something very weird happened
!		 open(OUT_UNIT,file=outputfname)
	     write(OUT_UNIT,*)'tab_cube: '
	     write(OUT_UNIT,*)'n_in = ', n_in,'   n_uz = ',n_uz
!	     close(OUT_UNIT)
	     call mystop
	 endif
	uzli_cube(:,site) = uzli(1:n_uz)


	enddo


	end subroutine tab_cube


!-----------------------------
!--- goodness-of-kink, GOK, used in 'leap' updates
!-----------------------------
      real*8 function GOK(x0,t0,x1,t1)
	integer, dimension(1:d) :: x0,x1  !site0,site1
	real*8  :: t0,t1

	real*8 :: dt, ddx
	integer :: dx, i

	GOK = 1.d200      ! default value: \infty if outside of the box


!	dt = abs(t1-t0); dt = min(dt, beta-dt)
!	if( dt<=rt_le )then


	dt=t1-t0; if(dt<0)dt=dt+beta

	if((dt<=rt_le).and.(dt>0))then

	   dt = dt/beta

	   ddx = 0.d0
	   do i=1,d
	      dx = abs(x0(i)-x1(i))
	      dx = min(dx,N(i)-dx)
	      ddx = ddx + ( 1.d0*dx/N(i) )**2
	   enddo

	   GOK = sqrt( ddx**2 + dt**2 )

	endif

	end function GOK



!--------------------------------------------
!--- Prepare weights for the 'leaps'
!--------------------------------------------
      subroutine tab_leaps 

	integer :: n_ev_total,site, ev,du3(d),ii,x1(d),x2(d),dx(d),j,i
	integer :: site0
	real*8 :: du, t1,w1, dovesok

	
	write(OUT_UNIT,*)'prepare mega-weights....'

	site=1; x1(:)=x(:,site); t1=0.d0
	allocate( uzli(1:Nsite),yes(1:Nsite) ); n_le=cube(site,rx_le,uzli)
	site0 = site

	n_ev_total=nt_same*n_le
	allocate(dx_le(d,n_ev_total),dt_le(n_ev_total),w_le(n_ev_total))

! get weights
	ev=0
	do i=1,n_le; site = uzli(i); x2=x(:,site) 
	 
	   do j=0,nt_same-1
			ev=ev+1
			dt_le(ev)=-ddt_same*(j+0.5d0)
			dx=x2-x1
			do ii=1,d;  if(dx(ii)>N2(ii))dx(ii)=dx(ii)-N(ii)   ! PBC
			enddo
	        dx_le(:,ev)=dx
						! "break the symmetry between the sites in a coordination sphere"
						! This is to say that the leaps, say, 'one site to the left' and
						! 'one site to the right' have the same  weights --- apart from the
						! round-off errors. Since different machines round off in 
						! a slightly different manner, the order of the events will be 
						! machine-dependent. While this is no problem, it causes 
						! the update sequences be machine-dependent, which leads
						! slightly different outputs --- and hence bug suspicion 
						! [the author spend three days looking for such a spurious bug].
			dovesok = ( dx(1)+dx(2)*N(1)+dx(3)*N(1)*N(2) ) * 1.d-14    
		    w_le(ev) = GREENFUN(site0,t1,site,t1+dt_le(ev))**2 + dovesok
	   enddo
	enddo

	deallocate(uzli,yes)

	write(OUT_UNIT,*)'...done'

! sort weights [bubble, no ne gum :)]
	write(OUT_UNIT,*)'sort weights....'
	do i=1,n_ev_total; do j=1,n_ev_total-1
	    if( w_le(j)<w_le(j+1) )then
	       du=w_le(j);  w_le(j)=w_le(j+1);   w_le(j+1)=du
	       du=dt_le(j); dt_le(j)=dt_le(j+1); dt_le(j+1)=du
	       du3(:)=dx_le(:,j); dx_le(:,j)=dx_le(:,j+1); 
		   dx_le(:,j+1)=du3(:)
	    endif
	enddo; enddo	
	write(OUT_UNIT,*)'.... done'


! cut off weights less than 10^{-4} of the max.
	w1=1.d0/w_le(1)
	do i=1,n_ev_total; if(w_le(i)*w1<1.d-4)exit
	enddo
	n_ev_cut=i-1
	write(OUT_UNIT,*)'cutoff done, n_ev_cut = ', n_ev_cut
	

! normalize weights
	w_norm = sum(w_le(1:n_ev_cut))
	w_le(1:n_ev_cut)=w_le(1:n_ev_cut)/w_norm

! build the search tree for the events based on the weigts
	call init_tree(n_ev_cut,w_le(1:n_ev_cut))


	write(OUT_UNIT,*)'....tab_leaps done'


	end subroutine tab_leaps



!----------------------
!---- Getname function for the Name Manager
!----------------------
      SUBROUTINE GetName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nick

      IF(nmnm<nmnm_max)THEN; 
	  nmnm=nmnm+1
        nick=namelist(nmnm); numnam(nick)=nmnm
      ELSE
	    ! has been moved into CheckGetName()
      ENDIF   

      END SUBROUTINE GetName

!------------------------------
!--- Check for GetName: nmnm>nmnm_max
!------------------------------
	subroutine CheckGetName

	if(nmnm>=nmnm_max)then; PRINT*,'GetName-> list is over!'
	   call mystop
	endif

	end subroutine CheckGetName

       
!----------------------
!---- DropName function for Name Manager
!----------------------
      SUBROUTINE DropName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nick
      INTEGER :: nm, lastname

      nm=numnam(nick)
      IF(nm<nmnm)THEN
         lastname=namelist(nmnm); namelist(nm)=lastname
         numnam(lastname)=nm; namelist(nmnm)=nick
         numnam(nick)=nmnm; nmnm=nmnm-1
      ELSE IF(nm==nmnm)THEN
		 nmnm=nmnm-1
      ELSE
	   ! has been moved into the CheckDropName()
      ENDIF
      RETURN 

      END SUBROUTINE DropName 


!------------------------------
!--- Check whether the name to be dropped exists
!------------------------------
	subroutine CheckDropName(nick)
	integer :: nick

	if(numnam(nick)>nmnm)then; 
	PRINT*,'DropName-> No such name:',nick; call mystop
	endif

	end subroutine CheckDropName



!---------------------------------
!--- Calculate the det from scratch & compare
!---------------------------------
	subroutine check_recalc

	real*8, allocatable :: a(:,:)
	real*8 :: det_1
	integer :: i,j

	if(pm<2)return
	
      allocate(a(lda,lda))
      det_1 = recalc_matrix(pm,lda,a)
	
      do i=1,pm; do j=1,pm
	    if(abs(a(i,j)-matr(i,j))>1d-8)then 
	       print*; print*,'check_recalc failed'
	       print*,i,j
		   print*,matr(i,j),a(i,j),abs(a(i,j)-matr(i,j))
	       print*,status,step,pm
	       print*,det_1
	       call mystop
	    endif
	enddo; enddo


	deallocate(a)


	end subroutine check_recalc


!---------------------------------
!--- Thermalize from scratch:
!      - diag updates only
!      - tune U
!---------------------------------
	subroutine therm1
	logical :: acpt
	real*8 :: r, pex, bun_in, bun_fin, det

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [diag] will be done by ', n_sw,' sweeps'
	close(OUT_UNIT)

	bun_in = beta*U_in*Nsite; bun_fin=bun

	step_t = n_sw*bun
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_t=i_t+1.d0
	   if(step>=step_t)exit

	   pex=step_t/i_t-1.d0

	   if (pex < 0.05) then; bun=bun_fin
	   else
            bun = bun_fin + (bun_in - bun_fin)*exp(-1.d0/pex)
	   end if

!--- Update: diagonal add/drop ONLY
	   r = rndm()
	   if     (r<0.5)then;    call add(acpt,det)
	   else if(r<=1.d0)then;  call drop(acpt,det)
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	enddo

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6, ' mln steps'
	write(OUT_UNIT,*)'  '
	close(OUT_UNIT)

	if(n_sw>0) call wrt   

	end subroutine therm1


!---------------------------------
!--- Thermalize via worm scheme
!---------------------------------
	subroutine therm2
	logical :: acpt
	real*8 :: r, det


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [worm] will be done by ',  n_sw,' sweeps'
	close(OUT_UNIT)

	step_t = n_sw*bun
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_t=i_t+1.d0; i_r = i_r + 1.d0
	   if(step>=step_t)exit


!--- Update :)
	   r = rndm()
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call create(acpt,det)
	   else if(r<prob(4))then;    call annihilate(acpt,det)
	   else if(r<prob(5))then;    call leap_add(acpt,det)
	   else if(r<prob(6))then;    call leap_drop(acpt,det)
	   else if(r<=prob(7))then;   call hop(acpt,det)
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	enddo

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6, ' mln steps'
	write(OUT_UNIT,*)'  '

	if(n_sw>0) call wrt   

	end subroutine therm2


!--------------------------
!--- clean up & stop 
!--------------------------
	subroutine mystop

	if(allocated(ass))deallocate(ass)
	if(allocated(back))deallocate(back); if(allocated(x))deallocate(x)
	if(allocated(uzli))deallocate(uzli)
	if(allocated(yes))deallocate(yes)
	if(allocated(uzli_cube_le))deallocate(uzli_cube_le)
	if(allocated(uzli_cube_ca))deallocate(uzli_cube_ca)
	if(allocated(dx_le))deallocate(dx_le)
	if(allocated(dt_le))deallocate(dt_le)
	if(allocated(w_le))deallocate(w_le)
	if(allocated(name_best))deallocate(name_best)
	if(allocated(dist_best))deallocate(dist_best)
	if(allocated(distance))deallocate(distance)
	if(allocated(dx_le))deallocate(dx_le)
	if(allocated(nkink))deallocate(nkink)
	if(allocated(kname))deallocate(kname)
	if(allocated(ksite))deallocate(ksite)
	if(allocated(ktau))deallocate(ktau)
	if(allocated(row))deallocate(row)
	if(allocated(clmn))deallocate(clmn)
	if(allocated(nm_row))deallocate(nm_row)
	if(allocated(nm_clmn))deallocate(nm_clmn)
	if(allocated(matr))deallocate(matr)
	if(allocated(m_u))deallocate(m_u)
	if(allocated(m_v))deallocate(m_v)
	if(allocated(m_w))deallocate(m_w)
	if(allocated(m_z))deallocate(m_z)
	if(allocated(m_u2))deallocate(m_u2)
	if(allocated(m_v2))deallocate(m_v2)
	if(allocated(m_s2))deallocate(m_s2)
	if(allocated(m_c2))deallocate(m_c2)


!	call MPI_FINALIZE( ierr )
	stop

	end subroutine mystop


	end program MAIN
