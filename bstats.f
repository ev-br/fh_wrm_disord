!-----------------------------------------------------------
!  Block statistic helper routines
!-----------------------------------------------------------
	module bstats
	implicit none; save
	
	public :: bSTAT, collate
	public :: mrg, mrg_save, mrg_conv
	public :: bool_to_str, bool_to_char
	
!*****************************************************************
	contains


!-------------------------------
!--- Analyze block statistics: average and dispersion
!-------------------------------
	subroutine bSTAT(arr,n,Zb,av,err)
	real*8, dimension(*) :: arr
	integer              :: n
	real*8               :: Zb
	real*8               :: av, err

	real*8 :: av2, dif

	av  = sum( arr(1:n) )/Zb/n
	av2 = sum( arr(1:n)**2 )/Zb/Zb/n

				!av2 = av2 + (arr(j)/Zb)**2

	dif = av2 - av*av; 	if(dif<0.d0)dif=0.d0

	err = sqrt( dif ) / sqrt(1.d0*n)


	end subroutine bSTAT



!------------------------
!--- Collate the array
!------------------------
      subroutine collate(arr,n)
      real*8, dimension(*) :: arr
      integer :: n, Zb        ! array length & # of measurements per array entry

      integer :: i

      do i=1,n/2
          arr(i)=arr(2*i)+arr(2*i-1)
      enddo

      end subroutine collate



!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
      subroutine mrg(arr,n,Zb, OU)
      integer, intent(in)              :: n, Zb
      real*8, dimension(1:n), intent(in) :: arr
      integer, intent(in) :: OU

      real*8  :: av, err, arr1(1:n), zb1
      integer :: i, n1


!
!  2 Sept 2007 :   Calls bSTAT(..) w/ real*8 Zb1
!

      zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n            

      write(OU, *)'-----------' !, int(log(1.d0*n)/log(2.d0))+1

      do;

! emit
        call bSTAT(arr1,n1,zb1,av,err)
        write(OU, fmt=777) av, err,n1               
 777    format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
        if(n1<3)exit

! merge
        n1=INT(n1/2); zb1=zb1*2.d0
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

      enddo

      WRITE(OU,*) '------------'

      end subroutine mrg



!-------------------------------
!--- Merge blocks & RETURN av +/- err 
!-------------------------------
    subroutine mrg_save(arr,n,Zb,n_o,av_o,err_o)
    integer, intent(in)              :: n, Zb
    real*8, dimension(1:n), intent(in) :: arr

    integer :: n_o
    real*8, dimension(*), intent(out) :: av_o,err_o

    real*8  :: av, err, arr1(1:n), zb1
    integer :: i, n1, jj
!
!  2 Sept 2007 :   Calls bSTAT(..) w/ real*8 Zb1
!
    zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n;  jj=1

    do;

! emit
       call bSTAT(arr1,n1,zb1,av,err)
       av_o(jj) = av; err_o(jj)=err
	   jj=jj+1

! enough?
       if(n1<3)then; exit
	   endif

! merge
       n1=INT(n1/2); zb1=zb1*2.d0
       do i=1,n1
           arr1(i) =  arr1(2*i-1) + arr1(2*i)
       enddo

    enddo

    end subroutine mrg_save



!-------------------------------
!--- Merge blocks, check convergence & RETURN av +/- err [if conv]
!-------------------------------
	subroutine mrg_conv(arr,n,Zb,av_o,err_o, conv)
	implicit none
	integer, intent(in)              :: n, Zb
	real*8, dimension(1:n), intent(in) :: arr
	real*8, dimension(1:n)  :: av_t,err_t
	real*8 :: av_o,err_o
	logical :: conv

	real*8, parameter :: LIMIT = 0.05    ! convergence limit: if an errorbar fluctuates within 5%, then it has probably converged :)

	real*8  :: av, err, arr1(1:n), zb1
	integer :: i, n1, jj
    real*8 :: prev, dummy
	real*8, allocatable :: delta(:)

    zb1 = 1.d0*zb;       arr1(1:n) = arr(1:n); n1=n;  jj=1

    do;

! emit
       call bSTAT(arr1,n1,zb1,av,err)
	   av_t(jj) = av; err_t(jj)=err
	   jj=jj+1

! enough?
       if(n1<3)then; exit
	   endif

! merge
       n1=INT(n1/2); zb1=zb1*2.d0
       do i=1,n1
           arr1(i) =  arr1(2*i-1) + arr1(2*i)
       enddo

    enddo

! av and err: av is 2nd last
	jj = jj-1              ! jj is now the # of entry in the blocking array
	av_o = av_t(jj-1)
	err_o = maxval( err_t(1:jj) )

! check convergence
	conv=.false.

	if( jj>4 )then          ! it makes sense checking
		allocate(delta(1:jj-1))
 		prev=err_t(jj)	
		do i=jj-1,1,-1
		    dummy = 0.5*( prev+err_t(i) )
		    delta(i) = (prev-err_t(i))/dummy  
		enddo
		
		do i=jj-3,1,-1
			if( all( (/delta(i),delta(i+1),delta(i+2)/) <LIMIT ) )then
				conv=.true.
				exit
			endif
		enddo

		deallocate(delta)
	endif

    end subroutine mrg_conv



!--------------------------------------
!--- Utility: convert bool to char for printing
!--------------------------------------
    character*1 function bool_to_char(conv)
    implicit none
    logical, intent(in) :: conv
    
    if(conv)then
        bool_to_char = 'T'
    else
        bool_to_char = 'F'
    endif
    end function bool_to_char


!--------------------------------------
!--- Utility: convert bool to char for printing
!--------------------------------------
    character*8 function bool_to_str(conv)
    implicit none
    logical, intent(in) :: conv
    
    if(conv)then
        bool_to_str = '/conv/  '
    else
        bool_to_str = '/unconv/'
    endif
    end function bool_to_str


	end module bstats
