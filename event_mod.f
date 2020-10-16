!-----------------------
!      Tree search 
!-----------------------
	module tree
	implicit none


	public :: make_tree, event, init_tree

	private
!	public

	integer :: nms                  ! max number of events
      integer :: namelast             ! last name in the list (real number of events)
      real*8, allocatable :: sp(:)    ! sp(name) is the sum of probabilities
	                                ! from the right end to the name 
      real*8, allocatable :: sprob(:) ! sprob (name) is the sum of probabilities
	                                ! from the left end to the name 	

	integer, allocatable :: left(:), right(:), up(:) ! tree links
      integer :: ev0                  ! starting point for the tree




!========================================================
	CONTAINS

!--------------------------------
!--- Create a tree given weight(..) & nms
!--------------------------------
      subroutine init_tree(num,weight)
	integer :: num
	real*8 :: weight(1:num)
	integer :: i


	nms=num; namelast=num     
	allocate(sp(0:nms),sprob(nms),left(nms),right(nms),up(nms))


	sprob(1)=weight(1)
	do i=2,namelast; 
	   sprob(i)=sprob(i-1)+weight(i)
	enddo

	sp(namelast)=0.d0
	do i=namelast-1,0,-1
	  sp(i)=sp(i+1)+weight(i+1)
	enddo


      call MAKE_TREE

!	print*,'ev0 = ',ev0


	deallocate(sp)


	end subroutine init_tree


!--------------------------------
!--- Create the tree out of known sp(..) & sprob(...)
!--------------------------------
      subroutine MAKE_TREE
	integer :: name, na, nb, n1  
	double precision :: d, d2, spa, spb


	na=0
      spa=sp(0)
      name= namelast
!---------- Left dynasty for name ------------------

 1    continue
	n1=na+1

      if(n1==name) then ! means left dead end for name
        left(name)=-1
        go to 2  ! go to create right daughter for name
      end if	  
	 
	spb=sp(name)
      d2=0.5000000001d0*(spa-spb)     ! exactly 0.5 => stucks for equal weights
      do
	  d=sp(n1)-spb
        if(d < d2) then
          left(name)=n1; up(n1)=name
          name=n1
	    go to 1
        end if
        n1=n1+1
	end do    
!-----------end left dynasty for name -----------------


!---------- Right dauter for name ------------------
 2    continue

      spa=sp(name)    
 
      ! looking for nb (right boundary of the interval)
      nb=name
	do
	   nb=up(nb)
	   spb=sp(nb)
         if (spb < spa) exit
      end do
	! end looking for nb

	n1=name+1

      if(n1==nb) then ! means right dead end for name 
        right(name)=-1
        go to 3  ! go to process right dead end for name 


      end if

      d2=0.500000001d0*(spa-spb) ! exactly 0.5 => stucks for equal weights
      do
	  d=sp(n1)-spb
        if(d < d2) then
          right(name)=n1; up(n1)=name
	    na=name
          name=n1
	    spa=sp(na)
	    go to 1 ! go to left dynasty for name
        end if
        n1=n1+1
      end do
!---------- end right dauter for name --------------



!---------- Processing right dead end for name ------------------
 3    continue  

 
      ! looking for new name 
	do
	   name=up(name)
	   if(name == namelast) then ! means the tree is constructed
	      ev0=left(namelast)
	      return 
	   end if
	   spb=sp(name)
         if (spb < spa) exit
      end do
	! end looking for new name

      go to 2


!---------- end of processing right dead end for name --------------

      end subroutine MAKE_TREE



!--------------------------------
!--- Return the event #
!--------------------------------
      integer function EVENT(x)

      double precision :: x !random number
      integer :: next

      EVENT=ev0 !starting point

	DO

	if (x <= sprob(EVENT)) then
	    next= left(EVENT)
          if (next < 0) return
      else
	    next=right(EVENT)
	    if (next < 0) then; EVENT=EVENT+1; return; end if  
      end if

	EVENT=next

      END DO
	end function EVENT

	end module tree

