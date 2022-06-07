subroutine writepart
!
! Writing particles part.dat TODO: binary output
!
 use mod_streams
 implicit none
!
 integer :: i
!
 character(4) :: nastore
!
 if (masterproc) write(*,*) 'Storing part', istore,'at time', telaps
 write(nastore,1004) istore
 1004 format(I4.4)
!
 if (masterproc) then
!
  write(*,*)  'Writing part'
!
  open(unit=10,file='plotpart_'//nastore//'.dat',form='formatted')
  write(10,*) 'variables = x y z up vp wp'
  write(10,*) 'zone i=',npart
  write(10,*) 'solutiontime=',telaps
  do i=1,npart
   write(10,100) xpart(i),ypart(i),zpart(i),upart(i),vpart(i),wpart(i)
  100     format(20ES20.10)
  enddo
  close(10)
!
 endif
!
end subroutine writepart