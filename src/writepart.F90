subroutine writepart
!
! Writing particle positions and velocities
!
 use mod_streams
 implicit none
!
 integer :: ind
!
 character(4) :: nastore
!
 if (masterproc) write(*,*) 'Storing part', istore_part, 'at time', telaps
 write(nastore,1004) istore_part
 1004 format(I4.4)
!
 if (masterproc) then
!
! write binary file
  open(unit=10,file='part_'//nastore//'.bin',form='unformatted')
  write(10) telaps ! solution time
  write(10) npart ! number of particles
  do ind=1,npart
   write(10) xpart(ind),ypart(ind),zpart(ind),upart(ind),vpart(ind),wpart(ind)
  enddo
  close(10)
!
! write tecplot ascii file
  open(unit=10,file='part_'//nastore//'.dat',form='formatted')
  write(10,*) 'variables = x y z up vp wp'
  write(10,*) 'zone i = ',npart
  write(10,*) 'solutiontime = ',telaps
  do ind=1,npart
   write(10,100) xpart(ind),ypart(ind),zpart(ind),upart(ind),vpart(ind),wpart(ind)
  100 format(20ES20.10)
  enddo 
  close(10)
!
 endif
!
end subroutine writepart