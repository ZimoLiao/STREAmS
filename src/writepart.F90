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
  write(*,*)  'Writing part'
!
! 1D function file (PLOT3D)
  open(unit=10,file='part_'//nastore//'.bin',form='unformatted')
  write(10) npart,6
  do ind=1,npart
   write(10) xpart(ind)
  enddo
  do ind=1,npart
   write(10) ypart(ind)
  enddo
  do ind=1,npart
   write(10) zpart(ind)
  enddo
  do ind=1,npart
   write(10) upart(ind)
  enddo
  do ind=1,npart
   write(10) vpart(ind)
  enddo
  do ind=1,npart
   write(10) wpart(ind)
  enddo
  close(10)
!
 endif
!
end subroutine writepart