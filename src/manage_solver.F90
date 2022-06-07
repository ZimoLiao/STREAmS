subroutine manage_solver
!
 use mod_streams
 implicit none
!
 logical :: updatestat,savefield,saverst,savepart
!
!call write_wallpressure
!
 updatestat = .false.
 savefield  = .false.
 saverst    = .false.
 savepart   = .false.
!
 if (mod(icyc,istat)==0) updatestat = .true.
 if (telaps>tsol(istore)) savefield = .true.
 if (telaps>tsol_restart(istore_restart)) saverst = .true.
 if (telaps>tsol_part(istore_part)) savepart = .true.
!
 if (updatestat.or.savefield.or.saverst.or.savepart) then ! TODO: optimize (only copy particles)
  if (xrecyc>0._mykind) call recyc
  call updateghost()
  call prims()
  call copy_gpu_to_cpu()
 endif
!
!Statistics
 if (updatestat) then
  if (iflow==-1) then
  elseif (iflow==0) then
   call stats1d()
  else
   call stats2d()
  endif
 endif
!
!Writing fields
 if (savefield) then
  if (enable_plot3d>0) call writefield()
  if (enable_vtk>0) call writefieldvtk()
  if (iflow>0) call writestatzbl()
  istore = istore+1
 endif
!
 if (savepart) then
  call writepart() ! ADD(lzmo): TODO: enable option
  istore_part = istore_part+1
 endif
!
 if (saverst) then
  if (io_type==1) call writerst_serial()
  if (io_type==2) call writerst()
  if (iflow==-1) then
  elseif (iflow==0) then
   call writestat1d()
  else
   if (io_type==1) then
    call writestat2d_serial()
    call writedf_serial()
   else
    call writestat2d()
    call writedf()
   endif
  endif
  istore_restart = istore_restart+1
 endif
!
 if (updatestat.or.savefield.or.saverst.or.savepart) then
  call reset_cpu_gpu()
 endif
!
end subroutine manage_solver
