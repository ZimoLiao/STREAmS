subroutine part
!
! Lagrangian particle solver (explicit Euler method)
! TODO: TVD-RK3
!
 use mod_streams
 implicit none

 integer :: i,j,k
 integer :: il,jl,kl,ir,jr,kr
 real(mykind) :: ip,jp,kp,iq,jq,kq
 real(mykind) :: uatpart,vatpart,watpart
 real(mykind) :: dt
 
! TODO: consider time step for particle advancing carefully
 if (cfl>0._mykind) then
  dt = dtmin*cfl 
 else
  dt = dtmin
 endif

 !$cuf kernel do(1) <<<*,*>>> 
  do i=1,npart
   ! nearby 
   do ir=1,nx
    if (xpart_gpu(i)<x_gpu(ir)) exit
   enddo
   do jr=1,ny
    if (ypart_gpu(i)<y_gpu(jr)) exit
   enddo
   do kr=1,nz
    if (zpart_gpu(i)<z_gpu(kr)) exit
   enddo
   il=ir-1
   jl=jr-1
   kl=kr-1

   ! flow velocity interpolation (tri-linear)
   ip=xpart_gpu(i)-x_gpu(il)
   jp=xpart_gpu(j)-x_gpu(jl)
   kp=xpart_gpu(k)-x_gpu(kl)
   iq=1.0_mykind-ip
   jq=1.0_mykind-jp
   kq=1.0_mykind-kp
   
   uatpart=ip*jp*kp*wv_gpu(ir,jr,kr,2)+iq*jp*kp*wv_gpu(il,jr,kr,2)+ip*jq*kp*wv_gpu(ir,jl,kr,2)+iq*jq*kp*wv_gpu(il,jl,kr,2)+ip*jp*kq*wv_gpu(ir,jr,kl,2)+iq*jp*kq*wv_gpu(il,jr,kl,2)+ip*jq*kq*wv_gpu(ir,jl,kl,2)+iq*jq*kq*wv_gpu(il,jl,kl,2)
   vatpart=ip*jp*kp*wv_gpu(ir,jr,kr,3)+iq*jp*kp*wv_gpu(il,jr,kr,3)+ip*jq*kp*wv_gpu(ir,jl,kr,3)+iq*jq*kp*wv_gpu(il,jl,kr,3)+ip*jp*kq*wv_gpu(ir,jr,kl,3)+iq*jp*kq*wv_gpu(il,jr,kl,3)+ip*jq*kq*wv_gpu(ir,jl,kl,3)+iq*jq*kq*wv_gpu(il,jl,kl,3)
   watpart=ip*jp*kp*wv_gpu(ir,jr,kr,4)+iq*jp*kp*wv_gpu(il,jr,kr,4)+ip*jq*kp*wv_gpu(ir,jl,kr,4)+iq*jq*kp*wv_gpu(il,jl,kr,4)+ip*jp*kq*wv_gpu(ir,jr,kl,4)+iq*jp*kq*wv_gpu(il,jr,kl,4)+ip*jq*kq*wv_gpu(ir,jl,kl,4)+iq*jq*kq*wv_gpu(il,jl,kl,4)

   ! time advancing
   xpart_gpu(i) = xpart_gpu(i)+dt*upart_gpu(i)
   ypart_gpu(i) = ypart_gpu(i)+dt*vpart_gpu(i)
   zpart_gpu(i) = zpart_gpu(i)+dt*wpart_gpu(i)
   
   upart_gpu(i) = upart_gpu(i)+dt*(uatpart-upart_gpu(i))/stpart
   vpart_gpu(i) = vpart_gpu(i)+dt*(vatpart-vpart_gpu(i))/stpart
   wpart_gpu(i) = wpart_gpu(i)+dt*(watpart-wpart_gpu(i))/stpart
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
end subroutine part