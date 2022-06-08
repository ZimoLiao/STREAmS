subroutine part
!
! Lagrangian particle solver (explicit Euler method)
! TODO: TVD-RK3
!
 use mod_streams
 implicit none

 integer :: ind
 integer :: il,jl,kl,ir,jr,kr
 real(mykind) :: ip,jp,kp,iq,jq,kq
 real(mykind) :: uatpart,vatpart,watpart
 real(mykind) :: dt
 real(mykind) :: ry,rz
 
! TODO: consider time step for particle advancing carefully
 if (cfl>0._mykind) then
  dt = dtmin*cfl 
 else
  dt = dtmin
 endif
!
! integration
!
 !$cuf kernel do(1) <<<*,*>>> 
  do ind=1,npart
   ! nearby 
   do ir=1,nx
    if (xpart_gpu(ind)<x_gpu(ir)) exit
   enddo
   do jr=1,ny
    if (ypart_gpu(ind)<y_gpu(jr)) exit
   enddo
   do kr=1,nz
    if (zpart_gpu(ind)<z_gpu(kr)) exit
   enddo
   il=ir-1
   jl=jr-1
   kl=kr-1

   ! flow velocity interpolation (tri-linear)
   ip=xpart_gpu(ind)-x_gpu(il)
   jp=ypart_gpu(ind)-y_gpu(jl)
   kp=zpart_gpu(ind)-z_gpu(kl)
   iq=1.0_mykind-ip
   jq=1.0_mykind-jp
   kq=1.0_mykind-kp
   
   uatpart=ip*jp*kp*wv_gpu(ir,jr,kr,2)+iq*jp*kp*wv_gpu(il,jr,kr,2)+ip*jq*kp*wv_gpu(ir,jl,kr,2)+iq*jq*kp*wv_gpu(il,jl,kr,2)+ip*jp*kq*wv_gpu(ir,jr,kl,2)+iq*jp*kq*wv_gpu(il,jr,kl,2)+ip*jq*kq*wv_gpu(ir,jl,kl,2)+iq*jq*kq*wv_gpu(il,jl,kl,2)
   vatpart=ip*jp*kp*wv_gpu(ir,jr,kr,3)+iq*jp*kp*wv_gpu(il,jr,kr,3)+ip*jq*kp*wv_gpu(ir,jl,kr,3)+iq*jq*kp*wv_gpu(il,jl,kr,3)+ip*jp*kq*wv_gpu(ir,jr,kl,3)+iq*jp*kq*wv_gpu(il,jr,kl,3)+ip*jq*kq*wv_gpu(ir,jl,kl,3)+iq*jq*kq*wv_gpu(il,jl,kl,3)
   watpart=ip*jp*kp*wv_gpu(ir,jr,kr,4)+iq*jp*kp*wv_gpu(il,jr,kr,4)+ip*jq*kp*wv_gpu(ir,jl,kr,4)+iq*jq*kp*wv_gpu(il,jl,kr,4)+ip*jp*kq*wv_gpu(ir,jr,kl,4)+iq*jp*kq*wv_gpu(il,jr,kl,4)+ip*jq*kq*wv_gpu(ir,jl,kl,4)+iq*jq*kq*wv_gpu(il,jl,kl,4)

   ! time advancing
   xpart_gpu(ind) = xpart_gpu(ind)+dt*upart_gpu(ind)
   ypart_gpu(ind) = ypart_gpu(ind)+dt*vpart_gpu(ind)
   zpart_gpu(ind) = zpart_gpu(ind)+dt*wpart_gpu(ind)
   
   upart_gpu(ind) = upart_gpu(ind)+dt*(uatpart-upart_gpu(ind))/stpart
   vpart_gpu(ind) = vpart_gpu(ind)+dt*(vatpart-vpart_gpu(ind))/stpart
   wpart_gpu(ind) = wpart_gpu(ind)+dt*(watpart-wpart_gpu(ind))/stpart
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
! boundary conditions
!
 if (iflow==0) then ! channel
 !$cuf kernel do(1) <<<*,*>>> 
  do ind=1,npart
    
    xpart_gpu(ind) = mod(xpart_gpu(ind),rlx) ! periodic in x-/z-directions
    if (xpart_gpu(ind)<0) xpart_gpu(ind) = xpart_gpu(ind)+rlx
    zpart_gpu(ind) = mod(zpart_gpu(ind),rlz)
    if (zpart_gpu(ind)<0) zpart_gpu(ind) = zpart_gpu(ind)+rlz
    
    if (ypart_gpu(ind)<-rly*0.5) then ! reflective boundary (non physical)
     ypart_gpu(ind) = -rly-ypart_gpu(ind)
     vpart_gpu(ind) = -vpart_gpu(ind)
    elseif (ypart_gpu(ind)>rly*0.5) then
     ypart_gpu(ind) = rly-ypart_gpu(ind)
     vpart_gpu(ind) = -vpart_gpu(ind)
    endif

  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (iflow==1) then ! boundary-layer
 !$cuf kernel do(1) <<<*,*>>> 
  do ind=1,npart

    zpart_gpu(ind) = mod(zpart_gpu(ind),rlz) ! periodic in z-directions
    if (zpart_gpu(ind)<0) zpart_gpu(ind) = zpart_gpu(ind)+rlz

    if (ypart_gpu(ind)<0) then ! reflective boundary (non physical)
      ypart_gpu(ind) = -ypart_gpu(ind)
      vpart_gpu(ind) = -vpart_gpu(ind)
    endif

    ! recycling boundary
    if (ypart_gpu(ind)>rly.or.xpart_gpu(ind)>rlx.or.xpart_gpu(ind)<0) then
      xpart_gpu(ind)=x_gpu(1)
      ypart_gpu(ind)=randy_gpu(mod(irand,nrand)+1)*rly ! TODO: check ramdom number
      zpart_gpu(ind)=randz_gpu(mod(irand,nrand)+1)*rlz

      ! set particle velocity equal to nearby flow velocity approximately
      do jr=1,ny
        if (ypart_gpu(ind)<y_gpu(jr)) exit
      enddo
      do kr=1,nz
        if (zpart_gpu(ind)<z_gpu(kr)) exit
      enddo

      upart_gpu(ind)=wv_gpu(1,jr,kr,2)
      vpart_gpu(ind)=wv_gpu(1,jr,kr,3)
      wpart_gpu(ind)=wv_gpu(1,jr,kr,4)

      irand = irand+1
    endif 

  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 endif
end subroutine part