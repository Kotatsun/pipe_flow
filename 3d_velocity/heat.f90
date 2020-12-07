
program heat
!  SMAC (original version by Dr. T. Ushijima)
!  Simplified Marker and Cell method
!  Solving Heat Transfer in flow between parallel walls.
implicit double precision(a-h,o-z)
integer, parameter:: n=60, nx=n*2, ny=n, nz=n                                          ! 格子数n
real(8) :: p(0:nx+1,0:ny+1,0:nz+1)                                                    & ! p:圧力
          ,u(0:nx+1,0:ny+1,0:nz+1),v(0:nx+1,0:ny+1,0:nz+1),w(0:nx+1,0:ny+1,0:nz+1)    & ! u,v:速度
          ,phi(0:nx+1,0:ny+1,0:nz+1),divup(1:nx,1:ny,0:nz+1)                          & ! phi: 補正圧力, divup: 予測速度の発散
          ,up(0:nx+1,0:ny+1,0:nz+1),vp(0:nx+1,0:ny+1,0:nz+1),wp(0:nx+1,0:ny+1,0:nz+1)   ! up, vp: 予測速度
real(8), parameter:: re=750.d0                                                        & ! レイノルズ数
                    ,r02=900d0                                                        & ! 血管径
                    ,r12=400d0                                                          ! カテーテル径
character*6:: num                                                                       ! ファイル番号用文字列
integer,parameter:: stepmax=20000                                                       ! 時間ステップ数
integer:: step                                                                          ! 時間ステップ
! grid size 
dy=0.1d0                                    
dx=dy
dz=dy
ddx=1.d0/dx
ddy=1.d0/dy
ddz=1.d0/dz
ddx2=ddx*ddx
ddy2=ddy*ddy
ddz2=ddz*ddz
! time step
dt=0.002d0                                          ! time step
! dt=min(dt,0.25*dx) ! 1時間ステップで流体が移流によって飛び出さないようにする
! dt=min(dt,0.2*re*dx*dx) ! 拡散の影響の考慮 (ここではプラントル数の影響は考慮していない)
! write(6,*) 'dt = ',dt
ddt=1.d0/dt

! initial condition
! icont=1 のときは，既存のデータを用いる. それ以外は，初期化する.
icont=0
if (icont.eq.1) then
  open(unit=9, file='flow000001.dat', form='unformatted', access="stream")
  read(9) u,v,w,p
  close(9)
else
  do k=0,nz+1
    do j=0,ny+1
      do i=0,nx+1
        u(i,j,k)=0.d0
        r2 = dble(k-225)**2+dble(j-225)**2
        if (r2<=r02) then
          u(i,j,k)=2*(1-r2/r02)
        endif
        v(i,j,k)=0.d0
        w(i,j,k)=0.d0
        p(i,j,k)=0.d0
      enddo
    enddo
  enddo
end if

! boundary condition 1
un=0.d0
uw=1.d0
us=0.d0
ue=0.d0
uf=0.d0
ub=0.d0
vn=0.d0
vw=0.d0
vs=0.d0
ve=0.d0
vf=0.d0
vb=0.d0
wn=0.d0
ww=0.d0
ws=0.d0
we=0.d0
wf=0.d0
wb=0.d0

!============== time marching ================
do step=1,stepmax
  if(mod(step,stepmax/2000)==0) write(6,*) 'step =',step
! boundary condition 2
!     right wall (east) or outlet
  do j=0,ny
    do k=0,nz
      w(nx+1,j,k)=w(nx,j,k)
      v(nx+1,j,k)=v(nx,j,k)
      u(nx,j,k)=u(nx-1,j,k)
      p(nx+1,j,k)=0.d0
!     left wall (west) or inlet
      w(0,j,k)=ww
      v(0,j,k)=vw
      u(0,j,k)=uw
      p(0,j,k)=p(1,j,k)
    enddo
  enddo
!         u(0,ny/2+1)=uw
!     lower wall (south)
  do i=0,nx
    do k=0,nz
      w(i,0,k)=2.*ws-w(i,1,k)
      u(i,0,k)=2.*us-u(i,1,k)
      v(i,0,k)=vs
      p(i,0,k)=p(i,1,k)
!     upper wall (north)
      w(i,ny+1,k)=2.*wn-w(i,ny,k)
      u(i,ny+1,k)=2.*un-u(i,ny,k)
      v(i,ny,k)=vn
      p(i,ny+1,k)=p(i,ny,k)
    enddo
  enddo
!     forward wall (forward)
  do i=0,nx
    do j=0,ny
      w(i,j,0)=2.*wf-w(i,j,1)
      u(i,j,0)=2.*uf-u(i,j,1)
      v(i,j,0)=vf
      p(i,j,0)=p(i,j,1)
!     backward wall (backward)
      w(i,j,nz)=wb
      u(i,j,nz+1)=2.*ub-u(i,j,nz)
      v(i,j,nz+1)=vb
      p(i,j,nz+1)=p(i,j,nz)
    enddo
  enddo

! compute the thermal field
  do k=1,nz
    do j=1,ny
      do i=1,nx
        cnvt=( ddx*((the(i,j,k)+the(i-1,j,k))*u(i-1,j,k)         &
                   -(the(i,j,k)+the(i+1,j,k))*u(i,j,k))          &
              +ddy*((the(i,j,k)+the(i,j-1,k))*v(i,j-1,k)         &
                   -(the(i,j,k)+the(i,j+1,k))*v(i,j,k))          &
              +ddz*((the(i,j,k)+the(i,j,k-1))*w(i,j,k-1)         &
                   -(the(i,j,k)+the(i,j,k+1))*w(i,j,k)) )/2.d0    

        dift=( ddx2*(the(i+1,j,k)-2.*the(i,j,k)+the(i-1,j,k))           &
              +ddy2*(the(i,j+1,k)-2.*the(i,j,k)+the(i,j-1,k))           &
              +ddz2*(the(i,j,k+1)-2.*the(i,j,k)+the(i,j,k-1)) )/(pr*re)

        the0(i,j,k)=the(i,j,k)+dt*(cnvt+dift)
      enddo
    enddo
  enddo
  do k=1,nz
    do j=1,ny
      do i=1,nx
        the(i,j,k)=the0(i,j,k)
      enddo
    enddo
  enddo

! =====================================
! predictor step
! for u_ijk
  do k=1,nz
    do j=1,ny
      do i=1,nx-1
        cnvu=ddx*((u(i+1,j,k)+u(i,j,k))**2                                &
                 -(u(i-1,j,k)+u(i,j,k))**2)/4.d0                          &
            +ddy*((u(i,j+1,k)+u(i,j,k))*(v(i+1,j,k)+v(i,j,k))             &
                 -(u(i,j,k)+u(i,j-1,k))*(v(i,j-1,k)+v(i+1,j-1,k)))/4.d0   &
            +ddz*((u(i,j,k+1)+u(i,j,k))*(w(i+1,j,k)+w(i,j,k))             &
                 -(u(i,j,k)+u(i,j,k-1))*(w(i,j,k-1)+w(i+1,j,k-1)))/4.d0   

        fijk=-ddx*(p(i+1,j,k)-p(i,j,k))-cnvu                      &
            +ddx2*(u(i+1,j,k)-2.d0*u(i,j,k)+u(i-1,j,k))/re        &
            +ddy2*(u(i,j+1,k)-2.d0*u(i,j,k)+u(i,j-1,k))/re        &
            +ddz2*(u(i,j,k+1)-2.d0*u(i,j,k)+u(i,j,k-1))/re  
        
        up(i,j,k)=u(i,j,k)+dt*fijk
      enddo
    enddo
  enddo
  ! right and left wall
  do j=0,ny
    do k=0,nz
      up(0,j,k)=uw
      up(nx,j,k)=up(nx-1,j,k)
    enddo
  enddo
  ! upper and lower wall
  do i=0,nx-1
    do k=0,nz
      up(i,0,k)=2.*us-up(i,1,k)
      up(i,ny+1,k)=2.*un-up(i,ny,k)
    enddo
  enddo
  ! forward and backward wall
  do i=0,nx-1
    do j=0,ny
      up(i,j,0)=2.*uf-up(i,j,1)
      up(i,j,nz+1)=2.*ub-up(i,j,nz)
    enddo
  enddo
! --------------------------------------------
!   for v_ijk
  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        cnvv=ddx*((u(i,j+1,k)+u(i,j,k))*(v(i+1,j,k)+v(i,j,k))               &
                 -(u(i-1,j+1,k)+u(i-1,j,k))*(v(i-1,j,k)+v(i,j,k)))/4.d0     &
            +ddy*((v(i,j+1,k)+v(i,j,k))**2                                  &
                 -(v(i,j,k)+v(i,j-1,k))**2)/4.d0                            &
            +ddz*((v(i,j,k+1)+v(i,j,k))*(w(i,j+1,k)+w(i,j,k))               &
                 -(v(i,j,k)+v(i,j,k-1))*(w(i,j,k-1)+w(i,j+1,k-1)))/4.d0    

        gijk=-ddy*(p(i,j+1,k)-p(i,j,k))-cnvv                            &
            +ddx2*(v(i+1,j,k)-2.d0*v(i,j,k)+v(i-1,j,k))/re              &
            +ddy2*(v(i,j+1,k)-2.d0*v(i,j,k)+v(i,j-1,k))/re              &
            +ddz2*(v(i,j,k+1)-2.d0*v(i,j,k)+v(i,j,k-1))/re  

        vp(i,j,k)=v(i,j,k)+dt*gijk
      enddo
    enddo
  enddo
  ! right and left wall
  do j=0,ny-1
    do k=0,nz
      vp(nx+1,j,k)=vp(nx,j,k)
      vp(0,j,k)=vw
    enddo
  enddo
  ! upper and lower wall
  do i=0,nx
    do k=0,nz
      vp(i,0,k)=vs
      vp(i,ny,k)=vn
    enddo
  enddo
  ! forward and backward wall
  do i=0,nx
    do j=0,ny-1
      vp(i,j,0)=vf
      vp(i,j,nz+1)=vb
    enddo
  enddo
! --------------------------------------------
!   for w_ijk
  do k=1,nz-1
    do j=1,ny
      do i=1,nx
        cnvw=ddx*((u(i,j,k+1)+u(i,j,k))*(w(i+1,j,k)+w(i,j,k))               &
                 -(u(i-1,j,k+1)+u(i-1,j,k))*(w(i-1,j,k)+w(i,j,k)))/4.d0     &
            +ddy*((v(i,j,k+1)+v(i,j,k))*(w(i,j+1,k)+w(i,j,k))               &
                 -(v(i,j-1,k+1)+v(i,j-1,k))*(w(i,j-1,k)+w(i,j,k)))/4.d0     &
            +ddz*((w(i,j,k+1)+w(i,j,k))**2                                  &
                 -(w(i,j,k)+w(i,j,k-1))**2)/4.d0                            

        hijk=-ddz*(p(i,j,k+1)-p(i,j,k))-cnvw                            &
            +ddx2*(w(i+1,j,k)-2.d0*w(i,j,k)+w(i-1,j,k))/re              &
            +ddy2*(w(i,j+1,k)-2.d0*w(i,j,k)+w(i,j-1,k))/re              &
            +ddz2*(w(i,j,k+1)-2.d0*w(i,j,k)+w(i,j,k-1))/re  

        wp(i,j,k)=w(i,j,k)+dt*hijk
      enddo
    enddo
  enddo
  ! right and left wall
  do j=0,ny
    do k=0,nz-1
      wp(nx+1,j,k)=wp(nx,j,k)
      wp(0,j,k)=ww
    enddo
  enddo
  ! upper and lower wall
  do i=0,nx
    do k=0,nz-1
      wp(i,0,k)=ws
      wp(i,ny+1,k)=wn
    enddo
  enddo
  ! forward and backward wall
  do i=0,nx
    do j=0,ny
      wp(i,j,0)=wf
      wp(i,j,nz)=wb
    enddo
  enddo
! =====================================
! =====================================
  ! write(6,*) 'evaluate continuity'
  ic=0
  div=0.0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        divup(i,j,k)=ddx*(up(i,j,k)-up(i-1,j,k))     &
                    +ddy*(vp(i,j,k)-vp(i,j-1,k))     &
                    +ddz*(wp(i,j,k)-wp(i,j,k-1))  
        div=div+divup(i,j,k)**2
        ic=ic+1
      enddo
    enddo
  enddo
! =====================================
! =====================================
! Solve the poisson equation (nabla)2 p=(nabla)up/dt by SOR
  ! write(6,*) 'solve the poisson equation for pressure'
! initialization
  do k=0,nz+1
    do i=0,nx+1
      do j=0,ny+1
        phi(i,j,k)=0.d0
      enddo
    enddo
  enddo
  eps=1.D-6
! 最大反復数 maxitr
! 収束させるためには maxitr=nx*ny 程度必要 
! 最大反復数は小さめにとっているので，計算の初めは 収束しないが，定常計算では影響なし.
  maxitr=nx*ny*nz/2
! 緩和係数 alpha=1.5~1.7 程度に設定
  alpha=1.7
! Start iteration
  do iter=1,maxitr
    error=0.d0
!$omp parallel do private(i,j,k,rhs,den,resid,dphi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhs=ddt*divup(i,j,k)
          resid=ddx2*(phi(i-1,j,k)-2.d0*phi(i,j,k)+phi(i+1,j,k))  &
               +ddy2*(phi(i,j-1,k)-2.d0*phi(i,j,k)+phi(i,j+1,k))  &
               +ddz2*(phi(i,j,k-1)-2.d0*phi(i,j,k)+phi(i,j,k+1))  &
               -rhs
          den=2.d0*(ddx2+ddy2+ddz2)
          dphi=alpha*resid/den
          error=max(abs(dphi),error)
          phi(i,j,k)=phi(i,j,k)+dphi
        enddo
      enddo
    enddo
!$omp end parallel do
    ! right and left wall
    do j=0,ny
      do k=0,nz
        phi(nx+1,j,k)=0.d0
        phi(0,j,k)=phi(1,j,k)
      enddo
    enddo
    ! upper and lower wall
    do i=0,nx
      do k=0,nz
        phi(i,0,k)=phi(i,1,k)
        phi(i,ny+1,k)=phi(i,ny,k)
      enddo
    enddo
    ! forward and backward wall
    do i=0,nx
      do j=0,ny
        phi(i,j,0)=phi(i,j,1)
        phi(i,j,nz+1)=phi(i,j,nz)
      enddo
    enddo
  ! 収束の判定
    if (error.lt.eps) exit
  enddo
  write(6,*) 'iter=  ', iter, error
! =====================================
! =====================================
! correct step
  do k=1,nz
    do j=1,ny
      do i=1,nx-1
        u(i,j,k)=up(i,j,k)-dt*ddx*(phi(i+1,j,k)-phi(i,j,k))
      enddo  
    enddo
  enddo
  do k=1,nz
    do j=1,ny-1
      do i=1,nx
        v(i,j,k)=vp(i,j,k)-dt*ddy*(phi(i,j+1,k)-phi(i,j,k))
      enddo  
    enddo
  enddo
  do k=1,nz-1
    do j=1,ny
      do i=1,nx
        w(i,j,k)=wp(i,j,k)-dt*ddz*(phi(i,j,k+1)-phi(i,j,k))
      enddo  
    enddo
  enddo
  do k=1,nz
    do j=1,ny
      do i=1,nx
        p(i,j,k)=p(i,j,k)+phi(i,j,k)
      enddo
    enddo
  enddo
! =====================================
! check the continuity for n+1 th step
! write(6,*) 'check the continuity for n+1 th step'
  ic=0
  div=0.0
  do k=1,nz
    do j=1,ny
      do i=1,nx
          divup(i,j,k)=ddx*(u(i,j,k)-u(i-1,j,k))    &
                      +ddy*(v(i,j,k)-v(i,j-1,k))    &
                      +ddz*(w(i,j,k)-w(i,j,k-1))
          ic=ic+1
          div=div+divup(i,j,k)**2
      enddo
    enddo
  enddo

! ********************************************* 
! 混合平均温度の計算(i=nx における値)
! *********************************************
  ! iN=nx
  ! tm=0.d0
  ! ! um=0.d0
  ! do j=1,ny
  !   do k=0,nz
  !     tm=tm+0.5d0*(u(iN,j,k)+u(iN-1,j,k))*the(iN,j,k)*dy
  !     !  um=um+0.5d0*(u(iN,j)+u(iN-1,j))*dy
  !   enddo
  ! enddo
  ! ! tm=tm/um
  ! tw=0.5d0*(the(iN,0,k)+the(iN,1,k)) むりっぽい
  ! s_Nu=q_w/(tw-tm)
  ! ! write(6,*) s_Nu

  ! =============== output =================
  if(mod(step, stepmax/stepmax) == 0) then
    ! ! calculate stream function
    ! do i=0,nx
    !   psi(i,0,int(nz/2))=0.d0
    !   do j=1,ny+1
    !     psi(i,j,int(nz/2))=psi(i,j-1,int(nz/2))+0.5*dy*(u(i,j-1,int(nz/2))+u(i,j,int(nz/2)))
    !   enddo
    ! enddo
    ! output file
    num='      '
    write(num,'(1i6)') step
    do j=1,6
      if(num(j:j)==' ') num(j:j)='0'
    end do
    open(unit=11, iostat=ios, file='output/flow'//num//'.txt')
    ! if (ios /= 0) then
    !   write(*,*) 'Failed to open file for output'
    !   stop
    ! end if
    do j=1, ny
      do i=1, nx
        x0=dx*(dble(i)-0.5)
        y0=dy*(dble(j)-0.5)
        u0=0.5*(u(i,j,int(nz/2))+u(i-1,j,int(nz/2)))
        v0=0.5*(v(i,j,int(nz/2))+v(i,j-1,int(nz/2)))
        p0=p(i,j,int(nz/2))
        divu0=divup(i,j,int(nz/2))
        ! psi0=0.5*(psi(i,j,int(nz/2))+psi(i-1,j,int(nz/2)))
        ! x座標, y座標, 速度u, 速度v, 圧力P, 連続の式, 流れ関数, 温度
        ! write(11,'(2f20.15, 1x, 6f20.15)') x0, y0, u0, v0, p0, divu0, psi0, the(i,j,int(nz/2))
        write(11,'(2f20.15, 1x, 5f20.15)') x0, y0, u0, v0, p0, divu0, the(i,j,int(nz/2))
      end do
      write(11,*)
    end do
    close(11)
    open(unit=31, file='output/dat/flow'//num//'.dat', access='stream', form='unformatted')
    write(31) u,v,w,p
    close(31)
  end if
  ! =======================================
enddo
! =========================================

write(21) u,v,w,p
end program heat