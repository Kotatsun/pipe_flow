program heat
!  SMAC (original version by Dr. T. Ushijima)
!  Simplified Marker and Cell method
!  Solving Heat Transfer in flow between parallel walls.
implicit double precision(a-h,o-z)
integer, parameter:: n=40, nx=40, ny=60, nz=60                                          & ! 格子数n
                    ,ix = 40, jx = 60, kx = 60
real(8) :: u(0:nx+1,0:ny+1,0:nz+1),v(0:nx+1,0:ny+1,0:nz+1),w(0:nx+1,0:ny+1,0:nz+1)      & ! u,v:速度
          ,the(0:nx+1,0:ny+1,0:nz+1),the0(0:nx+1,0:ny+1,0:nz+1)                         & ! the: スカラー(e.g. 温度，濃度)
          ,pr(0:nx+1,0:ny+1,0:nz+1),qr(0:nx+1,0:ny+1,0:nz+1)                              ! pr: プラントル数 qr: 熱源の補正係数
real(KIND(0.d0)), dimension(ix,jx,kx) :: q                                                ! q: 熱源分布
real(8), parameter:: mu  =4.7e-3                                                        & ! 粘度
                    ,Cp_b=3.617e3                                                       & ! 比熱 [kJ/K/kg]（血液）
                    ,Cp_w=3.306e3                                                       & ! 比熱 [kJ/K/kg]（血管壁）
                    ,Cp_t=2.3564e3                                                      & ! 比熱 [kJ/K/kg]（周囲組織）
                    ,k_b =0.52d0                                                        & ! 熱伝導率 [W/K/m]（血液）
                    ,k_w =0.46d0                                                        & ! 熱伝導率 [W/K/m]（血管壁）
                    ,k_t =0.273d0                                                       & ! 熱伝導率 [W/K/m]（周囲組織）
                    ,U_r =0.059d0                                                       & ! 基準速度 [m/s]
                    ,T_in=310.d0                                                        & ! 基準温度 [K]
                    ,rho =1.05e3                                                        & ! 密度 [kg/m^3]
                    ,D   =6.0e-3                                                        & ! 基準長さ（血管内径） [m]
                    ,D_w =8.0e-3                                                        & ! 基準長さ（血管外径） [m]
                    ,d_c =3.d0                                                            ! カテーテル径 [m]
                    ! ,q_r=46.08e-12                                                 
                    ! ,q_r=0.92165899e-9                                                  ! 熱源の補正係数

character filename*128
integer,parameter:: stepmax=590000                                                        ! 時間ステップ数
integer:: step                                                                            ! 時間ステップ
integer :: idf

re=rho*U_r*D/mu                                                                           ! Reynolds Number re
dy=(0.5e-3)/D                                                                             ! grid size 
! dy=0.1d0/6.d0
! dy=0.1d0                                  
dx=dy
dz=dy
ddx=1.d0/dx
ddy=1.d0/dy
ddz=1.d0/dz
ddx2=ddx*ddx
ddy2=ddy*ddy
ddz2=ddz*ddz
! time step
dt=0.001d0                                          ! time step
! dt=0.0005d0                                          ! time step
! dt=min(dt,0.25*dx) ! 1時間ステップで流体が移流によって飛び出さないようにする
! dt=min(dt,0.2*re*dx*dx) ! 拡散の影響の考慮 (ここではプラントル数の影響は考慮していない)
ddt=1.d0/dt
write(6,*) 'dt=', dt
write(6,*) 're=', re

! initial condition
! icont=1 : use the input data to initialize
icont=0
! read the heat field condition
open(newunit=idf, file='../resources/helical_5.5mm_n8_resize.dat',form='unformatted',access='stream')
  read(idf) q
  close(idf)
if (icont.eq.1) then
  open(unit=9, file='fort.21', form='unformatted', access="stream")
  read(9) u,v,w,the,step
  close(9)
else
  write(6,*) 'initialization'
  do k=0,nz+1
    do j=0,ny+1
      do i=0,nx+1
        u(i,j,k)=0.d0
        ! r2 = sqrt(dble(k-nz/2.d0)**2+dble(j-ny/2.d0)**2)*(dy*D)
        ! r3 = sqrt(dble(k-nz/2.d0)**2+dble(j-(ny/2.d0+D_w/(2.d0*dy*D)))**2)*(dy*D)
        ! if (r2<=D/2.d0) then
        !   u(i,j,k)=2.d0*(1.d0-(r2/(D/2.d0))**2)
        ! endif
        ! if (r3<=D_w/2.d0) then
        !   u(i,j,k)=0.d0
        ! endif
        v(i,j,k)=0.d0
        w(i,j,k)=0.d0
        the(i,j,k)=0.d0
        the0(i,j,k)=0.d0
        
        ! biological parameter
        pr(i,j,k)=mu*Cp_t/k_t
        qr(i,j,k)=D/(rho*Cp_t*T_in*U_r)*20.d0
        if (r2<=D/2.d0) then
          ! Blood
          pr(i,j,k)=mu*Cp_b/k_b
          qr(i,j,k)=D/(rho*Cp_b*T_in*U_r)*20.d0
        elseif (D/2.d0<r2 .and. r2<=D_w/2.d0) then
          ! Blood vessel wall
          pr(i,j,k)=mu*Cp_w/k_w
          qr(i,j,k)=D/(rho*Cp_w*T_in*U_r)*20.d0
        endif
      enddo
    enddo
  enddo
  write(6,*) 'end of initialization'
end if
write(6,*) 'pr=', pr(1,1,1)

!============== time marching ================
do step=1,stepmax
  ! if(mod(step,stepmax/stepmax)==0) write(6,*) 'step =',step
!     boundary condition
!     right wall (east) or outlet
  do j=0,ny
    do k=0,nz
      the(nx+1,j,k)=2.*the(nx,j,k)-the(nx-1,j,k)
!     left wall (west) or inlet
      the(0,j,k)=0.d0
    enddo
  enddo
!     lower wall (south)
  do i=0,nx
    do k=0,nz
      the(i,0,k)=the(i,1,k) 
!     upper wall (north)
      the(i,ny+1,k)=the(i,ny,k) 
    enddo
  enddo
!     forward wall (forward)
  do i=0,nx
    do j=0,ny
      the(i,j,0)=the(i,j,1) 
!     backward wall (backward)
      the(i,j,nz+1)=the(i,j,nz) 
    enddo
  enddo

! compute the thermal field
!$OMP parallel do private(i,j,k,cnvt,dift,dq) shared(ddx,ddy,ddz,ddx2,ddy2,ddz2,dt,the,the0,u,v,w,q)
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
              +ddz2*(the(i,j,k+1)-2.*the(i,j,k)+the(i,j,k-1)) )/(pr(i,j,k)*re)

        dq = q(i,j,k)*qr(i,j,k)

        the0(i,j,k)=the(i,j,k)+dt*(cnvt+dift+dq)
      enddo
    enddo
  enddo
!$OMP end parallel do
  do k=1,nz
    do j=1,ny
      do i=1,nx
        the(i,j,k)=the0(i,j,k)
      enddo
    enddo
  enddo
  ! write(6,*) 'end of computing thermal field'
  ! =====================================
  ! =============== output =================
  ! if(mod(step, 5000) == 0) then
  if(mod(step, 10000) == 0) then
    write(6,*) 'step =',step
    ! write(6,*) 'time = ',step*dt*0.06d0
    write(6,*) 'time = ',step*dt*D/U_r,' [s]'
    write(6,*) 'Q =    ',rho*Cp_b*(D**3)*dx*dy*dz*sum(the)*T_in,' [J]'
    write(6,*) 'w =    ',rho*Cp_b*(D**3)*dx*dy*dz*sum(the)*T_in/(step*dt*D/U_r),' [W]'
    write(6,*) 'max T =',(maxval(the)+1.d0)*T_in,' [K]'
    ! write(6,*) 'q =',(sum(q))/(80*90*90)

    write (filename, '("output/flow", i9.9, ".txt")') step
    open (11, file=filename, status='replace')
    do j=1, ny
      do i=1, nx
        x0=1e3*D*dx*(dble(i)-0.5)
        y0=1e3*D*dy*(dble(j)-0.5)
        u0=0.5*(u(i,j,int(nz/2))+u(i-1,j,int(nz/2)))
        thermo=(the(i,j,int(nz/2))+1.d0)*T_in
        ! x座標, y座標, 速度u, 温度, 時間ステップ, 最大温度
        write(11,'(2f20.15, 1x, 4f20.15)') x0, y0, u0, thermo, step*dt*D/U_r, (maxval(the)+1.d0)*T_in
      end do
      write(11,*)
    end do
    close(11)
    write(21) u,v,w,the,step
  end if
  ! =======================================
enddo
! =========================================

write(21) u,v,w,the,step
end program heat