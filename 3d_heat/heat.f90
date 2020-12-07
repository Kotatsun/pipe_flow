program heat
!  SMAC (original version by Dr. T. Ushijima)
!  Simplified Marker and Cell method
!  Solving Heat Transfer in flow between parallel walls.
implicit double precision(a-h,o-z)
integer, parameter:: n=80, nx=80, ny=90, nz=90                                    & ! 格子数n
                    ,ix = 80, jx = 90, kx = 90
! integer, parameter:: n=451, nx=401, ny=451, nz=451                                    & ! 格子数n
!                     ,ix = 401, jx = 451, kx = 451
real(8) :: u(0:nx+1,0:ny+1,0:nz+1),v(0:nx+1,0:ny+1,0:nz+1),w(0:nx+1,0:ny+1,0:nz+1)    & ! u,v:速度
          ,the(0:nx+1,0:ny+1,0:nz+1),the0(0:nx+1,0:ny+1,0:nz+1)                         ! the: スカラー(e.g. 温度，濃度)
real(KIND(0.d0)), dimension(ix,jx,kx) :: q                                              ! q: 熱源分布
real(8), parameter:: re=7500.d0                                                        & ! レイノルズ数
                    ,pr=5.6d0                                                         & ! プラントル数 pr (シュミット数)
                    ,r02=36.d0                                                        & ! 血管径
                    ,r12=4.d0                                                         & ! カテーテル径
                    ,q_r=0.92165899e-10                                                  ! 熱源の補正係数

! character*6:: num                                                                     ! ファイル番号用文字列
character filename*128
! integer,parameter:: stepmax=9600000                                                   ! 時間ステップ数
! integer,parameter:: stepmax=1920000                                                   ! 時間ステップ数
! integer,parameter:: stepmax=2000000                                                      ! 時間ステップ数
! integer,parameter:: stepmax=10000000                                                      ! 時間ステップ数
integer,parameter:: stepmax=20000000                                                      ! 時間ステップ数
integer:: step                                                                          ! 時間ステップ
integer :: idf
! grid size 
dy=0.5d0/6.d0 
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
! dt=0.0001d0                                          ! time step
! dt=0.0005d0                                          ! time step
dt=0.0005d0                                          ! time step
! dt=min(dt,0.25*dx) ! 1時間ステップで流体が移流によって飛び出さないようにする
! dt=min(dt,0.2*re*dx*dx) ! 拡散の影響の考慮 (ここではプラントル数の影響は考慮していない)
! write(6,*) 'dt = ',dt
ddt=1.d0/dt

write(6,*) 'dt=', dt
! initial condition
! icont=1 のときは，既存のデータを用いる. それ以外は，初期化する.
icont=0
open(newunit=idf, file='loss_arr_resize.dat',form='unformatted',access='stream')
! open(newunit=idf, file='test.dat',form='unformatted',access='stream')
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
        r2 = dble(k-45)**2+dble(j-45)**2
        if (r2<=r02) then
          u(i,j,k)=2.d0*(1.d0-r2/r02)
        endif
        r3 = dble(k-45)**2+dble(j-49)**2
        if ((r3<=r12 .and. i<=50)) then
          u(i,j,k)=0.d0
        endif
        v(i,j,k)=0.d0
        w(i,j,k)=0.d0
        the(i,j,k)=0.d0
        the0(i,j,k)=0.d0
      enddo
    enddo
  enddo
  write(6,*) 'end of initialization'
end if

!============== time marching ================
do step=1,stepmax
  ! if(mod(step,stepmax/stepmax)==0) write(6,*) 'step =',step
! boundary condition
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
              +ddz2*(the(i,j,k+1)-2.*the(i,j,k)+the(i,j,k-1)) )/(pr*re)

        dq = q(i,j,k)*q_r

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
  if(mod(step, 100000) == 0) then
    ! num='      '
    ! write(num,'(1i6)') step
    ! do j=1,6
    !   if(num(j:j)==' ') num(j:j)='0'
    ! end do
    ! Q = 1e3*4200.d0*216.d0*dx*dy*dz*1e-9*sum(the)*310.d0
    ! time = step*dt*0.06d0
    write(6,*) 'step =',step
    ! write(6,*) 'time = ',step*dt*0.06d0
    write(6,*) 'time = ',step*dt*0.006d0
    write(6,*) 'Q =    ',1e3*4200.d0*216.d0*dx*dy*dz*1e-9*sum(the)*310.d0
    write(6,*) 'w =   ',1e3*4200.d0*216.d0*dx*dy*dz*1e-9*sum(the)*310.d0/(step*dt*0.06d0)
    write(6,*) 'max T =',(maxval(the)+1.d0)*310.d0
    ! write(6,*) 'q =',(sum(q))/(401*451*451)
    ! write(6,*) 'q =',(sum(q))/(80*90*90)

    write (filename, '("output/flow", i9.9, ".txt")') step
    open (11, file=filename, status='replace')
    ! open(unit=11, iostat=ios, file='output/flow'//num//'.txt')
    ! if (ios /= 0) then
    !   write(*,*) 'Failed to open file for output'
    !   stop
    ! end if
    do j=1, ny
      do i=1, nx
        x0=6.d0*dx*(dble(i)-0.5)
        y0=6.d0*dy*(dble(j)-0.5)
        u0=0.5*(u(i,j,int(nz/2))+u(i-1,j,int(nz/2)))
        thermo=(the(i,j,int(nz/2))+1.d0)*310.d0
        ! v0=0.5*(v(i,j,int(nz/2))+v(i,j-1,int(nz/2)))
        ! psi0=0.5*(psi(i,j,int(nz/2))+psi(i-1,j,int(nz/2)))
        ! x座標, y座標, 速度u, 速度v, 圧力P, 連続の式, 流れ関数, 温度
        ! write(11,'(2f20.15, 1x, 6f20.15)') x0, y0, u0, v0, p0, divu0, psi0, the(i,j,int(nz/2))
        ! write(11,'(2f20.15, 1x, 2f20.15)') x0, y0, u0, the(i,j,int(nz/2))
        write(11,'(2f20.15, 1x, 4f20.15)') x0, y0, u0, thermo, step*dt*0.06d0, (maxval(the)+1.d0)*310.d0
      end do
      write(11,*)
    end do
    close(11)
    ! open(unit=31, file='output/dat/flow'//num//'.dat', access='stream', form='unformatted')
    ! write(31) u,v,w
    ! close(31)
    write(21) u,v,w,the,step
  end if
  ! =======================================
enddo
! =========================================

write(21) u,v,w,the,step
end program heat