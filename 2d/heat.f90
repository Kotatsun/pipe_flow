
program heat
!  SMAC (original version by Dr. T. Ushijima)
!  Simplified Marker and Cell method
!  Solving Heat Transfer in flow between parallel walls.
implicit double precision(a-h,o-z)
integer, parameter:: n=20, nx=n*10, ny=n            ! 格子数n
real(8) :: p(0:nx+1,0:ny+1)                       & ! p:圧力
          ,u(0:nx,0:ny+1),v(0:nx+1,0:ny)          & ! u,v:速度
          ,phi(0:nx+1,0:ny+1),divup(1:nx,1:ny)    & ! phi: 補正圧力, divup: 予測速度の発散
          ,up(0:nx,0:ny+1),vp(0:nx+1,0:ny)        & ! up, vp: 予測速度
          ,psi(0:nx,0:ny+1)                       & ! psi: 流れ関数
          ,the(0:nx+1,0:ny+1),the0(0:nx+1,0:ny+1)   ! the: スカラー(e.g. 温度，濃度) 
real(8), parameter:: re=50.d0                     & ! レイノルズ数
                    ,pr=0.7d0                     & ! プラントル数 pr (シュミット数) 
                    ,q_w=1.d0                       ! 壁面熱流束の値
character*6:: num                                   ! ファイル番号用文字列
integer,parameter:: stepmax=20000                   ! 時間ステップ数
integer:: step                                      ! 時間ステップ
! grid size 
dy=1.d0/dble(n)                                     
dx=dy
ddx=1.d0/dx
ddy=1.d0/dy
ddx2=ddx*ddx
ddy2=ddy*ddy
! time step
dt=0.002d0                                          ! time step
dt=min(dt,0.25*dx) ! 1時間ステップで流体が移流によって飛び出さないようにする
dt=min(dt,0.2*re*dx*dx) ! 拡散の影響の考慮 (ここではプラントル数の影響は考慮していない)
! write(6,*) 'dt = ',dt
ddt=1.d0/dt

! initial condition
! icont=1 のときは，既存のデータを用いる. それ以外は，初期化する.
icont=0
if (icont.eq.1) then
  open(unit=9,file='fort.21',form='unformatted', status='unknown')
  read(9) u,v,p
  close(9)
else
  do j=0,ny
      do i=0,nx
        u(i,j)=1.d0
        v(i,j)=0.d0
        p(i,j)=0.d0
        the(i,j)=0.d0
      enddo
  enddo
end if
do i=0,nx
  p(i,ny+1)=0.d0
  u(i,ny+1)=0.d0
  the(i,ny+1)=0.d0
enddo
do j=0,ny
  p(nx+1,j)=0.d0
  v(nx+1,j)=0.d0
  the(nx+1,j)=0.d0
enddo

! boundary condition 1
un=0.d0
uw=1.d0
us=0.d0
ue=0.d0
vn=0.d0
vw=0.d0
vs=0.d0
ve=0.d0

!============== time marching ================
do step=1,stepmax
  if(mod(step,stepmax/1000)==0) write(6,*) 'step =',step
! boundary condition 2
  do j=0,ny
!     right wall (east) or outlet
    v(nx+1,j)=v(nx,j)
    u(nx,j)=u(nx-1,j)
    p(nx+1,j)=0.0
    the(nx+1,j)=2.*the(nx,j)-the(nx-1,j)
!     left wall (west) or inlet
    v(0,j)=vw
    u(0,j)=uw
    p(0,j)=p(1,j)
    the(0,j)=0.d0
  enddo
!         u(0,ny/2+1)=uw
  do i=0,nx
!     lower wall (south)
    u(i,0)=2.*us-u(i,1)
    v(i,0)=vs
    p(i,0)=p(i,1)
    the(i,0)=the(i,1)+q_w*dy !heat flux constant
!     upper wall (north)
    u(i,ny+1)=2.*un-u(i,ny)
    v(i,ny)=vn
    p(i,ny+1)=p(i,ny)
    the(i,ny+1)=the(i,ny)+q_w*dy !heat flux constant
  enddo
! compute the thermal field
  do j=1,ny
    do i=1,nx
      cnvt=(ddx*((the(i,j)+the(i-1,j))*u(i-1,j)     &
            -(the(i,j)+the(i+1,j))*u(i,j))          &
            +ddy*((the(i,j)+the(i,j-1))*v(i,j-1)    &
            -(the(i,j)+the(i,j+1))*v(i,j)))/2.d0    
      
      dift=(ddx2*(the(i+1,j)-2.*the(i,j)+the(i-1,j)) &
            +ddy2*(the(i,j+1)-2.*the(i,j)+the(i,j-1)))/(pr*re)

      the0(i,j)=the(i,j)+dt*(cnvt+dift)
    enddo
  enddo
  do j=1,ny
    do i=1,nx
      the(i,j)=the0(i,j)
    enddo
  enddo


! 予測段，既知の速度・圧力から速度を予測する
! predictor step
! for u_ij
! (up-u)/dt=-dp/dx-duu/dx-duv/dy+(nabla)ˆ2 u
  ! write(6,*)'up'
  do j=1,ny
    do i=1,nx-1
        ! vij=(v(i,j)+v(i,j-1)+v(i+1,j)+v(i+1,j+1))*0.25
        ! cnvu=0.5*ddx*(u(i,j)*(u(i+1,j)-u(i-1,j))
        !      -abs(u(i,j))*(u(i-1,j)-2.*u(i,j)+u(i+1,j)))
        !      +0.5*ddy*(vij*(u(i,j+1)-u(i,j-1))
        !      -abs(vij)*(u(i,j-1)-2.*u(i,j)+u(i,j+1)))
        ! 移流項の離散化，上記コメントアウトされたものは一次精度
      cnvu=ddx*((u(i+1,j)+u(i,j))**2                      &
          -(u(i-1,j)+u(i,j))**2)/4.d0                     &
          +ddy*((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))       &
          -(u(i,j)+u(i,j-1))*(v(i,j-1)+v(i+1,j-1)))/4.d0

      fij=-ddx*(p(i+1,j)-p(i,j))-cnvu                     &
          +ddx2*(u(i+1,j)-2.d0*u(i,j)+u(i-1,j))/re        &
          +ddy2*(u(i,j+1)-2.d0*u(i,j)+u(i,j-1))/re
      
      up(i,j)=u(i,j)+dt*fij
    enddo
  enddo
!   write(6,603) (up(i,j),i=1,nx-1)
  do j=1,ny
    up(0,j)=uw
    up(nx,j)=up(nx-1,j)
  enddo
!   up(0,ny/2+1)=uw*1.01
  do i=1,nx-1
      up(i,0)=2.*us-up(i,1)
      up(i,ny+1)=2.*un-up(i,ny)
  enddo
!   for v_ij
!   (vp-v)/dt=-dp/dy-duv/dx-dvv/dy+(nabla)ˆ2 v
  ! write(6,*)'vp'
  do j=1,ny-1
    do i=1,nx
!          uij=0.25*(u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))
!          cnvv=0.5*ddx*(uij*(v(i+1,j)-v(i-1,j))
!               -abs(uij)*(v(i-1,j)-2.*v(i,j)+v(i+1,j)))
!               +0.5*ddy*(v(i,j)*(v(i,j+1)-v(i,j-1))
!               -abs(v(i,j))*(v(i,j-1)-2.*v(i,j)+v(i,j+1)))
!  cnvu: 移流項の離散化，上記コメントアウトされたものは一次精度
      cnvv=ddx*((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))             &
              -(u(i-1,j+1)+u(i-1,j))*(v(i-1,j)+v(i,j)))/4.d0    &
              +ddy*((v(i,j+1)+v(i,j))**2                        &
              -(v(i,j)+v(i,j-1))**2)/4.d0

      gij=-ddy*(p(i,j+1)-p(i,j))-cnvv                           &
              +ddx2*(v(i+1,j)-2.d0*v(i,j)+v(i-1,j))/re          &
              +ddy2*(v(i,j+1)-2.d0*v(i,j)+v(i,j-1))/re

      vp(i,j)=v(i,j)+dt*gij
    enddo
  enddo

! 境界条件を合わせる
! evaluate continuity
  do i=1,nx
    vp(i,0)=vs
    vp(i,ny)=vn
  enddo
  do j=1,ny-1
    vp(0,j)=vw
    vp(nx+1,j)=vp(nx,j)
  enddo
  ! write(6,*) 'evaluate continuity'
  ic=0
  div=0.0
  do j=1,ny
    do i=1,nx
      divup(i,j)=ddx*(up(i,j)-up(i-1,j)) &
                +ddy*(vp(i,j)-vp(i,j-1))
      div=div+divup(i,j)**2
      ic=ic+1
    enddo
!   write(6,603) (divup(i,j),i=1,nx)
  enddo
  ! write(6,*) sqrt(div/dble(ic))


! Solve the poisson equation (nabla)2 p=(nabla)up/dt by SOR
  ! write(6,*) 'solve the poisson equation for pressure'
! initialization
  do i=0,nx+1
    do j=0,ny+1
      phi(i,j)=0.d0
    enddo
  enddo
  eps=1.D-6
! 最大反復数 maxitr
! 収束させるためには maxitr=nx*ny 程度必要 
! 最大反復数は小さめにとっているので，計算の初めは 収束しないが，定常計算では影響なし.
  maxitr=nx*ny/10
! 緩和係数 alpha=1.5~1.7 程度に設定
  alpha=1.7
! Start iteration
  do iter=1,maxitr
    error=0.d0
!$omp parallel do private(i,j,resid,dphi)
    do j=1,ny
      do i=1,nx
        rhs=ddt*divup(i,j)
        resid=ddx2*(phi(i-1,j)-2.d0*phi(i,j)+phi(i+1,j))  &
              +ddy2*(phi(i,j-1)-2.d0*phi(i,j)+phi(i,j+1)) &
              -rhs
        den=2.d0*(ddx2+ddy2)
        dphi=alpha*resid/den
        error=max(abs(dphi),error)
        phi(i,j)=phi(i,j)+dphi
      enddo
    enddo
!$omp end parallel do
    do j=1,ny
      phi(0,j)=phi(1,j)
      phi(nx+1,j)=0.0
    enddo
    do i=1,nx
      phi(i,0)=phi(i,1)
      phi(i,ny+1)=phi(i,ny)
    enddo
! 収束の判定
    if (error.lt.eps) exit
  enddo
  write(6,*) 'iter=  ', iter, it, error
! correct step
  do j=1,ny
    do i=1,nx-1
      u(i,j)=up(i,j)-dt*ddx*(phi(i+1,j)-phi(i,j))
    enddo  
  enddo
  do j=1,ny-1
    do i=1,nx
      v(i,j)=vp(i,j)-dt*ddy*(phi(i,j+1)-phi(i,j))
    enddo  
  enddo
  do j=1,ny
    do i=1,nx
      p(i,j)=p(i,j)+phi(i,j)
    enddo
  enddo
! check the continuity for n+1 th step
  ! write(6,*) 'check the continuity for n+1 th step'
  ic=0
  div=0.0
  do j=1,ny
    do i=1,nx
        divup(i,j)=ddx*(u(i,j)-u(i-1,j))+ddy*(v(i,j)-v(i,j-1))
        ic=ic+1
        div=div+divup(i,j)**2
    enddo
    ! write(6,603)(divup(i,j),i=1,nx)
! 603  format(20(1X,E10.2))
  enddo
  ! write(6,*) sqrt(div/dble(ic))

! ********************************************* 
! 混合平均温度の計算(i=nx における値)
! *********************************************
  iN=nx
  tm=0.d0
  ! um=0.d0
  do j=1,ny
    tm=tm+0.5d0*(u(iN,j)+u(iN-1,j))*the(iN,j)*dy
    !  um=um+0.5d0*(u(iN,j)+u(iN-1,j))*dy
  enddo
  ! tm=tm/um
  tw=0.5d0*(the(iN,0)+the(iN,1))
  s_Nu=q_w/(tw-tm)
  ! write(6,*) s_Nu

  ! =============== output =================
  if(mod(step, stepmax/50) == 0) then
    ! calculate stream function
    do i=0,nx
      psi(i,0)=0.d0
      do j=1,ny+1
        psi(i,j)=psi(i,j-1)+0.5*dy*(u(i,j-1)+u(i,j))
      enddo
    enddo
    ! output file
    num='      '
    write(num,'(1i6)') step
    do j=1,6
      if(num(j:j)==' ') num(j:j)='0'
    end do
    open(11,file='flow'//num//'.txt')
    do j=0, ny
      do i=0, nx
        x0=dx*(dble(i)-0.5)
        y0=dy*(dble(j)-0.5)
        u0=0.5*(u(i,j)+u(i-1,j))
        v0=0.5*(v(i,j)+v(i,j-1))
        p0=p(i,j)
        divu0=divup(i,j)
        psi0=0.5*(psi(i,j)+psi(i-1,j))
        ! x座標, y座標, 速度u, 速度v, 圧力P, 連続の式, 流れ関数, 温度
        write(11,'(2f20.15, 1x, 6f20.15)') x0, y0, u0, v0, p0, divu0, psi0, the(i,j)
      end do
      write(11,*)
    end do
    close(11)
  end if
  ! =======================================
enddo
! =========================================

! output 結果の出力 c 流れ関数の計算
! do i=0,nx
!   psi(i,0)=0.d0
!   do j=1,ny+1
!     psi(i,j)=psi(i,j-1)+0.5*dy*(u(i,j-1)+u(i,j))
!   enddo
! enddo
! open (10, file='result.csv', status='replace')
! do j=1,ny
!   do i=1,nx
!     x0=dx*(dble(i)-0.5)
!     y0=dy*(dble(j)-0.5)
!     u0=0.5*(u(i,j)+u(i-1,j))
!     v0=0.5*(v(i,j)+v(i,j-1))
!     p0=p(i,j)
!     divu0=divup(i,j)
!     psi0=0.5*(psi(i,j)+psi(i-1,j))
! ! x座標，y座標，速度u,v 圧力P，連続の式，流れ関数，温度 
!     write(10,699) x0,y0,u0,v0,p0,divu0,psi0,the(i,j)
! 699 format (8(1X,E12.5))
!   enddo
!   write(10,*)
! enddo
! do j=1,ny-1
!   do i=1,nx-1
!     x0=dx*dble(i)
!     y0=dy*dble(j)
!     u0=0.5*(u(i,j)+u(i,j+1))
!     v0=0.5*(v(i,j)+v(i+1,j))
!     omega=dx*(v(i+1,j)-v(i,j))-dy*(u(i,j+1)-u(i,j))
!     psi0=0.5*(psi(i,j)+psi(i,j+1))
! ! x 座標，y 座標，速度 u,v，渦度，流れ関数 
!     write(11,699) x0,y0,u0,v0,omega,psi0
!   enddo
!   write(11,699)
! enddo
! 次回の計算のために速度圧力データを保存
write(21) u,v,p
end program heat