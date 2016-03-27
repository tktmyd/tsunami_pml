!! ------------------------------------------------------------------------- !!
!>
!! An example code of an implementation of the absorbing boundary condition
!! for tsunami wave propagation: Liner Dispersive Wave (LDW) case
!!
!! @license
!!   (c) 2015 Takuto Maeda, Hiroaki Tsushima, and Takashi Furumura
!!   License: MIT
!!
!!   In addition, the users are suggested to cite the following accompanying
!!   paper when you use this code.
!!
!!   Maeda, T., H. Tsushima, and T. Furumura,
!!   An effective absorbing boundary condition for linear long-wave and
!!   linear dispersive wave tsunami simulations,
!!   in press, 2016. 
!!
!<
!! ------------------------------------------------------------------------- !!
program ldw

  implicit none

  !! fixed parameters
  integer, parameter :: STDERR  = 0            !< Standard Error
  integer, parameter :: STDOUT  = 6            !< Standard Output
  real,    parameter :: PI      = atan(1.0)*4
  real,    parameter :: g0      = 9.80665      !< gravity
  integer, parameter :: na=20                  !< absorber thickness

  !! control parameters
  integer            :: nx,  ny                !< model size
  integer            :: nt                     !< time step size
  real               :: dx,  dy                !< grid width
  real               :: dt                     !< time step
  character(80)      :: title                  !< model title

  !! derived global variables
  real               :: dxi, dyi               !< inverse grid width
  real               :: xbeg, ybeg             !< model edge
  integer            :: nst                    !< number of stations

  !! arrays
  real,    allocatable :: eta(:,:)             !< tsunami height
  real,    allocatable :: mm(:,:)              !< tsunami vel
  real,    allocatable :: nn(:,:)              !< tsunami vel
  real,    allocatable :: hh(:,:)              !< bathymetry
  real,    allocatable :: hm(:,:)              !< x-averaged bathymetry
  real,    allocatable :: hn(:,:)              !< y-averaged bathymetry
  real,    allocatable :: fh(:,:)              !< land filter
  real,    allocatable :: fm(:,:)              !< land filter
  real,    allocatable :: fn(:,:)              !< land filter
  real,    allocatable :: x(:), y(:)           !< coordinate
  real,    allocatable :: g1x(:,:), g2x(:,:)   !< PML absorber
  real,    allocatable :: g1y(:,:), g2y(:,:)   !< PML absorber
  real,    allocatable :: eta_x(:,:)           !< PML eta(x)
  real,    allocatable :: dtm(:,:)             !< time deriv. of velocity
  real,    allocatable :: dtn(:,:)             !< time deriv. of velocity
  real,    allocatable :: dtm0(:,:)            !< time deriv. of velocity
  real,    allocatable :: dtn0(:,:)            !< time deriv. of velocity
  real,    allocatable :: wav(:,:)             !< waveform at stations
  integer, allocatable :: ist(:), jst(:)       !< station location (i,j)
  real,    allocatable :: xst(:), yst(:)       !< station location (x,y)
  real,    allocatable :: c1x(:,:), c1y(:,:)   !< coefs for LDW calculation
  real,    allocatable :: c2x(:,:), c2y(:,:)   !< coefs for LDW calculation
  real,    allocatable :: c3x(:,:), c3y(:,:)   !< coefs for LDW calculation
  real,    allocatable :: xw(:), yw(:)         !< PML weight function

  !! ----------------------------------------------------------------------- !!

  !! parameter input
  include "blk_param.f90"

  !! ----------------------------------------------------------------------- !!
  !>
  !! memory allocation
  !<
  !! --
  block

    integer :: i, j

    allocate( eta(1:nx,1:ny) )
    allocate( mm (1:nx,1:ny) )
    allocate( nn (1:nx,1:ny) )
    allocate( hh (1:nx,1:ny) )
    allocate( hm (1:nx,1:ny) )
    allocate( hn (1:nx,1:ny) )
    allocate( fh (1:nx,1:ny) )
    allocate( fm (1:nx,1:ny) )
    allocate( fn (1:nx,1:ny) )
    allocate( dtm(0:nx+1,0:ny+1) )
    allocate( dtn(0:nx+1,0:ny+1) )
    allocate( dtm0(0:nx+1,0:ny+1) )
    allocate( dtn0(0:nx+1,0:ny+1) )
    allocate( c1x(nx,ny), c2x(nx,ny), c3x(nx,ny) )
    allocate( c1y(nx,ny), c2y(nx,ny), c3y(nx,ny) )
    allocate( xw(nx), yw(ny) )
    allocate( g1x(2,nx), g2x(2,nx), g1y(2,ny), g2y(2,ny) )
    allocate( eta_x(nx,ny) )


    eta(:,:) = 0.0
    mm (:,:) = 0.0
    nn (:,:) = 0.0
    dtm(0:nx+1,0:ny+1) = 0
    dtn(0:nx+1,0:ny+1) = 0

    allocate( x(nx), y(ny) )
    do i=1, nx
       x(i) = (i-1) * dx + xbeg
    end do
    do j=1, ny
       y(j) = (j-1) * dy + ybeg
    end do

  end block


  include 'blk_initheight.f90'

  include 'blk_bathymetry.f90'

  include 'blk_station.f90'

  !! ----------------------------------------------------------------------- !!
  !>
  !! PML absorber settings
  !<
  !! --
  block

    real :: d1, d2, b1, b2
    integer :: i, j
    real :: d0, b0, R0
    real :: xxb, xxc
    real :: hx
    integer, parameter :: pd = 2
    integer, parameter :: pb = 2
    real,    parameter :: c0 = 500.0 !! assumed phase speed
    !! ----


    !! damping profile constants
    hx = na * dx
    R0 = 10**( - ( log10( real(na) ) - 1 ) / log10( 2.0 )  - 3.0 )
    d0 = - ( 1.0 / (2.0*hx) ) * ( pd +1 ) * c0 * log( R0 )
    b0 = 2.0

    !! initialize eta: assume initial amplitude wihtin PML layer is zero
    eta_x(1:nx,1:ny) = 0.0

    !! initialize
    do i=1, 2
       do j=1, nx
          g1x(i,j) = 1.0
          g2x(i,j) = 1.0
       end do
       do j=1, ny
          g1y(i,j) = 1.0
          g2y(i,j) = 1.0
       end do
    end do

    !! set absorber
    do i=1, na


       !! distance from PML boundary
       xxb = max(hx  - (i    ) * dx,0.)
       xxc = max(hx  - (i+0.5) * dx,0.)

       d1 = d0 * ( xxb / hx )**pd
       b1 = 1.0 + ( b0 - 1.0 ) * ( xxb / hx )**pb

       d2 = d0 * ( xxc / hx )**pd
       b2 = 1.0 + ( b0 - 1.0 ) * ( xxc / hx )**pb

       g1x(1,i) = (b1 - d1/2.0*dt) / (b1 + d1/2*dt)
       g1x(2,i) = (b1 - d2/2.0*dt) / (b1 + d2/2*dt)
       g2x(1,i) =  1.0             / (b1 + d1/2*dt)
       g2x(2,i) =  1.0             / (b1 + d2/2*dt)

       g1y(:,i) = g1x(:,i)
       g2y(:,i) = g2x(:,i)

       ! Opposite Side; Exchange Staggered Grid
       g1x(1,nx-i+1) = g1x(2,i)
       g1x(2,nx-i+1) = g1x(1,i)
       g2x(1,nx-i+1) = g2x(2,i)
       g2x(2,nx-i+1) = g2x(1,i)

       g1y(1,ny-i+1) = g1x(2,i)
       g1y(2,ny-i+1) = g1x(1,i)
       g2y(1,ny-i+1) = g2x(2,i)
       g2y(2,ny-i+1) = g2x(1,i)

    end do

  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! Use 10-m cut off for linear long wave
  !<
  block

    integer :: i, j
    real, parameter :: CUTOFF_DEPTH = 10.0

    do j=1, ny
       do i=1, nx
          if( hh(i,j) < 0.0 ) then
             eta(i,j) = 0.0
          else if ( hh(i,j) < CUTOFF_DEPTH ) then
             hh(i,j) = 10.0
          end if
       end do
    end do

  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! bathymetry averaging for staggered grid calculation
  !<
  block

    integer :: i, j

    ! x-averaged bathymetry
    do j=1, ny
       hm(1,j) = max( hh(1,j), 0.0 )
       do i=2, nx
          hm(i,j) = ( hh(i,j)+hh(i-1,j) ) / 2.0
          if( hh(i,j) <= 0.0  .or.  hh(i-1,j) <= 0.0 ) hm(i,j) = 0.0
       end do
    end do

    ! y-averaged bathymetry
    do i=1, nx
       do j=2, ny
          hn(i,j) = ( hh(i,j) + hh(i,j-1) ) / 2.0
          if( hh(i,j) <= 0.0  .or.  hh(i,j-1) <= 0.0 ) hn(i,j) = 0.0
       end do
       hn(i,1) = max( hh(i,1), 0.0 )
    end do


  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! define land-filter matrix for reflection boundaries
  !<
  block

    integer :: i, j

    fm(:,:) = 1.0
    fn(:,:) = 1.0
    fh(:,:) = 1.0

    do j=1, ny
       do i=1, nx
          if( hm(i,j) < 0.0 ) fm(i,j) = 0.0
          if( hn(i,j) < 0.0 ) fn(i,j) = 0.0
          if( hh(i,j) < 0.0 ) fh(i,j) = 0.0
       end do
    end do

  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! CFL stability condition check
  !<
  block

    if( min(dx,dy)/dt < sqrt( 2*g0*maxval(hh(:,:)) ) ) then
       write(STDERR,*) "ERROR: stability condition is violated: ",  &
            min(dx,dy)/dt,  sqrt( 2*g0*maxval(hh(:,:)) )
       stop
    end if

  end block

  !! ----------------------------------------------------------------------- !!
  !>
  !! Dispersive tsunami coefficieints
  !<
  !! --
  block
    real :: xdiv, ydiv
    integer :: i, j

    do j=1, ny
       do i=1, nx
          xdiv = 3*dx*dx / (3*dx*dx + 2*hm(i,j)*hm(i,j) )
          ydiv = 3*dy*dy / (3*dy*dy + 2*hn(i,j)*hn(i,j) )
          c1x(i,j) = - xdiv  * g0 * hm(i,j)
          c2x(i,j) = xdiv * hm(i,j)*hm(i,j) / ( 3 * dx*dx )
          c3x(i,j) = xdiv * hm(i,j)*hm(i,j) / ( 3 * dx*dy )
          c1y(i,j) = - ydiv  * g0 * hn(i,j)
          c2y(i,j) = ydiv * hn(i,j)*hn(i,j) / ( 3 * dy*dy )
          c3y(i,j) = ydiv * hn(i,j)*hn(i,j) / ( 3 * dx*dy )
       end do
    end do

    ! PML weight
    xw(:) = 1.0
    yw(:) = 1.0
    do i=1, na
       xw(i) = sin( PI * i / na / 2 )

       xw(nx-i+1) = xw(i)
       yw(i)      = xw(i)
       yw(ny-i+1) = yw(i)
    end do

  end block


  !! ----------------------------------------------------------------------- !!
  !>
  !! time stepping
  !<
  block

    integer :: it
    real :: dxeta(nx,ny), dyeta(nx,ny)
    real :: dxm(nx,ny), dyn(nx,ny)
    !! --

    do it= 1, nt

       if( mod(it,100) == 0 ) then
          write(STDERR,'(A,I4.4,A,I4.4)') trim(title) // &
               ": cycle = ", it, " / ", nt
       end if

       !! ------------------------------------------------------------------ !!
       !>
       !! difference of eta with respect to x, y
       !<
       !! --
       block

         integer :: i, j
         !! --

         do j=1, ny
            do i=2, nx
               dxeta(i,j) = ( eta(i,j) - eta(i-1,j) ) * dxi
            end do
            dxeta(1,j) = ( eta(1,j) -        0.0 ) * dxi
         end do

         do i=1, nx
            do j=2, ny
               dyeta(i,j) = ( eta(i,j) - eta(i,j-1) ) * dyi
            end do
            dyeta(i,1) = ( eta(i,1) -        0.0 ) * dyi
         end do

       end block


       !! ------------------------------------------------------------------ !!
       block

         integer, parameter :: ITER_MAX = 1000
         real,    parameter :: TOR = 1e-5
         integer :: k
         integer :: i, j
         real :: err, amax
         real :: mm1, nn1, div
         integer :: icol, jcol
         real :: tmp
         !! --


         !! initial values of time-derivative of vel are
         !! approximated by long-wave theory
         if( it == 1 ) then
            do j=1, ny
               do i=1, nx

                  dtm(i,j) = - g0*hm(i,j)*dxeta(i,j)*dt*fm(i,j)
                  dtn(i,j) = - g0*hn(i,j)*dyeta(i,j)*dt*fn(i,j)

               end do
            end do
         end if


         !! maximum value for conversion check
         amax =max( maxval( abs(dtm(2:nx-1,2:ny-1)) ), &
                    maxval( abs(dtn(2:nx-1,2:ny-1)) ) )
         if( amax < epsilon(1.0 ) ) amax = 1.0



         !! iteration for time-derivative of vel
         do k=1, ITER_MAX

            !! copy previous vel velocities
            do j=1, ny
               do i=1, nx
                  dtm0(i,j) = dtm(i,j)
                  dtn0(i,j) = dtn(i,j)
               end do
            end do


            err = -1e30
            do j=2, ny-1
               do i=2, nx-1
                  dtm(i,j) = c1x(i,j) * dxeta(i,j) &
                           + c2x(i,j) * ( dtm(i+1,j)+dtm(i-1,j) ) &
                           + c3x(i,j) * ( dtn(i,j+1)-dtn(i,j)-dtn(i-1,j+1)+dtn(i-1,j) )
                  err = max( err,  abs( dtm(i,j) - dtm0(i,j) ) / amax )
               end do
            end do

            do j=2, ny-1
               do i=2, nx-1
                  dtn(i,j) = c1y(i,j) *  dyeta(i,j) &
                           + c2y(i,j) * ( dtn(i,j+1)+dtn(i,j-1) ) &
                           + c3y(i,j) * ( dtm(i+1,j)-dtm(i,j)-dtm(i+1,j-1)+dtm(i,j-1))
                  err = max( err,  abs( dtn(i,j) - dtn0(i,j) ) / amax )
               end do
            end do

            if( err < TOR ) exit

         end do

         !! update tsunami velocity using dispersive terms
         do j=2, ny-1
            do i=2, nx-1
               mm(i,j) = g1x(1,i) * mm(i,j) &
                    + g2x(1,i) * ( ( xw(i) - 1 ) * g0 * hm(i,j) * dxeta(i,j) &
                                 + xw(i)   * dtm(i,j)  ) * dt * fm(i,j)
               nn(i,j) = g1y(1,j) * nn(i,j) &
                    + g2y(1,j) * ( ( yw(j) - 1 ) * g0 * hn(i,j) * dyeta(i,j) &
                                 + yw(j)   * dtn(i,j) ) * dt * fn(i,j)
            end do
         end do

       end block



       !! ------------------------------------------------------------------ !!
       !>
       !! difference of velocities with respect to x, y
       !<
       !! --
       block

         integer :: i, j
         !! --

         do j=1, ny
            dxm(nx,j) = ( 0.0 - mm(nx,j)  ) * dxi
            do i=1, nx-1
               dxm(i,j) = ( mm(i+1,j) - mm(i,j) ) * dxi
            end do
         end do

         do i=1, nx
            dyn(i,ny) = ( 0.0 - nn(i,ny) ) * dyi
            do j=1, ny-1
               dyn(i,j) = ( nn(i,j+1) - nn(i,j) ) * dyi
            end do
         end do

       end block



       !! ------------------------------------------------------------------ !!
       block

         integer i, j
         real :: dtxi, dtyi
         !! --

         dtxi = dt * dxi
         dtyi = dt * dyi

         do j=na+1, ny-na
            do i=na+1, nx-na
               eta(i,j) =  eta(i,j) - ( dxm(i,j) +  dyn(i,j) )*dt * fh(i,j)
            end do
         end do

       end block



       !! ------------------------------------------------------------------ !!
       block

         integer :: i, j
         real :: eta_y
         !! --

         do j=1, na
            do i=1, nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)

            end do
         end do

         do j=ny-na+1,ny
            do i=1, nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)
            end do
         end do

         do j=na+1, ny-na
            do i=1, na
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)
            end do
         end do

         do j=na+1, ny-na
            do i=nx-na+1,nx
               eta_y = eta(i,j) - eta_x(i,j)

               eta_x(i,j) = g1x(2,i) * eta_x(i,j) - g2x(2,i) * dxm(i,j) * dt
               eta_y      = g1y(2,j) * eta_y      - g2y(2,j) * dyn(i,j) * dt

               eta(i,j) = ( eta_x(i,j) + eta_y ) * fh(i,j)
            end do
         end do

       end block


       !! ------------------------------------------------------------------ !!
       block

         integer :: i

         do i=1, nst
            wav(i,it) = eta( ist(i), jst(i) )
         end do

       end block



    end do

  end block

  !! ----------------------------------------------------------------------- !!
  block

    integer :: i, it
    character(256) :: fn
    integer, parameter :: io = 100
    !! --

    do i=1, nst
       write(fn,'(A,I3.3,A)') trim(adjustl(title))//'_', i, '.wav'
       open(io, file=fn, action='write', status='unknown' )
       do it=0, nt
          write(io,*) it*dt, wav(i,it)
       end do
       close(io)
    end do

  end block

  !! ----------------------------------------------------------------------- !!
  block

    deallocate( eta )
    deallocate( mm  )
    deallocate( nn  )
    deallocate( hh  )
    deallocate( hm  )
    deallocate( hn  )
    deallocate( fh  )
    deallocate( fm  )
    deallocate( fn  )
    deallocate( x, y )
    deallocate( ist, jst, xst, yst )
    deallocate( wav )
    deallocate( g1x, g2x, g1y, g2y )
    deallocate( eta_x )
    deallocate( c1x, c2x, c3x )
    deallocate( c1y, c2y, c3y )
  end block

end program ldw
!! ------------------------------------------------------------------------- !!
