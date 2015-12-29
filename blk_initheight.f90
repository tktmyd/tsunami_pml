  !! ----------------------------------------------------------------------- !!
  !>
  !! Tsunami initial height.
  !! Tsunami height eta(:,:) ( in m-unit ) must be defined here.
  !<
  !! --
  block

    real, parameter :: aa = 16000.
    real, parameter :: bb = 16000.
    integer :: i, j, i0, j0
    real :: hx, hy
    !! ----

    i0 = nx/2
    j0 = ny/2
    eta(:,:) = 0.0
    do j=1, ny
       if( -bb <= (j-j0)*dy  .and. (j-j0) * dy <= bb ) then
          hy = ( 1 + cos( pi * ( j-j0 ) * dy / bb ) ) / 2.0
       else
          hy = 0.0
       end if
       do i=1, nx

          if( -aa <= (i-i0)*dx  .and. (i-i0) * dx <= aa ) then
             hx = ( 1 + cos( pi * ( i-i0 ) * dx / aa ) ) / 2.0
          else
             hx = 0.0
          end if

          eta(i,j) = hx * hy

       end do

    end do

  end block
