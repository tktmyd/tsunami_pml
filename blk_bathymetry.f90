  !! ----------------------------------------------------------------------- !!
  !>
  !! Bathymetry depth hh(:,:) ( in m-unit ) must be filled here. 
  !<
  block

    integer :: i, j
    !! --

    do i=1, nx
       do j=1, 60
          hh(i,ny-j+1) = 500.0
       end do
       do j=61, 400
          hh(i,ny-j+1) = 500 + (j-60) * dy * tan( 2.0 * PI / 180 )
       end do
    end do

  end block
