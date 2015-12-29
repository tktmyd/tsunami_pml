  !! ----------------------------------------------------------------------- !!
  !>
  !! Set station location.
  !! nst, xst(:), yst(:), ist(:), jst(:) should be set here with memory alloc.
  !<
  block

    integer :: i

    nst = 5

    allocate( ist(nst), jst(nst), xst(nst), yst(nst) )
    allocate( wav(nst,0:nt) )

    xst(:) = (/ 50, 150, 150,  50, 100 /) * 1000 ! m -> km
    yst(:) = (/ 50,  50, 150, 150, 150 /) * 1000 ! m -> km

    !! choose nearest station grid
    !!   xbeg + (ist-1)*dx < xst <= xbeg + (ist) * dx
    !!   ybeg + (jst-1)*dy < yst <= ybeg + (jst) * dy
    do i=1, nst
       ist(i) = ceiling( ( xst(i) - xbeg ) / dx + 0.5 )
       jst(i) = ceiling( ( yst(i) - ybeg ) / dy + 0.5 )
    end do

    wav(:,:) = 0.0
    
  end block
  
