  !! ----------------------------------------------------------------------- !!
  !>
  !! model size parameter setting
  !<
  !! --
  block

    nx = 400
    ny = 400
    nt = 2000

    title = 'tsunami'

    dx = 500.0
    dy = 500.0
    dt = 1.0

    dxi = 1.0 / dx
    dyi = 1.0 / dy

    xbeg = 0.0
    ybeg = 0.0

  end block
