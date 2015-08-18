subroutine trackf90(nt,sn,we,nx,ny,nxy,ba,bb,bc,bd,np,nz2,dxy,zt1,zt2,rt,pi,ui,vi,hfi,po,zo)
implicit none

integer,parameter :: nz1=200
real,parameter :: dp=10.,dz=10.          ! distance thresholds for local minima

integer :: nt,sn,we                      ! input var dimensions
integer :: nx,ny,nxy                     ! number of boxes, tied to box size
integer :: ba,bb,bc,bd                   ! box sizes
integer :: np,nz2
real :: dxy                              ! box grid spacing, for accurate calc of vorticity
real :: zt1,zt2,rt                       ! threshold values used to select whirls
real,dimension(we,sn,nt) :: pi,ui,vi,hfi ! input vars for pressure, u+v winds, and heat flux
integer,dimension(2,np,nt) :: po         ! output, pressure minima locations
integer,dimension(2,nz2,nt) :: zo        ! output, vorticity maxima locations

integer :: swapx,swapy,sa,sb,sc          ! used in sorting algorithm
integer :: xs,xe,ys,ye                   ! box indices
integer :: tmpi,tmpx,tmpy                ! temporary integers
integer :: f,missing,nb,ct1,ct2
integer :: t,i,j,k,x,y                   ! array indices used in loops
real :: tmp,r,swap,zeta,high
logical :: urev,vrev,hflux               ! logic tests for wind reversal and heat flux

integer,dimension(nx,ny) :: pres_x,pres_y
integer,dimension(nxy) :: pres_x1,pres_y1
integer,dimension(nz1) :: zeta_x,zeta_y
integer,dimension(8) :: nbx,nby          ! used in pressure box location calc
real,dimension(nx,ny) :: pres_r
real,dimension(nxy) :: pres_r1
real,dimension(nz1) :: zeta_r

! set other parameters and initial values
high = 1.e6
missing = -999
po = missing
zo = missing
if (bd.gt.bc) then
  print *, "bd>bc so bd set to bc"
  bd = bc
end if

do t=1,nt
  ! find local pressure minima in each box
  do x=1,nx
  do y=1,ny
    tmp = high
    xs = 1+ba*(x-1)
    xe = xs+(ba-1)
    ys = 1+ba*(y-1)
    ye = ys+(ba-1)
    if (xe.gt.we) xe = we
    if (ye.gt.sn) ye = sn
    do i=xs,xe
    do j=ys,ye
      if (pi(i,j,t).lt.tmp) then
        tmp = pi(i,j,t)
        tmpx = i
        tmpy = j
      end if
    end do
    end do
    pres_r(x,y) = tmp
    pres_x(x,y) = tmpx
    pres_y(x,y) = tmpy
  end do
  end do

  ! determine if corner, edge or central box
  f = missing
  do x=1,nx
  do y=1,ny
    nb = 8
    if (x.eq.1.and.y.eq.1) then
      nb = 3
      nbx = (/0,1,1,f,f,f,f,f/)
      nby = (/1,0,1,f,f,f,f,f/)
    else if (x.eq.1.and.y.eq.ny) then
      nb = 3
      nbx = (/0,1,1,f,f,f,f,f/)
      nby = (/-1,0,-1,f,f,f,f,f/)
    else if (x.eq.nx.and.y.eq.1) then
      nb = 3
      nbx = (/0,-1,-1,f,f,f,f,f/)
      nby = (/1,0,1,f,f,f,f,f/)
    else if (x.eq.nx.and.y.eq.ny) then
      nb = 3
      nbx = (/-1,0,-1,f,f,f,f,f/)
      nby = (/0,-1,-1,f,f,f,f,f/)
    end if  
    if (x.eq.1.and.nb.ne.3) then
      nb = 5
      nbx = (/0,1,1,1,0,f,f,f/)
      nby = (/-1,-1,0,1,1,f,f,f/)
    else if (x.eq.nx.and.nb.ne.3) then
      nb = 5
      nbx = (/0,-1,-1,-1,0,f,f,f/)
      nby = (/1,1,0,-1,-1,f,f,f/)
    else if (y.eq.1.and.nb.ne.3) then
      nb = 5
      nbx = (/-1,-1,0,1,1,f,f,f/)
      nby = (/0,1,1,1,0,f,f,f/)
    else if (y.eq.ny.and.nb.ne.3) then
      nb = 5
      nbx = (/-1,-1,0,1,1,f,f,f/)
      nby = (/0,-1,-1,-1,0,f,f,f/)
    end if
    if (nb.eq.8) then
      nbx = (/0,0,-1,-1,-1,1,1,1/)
      nby = (/-1,1,-1,0,1,-1,0,1/)
    end if
    
    ! determine if local minima are within close proximity
    ct1 = 0
    do i=1,nb
      tmpx = pres_x(x,y)-pres_x(x+nbx(i),y+nby(i))
      tmpy = pres_y(x,y)-pres_y(x+nbx(i),y+nby(i))
      if (sqrt(1.*(tmpx*tmpx+tmpy*tmpy)).lt.dp) then
        ct1 = ct1+1
        if (pres_r(x,y).lt.pres_r(x+nbx(i),y+nby(i))) then
          pres_r(x+nbx(i),y+nby(i)) = high
        else
          pres_r(x,y) = high
        end if
      end if
    end do
  end do
  end do

  ! convert 2d to 1d array
  k = 1
  do x=1,nx
  do y=1,ny  
    pres_r1(k) = pres_r(x,y)
    pres_x1(k) = pres_x(x,y)
    pres_y1(k) = pres_y(x,y)
    k = k+1
  end do
  end do

  ! sort all pressure points, find nrp smallest local minima
  sa = 0
  sb = 0
  sc = 1
  do while(sc.ne.0)
    do i=1,nxy-1
      if(pres_r1(i).gt.pres_r1(i+1)) then 
        swap = pres_r1(i)
        swapx = pres_x1(i)
        swapy = pres_y1(i)
        pres_r1(i) = pres_r1(i+1)
        pres_x1(i) = pres_x1(i+1)
        pres_y1(i) = pres_y1(i+1)
        pres_r1(i+1) = swap
        pres_x1(i+1) = swapx
        pres_y1(i+1) = swapy
        sb = sa+1
      end if
    end do
    sc = sb-sa
    sa = sb
  end do

  ! record the lowest pressure minima into output variable
  do i=1,np
    if (t.eq.120) print *, i,pres_x1(i),pres_y1(i)
    po(1,i,t) = pres_x1(i)
    po(2,i,t) = pres_y1(i)
  end do
  
  ! find local vorticity maxima in box centred on pressure minima
  ct1 = 0
  zeta_r = f*1.
  zeta_x = f
  zeta_y = f
  do k=1,np
    if (pres_r1(k).lt.high) then
      if (ct1.lt.nz1) then
        xs = pres_x1(k)-bb
        xe = pres_x1(k)+bb
        ys = pres_y1(k)-bb
        ye = pres_y1(k)+bb
        if (xs.lt.2) xs = 2
        if (ys.lt.2) ys = 2
        if (xe.gt.(we-1)) xe = we-1
        if (ye.gt.(sn-1)) ye = sn-1
        do i=xs,xe
        do j=ys,ye
          zeta = ((vi(i+1,j,t)-vi(i-1,j,t))-(ui(i,j+1,t)-ui(i,j-1,t)))/(2.*dxy)
          if (zeta.gt.zt1) then
            ct1 = ct1+1
            zeta_r(ct1) = zeta
            zeta_x(ct1) = i
            zeta_y(ct1) = j
          end if 
        end do
        end do
      else
        print *, "ct1 exceeding maximum allowed by nz1"
      end if
    end if
  end do

  ! eliminate the close together locations, keep only highest vorticity value
  do i=1,ct1-1
    if (abs(zeta_r(i)-(f*1.)).gt.1.) then
      do j=k,ct1
        tmpx = zeta_x(i)-zeta_x(j)
        tmpy = zeta_y(i)-zeta_y(j)
        if (sqrt(1.*(tmpx*tmpx+tmpy*tmpy)).lt.dz) then
          if (abs(zeta_r(i)).gt.abs(zeta_r(j))) tmpi = j
          if (abs(zeta_r(i)).le.abs(zeta_r(j))) tmpi = i
          zeta_r(tmpi) = f*1.
          zeta_x(tmpi) = f
          zeta_y(tmpi) = f
        end if
      end do
    end if
  end do

  ! check local winds, vorticity and ground heat flux in box around vorticity max centre
  k = 0
  do i=1,ct1
    if (abs(zeta_r(i)-f*1.).gt.1.) then
      tmpx = zeta_x(i)
      tmpy = zeta_y(i)
      urev = .true.
      vrev = .true.
      ! check for wind reversal around central point
      do j=1,bd
        if ((ui(tmpx,tmpy+j,t)*ui(tmpx,tmpy-j,t)).ge.0) urev = .false.
        if ((vi(tmpx+j,tmpy,t)*vi(tmpx-j,tmpy,t)).ge.0) vrev = .false.
      end do
      ! check to see if the box has enough high vorticity points
      if (urev.and.vrev) then
        xs = tmpx-bc
        xe = tmpx+bc
        ys = tmpy-bc
        ye = tmpy+bc
        if ((xs-bc).lt.1) xs = bc+1
        if ((ys-bc).lt.1) ys = bc+1
        if ((xe+bc).gt.we) xe = (we-bc)-1
        if ((ye+bc).gt.sn) ye = (sn-bc)-1
        if (xe.gt.xs) xe=xs+1
        if (ye.gt.ys) ye=ys+1
        ct2 = 0
        hflux = .false.
        do x=xs,xe
        do y=ys,ye
          zeta = ((vi(x+1,y,t)-vi(x-1,y,t))-(ui(x,y+1,t)-ui(x,y-1,t)))/(2.*dxy)
          !if (hfi(x,y,t).gt.0.) hflux = .true.
          hflux = .true.
          if (abs(zeta).ge.zt2) ct2 = ct2+1
        end do
        end do
        r = (ct2*1.)/((xe-xs+1)*(ye-ys+1)*1.)
        ! check to see that a non-zero flux due to fire is nearby
        if (r.gt.rt.and.hflux.and.k.lt.nz2) then
          k = k+1
          zo(1,k,t) = tmpx
          zo(2,k,t) = tmpy
        end if
      end if
    end if
  end do
end do

end subroutine trackf90
