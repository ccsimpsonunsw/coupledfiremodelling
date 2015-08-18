subroutine tkef90(nt,bt,sn,we,hgt,height,u,v,w,tke)
implicit none

integer :: nt,bt,sn,we
integer :: x,y,z,t,top
real :: agl,utmp,vtmp,wtmp,tke_tmp
real,dimension(bt) :: uav,vav,wav
real,dimension(we,sn) :: hgt
real,dimension(we,sn,bt) :: height
real,dimension(we,sn,bt,nt) :: u,v,w
real,dimension(we,sn) :: tke

do x=1,we
do y=1,sn

  ! determine num vert levels below 100 m agl at this (x,y,t)
  agl = 0.
  z = 1
  top = 1
  do while (agl.le.100.)
    agl = height(x,y,z)-hgt(x,y)
    top = z
    z = z+1
  end do

  ! determine average wind at each (x,y,z)
  do z=1,top
    utmp = 0.
    vtmp = 0.
    wtmp = 0.
    do t=1,nt
      utmp = utmp+u(x,y,z,t)
      vtmp = vtmp+v(x,y,z,t)
      wtmp = wtmp+w(x,y,z,t)
    end do
    uav(z) = utmp/(nt*1.)
    vav(z) = vtmp/(nt*1.)
    wav(z) = wtmp/(nt*1.)
  end do

  ! determine average square deviation at each (x,y,z)
  tke_tmp = 0.
  do z=1,top
    utmp = 0.
    vtmp = 0.
    wtmp = 0.
    do t=1,nt
      utmp = utmp+(u(x,y,z,t)-uav(z))*(u(x,y,z,t)-uav(z))
      vtmp = vtmp+(v(x,y,z,t)-vav(z))*(v(x,y,z,t)-vav(z))
      wtmp = wtmp+(w(x,y,z,t)-wav(z))*(w(x,y,z,t)-wav(z))
    end do
    tke_tmp = tke_tmp+0.5*(utmp+vtmp+wtmp)/(nt*1.)
  end do
  tke(x,y) = tke_tmp/(top*1.)

end do
end do


end subroutine tkef90
