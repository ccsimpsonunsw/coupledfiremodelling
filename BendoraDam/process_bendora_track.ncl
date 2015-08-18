;===============================================;
;===============================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
external EX01 "./fort_bendora_track.so"

begin

inname = asciiread("list_bendora_all.txt",-1,"string")
outname = "data_bendora_track.nc"

nsim = dimsizes(inname)
ts = 30
te = 149
nt = te-ts+1

do i=0,nsim-1
  print(i+"")
  
  in = addfile("wrfout_bendora_"+inname(i)+".nc","r")
  p = in->P(ts:te,0,:,:)
  u = in->UAH(ts:te,:,:)
  v = in->VAH(ts:te,:,:)
  hf = in->GRNHFX(ts:te,:,:)
  ;printVarSummary(u)
  ;printVarSummary(hf)

  if (i.eq.0) then
    we = dimsizes(p(0,0,:))
    sn = dimsizes(p(0,:,0))
    np = 80
    nz = 40
    nx = 15
    ny = 15
    nxy = nx*ny
    ba = 20
    bb = 5
    bc = 3
    bd = 2
    dxy = 30.
    zt1 = 0.05
    zt2 = 0.02
    rt = 0.50
    nvar = 2
    po = new((/nt,np,nvar/),integer)
    zo = new((/nt,nz,nvar/),integer)
    po@_FillValue = -999
    zo@_FillValue = -999
    trackp = new((/nsim,nt,np,nvar/),integer)
    trackz = new((/nsim,nt,nz,nvar/),integer)
    trackp@_FillValue = -999
    trackz@_FillValue = -999
  end if

  po = 0
  zo = 0
  EX01::trackf90(nt,sn,we,nx,ny,nxy,ba,bb,bc,bd,np,nz,dxy,zt1,zt2,rt,p,u,v,hf,po,zo)
  trackp(i,:,:,:) = po(:,:,:)
  trackz(i,:,:,:) = zo(:,:,:)

  delete(in)
  delete(p)
  delete(u)
  delete(v)
end do

;===============================================;
;===============================================;

system("/bin/rm -f "+outname)
out = addfile(outname,"c")
setfileoption(out,"DefineMode",True)

dim_names = (/"nsim","nt","np","nz","nvar"/)
dim_sizes = (/nsim,nt,np,nz,nvar/)
dim_unlim = (/False,False,False,False,False/)
filedimdef(out,dim_names,dim_sizes,dim_unlim)

filevardef(out,"TRACKP","float",(/"nsim","nt","np","nvar"/))
filevardef(out,"TRACKZ","float",(/"nsim","nt","nz","nvar"/))

out->TRACKP = (/trackp/)
out->TRACKZ = (/trackz/)

end

;============================================================================;
;============================================================================;
