;===============================================;
;===============================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
external EX01 "./fort_bendora_tke.so"

begin

inname = asciiread("list_bendora_all.txt",-1,"string")
outname = "data_bendora_tke.nc"

nsim = dimsizes(inname)
bts = 0
bte = 20
bt = bte-bts+1
ts = 30
te = 149
nt = te-ts+1

do i=0,nsim-1
  print(i+"")
  
  in = addfile("wrfout_bendora_"+inname(i)+".nc","r")
  hgt = in->HGT(:,:)
  phb = in->PHB(bts:bte,:,:)
  ph = in->PH(ts:te,bts:bte,:,:)
  u = in->UA(ts:te,bts:bte,:,:)
  v = in->VA(ts:te,bts:bte,:,:)
  w = in->WA(ts:te,bts:bte,:,:)

  if (i.eq.0) then
    we = dimsizes(u(0,0,0,:))
    sn = dimsizes(u(0,0,:,0))
    ;bt = dimsizes(u(0,:,0,0))
    print (nt+" "+bt+" "+sn+" "+we)
    height = new((/bt,sn,we/),float)
    tke = new((/sn,we/),float)
    tkeall = new((/nsim,sn,we/),float)
    tkeall@_FillValue = -999
  end if

  ;do t=0,nt-1
  height(:,:,:) = (phb(:,:,:)+ph(ts,:,:,:))/9.81
  ;end do

  tke = 0.
  EX01::tkef90(nt,bt,sn,we,hgt,height,u,v,w,tke)
  tkeall(i,:,:) = tke(:,:)

  delete(in)
  delete(hgt)
  delete(u)
  delete(v)
  delete(w)
end do

;===============================================;
;===============================================;

system("/bin/rm -f "+outname)
out = addfile(outname,"c")
setfileoption(out,"DefineMode",True)

dim_names = (/"nsim","sn","we"/)
dim_sizes = (/nsim,sn,we/)
dim_unlim = (/False,False,False/)
filedimdef(out,dim_names,dim_sizes,dim_unlim)

filevardef(out,"TKE","float",(/"nsim","sn","we"/))

out->TKE = (/tkeall/)

end

;============================================================================;
;============================================================================;
