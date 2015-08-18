;===============================================;
;===============================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

groupname = "bendora"
listname = "list_"+groupname+"_all.txt"
inname = asciiread(listname,-1,"string")
nsim = dimsizes(inname)
outname = "data_"+groupname+"_uvav.nc"

ts = 31
te = 149
nt = te-ts+1

do i=0,nsim-1
  print(i+" ")
  in = addfile("wrfout_"+groupname+"_"+inname(i)+".nc","r")
  if (i.eq.0) then
    dm = getfiledimsizes(in)
    sn  = dm(2)
    we  = dm(3)
    uah_av = new((/nsim,sn,we/),float)
    vah_av = new((/nsim,sn,we/),float)
    wspd_av = new((/nsim,sn,we/),float)
  end if
  uah = in->UAH(ts:te,:,:)
  uah_av(i,:,:) = dim_avg_n(uah,0)
  delete(uah)
  vah = in->VAH(ts:te,:,:)
  vah_av(i,:,:) = dim_avg_n(vah,0)
  delete(vah)
  delete(in)
  wspd_av(i,:,:) = sqrt(uah_av(i,:,:)*uah_av(i,:,:)+vah_av(i,:,:)*vah_av(i,:,:))
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

filevardef(out,"UAH_AV","float",(/"nsim","sn","we"/))
filevardef(out,"VAH_AV","float",(/"nsim","sn","we"/))
filevardef(out,"WSPD_AV","float",(/"nsim","sn","we"/))

out->UAH_AV = (/uah_av/)
out->VAH_AV = (/vah_av/)
out->WSPD_AV = (/wspd_av/)

end

;============================================================================;
;============================================================================;
