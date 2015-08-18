;===============================================;
;===============================================;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

inname = "wrfout_bendora_wnw15c_nt.nc"
outname = "wrfout_bendora_wnw15c.nc"

in1 = addfile(inname,"r")
dm1 = getfiledimsizes(in1)
nt = dm1(0)
bt  = dm1(4)
we  = dm1(5)
sn  = dm1(3)
bts = bt+1
sns = sn+1
wes = we+1
wef = dm1(8)
snf = dm1(9)

hgt           = in1->HGT(0,:,:)
zsf           = in1->ZSF(0,:,:)
tign_g        = in1->TIGN_G(nt-1,:,:)
grnhfx        = in1->GRNHFX(:,:,:)
grnqfx        = in1->GRNQFX(:,:,:)
uah           = in1->UAH(:,:,:)
vah           = in1->VAH(:,:,:)
uah2 = wrf_user_unstagger(uah,"X")
vah2 = wrf_user_unstagger(vah,"Y")
delete(uah)
delete(vah)
phb  = in1->PHB(0,:,:,:)
ph   = in1->PH(:,:,:,:)
phb2 = wrf_user_unstagger(phb,"Z")
ph2  = wrf_user_unstagger(ph,"Z")
delete(phb)
delete(ph)
pb   = in1->PB(0,:,:,:)
p    = in1->P(:,:,:,:)
ua   = wrf_user_getvar(in1,"ua",-1)
va   = wrf_user_getvar(in1,"va",-1)
wa   = wrf_user_getvar(in1,"wa",-1)

;===============================================;
;===============================================;

system("/bin/rm -f "+outname)
out = addfile(outname,"c")
setfileoption(out,"DefineMode",True)

dim_names = (/"nt","bt","sn","we","snf","wef"/)
dim_sizes = (/nt,bt,sn,we,snf,wef/)
dim_unlim = (/True,False,False,False,False,False/)
filedimdef(out,dim_names,dim_sizes,dim_unlim)

copy_VarAtts(in1,out)

;2D variables
filevardef(out,"HGT","float",(/"sn","we"/))
filevardef(out,"TIGN_G","float",(/"snf","wef"/))
filevardef(out,"ZSF","float",(/"snf","wef"/))

;3D variables
filevardef(out,"GRNHFX","float",(/"nt","sn","we"/))
filevardef(out,"GRNQFX","float",(/"nt","sn","we"/))
filevardef(out,"UAH","float",(/"nt","sn","we"/))
filevardef(out,"VAH","float",(/"nt","sn","we"/))
filevardef(out,"PHB","float",(/"bt","sn","we"/))
filevardef(out,"PB","float",(/"bt","sn","we"/))

;4D variables
filevardef(out,"P","float",(/"nt","bt","sn","we"/))
filevardef(out,"PH","float",(/"nt","bt","sn","we"/))
filevardef(out,"UA","float",(/"nt","bt","sn","we"/))
filevardef(out,"VA","float",(/"nt","bt","sn","we"/))
filevardef(out,"WA","float",(/"nt","bt","sn","we"/))

;2D variables
out->HGT           = (/hgt/)
out->TIGN_G        = (/tign_g/)
out->ZSF           = (/zsf/)

;3D variables
out->GRNHFX        = (/grnhfx/)
out->GRNQFX        = (/grnqfx/)
out->UAH           = (/uah2/)
out->VAH           = (/vah2/)
out->PHB           = (/phb2/)
out->PB            = (/pb/)

;4D variables
out->P             = (/p/)
out->PH            = (/ph2/)
out->UA            = (/ua/)
out->VA            = (/va/)
out->WA            = (/wa/)

end

;============================================================================;
;============================================================================;
