	subroutine fitx0(i,orb,rcut,njrc,e,l,xj,n,jrt,xideal,phi,
     1                   zeff,v,q0,xm1,xm2,nr,r,dr,r2,dl,rel,factor)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax)
	dimension njrc(4),orb(nrmax,iorbs)
	dimension q0(nrmax),xm1(nrmax),xm2(nrmax)
	dimension phi(nrmax),v(nrmax)
	vl=-1000000.d0
	vh=+1000000.d0
 115    idoflag=2
	ns=1
	xkappa=-1.d0
	call setqmm(i,orb,l,ns,idoflag,v,zeff,dummy,rel,
     1              nr,r,r2,dl,q0,xm1,xm2,njrc,dummy)
	call integ(e,l,xkappa,n,nn,jrt,ief,xactual,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
	if (nn.ne.0) then
	  vl=v(1)
	  xla=1.d0
	else
	  if (xactual.gt.xideal) then
	    vh=v(1)
	  else
	    vl=v(1)
	  endif
	  xerror=xideal-xactual
	  if (dabs(xerror).lt.0.000000001d0) return
	  dxdla=0.d0
	  do 120 ii=1,jrt
	  dxdla=dxdla+dr(ii)*phi(ii)*phi(ii)*hb(r(ii)/rcut,factor)
 120      continue
	  dxdla=2.d0*dxdla/(phi(jrt)*phi(jrt))
	  xla=xerror/dxdla
	endif
	vmaybe=v(1)+xla
	if ((vmaybe.gt.vh).or.(vmaybe.lt.vl)) xla=(vl+vh)/2.d0-v(1)
	do 130 ii=1,jrt-1
	v(ii)=v(ii)+xla*hb(r(ii)/rcut,factor)
 130    continue
	goto 115
	end
