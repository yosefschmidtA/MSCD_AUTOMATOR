	subroutine augment(e,l,xj,phi,v,nr,r,dl)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension phi(nrmax),phi2(nrmax),v(nrmax),r(nrmax)
	c=137.038d0
	cc=c*c
	c2=cc+cc
	xkappa=-1
	if (dabs(xj).gt.dfloat(l)+0.25d0) xkappa=-l-1
	if (dabs(xj).lt.dfloat(l)-0.25d0) xkappa= l
	do 6110 j=4,nr-3
	if (phi(j).ne.0.d0) then
	  g0=phi(j)
	  ga=(phi(j+1)-phi(j-1))
	  gb=(phi(j+2)-phi(j-2))/2.d0
	  gc=(phi(j+3)-phi(j-3))/3.d0
	  gg=((1.5d0*ga-0.6d0*gb+0.1d0*gc)/(2.d0*dl)+xkappa*g0)/r(j)
	  f0=c*gg/(e-v(j)+c2)
	  phi2(j)=dsqrt(g0*g0+f0*f0)
	  if (g0.lt.0.d0) phi2(j)=-phi2(j)
	else
	  phi2(j)=phi(j)
	endif
 6110   continue
	do 6115 j=1,3
	phi2(j)=phi(j)*phi(4)/phi2(4)
 6115   continue
	do 6120 j=1,nr
	phi(j)=phi2(j)
 6120   continue
	return
	end
