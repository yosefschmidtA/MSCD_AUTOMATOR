	subroutine setgrid(nr,rmin,rmax,r,dr,r2,dl)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax)
	ratio=rmax/rmin
	dl=dlog(ratio)/dfloat(nr)
	xratio=dexp(dl)
	xr1=dsqrt(xratio)-dsqrt(1.d0/xratio)
	do 2010 i=1,nr
	r(i)=rmin*xratio**dfloat(i)
	dr(i)=r(i)*xr1
	r2(i)=r(i)*r(i)
 2010   continue
	return
	end
