	subroutine fourier(nr,r,dr,r2,vi)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax),vi(nrmax,7)
	dimension a(nrmax),v1(nrmax),v2(nrmax)
	do 350 l=0,2
	lp2=l+l+1
	dl=dlog(r(2)/r(1))
	dl1=12.d0*dl
	dl2=12.d0*dl*dl
	do 220 i=1,nr
	a(i)=r(i)*vi(i,lp2)
 220    continue
	do 225 i=3,nr-2
	al =(-(a(i+2)-a(i-2))+ 8.d0*(a(i+1)-a(i-1))           )/dl1
c       all=(-(a(i+2)+a(i-2))+16.d0*(a(i+1)+a(i-1))-30.d0*a(i))/dl2
	ar =al/r(i)
	v1(i)=ar
 225    continue
	open (unit=20+l,status='unknown')
	do 300 ii=1,200
	q=dfloat(ii)/10.d0
	vq=0.d0
	do 250 i=3,nr-2
	vq=vq+dr(i)*dcos(q*r(i))*v1(i)
 250    continue
	write (20+l,*) q,vq
 300    continue
	close(unit=1)
 350    continue
	return
	end
