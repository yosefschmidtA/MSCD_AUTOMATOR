	subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,
     1    xntot,nel)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax),njrc(4)
	write (6,*) 'ENTER Z, NR'
	read (5,*) zorig,nr
	rmin=0.0001d0/zorig
	rmax=800.d0/dsqrt(zorig)
	call setgrid(nr,rmin,rmax,r,dr,r2,dl)
	do 5 j=1,4
	njrc(j)=0
 5      continue
	xntot=0.d0
	nel=0
	return
	end
