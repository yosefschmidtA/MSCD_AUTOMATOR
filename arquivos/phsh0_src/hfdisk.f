	subroutine hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1                    zorig,xntot,ixflag,nel,
     1                    no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension no(iorbs),nl(iorbs),xnj(iorbs),is(iorbs)
	dimension ev(iorbs),ek(iorbs),occ(iorbs),r(nrmax)
	dimension phe(nrmax,iorbs),orb(nrmax,iorbs)
	dimension njrc(4),vi(nrmax,7),rho(nrmax)
	character*10 filename
	pi=4.*atan(1.)
	pi4=4.*pi
	rden=3.0
	if (iu.lt.0) then
	  write (6,*) 'PLEASE ENTER FULL FILENAME.'
	  read (5,52) filename
 52       format (a10)
	  iu1=1
	  open (unit=iu1,status='unknown',file=filename)
	else
	  iu1=iu
	  open (unit=iu1,status='unknown')
	endif
c
c define the logarithmic grid
c
	do 5 i=1,nr
	   r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
5       continue
c
c obtain the charge density on the logarithmic grid
c
	do 20 i=1,nr
	  rho(i)=.0
	  do 10 ii=1,nel
	      rho(i)=rho(i)+occ(ii)*phe(i,ii)**2
10        continue
20      continue
c
c write output file
c
	iprint=0
	write(iu1,550)iprint
550     format('RELA'/'RELAT. ATOMIC CHARGE DENSITY'/I1)
	write(iu1,54) rmin,rmax,nr,zorig
54      format (d15.8,d15.8,i5,f5.2)
	 nden=nr*(log(rden/rmin)/log(rmax/rmin))
	 write(iu1,56) (rho(j),j=1,nr)
56      format (f15.11)
	close (unit=iu1)
	return
	end
