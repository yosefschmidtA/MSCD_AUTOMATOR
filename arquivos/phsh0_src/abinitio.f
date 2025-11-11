	subroutine abinitio(etot,nst,rel,alfa,nr,r,dr,r2,dl,
     1    phe,njrc,vi,zorig,xntot,nel,no,nl,xnj,
     1    ev,occ,is,ek,orb,iuflag)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax),v(nrmax)
	dimension no(iorbs),nl(iorbs),nm(iorbs),xnj(iorbs)
	dimension ev(iorbs),occ(iorbs),is(iorbs),ek(iorbs)
	dimension phe(nrmax,iorbs),njrc(4),vi(nrmax,7)
	dimension orb(nrmax,iorbs),rpower(nrmax,0:15)
c  this will be good for going up to and including l=3...
	do 10 i=0,7
	xi=i
	do 10 k=1,nr
	rpower(k,i)=r(k)**xi
 10     continue        

c  read in nfc, nel.  refer to the documentation for their meanings.

 168    write (6,*) 'PLEASE ENTER NFC, NEL, RATIO, ETOL, XNUM'
	read (5,*) nfc,nel,ratio,etol,xnum

c  for all of the electrons, read in the quantum numbers.
c  get the total Hartree-active charge.  initialize eigenvalues.

	xntot=0.d0

	write (6,*) 'PLEASE ENTER N L M J S OCC.'

	do 100 i=nfc+1,nel
	read (5,*) no(i),nl(i),nm(i),xnj(i),is(i),occ(i)
	ev(i)=0.d0
	xntot=xntot+occ(i)
	do 100 j=1,nr
	phe(j,i)=0.d0
	orb(j,i)=0.d0
 100    continue

c  initialize the parameters for self-consistency loop.
c  ratio is the mixture of old and new field mixing.

 110    call atsolve(etot,nst,rel,alfa,
     1    eerror,nfc,nr,r,dr,r2,dl,phe,
     1    njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,
     1    ratio,orb,rpower,xnum,etot2,iuflag)

	eerror=eerror*(1.d0-ratio)/ratio
	write (6,112) eerror,etot
112     format (1x,3f14.6)
	if (eerror.gt.etol) goto 110

c  write out information about the atom.

 120    do 130 i=1,nel
	nj=xnj(i)+xnj(i)
	write (6,122) no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
 122    format(1x,2i4,i2,i4,a2,i4,f10.4,f18.6)
 130    continue

	write (6,132) 'TOTAL ENERGY =  ',etot,etot*27.2116d0
 132    format (1x,a16,2f14.6)

	return

	end

