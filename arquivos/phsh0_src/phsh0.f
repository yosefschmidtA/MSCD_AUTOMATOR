C  file PSPROG.AB3  Feb. 5, 1995
C  containing programs PHSH0.FOR and PHSH1.FOR
C
C---------------------------------------------------------------------
C  program PHSH0.FOR
C---------------------------------------------------------------------
c
c  there are nr grid points, and distances are in bohr radii...
c
c  r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr)) , i=1,2,3,...nr-1,nr
c
c
c
c  the orbitals are store in phe(), first index goes 1...nr, the
c  second index is the orbital index (i...nel)
c
c  look at the atomic files after printing this out to see everything...
c
c  suffice it to say, that the charge density at radius r(i)
c  in units of electrons per cubic bohr radius is given by 
c
c  sum of j=1...nel, 
c  occ(j)*phe(i,j)*phe(i,j)/(4.d0*3.14159265....*r(i)*r(i))... 
c
c  think of the phe functions as plotting the radial wave-functions
c  as a function of radius...on our logarithmic mesh...
c
c  final note:  
c
c  the Dirac equation is solved for the orbitals, whereas their density
c  is treated by setting phe to the square root of Dirac's F*F+G*G
c  times the sign of G...
c  
c  so we are doing Dirac-Fock, except that we are not treating exchange 
c  exactly, in terms of working with major and minor components of the 
c  orbitals, and the phe's give the CORRECT CHARGE DENSITY...
c
c  the above approximation ought to be very small for valence states,
c  so you need not worry about it...
c
c  the Breit interaction has been neglected altogether...it should not 
c  have a huge effect on the charge density you are concerned with...
C
C author: Eric Shirley
C
	program hartfock
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax)
	dimension no(iorbs),nl(iorbs),xnj(iorbs)
	dimension ev(iorbs),occ(iorbs),is(iorbs)
	dimension ek(iorbs),phe(nrmax,iorbs),orb(nrmax,iorbs)
	dimension njrc(4),vi(nrmax,7),rho(nrmax)
	dimension v(nrmax),q0(nrmax),xm1(nrmax),xm2(nrmax)
	dimension w(33,33),wi(33,33),rhs(33),co(33)
	dimension xint(0:12),vav(11),rint(0:12)
	dimension pin(0:11),sig(0:11),vctab(nrmax,0:3) 
	character*1 ichar
	character*11 jive
	character*60 jive2
	character*70 jive3
c	open(unit=5, file='atorb', status='old')
	rel=0.d0
 10     read (5,20) ichar
 20     format (1a1)
	if (ichar.eq.'d') then
	  write (6,*) 'PLEASE ENTER RELATIVITY FACTOR.  (0=NR, 1=REL.)'
	  read (5,*) rel
	endif
	if (ichar.eq.'x') then
	  write (6,*) 'PLEASE ENTER EXCHANGE CORRELATION (ALPHA).'
	  write (6,*) '0=HARTREE-FOCK, POSITIVE=LDA, NEGATIVE=XALPHA.'
	  read  (5,*) alfa
	endif
	if (ichar.eq.'a') then
	  call abinitio(etot,nst,rel,alfa,nr,r,
     1    dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,
     1    no,nl,xnj,ev,occ,is,ek,orb,iuflag)
	endif
	if (ichar.eq.'i') call initiali(zorig,nr,rmin,rmax,
     1    r,dr,r2,dl,njrc,xntot,nel)
	if (ichar.eq.'q') stop
	if (ichar.eq.'w') then
	  ixflag=1
	  iu=-1
	  ir=0
	  call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
	endif
	if (ichar.eq.'r') then
	  iu=-1
	  ir=1
	  call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
	  call setgrid(nr,rmin,rmax,r,dr,r2,dl)
	endif
	if (ichar.eq.'u') then
	  write (6,*) 'PLEASE ENTER IUFLAG. (0=U, 1=SU, 2=R).'
	  read (5,*) iuflag
	endif
	if (ichar.eq.'c') then
	  write (6,*) 'PLEASE ENTER ALPHA,RS,RP,RD.'
	  read (5,*) corpol,rs,rp,rd
	  do 100 k=1,nr
	  fs=(1.d0-dexp(-(r(k)/rs)**2.d0))**2.d0
	  fp=(1.d0-dexp(-(r(k)/rp)**2.d0))**2.d0
	  fd=(1.d0-dexp(-(r(k)/rd)**2.d0))**2.d0
	  vctab(k,0)=-corpol/2.d0*fs*fs/r(k)**4.d0
	  vctab(k,1)=-corpol/2.d0*fp*fp/r(k)**4.d0
	  vctab(k,2)=-corpol/2.d0*fd*fd/r(k)**4.d0
 100      continue
	endif
	if (ichar.eq.'f') then
	  write (6,*) 'PLEASE ENTER IUNIT,CORPOL'
	  read  (5,*) iunit,corpol
	  write (6,*) 'PLEASE ENTER ILEV,INUM,EOLD'
	  read  (5,*) ilev,inum,eold
	  xl=nl(ilev)
	  if (inum.eq.1) then
	    read (5,*) eav
	  else
	    read (5,*) e1,e2
	    eav=(e1*xl+e2*(xl+1.d0))
     1         /(   xl+    xl+1.d0 )
	  endif
	  if (eav.lt.0.d0) eav=eold+eav
	  if (iunit.eq.2) eav=eav/2.d0
	  if (iunit.eq.3) eav=eav/27.2116d0
	  if (iunit.eq.4) eav=eav*0.000123985d0/27.2116d0
	  sd=dabs(dabs(eav)-dabs(ev(ilev)))
	  rl= 0.d0
	  rh=10.d0
	  sl= 0.d0
	  sh= 0.d0
 300      if (sl*sh.le.0.00000001d0) rc=rl+(rh-rl)/2.d0
	  if (sl*sh.gt.0.00000001d0) rc=rl+(rh-rl)*(sd-sl)/(sh-sl)
	  sc=0.d0
	  do 320 i=1,nr
	  f=(1.d0-dexp(-(r(i)/rc)**2.d0))**2.d0
	  vcpp=corpol/(2.d0*r(i)**4.d0)*f*f
	  sc=sc+dr(i)*phe(i,ilev)*phe(i,ilev)*vcpp
 320      continue
	  if (sc.gt.sd) rl=rc
	  if (sc.gt.sd) sl=sc
	  if (sc.lt.sd) rh=rc
	  if (sc.lt.sd) sh=sc
	  write (6,*) rc,sc
	  if (dabs(sc-sd).gt.0.000001d0) goto 300
	endif
	if (ichar.eq.'p') then
	  call pseudo(etot,nst,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,
     1                phe,orb,njrc,vi,zorig,xntot,nel,
     1                no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
	endif
	if (ichar.eq.'g') then
	  read(5,*) iu
	  read(5,2202) jive
	  read(5,2212) jive2
	  read(5,2222) jive3
 2202     format(1x,1a11)
 2212     format(1x,1a60)
 2222     format(1x,1a70)
	  zizv=dabs(r(nr-1)*vi(nr-1,1))
	  write (iu,2202) jive
	  write (iu,2212) jive2
	  write (iu,2222) jive3
	  write (iu,*) 3,nr,zizv
	  write (iu,*) (r(i),i=1,nr)
	  write (iu,*) 0,(vi(k,1),k=1,nr)
	  write (iu,*) 1,(vi(k,3),k=1,nr)
	  write (iu,*) 2,(vi(k,5),k=1,nr)
	  write (iu,*) (0.d0,k=1,nr)
	  do 500 j=1,nr
	  rh=0.d0
	  do 480 k=1,nel
	  rh=rh+phe(j,k)*phe(j,k)*occ(k)
 480      continue
	  write (iu,*) rh
 500      continue
	endif
	if (ichar.eq.'v') then
	  rold=0.d0
	  open(unit=10)
	  open(unit=11)
	  open(unit=12)
	  do 600 k=1,nr
c         if ((r(k).lt.10.d0).and.((r(k)-rold).gt.0.05d0)) then
	    write (10,*) r(k),vi(k,1)*r(k)
	    write (11,*) r(k),vi(k,3)*r(k)
	    write (12,*) r(k),vi(k,5)*r(k)
	    rold=r(k)
c         endif
 600      continue
	  close(unit=10)
	  close(unit=11)
	  close(unit=12)
	endif
	if (ichar.eq.'V') call fourier(nr,r,dr,r2,vi)
	goto 10
	end
