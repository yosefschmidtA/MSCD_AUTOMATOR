	subroutine pseudo(etot,nst,rel,alfa,
     1  nr,rmin,rmax,r,dr,r2,dl,
     1  phe,orb,njrc,vi,zorig,xntot,nel,
     1  no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax)
	dimension q0(nrmax),xm1(nrmax),xm2(nrmax),phi(nrmax),v(nrmax)
	dimension njrc(4),njrcdummy(4),vi(nrmax,7),phe(nrmax,iorbs)
	dimension no(iorbs),nl(iorbs),nm(iorbs),xnj(iorbs)
	dimension ek(iorbs),ev(iorbs),occ(iorbs),is(iorbs)
	dimension rpower(nrmax,0:7),orb(nrmax,iorbs)
	dimension vctab(nrmax,0:3)
	do 10 i=1,nel
	nm(i)=0
 10     continue
	njrcdummy(1)=njrc(1)
	njrcdummy(2)=njrc(2)
	njrcdummy(3)=njrc(3)
	njrcdummy(4)=njrc(4)
	write (6,*) 'PLEASE ENTER NP,CORPOL,RNORM'
	read (5,*) np,corpol,rnorm
	xntot=0.d0
	do 1200 i=np,nel
	write (6,42) 'l=',nl(i),' ...'
 42     format(1x,1a2,1i1,1a4)
	lp2=nl(i)+nl(i)+1
	e=ev(i)
	do 57 j=1,nr
	orb(j,i)=orb(j,i)+vctab(j,nl(i))
 57     continue
	idoflag=1
	ns=1
	call setqmm(i,orb,nl(i),ns,idoflag,vi(1,lp2),zeff,zorig,rel,
     1              nr,r,r2,dl,q0,xm1,xm2,njrcdummy,vi)
	do 60 j=1,nr
	orb(j,i)=0.d0
 60     continue
c
c  you can replacing subroutine pseudize with any type of PP generation
c  you want...however, kleinman-bylanderization would take more coding...
c
	call pseudize(i,orb,e,nl(i),xnj(i),no(i),njrc,zeff,vi(1,lp2),
     1                q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
	write (6,*) 'WE HAVE GOT THUS FAR...'
	no(i)=nl(i)+1
	ruse=0.d0
	xkappa=-1.d0
	call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),
     1               zorig,zeff,ev(i),phe(1,i),vi(1,lp2),
     1               q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
	write (6,*) nl(i),ev(i)
	xntot=xntot+occ(i)
	if (lp2.eq.1) goto 1200
	do 1170 j=1,nr
	vi(j,lp2-1)=vi(j,lp2)
 1170   continue
 1200   continue
	write (6,*) 'everything is pseudized'
	do 1210 i=np,nel
	inew=1+i-np
	no (inew)=no (i)
	nl (inew)=nl (i)
	nm (inew)=nm (i)
	xnj(inew)=xnj(i)
	is (inew)=1
	ev (inew)=ev (i)
	occ(inew)=occ(i)
	do 1210 j=1,nr
	phe(j,inew)=phe(j,i)
 1210   continue
	nel=1+nel-np
	do 1212 i=0,7
	xi=i
	do 1212 k=1,nr
	rpower(k,i)=r(k)**xi
 1212   continue
	write (6,*) 'everything is scaled down...ready for unscreening'
	xnum=100.d0
	ratio=1.d0
	call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,
     1              xntot,phe,ratio,orb,occ,is,
     1              nel,nl,nm,no,xnj,rpower,xnum,etot2,iuflag)
	write (6,*) 'screening effects in pseudo atom computed...'
	do 1250 k=1,nel
	lp2=nl(k)+nl(k)+1
	do 1250 j=1,nr
		      vi(j,lp2  )=vi(j,lp2  )-orb(j,k)
	if (lp2.gt.1) vi(j,lp2-1)=vi(j,lp2-1)-orb(j,k)
 1250   continue
	write (6,*) 'we got past the unscreening...'
	do 1300 j=1,nr
	vl =     (     vi(j,2)+2.d0*vi(j,3))/3.d0
	vso=2.d0*(     vi(j,3)-     vi(j,2))/3.d0
	vi(j,2)=vso
	vi(j,3)=vl
	vl =     (2.d0*vi(j,4)+3.d0*vi(j,5))/5.d0
	vso=2.d0*(     vi(j,5)-     vi(j,4))/5.d0
	vi(j,4)=vso
	vi(j,5)=vl
 2222   format (5f8.4)
	vl =     (3.d0*vi(j,6)+4.d0*vi(j,7))/7.d0
	vso=2.d0*(     vi(j,7)-     vi(j,6))/7.d0
	vi(j,6)=vso
	vi(j,7)=vl
 1300   continue
	rel=0.d0
	write (6,*) 'we got past the spin-orbit jazz'
	izuse=dabs(vi(nr-2,1)*r(nr-2))+0.5d0
	 zuse=izuse
	do 2100 k=1,nr
	if (r(k).gt.rnorm) then
	  videal=-zuse/r(k)-corpol/(2.d0*r(k)**4.d0)
	  vi(k,1)=videal
	  vi(k,3)=videal
	  vi(k,5)=videal
	  vi(k,7)=videal
	  vi(k,2)=0.d0
	  vi(k,4)=0.d0
	  vi(k,6)=0.d0
	endif
 2100   continue
	write (6,*) 'we got to the end'
	return
	end
