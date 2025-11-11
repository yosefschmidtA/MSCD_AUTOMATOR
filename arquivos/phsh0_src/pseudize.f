	subroutine pseudize(i,orb,ev,l,xj,n,njrc,zeff,v,
     1                      q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
	implicit real*8 (a-h,o-z)
	parameter (iorbs=33,iside=600)
	parameter (io2=iorbs*(iorbs+1)/2)
	parameter (ijive=io2*(io2+1)/2)
	parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
	dimension r(nrmax),dr(nrmax),r2(nrmax)
	dimension phi(nrmax),v(nrmax),rf(3),vf(3)
	dimension phi0(nrmax),yl(nrmax),vraw(nrmax)
	dimension q0(nrmax),xm1(nrmax),xm2(nrmax)
	dimension njrc(4),orb(nrmax,iorbs)
	lp=l+1
	xkappa=-1.d0
	istop=nr
 40     istop=istop-1
	if (ev.le.q0(istop)) goto 40
	call integ(ev,l,xkappa,n,nn,istop,ief,xdummy,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
 50     write (6,*) 'PLEASE ENTER THE CUTOFF RADIUS, AND FACTOR.'
	read (5,*) rcut,factor
	if (rcut.lt.0.d0) then
	  xnodefrac=-rcut
	  j=istop
 55       j=j-1
	  if (phi(j-1)/phi(j).gt.1.d0) goto 55
	  if (n.gt.l+1) then
	    k=j
 60         k=k-1
	    if (phi(k-1)/phi(k).gt.0.d0) goto 60
	  else
	    k=1
	  endif
	  rcut=r(k)+xnodefrac*(r(j)-r(k))
	endif
	jrc=1.d0+dfloat(nr-1)*dlog(rcut /rmin)/dlog(rmax/rmin)
	rcut=r(jrc)
	rtest=2.d0*rcut
	jrt=1.d0+dfloat(nr-1)*dlog(rtest/rmin)/dlog(rmax/rmin)
	njrc(lp)=jrt
	rtest=r(jrt)
	switch=phi(jrt)/dabs(phi(jrt))
	write (6,92) 'RCUTOFF = ',rcut,'  JRC = ',jrc
	write (6,92) 'RTEST   = ',rtest,  '  JRT = ',jrt
 92     format (1x,1a10,1f8.4,1a8,1i5)
 94     format (1x,2d15.8)
	call integ(ev,l,xkappa,n,nn,jrt,ief,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
	do 8000 ii=1,jrt
	phi(ii)=phi(ii)/phi(jrt)
 8000   continue
	xn00=0.d0
	do 8010 ii=1,jrt-1
	xn00=xn00+dr( ii)*phi( ii)*phi( ii)
 8010   continue
	xn00=xn00+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
	de=0.0001d0
	ee=ev+de/2.d0
	call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
	ee=ev-de/2.d0
	call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
	c00=(xm-xp)/(2.d0*de)
	write (6,94)  c00,x00
	write (6,94) xn00
	ruse=0.d0
	v0=v(jrc)
	dvdl  =(8.d0*(v(jrc+1)-v(jrc-1))-(v(jrc+2)-v(jrc-2)))
     1         /(12.d0*dl)
	ddvdll=(16.d0*(v(jrc+1)+v(jrc-1))
     1         -30.d0*v(jrc)-v(jrc+2)-v(jrc-2))
     1         /(12.d0*dl*dl)
	dldr=1.d0/r(jrc)
	ddldrr=-1.d0/r2(jrc)
	v1=dvdl*dldr
	v2=dvdl*ddldrr+ddvdll*dldr*dldr
	b4=(v2*rcut-v1)/(8.d0*rcut**3.d0)
	b2=(v1-4.d0*b4*rcut**3.d0)/(2.d0*rcut)
	b0=v0-b4*rcut**4.d0-b2*rcut**2.d0
	do 110 ii=1,jrc
	rr=r(ii)
	v(ii)=b0+b2*rr**2.d0+b4*rr**4.d0
 110    continue
	call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
 180    do 200 ii=1,jrt
	phi0(ii)=phi(ii)
	vraw(ii)=v(ii)
 200    continue
 210    xi0=0.d0
	xi1=0.d0
	xi2=0.d0
	do 220 ii=1,jrt
	f=hb(r(ii)/rcut,factor)
	ph2=dr(ii)*phi0(ii)*phi0(ii)
	xi0=xi0+ph2
	if (ii.le.jrt) then
	  xi1=xi1+ph2*f
	  xi2=xi2+ph2*f*f
	endif
 220    continue
	ph2=phi0(jrt)*phi0(jrt)
	xi0=xi0/ph2
	xi1=xi1/ph2
	xi2=xi2/ph2
	quant=xi1*xi1+xi2*(c00-xi0)
	if (quant.gt.0.d0) then
	  deltal=(dsqrt(xi1*xi1+xi2*(c00-xi0))-xi1)/xi2
	else
	  deltal=(c00-xi0)/(2.d0*xi1)
	endif
	write (6,222) 'DELTAL = ',deltal
 222    format (1x,1a9,1f11.8)
 225    do 230 ii=1,jrt
	yl (ii)=phi0(ii)*hb(r(ii)/rcut,factor)
	phi(ii)=phi0(ii)+deltal*yl(ii)
	if (phi(ii).lt.0.d0) then
	  write (6,*) 'BIG TROUBLE!!! CROSS AXIS!!!'
	  stop
	endif
 230    continue
	do 300 ii=1,jrt-1
	if ((phi(ii).eq.0.).or.(yl(ii).eq.0.)) goto 1170
	jj=ii
	if (ii.eq.1) jj=2
	do 240 j=jj-1,jj+1
	rf(2+j-jj)=r(j)
	vf(2+j-jj)=hb(r(j)/rcut,factor)
 240    continue
	call parabreg(f,fp,fpp,rf,vf)
	do 242 j=jj-1,jj+1
	vf(2+j-jj)=phi0(j)
 242    continue
	call parabreg(psi,psip,psipp,rf,vf)
	v(ii)=vraw(ii)+
     1        (1.d0-phi0(ii)/phi(ii))*(2.d0*psip/psi*fp/f+fpp/f)/2.d0
 300    continue
 1170   call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
	call integ(ev,l,xkappa,n,nn,jrt,ief,x0,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
	do 8015 ii=1,jrt
	phi(ii)=phi(ii)/phi(jrt)
 8015   continue
	xn0=0.d0
	do 8020 ii=1,jrt-1
	xn0=xn0+dr( ii)*phi( ii)*phi( ii)
 8020   continue
	xn0=xn0+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
	de=0.0001d0
	ee=ev+de/2.d0
	call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
	ee=ev-de/2.d0
	call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
	c0=(xm-xp)/(2.d0*de)
	write (6,94)  c0,x0
	write (6,94) xn0
	if (dabs(c0-c00).ge.0.000000001d0) then
	  dqddel=2.*(xi1+deltal*xi2)
	  deltal=deltal+(c00-c0)/dqddel
	  goto 225
	endif
	write (6,94)  c0,x0
	write (6,*) 'NCPP ACHIEVED !!!'
	return
	end
