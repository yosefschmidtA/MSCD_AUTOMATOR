	subroutine parabreg(f,fp,fpp,rf,vf)
	implicit real*8 (a-h,o-z)
	dimension rf(3),vf(3)
	f=vf(2)
	r21=rf(2)-rf(1)
	r32=rf(3)-rf(2)
	v21=vf(2)-vf(1)
	v32=vf(3)-vf(2)
	fp=(v21+v32)/(r21+r32)
	fpp=(v32/r32-v21/r21)/((r21+r32)/2.d0)
	return
	end
