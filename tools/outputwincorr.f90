program outputwincorr
use myomp
use head
implicit none

integer(4)::LL(3),LLk(3)
real(4),allocatable::vc3(:,:,:)

real(4)::boxinfo(2,3)
real(4)::fsvec(3)

integer(4)::i,j,k,dir

real(8),allocatable::pk8(:,:)
real(8),allocatable::pk2d(:,:,:)
real(4)::fs
real(4)::npeff
integer(4)::pp
logical(4)::debug=.false.

real(4),parameter::pi=4.*atan(1.)
real(4)::kbasic

LL=[L,L,L]
LLk=[L+2,L,L]
boxinfo=reshape([0,LL(1),0,LL(2),0,LL(3)],[2,3])

call binscheme 
call memo(0)


if (sncorr) then
  pp=0
  if (second) then
    if (wincorr) pp=pp+painter
    if (wincorr2) pp=pp+painter2
  else
    pp=pp+painter*2
  endif
  if (interlace) then
    call getshotnoise_interlace_new(vc3,LL,1.,pp)
  else
    call getshotnoise(vc3,LL,1.,pp)
  endif
  call calcpkell(vc3,LL,pk8,kbin,outsnname)
else
  write(*,*) 'shot noise is not corrected'
endif

if (wincorr.or.(wincorr2.and.second)) then
  pp=0
  if (second) then
    if (wincorr) pp=pp+painter
    if (wincorr2) pp=pp+painter2
  else
    if (wincorr) pp=pp+painter*2
  endif

  call getnormalization(vc3,LL,pp)
  call calcpkell(vc3,LL,pk8,kbin,outwinname)
endif

call memo(1)

contains

subroutine binscheme
  implicit none
  kbasic=2*pi/box
  write(*,'(a,f8.3,a,f8.3)') 'kbasic=',kbasic,'  knyquist=',kbasic*L/2
  if (kbin.eq.0.and.dk.eq.0.) then
    write(*,*) 'ERROR: kbin and dk cannot to 0 in the same time'
    stop
  endif

  if (usedefaulkminkmax) then
    if (logbin) then
      kmax=alog10(kbasic*L/2)
      kmin=alog10(kbasic*1.28)
    else
      kmax=kbasic*(L+1)/2
      kmin=kbasic/2
    endif
  else
    if (logbin) then
      kmax=alog10(kmax)
      if (kmin.eq.0) then
        write(*,*) 'ERROR: kmin=0, set kmin=1.28*kbasic'
        kmin=alog10(kbasic*1.28)
      else
        kmin=alog10(kmin)
      endif
    endif
  endif
  
  if (kbin.ne.0) then
    dk=(kmax-kmin)/kbin
  else
    kbin=floor((kmax-kmin)/dk)
    kmax=kmin+dk*kbin
  endif
  if (logbin) then
    write(*,'(a,f8.3,a,f8.3)') 'kmin=',10.**kmin,'  kmax=',10.**kmax
    write(*,'(a,f8.3,a,i8)') 'dlogk=',dk,'  kbin=',kbin
  else
    write(*,'(a,f8.3,a,f8.3)') 'kmin=',kmin,'  kmax=',kmax
    write(*,'(a,f8.3,a,i8)') 'dk=',dk,'  kbin=',kbin
  endif
endsubroutine binscheme

subroutine memo(command)
  implicit none
  integer(4)::command
  if (command.eq.0) then
    write(*,*) 'memo for vc3:',(float(LL(1)+2)*LL(2)*LL(3))*4/1024.**3,'G'
    allocate(vc3(L+2,L,L))
    allocate(pk8(9,kbin))
  elseif(command.eq.1) then
    deallocate(vc3)
    deallocate(pk8)
  endif
endsubroutine memo

subroutine getnormalization(vc,LL,pp)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  real(4)::hpi(3),a1,a2,a3,d2,ri,rj,rk
  real(4)::wx(LL(1)/2+1),wy(LL(2)),wz(LL(3))
  integer(4)::pp
  if (pp.eq.0) then
    write(*,*) 'no window function correction'
    return
  endif
  write(*,*) 'window function is deconvolved, pp=',pp
  hpi=pi/float(LL)

  wz(0)=1
  !$omp parallel do default(private) shared(hpi,LL,wz) schedule(static)
  do k=2,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    a3=hpi(3)*rk
    wz(k)=sin(a3)/a3
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(hpi,LL,wy) schedule(static)
  do j=1,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    a2=hpi(2)*rj
    wy(j)=sin(a2)/a2
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(hpi,LL,wx) schedule(static)
  do i=1,LL(1)+2,2
    ri=i/2
    a1=hpi(1)*ri
    wx(i/2+1)=sin(a1)/a1
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(vc,wx,wy,wz,pp) schedule(static)
  do k=1,LL(3)
    do j=1,LL(2)
      do i=1,LL(1)+2,2
        d2=(wx(i/2+1)*wy(j)*wz(k))**pp
        vc(i:i+1,j,k)=d2
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine getnormalization

subroutine getshotnoise(vc,LL,npeff,pp)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  real(4)::npeff
  integer(4)::pp
  real(4)::hpi(3),a1,a2,a3,d2,ri,rj,rk
  real(4)::wx(LL(1)/2+1),wy(LL(2)),wz(LL(3))
  if (npeff.eq.0) then
    write(*,*) 'can not correct shot noise with unknown particle number'
    return
  endif
  write(*,*) 'shot noise is corrected'
  hpi=pi/float(LL)

  !$omp parallel do default(private) shared(hpi,LL,pp,wz) schedule(static)
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    a3=sin(hpi(3)*rk)
    wz(k)=wcorr(a3,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(hpi,LL,pp,wy) schedule(static)
  do j=1,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    a2=sin(hpi(2)*rj)
    wy(j)=wcorr(a2,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(hpi,LL,pp,wx) schedule(static)
  do i=1,LL(1)+2,2
    ri=i/2
    a1=sin(hpi(1)*ri)
    wx(i/2+1)=wcorr(a1,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(vc,wx,wy,wz,LL,npeff) schedule(static)
  do k=1,LL(3)
    do j=1,LL(2)
      do i=1,LL(1)+2,2
        d2 = wx(i/2+1)*wy(j)*wz(k)
        vc(i:i+1,j,k)=d2/npeff
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine getshotnoise

function wcorr(aa,pp)
  implicit none
  real(4)::wcorr
  real(4)::aa  ! sinx
  integer(4)::pp
  real(4)::bb  ! cosx
  if (mod(pp,2).ne.0) bb=sqrt(1.-aa**2)
  if (pp.eq.2) then
    wcorr=1
  elseif (pp.eq.3) then
    wcorr=bb
  elseif (pp.eq.4) then
    wcorr=1-2./3*aa**2
  elseif (pp.eq.5) then
    wcorr=bb*(1-1./3*aa**2)
  elseif (pp.eq.6) then
    wcorr=1-aa**2+2./15*aa**4
  elseif (pp.eq.7) then
    wcorr=bb*(1-2./3*aa**2+2./45*aa**4)
  elseif(pp.eq.8) then
    wcorr=1-4./3*aa**2+2./5*aa**4-4./315*aa**6
  else
    wcorr=1
  endif
endfunction wcorr

subroutine getshotnoise_interlace(vc,LL,npeff,pp)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  real(4)::npeff
  integer(4)::pp
  real(4)::hpi(3),a1,a2,a3,d2,ri,rj,rk
  real(4)::wx(LL(1)/2+1),wy(LL(2)),wz(LL(3))
  if (npeff.eq.0) then
    write(*,*) 'can not correct shot noise with unknown particle number'
    return
  endif
  write(*,*) 'shot noise is corrected (interlaced)'
  write(*,*) 'pp=',pp
  write(*,*) 'npeff=',npeff
  hpi=pi/float(LL)/2
  !$omp parallel do default(private) shared(vc,hpi,LL,npeff,painter,pp) schedule(static)
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    a3=sin(hpi(3)*rk)
    a3=wcorr_interlace(a3,pp)

    do j=1,LL(2)
      if (j.le.LL(2)/2+1) then
        rj=j-1
      else
        rj=j-1-LL(2)
      endif
      a2=sin(hpi(2)*rj)
      a2=wcorr_interlace(a2,pp)

      do i=1,LL(1)+2,2
        ri=i/2
        a1=sin(hpi(1)*ri)
        a1=wcorr_interlace(a1,pp)

        vc(i:i+1,j,k)=a1*a2*a3/npeff
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine getshotnoise_interlace

function wcorr_interlace(aa,pp)
  implicit none
  real(4)::wcorr_interlace
  real(4)::aa  ! sinx
  integer(4)::pp
  real(4)::bb  ! cosx
  bb=sqrt(1.-aa**2)
  if (pp.eq.2) then
    wcorr_interlace=bb**2
  elseif (pp.eq.3) then
    wcorr_interlace=bb**4
  elseif (pp.eq.4) then
    wcorr_interlace=(1-2./3*aa**2)*bb**4
  elseif (pp.eq.5) then
    wcorr_interlace=(1-1./3*aa**2)*bb**6
  elseif (pp.eq.6) then
    wcorr_interlace=(1-aa**2+2./15*aa**4)*bb**6
  elseif (pp.eq.7) then
    wcorr_interlace=(1-2./3*aa**2+2./45*aa**4)*bb**8
  elseif(pp.eq.8) then
    wcorr_interlace=(1-4./3*aa**2+2./5*aa**4-4./315*aa**6)*bb**8
  else
    wcorr_interlace=1
  endif
endfunction wcorr_interlace

subroutine getshotnoise_interlace_new(vc,LL,npeff,pp)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  real(4)::npeff
  integer(4)::pp
  real(4)::hpi(3),a1,a2,a3,d2,ri,rj,rk
  real(4)::wx(LL(1)/2+1,2),wy(LL(2),2),wz(LL(3),2)
  if (npeff.eq.0) then
    write(*,*) 'can not correct shot noise with unknown particle number'
    return
  endif
  write(*,*) 'shot noise is corrected (interlaced)'
  write(*,*) 'pp=',pp
  write(*,*) 'npeff=',npeff
  hpi=pi/float(LL)

  !$omp parallel do default(private) shared(LL,hpi,wz,pp) schedule(static)
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    a3=hpi(3)*rk/2
    wz(k,1:2)=wcorr_interlace_new(a3,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(LL,hpi,wy,pp) schedule(static)
  do j=1,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    a2=hpi(2)*rj/2
    wy(j,1:2)=wcorr_interlace_new(a2,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(LL,hpi,wx,pp) schedule(static)
  do i=1,LL(1)+2,2
    ri=i/2
    a1=hpi(1)*ri/2
    wx(i/2+1,1:2)=wcorr_interlace_new(a1,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(vc,LL,wx,wy,wz,npeff) schedule(static)
  do k=1,LL(3)
    do j=1,LL(2)
      do i=1,LL(1)+2,2
        d2= wx(i/2+1,1)*wy(j,1)*wz(k,1) + wx(i/2+1,1)*wy(j,2)*wz(k,2) + wx(i/2+1,2)*wy(j,1)*wz(k,2) + wx(i/2+1,2)*wy(j,2)*wz(k,1)
        vc(i:i+1,j,k)=d2/npeff
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine getshotnoise_interlace_new

function wcorr_interlace_new(aa,pp)
  implicit none
  real(4)::wcorr_interlace_new(2)
  real(4)::aa
  integer(4)::pp
  real(4)::sinx,cosx
  sinx=sin(aa)
  cosx=cos(aa)

  if (pp.eq.2) then
    wcorr_interlace_new(1)=cosx**2
    wcorr_interlace_new(2)=sinx**2
  !elseif (pp.eq.3) then
  !  wcorr_interlace_new(1)=coshx**4
  elseif (pp.eq.4) then
    wcorr_interlace_new(1)=(1-2./3*sinx**2)*cosx**4
    wcorr_interlace_new(2)=(1-2./3*cosx**2)*sinx**4
  !elseif (pp.eq.5) then
  !  wcorr_interlace_new(1)=(1-1./3*sinhx**2)*coshx**6
  elseif (pp.eq.6) then
    wcorr_interlace_new(1)=(1-sinx**2+2./15*sinx**4)*cosx**6
    wcorr_interlace_new(2)=(1-cosx**2+2./15*cosx**4)*sinx**6
  !elseif (pp.eq.7) then
  !  wcorr_interlace_new(1)=(1-2./3*sinhx**2+2./45*sinhx**4)*coshx**8
  elseif(pp.eq.8) then
    wcorr_interlace_new(1)=(1-4./3*sinx**2+2./5*sinx**4-4./315*sinx**6)*cosx**8
    wcorr_interlace_new(2)=(1-4./3*cosx**2+2./5*cosx**4-4./315*cosx**6)*sinx**8
  else
    wcorr_interlace_new(1:2)=1
  endif
endfunction wcorr_interlace_new



subroutine calcpkell(vc,LL,pk8,kbin,filename)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  integer(4)::kbin
  character(*)::filename
  real(8)::pk8(9,kbin)
  real(4)::ri,rj,rk,kmode,mu
  real(4)::pl0,pl2,pl4
  integer(4)::bin
  real(4)::dr
  real(4)::factor
    
  pk8=0.
  write(*,*) 'calculating power spectrum'
    
  !$omp parallel do default(private) shared(vc,LL,kbasic,kmin,dk,kbin,logbin) &
  !$omp reduction(+:pk8) schedule(guided)
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    rk=kbasic*rk
    do j=1,LL(2)
      if (j.le.LL(2)/2+1) then
        rj=j-1
      else
        rj=j-1-LL(2)
      endif
      rj=kbasic*rj
      do i=1,LL(1)+2,2
        ri=i/2
        ri=kbasic*ri
    
        if (i.eq.1.and.j.eq.1.and.k.eq.1) cycle
        if (i.eq.1) then
          factor=1
        else
          factor=2
        endif
        kmode=sqrt(ri*ri+rj*rj+rk*rk)
        mu=rk/kmode
        if (logbin) kmode=alog10(kmode)
        bin=floor((kmode-kmin)/dk)+1
        if (bin.lt.1.or.bin.gt.kbin) cycle
        pl0=1.*vc(i,j,k)
        pl2=0.5*(3*mu**2-1)*vc(i,j,k)
        pl4=0.125*(35*mu**4-30*mu**2+3)*vc(i,j,k)
        pk8(1,bin)=pk8(1,bin)+kmode*factor
        pk8(2,bin)=pk8(2,bin)+kmode**2*factor
        pk8(3,bin)=pk8(3,bin)+pl0*factor
        pk8(4,bin)=pk8(4,bin)+pl0**2*factor
        pk8(5,bin)=pk8(5,bin)+pl2*factor
        pk8(6,bin)=pk8(6,bin)+pl2**2*factor
        pk8(7,bin)=pk8(7,bin)+pl4*factor
        pk8(8,bin)=pk8(8,bin)+pl4**2*factor
        pk8(9,bin)=pk8(9,bin)+1.*factor
      enddo
    enddo
  enddo
  !$omp end parallel do
    
  pk8(1,:)=pk8(1,:)/pk8(9,:)
  pk8(2,:)=pk8(2,:)/pk8(9,:)-pk8(1,:)**2
  pk8(3,:)=pk8(3,:)/pk8(9,:)
  pk8(4,:)=pk8(4,:)/pk8(9,:)-pk8(3,:)**2
  pk8(5,:)=pk8(5,:)/pk8(9,:)
  pk8(6,:)=pk8(6,:)/pk8(9,:)-pk8(5,:)**2
  pk8(7,:)=pk8(7,:)/pk8(9,:)
  pk8(8,:)=pk8(8,:)/pk8(9,:)-pk8(7,:)**2

    
  if (usecenter) then
    do i=1,kbin
      pk8(1,i)=kmin+(i-0.5)*dk
    enddo
    if (logbin) pk8(1,:)=10.**pk8(1,:)
    pk8(2,:)=0.
  else
    if (logbin) then
      pk8(1,:)=10.**pk8(1,:)
      pk8(2,:)=(10.**sqrt(pk8(2,:))-1.)*pk8(1,:)
    endif
  endif
    
  write(*,*) 'writing: ',trim(filename)
  open(32,file=filename,form='formatted',status='replace')
  do i=1,kbin
    write(32,'(9e14.6)') pk8(1:9:2,i)
    if (debug) write(*,*) real(pk8(1:9:2,i))
  enddo
  close(32)
endsubroutine calcpkell



endprogram outputwincorr
