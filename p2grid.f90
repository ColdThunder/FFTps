module p2grid
implicit none

contains

function index2lagpos(id,ng)
implicit none
integer(8)::id
integer(4)::ng
real(4)::index2lagpos(3)
! index range from 0 - ng^3-1
! output range is (0,ng]
integer(8)::id2
integer(4)::i,j,k
id2=id
k=id2/ng**2
id2=id2-k*ng**2
j=id2/ng
i=id2-j*ng
index2lagpos=[k,j,i]
endfunction index2lagpos

! ===== density =====
subroutine massassign(np,pos,ng,den8,cmd)
use omp_lib
implicit none
integer(4)::np
real(4)::pos(3,np)
integer(4)::ng
real(8)::den8(ng,ng,ng)
integer(4)::cmd
integer(4)::ibin,jbin,kbin,ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
real(4)::hx,hy,hz,hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1
integer(4)::pid
real(8)::tmp1,tmp2
integer(4)::i,j,k
!$omp parallel do default(shared)
do k=1,ng
  den8(:,:,k)=0.d0
enddo
!$omp end parallel do
write(*,*) 'begin mass assignment by method',cmd
if (cmd.eq.1) then
  ! NGP
  !$omp parallel do default(private) shared(np,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    !$omp atomic
    den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+1.d0
  enddo
  !$omp end parallel do
elseif (cmd.eq.2) then
  ! CIC
  !$omp parallel do default(private) shared(np,ng,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=1.-abs(hx)
    if (hx.gt.0.) then
      hxp1=hx
      hxm1=0.
    else
      hxp1=0.
      hxm1=-hx
    endif
    hy0=1.-abs(hy)
    if (hy.gt.0.) then
      hyp1=hy
      hym1=0.
    else
      hyp1=0.
      hym1=-hy
    endif
    hz0=1.-abs(hz)
    if (hz.gt.0.) then
      hzp1=hz
      hzm1=0.
    else
      hzp1=0.
      hzm1=-hz
    endif
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
elseif (cmd.eq.3) then
  ! TSC
  !$omp parallel do default(private) shared(np,ng,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=0.75-hx**2
    hxp1=0.5*(0.5+hx)**2
    hxm1=0.5*(0.5-hx)**2
    hy0=0.75-hy**2
    hyp1=0.5*(0.5+hy)**2
    hym1=0.5*(0.5-hy)**2
    hz0=0.75-hz**2
    hzp1=0.5*(0.5+hz)**2
    hzm1=0.5*(0.5-hz)**2
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
else
  write(*,*) 'ERROR: mass assignment not specified'
endif
write(*,*) 'statistics:'
tmp1=0.d0
tmp2=0.d0
!$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
do k=1,ng
  den8(:,:,k)=den8(:,:,k)/np*real(ng)**3-1.d0
  tmp1=tmp1+sum(den8(:,:,k))/real(ng)**3
  tmp2=tmp2+sum(den8(:,:,k)**2)/real(ng)**3
enddo
!$omp end parallel do
tmp2=sqrt(tmp2)
write(*,*) 'mean:',real(tmp1)
write(*,*) 'sigma:',real(tmp2)
endsubroutine massassign

! ===== density =====
subroutine massassigndown(np,pos,ng,den8,cmd,nc,coarse)
use omp_lib
implicit none
integer(8)::np
real(4)::pos(3,np)
integer(4)::ng
real(8)::den8(ng,ng,ng)
integer(4)::cmd
integer(4)::nc,coarse
integer(4)::ibin,jbin,kbin,ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
real(4)::hx,hy,hz,hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1
integer(8)::pid
integer(8)::i8,j8,k8
real(8)::tmp1,tmp2
integer(4)::i,j,k
!$omp parallel do default(shared)
do k=1,ng
  den8(:,:,k)=0.d0
enddo
!$omp end parallel do
write(*,*) 'begin mass assignment by method',cmd
if (cmd.eq.1) then
  ! NGP
  !$omp parallel do default(private) shared(np,pos,den8,nc,coarse) schedule(static)
  do pid=1,np
    k8=(pid-1)/nc/nc
    i8=(pid-1)-k8*nc*nc
    j8=i8/nc
    i8=i8-j8*nc
    if (any( mod([i8,j8,k8],coarse).ne.0 )) cycle
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    !$omp atomic
    den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+1.d0
  enddo
  !$omp end parallel do
elseif (cmd.eq.2) then
  ! CIC
  !$omp parallel do default(private) shared(np,ng,pos,den8,nc,coarse) schedule(static)
  do pid=1,np
    k8=(pid-1)/nc/nc
    i8=(pid-1)-k8*nc*nc
    j8=i8/nc
    i8=i8-j8*nc
    if (any( mod([i8,j8,k8],coarse).ne.0 )) cycle
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=1.-abs(hx)
    if (hx.gt.0.) then
      hxp1=hx
      hxm1=0.
    else
      hxp1=0.
      hxm1=-hx
    endif
    hy0=1.-abs(hy)
    if (hy.gt.0.) then
      hyp1=hy
      hym1=0.
    else
      hyp1=0.
      hym1=-hy
    endif
    hz0=1.-abs(hz)
    if (hz.gt.0.) then
      hzp1=hz
      hzm1=0.
    else
      hzp1=0.
      hzm1=-hz
    endif
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
elseif (cmd.eq.3) then
  ! TSC
  !$omp parallel do default(private) shared(np,ng,pos,den8,nc,coarse) schedule(static)
  do pid=1,np
    k8=(pid-1)/nc/nc
    i8=(pid-1)-k8*nc*nc
    j8=i8/nc
    i8=i8-j8*nc
    if (any( mod([i8,j8,k8],coarse).ne.0 )) cycle
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=0.75-hx**2
    hxp1=0.5*(0.5+hx)**2
    hxm1=0.5*(0.5-hx)**2
    hy0=0.75-hy**2
    hyp1=0.5*(0.5+hy)**2
    hym1=0.5*(0.5-hy)**2
    hz0=0.75-hz**2
    hzp1=0.5*(0.5+hz)**2
    hzm1=0.5*(0.5-hz)**2
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
else
  write(*,*) 'ERROR: mass assignment not specified'
endif
write(*,*) 'statistics:'
tmp1=0.d0
tmp2=0.d0
!$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
do k=1,ng
  den8(:,:,k)=den8(:,:,k)/np*coarse**3*real(ng)**3-1.d0
  tmp1=tmp1+sum(den8(:,:,k))/real(ng)**3
  tmp2=tmp2+sum(den8(:,:,k)**2)/real(ng)**3
enddo
!$omp end parallel do
tmp2=sqrt(tmp2)
write(*,*) 'mean:',real(tmp1)
write(*,*) 'sigma:',real(tmp2)
endsubroutine massassigndown


subroutine massassign_double(np,pos,ng,den8,cmd)
use omp_lib
implicit none
integer(4)::np
real(8)::pos(3,np)
integer(4)::ng
real(8)::den8(ng,ng,ng)
integer(4)::cmd
integer(4)::ibin,jbin,kbin,ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
real(8)::hx,hy,hz,hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1
integer(4)::pid
real(8)::tmp1,tmp2
integer(4)::i,j,k
!$omp parallel do default(shared)
do k=1,ng
  den8(:,:,k)=0.d0
enddo
!$omp end parallel do
write(*,*) 'begin mass assignment by method',cmd
if (cmd.eq.1) then
  ! NGP
  !$omp parallel do default(private) shared(np,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    !$omp atomic
    den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+1.d0
  enddo
  !$omp end parallel do
elseif (cmd.eq.2) then
  ! CIC
  !$omp parallel do default(private) shared(np,ng,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=1.-abs(hx)
    if (hx.gt.0.) then
      hxp1=hx
      hxm1=0.
    else
      hxp1=0.
      hxm1=-hx
    endif
    hy0=1.-abs(hy)
    if (hy.gt.0.) then
      hyp1=hy
      hym1=0.
    else
      hyp1=0.
      hym1=-hy
    endif
    hz0=1.-abs(hz)
    if (hz.gt.0.) then
      hzp1=hz
      hzm1=0.
    else
      hzp1=0.
      hzm1=-hz
    endif
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
elseif (cmd.eq.3) then
  ! TSC
  !$omp parallel do default(private) shared(np,ng,pos,den8) schedule(static)
  do pid=1,np
    ibin=ceiling(pos(1,pid))
    jbin=ceiling(pos(2,pid))
    kbin=ceiling(pos(3,pid))
    hx=pos(1,pid)-ibin+0.5
    hy=pos(2,pid)-jbin+0.5
    hz=pos(3,pid)-kbin+0.5
    hx0=0.75-hx**2
    hxp1=0.5*(0.5+hx)**2
    hxm1=0.5*(0.5-hx)**2
    hy0=0.75-hy**2
    hyp1=0.5*(0.5+hy)**2
    hym1=0.5*(0.5-hy)**2
    hz0=0.75-hz**2
    hzp1=0.5*(0.5+hz)**2
    hzm1=0.5*(0.5-hz)**2
    ibinm1=mod(ibin-2+ng,ng)+1
    jbinm1=mod(jbin-2+ng,ng)+1
    kbinm1=mod(kbin-2+ng,ng)+1
    ibinp1=mod(ibin+ng,ng)+1
    jbinp1=mod(jbin+ng,ng)+1
    kbinp1=mod(kbin+ng,ng)+1
    !$omp atomic
    den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1
    !$omp atomic
    den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0
    !$omp atomic
    den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0
    !$omp atomic
    den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0
    !$omp atomic
    den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0
    !$omp atomic
    den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0
    !$omp atomic
    den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0
    !$omp atomic
    den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0
    !$omp atomic
    den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0
    !$omp atomic
    den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0
    !$omp atomic
    den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1
    !$omp atomic
    den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1
    !$omp atomic
    den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1
    !$omp atomic
    den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1
    !$omp atomic
    den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1
    !$omp atomic
    den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1
    !$omp atomic
    den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1
    !$omp atomic
    den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1
    !$omp atomic
    den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1
  enddo
  !$omp end parallel do
else
  write(*,*) 'ERROR: mass assignment not specified'
endif
write(*,*) 'statistics:'
tmp1=0.d0
tmp2=0.d0
!$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
do k=1,ng
  den8(:,:,k)=den8(:,:,k)/np*real(ng)**3-1.d0
  tmp1=tmp1+sum(den8(:,:,k))/real(ng)**3
  tmp2=tmp2+sum(den8(:,:,k)**2)/real(ng)**3
enddo
!$omp end parallel do
tmp2=sqrt(tmp2)
write(*,*) 'mean:',real(tmp1)
write(*,*) 'sigma:',real(tmp2)
endsubroutine massassign_double




endmodule p2grid
