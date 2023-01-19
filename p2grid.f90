module p2grid
  implicit none
  
  interface massassign
    module procedure massassign_r4,massassign_r8,massassign_long_r4,massassign_long_r8
    module procedure massassignw_r4,massassignw_r8,massassignw_long_r4,massassignw_long_r8
  endinterface massassign
  
  
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
  subroutine massassign_r4(np,pos,LL,den8,cmd)
    use omp_lib
    implicit none
    integer(4)::np
    real(4)::pos(3,np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
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
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr
      enddo
      !$omp end parallel do
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
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
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
  endsubroutine massassign_r4
  
  ! ===== density =====
    ! ===== density =====
  subroutine massassign_r8(np,pos,LL,den8,cmd)
    use omp_lib
    implicit none
    integer(4)::np
    real(8)::pos(3,np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
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
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr
      enddo
      !$omp end parallel do
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
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
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
  endsubroutine massassign_r8
  
  ! ===== density =====
  subroutine massassign_long_r4(np,pos,LL,den8,cmd)
    use omp_lib
    implicit none
    integer(8)::np
    real(4)::pos(3,np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
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
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr
      enddo
      !$omp end parallel do
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
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
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
  endsubroutine massassign_long_r4
  
  ! ===== density =====
    ! ===== density =====
  subroutine massassign_long_r8(np,pos,LL,den8,cmd)
    use omp_lib
    implicit none
    integer(8)::np
    real(8)::pos(3,np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
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
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr
      enddo
      !$omp end parallel do
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,den8) schedule(static)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
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
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
  endsubroutine massassign_long_r8
  

  !!!!!! with weights
    ! ===== density =====
  subroutine massassignw_r4(np,pos,mass,LL,den8,cmd,npeff)
    use omp_lib
    implicit none
    integer(4)::np
    real(4)::pos(3,np)
    real(4)::mass(np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    real(kind(mass))::npeff
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    real(8)::w1,w2
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
      den8(:,:,k)=0.d0
    enddo
    !$omp end parallel do
    write(*,*) 'begin mass assignment by method',cmd
    w1=0.d0
    w2=0.d0
    if (cmd.eq.1) then
      ! NGP
      !$omp parallel do default(private) shared(np,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        !$omp atomic
        den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    elseif (cmd.eq.2) then
      ! CIC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl*mass(pid)
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl*mass(pid)
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl*mass(pid)
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl*mass(pid)
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr*mass(pid)
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr*mass(pid)
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr*mass(pid)
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
        !$omp atomic
        den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    else
      write(*,*) 'ERROR: mass assignment not specified'
    endif
    write(*,*) 'statistics:'
    tmp1=0.d0
    tmp2=0.d0
    !$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
    npeff=w1**2/w2
  endsubroutine massassignw_r4
  
  ! ===== density =====
  subroutine massassignw_r8(np,pos,mass,LL,den8,cmd,npeff)
    use omp_lib
    implicit none
    integer(4)::np
    real(8)::pos(3,np)
    real(8)::mass(np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    real(kind(mass))::npeff
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    real(8)::w1,w2
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
      den8(:,:,k)=0.d0
    enddo
    !$omp end parallel do
    write(*,*) 'begin mass assignment by method',cmd
    w1=0.d0
    w2=0.d0
    if (cmd.eq.1) then
      ! NGP
      !$omp parallel do default(private) shared(np,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        !$omp atomic
        den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    elseif (cmd.eq.2) then
      ! CIC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl*mass(pid)
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl*mass(pid)
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl*mass(pid)
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl*mass(pid)
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr*mass(pid)
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr*mass(pid)
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr*mass(pid)
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
        !$omp atomic
        den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    else
      write(*,*) 'ERROR: mass assignment not specified'
    endif
    write(*,*) 'statistics:'
    tmp1=0.d0
    tmp2=0.d0
    !$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
    npeff=w1**2/w2
  endsubroutine massassignw_r8
  
  ! ===== density =====
  subroutine massassignw_long_r4(np,pos,mass,LL,den8,cmd,npeff)
    use omp_lib
    implicit none
    integer(8)::np
    real(4)::pos(3,np)
    real(4)::mass(np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    real(kind(mass))::npeff
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    real(8)::w1,w2
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
      den8(:,:,k)=0.d0
    enddo
    !$omp end parallel do
    write(*,*) 'begin mass assignment by method',cmd
    w1=0.d0
    w2=0.d0
    if (cmd.eq.1) then
      ! NGP
      !$omp parallel do default(private) shared(np,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        !$omp atomic
        den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    elseif (cmd.eq.2) then
      ! CIC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl*mass(pid)
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl*mass(pid)
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl*mass(pid)
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl*mass(pid)
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr*mass(pid)
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr*mass(pid)
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr*mass(pid)
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
        !$omp atomic
        den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    else
      write(*,*) 'ERROR: mass assignment not specified'
    endif
    write(*,*) 'statistics:'
    tmp1=0.d0
    tmp2=0.d0
    !$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
    npeff=w1**2/w2
  endsubroutine massassignw_long_r4
  
  ! ===== density =====
  subroutine massassignw_long_r8(np,pos,mass,LL,den8,cmd,npeff)
    use omp_lib
    implicit none
    integer(8)::np
    real(8)::pos(3,np)
    real(8)::mass(np)
    integer(4)::LL(3)
    real(8)::den8(LL(1),LL(2),LL(3))
    integer(4)::cmd
    real(kind(mass))::npeff
    ! variables for ngp
    integer(4)::ibin,jbin,kbin
    ! variables for cic
    integer(4)::il,ir,jl,jr,kl,kr
    real(kind(pos))::hx,hy,hz,hxl,hxr,hyl,hyr,hzl,hzr
    ! variables for tsc
    integer(4)::ibinm1,jbinm1,kbinm1,ibinp1,jbinp1,kbinp1
    real(kind(pos))::hx0,hy0,hz0,hxm1,hym1,hzm1,hxp1,hyp1,hzp1

    integer(kind(np))::pid
    real(8)::tmp1,tmp2,vol
    integer(4)::i,j,k
    real(8)::w1,w2
    vol=product(real(LL))
    !$omp parallel do default(shared)
    do k=1,LL(3)
      den8(:,:,k)=0.d0
    enddo
    !$omp end parallel do
    write(*,*) 'begin mass assignment by method',cmd
    w1=0.d0
    w2=0.d0
    if (cmd.eq.1) then
      ! NGP
      !$omp parallel do default(private) shared(np,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        !$omp atomic
        den8(ibin,jbin,kbin)=den8(ibin,jbin,kbin)+mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    elseif (cmd.eq.2) then
      ! CIC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
      do pid=1,np
        ibin=ceiling(pos(1,pid))
        jbin=ceiling(pos(2,pid))
        kbin=ceiling(pos(3,pid))
        hx=pos(1,pid)-ibin+0.5
        hy=pos(2,pid)-jbin+0.5
        hz=pos(3,pid)-kbin+0.5
        if (hx.gt.0.) then
          il=ibin
          ir=mod(ibin+LL(1),LL(1))+1
          hxl=1-hx
          hxr=hx
        else
          il=mod(ibin-2+LL(1),LL(1))+1
          ir=ibin
          hxl=-hx
          hxr=1+hx
        endif
        if (hy.gt.0.) then
          jl=jbin
          jr=mod(jbin+LL(2),LL(2))+1
          hyl=1-hy
          hyr=hy
        else
          jl=mod(jbin-2+LL(2),LL(2))+1
          jr=jbin
          hyl=-hy
          hyr=1+hy
        endif
        if (hz.gt.0.) then
          kl=kbin
          kr=mod(kbin+LL(3),LL(3))+1
          hzl=1-hz
          hzr=hz
        else
          kl=mod(kbin-2+LL(3),LL(3))+1
          kr=kbin
          hzl=-hz
          hzr=1+hz
        endif
        !$omp atomic
        den8(il,jl,kl)=den8(il,jl,kl)+hxl*hyl*hzl*mass(pid)
        !$omp atomic
        den8(ir,jl,kl)=den8(ir,jl,kl)+hxr*hyl*hzl*mass(pid)
        !$omp atomic
        den8(il,jr,kl)=den8(il,jr,kl)+hxl*hyr*hzl*mass(pid)
        !$omp atomic
        den8(ir,jr,kl)=den8(ir,jr,kl)+hxr*hyr*hzl*mass(pid)
        !$omp atomic
        den8(il,jl,kr)=den8(il,jl,kr)+hxl*hyl*hzr*mass(pid)
        !$omp atomic
        den8(ir,jl,kr)=den8(ir,jl,kr)+hxr*hyl*hzr*mass(pid)
        !$omp atomic
        den8(il,jr,kr)=den8(il,jr,kr)+hxl*hyr*hzr*mass(pid)
        !$omp atomic
        den8(ir,jr,kr)=den8(ir,jr,kr)+hxr*hyr*hzr*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
    elseif (cmd.eq.3) then
      ! TSC
      !$omp parallel do default(private) shared(np,LL,pos,mass,den8) schedule(static) reduction(+:w1,w2)
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
        ibinm1=mod(ibin-2+LL(1),LL(1))+1
        jbinm1=mod(jbin-2+LL(2),LL(2))+1
        kbinm1=mod(kbin-2+LL(3),LL(3))+1
        ibinp1=mod(ibin+LL(1),LL(1))+1
        jbinp1=mod(jbin+LL(2),LL(2))+1
        kbinp1=mod(kbin+LL(3),LL(3))+1
        !$omp atomic
        den8(ibinm1,jbinm1,kbinm1)=den8(ibinm1,jbinm1,kbinm1)+hxm1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinm1)=den8(ibin  ,jbinm1,kbinm1)+hx0 *hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinm1)=den8(ibinp1,jbinm1,kbinm1)+hxp1*hym1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinm1)=den8(ibinm1,jbin  ,kbinm1)+hxm1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinm1)=den8(ibin  ,jbin  ,kbinm1)+hx0 *hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinm1)=den8(ibinp1,jbin  ,kbinm1)+hxp1*hy0 *hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinm1)=den8(ibinm1,jbinp1,kbinm1)+hxm1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinm1)=den8(ibin  ,jbinp1,kbinm1)+hx0 *hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinm1)=den8(ibinp1,jbinp1,kbinm1)+hxp1*hyp1*hzm1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbin  )=den8(ibinm1,jbinm1,kbin  )+hxm1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbin  )=den8(ibin  ,jbinm1,kbin  )+hx0 *hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbin  )=den8(ibinp1,jbinm1,kbin  )+hxp1*hym1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbin  )=den8(ibinm1,jbin  ,kbin  )+hxm1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbin  )=den8(ibin  ,jbin  ,kbin  )+hx0 *hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbin  )=den8(ibinp1,jbin  ,kbin  )+hxp1*hy0 *hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbin  )=den8(ibinm1,jbinp1,kbin  )+hxm1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbin  )=den8(ibin  ,jbinp1,kbin  )+hx0 *hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbin  )=den8(ibinp1,jbinp1,kbin  )+hxp1*hyp1*hz0 *mass(pid)
        !$omp atomic
        den8(ibinm1,jbinm1,kbinp1)=den8(ibinm1,jbinm1,kbinp1)+hxm1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinm1,kbinp1)=den8(ibin  ,jbinm1,kbinp1)+hx0 *hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinm1,kbinp1)=den8(ibinp1,jbinm1,kbinp1)+hxp1*hym1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbin  ,kbinp1)=den8(ibinm1,jbin  ,kbinp1)+hxm1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbin  ,kbinp1)=den8(ibin  ,jbin  ,kbinp1)+hx0 *hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbin  ,kbinp1)=den8(ibinp1,jbin  ,kbinp1)+hxp1*hy0 *hzp1*mass(pid)
        !$omp atomic
        den8(ibinm1,jbinp1,kbinp1)=den8(ibinm1,jbinp1,kbinp1)+hxm1*hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibin  ,jbinp1,kbinp1)=den8(ibin  ,jbinp1,kbinp1)+hx0 *hyp1*hzp1*mass(pid)
        !$omp atomic
        den8(ibinp1,jbinp1,kbinp1)=den8(ibinp1,jbinp1,kbinp1)+hxp1*hyp1*hzp1*mass(pid)
        w1=w1+mass(pid)
        w2=w2+mass(pid)**2
      enddo
      !$omp end parallel do
    else
      write(*,*) 'ERROR: mass assignment not specified'
    endif
    write(*,*) 'statistics:'
    tmp1=0.d0
    tmp2=0.d0
    !$omp parallel do default(shared) reduction(+:tmp1,tmp2) 
    do k=1,LL(3)
      den8(:,:,k)=den8(:,:,k)/np*vol-1.d0
      tmp1=tmp1+sum(den8(:,:,k))/vol
      tmp2=tmp2+sum(den8(:,:,k)**2)/vol
    enddo
    !$omp end parallel do
    tmp2=sqrt(tmp2)
    write(*,*) 'mean:',real(tmp1)
    write(*,*) 'sigma:',real(tmp2)
    npeff=w1**2/w2
  endsubroutine massassignw_long_r8
  
endmodule p2grid
  
