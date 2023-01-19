program FFTps
use omp_lib
use myomp
use p2grid
use head
implicit none

integer(4)::LLk(3)
real(4),allocatable::pos(:,:)
real(8),allocatable::den8(:,:,:)
real(4),allocatable::vc1(:,:,:)
real(4),allocatable::vc2(:,:,:)
real(4),allocatable::vc3(:,:,:)

real(4)::boxinfo(2,3)
real(4)::fsvec(3)

integer(4)::i,j,k,dir

integer(kind(ngal))::pid

real(8),allocatable::pk8(:,:)
real(8),allocatable::pk2d(:,:,:)
real(4)::fs
real(4)::npeff
integer(4)::pp
logical(4)::debug=.false.

real(4),parameter::pi=4.*atan(1.)
real(4)::kbasic(3)

write(*,*) 'calculating power spectrum using FFT'

call omp_set_num_threads(NUM_THREADS)

call basiccheck
call binscheme 

LLk=[LL(1)+2,LL(2),LL(3)]
boxinfo=reshape([0,LL(1),0,LL(2),0,LL(3)],[2,3])

call memo(0)

!! galaxy catalog/field 1
if (flag_gcata) then
  allocate(pos(3,ngal))
  call readcata(infg,nghead,ngal,pos,gxyz)
  fsvec=LL/box
  call ompscale(pos,3,ngal,fsvec)

  fsvec=[0.,0.,0.]
  call ompshift(pos,3,ngal,fsvec,boxinfo)
  call massassign(ngal,pos,LL,den8,painter)  ! den8 -> delta_g
  call ompcopy3d(den8,vc1(1:LL(1),:,:),LL)  ! vc1 -> delta_g

  if (interlace) then
    fsvec=[0.5,0.5,0.5]
    call ompshift(pos,3,ngal,fsvec,boxinfo)
    call massassign(ngal,pos,LL,den8,painter)  ! den8 -> delta_g(interlace)
    call ompcopy3d(den8,vc3(1:LL(1),:,:),LL)  ! vc3 -> delta_g(interlace)
  endif

  deallocate(pos)
else
  call readfield(infg,LL,vc1(1:LL(1),:,:))  ! vc1 -> delta_g, if the input is field
endif
! anyway (input catalog or field), vc1 -> delta_g, vc3 -> delta_g (interlace)

!! random catalog/field 1
if (flag_r) then
  if (flag_rcata) then
    allocate(pos(3,nran))
    call readcata(infr,nrhead,nran,pos,rxyz)
    fsvec=LL/box
    call ompscale(pos,3,nran,fsvec)

    fsvec=[0.,0.,0.]
    call ompshift(pos,3,nran,fsvec,boxinfo)
    call massassign(nran,pos,LL,den8,painter)  ! den8 -> delta_r
    call ompminus3d(vc1(1:LL(1),:,:),den8,LL)  ! vc1 -> delta_g-delta_r

    if (interlace) then
      fsvec=[0.5,0.5,0.5]
      call ompshift(pos,3,nran,fsvec,boxinfo)
      call massassign(nran,pos,LL,den8,painter)  ! den8 -> delta_r (interlace)
      call ompminus3d(vc3(1:LL(1),:,:),den8,LL)  ! vc3 -> delta_g-delta_r (interlace)
    endif

    deallocate(pos)
  else
    call readfield(infr,LL,vc3(1:LL(1),:,:))  ! vc3 -> delta_r, if the input is field
    call ompminus3d(vc1(1:LL(1),:,:),vc3(1:LL(1),:,:),LL) ! vc1 -> delta_g-delta_r
  endif
  ! anyway (input catalog or field), vc1 -> delta_g-delta_r, vc3 -> delta_g-delta_r (interlace)

endif

call fftr2c_inplace(vc1,LL)  ! vc1 -> delta(k)
if (interlace) then
  call fftr2c_inplace(vc3,LL)  ! vc3 -> delta(k) (interlace)
  fsvec=box/LL/2
  call phaseshift(vc3,LL,fsvec)  ! correct interlace phase shift
  call ompplus3d(vc1,vc3,LLk)
  call omptimes3d(vc1,0.5,LLk)  ! vc1 -> delta(k) (interlaced)
endif
! anyway (interlace or not), vc1 -> delta(k), vc3 is free

!!! second field
if (second) then
  !! galaxy catalogue/field 2
  if (flag_gcata2) then
    allocate(pos(3,ngal2))
    call readcata(infg2,nghead2,ngal2,pos,gxyz2)
    fsvec=LL/box
    call ompscale(pos,3,ngal2,fsvec)
  
    fsvec=[0.,0.,0.]
    call ompshift(pos,3,ngal2,fsvec,boxinfo)
    call massassign(ngal2,pos,LL,den8,painter2)  ! den8 -> delta_g (2nd)
    call ompcopy3d(den8,vc2(1:LL(1),:,:),LL)  ! vc2 -> delta_g (2nd)
  
    if (interlace2) then
      fsvec=[0.5,0.5,0.5]
      call ompshift(pos,3,ngal2,fsvec,boxinfo)
      call massassign(ngal2,pos,LL,den8,painter2)  ! den8 -> delta_g (2nd, interlace)
      call ompcopy3d(den8,vc3(1:LL(1),:,:),LL)  ! vc3 -> delta_g (2nd, interlace)
    endif
  
    deallocate(pos)
  else
    call readfield(infg2,LL,vc2(1:LL(1),:,:)) ! vc2 -> delta_g (2nd), if the input is field
  endif
  ! anyway (input catalog or field), vc2 -> delta_g (2nd), vc3 -> delta_g (2nd, interlace)
  
  !! random catalog/field 2
  if (flag_r2) then
    if (flag_rcata2) then
      allocate(pos(3,nran2))
      call readcata(infr2,nrhead2,nran2,pos,rxyz2)
      fsvec=LL/box
      call ompscale(pos,3,nran2,fsvec)
  
      fsvec=[0.,0.,0.]
      call ompshift(pos,3,nran2,fsvec,boxinfo)
      call massassign(nran2,pos,LL,den8,painter2)  ! den8 -> delta_r (2nd)
      call ompminus3d(vc2(1:LL(1),:,:),den8,LL)  ! vc2 -> delta_g-delta_r (2nd)
  
      if (interlace2) then
        fsvec=[0.5,0.5,0.5]
        call ompshift(pos,3,nran2,fsvec,boxinfo)
        call massassign(nran2,pos,LL,den8,painter2)  ! den8 -> delta_r (2nd, interlace)
        call ompminus3d(vc3(1:LL(1),:,:),den8,LL) ! vc3 -> delta_g-delta_r (2nd, interlace)
      endif
  
      deallocate(pos)
    else
      call readfield(infr2,LL,vc3(1:LL(1),:,:))
      call ompminus3d(vc2(1:LL(1),:,:),vc3(1:LL(1),:,:),LL)
    endif
    ! anyway (input catalog or field), vc2 -> delta_g-delta_r (2nd), vc3 -> delta_g-delta_r (2nd, interlace)
  
  endif

  call fftr2c_inplace(vc2,LL)  ! vc2 -> delta(k) (2nd)
  if (interlace2) then
    call fftr2c_inplace(vc3,LL)  ! vc3 -> delta(k) (2nd, interlace)
    fsvec=box/LL/2
    call phaseshift(vc3,LL,fsvec)  ! correct interlace phase shift
    call ompplus3d(vc2,vc3,LLk)
    call omptimes3d(vc2,0.5,LLk)  ! vc2 -> delta(k) (2nd, interlaced)
  endif
  ! anyway (interlace or not), vc2 -> delta(k) (2nd), vc3 is free
endif

!! all the input is done
!! galaxy catalogue or galaxy+random catalogyes
!! input catalogue or input fields
!! 1 input for 2 inputs
!! interlace or not


!! 3D power spectrum
if (second) then
  call adjustdata(vc1,vc2,vc3,LL)
else
  call adjustdata(vc1,vc1,vc3,LL)
endif


if (sncorr) then
  pp=0
  if (second) then
    if (flag_r) then
      npeff=real(ngal)*nran/(ngal+nran)
    else
      npeff=real(ngal)
    endif
    if (othernpeff) npeff=npeffinput
    if (wincorr) pp=pp+painter
    if (wincorr2) pp=pp+painter2
  else
    if (flag_r) then
      npeff=real(ngal)*nran/(ngal+nran)
    else
      npeff=real(ngal)
    endif
    pp=pp+painter*2
  endif
  if (interlace) then
    call subtractnoise_interlace(vc3,LL,npeff,pp)
  else
    call subtractnoise(vc3,LL,npeff,pp)
  endif
else
  write(*,*) 'shot noise is not corrected'
endif

pp=0
if (second) then
  if (wincorr) pp=pp+painter
  if (wincorr2) pp=pp+painter2
else
  if (wincorr) pp=pp+painter*2
endif
call normalization(vc3,LL,pp)

call calcpkell(vc3,LL,pk8,kbin,outname)

if (flag_pk2d) call calcpkmu(vc3,LL,pk2d,kbin,ubin,outname2)

call memo(1)

contains

subroutine basiccheck
  implicit none
  character(512)::message

  ! interlace conflicts with input density field
  message='ERROR: interlace can only be used on catalogue'
  if (interlace.and. (.not.flag_gcata)) then
    write(*,*) trim(message)//' (The input is a galaxy density FIELD 1)'
    stop
  endif
  if (interlace.and. (.not.flag_rcata)) then
    write(*,*) trim(message)//' (The input is a random density FIELD 1)'
    stop
  endif
  if (second.and.interlace2.and. (.not.flag_gcata2)) then
    write(*,*) trim(message)//' (The input is a galaxy density FIELD 2)'
    stop
  endif
  if (second.and.interlace2.and. (.not.flag_rcata2)) then
    write(*,*) trim(message)//' (The input is a random density FIELD 2)'
    stop
  endif
endsubroutine basiccheck

subroutine binscheme
  implicit none
  kbasic=2*pi/box
  write(*,'(a,3f8.3)') 'kbasic=',kbasic
  write(*,'(a,3f8.3)') 'knyquist=',kbasic*LL/2
  if (kbin.eq.0.and.dk.eq.0.) then
    write(*,*) 'ERROR: kbin and dk cannot to 0 in the same time'
    stop
  endif

  if (usedefaultkminkmax) then
    if (logbin) then
      kmax=alog10(maxval(kbasic*LL/2))
      kmin=alog10(minval(kbasic)*1.28)
    else
      kmax=maxval(kbasic*(LL+1)/2)
      kmin=minval(kbasic/2)
    endif
  else
    if (logbin) then
      kmax=alog10(kmax)
      if (kmin.le.0) then
        write(*,*) 'ERROR: kmin=0, set kmin=1.28*kbasic'
        kmin=alog10(minval(kbasic)*1.28)
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
integer(4)::npmax
if (command.eq.0) then

  if (flag_gcata.or.flag_gcata2.or.flag_rcata.or.flag_rcata2) then
    npmax=0
    if (flag_gcata) npmax=max(npmax,ngal)
    if (flag_rcata) npmax=max(npmax,nran)
    if (flag_gcata2) npmax=max(npmax,ngal2)
    if (flag_rcata2) npmax=max(npmax,nran2)
    write(*,*) 'memo for pos:',(float(npmax)*3)*4/1024.**3,'G'
  endif

  write(*,*) 'memo for den8:',product(real(LL))*8/1024.**3,'G'
  allocate(den8(LL(1),LL(2),LL(3)))

  write(*,*) 'memo for vc1:',(float(LL(1)+2)*LL(2)*LL(3))*4/1024.**3,'G'
  allocate(vc1(LL(1)+2,LL(2),LL(3)))

  if (second) then
    write(*,*) 'memo for vc2:',(float(LL(1)+2)*LL(2)*LL(3))*4/1024.**3,'G'
    allocate(vc2(LL(1)+2,LL(2),LL(3)))
  endif

  write(*,*) 'memo for vc3:',(float(LL(1)+2)*LL(2)*LL(3))*4/1024.**3,'G'
  allocate(vc3(LL(1)+2,LL(2),LL(3)))

  allocate(pk8(9,kbin))
  if (flag_pk2d) allocate(pk2d(7,kbin,ubin))
elseif(command.eq.1) then
  deallocate(den8)
  deallocate(vc1)
  if (second) deallocate(vc2)
  deallocate(vc3)

  deallocate(pk8)
  if (flag_pk2d) deallocate(pk2d)
endif
endsubroutine memo

subroutine readcata(filename,nhead,np,pos,col)
  implicit none
  character(*)::filename
  integer(4)::nhead
  integer(4)::np
  real(4)::pos(3,np)
  integer(4)::col(3)

  integer(4)::colmax
  real(4),allocatable::line(:)
  ! read in galaxy
  colmax=maxval(col)
  allocate(line(colmax))
  write(*,*) 'reading: ',trim(filename)
  open(31,file=filename,status='old',form='formatted')
  do pid=1,nhead
    read(31,*)
  enddo
  do pid=1,np
    read(31,*) line
    pos(1:3,pid)=line(col(1:3))
  enddo
  close(31)
  deallocate(line)
endsubroutine readcata

subroutine readfield(filename,LL,den)
  implicit none
  character(*)::filename
  integer(4)::LL(3)
  real(4)::den(LL(1),LL(2),LL(3))

  write(*,*) 'reading: ',trim(filename)
  open(31,file=filename,status='old',access='stream')
  read(31) den(1:LL(1),1:LL(2),1:LL(3))
  close(31)
endsubroutine readfield

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
  rk=kbasic(3)*rk
  do j=1,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    rj=kbasic(2)*rj
    do i=1,LL(1)+2,2
      ri=i/2
      ri=kbasic(1)*ri

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

pk8(3:8,:)=pk8(3:8,:)*product(box)

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


subroutine calcpkmu(vc,LL,pk2d,kbin,ubin,filename)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  integer(4)::kbin
  integer(4)::ubin
  character(*)::filename
  real(8)::pk2d(7,kbin,ubin)
  real(4)::ri,rj,rk,kmode,mu,du
  integer(4)::bin1,bin2
  real(4)::factor
  real(4)::kvalue(kbin)
  real(4)::uvalue(ubin)

  du=(umax-umin)/ubin
  pk2d=0.
  write(*,*) 'calculating 2D power spectrum'
  write(*,*) 'kmin,kmax=',kmin,kmax
  write(*,*) 'kbin=',kbin
  write(*,*) 'dk=',dk
  write(*,*) 'umin,umax=',umin,umax
  write(*,*) 'ubin=',ubin
  write(*,*) 'du=',du
  !$omp parallel do default(private) shared(vc,LL,kbasic,kmin,umin,dk,du,kbin,ubin,logbin) &
  !$omp reduction(+:pk2d) schedule(guided)
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    rk=kbasic(3)*rk
    do j=1,LL(2)
      if (j.le.LL(2)/2+1) then
        rj=j-1
      else
        rj=j-1-LL(2)
      endif
      rj=kbasic(2)*rj
      do i=1,LL(1)+2,2
        ri=i/2
        ri=kbasic(1)*ri

        if (i.eq.1.and.j.eq.1.and.k.eq.1) cycle
        if (i.eq.1) then
          factor=1
        else
          factor=2
        endif
        kmode=sqrt(ri*ri+rj*rj+rk*rk)
        if (logbin) kmode=alog10(kmode)
        mu=rk/kmode
        bin1=floor((kmode-kmin)/dk)+1
        bin2=floor((mu-umin)/du)+1
        if (bin1.lt.1.or.bin1.gt.kbin) cycle
        if (bin2.lt.1.or.bin2.gt.ubin) cycle
        pk2d(1,bin1,bin2)=pk2d(1,bin1,bin2)+kmode*factor
        pk2d(2,bin1,bin2)=pk2d(2,bin1,bin2)+kmode**2*factor
        pk2d(3,bin1,bin2)=pk2d(3,bin1,bin2)+mu*factor
        pk2d(4,bin1,bin2)=pk2d(4,bin1,bin2)+mu**2*factor
        pk2d(5,bin1,bin2)=pk2d(5,bin1,bin2)+vc(i,j,k)*factor
        pk2d(6,bin1,bin2)=pk2d(6,bin1,bin2)+vc(i,j,k)**2*factor
        pk2d(7,bin1,bin2)=pk2d(7,bin1,bin2)+1.*factor
      enddo
    enddo
  enddo
  !$omp end parallel do
  
  pk2d(1,:,:)=pk2d(1,:,:)/pk2d(7,:,:)
  pk2d(2,:,:)=pk2d(2,:,:)/pk2d(7,:,:)-pk2d(1,:,:)**2
  pk2d(3,:,:)=pk2d(3,:,:)/pk2d(7,:,:)
  pk2d(4,:,:)=pk2d(4,:,:)/pk2d(7,:,:)-pk2d(3,:,:)**2
  pk2d(5,:,:)=pk2d(5,:,:)/pk2d(7,:,:)
  pk2d(6,:,:)=pk2d(6,:,:)/pk2d(7,:,:)-pk2d(5,:,:)**2

  do i=1,kbin
    kvalue(i)=kmin+(i-0.5)*dk
  enddo
  do i=1,ubin
    uvalue(i)=umin+(i-0.5)*du
  enddo

  if (logbin) kvalue=10.**kvalue
  
  write(*,*) 'writing: ',trim(filename)
  open(32,file=filename,form='formatted',status='replace')
  do i=1,kbin
    do j=1,ubin
      write(32,'(4e14.6)') kvalue(i),uvalue(j),pk2d(5,i,j),pk2d(7,i,j)
      if (debug) write(*,*) kvalue(i),uvalue(j),real(pk2d(5,i,j)),pk2d(7,i,j)
    enddo
  enddo
  close(32)

endsubroutine calcpkmu


subroutine adjustdata(vc1,vc2,vc3,LL)
  implicit none
  integer(4)::LL(3)
  real(4)::vc1(LL(1)+2,LL(2),LL(3))
  real(4)::vc2(LL(1)+2,LL(2),LL(3))
  real(4)::vc3(LL(1)+2,LL(2),LL(3))
  real(4)::vctemp(2,LL(2))
  real(4)::factor
  factor=1./LL(1)/LL(2)/LL(3)
  !$omp parallel do default(shared) schedule(static) private(i,k,vctemp)
  do k=1,LL(3)
    do i=1,LL(1)+2,2
      vctemp(1,:)=vc1(i,:,k)*vc2(i,:,k)+vc1(i+1,:,k)*vc2(i+1,:,k)
      vctemp(2,:)=vc2(i,:,k)*vc1(i+1,:,k)-vc1(i,:,k)*vc2(i+1,:,k)
      vc3(i:i+1,:,k)=vctemp(1:2,:)*factor**2
    enddo
  enddo
  !$omp end parallel do
endsubroutine adjustdata

subroutine subtractnoise(vc,LL,npeff,pp)
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
        vc(i:i+1,j,k)=vc(i:i+1,j,k)-d2/npeff
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine subtractnoise

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

subroutine subtractnoise_interlace(vc,LL,npeff,pp)
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
    wz(k,1:2)=wcorr_interlace(a3,pp)
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
    wy(j,1:2)=wcorr_interlace(a2,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(LL,hpi,wx,pp) schedule(static)
  do i=1,LL(1)+2,2
    ri=i/2
    a1=hpi(1)*ri/2
    wx(i/2+1,1:2)=wcorr_interlace(a1,pp)
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(vc,LL,wx,wy,wz,npeff) schedule(static)
  do k=1,LL(3)
    do j=1,LL(2)
      do i=1,LL(1)+2,2
        d2= wx(i/2+1,1)*wy(j,1)*wz(k,1) + wx(i/2+1,1)*wy(j,2)*wz(k,2) + wx(i/2+1,2)*wy(j,1)*wz(k,2) + wx(i/2+1,2)*wy(j,2)*wz(k,1)
        vc(i:i+1,j,k)=vc(i:i+1,j,k)-d2/npeff
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine subtractnoise_interlace

function wcorr_interlace(aa,pp)
  implicit none
  real(4)::wcorr_interlace(2)
  real(4)::aa
  integer(4)::pp
  real(4)::sinx,cosx
  sinx=sin(aa)
  cosx=cos(aa)

  if (pp.eq.2) then
    wcorr_interlace(1)=cosx**2
    wcorr_interlace(2)=sinx**2
  !elseif (pp.eq.3) then
  !  wcorr_interlace_new(1)=coshx**4
  elseif (pp.eq.4) then
    wcorr_interlace(1)=(1-2./3*sinx**2)*cosx**4
    wcorr_interlace(2)=(1-2./3*cosx**2)*sinx**4
  !elseif (pp.eq.5) then
  !  wcorr_interlace_new(1)=(1-1./3*sinhx**2)*coshx**6
  elseif (pp.eq.6) then
    wcorr_interlace(1)=(1-sinx**2+2./15*sinx**4)*cosx**6
    wcorr_interlace(2)=(1-cosx**2+2./15*cosx**4)*sinx**6
  !elseif (pp.eq.7) then
  !  wcorr_interlace_new(1)=(1-2./3*sinhx**2+2./45*sinhx**4)*coshx**8
  elseif(pp.eq.8) then
    wcorr_interlace(1)=(1-4./3*sinx**2+2./5*sinx**4-4./315*sinx**6)*cosx**8
    wcorr_interlace(2)=(1-4./3*cosx**2+2./5*cosx**4-4./315*cosx**6)*sinx**8
  else
    wcorr_interlace(1:2)=1
  endif
endfunction wcorr_interlace

subroutine normalization(vc,LL,pp)
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

  wz(1)=1
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

  wy(1)=1
  !$omp parallel do default(private) shared(hpi,LL,wy) schedule(static)
  do j=2,LL(2)
    if (j.le.LL(2)/2+1) then
      rj=j-1
    else
      rj=j-1-LL(2)
    endif
    a2=hpi(2)*rj
    wy(j)=sin(a2)/a2
  enddo
  !$omp end parallel do

  wx(1)=1
  !$omp parallel do default(private) shared(hpi,LL,wx) schedule(static)
  do i=3,LL(1)+2,2
    ri=i/2
    a1=hpi(1)*ri
    wx(i/2+1)=sin(a1)/a1
  enddo
  !$omp end parallel do

  !$omp parallel do default(private) shared(vc,wx,wy,wz,pp,LL) schedule(static)
  do k=1,LL(3)
    do j=1,LL(2)
      do i=1,LL(1)+2,2
        d2=(wx(i/2+1)*wy(j)*wz(k))**pp
        vc(i:i+1,j,k)=vc(i:i+1,j,k)/d2
      enddo
    enddo
  enddo
  !$omp end parallel do
endsubroutine normalization

subroutine phaseshift(vc,LL,x0)
  implicit none
  integer(4)::LL(3)
  real(4)::vc(LL(1)+2,LL(2),LL(3))
  complex(4)::cc
  real(4)::x0(3)
  real(4)::ri,rj,rk
  integer(4)::i,j,k
  do k=1,LL(3)
    if (k.le.LL(3)/2+1) then
      rk=k-1
    else
      rk=k-1-LL(3)
    endif
    rk=kbasic(3)*rk
    do j=1,LL(2)
      if (j.le.LL(2)/2+1) then
        rj=j-1
      else
        rj=j-1-LL(2)
      endif
      rj=kbasic(2)*rj
      do i=1,LL(1)+2,2
        ri=i/2
        ri=kbasic(1)*ri
  
        cc=cmplx(vc(i,j,k),vc(i+1,j,k))
        cc=cc*exp(sum(x0*[ri,rj,rk])*cmplx(0.,1.))
        vc(i:i+1,j,k)=[real(cc),aimag(cc)]
      enddo
    enddo
  enddo
  endsubroutine phaseshift


endprogram FFTps

include 'mkl_dfti.f90'
include 'intel_3d_fft.f90'
