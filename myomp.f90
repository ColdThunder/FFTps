module myomp
implicit none

interface ompcopy3d
  module procedure ompcopy3dr4r4,ompcopy3dr8r4,ompcopy3dr4r8
  module procedure ompcopy3dc4c4
endinterface ompcopy3d

interface ompplus3d
  module procedure ompplus3dr4r4,ompplus3dr4r8
endinterface ompplus3d

interface ompminus3d
  module procedure ompminus3dr4r4
endinterface ompminus3d

interface omptimes3d
  module procedure omptimes3dr4r4f
  module procedure omptimes3dr4r4
endinterface omptimes3d

interface ompscale
  module procedure ompscalelongr4,ompscaleshortr4
endinterface ompscale

interface ompshift
  module procedure ompshiftlongr4,ompshiftshortr4
endinterface ompshift

interface ompsigma
  module procedure ompsigmalongr4
endinterface ompsigma

interface ompminus
  module procedure ompminuslongr4
endinterface ompminus

interface ompcopy
  module procedure ompcopylongr4
endinterface ompcopy

interface ompset
  module procedure ompsetlongr4
endinterface ompset

contains


!!!!! =====================================
!!!!! OMPPLUS3D
!!!!! =====================================
subroutine ompplus3dr4r8(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(8)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared)
do k=1,LL(3)
  den1(:,:,k)=den1(:,:,k)+den2(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompplus3dr4r8

subroutine ompplus3dr4r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared)
do k=1,LL(3)
  den1(:,:,k)=den1(:,:,k)+den2(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompplus3dr4r4


!!!!! =====================================
!!!!! OMPMINUS3D
!!!!! =====================================
subroutine ompminus3dr4r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared)
do k=1,LL(3)
  den1(:,:,k)=den1(:,:,k)-den2(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompminus3dr4r4

!!!!! =====================================
!!!!! OMPTIMES3D
!!!!! =====================================
subroutine omptimes3dr4r4f(den1,a,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(4)::a
integer(4)::k
!$omp parallel do default(shared)
do k=1,LL(3)
  den1(:,:,k)=den1(:,:,k)*a
enddo
!$omp end parallel do
endsubroutine omptimes3dr4r4f

subroutine omptimes3dr4r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared)
do k=1,LL(3)
  den1(:,:,k)=den1(:,:,k)*den2(:,:,k)
enddo
!$omp end parallel do
endsubroutine omptimes3dr4r4


!!!!! =====================================
!!!!! OMPSHIFT
!!!!! =====================================
subroutine ompshiftshortr4(pos,col,row,a,arange)
implicit none
integer(4)::col
integer(4)::row
real(4)::pos(col,row)
real(4)::a(col)
real(4)::arange(2,col)
real(4)::asize(col)
integer(4)::pid
asize=arange(2,:)-arange(1,:)
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos(1:col,pid)=pos(1:col,pid)+a(1:col)
  where (pos(1:col,pid).le.arange(1,1:col)) pos(1:col,pid)=pos(1:col,pid)+asize(1:col)
  where (pos(1:col,pid).gt.arange(2,1:col)) pos(1:col,pid)=pos(1:col,pid)-asize(1:col)
enddo
!$omp end parallel do
endsubroutine ompshiftshortr4

subroutine ompshiftlongr4(pos,col,row,a,arange)
implicit none
integer(4)::col
integer(8)::row
real(4)::pos(col,row)
real(4)::a(col)
real(4)::arange(2,col)
real(4)::asize(col)
integer(8)::pid
asize=arange(2,:)-arange(1,:)
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos(1:col,pid)=pos(1:col,pid)+a(1:col)
  where (pos(1:col,pid).le.arange(1,1:col)) pos(1:col,pid)=pos(1:col,pid)+asize(1:col)
  where (pos(1:col,pid).gt.arange(2,1:col)) pos(1:col,pid)=pos(1:col,pid)-asize(1:col)
enddo
!$omp end parallel do
endsubroutine ompshiftlongr4

!!!!! =====================================
!!!!! OMPSCALE
!!!!! =====================================
subroutine ompscaleshortr4(pos,col,row,a)
implicit none
integer(4)::col
integer(4)::row
real(4)::pos(col,row)
real(4)::a(col)
integer(4)::pid
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos(1:col,pid)=pos(1:col,pid)*a(1:col)
enddo
!$omp end parallel do
endsubroutine ompscaleshortr4

subroutine ompscalelongr4(pos,col,row,a)
implicit none
integer(4)::col
integer(8)::row
real(4)::pos(col,row)
real(4)::a(col)
integer(8)::pid
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos(1:col,pid)=pos(1:col,pid)*a(1:col)
enddo
!$omp end parallel do
endsubroutine ompscalelongr4


!!!!! =====================================
!!!!! OMPCOPY3D
!!!!! =====================================
subroutine ompcopy3dr8r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(8)::den1(LL(1),LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(:,:,k)=den1(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompcopy3dr8r4

subroutine ompcopy3dr4r8(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(8)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(:,:,k)=den1(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompcopy3dr4r8

subroutine ompcopy3dr4r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(:,:,k)=den1(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompcopy3dr4r4

subroutine ompcopy3dc4c4(den1,den2,LL)
implicit none
integer(4)::LL(3)
complex(4)::den1(LL(1),LL(2),LL(3))
complex(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(:,:,k)=den1(:,:,k)
enddo
!$omp end parallel do
endsubroutine ompcopy3dc4c4

subroutine ompcopy3dc4r4(den1,den2,LL)
implicit none
integer(4)::LL(3)
complex(4)::den1(LL(1)/2,LL(2),LL(3))
real(4)::den2(LL(1),LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(1:LL(1):2,:,k)=real(den1(:,:,k))
  den2(2:LL(1):2,:,k)=aimag(den1(:,:,k))
enddo
!$omp end parallel do
endsubroutine ompcopy3dc4r4

subroutine ompcopy3dr4c4(den1,den2,LL)
implicit none
integer(4)::LL(3)
real(4)::den1(LL(1),LL(2),LL(3))
complex(4)::den2(LL(1)/2,LL(2),LL(3))
integer(4)::k
!$omp parallel do default(shared) schedule(static)
do k=1,LL(3)
  den2(:,:,k)=cmplx(den1(1:LL(1):2,:,k),den1(2:LL(1):2,:,k))
enddo
!$omp end parallel do
endsubroutine ompcopy3dr4c4


!!!!! =====================================
!!!!! OMPSIGMA
!!!!! =====================================

subroutine ompsigmalongr4(pos,col,row,sigma)
implicit none
integer(4)::col
integer(8)::row
real(4)::pos(col,row)
real(4)::sigma
real(8)::sigma8
integer(8)::pid
sigma8=0.
!$omp parallel do default(shared) reduction(+:sigma8) schedule(static)
do pid=1,row
  sigma8=sigma8+sum(pos(1:col,pid)**2)
enddo
!$omp end parallel do
sigma=sqrt(sigma8/row)
endsubroutine ompsigmalongr4


!!!!! =====================================
!!!!! OMPMINUS
!!!!! =====================================

subroutine ompminuslongr4(pos1,pos2,col,row)
implicit none
integer(4)::col
integer(8)::row
real(4)::pos1(col,row)
real(4)::pos2(col,row)
integer(8)::pid
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos1(1:col,pid)=pos1(1:col,pid)-pos2(1:col,pid)
enddo
!$omp end parallel do
endsubroutine ompminuslongr4

!!!!! =====================================
!!!!! OMPCOPY
!!!!! =====================================

subroutine ompcopylongr4(pos1,pos2,col,row)
implicit none
integer(4)::col
integer(8)::row
real(4)::pos1(col,row)
real(4)::pos2(col,row)
integer(8)::pid
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos2(1:col,pid)=pos1(1:col,pid)
enddo
!$omp end parallel do
endsubroutine ompcopylongr4

!!!!! =====================================
!!!!! OMPSET
!!!!! =====================================

subroutine ompsetlongr4(pos,row,ff)
implicit none
integer(8)::row
real(4)::pos(row)
real(4)::ff
integer(8)::pid
!$omp parallel do default(shared) schedule(static)
do pid=1,row
  pos(pid)=ff
enddo
!$omp end parallel do
endsubroutine ompsetlongr4

endmodule myomp
