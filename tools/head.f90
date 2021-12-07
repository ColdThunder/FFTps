
module head
implicit none

real(4),parameter::box=100.  ! box size in Mpc/h
integer(4),parameter::L=256   ! cells per dimension

character(512)::outsnname='./sncorr256-cic.txt'
character(512)::outwinname='./wincorr256-cic.txt'

logical(4)::sncorr=.true.   ! whether correct shot noise

logical(4)::interlace=.false.
integer(4)::painter=2  ! 1 for cic
logical(4)::wincorr=.true.   ! whether deconvolve window for field1

logical(4)::usedefaulkminkmax=.true.
real(4)::kmin=0.    ! input kmin kmax if usedefault=false
real(4)::kmax=0.    ! input kmin kmax if usedefault=false

logical(4)::logbin=.false.
integer(4)::kbin=64   ! kbin and dk can not be zero simuteneously
real(4)::dk=0.         ! kbin and dk can not be zero simuteneously
logical(4)::usecenter=.true.

integer(4)::NUM_THREADS=32

logical(4)::flag_pk2d=.false.
character(512)::outname2='./output/pk2D-rsd.txt'
real(4)::umin=0.
real(4)::umax=1.
integer(4)::ubin=5

!!!!!! second field
logical(4)::second=.false.

logical(4)::interlace2=.false.
integer(4)::painter2=2   ! 2 for cic
logical(4)::wincorr2=.true.   ! whether deconvolve window for field2, for auto, field2=field1

!!!!! other shot noise choice
logical(4)::othernpeff=.false.  ! set effective particle number to correct shot noise
real(4)::npeffinput=0.          ! 


endmodule head
