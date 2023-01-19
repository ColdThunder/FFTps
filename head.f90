module head
implicit none

real(4)::box(3)=[1000.,1000.,1000. ] ! box size in Mpc/h
integer(4)::LL(3)=[512,512,512]   ! cells per dimension

character(512)::outname='./output/pkl-rsd.txt'

character(512)::infg='./data/ELG-wpmax-v3-snap97-redshift0.99_dens2.dat'  ! galaxy catalogue name
logical(4)::flag_gcata=.true.  ! input data is a catalogue, otherwise, a field
integer(4)::gxyz(3)=[1,2,4]     ! columns for x y z
integer(4)::nghead=0            ! number of lines for header
integer(4)::ngal=2500000          ! number of lines for galaxy

logical(4)::flag_r=.true.      ! use random or not
character(512)::infr='./data/randomsx10_seed0.txt'    ! random catalgue name
logical(4)::flag_rcata=.true.  ! input random is a catalogue, otherwise, a field
integer(4)::rxyz(3)=[1,2,3]       ! columns for x y z
integer(4)::nrhead=0              ! number of lines for header
integer(4)::nran=32340460        ! number of lines for random

logical(4)::sncorr=.true.   ! whether correct shot noise

logical(4)::interlace=.false.
integer(4)::painter=2   ! 2 for cic
logical(4)::wincorr=.true.   ! whether deconvolve window for field1

logical(4)::usedefaultkminkmax=.true.
real(4)::kmin=0.    ! input kmin kmax if usedefault=false
real(4)::kmax=0.    ! input kmin kmax if usedefault=false

logical(4)::logbin=.true.
integer(4)::kbin=20   ! kbin and dk can not be zero simuteneously
real(4)::dk=0.         ! kbin and dk can not be zero simuteneously
logical(4)::usecenter=.true.

integer(4)::NUM_THREADS=32

logical(4)::flag_pk2d=.true.
character(512)::outname2='./output/pk2D-rsd.txt'
real(4)::umin=0.
real(4)::umax=1.
integer(4)::ubin=5

!!!!!! second field
logical(4)::second=.true.
logical(4)::flag_gcata2=.true.  ! input data is a catalogue, otherwise, a field
character(512)::infg2='./data/ELG-wpmax-v3-snap97-redshift0.99_dens2.dat'  ! galaxy catalogue name
integer(4)::gxyz2(3)=[1,2,4]     ! columns for x y z
integer(4)::nghead2=0            ! number of lines for header
integer(4)::ngal2=2500000          ! number of lines for galaxy

logical(4)::flag_r2=.true.      ! use random or not
logical(4)::flag_rcata2=.true.  ! input random is a catalogue, otherwise, a field
character(512)::infr2='./data/randomsx10_seed0.txt'    ! random catalgue name
integer(4)::rxyz2(3)=[1,2,3]       ! columns for x y z
integer(4)::nrhead2=0              ! number of lines for header
integer(4)::nran2=32340460        ! number of lines for random

logical(4)::interlace2=.false.
integer(4)::painter2=2   ! 2 for cic
logical(4)::wincorr2=.true.   ! whether deconvolve window for field2, for auto, field2=field1

!!!!! other shot noise choice
logical(4)::othernpeff=.false.  ! set effective particle number to correct shot noise
real(4)::npeffinput=0.          ! 


endmodule head
