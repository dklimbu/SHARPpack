      module global_module
      !integer,parameter   :: nlit=20
      !real*8,parameter    :: rnlit = 1.0/REAL(nlit,8)
      real*8,parameter    :: hbar=1.d0
      complex*16,parameter:: eye=(0.,1.)

      integer             :: np     !20
      real*8              :: dt     !=0.01d0
      real*8              :: dtq    != dt*rnlit

      integer             :: model
      integer             :: nstates
      integer             :: nsteps
      integer             :: nsample
      integer             :: ntraj
      integer             :: ntrajR
      integer             :: nb
      integer             :: iskip 
      integer             :: vinit
      integer             :: Rinit
      real*8              :: P0
      real*8              :: R0 
      real*8              :: omega
      character(len=36)   :: modelname
      character(len=8)    :: method

      logical             :: llan,lnhc,lpcet
      logical             :: ldtl
  
      integer,parameter   :: nrite_hopp=110
      integer,parameter   :: nrite_dcoup=120
      integer,parameter   :: nrite_therm=130
      integer,parameter   :: nrite_bath=140
 
      integer             :: vrkey
      integer             :: nfrust_hop 
      integer             :: nfrust_hop2
      integer             :: ncpu 
      integer             :: vckey

! Conversion factors to the atomic units
! Energy from eV, Time from ps, Distance from angstrom
! Temperature from K, Frequency from cm^-1

      real*8,parameter    :: energy=27.21d0,tim=0.00002419d0,dist=0.5291d0
      real*8,parameter    :: temp=315774d0,freq=219463.343d0

! Classical environment parameters
      real*8              :: beta
      real*8              :: mp 
      real*8,parameter    :: pi=3.14159265d0  !DACOS(-1.d0)
!      real*8,parameter :: mp=2000.d0  !1.d0

!! Spin-Boson parameters
!      real*8, parameter :: zxi = 0.13d0
!      real*8, parameter :: wc = 2.5d0
      real*8, parameter :: wmb = 3.0d0   ! wm = wm/wc ~ 3 
      real*8            :: wmax          != wc * (1.0 - exp(-wmb))/np
!      real*8, parameter :: delta = 1.0d0
!      real*8, parameter :: eps = 0.0d0
      real*8             :: zxi, wc, delta, eps

!! Spin-Boson2-Debye parameters
! Diabatic matrix parameters
      real*8, parameter :: enu=208.5d0/freq
      real*8            :: E_r, KT
!      real*8, parameter :: eps=0.d0*enu/freq
!      real*8, parameter :: delta=0.2d0*enu/freq
! Classical environment parameters
!      real*8,parameter :: beta=1.d0/(2.d0*enu/freq)
!      real*8,parameter :: mp=1.d0
!      real*8,parameter :: wc=0.2d0*enu/freq
!      real*8,parameter :: wmax=20.d0*wc
!      real*8,parameter :: E_r=0.02d0*enu/freq

! NVT_PILE parameters
      real*8            :: tau0

! NVT_NOSE_HOOVER_CHAIN parameters
      integer,parameter :: nsy=7
      integer           :: nrespa
      integer           :: nchain
      real*8            :: weight(nsy)
      real*8            :: qmass_t, qmass_part
      real*8            :: taut

!!!!!!PCET MODEL !!!!
      ! Classical environment parameters
      real*8,parameter::f0=55.7d0,eps_inf=4.2d0,eps0=79.2d0
      real*8,parameter :: tauo=0.0103d0/tim
      real*8,parameter :: taud=8.72d0/tim
      real*8,parameter :: lambda=0.65d0/energy
!      real*8,parameter :: delta=0.0d0/energy
!      real*8,parameter :: beta=1.d0/298.d0*temp   
!      real*8,parameter :: mass(1:NDIM*NPART)=0.265/tim**2   !<=mp(water) in sysdef
      real*8,parameter :: taul=eps_inf*(tauo+taud)/eps0
      real*8           :: sigma         ! width of the random noise
!      real*8,parameter :: sigma=dsqrt(2.d0*f0*taul/(dt*beta))  ! width of the random noise
! Quantum subsystem parameters
!      real*8,parameter :: m=1836.d0
!      real*8,parameter :: omega=3000.d0/freq
      real*8,parameter :: qpA=0.5/dist
      real*8,parameter :: Vda=0.03/energy
      integer          :: nbasis         ! nstates=2*nbasis

!!! LinearChainModel !!!
      integer :: nprint
      real*8  :: gamaLC
      real*8  :: sigmaLC
      real*8  :: d_ab(3,3)     
      real*8  :: v11
      real*8  :: v22
      real*8,parameter :: kJ_mol2au=0.00038125d0
      real*8,parameter :: amu2au=1822.4d0

!!! with 3-states LinearChainModel !!!
      real*8  :: v33

      end module global_module
