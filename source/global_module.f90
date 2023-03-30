      module global_module
!**********************************************************************
!     SHARP Pack module for defining global variables
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      real*8,parameter    :: hbar=1.d0
      complex*16,parameter:: eye=(0.,1.)

      integer             :: np     
      real*8              :: dt    
      real*8              :: dtq

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

      logical             :: llan,lnhc
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
      real*8,parameter    :: pi=3.14159265d0  !DACOS(-1.d0)
      real*8              :: beta
      real*8              :: mp 
      real*8              :: KT

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
