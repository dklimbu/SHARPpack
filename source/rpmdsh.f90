PROGRAM ring_polymer_surface_hopping

  use global_module
  use sysdef_module
  use modelvar_module
  use models_module
  use initial_module
  use nonadiabatic_module
  use propagation_module
  use runtraj_module
  use print_module

  implicit none

  integer  :: i, icpu
  integer  :: iseed,seed_dimension,time
 
  real*8   :: RANF
  real*8   :: tt0,tt1

  logical  :: lkval

  integer, allocatable :: seed(:)

  call CPU_TIME(tt0)

! read input parameters from "input.in" file
  call sysdef(lkval)

! set default electronic time step if not assigned as input
  if(dtq.eq.0.d0)dtq=dt/10.d0

! adjust nTraj based on nCPUs
  ntraj = int(real(ntraj)/real(ncpu))

! Initialize RNG seed
  call RANDOM_SEED(size=seed_dimension)
  allocate(seed(seed_dimension))
  do i=1,seed_dimension
     seed(i) =time()+3*i-1
  enddo
  if(ncpu >1)then
    open(1,file='icpu.in',status='old')
    read(1,*) icpu
    close(1)
    seed = seed + icpu*10
  endif

  call RANDOM_SEED(PUT=seed)

! assign some model specific parameters like mass, etc
  call modelParam(model)

! allocate variables and initialize values
  call modelallocat()

! open file for writing purpose into the files
  call openfile()

!=======================initialization==================================
 
! print input parameters as "param.out"
  call printin()

! main trajectory loop
  call runTraj()

! writing hopping statistics
  call hopping_stat()

  call CPU_TIME(tt1)
  if(ldtl)then
    write(nrite_dcoup,*) '#RUN TIME: ',(tt1-tt0), 'seconds'
    write(nrite_dcoup,*) '#RUN TIME: ',(tt1-tt0)/ntraj, 'seconds/traj'
  endif

! close all open files
  call close_file()

! print output results
  call printresults(lkval)

! deallocate the arrays
  call modeldeallocat()

END PROGRAM ring_polymer_surface_hopping

