  program average2db
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PROGRAM TO CALCULATE RUNNING AVERAGE OF POPULATION IN DETAILED
!  BALANCE MODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  integer  :: i,io
  integer  :: df,ncount,narg

  real(8)  :: time,tstart
  real(8),allocatable :: sumx(:),x(:)

  character(len=160) :: dfile
  character(len=30) :: argv2,argv3

  !write(0,*), 'Data file name?'
  !read(*,*) dfile
  !write(0,*), 'No. of data column(s) (except first column)?'
  !read(*,*) df
  !write(0,*), 'Skip first n ps step data?'
  !read(*,*) tstart

  !dfile='adiabat1.txt'
  !df = 2
  !tstart = 30.0

  call get_command_argument(1,dfile)
  call get_command_argument(2,argv2)
  call get_command_argument(3,argv3)
  
  narg = command_argument_count()

  if(narg .ne. 3)then
    write(0,*)
    write(0,*) ' USAGE:: ./average2db.x dfile df tstart > outfile.xyz'
    write(0,*) 
    stop
  endif

  read(argv2,*) df
  read(argv3,*) tstart

  write(0,'(3A,f8.2)') 'file: ',trim(dfile),' tstart:',tstart

  tstart = tstart * 41340 !! (ps into a.u.)

  allocate(sumx(df),x(df))

  open(2,file=dfile,status='old')

  sumx = 0.d0
  ncount = 0
  read(2,*)
  read(2,*)
  do !istep = 1, tstep
    read(2,*,End=1) time,x(1:df)
    if(io > 0) stop
    if(time > tstart)then
      sumx = sumx + x
      ncount = ncount + 1
    endif
  enddo

 1 write(0,*) 'Completed!!!'
   
   close(2)

   sumx = sumx/ncount
   write(*,'(A,20f15.6)') 'adiabat1.txt',sumx

    end program
