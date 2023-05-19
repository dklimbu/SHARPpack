!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PROGRAM TO CALCULATE AVERAGE FROM PARALLEL SHARP RUN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program average
  implicit none

  integer  :: i,j,io, iformat
  integer  :: df,ndir,headskip

  real(8)  :: time
  real(8),allocatable :: sumx(:),x(:)

  character(len=80) :: infile,root_dir,dfile,outfil
  character(len=4)  :: idir
  character(len=150)  ::header1,header2

  logical  :: file_exist

  !write(0,*), 'Root directory (without numeric)?'
  !read(*,*) root_dir
  !write(0,*), 'Data file name?'
  !read(*,*) dfile
  !write(0,*), 'No. of directory(s)?'
  !read(*,*) ndir
  !write(0,*), 'No. of data column(s) (except first column)?'
  !read(*,*) df
  !write(0,*), 'No.of heading row(s) skip ?'
  !read(*,*) headskip
  !write(0,*), 'numeric format for root_dir (x:1, xx:2, xxx:3), enter (1/2/3) ?'
  !read(*,*) iformat

  open(1,file='input', status='old')
    read(1,*) root_dir
    read(1,*) dfile
    read(1,*) ndir
    read(1,*) df
    read(1,*) iformat
  close(1)

  outfil = dfile(1:len_trim(dfile)-4)//'_ave.out'

  headskip=2
  allocate(sumx(df),x(df))

  write(0,*)
  write(0,*) '************************************'
  write(0,*), 'Averaging from File(s):'
  do i = 1, ndir

    if(iformat .eq. 1)then
      write(idir,'(I1.1)') i
    elseif(iformat .eq. 2)then
      write(idir,'(I2.2)') i
    elseif(iformat .eq. 3)then
      write(idir,'(I3.3)') i
    else
      write(0,*) 'Wrong numeric format!!'
      stop
    endif

    infile = trim(root_dir)//trim(idir)// '/' //trim(dfile)

    inquire(file=infile,exist=file_exist)
    if(file_exist)then
      open(unit=i,file=infile)
    else
      write(0,*) 'File does not exist!'
      write(0,*) 'Check input file name!'
      write(0,*) 
      write(0,*) '************************************'
      stop
    endif

    write(0,*),i,trim(infile)
    
    !do j = 1, headskip
     read(i,'(A)') header1
     read(i,'(A)') header2
    !enddo

  enddo

  open(1000,file=outfil,status='unknown')
  write(1000,'(A)') trim(header1)
  write(1000,'(A)') trim(header2)
  do !istep = 1, tstep
    sumx = 0.d0
    do i = 1, ndir
      read(i,*,End=1) time,x(1:df)
      if(io > 0) stop
      sumx = sumx + x
    enddo
    sumx = sumx/ndir

    write(1000,'(f13.4,9f15.6)') time, sumx
  enddo

 1 write(0,*) 'Completed!!!'
  close(2)
  write(0,*) 'See Averaged File: ',trim(outfil)
  write(0,*) '************************************'
  write(0,*)

    end program
