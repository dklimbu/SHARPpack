      module sysdef_module
!**********************************************************************
!     SHARP PACK module for defining simulation system     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use parse_module
      use global_module

      contains

      subroutine sysdef(lkval)
!**********************************************************************
!     SHARP PACK routine for reading simulation control input
!     parameters     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
      
      character*1        :: directive(lenrec)
      logical            :: safe, loop, kill
      integer            :: i, idum, idnode, stat

      logical            :: lkval

      integer, parameter :: nread = 1, nrite=2
      character(len=24)  :: datentime

      idnode = 0

      loop = .true.
      ldtl = .false.
      lkval = .false.
      lfft = .false.
      
      np = 1
      dt = 1.0
      dtq = 0.1d0

      keymodel = 1
      nstates = 0
      nsteps = 0
      nsample = 0
      ntraj = 0
      ntrajR = 0
      nb = 1
      iskip = 100
      v0key = 0
      R0key = 1
      P0 = 0
      R0 = 0
      omega = 0
      modelname = 'Model'
      method = 'FSSH'

      vrkey = 0  !velocity reversal key
      nfrust_hop = 0  ! number of frustrated hop
      nfrust_hop2 = 0  ! number of frustrated hop
      ncpu = 1

      vckey = 1 ! velocity rescaling scheme: CA, 2: BA

      beta = 0.d0
      mp = 1.d0

      open(nread,file='param.in',status='old',IOSTAT=stat)
      if(stat .ne. 0) then
         write(0,*)  
         write(0,*) " ERROR!! "
         write(0,*) " 'param.in' file NOT FOUND! "
         write(0,*) " Please, provide necessary input file! "
         write(0,*)  
         stop
      end if

      call fdate(datentime)

      do while(loop)
        
         call getrec(safe,idnode,nread)
       
!c     convert to lowercase and strip out leading blanks
         call lowcase(record,lenrec)
         call strip(record,lenrec)
         call copystring(record,directive,lenrec)

         if(record(1).eq.'#'.or.record(1).eq.' ')then
!c     record is commented out
            cycle

         elseif(findstring('npart',directive,idum))then
!c     number of particle(s)
           np = intstr(directive,lenrec,idum)

         elseif(findstring('nstep',directive,idum))then
!c     number of simulation steps
           nsteps = intstr(directive,lenrec,idum)

         elseif(findstring('nsamp',directive,idum))then
!c     number of rpmd sample
           nsample = intstr(directive,lenrec,idum)

         elseif(findstring('ntraj',directive,idum))then
!c     number of of trajectories
           ntraj = intstr(directive,lenrec,idum)

         elseif(findstring('model',directive,idum))then
!c     selection of model
           if(findstring('tully1',directive,idum))then
             keymodel = 1
             modelname = 'TULLY MODEL 1'
           elseif(findstring('tully2',directive,idum))then
             keymodel = 2
             modelname = 'TULLY MODEL 2'
           elseif(findstring('tully3',directive,idum))then
             keymodel = 3
             modelname = 'TULLY MODEL 3'
           elseif(findstring('morse1',directive,idum))then
             keymodel = 4
             modelname = 'MORSE MODEL 1'
           elseif(findstring('morse2',directive,idum))then
             keymodel = 5
             modelname = 'MORSE MODEL 2'
           elseif(findstring('morse3',directive,idum))then
             keymodel = 6
             modelname = 'MORSE MODEL 3'
           elseif(findstring('db2lchain',directive,idum))then
             keymodel = 12
             modelname = '2-State with Linear Chain MODEL'
           elseif(findstring('db3lchain',directive,idum))then
             keymodel = 13
             modelname = '3-State with Linear Chain MODEL'
           endif

         elseif(findstring('nmode',directive,idum))then
!c     fast-fourier transform for nomral mode
           if(findstring('fft',directive,idum))then
            lfft = .true.
!c     matrix transform for nomral mode
           elseif(findstring('mat',directive,idum))then
            lfft = .false.
           endif

         elseif(findstring('nstate',directive,idum))then
!c     number of nstates
           nstates = intstr(directive,lenrec,idum)

         elseif(findstring('nbead',directive,idum))then
!c     number of nbeads
           nb = intstr(directive,lenrec,idum)
           if(nb .eq. 1)then
             method = 'FSSH'
           elseif(nb .gt. 1)then
             method='RPSH'
           else
!c     default
             nb = 1
             method = 'FSSH'
           endif

         elseif(findstring('ncpu',directive,idum))then
!c     number of cpus for parallel job
           ncpu = intstr(directive,lenrec,idum)

         elseif(findstring('timestep',directive,idum))then
!c     classical timestep (dt_c)
           dt = dblstr(directive,lenrec,idum)
           dtq=dt/10.d0

!         elseif(findstring('electime',directive,idum))then
!c     electronic timestep (dtq)
!           dtq = dblstr(directive,lenrec,idum)

         elseif(findstring('temp',directive,idum))then
!c     number of initial momentum
           kT = dblstr(directive,lenrec,idum)

         elseif(findstring('kinitial',directive,idum))then
!c     number of initial momentum
           P0 = dblstr(directive,lenrec,idum)

         elseif(findstring('beadpos',directive,idum))then
!c     position initialization for nbeads
           if(findstring('same',directive,idum))then
             R0key = 0
           elseif(findstring('gauss',directive,idum))then
             R0key = 1
           elseif(findstring('wigner',directive,idum))then
             R0key = 2
           else
             write(0,*) 'Default gaussian bead position initilazation!!'
             R0key = 1
           endif

         elseif(findstring('vinitial',directive,idum))then
!c     number of initial momentum
           if(findstring('fixed',directive,idum))then
             v0key = 0
           elseif(findstring('gauss',directive,idum))then
             v0key = 1
           elseif(findstring('wigner',directive,idum))then
             v0key = 2
           else
             write(0,*)'Invalid momentum initilazation,setting defacult'
             v0key = 1
             !!stop 
           endif

         elseif(findstring('vreverse',directive,idum))then
!c     velocity reversal scheme for frustrated hop
           if(findstring('never',directive,idum))then
             vrkey = 0
           elseif(findstring('alway',directive,idum))then
             vrkey = 1
           elseif(findstring('delv1',directive,idum))then
!c          Truhlar's delV scheme for frustrated hop
             vrkey = 2
           elseif(findstring('delv2',directive,idum))then
!c          Subotnik's delV2 scheme for frustrated hop
             vrkey = 3
           endif

         elseif(findstring('vrescale',directive,idum))then
!c     velocity  rescale to conserve energy
           if(findstring('cent',directive,idum))then
             vckey = 1
           elseif(findstring('bead',directive,idum))then
             vckey = 2
           endif

         elseif(findstring('usrkval',directive,idum))then
!c     number of initial momentum
           if(findstring('yes',directive,idum))then
             write(0,*) 'Enter k value:'
             read(*,*) P0
           endif

         elseif(findstring('acckval',directive,idum))then
!c     number of initial momentum
           if(findstring('yes',directive,idum))then
             lkval = .true.
           endif

         elseif(findstring('rundtail',directive,idum))then
!c     printing running detail
           if(findstring('yes',directive,idum))then
             ldtl = .true.
              dlevel=intstr(directive,lenrec,idum)
           endif

         elseif(findstring('iprint',directive,idum))then
!c     number of printing
           iskip = intstr(directive,lenrec,idum)

         elseif(findstring('finish',directive,idum))then
!c     safe termination of reading CONTROL file
            loop=.false.

         else
!c     unrecognised directive in control file
           kill=.true.
           if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
!           call error(idnode,-3)
         endif

      enddo
      close(nread)
! end reading input.in file
      end subroutine 

!**********************************************************************
      end module sysdef_module
