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
      
      np = 1
      dt = 1.0
      dtq = 0.d0

      model = 0
      nstates = 0
      nsteps = 0
      nsample = 0
      ntraj = 0
      ntrajR = 0
      nb = 1
      iskip = 100
      vinit = 0
      Rinit = 1
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

         elseif(findstring('nparticle',directive,idum))then
!c     number of timesteps
           np = intstr(directive,lenrec,idum)

         elseif(findstring('nsteps',directive,idum))then
!c     number of timesteps
           nsteps = intstr(directive,lenrec,idum)

         elseif(findstring('nsample',directive,idum))then
!c     number of timesteps
           nsample = intstr(directive,lenrec,idum)

         elseif(findstring('ntraj',directive,idum))then
!c     number of timesteps
           ntraj = intstr(directive,lenrec,idum)

         elseif(findstring('model',directive,idum))then
!c     selection of model
           if(findstring('tully1',directive,idum))then
              model = intstr(directive,lenrec,idum)
              modelname = 'TULLY MODEL 1'
           elseif(findstring('tully2',directive,idum))then
              model = intstr(directive,lenrec,idum)
              modelname = 'TULLY MODEL 2'
           elseif(findstring('tully3',directive,idum))then
              model = intstr(directive,lenrec,idum)
              modelname = 'TULLY MODEL 3'
           elseif(findstring('morse1',directive,idum))then
              model = intstr(directive,lenrec,idum)+3
              modelname = 'MORSE MODEL 1'
           elseif(findstring('morse2',directive,idum))then
              model = intstr(directive,lenrec,idum)+3
              modelname = 'MORSE MODEL 2'
           elseif(findstring('morse3',directive,idum))then
              model = intstr(directive,lenrec,idum)+3
              modelname = 'MORSE MODEL 3'
           elseif(findstring('db2lchain',directive,idum))then
              model = 12
              modelname = '2-State with Linear Chain MODEL'
           elseif(findstring('db3lchain',directive,idum))then
              model = 13
              modelname = '3-State with Linear Chain MODEL'
           endif

         elseif(findstring('nstates',directive,idum))then
!c     number of nstates
           nstates = intstr(directive,lenrec,idum)

         elseif(findstring('method',directive,idum))then
!c     method and number of nbeads
           if(findstring('fssh',directive,idum))then
              nb = 1
              method = 'FSSH'
           elseif(findstring('rpsh',directive,idum))then
              nb=intstr(directive,lenrec,idum)
              method='RPSH'
           endif

         elseif(findstring('ncpu',directive,idum))then
!c     number of cpus for parallel job
           ncpu = intstr(directive,lenrec,idum)

         elseif(findstring('timestep',directive,idum))then
!c     classical timestep (dt_c)
           dt = dblstr(directive,lenrec,idum)

         elseif(findstring('etimestep',directive,idum))then
!c     electronic timestep (dtq)
           dtq = dblstr(directive,lenrec,idum)

         elseif(findstring('temp',directive,idum))then
!c     number of initial momentum
           kT = dblstr(directive,lenrec,idum)

         elseif(findstring('kinitial',directive,idum))then
!c     number of initial momentum
           P0 = dblstr(directive,lenrec,idum)

         elseif(findstring('beadpos',directive,idum))then
!c     position initialization for nbeads
           if(findstring('same',directive,idum))then
             Rinit = 0
           elseif(findstring('gaussian',directive,idum))then
             Rinit = 1
           elseif(findstring('wigner',directive,idum))then
             Rinit = 2
           else
             write(0,*) 'Default gaussian bead position initilazation!!'
             Rinit = 1
           endif

         elseif(findstring('vinitial',directive,idum))then
!c     number of initial momentum
           if(findstring('fixed',directive,idum))then
             Vinit = 0
           elseif(findstring('gaussian',directive,idum))then
             Vinit = 1
           elseif(findstring('wigner',directive,idum))then
             Vinit = 2
           else
             write(0,*)'Invalid momentum initilazation,setting defacult'
             Vinit = 1
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
!c     velocity  scheme for frustrated hop
           if(findstring('ca',directive,idum))then
             vckey = 1
           elseif(findstring('ba',directive,idum))then
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
