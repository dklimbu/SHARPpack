      module print_module
!**********************************************************************
!     SHARP PACK module that contains printing variables     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module
      use modelvar_module
      implicit none

      contains

      subroutine printin()
!**********************************************************************
!     SHARP PACK subroutine to print system model parameters     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none

      integer, parameter :: nrite=2
      character(len=24)  :: datentime

      call fdate(datentime)

      open(nrite,file='param.out',status='unknown')

      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,*) '    *******      ***      ***        **         ********       ******        '
      write(nrite,*) '   **********    ***      ***       ****        **********     *********     '
      write(nrite,*) '  ***      ***   ***      ***      ******       ***     ***    ***     ***   '
      write(nrite,*) '  ***            ***      ***      *** ***      ***      ***   ***      ***  '
      write(nrite,*) '   *********     ************     ***   ***     ***    ****    ***     ***   '
      write(nrite,*) '     *********   ************     *********     *** *****      *********     '
      write(nrite,*) '            ***  ***      ***    ***********    *** ****       ********      '
      write(nrite,*) '  ***      ***   ***      ***    ***     ***    ***    ***     ***           '
      write(nrite,*) '   **********    ***      ***   ***       ***   ***     ***    ***           '
      write(nrite,*) '     *******     ***      ***  ****       ****  ***      ****  ***           '
      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,'(1x,A,2x,A)') '                              ',trim(method)                             
      write(nrite,*) '                    SIMULATION CONTROL PARAMETERS                            '
      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                 '
      if(ncpu == 1)then
       write(nrite,*)'   Serial Execution on 1 CPUs                    '
      elseif(ncpu > 1)then
       write(nrite,*)'   Parallel Execution: Running on',ncpu,' CPUs   '
      endif
      write(nrite,*) '                                                 '
      write(nrite,*) "JOB STARTED AT: ", datentime   
      write(nrite,*) '                                                 '

      write(nrite,'(/,1x,A,2x,A)')     'MODEL                               ', modelname

      write(nrite,'(/,1x,A,2x,I8)')    'number of states                    ', nstates

      write(nrite,'(/,1x,A,2x,I8)')    'number of particle(s)               ', np

      write(nrite,'(/,1x,A,2x,I8)')    'number of nbead(s)                  ', nb

      write(nrite,'(/,1x,A,2x,I8)')    'number of ntraj (in each CPU)       ', ntraj

      if(nsample .ne. 0)then
         write(nrite,'(/,1x,A,2x,I8)') 'number of nsample                   ', nsample
      endif

      write(nrite,'(/,1x,A,2x,I8)')    'number of nsteps                    ', nsteps

      write(nrite,'(/,1x,A,2x,f8.2,A)')'R0 at which wave-function excited   ', R0,' (a.u.)'

      write(nrite,'(/,1x,A,2x,f8.2,A)')'P0, initial momentum, k             ', P0,' (a.u.)'

      if(Rinit == 0)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'Same'"
      elseif(Rinit == 1)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'Gaussian'"
      elseif(Rinit == 2)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'Wigner distribution'"
      endif

      if(vinit == 0)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initializatoin scheme            ', "'Deterministic'"
      elseif(vinit == 1)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initialization scheme            ', "'Gaussian distribution'"
      elseif(vinit == 2)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initialization scheme            ', "'Wigner distribution'"
      endif

      if(vrkey == 0)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Never' "
      elseif(vrkey == 1)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Always' "
      elseif(vrkey == 2)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Truhlar' "
      elseif(vrkey == 3)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Subotnik' "
      endif

      if(vckey == 1)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Rescaling Scheme               'Centroid' "
      elseif(vckey == 2)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Rescaling Scheme               'Beads-Average' "
      endif

      if(lpcet)then
        write(nrite,'(/,1x,A,2x,E12.4)')  'mass of particle, mp (a.u.)         ', mp
      else
        write(nrite,'(/,1x,A,2x,f8.1)')  'mass of particle, mp (a.u.)         ', mp
      endif

      write(nrite,'(/,1x,A,2x,e13.6)') 'Inverse temperature of system, β    ', beta

      if(model >3 .and. model <=6)then
         write(nrite,'(/,1x,A,2x,e13.6)') 'Vibrational frequency, ω            ', omega
      endif

      write(nrite,'(/,1x,A,2x,f8.4)')  'Simulation time step, dt (a.u.)     ', dt

      write(nrite,'(/,1x,A,2x,f8.4)')  'Electronic time step, dtE (a.u.)    ', dtq

      write(nrite,'(/,1x,A,2x,I8)')    'print skip                          ', iskip

      if(model==7)then
         write(nrite,'(/,1x,A,2x,e13.6)') 'SPIN-BOSON PARAMATERS:'
         write(nrite,'(/,1x,A,2x,f8.4)')    '   Δ             ', delta
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ϵ             ', eps
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ξ             ', zxi
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωc            ', wc
         write(nrite,'(/,1x,A,2x,f8.4)')    '   β             ', beta
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH SPIN-BOSON PARAMATERS.'
      endif

      if(model==8)then
         write(nrite,'(/,1x,A,2x,e13.6)')   'SPIN-BOSON (DEBYE) PARAMATERS:'
         write(nrite,'(/,1x,A,2x,f8.4)')    '   Δ             ', delta/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ϵ             ', eps/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   E_r           ', E_r/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωc            ', wc/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωmax          ', wmax/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   KT(1/β)       ',1.d0/(beta*enu)
         write(nrite,'(/,1x,A,2x,f8.4,2x,e13.6)')    '   enu           ', enu*freq, enu
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH SPIN-BOSON PARAMATERS.'
      endif

      if(model==9)then
         write(nrite,'(/,1x,A,2x,e13.6)')   'FMO MODEL (DEBYE) PARAMATERS:'
         write(nrite,'(/,1x,A,2x,f8.4)')    '   E_r           ', E_r*freq
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωc            ', wc*freq
         write(nrite,'(/,1x,A,2x,e13.4)')   '   ωmax          ', wmax*freq
         write(nrite,'(/,1x,A,2x,f8.4)')    '   KT(1/β)       ', temp/beta
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH FMO PARAMATERS.'
      endif

      if(model==10)then
         write(nrite,'(/,1x,A,2x,e13.6)')   'PCET MODEL PARAMATERS:'
         write(nrite,'(/,1x,A,2x,I8)')      '   nbasis        ', nbasis
         write(nrite,'(/,1x,A,2x,f8.4)')    '   Δ             ', delta*energy
         write(nrite,'(/,1x,A,2x,f8.4)')    '   V_DA          ', Vda*energy
         write(nrite,'(/,1x,A,2x,f8.1)')    '   ωp            ', omega*freq
         write(nrite,'(/,1x,A,2x,f8.4)')    '   Temp(T)       ', temp/beta
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ε_o           ', eps0
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ε_∞           ', eps_inf
         write(nrite,'(/,1x,A,2x,f8.4)')    '   fo            ', f0
         write(nrite,'(/,1x,A,2x,f8.4)')    '   τ_o           ', tauo*tim
         write(nrite,'(/,1x,A,2x,f8.4)')    '   τ_d           ', taud*tim
         write(nrite,'(/,1x,A,2x,f8.4)')    '   τ_l           ', taul*tim
         write(nrite,'(/,1x,A,2x,f8.4)')    '   m             ', mp*tim*tim
         write(nrite,'(/,1x,A,2x,f8.4)')    '   λ             ', lambda*energy
         write(nrite,'(/,1x,A,2x,f8.4)')    '   qpA           ', qpA*dist
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH PCET PARAMATERS.'
      endif

      if(model==12)then
         write(nrite,'(/,1x,A,2x,e13.6)') '2-State LinearChain MODEL PARAMATERS:'
         write(nrite,'(/,1x,A,2x,e13.6)') '   T (in Kelvin)   ', kT
         write(nrite,'(/,1x,A,2x,e13.6)') '   V11 (kJ/mol)    ', v11/kJ_mol2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   V22 (kJ/mol)    ', v22/kJ_mol2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_12 (A^-1))    ', d_ab(1,2)/0.52918d0
         write(nrite,'(/,1x,A,2x,I4)')    '   N               ', np
         write(nrite,'(/,1x,A,2x,e13.6)') '   m (amu)         ', mp/amu2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   Vo (kJ/mol)     ', 175.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   a (A^-1)        ', 4.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   γ (s^-1)        ', 1.0e14
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH 2-State LinearChain PARAMATERS.'
      endif
      if(model==13)then
         write(nrite,'(/,1x,A,2x,e13.6)') '3-State LinearChain MODEL PARAMATERS:'
         write(nrite,'(/,1x,A,2x,e13.6)') '   T (in Kelvin)   ', kT
         write(nrite,'(/,1x,A,2x,e13.6)') '   V11 (a.u.)      ', v11
         write(nrite,'(/,1x,A,2x,e13.6)') '   V22 (a.u.)      ', v22
         write(nrite,'(/,1x,A,2x,e13.6)') '   V33 (a.u.)      ', v33
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_12 (a.u.)     ', d_ab(1,2)
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_23 (a.u.)     ', d_ab(2,3)
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_13 (a.u.)     ', d_ab(1,3)
         write(nrite,'(/,1x,A,2x,I4)')    '   N               ', np
         write(nrite,'(/,1x,A,2x,e13.6)') '   m (amu)         ', mp/amu2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   Vo (kJ/mol      ', 175.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   a (A^-1)        ', 4.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   γ (s^-1)        ', 1.0e14
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH 3-State LinearChain PARAMATERS.'
      endif

      write(nrite,*) '                                                 '
      write(nrite,*) '*************************************************************************'

      if(llan.and.(nsample>0))then
         write(nrite,"(/,1x,'Path Integral Langevin Equation (PILE)', &
           /,1x,'thermostat relaxation time',1p,e12.4)")tau0
         write(nrite,'(/,1x,A,2x,100e13.6)') 'Bath frequency, ωb[i]            ', wb
      endif

      if(lnhc.and.(nsample>0))then
         write(nrite,"(/,1x,'Nose-Hoover Chain', &
           /,1x,'thermostat relaxation time',1p,e12.4,&
           /,1x,'number of RESPA steps             ',1p,i6, &
           /,1x,'number of chains     ',1p,i6)") taut,nrespa,nchain
      endif

      close(nrite)

      return

      end subroutine


      subroutine printresults(lkval)
!**********************************************************************
!     SHARP PACK subroutine to print detail results of populations 
!     by different methods
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none
      
      integer :: i, j, itime, iiskip
      logical :: lkval,lfile

      iiskip = int(1.0/dt)

        open(56,file='adiabat1.out',status='unknown')
        open(66,file='adiabat2.out',status='unknown')
        open(76,file='diabat1.out',status='unknown')
        open(77,file='diabat2.out',status='unknown')
        open(78,file='diabat3.out',status='unknown')

        write(56,'(A13,2x,A)') '# MODEL:    ', modelname
        write(66,'(A13,2x,A)') '# MODEL:    ', modelname
        write(76,'(A13,2x,A)') '# MODEL:    ', modelname
        write(77,'(A13,2x,A)') '# MODEL:    ', modelname
        write(78,'(A13,2x,A)') '# MODEL:    ', modelname

        write(56,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC POPULATIONS: METHOD I'
        write(66,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC POPULATIONS: METHOD II'
        write(76,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD I'
        write(77,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD II'
        write(78,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD III'
        
       if(ntrajR > 0)then 
        open(86,file='adiabatR1.out',status='unknown')
        open(87,file='adiabatR2.out',status='unknown')
        open(88,file='diabatR3.out',status='unknown')

        write(86,'(A13,2x,A)') '# MODEL:    ', modelname
        write(87,'(A13,2x,A)') '# MODEL:    ', modelname
        write(88,'(A13,2x,A)') '# MODEL:    ', modelname
        
        write(86,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC REFLECTED POPULATIONS: METHOD I'
        write(87,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC REFLECTED POPULATIONS: METHOD II'
        write(88,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC REFLECTED POPULATIONS: METHOD III'
       endif

        if(dt.ge.1.0)then
        do itime=0,nprint
           write(56,222) itime*iskip*dt,(redmat(i,itime)/real(ntraj),i=1,nstates)
           write(66,222) itime*iskip*dt,((redmat_ec(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
           write(76,222) itime*iskip*dt,(diabat1(i,itime)/real(ntraj),i=1,nstates)
           write(77,222) itime*iskip*dt,(diabat2(i,itime)/real(ntraj),i=1,nstates)
           write(78,222) itime*iskip*dt,(diabat3(i,itime)/real(ntraj),i=1,nstates)
         
          if(ntrajR > 0)then 
           write(86,222) itime*iskip*dt,(redmatR(i,itime)/real(ntraj),i=1,nstates)
           write(87,222) itime*iskip*dt,((redmat_ecR(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
           write(88,222) itime*iskip*dt,(diabat3R(i,itime)/real(ntraj),i=1,nstates)
          endif
        end do

        else
        do itime=iiskip,nprint,iiskip
           write(56,222) itime*iskip*dt,(redmat(i,itime)/real(ntraj),i=1,nstates)
           write(66,222) itime*iskip*dt,((redmat_ec(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
           write(76,222) itime*iskip*dt,(diabat1(i,itime)/real(ntraj),i=1,nstates)
           write(77,222) itime*iskip*dt,(diabat2(i,itime)/real(ntraj),i=1,nstates)
           write(78,222) itime*iskip*dt,(diabat3(i,itime)/real(ntraj),i=1,nstates)

          if(ntrajR > 0)then 
           write(86,222) itime*iskip*dt,(redmatR(i,itime)/real(ntraj),i=1,nstates)
           write(87,222) itime*iskip*dt,((redmat_ecR(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
           write(88,222) itime*iskip*dt,(diabat3R(i,itime)/real(ntraj),i=1,nstates)
          endif
        end do
        endif 
        close(56)
        close(66)
        close(76)
        close(77)
        close(78)

        if(ntrajR > 0)then 
        close(86)
        close(87)
        close(88)
        endif

      if(lkval)then
        
        inquire(file='pkval.out', exist=lfile)
        if(lfile)then
           open(3,file='pkval.out',status='old',access='append')
        else
           open(3,file='pkval.out',status='new')
           write(3,'(A36,2x,A)') '# BRANCHING PROBABILTY OF MODEL:    ', modelname
           write(3,'(A)')        '# k              T1             T2             R1            R2            nR          nFrust      nFrust2'
        endif

        write(3,222) P0, ((redmat(i,nprint)-redmatr(i,nprint))/(ntraj), i=1,nstates),(redmatR(i,nprint)/(ntraj),i=1,nstates),&
                real(ntrajR)/ntraj,real(nfrust_hop)/ntraj,real(nfrust_hop2)/ntraj

        close(3)
      endif

222 format(120(e13.5E3,2x))
!222 format(f15.6,100(2x,f13.8))
      return

      end subroutine


      subroutine printlogo(nrite)
!**********************************************************************
!     SHARP PACK subroutine to print logo
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none

      integer            :: nrite
      !character(len=*)  :: ofile

      !open(nrite,file=ofile,status='unknown')

      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,*) '    *******      ***      ***        **         ********       ******        '
      write(nrite,*) '   **********    ***      ***       ****        **********     *********     '
      write(nrite,*) '  ***      ***   ***      ***      ******       ***     ***    ***     ***   '
      write(nrite,*) '  ***            ***      ***      *** ***      ***      ***   ***      ***  '
      write(nrite,*) '   *********     ************     ***   ***     ***    ****    ***     ***   '
      write(nrite,*) '     *********   ************     *********     *** *****      *********     '
      write(nrite,*) '            ***  ***      ***    ***********    *** ****       ********      '
      write(nrite,*) '  ***      ***   ***      ***    ***     ***    ***    ***     ***           '
      write(nrite,*) '   **********    ***      ***   ***       ***   ***     ***    ***           '
      write(nrite,*) '     *******     ***      ***  ****       ****  ***      ****  ***           '
      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '

      return 
      end subroutine


      subroutine openfile()
!**********************************************************************
!     SHARP PACK subroutine to print open files
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
  
      character(len=24) :: file_hopp, file_dcoup, file_therm

      file_hopp = 'hoppinghist.out'
      file_dcoup = 'dcoupling.out'
      file_therm = 'thermostate.out'

      open(nrite_hopp,file=file_hopp,status='unknown')
      call printlogo(nrite_hopp)
      write(nrite_hopp,'(A13,2x,A)') '# MODEL:    ', modelname
      
      if(ldtl)then
         open(nrite_dcoup,file=file_dcoup,status='unknown')
         write(nrite_dcoup,'(A13,2x,A)') '# MODEL:    ', modelname
         write(nrite_dcoup,'(1x,A)') '# nTraj      nSteps      R_1 (a.u.)    Velocity_1 (a.u.)     &
           KE    PE     Ering     TotalEnergy    Energy (nstates)           dCoupling (nstates,nstates)  &
           Diabatic Potentials Vij(nstates,nstates),         istate   Energy(istate)'
      endif

      if((lnhc.or.llan).and.(nsample>0))then
        open(nrite_therm,file=file_therm,status='unknown')
        write(nrite_therm,'(1x,A)') '# nTraj      nSteps      Temperature(K)   KEO(a.u)      PEO(a.u.)     Thermo_tot,  eta    peta'
      endif

      return

      end subroutine openfile


      subroutine hopping_stat()
!**********************************************************************
!     SHARP PACK subroutine for writing hopping statistics
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer :: i,j

      write(nrite_hopp,'("============= Hopping Statistics ==============")')
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
            write(nrite_hopp,'(" # Accepted jump from ",I3," to ",I3," : ",I8)') i,j,nJump(i,j)
            write(nrite_hopp,'(" # Rejected jump from ",I3," to ",I3," : ",I8)') i,j,nJumpFail(i,j)
          endif
        enddo
      enddo
      write(nrite_hopp,'(" #Total Successful Jump:",I8)') sum(nJump)
      write(nrite_hopp,'(" #Total Attempted Jump:",I8)') sum(nJump)+sum(nJumpFail)
      write(nrite_hopp,*)
      write(nrite_hopp,*) "#Total Initial states:", nIniStat
      write(nrite_hopp,*)
      write(nrite_hopp,*) "#Total Frustrated Hop:", nfrust_hop
      write(nrite_hopp,*) "#Reversed Frustrated Hop:", nfrust_hop2
      write(nrite_hopp,'("===============================================")')

      return

      end subroutine hopping_stat


      subroutine close_file()
!**********************************************************************
!     SHARP PACK subroutine to close files
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      close(nrite_hopp)
      if(ldtl)close(nrite_dcoup)
      close(nrite_therm)
      
      return 
      end subroutine close_file

!***********************************************************************
     subroutine printdetail()

     implicit none

     integer  :: ibd
     integer  :: itraj,itime


     end subroutine


!***********************************************************************

      end module print_module
