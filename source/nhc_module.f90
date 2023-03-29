      module nhc_module
!c**********************************************************************
!c     SHARP Pack module for defining Nose-Hoover Chain arrays for
!c     different ensembles
!c     
!c     authors    - D.K. Limbu & F.A. Shakib     
!c     copyright  - D.K. Limbu & F.A. Shakib
!c
!c     Method Development and Materials Simulation Laboratory
!c     New Jersey Institute of Technology
!c**********************************************************************

      implicit none

      real(8), allocatable, save :: eta_nhc(:)
      real(8), allocatable, save :: peta(:)
      real*8  :: thermo_tot 
      
      contains
      
      subroutine alloc_nhc_arrays(nchain)
!c**********************************************************************
!c     SHARP Pack routine to allocate arrays for Nose-Hoover Chain
!c     
!c     authors    - D.K. Limbu & F.A. Shakib     
!c     copyright  - D.K. Limbu & F.A. Shakib
!c
!c     Method Development and Materials Simulation Laboratory
!c**********************************************************************

      implicit none

      logical safe
      integer, intent(in) :: nchain
      integer, dimension(1:4) :: fail

      safe=.true.

      fail(:)=0
      
      allocate (eta_nhc(1:nchain),stat=fail(1))
      allocate (peta(1:nchain),stat=fail(2))
      
      if(any(fail.gt.0))safe=.false.
      
      end subroutine alloc_nhc_arrays
      
      subroutine dealloc_nhc_arrays()
      
!c**********************************************************************
!c     SHARP Pack routine to deallocate arrays for Nose-Hoover Chain
!c     
!c     authors    - D.K. Limbu & F.A. Shakib     
!c     copyright  - D.K. Limbu & F.A. Shakib
!c
!c     Method Development and Materials Simulation Laboratory
!c**********************************************************************
      
      implicit none

      logical safe
      integer, dimension(2) :: fail
      
      fail(:)=0
      safe=.true.
      
      deallocate(eta_nhc,peta,stat=fail(1))
      
      if(any(fail.gt.0))safe=.false.
        
      end subroutine dealloc_nhc_arrays


      subroutine nhc_init(nchain)
      
!c**********************************************************************
!c     SHARP Pack routine to initialise NVT-NHC thermostat
!c     
!c     authors    - D.K. Limbu & F.A. Shakib     
!c     copyright  - D.K. Limbu & F.A. Shakib
!c
!c     Method Development and Materials Simulation Laboratory
!c**********************************************************************
      
      implicit none
      
      integer i
      integer, intent(in) :: nchain

        do i=1,nchain
          eta_nhc(i)=0.d0
          peta(i)=0.d0
        enddo


      end subroutine nhc_init


      subroutine nhc_part(vp,KE)

!c*********************************************************************
!c     SHARP Pack routine to integrate and apply NHC thermostat
!c     together with NVT ensemble 
!c
!c     authors    - D.K. Limbu & F.A. Shakib     
!c     copyright  - D.K. Limbu & F.A. Shakib
!c
!c     Method Development and Materials Simulation Laboratory
!c*********************************************************************

      use global_module
      implicit none

      integer :: idnode,mxnode
      integer :: i,j,nw,respa,natm

      real*8  :: rp(np,nb), vp(np,nb)
      real*8  :: KE,scale_fac
      real*8  :: dt2,dt4,dt8,kbT,dNkT
      real*8  :: qmass(nchain)

      integer,save :: ic=0

!c     define block indices

      dt2=0.50d0*dt
      dt4=0.25d0*dt
      dt8=0.125d0*dt
!      kT=sigma_nhc
!      dNkT=sigma
      kbT= real(nb)/beta   ! 1.d0/beta gives NaN for higher nb >4 (1/26/2022)
      dNkT=1.d0*np*kbT

      qmass(1) = dNKT*taut*taut
      do i = 2, nchain
        qmass(i) = kbT*taut*taut
      enddo

      if(ic .eq. 0)then
!c     Assign Suzuki-Yoshida weights
         weight(1)=0.784513610477560d0
         weight(2)=0.235573213359357d0
         weight(3)=-1.17767998417887d0
         weight(4)=1.315186321d0
         weight(5)=-1.17767998417887d0
         weight(6)=0.235573213359357d0
         weight(7)=0.784513610477560d0

!         qmass_t = 2.d0*sigma*taut**2
!         qmass_part = 2.d0*sigma_nhc*taut**2

!         qmass_t = dNkT*taut*taut   !removed 2.0 as implement in Tuckerman
!         qmass_part = kbT*taut*taut

         ic=1
      endif

      call realft(vp,np,nb,+1)

!c     Start Suzuki-Yoshida scheme      
      do nw=nsy,1,-1

!c     Start RESPA loop
        do respa=1,nrespa

!c      Calculate kinetic energy
!          KE = mp*sum(vp*vp)   !removded 0.5 as in Tuckerman
          KE = mp*sum(vp(:,1)*vp(:,1))   !removded 0.5 as in Tuckerman

!c     Start 1st Suzuki-Yoshida scheme
          peta(nchain) = peta(nchain) + &
          (weight(nw)*dt4/nrespa)*((peta(nchain-1)**2/qmass(nchain-1))-kbT)
         
          do j=nchain-1,1,-1
            peta(j) = peta(j)* &
            exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass(j+1))

            if (j.eq.1) then
              peta(j) = peta(j) + (weight(nw)*dt4/nrespa)*(KE-dNkT)
            else
              peta(j) = peta(j) + &
              (weight(nw)*dt4/nrespa)*((peta(j-1)**2/qmass(j-1))-kbT)
            endif
            peta(j) = peta(j)* &
            exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass(j+1))
          enddo      
  
          do j=nchain,1,-1
            eta_nhc(j) = eta_nhc(j) - &
            (weight(nw)*dt2/nrespa)*peta(j)/qmass(j)
          enddo

!c     Thermostat the velocities
          scale_fac=exp((-weight(nw)*dt2/nrespa)*peta(1)/qmass(1))

!         do i=iatm0,iatm1
!           vxx(i)=scale*vxx(i)
!           vyy(i)=scale*vyy(i)
!           vzz(i)=scale*vzz(i)
!         enddo
          vp = scale_fac*vp

          KE=KE*scale_fac*scale_fac

!c     Start 2nd Suzuki-Yoshida scheme
          do j=1,nchain-1
            peta(j) = peta(j)* &
            exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass(j+1))
            if (j.eq.1) then
              peta(j) = peta(j) + (weight(nw)*dt4/nrespa)*(KE-dNkT)
            else
              peta(j) = peta(j) + &
              (weight(nw)*dt4/nrespa)*((peta(j-1)**2/qmass(j-1))-kbT)
            endif
            peta(j) = peta(j)* &
            exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass(j+1))
          enddo      

          peta(nchain) = peta(nchain) + &
          (weight(nw)*dt4/nrespa)*((peta(nchain-1)**2/qmass(nchain-1))-kbT)
      
        enddo
      enddo

      call realft(vp,np,nb,-1)

      thermo_tot = 0.5d0*peta(1)*peta(1)/qmass(1) + dNkT*eta_nhc(1)
      do i = 2, nchain 
        thermo_tot = thermo_tot + 0.5d0*peta(i)*peta(i)/qmass(i) &
                      + kbT*eta_nhc(i)
      enddo

      if(KE .ne. KE)then
        write(0,*) 'NaN detected in NHC THERMOSTAT !'
        write(0,100) 'peta:', peta
        write(0,100) 'eta:', eta_nhc
        write(0,100) 'KE:', KE
        stop
      endif
 100  format(A6,100e13.5)

      return

      end subroutine nhc_part

!c**********************************************************************
      end module nhc_module 
