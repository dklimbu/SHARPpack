      module initial_module
      use global_module
      implicit none

      contains

      subroutine sample_init_vp(vp,rp)
!**********************************************************************
!     SHARP PACK routine to sample initial velocity and position 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: sigR,sigP,alpha,sigRw,sigPw
      !real*8              :: gaussn

      if(model .le. 3)then
         alpha = 0.25d0
         sigR = dsqrt(1.d0/(2.d0*alpha))
         sigP = dsqrt(hbar*hbar*alpha/2.d0)
         sigRw = dsqrt(1.d0/(4.d0*alpha))
         sigPw = dsqrt(hbar*hbar*alpha)
         R0 = -15.d0
      elseif(model==4)then
         sigR = dsqrt(hbar/(2.d0*mp*omega))
         sigP = dsqrt(mp*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         R0 = 2.9d0
         Vinit = 1
      elseif(model==5)then
         sigR = dsqrt(hbar/(2.d0*mp*omega))
         sigP = dsqrt(mp*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         R0 = 3.3d0
         Vinit = 1
      elseif(model==6)then
         sigR = dsqrt(hbar/(2.d0*mp*omega))
         sigP = dsqrt(mp*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         R0 = 2.1d0
         Vinit = 1
      endif

      do ip=1,np
         
         do ibd=1,nb
            ! initialize bead positions
            if(Rinit==0)then
              ! assign each bead same position
              rp(ip,ibd) = sigR
            elseif(Rinit==1)then
              ! assign each bead gaussian distributed position
              rp(ip,ibd) = gaussn()*sigR +R0  !- c(ip,1)/Wb(ip)**2

            elseif(Rinit==2)then
              ! assign each bead wigner distributed position
              rp(ip,ibd) = gaussn()*sigRw +R0  !- c(ip,1)/Wb(ip)**2
            endif

            ! initialize bead momenta
            if(vinit == 0)then
               ! initialize deterministic bead momenta 
               vp(ip,ibd) = P0
            elseif(vinit == 1)then
               ! initialize Gaussian distributed bead momenta 
               vp(ip,ibd) = gaussn()*sigP + P0
            elseif(vinit == 2)then
               ! initialize Wigner distributed bead momenta 
               vp(ip,ibd) = gaussn()*sigPw + P0
            else
               write(0,*) 'INCORRECT MOMENTUM INIT SCHEME SELECTED!!!'
               stop
            endif
         end do

      end do

      vp = vp/mp

      return
      end subroutine sample_init_vp


      subroutine sample_init_vp_lchain(vp,rp)
!**********************************************************************
!     SHARP PACK routine to sample initial velocity and position for
!     N-Linear Chain Model
!
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: sigR,sigP
      !real*8              :: gaussn

      sigP = dsqrt(mp*real(nb)/beta)
      
      do ip=1,np

         do ibd=1,nb
           ! initialize bead positions
           rp(ip,ibd) = (ip-1)*2.d0/np
           !rp(ip,ibd) = gaussn()

           ! initialize Gaussian distributed bead momenta 
           vp(ip,ibd) = gaussn()*sigP
         enddo
      enddo

      vp = vp/mp

      end subroutine sample_init_vp_lchain
      

      FUNCTION gaussn()
!**********************************************************************
!     SHARP PACK function to calculate normal distribution center at
!     origin with unit standard deviation 
!     
!**********************************************************************
      IMPLICIT NONE

      REAL(8), PARAMETER  :: PI2=2.0*3.141592654
      REAL(8)             :: Z1,Z2,gaussn

      CALL RANDOM_NUMBER(z1)
      CALL RANDOM_NUMBER(z2)
      gaussn=dSQRT(-2.d0*dLOG(Z1))*dCOS(PI2*Z2)
!      gaussn=0.00193
      RETURN
      END FUNCTION gaussn


      subroutine modelParam(baseModel)
!**********************************************************************
!     SHARP PACK routine to set some specific  model parameters
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
      
      integer  :: ib,baseModel

!Mass, beta parameters for Tully's Models
      if(baseModel .le. 3)then  
        nstates = 2
        mp = 2000.0d0
        beta = mp/(P0*P0)   !1.0d0/P0            !mp/(P0**2)

! Parameters for Morse Models
      elseif(baseModel.gt.3 .and. model .le.6)then  
        nstates = 3
        mp = 20000.0d0
        beta = 1052.58d0   ! beta corresponding to 300K
        omega = 0.005d0

!Mass, beta parameters for 2-state linear chain Model
      elseif(baseModel .eq. 12)then  
        nstates = 2
        mp = 12.d0 * amu2au   !! amu --> a.u.
        beta = 1.d0/kT*temp    !! a.u.
        v11 = 0.d0 * kJ_mol2au  !! kJ/mol --> a.u.
        v22 = 8.0d0 * kJ_mol2au  !! kJ/mol --> a.u.
        d_ab = 0.d0
        d_ab(1,2)= -6.d0 * 0.52918d0  !-6.0 A^(-1)
        d_ab(2,1)=-d_ab(1,2)
        gamaLC(1) = 0.002418d0   !! 10^14 S^-1 --> a.u.
        sigmaLC(1) = sqrt(2.d0*gamaLC(1)*mp*nb/beta/dt) 
        do ib = 2, nb
          gamaLC(ib) = 2.d0*(2.d0*nb/(beta*hbar)*sin((ib-1)*pi/nb))  ! = 2w(k)
          sigmaLC(ib) = sqrt(2.d0*gamaLC(ib)*mp*nb/beta/dt)
        enddo 

!Mass, beta parameters for 3-state SuperExchange linear chain Model
      elseif(baseModel .eq. 13)then  
        nstates = 3
        mp = 12.d0 * amu2au   !! amu --> a.u.
        beta = 1.d0/kT*temp    !! a.u.
        v11 = 0.0   !!a.u.
        v22 = 0.01  !!a.u.
        v33 = 0.005 !!a.u.
        d_ab = 0.d0
        d_ab(1,2)= -6.5d0 !!a.u.
        d_ab(2,1)=-d_ab(1,2) 
        d_ab(2,3)= 8.d0  !!a.u.
        d_ab(3,2)=-d_ab(2,3)
        gamaLC(1) = 0.002418d0   !! 10^14 S^-1 --> a.u.
        sigmaLC(1) = sqrt(2.d0*gamaLC(1)*mp*nb/beta/dt) 
        do ib = 2, nb
          gamaLC(ib) = 2.d0*(2.d0*nb/(beta*hbar)*sin((ib-1)*pi/nb))  ! = 2w(k)
          sigmaLC(ib) = sqrt(2.d0*gamaLC(ib)*mp*nb/beta/dt)
        enddo

      endif

      
! set the point of excited state
      if(model .le. 3)then
         R0 = -15.d0
      elseif(model==4)then
         R0 = 2.9d0
      elseif(model==5)then
         R0 = 3.3d0
      elseif(model==6)then
         R0 = 2.1d0
      elseif(model==12)then  !LinearCahin MODEL
         R0 = 0.0d0
         P0 = 0.d0
      elseif(model==13)then  !LinearCahin MODEL
         R0 = 0.0d0
         P0 = 0.d0
      endif

      return
      end subroutine 

!**********************************************************************
    end module initial_module
