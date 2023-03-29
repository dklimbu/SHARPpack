      module models_module
      use global_module
      use modelvar_module, only : Wb, c

      implicit none

      real*8,private :: A, B, Ct, D, Eo
      real*8,private :: Dm(3),bm(3),Re(3), cm(3)
      real*8,private :: Aij(3), Rij(3), amj(3)
      real*8,private :: aa(3),bb(3)  ! for 3-state super exchange model

      contains

      subroutine gethel(rxyz,hel, dhel)
!**********************************************************************
!     
!     SHARP PACK subroutine to define Hamiltonian of models
!     
!           | V11   V12   ... |
!      H  = | V21   V22   ... |
!           | ...   ...   ... |
!  
!     model = 1-3 (Tully Model 1,2,3)
!     model = 4-6 (Morse Model 1,2,3)     
!     model = 7-8 (Spin Boson Model 1,2)  
!     model = 9   (FMO Model) -> not working??
!     model = 10  (PCET Model) -> FSSH only
!     modle = 11  (Super Exchange Model) -> reflection coefficient??
!     
!**********************************************************************

      implicit none

      real*8 :: rxyz(np)
      real*8 :: hel(NSTATES,NSTATES)
      real*8 :: dhel(NSTATES,NSTATES,np)

      real*8,parameter :: m=1836.d0
      real*8  :: Ho(nstates,nstates)

      integer :: i, j, ip

      real*8  :: ff(np),Vq

      hel = 0.d0
      dhel = 0.d0

! if(model == 1) <= Tully Model I
      if(model == 1)then
         A = 0.01d0
         B = 1.6d0
         Ct = 0.005d0
         D = 1.0d0
         if(rxyz(1) > 0.0)then
             hel(1,1) = A*(1.0-exp(-B*rxyz(1)))
             dhel(1,1,1) = A*B*exp(-B*rxyz(1))
         else
             hel(1,1) = -A*(1.0-exp(B*rxyz(1)))
             dhel(1,1,1) = A*B*exp(B*rxyz(1))
         endif

         hel(2,2) = -1.*hel(1,1)
         hel(1,2) = Ct*exp(-D*rxyz(1)*rxyz(1))
         hel(2,1) = hel(1,2)

         dhel(2,2,1) = -1.*dhel(1,1,1)
         dhel(1,2,1) = hel(1,2)*(-2.*D*rxyz(1))
         dhel(2,1,1) = dhel(1,2,1)

! if(model == 2) <= Tully Model II
      elseif(model == 2)then
         A = 0.10d0
         B = 0.28d0
         Ct = 0.015d0
         D = 0.06d0
         Eo = 0.05d0

         hel(2,2) = -A*exp(-B*rxyz(1)*rxyz(1)) + Eo
         hel(1,2) = Ct*exp(-D*rxyz(1)*rxyz(1))
         hel(2,1) = hel(1,2)

         dhel(2,2,1) = (hel(2,2)-Eo) * (-2.*B*rxyz(1))
         dhel(1,2,1) = hel(1,2) * (-2.*D*rxyz(1))
         dhel(2,1,1) = dhel(1,2,1)

! if(model == 3) <= Tully Model III
      elseif(model == 3)then
         A = 0.0006d0
         B = 0.10d0
         Ct = 0.90d0
      ! interchange sign of A here for currect initialazation of wavepacket
      ! into state1

         hel(1,1) = -A  
         hel(2,2) = A

         if(rxyz(1) < 0.0)then
            hel(1,2) = B*exp(Ct*rxyz(1))
            dhel(1,2,1) = B*Ct*exp(Ct*rxyz(1))
         else
            hel(1,2) = B*(2.d0-exp(-Ct*rxyz(1)))
            dhel(1,2,1) = B*Ct*exp(-Ct*rxyz(1))
         endif

         hel(2,1) = hel(1,2)
         dhel(2,1,1) = dhel(1,2,1)

! if(model == 4) <= Morse Model I
      elseif(model == 4)then
         Dm = (/0.003d0, 0.004d0, 0.003d0/)
         bm = (/0.650d0, 0.600d0, 0.650d0/)
         Re = (/5.000d0, 4.000d0, 6.000d0/)
         cm = (/0.000d0, 0.010d0, 0.006d0/)

         Aij = (/0.002d0, 0.000d0, 0.002d0/)
         Rij = (/3.400d0, 0.000d0, 4.800d0/)
         amj = (/16.00d0, 0.000d0, 16.00d0/)

         !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i) * (1.0d0 - dexp(-bm(i)*(rxyz(1)-Re(i))))**2 +cm(i)

            dhel(i,i,1) = 2.d0*Dm(i)*bm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))*exp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(2,3) = Aij(3)*dexp(-amj(3)*(rxyz(1)-Rij(3))**2)

         hel(2,1) = hel(1,2)
         hel(3,2) = hel(2,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(2,3,1) = -2*amj(3)*(rxyz(1)-Rij(3)) * hel(2,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,2,1) = dhel(2,3,1)

! if(model == 5) <=  Morse Model II      
      elseif(model == 5)then
         Dm = (/0.020d0, 0.010d0, 0.003d0/)
         bm = (/0.650d0, 0.400d0, 0.650d0/)
         Re = (/4.500d0, 4.000d0, 4.400d0/)
         cm = (/0.000d0, 0.010d0, 0.020d0/)

         Aij = (/0.005d0, 0.005d0, 0.000d0/)
         Rij = (/3.660d0, 3.340d0, 0.000d0/)
         amj = (/32.00d0, 32.00d0, 0.000d0/)

        !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i) * (1.0d0 - dexp(-bm(i)*(rxyz(1)-Re(i))))**2 +cm(i)

            dhel(i,i,1) = 2*Dm(i)*bm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))*exp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(1,3) = Aij(2)*dexp(-amj(2)*(rxyz(1)-Rij(2))**2)

         hel(2,1) = hel(1,2)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(1,3,1) = -2*amj(2)*(rxyz(1)-Rij(2)) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(model == 6) <=  Morse Model III      
      elseif(model == 6)then
         Dm = (/0.020d0, 0.020d0, 0.003d0/)
         bm = (/0.400d0, 0.650d0, 0.650d0/)
         Re = (/4.000d0, 4.500d0, 6.000d0/)
         cm = (/0.020d0, 0.000d0, 0.020d0/)

         Aij = (/0.005d0, 0.005d0, 0.000d0/)
         Rij = (/3.400d0, 4.970d0, 0.000d0/)
         amj = (/32.00d0, 32.00d0, 0.000d0/)

         !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i) * (1.0d0 - dexp(-bm(i)*(rxyz(1)-Re(i))))**2 +cm(i)

            dhel(i,i,1) = 2*Dm(i)*bm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))*exp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(1,3) = Aij(2)*dexp(-amj(2)*(rxyz(1)-Rij(2))**2)

         hel(2,1) = hel(1,2)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(1,3,1) = -2*amj(2)*(rxyz(1)-Rij(2)) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(model == 7) <= Spin-Boson Model I (from Ananth2013, JCP 139,
! 124102 and modified to get result of Mandal2018)
      elseif(model == 7)then
         hel(1,1) = eps
         hel(2,2) = -eps
         dhel=0.d0

         do i=1,nstates
            do j=1,nstates
               if(i.ne.j)then
                  hel(i,j)=delta
               else
                  do ip=1,np
                     hel(i,j) = hel(i,j) + 0.5d0*mp*Wb(ip)**2*rxyz(ip)**2 + c(ip,i)*rxyz(ip)
                     dhel(i,j,ip) = mp * Wb(ip)**2 * rxyz(ip) + c(ip,i)
                  enddo
               endif
            enddo
         enddo

! if(model == 8) <= Spin-Boson Model-II-Debye (J. Chem. Phys. 144, 094104 (2016)
      elseif(model == 8)then
         hel(1,1) = eps!*.5d0
         hel(2,2) = -eps!*.5d0

         do i=1,nstates
            do j=1,nstates
               if(i.ne.j)then
                  hel(i,j)=delta
               else
                  do ip=1,np
                     hel(i,j) = hel(i,j) + 0.5d0*mp*Wb(ip)**2*rxyz(ip)**2 + c(ip,i)*rxyz(ip)
                     dhel(i,j,ip) = mp * Wb(ip)**2 * rxyz(ip) + c(ip,i)
                  enddo
               endif
            enddo
         enddo

! if(model == 9) <= 7-states FMO model)
      elseif(model == 9)then
         Ho = 0.0 
         Ho(:,1)=(/12410.0, -87.7,   5.5,   -5.9,    6.7,   -13.7,  -9.9/)
         Ho(:,2)=(/-87.7,  12530.0,  30.8,   8.2,    0.7,    11.8,   4.3/)
         Ho(:,3)=(/5.5,   30.8,   12210.0,  -53.5,  -2.2,   -9.6,    6.0/)
         Ho(:,4)=(/-5.9,  8.2,   -53.5,   12320.0,  -70.7,  -17.0, -63.3/)
         Ho(:,5)=(/6.7,   0.7,    -2.2,   -70.7,   12480.0,  81.1,  -1.3/)
         Ho(:,6)=(/-13.7, 11.8,   -9.6,   -17.0,   81.1,   12630.0,  39.7/)
         Ho(:,7)=(/-9.9,   4.3,    6.0,   -63.3,   -1.3,    39.7,   12440.0/)

         Ho = Ho / freq

!         hel(1,1) = eps*.5d0
!         hel(2,2) = -eps*.5d0

         do i=1,nstates
            hel(i,i) = Ho(i,i)
            do j=1,nstates
               if(i.ne.j)then
                  hel(i,j)=Ho(i,j)    !delta
               else
                  do ip=1,np
                     hel(i,j) = hel(i,j) + 0.5d0*mp*Wb(ip)**2*rxyz(ip)**2 + c(ip,1)*rxyz(ip)
                     dhel(i,j,ip) = mp * Wb(ip)**2 * rxyz(ip) + c(ip,1)
                  enddo
               endif
            enddo
         enddo

! if(model == 10) <= PCET model
      elseif(model == 10)then
         if(nbasis .eq. 1)then
            do i=1,nbasis
               hel(i,i)=0.5*m*omega**2*(rxyz(1))**2 + 0.5d0*f0*rxyz(1)**2
               hel(i+nbasis,i)=Vda
               hel(i,i+nbasis)=Vda
               hel(i+nbasis,i+nbasis)=0.5*m*omega**2*(rxyz(1)-qpA)**2+ 0.5*f0*(rxyz(1)-sqrt(2*lambda/f0))**2-delta
            enddo

            do i=1,nbasis
               dhel(i,i,1)=m*omega**2*rxyz(1) + f0*rxyz(1)
               dhel(i+nbasis,i+nbasis,1)=m*omega**2*(rxyz(1)-qpA) + f0*(rxyz(1)-sqrt(2*lambda/f0))
            enddo
         else
            do i=1,nbasis
               hel(i,i)=hbar*omega*((i-1)+0.5)+0.5d0*f0*rxyz(1)**2
               hel(i+nbasis,i)=Vda
               hel(i,i+nbasis)=Vda
               hel(i+nbasis,i+nbasis)=hbar*omega*((i-1)+0.5)+0.5*m*omega**2*qpA**2+0.5*f0*(rxyz(1)-sqrt(2*lambda/f0))**2-delta
               hel(i+nbasis-1,i+nbasis)=-dsqrt((m*(omega**3)*(qpA**2)*hbar*(i-1))/2)
               hel(i+nbasis+1,i+nbasis)=-dsqrt((m*(omega**3)*(qpA**2)*hbar*i)/2)
            enddo

            do i=1,nbasis
               dhel(i,i,1)=f0*rxyz(1)
               dhel(i+nbasis,i+nbasis,1)=f0*(rxyz(1)-sqrt(2*lambda/f0))
            enddo
         endif

! if(model == 11) <=  3-State Super Exchange Model      
      elseif(model == 11)then
         aa = (/0.d0,0.01d0,0.005d0/)
         bb = (/0.001d0,0.01d0,0.d0/)

         !Diagonal elements
         do i = 1,3
             hel(i,i) = aa(i)

             dhel(i,i,1) = 0.d0
         enddo

         !Off-Diagonal elements
         hel(1,2) = bb(1)*exp(-rxyz(1)*rxyz(1)*0.5d0)
         hel(2,3) = bb(2)*exp(-rxyz(1)*rxyz(1)*0.5d0)
         hel(1,3) = 0.d0

         hel(2,1) = hel(1,2)
         hel(3,2) = hel(2,3)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -rxyz(1) * hel(1,2)
         dhel(2,3,1) = -rxyz(1) * hel(2,3)
         dhel(1,3,1) = -rxyz(1) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,2,1) = dhel(2,3,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(model == 12) <= detailed balance  2-State with N-Chain Model      
      elseif(model == 12)then

         !call pot_Lchain(rxyz,ff,Vq)

         hel(1,1) = v11
         hel(2,2) = v22
                  
        !do ip=1,1
        !  do i=1,nstates
        !    dhel(i,i,ip) =0.d0 ! ff(ip)
        !  enddo
        !enddo

! if(model == 13) <= detailed balance 3-State SuperExchange with N-Chain Model      
      elseif(model == 13)then

         hel(1,1) = v11 
         hel(2,2) = v22 
         hel(3,3) = v33

      else
         write(0,*) 'Wrong Model System'
         write(0,*) 'Choose Proper Model System!!!'
         stop
      endif

      return
      end subroutine gethel
     

      SUBROUTINE DIAG(EVALUES,EVECT,CRV)
!**********************************************************************
!     
!     CRV: HERMITIAN MATRIX (INPUT)
!     EVECT: EIGENVECTORS (OUTPUT)
!     EVALUES: EIGENVALUES (OUTPUT)
! 
!**********************************************************************

      use global_module, only: nstates
      implicit none

      character :: JOBZ,UPLO
      integer   :: INFO
      integer   :: N,NAP,LDZ
      integer   :: I,J,IND
      real*8    :: AP(nstates*(nstates+1)/2),WORK(3*nstates)
      real*8    :: EVALUES(nstates)
      real*8    :: CRV(nstates,nstates),EVECT(nstates,nstates)
       
      N=nstates
      NAP=N*(N+1)/2
      LDZ=N

      EVALUES=0.
      EVECT=0.
      JOBZ='V' ! calculate both eigenvalue and eigenvector

      UPLO='L' ! lower diagonal matrix

      IND=0

      DO J=1,N
         DO I=J,N
            IND=IND+1
            AP(IND)=CRV(I,J)
         END DO
      END DO

      CALL DSPEV(JOBZ,UPLO,N,AP,EVALUES,EVECT,LDZ,WORK,INFO)

      RETURN
      
      END SUBROUTINE DIAG

!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE pot_Lchain(r,ff,pot)
!**********************************************************************
!     
!     SHARP PACK subroutine to calculate force 
!     
!**********************************************************************

      use global_module, only : np
      implicit none

      integer             :: i,j,iL,iR
      real*8, intent(in)  :: r(np)
      real*8, intent(out) :: ff(np)
      real*8              :: rij,rr
      real*8              :: Vo,a,pot

      Vo = 175.d0 * kJ_mol2au  !!0.00038125d0  !! kJ/mol --> a.u.
      a = 4.d0 * 0.5292d0 !! A^(-1) --> a.u^(-1)
      pot=0.d0
      ff=0.d0

      do i=1,1

      ! potential calculation
        if(i==np)then
           rij=r(i)-10.d0
           rr = sqrt(rij*rij)
        else
           rij = r(i)-r(i+1)
           rr = sqrt(rij*rij)
        endif
        pot=pot+Vo*(a*a*rr*rr-a*a*a*rr*rr*rr+0.58d0*a*a*a*a*rr*rr*rr*rr)

      ! force calculation
        iL = i-1
        iR = i+1
        if(i==1)iL=1

        do j=iL,iR
          if(i .ne. j)then
            if((i==np).and.(j.gt.np))then
              rij = r(i) - 10.d0
              rr = sqrt(rij*rij)
            else
              rij = r(i)-r(j)
              rr = sqrt(rij*rij)
            endif

            ff(i) = ff(i)-Vo*(2.d0*a*a*rr-3.d0*a*a*a*rr*rr+2.32d0*a*a*a*a*rr*rr*rr)*rij/rr
            !ff(i) =ff(i)-Vo*(2.d0*a*a*rij-3.d0*a*a*a*rij*rij+2.32d0*a*a*a*a*rij*rij*rij)
          endif
        enddo

      enddo

      END SUBROUTINE pot_Lchain


!**********************************************************************
!     
!    bath initialization for spin-boson/FMO models 
! 
!**********************************************************************
      subroutine bath_initial()
      use global_module
      use modelvar_module, only: Wb,c
      implicit none
      
      integer  :: i
! first time dynamics step, initialize everything.
      if(model==7)then
        do i=1,np
          Wb(i) = -wc * log(1.0 - real(i)*wmax/wc)      !tan((real(i)*atan(wmax/wc))/real(np))*wc
          c(i,1) = Wb(i) * dsqrt(zxi*wmax)
          c(i,2) = -c(i,1)
        enddo
      endif

! Bath Discretization: Debye Spectral for Spin-Boson Model
      if(model==8)then
        open(nrite_bath,file='bathfrequency.out',status='unknown')
        do i=1,np
          Wb(i) = tan((real(i)*atan(wmax/wc))/real(np))*wc
          c(i,1) = Wb(i)*dsqrt((E_r*atan(wmax/wc))/(pi*real(np)))
!         Wb(i) = wc*tan((i-0.5d0)*pi/(2.d0*np))
!         c(i,1) = Wb(i)*dsqrt(E_r/(2.d0*np))
          c(i,2) = -c(i,1)
          write(nrite_bath,'(i4,2e15.6)') i, wb(i)/enu,c(i,1)/enu
        enddo
        write(nrite_bath,'(A,3d15.6)') '#KT (in K), KT (amu), beta (amu):', (KT*enu*freq/0.695), (1.0/beta/enu), beta
        write(nrite_bath,'(A,2d15.6)') '#w_min, w_max: ', minval(Wb)/enu,maxval(Wb)/enu
        write(nrite_bath,'(A,2d15.6)') '#c_i_min, c_i_max: ', minval(c(:,1))/enu,maxval(c(:,1))/enu
        close(nrite_bath)
      endif

! Bath Discretization: Debye Spectral for FMO Model
      if(model==9)then
        open(120,file='bathfrequency.out',status='unknown')
        do i=1,np
          Wb(i) = tan((real(i)*atan(wmax/wc))/real(np))*wc
          c(i,1) = Wb(i)*dsqrt((E_r*atan(wmax/wc))/(pi*real(np)))
          write(120,'(i4,2e15.6)') i, wb(i)/enu,c(i,1)/enu
       enddo
       close(120)
       endif

!   if(model .le. 6)Wb = omega
      end subroutine bath_initial

      end module models_module
