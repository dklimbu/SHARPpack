      module models_module
!**********************************************************************
!     SHARP Pack module for defining system models
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
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
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!
!           | V11   V12   ... |
!      H  = | V21   V22   ... |
!           | ...   ...   ... |
!  
!     model = 1-3 (Tully Model 1,2,3)
!     model = 4-6 (Morse Model 1,2,3)     
!     model = 12 (2-state coupled with N-linear chain)     
!     model = 13 (2-state coupled with N-linear chain)     
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

! if(model == 12) <= detailed balance  2-State with N-Chain Model      
      elseif(model == 12)then

         hel(1,1) = v11
         hel(2,2) = v22

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

!**********************************************************************
      end module models_module
