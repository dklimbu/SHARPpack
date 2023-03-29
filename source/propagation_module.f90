      module propagation_module
      use global_module
      use models_module

      contains


      SUBROUTINE ADVANCE_ADIA(eva,vdotd_old,vdotd_new,adia)
!**********************************************************************
!     
!     SHARP PACK subroutine to advance the electronic coefficients
!     by 4th order Runge-Kutta (RK) method
!     
!**********************************************************************
      implicit none

      integer    :: i
      real*8     :: eva(NSTATES,2)
      real*8     :: vdotd_old(NSTATES,NSTATES),vdotd_new(NSTATES,NSTATES)
      real*8     :: eva_half(NSTATES)
      real*8     :: vdotd_half(NSTATES,NSTATES)
      complex*16 :: adia(NSTATES)
      complex*16 :: k1(NSTATES),k2(NSTATES),k3(NSTATES),k4(NSTATES)

      eva_half = 0.5*(eva(:,1)+eva(:,2))

      vdotd_half = 0.5*(vdotd_old+vdotd_new)
      
! 4th order RK scheme

      do i=1,nstates
        k1(i) = adia(i)*eva(i,1)/eye - sum(vdotd_old(i,:)*adia(:))
      end do

      do i=1,nstates
        k2(i) = (adia(i)+k1(i)*0.5*dtq)*eva_half(i)/eye-sum(vdotd_half(i,:)*(adia(:)+k1(:)*0.5*dtq))
      end do

      do i=1,nstates
        k3(i) = (adia(i)+k2(i)*0.5*dtq)*eva_half(i)/eye-sum(vdotd_half(i,:)*(adia(:)+k2(:)*0.5*dtq))
      end do

      do i=1,nstates
        k4(i) = (adia(i)+k3(i)*dtq)*eva(i,2)/eye-sum(vdotd_new(i,:)*(adia(:)+k3(:)*dtq))
      end do

      adia = adia+dtq/6.0*(k1+2.0*k2+2.0*k3+k4)

      return
      END SUBROUTINE ADVANCE_ADIA

!---------------------------------------------------------------
      SUBROUTINE ADVANCE_MD(istate,rp,vp,fp)
!**********************************************************************
!     
!     SHARP PACK subroutine to advance the classical positions
!     and velocities by velocity verlet
!     
!**********************************************************************

      implicit none

      integer              :: i,ip,ibd
      integer, intent(in)  :: istate
      real*8, intent(inout):: rp(np,nb),vp(np,nb),fp(np,nb)
      real*8               :: fran(np,nb)
      real*8               :: dhel(NSTATES,NSTATES,NP)
      real*8               :: hel(NSTATES,NSTATES)
      real*8               :: ee(NSTATES)
      real*8               :: p(np,nb)
      real*8               :: temppsi(NSTATES,NSTATES)
      real*8               :: gaussn

      real*8               :: aa(np,nb)

! for LinearChainModel
      real*8    :: fR(nb)

! RP propogation, velocity-Verlet

     !!ifLangevine dynamics for PCET 
      if(lpcet)then
        do ibd = 1, nb
          do ip = 1, np
            fran(ip,ibd) = sigma*gaussn()
          enddo
        enddo
        aa = (fp - f0*taul*vp + fran)/mp

      !!Langevin dynamic for just N-th particle in LinearChain Model
      elseif((model==12).or.(model==13))then

        do ibd = 1, nb
           fR(ibd) = sigmaLC * gaussn()     
        enddo

        aa = fp/mp
        aa(np,:) = (fp(np,:) - gamaLC*mp*vp(np,:) + fR)/mp

      else

        aa = fp/mp

      endif

     !!! Velocity-Verlet
      vp=vp+0.5*dt*aa     !! fp/mp

      p=vp*mp

      CALL freerp(np,p,rp)
      !rp = rp + vp*dt

      vp=p/mp

      do ibd=1,nb
        if((model==12).or.(model==13))then
        ! LinearChainModel classical force calculation
          CALL FORCE_Lchain(rp(1:np,ibd), fp(1:np,ibd))

        else
         CALL gethel(rp(:,ibd),hel(:,:),dhel(:,:,:))
         CALL Diag(ee(:),temppsi(:,:),hel(:,:))
         CALL FORCE(temppsi(:,:),istate, dhel(:,:,1:np), fp(1:np,ibd))

        endif
      enddo

! second half of RP propogation

      if(lpcet)then
        aa = (fp - f0*taul*vp + fran)/mp

      elseif((model==12).or.(model==13))then
        aa = fp/mp
        aa(np,:) = (fp(np,:) - gamaLC*mp*vp(np,:) + fR)/mp

      else
        aa = fp/mp
      endif

      vp=vp+0.5*dt*aa    !!fp/mp

      return
      END SUBROUTINE ADVANCE_MD
      
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FORCE(psi,istate,dhel,fp)
!**********************************************************************
!     
!     SHARP PACK subroutine to calculate force 
!     
!**********************************************************************

      use global_module, only : np, NSTATES
      implicit none

      integer             :: i,j,ip
      integer, intent(in) :: istate
      real*8, intent(in)  :: psi(NSTATES,NSTATES)
      real*8, intent(in)  :: dhel(NSTATES,NSTATES,NP)
      real*8, intent(out) :: fp(np)

      fp=0.d0
      
      do ip=1,np
        do i=1,nstates   
          do j=1,nstates
            fp(ip)=fp(ip)-psi(i,istate)*dhel(i,j,ip)*psi(j,istate)
          enddo
        end do
      end do

      END SUBROUTINE FORCE


!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FORCE_Lchain(r,ff)
!**********************************************************************
!     
!     SHARP PACK subroutine to calculate force of LinearChainModel
!     
!**********************************************************************

      use global_module, only : np
      implicit none

      integer             :: i,j,iL,iR
      real*8, intent(in)  :: r(np)
      real*8, intent(out) :: ff(np)
      real*8              :: rij,rr
      real*8              :: Vo,a

      Vo = 175.d0 * kJ_mol2au  !!0.00038125d0  !! kJ/mol --> a.u.
      a = 4.d0 * 0.5292d0 !! A^(-1) --> a.u^(-1)

      ff=0.d0

      do i=1,np
        iL = i-1
        iR = i+1
        if(i==1)iL=1
       ! if(i==np)iR=np

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
          endif
        enddo

      enddo

      END SUBROUTINE FORCE_Lchain

!*********** END MODULE ******************

      end module propagation_module
