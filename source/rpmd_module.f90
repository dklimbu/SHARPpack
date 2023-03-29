      subroutine run_traj(Wb,rp,vp,KE)
!**********************************************************************
!     SHARP PACK routine to run RPMD trajectory     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module
      use nhc_module
      implicit none

      real*8, intent(in) :: Wb(np)
      real*8             :: rp(np,nb),vp(np,nb)
      real*8             :: fp(np,nb),p(np,nb)
      real*8             :: dt2,KE

      dt2 = dt*0.5d0

      if(llan)then  
!  Langevin pimd dynamics
         call Langevin(vp,KE)
      elseif(lnhc)then
!  NoseHooverChain pimd dynamics
         call nhc_part(vp,KE)
      endif

      call grhel(rp,fp)

      vp = vp + dt2*fp/mp

      p = vp*mp

! Evolve the ring polymer according to its ground state potential

      call freerp(np,p,rp)

      vp = p/mp

      call grhel(rp,fp)

      vp = vp + dt2*fp/mp

      if(llan)then  
!  Langevin pimd dynamics
         call Langevin(vp,KE)
      elseif(lnhc)then
!  NoseHooverChain pimd dynamics
         call nhc_part(vp,KE)
      endif

      end subroutine run_traj


      subroutine grhel(rp,fp)
!**********************************************************************
!     SHARP PACK subroutine to calculate force based on ground state 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module
      use modelvar_module, only: Wb,c
      implicit none

      integer             :: ip,ibd
      real*8              :: rp(np,nb)
      real*8, intent(out) :: fp(np,nb)
      real*8              :: hel(nstates,nstates)
      real*8              :: dhel(nstates,nstates,np)

      do ip=1,np
         do ibd=1,nb
            fp(ip,ibd) = -mp*(rp(ip,ibd)-R0)*Wb(ip)**2 !- c(ip,1) 
         enddo
      enddo
   
      !do ibd=1,nb
      !   CALL gethel(rp(:,ibd),hel(:,:),dhel(:,:,:))
      !   fp(1:np,ibd) = -dhel(1,1,:)
      !enddo

      return

      end subroutine grhel


      subroutine sample_init_vp(vp,rp)
!**********************************************************************
!     SHARP PACK routine to sample initial velocity and position 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module
      use modelvar_module, only : Wb,c
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: gaussn,sigR,sigP,alpha,sigRw,sigPw

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
      elseif((model>=7).and.(model<=9))then
         R0 = 0.0d0
      elseif(model==10)then  !PCET MODEL
         sigR = dsqrt(1.d0/(beta*f0))
         sigP = dsqrt(1.d0/(beta*mp))
         sigRw = sigR
         sigPw = sigP
         R0 = 0.0d0
         P0 = 0.d0
      elseif(model==11)then  !SuperExchange Model
         alpha = 0.25d0
         sigR = dsqrt(1.d0/(2.d0*alpha))
         sigP = dsqrt(alpha/2.d0)
         sigRw = dsqrt(1.d0/(4.d0*alpha))
         sigPw = dsqrt(hbar*hbar*alpha)
         R0 = -10.d0
      endif

      do ip=1,np
         ! if((model > 6).and.(Vinit == 1))then
         if((model .gt. 6).and.(model .le. 9))then
           if(Vinit == 1)then
              ! initialize Gaussian distributed bead momenta 
              sigP = dsqrt(mp*real(nb)/beta)
              sigR = sigP/(mp*Wb(ip))
           elseif(Vinit == 2)then
              ! initialize Wigner distributed bead momenta 
              sigPw = dsqrt(mp*hbar*Wb(ip)/(2.d0*tanh(0.5d0*beta/np*Wb(ip))))
              sigRw = sigPw/(mp*Wb(ip))
           else
              write(0,*) 'INCORRECT MOMENTUM INIT SCHEME SELECTED!!!'
              stop
           endif
         endif

         if((Rinit==0).and.(Vinit == 1))sigR = gaussn()*sigR + R0  
         if((Rinit==0).and.(Vinit == 2))sigR = gaussn()*sigRw + R0  
         do ibd=1,nb
            ! initialize bead positions
            if(Rinit==0)then
              ! assign each bead same position
              rp(ip,ibd) = sigR
              if(lpcet) rp(ip,ibd) = sigR + dsqrt(2.d0*lambda/f0)
            elseif(Rinit==1)then
              ! assign each bead gaussian distributed position
              rp(ip,ibd) = gaussn()*sigR +R0  !- c(ip,1)/Wb(ip)**2

              if(lpcet) rp(ip,ibd) = sigR + dsqrt(2.d0*lambda/f0)
            elseif(Rinit==2)then
              ! assign each bead wigner distributed position
              rp(ip,ibd) = gaussn()*sigRw +R0  !- c(ip,1)/Wb(ip)**2
              if(lpcet) rp(ip,ibd) = sigR + dsqrt(2.d0*lambda/f0)
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


      subroutine freerp (nf,p,q)
!**********************************************************************
!     SHARP PACK routine for 
!     -----------------------------------------------------------------
!     Free harmonic ring-polymer evolution through a time interval dt.
!     -----------------------------------------------------------------
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only: nb

      implicit none
      integer,parameter :: nbmax=128   !!1024
      integer      :: j,k,nf
      integer,save :: init = 0

      real*8       :: p(nf,nb),q(nf,nb)
      real*8       :: pjknew
      real*8,save  :: poly(4,nbmax)
!
      if (init .eq. 0) then
         if (nb .gt. nbmax) stop 'freerp 1'
         call ring(nf,poly)
         init = 1
      endif

      if (nb .eq. 1) then
         do j = 1,nf
            q(j,1) = q(j,1)+p(j,1)*poly(3,1)
         enddo
      else
         call realft(p,nf,nb,+1)
         call realft(q,nf,nb,+1)
         do k = 1,nb
            do j = 1,nf
               pjknew = p(j,k)*poly(1,k)+q(j,k)*poly(2,k)
               q(j,k) = p(j,k)*poly(3,k)+q(j,k)*poly(4,k)
               p(j,k) = pjknew
            enddo
         enddo
         call realft(p,nf,nb,-1)
         call realft(q,nf,nb,-1)
      endif

      return
      end subroutine

      subroutine ring(nf,poly)
!**********************************************************************
!     SHARP PACK routine to calculate 
!     -----------------------------------------------------------------
!     Monodromy matrix elements for free ring-polymer evolution.
!     -----------------------------------------------------------------
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only: nb,dt,mp,beta,hbar

      implicit none

      integer  :: k,nf
      real*8   :: poly(4,nb)
      real*8   :: betan,twown,pibyn
      real*8   :: wk,wt,wm,cwt,swt
!
      poly(1,1) = 1.d0
      poly(2,1) = 0.d0
      poly(3,1) = dt/mp
      poly(4,1) = 1.d0

      if (nb .gt. 1) then
         betan = beta/nb
         twown = 2.d0/(betan*hbar)
         pibyn = dacos(-1.d0)/nb

         do k = 1,nb/2
            wk = twown*dsin(k*pibyn)
            wt = wk*dt
            wm = wk*mp
            cwt = dcos(wt)
            swt = dsin(wt)
            poly(1,k+1) = cwt
            poly(2,k+1) = -wm*swt
            poly(3,k+1) = swt/wm
            poly(4,k+1) = cwt
         enddo

         do k = 1,(nb-1)/2
            poly(1,nb-k+1) = poly(1,k+1)
            poly(2,nb-k+1) = poly(2,k+1)
            poly(3,nb-k+1) = poly(3,k+1)
            poly(4,nb-k+1) = poly(4,k+1)
         enddo
      endif

      return
      end subroutine

      subroutine realft(data,m,n,mode)
!**********************************************************************
!     SHARP PACK routine to calculate 
!     -----------------------------------------------------------------
!     FFT of m real arrays (if mode = 1) or complex Hermitian
!     arrays in real storage (if mode = -1), using -lfftw3.
!     Works equally well with f77 and ifc.
!     -----------------------------------------------------------------
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
      integer,parameter    :: nmax=1024

      integer              :: m,n,mode
      integer              :: j,k,np
      integer*8            :: plana,planb

      real*8               :: scale
      real*8               :: copy(nmax)
      real*8,intent(inout) :: data(m,n)
      
      save copy,scale,plana,planb,np
!
      if (n .ne. np) then
         if (n .gt. nmax) stop 'realft 1'
         scale = dsqrt(1.d0/n)
         call dfftw_plan_r2r_1d(plana,n,copy,copy,0,64)
         call dfftw_plan_r2r_1d(planb,n,copy,copy,1,64)
         np = n
      endif

      do k = 1,m
         do j = 1,n
            copy(j) = data(k,j)
         enddo

         if (mode .eq. 1) then
            call dfftw_execute(plana)
         else if (mode .eq. -1) then
            call dfftw_execute(planb)
         else
            stop 'realft 2'
         endif

         do j = 1,n
            data(k,j) = scale*copy(j)
         enddo
      enddo

      return
      end


      subroutine get_energy(itraj,istep)
!**********************************************************************
!     SHARP PACK routine to calculate energy 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only : mp,np,nb
      use modelvar_module, only : rp,vp,Wb,rc,vc

      implicit none

      real*8            :: Ek, Ekc, Ering, Eringc

      integer  :: i, j, itraj,istep

      Ek = 0.0d0
      Ering = 0.0d0
      Ekc = 0.0d0
      Eringc = 0.0d0

!      Ek = 0.5 * mp * sum(vp*vp)
      Ekc = 0.5 * mp * sum(vc*vc)

      do i = 1, np
        do j = 1, nb
          Ek = Ek + 0.5 * mp * vp(i,j)*vp(i,j)
        enddo
      enddo


      do i = 1, np
        Ering = Ering + 0.5d0 * mp * Wb(i)**2 * (rp(i,1))**2
        Eringc = Eringc + 0.5d0 * mp * Wb(i)**2 * rc(i)**2
        do j = 2, np
          Ering = Ering + 0.5d0 * mp * Wb(i)**2 * (rp(i,j))**2
        enddo
      enddo

      write(100,'(2I8,30(e15.6E3))') itraj, istep, rc,vc, Ek, Ering, Ekc,Eringc,rp,vp

      end subroutine


      subroutine Langevin(vp,KE)
!**********************************************************************
!     SHARP PACK routine to calculate 
!     Langevin thermostat as implemented in Ceriotti2010
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only : mp,np,nb,beta,dt,pi,tau0
!      use modelvar_module, only : rp,vp,Wb,rc,vc

      implicit none

      integer :: ip, ib

      real*8,intent(inout)  :: vp(np,nb)
      real*8  :: p(np,nb)
      real*8  :: gama(nb), wk(nb), C1(nb), C2(nb), ranf(nb)
      real*8  :: betan, ph, dt2,KE
      real*8  :: gaussn

      p(:,:) = vp(:,:)*mp

      ph = pi/nb
      betan = beta/nb
      dt2 = dt/2.d0
!      tau0 = 0.7d0

      do ib = 1, nb
        wk(ib) = 2.d0/betan * sin((ib-1)*ph)
        gama(ib) = 2.d0 * wk(ib)
        if(ib==1) gama(ib) = 1.0/tau0
        C1(ib) = exp(-dt2*gama(ib))
        C2(ib) = dsqrt(1.d0 - C1(ib)*C1(ib))
!        write(111,*) ib, C1(ib),C2(ib)
      enddo

!  Langevin start here
!  normal mode transformation
      call realft(p,np,nb,+1)

      do ip = 1, np
         do ib = 1, nb
            ranf(ib) = gaussn()
         enddo

         p(ip,:) = C1*p(ip,:) + dsqrt(mp/betan)*C2 * ranf
      enddo

!  normal2bead transformation
      call realft(p,np,nb,-1)
      
      vp(:,:) = p(:,:)/mp

      KE = 0.5d0*sum(mp*vp*vp)

      end subroutine Langevin


      FUNCTION gaussn()
!**********************************************************************
!     SHARP PACK function to calculate normal distribution center at
!     origin with unit standard deviation 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
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
      use global_module
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: gaussn,sigR,sigP


      sigP = dsqrt(mp*real(nb)/beta)
      
      do ip=1,np

         do ibd=1,nb
           ! initialize bead positions
           rp(ip,ibd) = (ip-1)*10.d0/np

           ! initialize Gaussian distributed bead momenta 
           vp(ip,ibd) = gaussn()*sigP
         enddo
      enddo

      vp = vp/mp

      end subroutine sample_init_vp_lchain

