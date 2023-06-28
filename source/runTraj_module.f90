  module runTraj_module
!**********************************************************************
!     SHARP Pack module for running trajectories
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
  use global_module
  use modelvar_module
  use models_module
  use initial_module
  use rpmd_module
  use nonadiabatic_module
  use propagation_module
  use print_module

  implicit none

  contains

  subroutine runTraj()
!**********************************************************************
!     SHARP Pack routine to run trajectories
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer  :: i,j,k,l,itime,itraj,ibd,ip
  integer  :: istate, inext
  integer  :: nlit
  integer  :: iprint
  real*8   :: rnlit

  real*8   :: RANF
  real*8   :: crit,length,currentlength
  real*8   :: KE,PE,Ering,TotE,Vn

  integer,parameter  :: nrite_dtl1=1001
  integer,parameter  :: nrite_dtl2=1002

  ntrajR = 0
 
  nlit=int(dt/dtq)
  rnlit = 1.0/REAL(nlit,8)

  if(dlevel .eq. 1)then
     open(nrite_dtl1,file='HISTORYc',status='unknown')
     write(nrite_dtl1,'(1x,A)') '# nTraj       Time(fs)    rc     vc '
  endif
  if(dlevel .eq. 2)then
     open(nrite_dtl2,file='HISTORYb',status='unknown')
     write(nrite_dtl2,'(1x,A)') '# nTraj       Time(fs)    rp(np,nb)     vp(np,nb) '
  endif

! main trajectory loop
  DO itraj=1,ntraj
     
! randomly sample positon & velocity from gaussian distribution

    call sample_init_vp(vp_samp,rp_samp)

    if((model==12).or.(model==13))call sample_init_vp_lchain(vp_samp,rp_samp)

    rp = rp_samp
    vp = vp_samp

! Caculating the initial centroid variables

    rc=0.d0
    vc=0.d0

    do ibd=1,nb
      rc(:)=rc(:)+rp(:,ibd)/real(nb)
      vc(:)=vc(:)+vp(:,ibd)/real(nb)
    enddo

! Initialize position and momentum of solvent

! get the intial wavefunction 

    call gethel(rc(:),hel(:,:,1),dhel_rc(:,:,:))

    CALL Diag(eva(:,1),psi(:,:,1),hel(:,:,1))

    adia(:)=CMPLX(psi(1,:,1),0.)

    CALL RANDOM_NUMBER(RANF)

    crit=RANF

    currentlength=0.d0
    do i=1,nstates
      currentlength=currentlength+psi(1,i,1)**2
      if (crit.le.currentlength) then
        exit
      else
        continue
      endif
    enddo
    istate=i
    inext=i
    nIniStat(istate) = nIniStat(istate) + 1

    call compute_vdotd_old(vdotd_old,psi)

!=================================================================================
! Initialize the forces

    do ibd=1,nb
      CALL gethel(rp(:,ibd),hel(:,:,1),dhel(:,:,:))
      CALL Diag(eva(:,1),psi(:,:,1),hel(:,:,1))

      if((model==12).or.(model==13))then
        CALL FORCE_Lchain(rp(1:np,ibd), fp(1:np,ibd))

      else
        CALL FORCE(psi(:,:,1),istate, dhel(:,:,1:np), fp(1:np,ibd))

      endif
    enddo

!==================================================================================
! start a trajectory for n-steps     
    DO itime = 1, NSTEPS
        
! main nuclear propagation step
      CALL ADVANCE_MD(istate,rp,vp,fp)
       
! Caculating the centroid variables at each time step

      rc=0.d0
      vc=0.d0
   
      do ibd=1,nb
        rc(:)=rc(:)+rp(:,ibd) !/real(nb)
        vc(:)=vc(:)+vp(:,ibd) !/real(nb)
      enddo
      rc(:)=rc(:)/real(nb)
      vc(:)=vc(:)/real(nb)

      CALL gethel(rc(:),hel(:,:,2),dhel_rc(:,:,:))
      ! get the intial wavefunction 
      CALL Diag(eva(:,2),psi(:,:,2),hel(:,:,2))

      do i = 1, nstates
        psi(:,i,2) = psi(:,i,2) * dot_product(psi(:,i,1),psi(:,i,2))/abs(dot_product(psi(:,i,1),psi(:,i,2)))
      end do

      ! middle decoup based on SHS-Tully 94

      call compute_vdotd(vdotd_new,psi)

!      difference for ave_eva and vdotd
      deva(:) = eva(:,2) - eva(:,1)  
       
      vdotd = vdotd_new - vdotd_old
   
      ! electronic time steps
      DO j = 1, nlit
        eva(:,2) = eva(:,1) + deva * rnlit
        vdotd_new = vdotd_old + vdotd * rnlit 

        ! advance the adiabatic wave fucntion
        CALL ADVANCE_ADIA(eva(:,:),vdotd_old(:,:),vdotd_new(:,:),adia(:))

        CALL Hopping(istate,inext,vdotd_new(:,:),adia(:))         
           
        ! give the current to previous
        eva(:,1) = eva(:,2)      !eva(now) give to eva(old)
        vdotd_old = vdotd_new
      END DO

      psi(:,:,1) = psi(:,:,2)
            
      IF (istate.NE.inext) THEN
      ! total potential energy difference for the whole ring polymer hamiltonian 
 
       CALL EKINRESCALE(istate,inext,vc,vp,rp,eva(:,:),dhel_rc(:,:,:),psi(:,:,:))    
      END IF

!=============store the adiabatic/diabatic density matrix by differnet ways============
      if((mod(itime,iskip).eq.0).or.(itime==1))then
        iprint = int(itime/iskip)
        CALL  pop_estimator(istate,iprint,diabat1,diabat2,diabat3,redmat,redmat_ec)
      
! reflected probability calculation 
        if((vc(1)<0).and.(rc(1) < 0))then
          CALL  pop_estimator(istate,iprint,diabat1R,diabat2R,diabat3R,redmatR,redmat_ecR)
        endif
      endif

     ! printdetail()
     IF(ldtl)THEN
       ! rc=0.d0
       vc=0.d0
       do ibd=1,nb
       !  rc(:)=rc(:)+rp(:,ibd) !/real(nb)
         vc(:)=vc(:)+vp(:,ibd) !/real(nb)
       enddo
       !rc(:)=rc(:)/real(nb)
       vc(:)=vc(:)/real(nb)

       if(((ntraj .lt.10) .or. (mod(itraj,10).eq.0)).and.(mod(itime,iskip).eq.0))then
         call calEnergy(rp,vp,KE,Vn,Ering,TotE,istate)
         write(nrite_dcoup,101) itraj,itime*dt,rc(1),vc(1),KE,Vn,Ering,TotE,eva(:,1), &
                 vdotd_new,hel(:,:,2),real(istate),eva(istate,1)

         if(dlevel .eq. 1) write(nrite_dtl1,102) itraj,itime*dt,rc,vc
         if(dlevel .eq. 2) write(nrite_dtl2,102) itraj,itime*dt,(rp(ip,:),ip=1,np),(vp(ip,:),ip=1,np)

       endif
     ENDIF

! nuc step loop
    ENDDO
    if((vc(1)<0).and.(rc(1)<0)) ntrajR = ntrajR + 1

   ! traj loop
  ENDDO

 101 format(i10,1x,f15.3,1x,150(e18.8E3,2x))
 102 format(i10,1x,f15.3,1x,350(f15.6,2x))
   
  if(dlevel .eq. 1) close(nrite_dtl1)
  if(dlevel .eq. 2) close(nrite_dtl2)

  end subroutine


  subroutine compute_vdotd_old(vdotd_old,psi)
!**********************************************************************
!     SHARP Pack routine to compute vdotd_old
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: i,j,k,l,ip
  real*8,intent(out) :: vdotd_old(nstates,nstates)
  real*8,intent(in)  :: psi(nstates,nstates,2) ! old,current or new electronic step

  vdotd_old=0.d0

  if((model==12).or.(model==13))then
    do k=1,nstates
      do l=1,nstates
        if(k.ne.l)then
          do ip=1,1  !np
            vdotd_old(k,l) = vdotd_old(k,l)+vc(ip)*d_ab(k,l)
          enddo
        endif
      enddo
    enddo

  else
    do k=1,nstates
      do l=1,nstates
        if (k.ne.l) then
          !do i=1,nstates
          !  do j=1,nstates
          do ip=1,np
          !      vdotd_old(k,l) = vdotd_old(k,l)+psi(i,k,1)*dhel_rc(i,j,ip)*psi(j,l,1)*vc(ip)/(eva(l,1)-eva(k,1))
            vdotd_old(k,l) = vdotd_old(k,l)+sum(psi(:,k,1)*matmul(dhel_rc(:,:,ip),psi(:,l,1)))*vc(ip)/(eva(l,1)-eva(k,1))
          enddo
          !  enddo
          !end do
        endif
      end do
    end do
  endif

  end subroutine compute_vdotd_old


  subroutine compute_vdotd(vdotd_new,psi)
!**********************************************************************
!     SHARP Pack routine to compute vdotd_new
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: i,j,k,l,ip
  real*8,intent(out) :: vdotd_new(nstates,nstates)
  real*8,intent(in)  :: psi(nstates,nstates,2) ! old,current or new electronic step


  ! middle decoup based on SHS-Tully 94

  vdotd_new=0.d0

! CHECK using eq 32\times velocity Tully94 paper
  do k=1,nstates
    do l=1,nstates
      if(k.ne.l) then
        if((model==12).or.(model==13))then
        !!v*d_ij = p1*d_ij/m as for LinearChainModel
          vdotd_new(k,l) = vc(1)*d_ab(k,l)
        else
          vdotd_new(k,l) = sum(psi(:,k,1) * psi(:,l,2) - psi(:,k,2)*psi(:,l,1))/(2.*dt)
        endif
      end if
    enddo
  enddo

  end subroutine compute_vdotd
        

  subroutine calEnergy(rp,vp,KE,Vn,Ho,Hn,isurf)
!**********************************************************************
!     SHARP Pack routine to calculate energies at particular time
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
!  use global_module
  implicit none

  integer  :: ip,ibd, isurf

  real*8   :: rp(np,nb),vp(np,nb)
  real*8   :: omgn,omgn2
  real*8   :: eva_i,KE,Vn,Ho,Hn

  real*8   :: eva_b(nstates,2), psi_b(nstates,nstates,2)
  real*8   :: hel_b(nstates,nstates,2), dhel_b(nstates,nstates,np)

  omgn = nb/(beta*hbar)
  omgn2 = omgn*omgn

  KE = sum(vp*vp)*0.5d0*mp

!  PE = eva_i*nb

  ! calculate active energy of each bead
  
  Vn = 0.d0
  do ibd=1,nb

      CALL gethel(rp(:,ibd),hel_b(:,:,1),dhel_b(:,:,:))

      CALL Diag(eva_b(:,1),psi_b(:,:,1),hel_b(:,:,1))

      Vn = Vn + eva_b(isurf,1)
  enddo



  Ho = 0.d0
  do ip = 1, np
    do ibd = 1, nb
      if(ibd .eq. 1)then
        Ho = Ho + (rp(ip,ibd) - rp(ip,nb))**2
      else  
        Ho = Ho + (rp(ip,ibd) - rp(ip,ibd-1))**2
      endif
    enddo
  enddo

  Ho = 0.5d0*mp*omgn2*Ho

  Hn = KE+Vn+Ho

  end subroutine calEnergy
  
!**********************************************************************
  end module runTraj_module
