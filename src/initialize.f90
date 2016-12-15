! initialize all include parameters, kernel and fields

module initialize
  
  implicit none
  private

  public :: init_all 

contains

!=============================================================
subroutine init_all(pb)
  
  use problem_class
  use mesh, only : init_mesh
  use constants, only : PI, MPI_parallel 
  use my_mpi, only: MY_RANK, NPROCS
  use fault_stress, only : init_kernel,nnLocalfft_perproc,nnoffset_perproc,& 
                           nnLocal_perproc,nnoffset_glob_perproc,&
                           nwLocal_perproc,nwoffset_glob_perproc 
  use output, only : ot_init, ox_init
  use friction, only : set_theta_star, friction_mu
  use utils, only : save_vectorV

!!$  use omp_lib

  type(problem_type), intent(inout) :: pb

  integer :: TID, NTHREADS, iproc, nwLocal, nLocal, nx, &
             nxfft, nwGlobal, nnGlobal
! Reading mesh partitions.
  call init_mesh(pb%mesh)

! If MPI parallel then gather the whole fault to compute the Kernel(iloc,iglob,nxfft)
if (MPI_parallel) then
    nwLocal = pb%mesh%nw !nw along dip taken by processor i.
    nx = pb%mesh%nx 
    nxfft = 2*pb%mesh%nx !For the kernel
    nLocal=nx*nwLocal !For fault nodes
    if (MY_RANK==0) write(6,*) 'NPROCS:',NPROCS
!  Assamble the array of nwLocal_perproc taken from each processor.
    allocate(nwLocal_perproc(0:NPROCS-1))
    allocate(nnLocalfft_perproc(0:NPROCS-1))
    allocate(nnLocal_perproc(0:NPROCS-1))
!   noffset_perproc: Allocation of the vectors sent by each processor. 
!   This is needed to make each vector in different memory place to avoid superposition.
    allocate(nnoffset_perproc(0:NPROCS-1))
    allocate(nnoffset_glob_perproc(0:NPROCS-1))
    allocate(nwoffset_glob_perproc(0:NPROCS-1))
    nwLocal_perproc=0
    nnLocalfft_perproc=0
    nnLocal_perproc=0
    nnoffset_perproc=0
    nnoffset_glob_perproc=0
    nwoffset_glob_perproc=0
! Assambled array of number of points per processor and send to all processors.
    write(6,*) 'nwlocal*nxfft:',nwLocal*nxfft
    write(6,*) 'nwlocal:',nwLocal
    call synchronize_all()
    call gather_alli(nwLocal,nwLocal_perproc,NPROCS)
    call gather_alli(nwLocal*nxfft,nnLocalfft_perproc,NPROCS)
    call gather_alli(nwLocal*nx,nnLocal_perproc,NPROCS)
    call synchronize_all()
    do iproc=0,NPROCS-1
      nnoffset_perproc(iproc)=sum(nnLocalfft_perproc(0:iproc))-nnLocalfft_perproc(iproc)
      nnoffset_glob_perproc(iproc)=sum(nnLocal_perproc(0:iproc))-nnLocal_perproc(iproc)
      nwoffset_glob_perproc(iproc)=sum(nwLocal_perproc(0:iproc))-nwLocal_perproc(iproc)
    enddo

   if (MY_RANK==0) then 
      write(6,*) 'nwLocal_perproc:',nwLocal_perproc  
      write(6,*) 'nnLocal_perproc:',nnLocal_perproc  
      write(6,*) 'nnLocalfft_perproc:',nnLocalfft_perproc  
      write(6,*) 'nnoffset_perproc:',nnoffset_perproc 
      write(6,*) 'nnoffset_glob_perproc:',nnoffset_glob_perproc 
      write(6,*) 'nwoffset_glob_perproc:',nwoffset_glob_perproc 
   endif

!    k%nnLocal=nnoffset_glob_perproc(MY_RANK)
    nwGlobal=sum(nwLocal_perproc)
    nnGlobal=nwGlobal*nx     
    pb%mesh%nnglob = nnGlobal !Needed later in the code
    allocate(pb%mesh%xglob(nnGlobal),pb%mesh%yglob(nnGlobal),&
             pb%mesh%zglob(nnGlobal),pb%mesh%dwglob(nwGlobal),pb%mesh%dipglob(nnGlobal))
! Adding global mesh for computing the kernel(ilocal,global,nxfft)
    call gather_allvdouble(pb%mesh%x,nLocal,pb%mesh%xglob,nnLocal_perproc, & 
                           nnoffset_glob_perproc,nnGlobal,NPROCS)
    call gather_allvdouble(pb%mesh%y,nLocal,pb%mesh%yglob,nnLocal_perproc, & 
                           nnoffset_glob_perproc,nnGlobal,NPROCS)
    call gather_allvdouble(pb%mesh%z,nLocal,pb%mesh%zglob,nnLocal_perproc, & 
                           nnoffset_glob_perproc,nnGlobal,NPROCS)
    call gather_allvdouble(pb%mesh%dip,nLocal,pb%mesh%dipglob,nnLocal_perproc, & 
                           nnoffset_glob_perproc,nnGlobal,NPROCS)
    call gather_allvdouble(pb%mesh%dw,nwLocal,pb%mesh%dwglob,nwLocal_perproc, & 
                           nwoffset_glob_perproc,nwGlobal,NPROCS)
!    call save_vectorV(pb%mesh%xglob,pb%mesh%yglob,pb%mesh%zglob,pb%mesh%zglob,&
!                      MY_RANK,'fault_xyz_global',nwGlobal,nx)
!    call save_vectorV(pb%mesh%x,pb%mesh%y,pb%mesh%z,pb%mesh%z,&
!                      MY_RANK,'fault_xyz_ilocal',nwLocal,nx)
endif
  
  write(6,*) 'Initializing parameters: ...'
!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_all;  derivs_all/derivs
   
  pb%v_pre = 0.d0
  pb%v_pre2 = 0.d0
  pb%pot_pre = 0.d0
  pb%ox%dyn_stat = 0
  pb%ox%dyn_stat2 = 0
  pb%ox%dyn_count = 0
  pb%ox%dyn_count2 = 0
    !---------------------- dt_max & perturbation------------------
  if (pb%Aper /= 0.d0 .and. pb%Tper > 0.d0) then
    if (pb%dt_max > 0.d0) then
      pb%dt_max = min(pb%dt_max,0.2d0*pb%Tper)
    else
      pb%dt_max = 0.2d0*pb%Tper
    endif
  endif
  if (pb%Tper > 0.d0) then
    pb%Omper = 2.d0*PI/pb%Tper
  else
    pb%Omper = 0.d0
  endif
  !---------------------- dt_max & perturbation------------------

  !---------------------- impedance ------------------
  if (pb%beta > 0d0) then 
    pb%zimpedance = 0.5d0*pb%smu/pb%beta
  else
    pb%zimpedance = 0.d0
  endif
  write(6,*) 'impedance = ', pb%zimpedance
  !---------------------- impedance ------------------
     
  !---------------------- ref_value ------------------        
  call set_theta_star(pb)
  pb%tau_init = pb%sigma * friction_mu(pb%v,pb%theta,pb) + pb%coh

  pb%tau = pb%tau_init
  pb%slip = 0d0 
  !---------------------- ref_value ----------------- 

  !---------------------- init_value for solver ----------------- 
  pb%time = 0.d0
  pb%itstop = -1
  pb%it = 0
  !---------------------- init_value for solver ----------------- 
  call init_kernel(pb%lam,pb%smu,pb%mesh,pb%kernel)

  write(6,*) "MY_RANK=",MY_RANK,',begin init_kernel'
  call ot_init(pb)
  call ox_init(pb)
  write(6,*) "MY_RANK=",MY_RANK,',finish init_kernel'

  if (MY_RANK==0) write(6,*) 'Initialization completed'
  call synchronize_all()

  ! Info about threads 
!!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!!$  TID = OMP_GET_THREAD_NUM()
!!$  write(6,*) 'Thread index = ', TID
!!$OMP BARRIER
!!$  if (TID == 0) then
!!$    NTHREADS = OMP_GET_NUM_THREADS()
!!$    write(6,*) 'Total number of threads = ', NTHREADS
!!$  end if
!!$OMP END PARALLEL

end subroutine init_all

end module initialize
