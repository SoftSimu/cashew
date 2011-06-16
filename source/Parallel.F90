!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! CASHEW                                                            <br /> !
!  Coarse Approach Simulator for Hydrogen-bonding effects in Water  <br /> !
!                                                                   <br /> !
! Molecular dynamics simulator for the 3D Mercedec-Benz model of    <br /> !
! water and atomic particles.                                       <br /> !
!                                                                   <br /> !
! Teemu Hynninen 2009                                               <br /> !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <br />
! parallel.f90 is a module containing routines needed for the
! parallel implementation of Cashew. The module also includes the
! "include mpif.h" command so that one need not include it in the main program
! files, "use mpi_mod" is enough.
! <br /><br />Back to <a href="cashew.html">cashew</a>
!
! *master_cpu=0 the master cpu index
! *max_shift=50 maximum shift of molecules from one cpu to another in load balancing
! *min_parts=16 minimum load for a single cpu
! *loadout=9913 output index for a file listing the workloads of all processors, for debug and efficieny testing only
! *loadwrite=.false. switch for controlling whether the workload file is written
! *cpu_id identification index of the cpu assigned by the MPI interface
! *start_mpi_load the starting index for the particle range this cpu is supposed to handle in force calculation
! *end_mpi_load the ending index for the particle range this cpu is supposed to handle in force calculation
! *first_mpi_mb the index of the first mb molecule this cpu is supposed to handle in force calculation
! *last_mpi_mb the index of the last mb molecule this cpu is supposed to handle in force calculation
! *first_mpi_atom the index of the first atom this cpu is supposed to handle in force calculation
! *last_mpi_atom the index of the last atom this cpu is supposed to handle in force calculation
! *start_mpi_load the starting index for the particle range this cpu is supposed to handle in evenly split parallellization
! *end_mpi_load the ending index for the particle range this cpu is supposed to handle in evenly split parallellization
! *first_mpi_mb the index of the first mb molecule this cpu is supposed to handle in evenly split parallellization
! *last_mpi_mb the index of the last mb molecule this cpu is supposed to handle in evenly split parallellization
! *first_mpi_atom the index of the first atom this cpu is supposed to handle in evenly split parallellization
! *last_mpi_atom the index of the last atom this cpu is supposed to handle in evenly split parallellization
! *stopwatch timer for workload estimate
! *my_load the cpu time workload recorded by this cpu
! *n_cpus the total numver of cpus
! *n_parts the total number of particles
! *n_ms the number of mbs
! *n_at the number of atoms
! *mpistat mpi status variable
! *mpistatus mpi status array
! *mpi_qtrn mpi buffer
! *mpi_mb mpi buffer
! *mpi_atom mpi buffer
! *loadsteps=1 counter of md steps for timing load balancing
! *loadcount counter of md steps for timing load balancing
! *mpi_m_force mpi buffer for mb forces
! *mpi_m_torque mpi buffer for mb torques
! *mpi_a_force mpi buffer for atom forces
! *mpi_bonds mpi buffer for bonds
! *mpi_n_bonds mpi buffer for number of bonds
! *loads buffer for workload records
! *limits buffer for the starting indices of particle ranges in the workload splitting between cpus
! *new_limits buffer for new starting indices of ranges given by workload balancing
! *sendbuf sending buffer
! *bondcount numbesr of molecules the cpus needs to handle in parallellized close neighbor finding (update_bond_numbers_p())
! *bondpos indices of the first mbs (-1) the cpus need to handle in parallellized close neighbor finding (update_bond_numbers_p())
! *mpi_m_n_nbor mpi buffer for mb numbers of neighbors
! *mpi_a_n_nbor mpi buffer for atom numbers of neighbors

MODULE mpi_mod
  
  USE parameters
  IMPLICIT NONE
#ifdef MPI
  include 'mpif.h'
#endif

  ! constants
  INTEGER, PARAMETER :: master_cpu = 0, &
       max_shift = 50, &
       min_parts = 16, &
       loadout = 9913
  LOGICAL, PARAMETER :: loadwrite = .false.

  ! processor-specific variables
  INTEGER :: cpu_id
  INTEGER :: start_mpi_load, end_mpi_load, first_mpi_mb, last_mpi_mb, first_mpi_atom, last_mpi_atom
  INTEGER :: start_mpi_split, end_mpi_split, &
       first_mpi_mb_split, last_mpi_mb_split, first_mpi_atom_split, last_mpi_atom_split
  INTEGER :: loadcount
  REAL(KIND=dp) :: stopwatch, my_load

  ! general variables
  INTEGER :: n_cpus, n_parts, n_ms, n_as, mpistat
#ifdef MPI
  INTEGER :: mpistatus(MPI_STATUS_SIZE)
#endif
  INTEGER :: mpi_qtrn, mpi_mb, mpi_atom
  INTEGER :: loadsteps = 1

  REAL(KIND=dp), ALLOCATABLE :: mpi_m_force(:,:), mpi_torque(:,:), mpi_a_force(:,:)
  ! arrays for temporary storage of very small force terms (for numeric accuracy)
  REAL(KIND=dp), POINTER :: small_m_force(:,:), small_torque(:,:), small_a_force(:,:)
  REAL(KIND=dp), POINTER :: tiny_m_force(:,:), tiny_torque(:,:), tiny_a_force(:,:)

  REAL(KIND=dp), ALLOCATABLE :: loads(:), mpi_bonds(:), mpi_n_bond(:)
  INTEGER, ALLOCATABLE :: limits(:), new_limits(:), sendbuf(:), bondcount(:), bondpos(:)
  INTEGER, POINTER :: mpi_m_n_nbor(:), mpi_a_n_nbor(:)

CONTAINS

  ! Routine for writing force calculation workload analysis data to a file
  SUBROUTINE write_loadmonitor()
    IMPLICIT NONE
    INTEGER :: ii

    IF(loadwrite)THEN
       IF(cpu_id == master_cpu)THEN
          DO ii = 0, n_cpus-1
             WRITE(loadout,'(I4,I6,ES9.2,A)',advance='no') ii, limits(ii), loads(ii), ", "
          END DO
          WRITE(loadout,'(A,ES9.2,F6.1,A,I2)',advance='yes') " max: ", MAXVAL(loads(:)), &
               n_cpus*100*MAXVAL(loads(:))/SUM(loads(:)), " % - next ", loadcount
       END IF
    END IF

  END SUBROUTINE write_loadmonitor

  ! Opens the output for writing workload data to a file called "mpi_load.out"
  SUBROUTINE open_loadmonitor()
    IMPLICIT NONE

    loadcount = 1
    IF(cpu_id == master_cpu)THEN
       IF(loadwrite)THEN
          OPEN(loadout,file="mpi_load.out")
       END IF
       ALLOCATE(loads(0:n_cpus-1))
       ALLOCATE(limits(-1:n_cpus-1))
       ALLOCATE(new_limits(-1:n_cpus-1))
       ALLOCATE(sendbuf(0:n_cpus-1))
       loads = 0.0
       limits = 0
       new_limits = 0
    END IF

  END SUBROUTINE open_loadmonitor

  ! Closes the output for wirting workload data
  SUBROUTINE close_loadmonitor()
    IMPLICIT NONE

    IF(loadwrite)THEN
       IF(cpu_id == master_cpu)THEN
          CLOSE(loadout)
       END IF
    END IF

  END SUBROUTINE close_loadmonitor  

#ifdef MPI
  ! Adjusts the load between cpus in the parallel force calculation.
  SUBROUTINE adjust_load()
    IMPLICIT NONE
    INTEGER :: shifter, ii, extra
    REAL :: efficiency, target_load, shifts(-1:n_cpus-2), biggest, loadsum
 
    ! The force calculation split is handled by giving each proceccor some ranges of
    ! indices for mbs and atoms. For the indices i in the mb range, the cpu calculates 
    ! the mb-mb forces for pairs (i,j) : i<j and the mb-atom forces for all pairs (i,*).
    ! For the indices k in the atom range, the cpu calculates the atom-atom forces for pairs
    ! (k,l) : k < l.
    !
    ! Visually, the particle pairs needed to be calculated can be arranged like this:
    !
    !       MBs  atoms
    ! atoms *****aaaa
    !       *****aaa
    !       *****aa
    !       *****a
    ! MBs   mmmmm         m = mb-mb pair
    !       mmmm          * = mb-atom pair
    !       mmm           a = atom-atom pair
    !       mm
    !       m
    !
    !       column
    !       123456789
    !       MBs  atoms
    !       123451234
    !
    ! And the workload can be thought of as being split to vertical columns in this diagram.
    ! In the workload splitting, the range of columns each cpu should handle in this scheme is first 
    ! determined. These column indices are then translated into index ranges for mbs and atoms.


    ! loadsteps counts the number of MD steps from the start of the run
    loadsteps = loadsteps + 1
    ! loadcount counts the number of steps since last adjustment
    ! don't adjust every turn: if loadcount is positive, just decrease it by one.
    ! adjust once loadcount reaches zero
    IF(loadcount > 0)THEN
       loadcount = loadcount - 1
       RETURN
    END IF

    ! adjust loads now!

    ! gather the load estimates from all processors to the master cpu
    CALL MPI_GATHER(my_load,1,MPI_DOUBLE_PRECISION,loads,1,MPI_DOUBLE_PRECISION,master_cpu,MPI_COMM_WORLD,mpistat)
    ! gather the end indices in the load split from all processors to the master cpu
    CALL MPI_GATHER(end_mpi_load,1,MPI_INTEGER,sendbuf,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)

    ! determine the new index ranges on the master cpu
    IF(cpu_id == master_cpu)THEN

       ! copy the gathered limits from the sendbuf
       ! limit has a bigger index range than the number of cpus so sendbuf had to be used above
       ! in the MPI communications
       limits(0:n_cpus-1) = sendbuf(0:n_cpus-1)       
       limits(-1) = 0 ! limit(ii) is the end of range for cpu ii, limit(ii-1)+1 is the start. cpu 0 starts from index 1.

       ! copy the old limist to array new_limits
       new_limits(-1:n_cpus-1) = limits(-1:n_cpus-1)
       
       ! calculate the average load per cpu from previous steps.
       ! if the load stays the same (it likely won't...), this is the load
       ! all the cpus should have
       target_load = SUM(loads(:))/n_cpus

       IF(loadsteps < 20 .OR. loadcount < 0)THEN ! due to the loadcount < 0 conditions, this if-clause is always performed
          
          ! estimate how much should limits be shifted
          loadsum = 0.0
          DO ii = 0, n_cpus-2 ! for n_cpus there are n_cpus+1 limits, of which the first and last are fixed. loop over the other n_cpus-1

             ! loadsum :
             ! ii = 0: loadsum = target-loads(0), so it directly tells how much is cpu 0 off balance
             !   this will be balanced by shifting the range between cpus 0 and 1 -> load of cpu 0 increases by target-loads(0)
             !   and load of cpu 1 decreases by this amount, mark balancer(0) = target-loads(0) = loadsum
             ! ii = 1: loadsum = 2*target-loads(0)-loads(1) = target+balancer(0)-loads(1)
             !   cpu 1 was off balance by target-loads(1), but the previous correction already changed this by +balancer(0)
             !   the rest of the off balance will be corrected by a shift between cpus 1 and 2 -> load of 1 incr. by
             !   target+balancer(0)-loads(1) = balancer(1) = loadsum, and that of 2 decr. by balancer(1)
             !     ...
             ! ii = n: loadsum = (n+1)*target - sum_i=0^n loads(i) = balancer(n-1)+target-loads(n)
             loadsum = loadsum + target_load-loads(ii) 

             ! estimate the shift in indices from the needed balancing in workload:
             ! loadshift_estimate is an estimate on how much the load changes if the index is shifted by one.
             ! of course, the load is not really linearly dependent on the index
             shifts(ii) = loadsum/&
                  loadshift_estimate(loads(ii)+loads(ii+1),limits(ii),limits(ii-1)+1,limits(ii+1))

          END DO

          ! if we are far from balance, some index shift estimates may be too big.
          ! if this is the case, scale down the shifts
          biggest = MAXVAL(ABS(shifts(:)))
          shifts(:) = shifts(:)/biggest*max_shift

          DO ii = 0, n_cpus-2 ! adjust the index limit between cpus ii & ii+1
             IF(shifts(ii) < 0)THEN  ! transfer load from ii to ii+1, move limit down
                shifter = -MIN(FLOOR(-shifts(ii)),new_limits(ii)-new_limits(ii-1)-1) ! make sure the limits don't cross
                new_limits(ii) = new_limits(ii) + shifter
             ELSE IF(shifts(ii) > 0)THEN ! transfer load from ii+1 to ii, move limit up
                ! here, the new limit is not allowed to pass the current value of the next limit
                ! even if the new value of the next limit is going to be larger than that.
                ! this ensures no limit crossing, but may slow down the balancing if
                ! the particle ranges are small (lots of cpus are used).
                shifter = MIN(FLOOR(shifts(ii)),new_limits(ii+1)-new_limits(ii)-1)
                new_limits(ii) = new_limits(ii) + shifter
             END IF
          END DO

       ELSE ! a simplified load shift scheme. this in not used currently as the one above is much better

          DO ii = 0, n_cpus-2 ! adjust limit between cpus ii & ii+1
             IF(loads(ii)-loads(ii+1) > 0)THEN  ! load from ii to ii+1, move limit down
                shifter = shift( (loads(ii)-loads(ii+1))/&
                     loadshift_estimate(loads(ii)+loads(ii+1),limits(ii),limits(ii-1)+1,limits(ii+1)) )
                shifter = -MIN(shifter,new_limits(ii)-new_limits(ii-1)-1)
                new_limits(ii) = new_limits(ii) + shifter
             ELSE IF(loads(ii)-loads(ii+1) < 0)THEN ! load from ii+1 to ii, move limit up
                shifter = shift((loads(ii+1)-loads(ii))/&
                     loadshift_estimate(loads(ii)+loads(ii+1),limits(ii),limits(ii-1)+1,limits(ii+1)) )
                shifter = MIN(shifter,new_limits(ii+1)-new_limits(ii)-1)
                new_limits(ii) = new_limits(ii) + shifter
             END IF
          END DO

       END IF

       ! copy the new limits to sendbuf for MPI communications
       sendbuf(0:n_cpus-1) = new_limits(0:n_cpus-1)
       
       ! decide on when to do the next balancing
       IF(loadsteps < 10)THEN ! we are still in the beginning of the run, so do the balancing often
          loadcount = 1 ! do two MD steps before another balancing
       ELSE ! the balancing has been done for a few times. estimate the efficiency and decide based on that

          ! efficiency estimate is the ratio of the max load and the average load
          !efficiency = MAXVAL(loads(:))/target_load
          efficiency = n_cpus*MAXVAL(loads(:))/SUM(loads(:))

          IF(efficiency < 1.05)THEN ! unbelievable efficiency
             loadcount = 100
          ELSE IF(efficiency < 1.10)THEN ! awesome efficiency
             loadcount = 50
          ELSE IF(efficiency < 1.20)THEN ! great efficiency
             loadcount = 25
          ELSE IF(efficiency < 1.50)THEN ! ok efficiency
             loadcount = 10
          ELSE IF(efficiency > 2.00)THEN ! bad efficiency
             loadcount = -1 ! rebalance again after next step
          ELSE ! not very good efficiency (1.5 - 2.0)
             loadcount = 1 ! rebalance in two steps
          END IF
       END IF

    END IF
    
    ! send the new end limits to other cpus
    CALL MPI_SCATTER(sendbuf,1,MPI_INTEGER,end_mpi_load,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    
    ! shift the sendbuf by one index
    IF(cpu_id == master_cpu) sendbuf(0:n_cpus-1) = new_limits(-1:n_cpus-2)

    ! send the new starting limits (-1) to other cpus (end of cpu i = start of cpu (i+1) -1)
    CALL MPI_SCATTER(sendbuf,1,MPI_INTEGER,start_mpi_load,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)

    ! correct the starting limits
    start_mpi_load = start_mpi_load+1

    ! send the loadcount to all cpus
    CALL MPI_BCAST(loadcount,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)

    ! clear data
    CALL wipe_load()

    ! write workload data if required
    IF(loadwrite) CALL write_loadmonitor()

    ! translate the index ranges to particle index ranges (each cpu does this by itself)
    CALL load_to_particles()

  END SUBROUTINE adjust_load
#endif

  ! Calculates an estimate on how much workload is shifted from a cpu to the previous if the
  ! limiting index shifts by one
  ! *theload current load of the two cpus
  ! *midindex the current index where the load between the two cpus is split
  ! *leftlim the starting index for the first cpu
  ! *rightlim the ending index for the second cpu
  ! *loadshift_estimate the estimated shift in load from cpu i to cpu i-1
  FUNCTION loadshift_estimate(theload,midindex,leftlim,rightlim)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: leftlim, rightlim, midindex
    REAL(KIND=dp) :: loadshift_estimate, theload

    ! The situation is such:
    ! 
    ! --11222--  ->  --11122--
    ! --11222-       --11122-
    ! --11222        --11122
    ! --1122         --1112
    ! --112          --111
    ! --11           --11
    ! --1            --1
    ! --             --
    ! -              -
    ! 
    ! So, it is estimated how many particle pairs are being transferred, and
    ! this is scaled to workload using the known current load

    loadshift_estimate = theload*(n_parts-midindex) / &
         ( (rightlim-leftlim+1) * (n_parts-(leftlim+rightlim)/2) )

  END FUNCTION loadshift_estimate

  ! Obsolete function.
  FUNCTION shift(effect)
    IMPLICIT NONE
    REAL(KIND=dp) :: effect
    INTEGER :: shift

    shift = FLOOR( MIN( sqrt(1.d0+2.d0*effect)-1.0, MAX(1.d0,REAL(max_shift,KIND=dp)) ) )

  END FUNCTION shift

  ! Adds the given load to my_load counter
  ! *addload the amount of load to be added
  SUBROUTINE record_load(addload)
    IMPLICIT NONE
    REAL(KIND=dp) :: addload

    my_load = my_load + addload

  END SUBROUTINE record_load

  ! Sets my_load counter to zero
  SUBROUTINE wipe_load()
    IMPLICIT NONE
    
    my_load = 0.0
    
  END SUBROUTINE wipe_load

  ! Records the wall clock time to stopwatch
  SUBROUTINE start_timer()
    IMPLICIT NONE

#ifdef MPI
    stopwatch = MPI_WTIME()
#else
    CALL cpu_time(stopwatch)
#endif

  END SUBROUTINE start_timer

  ! Reads the elapsed wall clock time since the previous starting of the timer
  ! and then restarts the timer.
  ! *timer the elapsed real time
  FUNCTION timer()
    IMPLICIT NONE
    REAL(KIND=dp) :: timer

    timer = stopwatch       ! previous time
#ifdef MPI
    stopwatch = MPI_WTIME() ! current time
#else
    CALL cpu_time(stopwatch)
#endif
    timer = stopwatch-timer ! time difference

  END FUNCTION timer

  ! Initializes the particle counters and arrays needed by the mpi machinery
  ! *mbs number of mbs
  ! *ats number atoms
  SUBROUTINE set_particles(mbs,ats)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mbs, ats
    INTEGER :: ii

    n_ms = mbs
    n_as = ats
    n_parts = n_ms+n_as
    IF(n_parts/n_cpus < min_parts .AND. n_cpus > 1)THEN ! too many cpu's; end this madness!
       IF(cpu_id == master_cpu)THEN
          WRITE(*,'(A,I8)') "Too many cpus for a system of this size, try ", &
               MAX( 1,floor(2.0**floor( log(REAL(n_parts)/REAL(min_parts))/log(2.0) )) )
       END IF
       CALL finish_mpi()
       STOP
    END IF
    ALLOCATE(mpi_m_force(3,n_ms))
    ALLOCATE(mpi_torque(3,n_ms))
    ALLOCATE(mpi_a_force(3,n_as))
    ALLOCATE(small_m_force(3,n_ms))
    ALLOCATE(small_torque(3,n_ms))
    ALLOCATE(small_a_force(3,n_as))
    ALLOCATE(tiny_m_force(3,n_ms))
    ALLOCATE(tiny_torque(3,n_ms))
    ALLOCATE(tiny_a_force(3,n_as))
    ALLOCATE(mpi_bonds( 1 : ((cpu_id+1)*n_ms)/n_cpus - (cpu_id*n_ms)/n_cpus ))
    ALLOCATE(bondcount(n_cpus))
    ALLOCATE(bondpos(n_cpus))
    
    DO ii = 0, n_cpus-1
       bondcount(ii+1) = ((ii+1)*n_ms)/n_cpus - (ii*n_ms)/n_cpus
       bondpos(ii+1) = (ii*n_ms)/n_cpus
    END DO

    ALLOCATE(mpi_m_n_nbor(n_ms))
    ALLOCATE(mpi_a_n_nbor(n_as))

  END SUBROUTINE set_particles

  ! Picks initial particle ranges for the cpus
  SUBROUTINE initial_loadlimits()
    IMPLICIT NONE
    INTEGER :: shift, leftshift

    IF(cpu_id == 0)THEN
       start_mpi_load = 1
       start_mpi_split = 1
    ELSE
       start_mpi_load = FLOOR(REAL(n_parts)*(1.0-sqrt(1.0-REAL(cpu_id)/REAL(n_cpus))))+1
       start_mpi_split = cpu_id*n_parts/n_cpus +1
    END IF
    IF(cpu_id == n_cpus)THEN
       end_mpi_load = n_parts 
       end_mpi_split = n_parts 
    ELSE
       end_mpi_load = FLOOR(REAL(n_parts)*(1.0-sqrt(1.0-REAL(cpu_id+1)/REAL(n_cpus))))
       end_mpi_split = (cpu_id+1)*n_parts/n_cpus
    END IF

#ifdef MPI
    ! if there are size 0 loads, shift the limits
    leftshift = 0
    shift = 0
    IF(cpu_id > 0)THEN
       CALL MPI_RECV(leftshift,1,MPI_INTEGER,cpu_id-1,1010,MPI_COMM_WORLD,mpistatus,mpistat)
       start_mpi_load = start_mpi_load + leftshift
    END IF
    IF(end_mpi_load < start_mpi_load)THEN
       shift = start_mpi_load - end_mpi_load
       end_mpi_load = end_mpi_load + shift       
    END IF
    IF(cpu_id < n_cpus-1)THEN
       CALL MPI_SEND(shift,1,MPI_INTEGER,cpu_id+1,1010,MPI_COMM_WORLD,mpistat)
    END IF
#endif

    CALL load_to_particles()

  END SUBROUTINE initial_loadlimits

  ! Translates the general particle index ranges to
  ! mb and atom index ranges for load splitting
  SUBROUTINE load_to_particles()
    IMPLICIT NONE

    first_mpi_mb = 0
    last_mpi_mb = 0
    first_mpi_atom = 0
    last_mpi_atom = 0
    first_mpi_mb_split = 0
    last_mpi_mb_split = 0
    first_mpi_atom_split = 0
    last_mpi_atom_split = 0

    IF(start_mpi_load <= n_ms)THEN
       first_mpi_mb = start_mpi_load  
    ELSE
       first_mpi_atom = start_mpi_load-n_ms
    END IF
    IF(start_mpi_split <= n_ms)THEN
       first_mpi_mb_split = start_mpi_split
    ELSE
       first_mpi_atom_split = start_mpi_split-n_ms
    END IF

    IF(end_mpi_load <= n_ms)THEN
       last_mpi_mb = end_mpi_load
    ELSE
       IF(first_mpi_mb > 0)THEN
          last_mpi_mb = n_ms
          first_mpi_atom = 1
       END IF
       last_mpi_atom = end_mpi_load-n_ms
    END IF
    IF(end_mpi_split <= n_ms)THEN
       last_mpi_mb_split = end_mpi_split
    ELSE
       IF(first_mpi_mb_split > 0)THEN
          last_mpi_mb_split = n_ms
          first_mpi_atom_split = 1
       END IF
       last_mpi_atom_split = end_mpi_split-n_ms
    END IF

  END SUBROUTINE load_to_particles

  ! Initialializes MPI structures
  SUBROUTINE initialize_mpi()
    IMPLICIT NONE

#ifdef MPI
    CALL MPI_INIT(mpistat)
    IF(mpistat /= MPI_SUCCESS)THEN
       WRITE(*,*) "MPI error"
       STOP
    END IF
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpus, mpistat)
    IF(mpistat /= MPI_SUCCESS)THEN
       WRITE(*,*) "MPI error"
       STOP
    END IF
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, cpu_id, mpistat)
    IF(mpistat /= MPI_SUCCESS)THEN
       WRITE(*,*) "MPI error"
       STOP
    END IF
#else
    n_cpus = 1
    cpu_id = 0
#endif

    CALL open_loadmonitor()

  END SUBROUTINE initialize_mpi

  ! Finishes MPI run
  SUBROUTINE finish_mpi()
    IMPLICIT NONE

    CALL close_loadmonitor()

#ifdef MPI
    CALL MPI_FINALIZE(mpistat)
    IF(mpistat /= MPI_SUCCESS)THEN
       WRITE(*,*) "MPI error"
       STOP
    END IF
#endif

  END SUBROUTINE finish_mpi

#ifdef MPI
  ! Broadcasts force information from the master cpu to helpers
  SUBROUTINE sync_forces(n_mols,n_ats,mb_forces,mb_torques,atom_forces,virial)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_mols, n_ats
    REAL(KIND=dp), POINTER :: mb_forces(:,:), mb_torques(:,:), atom_forces(:,:)
    REAL(KIND=dp), INTENT(INOUT) :: virial
    
    IF(n_mols > 0)THEN
       CALL MPI_BCAST(mb_forces,SIZE(mb_forces),MPI_DOUBLE_PRECISION,master_cpu,MPI_COMM_WORLD,mpistat)
       CALL MPI_BCAST(mb_torques,SIZE(mb_torques),MPI_DOUBLE_PRECISION,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF
    IF(n_ats > 0)THEN
       CALL MPI_BCAST(atom_forces,SIZE(atom_forces),MPI_DOUBLE_PRECISION,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF
    CALL MPI_BCAST(virial,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpistat)

  END SUBROUTINE sync_forces
#endif

!!$  ! Gets and sets the custom mpi types, not used currently
!!$  SUBROUTINE set_mpi_qtrn(qq)
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: qq
!!$
!!$    mpi_qtrn = qq
!!$
!!$  END SUBROUTINE set_mpi_qtrn
!!$
!!$  SUBROUTINE set_mpi_mb(qq)
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: qq
!!$
!!$    mpi_mb = qq
!!$
!!$  END SUBROUTINE set_mpi_mb
!!$
!!$  SUBROUTINE set_mpi_atom(qq)
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: qq
!!$
!!$    mpi_atom = qq
!!$
!!$  END SUBROUTINE set_mpi_atom

END MODULE mpi_mod
