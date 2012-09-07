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
! cashew.f90 is the main program file for cashew.
! it contains the molecular dynamics operations that
! manipulate the positions and velocities of particles.
! for all other operations however, it calls routines in
! external modules.

! *system_name the name of the simulation
! *namelength length of the name
! *allostat used for checking allocation success
! *output index for the main output file
! *statout index for the statistics output file
! *supercell length of the simulation cell in x, y and z directions (this probably should have gone into "control_params")
! *periodic_boundary each component (x,y,z) is true if boundaries are periodic in the corresponding direction
! *boundary_type indices for boundary conditions (allows more than just periodic and free bounds)
! *boundary_value if the boundary conditions require a parameter, it is stored here
! *molecules list of MB molecules
! *atoms list of atomic particles
! *mb_dummy dummy molecule list used for storing the state in cg
! *at_dummy dummy atom list used for storing the state in cg
! *n_mols number of MB particles
! *n_ats number of atomic particles
! *n_elements number of different types of atomic particles
! *n_freedom number of degrees of freedom in the system
! *is_constrained true if constraints are being used
! *n_types list telling how many atoms of each type are there in the simulation (sum = n_ats)
! *mb_bonds list storing the number of MB-MB bonds seen by each MB molecule (for bond order calculations)
! *mbmb_neighbors list of MB-MB neighbors for each MB molecule
! *mbmb_n_nbors list of the number of MB-MB neighbors for each MB molecule
! *mbat_neighbors list of MB-atom neighbors for each MB molecule
! *mbat_n_nbors list of the number of MB-atom neighbors for each MB molecule
! *atat_neighbors list of atom-atom neighbors for each atom
! *atat_n_nbors list of the number of atom-atom neighbors for each atom
! *e_pot total potential energy
! *e_lenjon MB-MB lennard-jones energy
! *e_bonds MB-MB hyrogen bond energy
! *e_mbat MB-atom potential energy
! *e_atat atom-atom potential energy
! *e_constr potential energy due to constraining potentials/forces
! *e_kin kinetic energy
! *e_lin linear movement component of kinetic energy
! *e_rot rotational movement component of kinetic energy
! *last_e stores the previous total energy for cg
! *mb_forces forces acting on MB molecules
! *mb_torques torques acting on MB molecules
! *atom_forces forces acting on atoms
! *old_mb_forces list of forces acting on mb molecules at previous time, for cg
! *old_mb_torques list of torques acting on mb molecules at previous time, for cg
! *old_atom_forces list of forces acting on atoms at previous time, for cg
! *mb_conj_forces list of (transition) conjugate gradients for mb molecules
! *mb_conj_torques list of (rotation) conjugate gradients for mb molecules
! *atom_conj_forces list of (transition) conjugate gradients for atoms
! *n_bond bond count stored for bonding analysis
! *physical_params a custom type variable storing the parameters of the MB-model
! *control_params a custom type variable storing the control parameters for the simulation
! *simulation_time virtual time of the simulation
! *speaktimer determines the point in virtual time when the program should output information
! *xyztimer determines the virtual time for writing xyz file
! *bondtimer determines the virtual time for writing bonds file
! *rdftimer determines the virtual time for writing rdf file
! *adftimer determines the virtual time for writing adf file
! *inptimer determines the virtual time for wirting continuation input file
! *offset half the timstep, needed to account for numeric error in printing time evaluation
! *stattimer determines the virtual time for writing statistics
! *drift an upper limit estimate for how much any single particle has moved (for neighbors lists updating)
! *drift_range the value of drift for which the neighbors lists must be updated
! *timeA cpu time storage
! *timeB cpu time storage
! *time_io cpu time counter for input/output
! *time_force cpu time counter for force calculations
! *time_stat cpu time counter for statistics calculations
! *time_move cpu time counter for system updating
! *time_nbor cpu time counter for neighbor finding
! *time_tot total cpu time counter
! *cgbeta beta coefficient for conjugate gradients
! *barotimer determines the virtual time for applying barostat
! *barosteps the number of md steps between barostats, used for averaging the pressure
! *temps temperature
! *press (average) pressure
! *virial the inner virial, for calculating pressure
! *scaleratio the ratio of volumes for scaling (barostat)
! *finish marker for finishing cg loop
! *updated_stats true if statistics are up to date. marker for statistics output
! *rdf_ave the sum of measured rdfs
! *rdf_temp a temporary array for storing the rdf
! *rdf_sq the sum of squares of measured rdfs
! *adf_ave the sum of measured adfs
! *adf_temp a temporary array for storing the adf
! *adf_sq the sum of squares of measured adfs
! *rdf_count the number of times the rdf has been measured
! *adf_count the number of times the adf has been measured
! *tvector vector storing simulation starting time in year, month, day, UTC-correction in min, hour, min, sec, millisec - obtained using "date_and_time"
! *newt another time vector (for wall clock)
! *timediff vector in which the time difference between newt and tvector is written in day, hours, min, sec, millisec, control
! *timelabel string for date formatting
! *clocklabel string for time difference formatting
! *i index
! *j index
! *k index
! *boxes number of subcells in x, y and z directions
! *box_mbs list of MBs in subcell (x,y,z)
! *box_ats list of atoms in subcell (x,y,z)
! *box_mb_count number of MBs in subcell (x,y,z)
! *box_mb_count number of atoms in subcell (x,y,z)


PROGRAM cashew

  USE mb_model     ! Physics module
  USE file_handler ! I/O module
  USE functions    ! Functions module
  USE quaternions  ! Quaternion module
  USE parameters   ! Parameter module
  USE mpi_mod      ! Parallellization module
  USE mt95         ! Random number generator module
  IMPLICIT NONE

  ! i/o variables
  CHARACTER(LEN=100) :: system_name
  INTEGER :: namelength
  INTEGER :: allostat, output, statout

  ! system setup
  REAL(KIND=dp) :: supercell(3), boundary_value(3)
  LOGICAL :: periodic_boundary(3)
  INTEGER :: boundary_type(3)  
  TYPE(mb), POINTER :: molecules(:), mb_dummy(:)
  TYPE(atom), POINTER :: atoms(:), at_dummy(:)
  TYPE(gop), POINTER :: gos(:), go_dummy(:)
  INTEGER :: N_mols, N_ats, N_gos, N_elements, N_freedom
  LOGICAL :: is_constrained
  INTEGER, POINTER :: n_types(:)

  ! data storage
  REAL(KIND=dp), POINTER :: mb_bonds(:) ! bond numbers
  INTEGER, POINTER :: & ! neighbor lists
       mbmb_neighbors(:,:), mbmb_n_nbors(:), &
       mbat_neighbors(:,:), mbat_n_nbors(:), &
       mbgo_neighbors(:,:), mbgo_n_nbors(:), &
       atat_neighbors(:,:), atat_n_nbors(:)
  !LOGICAL, POINTER :: pair_done(:,:), atom_pair_done(:,:) ! to prevent double counting (obsolete)
  REAL(KIND=dp) :: e_pot,e_lenjon,e_bonds, e_mbat, e_atat, e_constr, e_kin, &
       e_lin, e_rot, last_e, & ! energy components
       e_mbgo, e_go,  e_stretch, e_bend, e_torsion, e_native, e_nonnat
  REAL(KIND=dp), POINTER :: mb_forces(:,:), mb_torques(:,:), atom_forces(:,:), & ! forces
       go_forces(:,:), old_go_forces(:,:), go_conj_forces(:,:), &
       old_mb_forces(:,:), old_mb_torques(:,:), old_atom_forces(:,:), & ! conjugate gradients
       mb_conj_forces(:,:), mb_conj_torques(:,:), atom_conj_forces(:,:), &
       n_bond(:) ! bond count

  ! special types storing input parameters (defined in parameters module)
  TYPE(mbps) :: physical_params
  TYPE(cps) :: control_params
  TYPE(gops) :: go_params

  ! virtual time, particle movement monitor
  REAL(KIND=dp) :: simulation_time, speaktimer, xyztimer, bondtimer, inptimer, rdftimer, offset, &
       adftimer, drift, drift_range, &
       timeA, timeB, time_io, time_force, time_stat, time_nbor, time_move, time_tot, &
       cgbeta, barotimer
  REAL(KIND=dp), ALLOCATABLE :: stattimer(:)
  INTEGER :: barosteps
  LOGICAL :: finish, updated_stats
  ! temperature, pressure
  REAL(KIND=dp) :: temps, press, virial, scaleratio
  ! rdf averaging
  REAL(KIND=dp) :: rdf_ave(rdfbins,4), rdf_sq(rdfbins,4), rdf_temp(rdfbins,4), &
       adf_ave(rdfbins,rdfbins,3), adf_sq(rdfbins,rdfbins,3), adf_temp(rdfbins,rdfbins,3)
  INTEGER :: rdf_count, adf_count

  ! time
  INTEGER :: tvector(8), newt(8), timediff(6)
  CHARACTER(LEN=10) :: timelabel
  CHARACTER(LEN=18) :: clocklabel

  ! box decomposition
  INTEGER :: boxes(3,2)
  INTEGER, POINTER :: box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), box_gos(:,:,:,:,:), &
       box_mb_count(:,:,:,:), box_at_count(:,:,:,:), box_go_count(:,:,:,:)
  LOGICAL :: boxed
  
  ! general
  INTEGER :: i, j, k

  REAL(KIND=dp) :: shifty, nep(3), nem(3)


  !**********************
  !
  ! Start Initialization
  !
  !**********************

  ! mpi
  CALL initialize_mpi()

  ! cpu time counters
  timeB = 0.d0
  timeA = 0.d0
  time_tot = 0.d0
  time_io = 0.d0
  time_force = 0.d0
  time_stat = 0.d0
  time_nbor = 0.d0
  time_move = 0.d0
  ! record the first cpu time checkpoint
  CALL checkpoint(time_tot,timeA,timeB)

  ! Read system name as a command line parameter
  IF( command_argument_count()==0)THEN
     WRITE(*,*) "Give the system name as command-line parameter."
     WRITE(*,*) "The input-files should be named 'system_name."//MB_IN//"'."
     STOP
  ELSE
     CALL get_command_argument(1,system_name)
     IF(INDEX(system_name,'.mb') == 0)THEN
         namelength = INDEX(system_name,' ')
     ELSE
         namelength = INDEX(system_name,'.mb')
     END IF
  END IF  


  !
  ! Extract a seed for the random number generator from time
  ! tvector holds: year, month, day, ?, h, min, s, ms
  !
  CALL DATE_AND_TIME(values=tvector)
  control_params%rnd_seed = SUM(tvector)+49*tvector(7)+537*tvector(8)  
  
  ! make sure all cpus have the same random number sequence
#ifdef MPI
  CALL MPI_BCAST(control_params%rnd_seed,1,MPI_DOUBLE_PRECISION,master_cpu,MPI_COMM_WORLD,mpistat)
#endif

  ! the RNG is initialized in input parsing, since if a seed is given
  ! in the input file, it overrides the one obtained from the clock

  !!!!!!!!!!!!!!!!!!!!!!!!!
  !  Read the input file  !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! We don't know yet if the user wants output (verbose hasn't been specified yet),
  ! so it's best to not write out yet...
  ! WRITE(*,*) "Reading input file "//system_name(1:namelength-1)//"."//MB_IN
  CALL parse_input(system_name(1:namelength-1)//"."//MB_IN,&
       control_params,physical_params,go_params,&
       supercell,periodic_boundary,boundary_type,boundary_value,&
       molecules,atoms,gos,&
       n_mols,n_ats,n_gos,n_elements,n_freedom,n_types,is_constrained)
  ! Notify of read input
  IF(control_params%verbose > 0.d0 .AND. cpu_id == master_cpu)THEN
     WRITE(*,*) "Read the input file "//system_name(1:namelength-1)//"."//MB_IN
  END IF
  ! cpu time checkpoint
  CALL checkpoint(time_io,timeA,timeB)

  ! Initialize data arrays
  IF(control_params%verbose > 0.d0 .AND. cpu_id == master_cpu) &
       WRITE(*,*) "Initializing data structures"  
  CALL init_data_arrays(molecules,atoms,gos,supercell,periodic_boundary,&
       mb_bonds,boxed,boxes,&
       box_mbs,box_ats,box_gos,&
       box_mb_count,box_at_count,box_go_count,&
       mbmb_neighbors,mbmb_n_nbors,&
       mbat_neighbors,mbat_n_nbors,&
       mbgo_neighbors,mbgo_n_nbors,&
       atat_neighbors,atat_n_nbors,&
       physical_params,control_params,go_params)  

  ! Allocate force and torque arrays
  ALLOCATE(mb_forces(3,n_mols), STAT = allostat)
  IF(allostat /= 0) CALL abort("allocating force array")
  ALLOCATE(mb_torques(3,n_mols), STAT = allostat)
  IF(allostat /= 0) CALL abort("allocating torque array")
  ALLOCATE(atom_forces(3,n_ats), STAT = allostat)
  IF(allostat /= 0) CALL abort("allocating force array")
  ALLOCATE(go_forces(3,n_gos), STAT = allostat)
  IF(allostat /= 0) CALL abort("allocating force array")



  ! if bond monitoring is requested, allocate arrays for that
  IF(control_params%bond_writer /= noxyz_index)THEN
     ALLOCATE(n_bond(n_mols))
     ALLOCATE(mpi_n_bond(n_mols))
  END IF

  ! allocate arrays for conjugate gradients (ignored for MD)
  IF(control_params%md_algo == cg_index)THEN
     ALLOCATE(mb_dummy(n_mols), STAT=allostat)
     IF(allostat /= 0) CALL abort("allocating particle array")
     ALLOCATE(at_dummy(n_ats), STAT=allostat)
     IF(allostat /= 0) CALL abort("allocating particle array")

     ALLOCATE(old_mb_forces(3,n_mols), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating force array")
     ALLOCATE(old_mb_torques(3,n_mols), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating torque array")
     ALLOCATE(old_atom_forces(3,n_ats), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating force array")
     ALLOCATE(mb_conj_forces(3,n_mols), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating force array")
     ALLOCATE(mb_conj_torques(3,n_mols), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating torque array")
     ALLOCATE(atom_conj_forces(3,n_ats), STAT = allostat)
     IF(allostat /= 0) CALL abort("allocating force array")

     ! velocities are meaningless in CG
     DO i = 1, n_mols
        molecules(i)%vel = 0.d0
        molecules(i)%angvel = 0.d0
     END DO
     DO i = 1, n_ats
        atoms(i)%vel = 0.d0
     END DO

  END IF

  ! Initialize arrays needed by the parallellization machinery
  CALL set_particles(n_mols,n_ats,n_gos)
  ! Determine the initial workload split. That is, find
  ! suitable sets of particles that each cpu should handle in
  ! the parallellized parts of the program.
  CALL initial_loadlimits()

  ! output
  IF(control_params%verbose > 0.d0 .AND. cpu_id == master_cpu) &
       WRITE(*,*) "Done"

  ! cpu time checkpoint
  CALL checkpoint(time_stat,timeA,timeB)

  !***************************
  !
  ! Initialization complete
  !
  !***************************




  !***********************************
  !
  ! Analyze system and write output
  !
  !***********************************

  ! Calculate different types of energy and deduce temperature, pressure etc.
  press = 0.d0
  CALL kin_energy(molecules,atoms,gos,n_freedom,e_kin,temps,e_lin,e_rot)
  CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
       mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
       mbgo_neighbors,mbgo_n_nbors, &
       mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
       physical_params,go_params,is_constrained,&
       e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
       e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)

  ! Calculate initial forces
  CALL calc_forces(molecules,atoms,gos,&
       mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
       mbgo_neighbors,mbgo_n_nbors, &
       mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
       control_params,physical_params,go_params,&
       is_constrained,.false.,&
       mb_forces,mb_torques,atom_forces,go_forces,virial,n_bond)

!!$  ! numerically calculate forces for debugging
!!$  shifty = 0.0001
!!$  do i = 1, 1
!!$     do j = 1, 3
!!$        molecules(i)%pos(j) = molecules(i)%pos(j)+shifty
!!$        CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
!!$             mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
!!$             mbgo_neighbors,mbgo_n_nbors, &
!!$             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
!!$             physical_params,go_params,is_constrained,&
!!$             e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
!!$             e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)
!!$        nep(j) = e_pot
!!$        molecules(i)%pos(j) = molecules(i)%pos(j)-2.d0*shifty
!!$        CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
!!$             mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
!!$             mbgo_neighbors,mbgo_n_nbors, &
!!$             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
!!$             physical_params,go_params,is_constrained,&
!!$             e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
!!$             e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)
!!$        nem(j) = e_pot
!!$        molecules(i)%pos(j) = molecules(i)%pos(j)+shifty
!!$     end do
!!$     write(*,*) "numeric and analytic force for mb ", i
!!$     write(*,'(F20.10,F20.10,F20.10)') &
!!$          -(nep(1)-nem(1))/(2.d0*shifty), &
!!$          -(nep(2)-nem(2))/(2.d0*shifty), &
!!$          -(nep(3)-nem(3))/(2.d0*shifty)
!!$     write(*,'(F20.10,F20.10,F20.10)') &
!!$          mb_forces(1:3,i)
!!$  end do
!!$  do i = 1, 5
!!$     do j = 1, 3
!!$        gos(i)%pos(j) = gos(i)%pos(j)+shifty
!!$        CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
!!$             mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
!!$             mbgo_neighbors,mbgo_n_nbors, &
!!$             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
!!$             physical_params,go_params,is_constrained,&
!!$             e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
!!$             e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)
!!$        nep(j) = e_pot
!!$        gos(i)%pos(j) = gos(i)%pos(j)-2.d0*shifty
!!$        CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
!!$             mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
!!$             mbgo_neighbors,mbgo_n_nbors, &
!!$             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
!!$             physical_params,go_params,is_constrained,&
!!$             e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
!!$             e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)
!!$        nem(j) = e_pot
!!$        gos(i)%pos(j) = gos(i)%pos(j)+shifty
!!$     end do
!!$     write(*,*) "numeric and analytic force for go ", i
!!$     write(*,'(F20.10,F20.10,F20.10)') &
!!$          -(nep(1)-nem(1))/(2.d0*shifty), &
!!$          -(nep(2)-nem(2))/(2.d0*shifty), &
!!$          -(nep(3)-nem(3))/(2.d0*shifty)
!!$     write(*,'(F20.10,F20.10,F20.10)') &
!!$          go_forces(1:3,i)
!!$  end do


  ! conjugate gradients
  IF(control_params%md_algo == cg_index)THEN
     CALL abort("conjugate gradient method has been disabled for this version of Cashew due to lack of support for go-particles")
     old_mb_forces = mb_forces
     old_mb_torques = mb_torques
     old_atom_forces = atom_forces
     mb_conj_forces = mb_forces
     mb_conj_torques = mb_torques
     atom_conj_forces = atom_forces
     last_e = e_pot + 2.d0*control_params%time_step

     ! copy particles to dummy lists
     DO j = 1, n_mols
        mb_dummy(j) = molecules(j)
     END DO
     DO j = 1, n_ats
        at_dummy(j) = atoms(j)
     END DO

  END IF

  ! cpu time checkpoint
  CALL checkpoint(time_stat,timeA,timeB)

  ! open outputs
  IF(cpu_id == master_cpu)THEN
     ! main putput
     output = 975075
     OPEN(output,FILE=system_name(1:namelength-1)//"."//MAIN_OUT)
     IF(control_params%n_statfiles > 0)THEN
        ! statistics outputs
        statout = 2515975
        DO i = 1, control_params%n_statfiles
           IF(control_params%n_statfiles == 1)THEN
              timelabel = real_to_string(REAL(i,KIND=dp))
              OPEN(statout+i,FILE=system_name(1:namelength-1)//STATS//&
                   "."//MAIN_OUT)
           ELSE
              timelabel = real_to_string(REAL(i,KIND=dp))
              OPEN(statout+i,FILE=system_name(1:namelength-1)//STATS//&
                   "_"//timelabel(8:10)//"."//MAIN_OUT)
           END IF
        END DO
     END IF
     
     ! write output: header
     CALL write_header(output,system_name(1:namelength-1),&
          molecules,atoms,gos,n_elements,n_types,tvector,&
          physical_params,control_params,go_params,&
          supercell,boundary_type,boundary_value,periodic_boundary,is_constrained,&
          e_kin,e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,temps,&
          e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat, &
          mb_forces,mb_torques,atom_forces,go_forces,&
          mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors)
     
     ! printout
     IF(control_params%verbose > 0.d0)THEN
        IF(control_params%runtype == full_index)THEN
           WRITE(*,*) ""
           WRITE(*,'(A)') "Starting simulation"
           WRITE(*,*) ""
        ELSE
           WRITE(*,*) ""
           WRITE(*,'(A)') "Configuration properties"
           WRITE(*,*) ""
        END IF
        WRITE(*,'(A,F20.4)') "Temperature:                ", temps
        WRITE(*,'(A,F20.4)') "Total energy:               ", e_kin+e_pot
        WRITE(*,'(A,F20.4)') "  kinetic:                  ", e_kin
        WRITE(*,'(A,F20.4)') "  potential:                ", e_pot
        WRITE(*,'(A,F20.4)') "     MB-MB Lennard-Jones:   ", e_lenjon
        WRITE(*,'(A,F20.4)') "     MB-MB hydrogen bond:   ", e_bonds
        IF(n_ats > 0)THEN
           WRITE(*,'(A,F20.4)') "     MB-atom potential:     ", e_mbat
           WRITE(*,'(A,F20.4)') "     atom-atom potential:   ", e_atat
        END IF
       IF(n_gos > 0)THEN
          WRITE(*,'(A,F20.4)') "     MB-GO potential:       ", e_mbgo
          WRITE(*,'(A,F20.4)') "     GO-GO potential:       ", e_go
          WRITE(*,'(A,F20.4)') "        GO stretching:      ", e_stretch
          WRITE(*,'(A,F20.4)') "        GO bending:         ", e_bend
          WRITE(*,'(A,F20.4)') "        GO torsion:         ", e_torsion
          WRITE(*,'(A,F20.4)') "        GO native:          ", e_native
          WRITE(*,'(A,F20.4)') "        GO non-native:      ", e_nonnat
       END IF
        IF(is_constrained)THEN
           WRITE(*,'(A,F20.4)') "     constraint potential:  ", e_constr
        END IF
        WRITE(*,*) ""
     END IF

     ! write stats
     DO k = 1, control_params%n_statfiles
        IF(control_params%stat_start(k) < control_params%time_step)THEN
           CALL write_stats(statout,molecules,atoms,temps,e_kin,e_pot,&
                e_lin,e_rot,e_lenjon,e_bonds,e_mbat,e_atat,virial,supercell,0.d0,"      00:00:00.000",&
                mb_forces,mb_torques,atom_forces,control_params,k)
        END IF
     END DO
     
  END IF

  ! cpu time checkpoint
  CALL checkpoint(time_io,timeA,timeB)





  !********************************************
  !
  ! The actual molecular dynamics part begins   
  !
  !********************************************

  ! In leapfrog, velocities and positions are calculated
  ! dt/2 apart. So, advance initial velocities by dt/2
  ! Note: leapfrog is not fully implemented, one should use velocity-Verlet
  IF(control_params%md_algo == leapfrog_index)THEN
     CALL update_velocities(molecules,atoms,gos,n_mols,n_ats,n_gos,&
          mb_forces,mb_torques,atom_forces,go_forces,&
          physical_params,control_params,go_params,0.5d0,is_constrained)       
  END IF

  ! cpu time checkpoint
  CALL checkpoint(time_move,timeA,timeB)

  ! set counters and timers
  simulation_time = 0.d0
  ! offset
  IF(control_params%md_algo == cg_index)THEN
     offset = 0.5d0
  ELSE
     offset = 0.5d0*control_params%time_step
  END IF
  ! screen/fileprint
  speaktimer = ABS(control_params%verbose)-offset
  ! stats
  IF(control_params%n_statfiles > 0)THEN
     ALLOCATE(stattimer(control_params%n_statfiles))
     DO k = 1, control_params%n_statfiles
        IF(control_params%stat_start(k) < offset)THEN
           stattimer(k) = control_params%stat_interval(k)-offset
        ELSE
           stattimer(k) = control_params%stat_start(k)-offset
        END IF
     END DO
  ELSE
     ALLOCATE(stattimer(1))
     stattimer(1) = control_params%time_max*2.d0
  END IF
  ! xyz
  IF(control_params%xyz_writer == ixyz_index)THEN
     xyztimer = control_params%xyz_interval-offset
  ELSE
     xyztimer = control_params%time_max*2.d0
  END IF
  ! bond
  IF(control_params%bond_writer == ixyz_index)THEN
     bondtimer = control_params%bond_interval-offset
  ELSE
     bondtimer = control_params%time_max*2.d0
  END IF
  ! continuation
  IF(control_params%inp_writer == ixyz_index)THEN
     inptimer = control_params%inp_interval-offset
  ELSE
     inptimer = control_params%time_max*2.d0
  END IF
  ! rdf
  IF(control_params%rdf_writer == ixyz_index)THEN
     rdftimer = control_params%rdf_interval-offset
  ELSE IF(control_params%rdf_writer == average_index)THEN
     rdftimer = control_params%rdf_start-offset
     rdf_ave = 0.d0
     rdf_sq = 0.d0
     rdf_count = 0
  ELSE
     rdftimer = control_params%time_max*2.d0
  END IF
  ! adf
  IF(control_params%adf_writer == ixyz_index)THEN
     adftimer = control_params%adf_interval-offset
  ELSE IF(control_params%adf_writer == average_index)THEN
     adftimer = control_params%adf_start-offset
     adf_ave = 0.d0
     adf_sq = 0.d0
     adf_count = 0
  ELSE
     adftimer = control_params%time_max*2.d0
  END IF
  ! barometer
  IF(control_params%md_baro /= microcanonical_index)THEN
     barotimer = (control_params%baro_interval*control_params%time_step-offset)
     barosteps = FLOOR((control_params%baro_interval*control_params%time_step+offset)/control_params%time_step)
  ELSE
     barotimer = control_params%time_max*2.d0
     barosteps = 1
  END IF
  ! at each system update, the max shift in any single coordinate is recorded
  ! and drift stores the sum of these max shifts
  drift = 0.d0
  ! drift_range is the maximum allowed value for drift before neighbor
  ! lists must be updated: neighbor lists are formed in spheres
  ! with a radius of r + cut_ver where r is the interaction range and
  ! cut_ver is an additional marginal. It is important that the neighbor lists
  ! always include all particles within r, so the list must be updated once a
  ! particle may have passed from the outside of r+cut_ver to inside of r.
  ! since drift is the sum of maximum shifts in single coordinates, any one
  ! particle can at max have shifted by drift in the direction of one coordinate axis,
  ! and by sqrt(3)*drift in total if all x, y, and z coordinates have changed this much.
  ! Since the neighbor lists involve two particles, both may have moved by this amount
  ! and the maximum possible change in their distance is 2*sqrt(3)*drift.
  ! a particle can be missing in a neighbor list if this distance is more than cut_ver, so
  ! the proper trigger for list update is 2*sqrt(3)*drift > cut_ver, i.e.,
  ! the maximum allowed drift is drift_range = cut_ver*0.5/sqrt(3).
  drift_range = physical_params%cut_ver*0.5d0/sqrt3

  ! cpu time checkpoint
  CALL checkpoint(time_stat,timeA,timeB)

  ! if required, write data based on the initial configuration;
  ! only the master cpu writes
  IF(cpu_id == master_cpu)THEN
     
     ! write the initial configuration to xyz
     IF(control_params%xyz_writer == sxyz_index .OR. &
          control_params%xyz_writer == sexyz_index .OR. &
          control_params%xyz_writer == ixyz_index )THEN
        CALL write_xyz(system_name(1:namelength-1)//".xyz",&
             molecules,atoms,gos,physical_params%l_oh,.true.)
     END IF
     ! write the initial configuration as bonds
     IF(control_params%bond_writer == sxyz_index .OR. &
          control_params%bond_writer == sexyz_index .OR. &
          control_params%bond_writer == ixyz_index )THEN
        CALL write_bonds(system_name(1:namelength-1)//MB_MID//".xyz",system_name(1:namelength-1)//BOND_MID//".out",n_bond,&
             molecules,mbmb_neighbors,mbmb_n_nbors,supercell,periodic_boundary,control_params,physical_params,.true.)
     END IF
     ! write the initial configuration as an inputfile
     IF(control_params%inp_writer == sxyz_index .OR. &
          control_params%inp_writer == sexyz_index .OR. &
          control_params%inp_writer == ixyz_index )THEN
        CALL write_contfile(system_name(1:namelength-1)//CONT_START//"."//MB_IN,&
             molecules,atoms,gos,n_elements,&
             physical_params,control_params,go_params,&
             supercell,boundary_type,boundary_value,is_constrained)
     END IF
     
     ! cpu time checkpoint
     CALL checkpoint(time_io,timeA,timeB)
     
     ! write the initial rdf
     IF(control_params%rdf_writer == sxyz_index .OR. &
          control_params%rdf_writer == sexyz_index .OR. &
          control_params%rdf_writer == ixyz_index )THEN
        CALL write_rdf(system_name(1:namelength-1)//RDF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params,.true.)
     ELSE IF(control_params%rdf_writer == sxyz_index .AND. &
          control_params%rdf_start <= norm_tolerance)THEN
        CALL radial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params%rdf_particles,control_params%rdf_range,rdf_temp)
        rdf_ave = rdf_ave + rdf_temp
        rdf_sq = rdf_sq + rdf_temp*rdf_temp
        rdf_count = rdf_count+1
     END IF
     ! write the initial adf
     IF(control_params%adf_writer == sxyz_index .OR. &
          control_params%adf_writer == sexyz_index .OR. &
          control_params%adf_writer == ixyz_index )THEN
        CALL write_adf(system_name(1:namelength-1)//ADF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,&
             physical_params,control_params,.true.)
     ELSE IF(control_params%adf_writer == sxyz_index .AND. &
          control_params%adf_start <= norm_tolerance)THEN
        CALL axial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,&
             physical_params,control_params,adf_temp)
        adf_ave = adf_ave + adf_temp
        adf_sq = adf_sq + adf_temp*adf_temp
        adf_count = adf_count+1
     END IF
  END IF

  ! cpu time checkpoint
  CALL checkpoint(time_stat,timeA,timeB)

  ! end here if this is just a test run
  IF(control_params%runtype == test_index)THEN
     time_tot = timeA-time_tot
     ! Close the main input file
     IF(cpu_id == master_cpu)THEN
        CALL write_end(output,tvector,&
             time_io, time_force, time_stat, time_nbor, time_move, time_tot)
        CLOSE(output)
        IF(control_params%n_statfiles > 0)THEN
           DO k = 1, control_params%n_statfiles
              CLOSE(statout+k)
           END DO
        END IF
     END IF
     ! close mpi
     CALL finish_mpi()
     STOP    
  END IF






  !**************************************!
  ! Time advancing, main simulation loop !
  !**************************************!

  ! switch for ending the loop
  finish = .false.

  ! run as long as the maximum simulation time has not been reached,
  ! and the cg has not converged
  DO WHILE(simulation_time < control_params%time_max-offset .AND. .NOT.finish)

     IF(control_params%md_algo == cg_index)THEN ! conjugate gradients

        !******************
        ! Do a CG step
        !******************
        CALL cg_step(molecules,atoms,mb_dummy,at_dummy,n_mols,n_ats,&
             supercell,periodic_boundary,boundary_type,boundary_value,&
             physical_params,control_params,&
             mb_bonds,is_constrained,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
             boxed,boxes,box_mbs,box_ats,box_mb_count,box_at_count,&
             mb_forces,mb_torques,atom_forces,&
             old_mb_forces,old_mb_torques,old_atom_forces,&
             mb_conj_forces,mb_conj_torques,atom_conj_forces,&
             last_e,simulation_time,finish,cgbeta,n_bond)

     ELSE ! molecular dynamics

        !*********************
        ! Advance time in MD
        !*********************
        CALL advance_time(molecules,atoms,gos,n_mols,n_ats,n_gos,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,&
             atat_neighbors,atat_n_nbors,mbgo_neighbors,mbgo_n_nbors,&
             boxed,boxes,box_mbs,box_ats,box_gos,box_mb_count,box_at_count,box_go_count,&
             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
             is_constrained,&
             mb_forces,mb_torques,atom_forces,go_forces,virial,n_bond,&
             physical_params,control_params,go_params,&
             drift,drift_range,simulation_time)
       
        ! if a barostat is in use, apply it
        IF(control_params%md_baro /= microcanonical_index)THEN
           CALL apply_barostat()
        END IF

    END IF

     !***************************
     !  Output routines
     !***************************

    CALL interval_outputs()
    
  END DO

  !**************************************!
  !   End of the main simulation loop    !
  !**************************************!





  !***************************
  !  Output routines
  !***************************


  ! write output based on the final configuration;
  ! only the master cpu writes
  IF(cpu_id == master_cpu)THEN
     ! write the final configuration to xyz
     IF(control_params%xyz_writer == exyz_index)THEN
        CALL write_xyz(system_name(1:namelength-1)//".xyz",&
             molecules,atoms,gos,physical_params%l_oh,.true.) ! only end xyz is written
     ELSE IF(control_params%xyz_writer == sexyz_index )THEN
        CALL write_xyz(system_name(1:namelength-1)//".xyz",&
             molecules,atoms,gos,physical_params%l_oh,.false.) ! start xyz already exists
     END IF

     ! write the final configuration as bonds
     IF(control_params%bond_writer == exyz_index)THEN
        CALL write_bonds(system_name(1:namelength-1)//MB_MID//".xyz",system_name(1:namelength-1)//BOND_MID//".out",n_bond,&
             molecules,mbmb_neighbors,mbmb_n_nbors,supercell,periodic_boundary,control_params,physical_params,.true.)
     ELSE IF(control_params%bond_writer == sexyz_index )THEN
        CALL write_bonds(system_name(1:namelength-1)//MB_MID//".xyz",system_name(1:namelength-1)//BOND_MID//".out",n_bond,&
             molecules,mbmb_neighbors,mbmb_n_nbors,supercell,periodic_boundary,control_params,physical_params,.false.)
     END IF
     
     ! write the final configuration as an inputfile
     IF(control_params%inp_writer == exyz_index .OR. &
          control_params%inp_writer == sexyz_index )THEN
        CALL write_contfile(system_name(1:namelength-1)//CONT_END//"."//MB_IN,&
             molecules,atoms,gos,n_elements,&
             physical_params,control_params,go_params,&
             supercell,boundary_type,boundary_value,is_constrained)
     END IF
     
     ! cpu time checkpoint
     CALL checkpoint(time_io,timeA,timeB)
     
     ! write the final configuration to rdf
     IF(control_params%rdf_writer == exyz_index)THEN
        CALL write_rdf(system_name(1:namelength-1)//RDF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params,.true.)
     ELSE IF(control_params%rdf_writer == sexyz_index)THEN
        CALL write_rdf(system_name(1:namelength-1)//RDF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params,.false.)
     ELSE IF(control_params%rdf_writer == average_index)THEN
        CALL radial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params%rdf_particles,control_params%rdf_range,rdf_temp)
        rdf_ave = rdf_ave + rdf_temp
        rdf_sq = rdf_sq + rdf_temp*rdf_temp
        rdf_count = rdf_count+1
        rdf_ave = rdf_ave/REAL(rdf_count,KIND=dp)
        rdf_sq = rdf_sq/REAL(rdf_count,KIND=dp)
        CALL write_rdf(system_name(1:namelength-1)//RDF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,&
             mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
             control_params,.true.,&
             rdf_ave,sqrt((rdf_sq-rdf_ave*rdf_ave)/REAL(MAX(1,rdf_count-1),KIND=dp)))
     END IF
     
     ! write the final configuration to adf
     IF(control_params%adf_writer == exyz_index)THEN
        CALL write_adf(system_name(1:namelength-1)//ADF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,physical_params,&
             control_params,.true.)
     ELSE IF(control_params%adf_writer == sexyz_index)THEN
        CALL write_adf(system_name(1:namelength-1)//ADF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,physical_params,&
             control_params,.false.)
     ELSE IF(control_params%adf_writer == average_index)THEN
        CALL axial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,physical_params,&
             control_params,adf_temp)
        adf_ave = adf_ave + adf_temp
        adf_sq = adf_sq + adf_temp*adf_temp
        adf_count = adf_count+1
        adf_ave = adf_ave/REAL(adf_count,KIND=dp)
        adf_sq = adf_sq/REAL(adf_count,KIND=dp)
        CALL write_adf(system_name(1:namelength-1)//ADF_MID//".out",&
             molecules,atoms,supercell,periodic_boundary,physical_params,&
             control_params,.true.,&
             adf_ave,sqrt((adf_sq-adf_ave*adf_ave)/REAL(MAX(1,adf_count-1),KIND=dp)))
     END IF
  
     ! cpu time checkpoint
     CALL checkpoint(time_stat,timeA,timeB)
     time_tot = timeA-time_tot
     
     !
     ! Close the main input file
     !
     CALL write_end(output,tvector,&
          time_io, time_force, time_stat, time_nbor, time_move, time_tot)
     CLOSE(output)
  END IF
  ! close mpi
  CALL finish_mpi()
  STOP

CONTAINS

  ! advances time by one timestep. this is done by calling
  ! the routines "move_particles" and "update_velocities".
  ! the forces are also updated within the routine by
  ! calling the routine "calc_forces" in the mb_model module.
  ! *mbs MB molecules
  ! *ats atomic particles
  ! *n_mbs number of MB molecules
  ! *n_ats number of atomic particles
  ! *n_elems number of different types of atoms
  ! *mm_nbors MB-MB list of neighbors
  ! *ma_nbors MB-atom list of neighbors
  ! *aa_nbors atom-atom list of neighbors
  ! *mm_n_ns MB-MB number of neighbors
  ! *ma_n_ns MB-atom number of neighbors
  ! *aa_n_ns atom-atom number of neighbors
  ! *bonds number of bonds for each MB molecule (z_i)
  ! *cell supercell dimensions (x,y,z)
  ! *pbc periodic boundary conditions (x,y,z)
  ! *btype boundary type
  ! *bval value associated with the boundary
  ! *pair_done logical list for MB-MB pairs - used for preventing double counting (obsolete)
  ! *atom_pair_done logical list for atom-atom pairs (obsolete)
  ! *is_constr true if constraints have been defined, false otherwise
  ! *m_force forces acting on MB molecules
  ! *m_torque torques acting on MB molecules
  ! *a_force forces acting on atomic particles
  ! *params physical parameters
  ! *control control parameters (including the timestep)
  ! *drift the maximum change in a position coordinate. needed for optimal updating frequency of neighbors lists
  ! *clock virtual time of the simulation
  SUBROUTINE advance_time(mbs,ats,gos,n_mbs,n_ats,n_gos,&
       mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
       boxed,boxes,box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
       bonds,cell,pbc,btype,bval,is_constr,&
       m_force,m_torque,a_force,g_force,virial,n_bond,params,control,go,drift,drift_range,clock)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:), g_force(:,:)
    REAL(KIND=dp), INTENT(INOUT) :: drift, clock, virial
    INTEGER, INTENT(IN) :: n_mbs, n_ats, n_gos, btype(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3), drift_range
    LOGICAL, INTENT(IN) :: pbc(3), is_constr
    !LOGICAL, POINTER :: pair_done(:,:), atom_pair_done(:,:)
    REAL(KIND=dp), POINTER :: bonds(:), n_bond(:)
    INTEGER, POINTER :: &
         mm_nbors(:,:), mm_n_ns(:), &
         ma_nbors(:,:), ma_n_ns(:), &
         aa_nbors(:,:), aa_n_ns(:), &
         mg_nbors(:,:), mg_n_ns(:), &
         box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), box_gos(:,:,:,:,:),&
         box_mc(:,:,:,:), box_ac(:,:,:,:), box_gc(:,:,:,:)
    INTEGER, INTENT(IN) :: boxes(3,2)
    LOGICAL, INTENT(IN) :: boxed
    REAL(KIND=dp) :: shift

    ! increase the elapsed virtual time
    clock = clock + control%time_step

    ! velocity verlet - the main algorithm of the program
    IF(control%md_algo == velocity_verlet_index)THEN

       ! move by dt (from t to t+dt)
       CALL move_particles(mbs,ats,gos,n_mbs,n_ats,n_gos,cell,pbc,&
            m_force,m_torque,a_force,g_force,params,control,go,is_constr,shift)
       ! keep account on how much particles move, for triggering neighbor list updating
       drift = drift+shift

       ! notify of very large shifts.
       ! this may happen if the timestep is too large or if 
       ! the particles are initially too close to each other
       IF(shift > 0.5d0*params%max_cut)THEN
          CALL abort("huge shifts in atomic placements")
       END IF
       IF(shift > params%cut_ver)THEN
          WRITE(*,*) "Warning: large shifts in atomic placements"
       END IF

       ! cpu time checkpoint
       CALL checkpoint(time_move,timeA,timeB)

       IF(drift > drift_range)THEN ! update neighbor lists if the particles have moved enough
          drift = 0.d0
          CALL update_neighbors(mbs,ats,gos,cell,pbc,bonds,&
               boxed,boxes,box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
               mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
               params,go)
       ELSE
          ! the bond numbers needed for calculating the bond-order
          ! term must be updated after each step
          CALL update_bond_numbers_p(mbs,cell,pbc,mm_nbors,mm_n_ns,bonds,params%R_b,params%D_b,params%inv_Db)
       END IF

       ! cpu time checkpoint
       CALL checkpoint(time_nbor,timeA,timeB) 

       ! update vel by dt/2 (from t to t+dt/2) with old forces (t)
       CALL update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
            m_force,m_torque,a_force,g_force,params,control,go,0.5d0,is_constr)

       ! calculate new forces (t+dt) without thermostat
       CALL calc_forces(mbs,ats,gos,&
            mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
            bonds,cell,pbc,btype,bval,control,params,go,is_constr,.false.,&
            m_force,m_torque,a_force,g_force,virial,n_bond)

       ! cpu time checkpoint       
       CALL checkpoint(time_force,timeA,timeB)

       IF(control%md_thermo == langevin_index)THEN ! Langevin thermostat

          ! add random forces but no friction (t+dt)
          CALL add_random_force(ats,n_mbs,n_ats,n_gos,control,params,go,m_force,m_torque,a_force,g_force)
          ! cpu time checkpoint
          CALL checkpoint(time_force,timeA,timeB)
          
          ! update vel by dt/2 (from t+dt/2 to t+dt) with new forces (t+dt) with "future friction"
          CALL update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
               m_force,m_torque,a_force,g_force,params,control,go,0.5d0,is_constr,.true.)
          ! cpu time checkpoint
          CALL checkpoint(time_move,timeA,timeB)
          
          ! add friction to the force (t+dt)
          CALL add_friction_force(mbs,ats,gos,n_mbs,n_ats,n_gos,control,go,m_force,m_torque,a_force,g_force)
          ! cpu time checkpoint
          CALL checkpoint(time_force,timeA,timeB)

       ELSE IF(control%md_thermo == cooler_index )THEN ! cooler "thermostat" - only friction is applied resulting in decreasing temperature

          ! update vel by dt/2 (from t+dt/2 to t+dt) with new forces (t+dt) including "future friction"
          CALL update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
               m_force,m_torque,a_force,g_force,params,control,go,0.5d0,is_constr,.true.)
          ! cpu time checkpoint
          CALL checkpoint(time_move,timeA,timeB)

          ! add friction to the force (t+dt)
          CALL add_friction_force(mbs,ats,gos,n_mbs,n_ats,n_gos,control,go,m_force,m_torque,a_force,g_force)
          ! cpu time checkpoint
          CALL checkpoint(time_force,timeA,timeB)

       ELSE ! no thermostat

          ! update vel by dt/2 (from t+dt/2 to t+dt) with new forces (t+dt)
          CALL update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
               m_force,m_torque,a_force,g_force,params,control,go,0.5d0,is_constr)
          ! cpu time checkpoint
          CALL checkpoint(time_move,timeA,timeB)

       END IF

    ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog algorithm: this is here for testing, e.g. termostats have not been implemented

       ! move by dt (from t to t+dt) according to velocities at t+dt/2
       CALL move_particles(mbs,ats,gos,n_mbs,n_ats,n_gos,cell,pbc,&
            m_force,m_torque,a_force,g_force,params,control,go,is_constr,shift)
       ! keep account on how much particles move for neighbor list updating
       drift = drift+shift

       IF(shift > 0.5d0*params%max_cut)THEN
          CALL abort("huge shifts in atomic placements")
       END IF
       IF(shift > params%cut_ver)THEN
          WRITE(*,*) "Warning: large shifts in atomic placements"
       END IF

       ! cpu time checkpoint
       CALL checkpoint(time_move,timeA,timeB) 

       IF(drift > drift_range)THEN ! update neighbor lists if the particles have moved enough
          drift = 0.d0
          CALL update_neighbors(mbs,ats,gos,cell,pbc,bonds,&
               boxed,boxes,box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
               mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
               params,go)
       ELSE
          CALL update_bond_numbers(mbs,cell,pbc,mm_nbors,mm_n_ns,bonds,params%R_b,params%D_b,params%inv_Db)
       END IF

       ! cpu time checkpoint
       CALL checkpoint(time_nbor,timeA,timeB)      
       
       ! calculate new forces (t+dt)
       CALL calc_forces(mbs,ats,gos,&
            mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
            bonds,cell,pbc,btype,bval,control,params,go,is_constr,.true.,&
            m_force,m_torque,a_force,g_force,virial,n_bond)

       ! cpu time checkpoint
       CALL checkpoint(time_force,timeA,timeB)
       
       ! update vel by dt (from t+dt/2 to t+dt+dt/2) with new forces (t+dt)
       CALL update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
            m_force,m_torque,a_force,g_force,params,control,go,1.0d0,is_constr)

       ! cpu time checkpoint
       CALL checkpoint(time_move,timeA,timeB)

    END IF

    RETURN
  END SUBROUTINE advance_time

  ! Updates the velocities of all MB molecules and atoms according to the
  ! forces given as input. As algorithms may require advancing the forces
  ! by either a full or a half timestep, the variable timestep_fraction
  ! has been defined to control the amount of time by which the forces
  ! are updated. (The actual timestep is defined in the control parameter
  ! type.)
  ! possible constraints are also taken into account.
  ! *mbs MB molecules
  ! *ats atomic particles
  ! *n_mbs number of MB molecules
  ! *n_ats number of atomic particles
  ! *n_elems number of different types of atoms
  ! *is_constr true if constraints have been defined, false otherwise
  ! *m_force forces acting on MB molecules
  ! *m_torque torques acting on MB molecules
  ! *a_force forces acting on atomic particles
  ! *params physical parameters
  ! *control control parameters (including the timestep)
  ! *timestep_fraction 1.0 for a full timestep, 0.5 for a half step (other values work similarly, but should not be needed)
  ! *future friction if given and true, the velocities are updated assuming forces from a future time including v-dependent friction: v(t+dt) = v(t) + a(t+dt)*dt = v(t) + a'(t+dt)*dt - g*v(t+dt)*dt -> v(t+dt) = (v(t) + a'(t+dt)*dt)/(1-g*dt)
  SUBROUTINE update_velocities(mbs,ats,gos,n_mbs,n_ats,n_gos,&
       m_force,m_torque,a_force,g_force,params,control,go,timestep_fraction,is_constr,future_friction)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:), g_force(:,:)
    REAL(KIND=dp), INTENT(IN) :: timestep_fraction
    LOGICAL, INTENT(IN) :: is_constr
    INTEGER, INTENT(IN) :: n_mbs, n_ats, n_gos
    INTEGER :: ii, jj
    REAL(KIND=dp) :: dt, friction_corr!, friction_corr2
    LOGICAL, INTENT(IN), OPTIONAL :: future_friction
    LOGICAL :: include_friction

    ! the delta t in this step (in practice dt or dt/2)
    dt = control%time_step * timestep_fraction
    friction_corr = 1.d0

    include_friction = .false.
    ! Future friction means that a frictional force is applied and we need to
    ! update the velocity according to F(t+dt) including v(t+dt).
    ! Here, v(t+dt) can be solved from the equation resulting in
    ! scaling coefficient that is less than one due to the friction.
    IF(present(future_friction))THEN
       IF(future_friction)THEN
          include_friction = .true.
          friction_corr = 1.d0/(1.d0 + dt*control%thermo_value)
       END IF
    END IF

    ! constrained update
    IF(is_constr)THEN
       DO ii = 1, n_mbs ! loop over molecules
          DO jj = 1, 3 ! loop over x, y, z

             ! select the type of constraint
             SELECT CASE (mbs(ii)%constrained(jj)) 
             CASE(frozen_pos_index) ! frozen position -> v = 0
                mbs(ii)%vel(jj) = 0.d0
                mbs(ii)%angvel(jj) = mbs(ii)%angvel(jj) + dt*params%inv_i*m_torque(jj,ii)

             CASE(frozen_vel_index) ! frozen velocity -> v = const.
                ! don't change the velocity component
                !mbs(ii)%vel(jj) = mbs(ii)%vel(jj)                
                mbs(ii)%angvel(jj) = mbs(ii)%angvel(jj) + dt*params%inv_i*m_torque(jj,ii)

             CASE(all_frozen_index) ! all frozen -> v = 0, w = 0
                mbs(ii)%vel(jj) = 0.d0
                mbs(ii)%angvel(jj) = 0.d0

             CASE DEFAULT ! no constraints on the velocity
                mbs(ii)%vel(jj) = mbs(ii)%vel(jj) + dt*params%inv_m*m_force(jj,ii)
                mbs(ii)%angvel(jj) = mbs(ii)%angvel(jj) + dt*params%inv_i*m_torque(jj,ii)

             END SELECT
          END DO

          ! apply the scaling coefficient for "future friction"
          IF(include_friction)THEN 
             mbs(ii)%vel = friction_corr*mbs(ii)%vel
             mbs(ii)%angvel = friction_corr*mbs(ii)%angvel
          END IF

       END DO
       DO ii = 1, n_ats ! loop atoms
          DO jj = 1, 3 ! loop x, y, z

             ! select constraint
             SELECT CASE (ats(ii)%constrained(jj)) 
             CASE(frozen_pos_index) ! frozen position -> v = 0
                ats(ii)%vel(jj) = 0.d0

             CASE(frozen_vel_index) ! frozen velocity -> v = const.
                ! don't change the velocity component
                !ats(ii)%vel(jj) = ats(ii)%vel(jj)

             CASE(all_frozen_index) ! frozen all -> v = 0
                ats(ii)%vel(jj) = 0.d0

             CASE DEFAULT ! none
                ats(ii)%vel(jj) = ats(ii)%vel(jj) + dt*params%invm_atoms(ats(ii)%type)*a_force(jj,ii)
             END SELECT
          END DO

          ! apply the scaling coefficient for "future friction"
          IF(include_friction)THEN
             ats(ii)%vel = friction_corr*ats(ii)%vel
             !ats(ii)%vel = ats(ii)%vel/(1.d0 + dt*control%thermo_value*params%invm_atoms(ats(ii)%type))
          END IF

       END DO



       DO ii = 1, n_gos ! loop gos
          DO jj = 1, 3 ! loop x, y, z

             ! select constraint
             SELECT CASE (gos(ii)%constrained(jj)) 
             CASE(frozen_pos_index) ! frozen position -> v = 0
                gos(ii)%vel(jj) = 0.d0

             CASE(frozen_vel_index) ! frozen velocity -> v = const.
                ! don't change the velocity component
                !ats(ii)%vel(jj) = ats(ii)%vel(jj)

             CASE(all_frozen_index) ! frozen all -> v = 0
                gos(ii)%vel(jj) = 0.d0

             CASE DEFAULT ! none
                gos(ii)%vel(jj) = gos(ii)%vel(jj) + dt*g_force(jj,ii)/gos(ii)%mass
             END SELECT
          END DO

          ! apply the scaling coefficient for "future friction"
          IF(include_friction)THEN
             gos(ii)%vel = friction_corr*gos(ii)%vel
          END IF

       END DO



    ELSE ! unconstrained system
       DO ii = 1, n_mbs ! loop molecules

          ! update velocity: v += F/m*dt
          mbs(ii)%vel = mbs(ii)%vel + dt*params%inv_m*m_force(1:3,ii)
          mbs(ii)%angvel = mbs(ii)%angvel + dt*params%inv_i*m_torque(1:3,ii)

          ! update velocity: v += F/m*dt
          IF(include_friction)THEN
             mbs(ii)%vel = friction_corr*mbs(ii)%vel
             mbs(ii)%angvel = friction_corr*mbs(ii)%angvel
          END IF

       END DO
       DO ii = 1, n_ats ! loop atoms

          ! update velocity: v += F/m*dt
          ats(ii)%vel = ats(ii)%vel + dt*params%invm_atoms(ats(ii)%type)*a_force(1:3,ii)

          ! update velocity: v += F/m*dt
          IF(include_friction)THEN
             ats(ii)%vel = friction_corr*ats(ii)%vel
          END IF

       END DO
       DO ii = 1, n_gos ! loop atoms

          ! update velocity: v += F/m*dt
          gos(ii)%vel = gos(ii)%vel + dt*g_force(1:3,ii)/gos(ii)%mass

          ! update velocity: v += F/m*dt
          IF(include_friction)THEN
             gos(ii)%vel = friction_corr*gos(ii)%vel
          END IF

       END DO
    END IF

    RETURN
  END SUBROUTINE update_velocities

  ! moves the particles according to the current velocities and
  ! the forces acting on them (if using the velocity-verlet algorithm).
  ! possible constraints are also taken into account.
  ! *mbs MB molecules
  ! *ats atomic particles
  ! *n_mbs number of MB molecules
  ! *n_ats number of atomic particles
  ! *n_elems number of different types of atoms
  ! *cell supercell dimensions (x,y,z)
  ! *pbc periodic boundary conditions (x,y,z)
  ! *is_constr true if constraints have been defined, false otherwise
  ! *m_force forces acting on MB molecules
  ! *m_torque torques acting on MB molecules
  ! *a_force forces acting on atomic particles
  ! *params physical parameters
  ! *control control parameters (including the timestep)
  ! *max_shift the maximum change in a position coordinate. needed for optimal updating frequency of neighbors lists 
  SUBROUTINE move_particles(mbs,ats,gos,n_mbs,n_ats,n_gos,cell,pbc,&
       m_force,m_torque,a_force,g_force,params,control,go,is_constr,max_shift)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:), g_force(:,:)
    REAL(KIND=dp), INTENT(OUT) :: max_shift
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    INTEGER, INTENT(IN) :: n_mbs, n_ats, n_gos
    LOGICAL, INTENT(IN) :: is_constr, pbc(3)
    INTEGER :: ii, jj
    REAL(KIND=dp) :: dtdt, spin(3), shift(3)

    ! measure for the maximum shift in a single coordinate of a particle
    max_shift = 0.d0
    ! store dt**2/2
    dtdt = 0.5d0*control%time_step**2

    ! update the positions and orientations
    IF(is_constr)THEN ! constrained system
       DO ii = 1, n_mbs ! loop molecules
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm r += v*dt + a*dt**2/2
             DO jj = 1, 3 ! loop x,y,z

                ! calculate the shifts and rotations
                SELECT CASE (mbs(ii)%constrained(jj)) ! select constraint
                CASE(frozen_pos_index) ! frozen position -> r = const
                   shift(jj) = 0.d0
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj) + dtdt*params%inv_i*m_torque(jj,ii)
                CASE(all_frozen_index) ! all frozen -> r = const, q = const
                   shift(jj) = 0.d0
                   spin(jj) = 0.d0
                CASE(frozen_vel_index) ! frozen velocity -> v = const
                   shift(jj) = control%time_step*mbs(ii)%vel(jj) ! no acceleration term
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj) + dtdt*params%inv_i*m_torque(jj,ii)
                CASE DEFAULT ! no constraint on positions
                   shift(jj) = control%time_step*mbs(ii)%vel(jj) + dtdt*params%inv_m*m_force(jj,ii)
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj) + dtdt*params%inv_i*m_torque(jj,ii)
                END SELECT
             END DO
             
             ! shift the molecule
             mbs(ii)%pos = mbs(ii)%pos + shift
             ! account for periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(mbs(ii)%pos(jj) > cell(jj))THEN
                      mbs(ii)%pos(jj) = mbs(ii)%pos(jj) - cell(jj)
                      mbs(ii)%inipos(jj) = mbs(ii)%inipos(jj) - cell(jj)
                   END IF
                   IF(mbs(ii)%pos(jj) < 0.d0)THEN
                      mbs(ii)%pos(jj) = mbs(ii)%pos(jj) + cell(jj)
                      mbs(ii)%inipos(jj) = mbs(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO

             ! rotate molecule
             CALL rotate_molecule(mbs(ii),vec2q(spin))

             ! store the largest shift in a single coordinate
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog algorithm (no acceleration term involved)
             DO jj = 1, 3
                SELECT CASE (mbs(ii)%constrained(jj)) 
                CASE(frozen_pos_index)
                   shift(jj) = 0.d0
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj)
                CASE(all_frozen_index)
                   shift(jj) = 0.d0
                   spin(jj) = 0.d0
                CASE(frozen_vel_index)
                   shift(jj) = control%time_step*mbs(ii)%vel(jj)
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj)
                CASE DEFAULT
                   shift(jj) = control%time_step*mbs(ii)%vel(jj)
                   spin(jj) = control%time_step*mbs(ii)%angvel(jj)
                END SELECT
             END DO
             mbs(ii)%pos = mbs(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(mbs(ii)%pos(jj) > cell(jj))THEN
                      mbs(ii)%pos(jj) = mbs(ii)%pos(jj) - cell(jj)
                      mbs(ii)%inipos(jj) = mbs(ii)%inipos(jj) - cell(jj)
                   END IF
                   IF(mbs(ii)%pos(jj) < 0.d0)THEN
                      mbs(ii)%pos(jj) = mbs(ii)%pos(jj) + cell(jj)
                      mbs(ii)%inipos(jj) = mbs(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO
             CALL rotate_molecule(mbs(ii),vec2q(spin))
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO
       
       DO ii = 1, n_ats ! loop atoms
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm
             DO jj = 1, 3 ! loop x, y, z

                SELECT CASE (ats(ii)%constrained(jj)) ! select constraint
                CASE(frozen_pos_index) ! frozen position
                   shift(jj) = 0.d0
                CASE(frozen_vel_index) ! frozen velocity
                   shift(jj) = control%time_step*ats(ii)%vel(jj) ! no acceleration term
                CASE(all_frozen_index) ! all frozen
                   shift(jj) = 0.d0
                CASE DEFAULT ! none
                   shift(jj) = control%time_step*ats(ii)%vel(jj) + dtdt*params%invm_atoms(ats(ii)%type)*a_force(jj,ii)
                END SELECT
             END DO

             ! shift the position
             ats(ii)%pos = ats(ii)%pos + shift
             ! account for periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(ats(ii)%pos(jj) > cell(jj))THEN
                      ats(ii)%pos(jj) = ats(ii)%pos(jj) - cell(jj)
                      ats(ii)%inipos(jj) = ats(ii)%inipos(jj) - cell(jj)                      
                   END IF
                   IF(ats(ii)%pos(jj) < 0.d0)THEN
                      ats(ii)%pos(jj) = ats(ii)%pos(jj) + cell(jj)
                      ats(ii)%inipos(jj) = ats(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog algorithm
             DO jj = 1, 3
                SELECT CASE (ats(ii)%constrained(jj)) 
                CASE(frozen_pos_index)
                   shift(jj) = 0.d0
                CASE(frozen_vel_index)
                   shift(jj) = control%time_step*ats(ii)%vel(jj)
                CASE(all_frozen_index)
                   shift(jj) = 0.d0
                CASE DEFAULT
                   shift(jj) = control%time_step*ats(ii)%vel(jj)
                END SELECT
             END DO
             shift = control%time_step*ats(ii)%vel
             ats(ii)%pos = ats(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(ats(ii)%pos(jj) > cell(jj))THEN
                      ats(ii)%pos(jj) = ats(ii)%pos(jj) - cell(jj)
                      ats(ii)%inipos(jj) = ats(ii)%inipos(jj) - cell(jj)                      
                   END IF
                   IF(ats(ii)%pos(jj) < 0.d0)THEN
                      ats(ii)%pos(jj) = ats(ii)%pos(jj) + cell(jj)
                      ats(ii)%inipos(jj) = ats(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO



       DO ii = 1, n_gos ! loop gos
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm
             DO jj = 1, 3 ! loop x, y, z

                SELECT CASE (gos(ii)%constrained(jj)) ! select constraint
                CASE(frozen_pos_index) ! frozen position
                   shift(jj) = 0.d0
                CASE(frozen_vel_index) ! frozen velocity
                   shift(jj) = control%time_step*gos(ii)%vel(jj) ! no acceleration term
                CASE(all_frozen_index) ! all frozen
                   shift(jj) = 0.d0
                CASE DEFAULT ! none
                   shift(jj) = control%time_step*gos(ii)%vel(jj) + &
                        dtdt*g_force(jj,ii)/gos(ii)%mass
                END SELECT
             END DO

             ! shift the position
             gos(ii)%pos = gos(ii)%pos + shift
             ! account for periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(gos(ii)%pos(jj) > cell(jj))THEN
                      gos(ii)%pos(jj) = gos(ii)%pos(jj) - cell(jj)
                      gos(ii)%inipos(jj) = gos(ii)%inipos(jj) - cell(jj)                      
                   END IF
                   IF(gos(ii)%pos(jj) < 0.d0)THEN
                      gos(ii)%pos(jj) = gos(ii)%pos(jj) + cell(jj)
                      gos(ii)%inipos(jj) = gos(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog algorithm
             DO jj = 1, 3
                SELECT CASE (gos(ii)%constrained(jj)) 
                CASE(frozen_pos_index)
                   shift(jj) = 0.d0
                CASE(frozen_vel_index)
                   shift(jj) = control%time_step*gos(ii)%vel(jj)
                CASE(all_frozen_index)
                   shift(jj) = 0.d0
                CASE DEFAULT
                   shift(jj) = control%time_step*gos(ii)%vel(jj)
                END SELECT
             END DO
             shift = control%time_step*gos(ii)%vel
             gos(ii)%pos = gos(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(gos(ii)%pos(jj) > cell(jj))THEN
                      gos(ii)%pos(jj) = gos(ii)%pos(jj) - cell(jj)
                      gos(ii)%inipos(jj) = gos(ii)%inipos(jj) - cell(jj)                      
                   END IF
                   IF(gos(ii)%pos(jj) < 0.d0)THEN
                      gos(ii)%pos(jj) = gos(ii)%pos(jj) + cell(jj)
                      gos(ii)%inipos(jj) = gos(ii)%inipos(jj) + cell(jj)
                   END IF
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO


    ELSE ! no constraints


       DO ii = 1, n_mbs ! loop molecules
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm
             ! shift: r += v*dt + a*dt**2/2
             shift = control%time_step*mbs(ii)%vel + dtdt*params%inv_m*m_force(1:3,ii)
             ! shift the molecular position
             mbs(ii)%pos = mbs(ii)%pos + shift
             ! account for periodicity
             DO jj = 1, 3 ! loop x,y,z
                IF(pbc(jj))THEN ! periodic boundary
                   ! shift the molecule by the length of the cell
                   IF(mbs(ii)%pos(jj) > cell(jj)) mbs(ii)%pos(jj) = mbs(ii)%pos(jj) - cell(jj)
                   IF(mbs(ii)%pos(jj) < 0.d0)     mbs(ii)%pos(jj) = mbs(ii)%pos(jj) + cell(jj)
                END IF
             END DO

             ! rotation: phi = omega*dt + alpha*dt**2/2
             spin = control%time_step*mbs(ii)%angvel + dtdt*params%inv_i*m_torque(1:3,ii)
             ! rotate the molecule
             CALL rotate_molecule(mbs(ii),vec2q(spin))
             
             ! store the largest shift in a single coordinate
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog
             shift = control%time_step*mbs(ii)%vel
             mbs(ii)%pos = mbs(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(mbs(ii)%pos(jj) > cell(jj)) mbs(ii)%pos(jj) = mbs(ii)%pos(jj) - cell(jj)
                   IF(mbs(ii)%pos(jj) < 0.d0)     mbs(ii)%pos(jj) = mbs(ii)%pos(jj) + cell(jj)
                END IF
             END DO
             spin = control%time_step*mbs(ii)%angvel
             CALL rotate_molecule(mbs(ii),vec2q(spin))
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO
       
       DO ii = 1, n_ats ! loop atoms
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm

             ! shift: r+ = v*dt + a*dt**2/2
             shift = control%time_step*ats(ii)%vel + dtdt*params%invm_atoms(ats(ii)%type)*a_force(1:3,ii)
             ats(ii)%pos = ats(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3 ! loop x, y, z
                IF(pbc(jj))THEN ! periodic boundary
                   IF(ats(ii)%pos(jj) > cell(jj)) ats(ii)%pos(jj) = ats(ii)%pos(jj) - cell(jj)
                   IF(ats(ii)%pos(jj) < 0.d0)     ats(ii)%pos(jj) = ats(ii)%pos(jj) + cell(jj)
                END IF
             END DO
             ! record the max shift in a single coordinate
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog
             shift = control%time_step*ats(ii)%vel
             ats(ii)%pos = ats(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(ats(ii)%pos(jj) > cell(jj)) ats(ii)%pos(jj) = ats(ii)%pos(jj) - cell(jj)
                   IF(ats(ii)%pos(jj) < 0.d0)     ats(ii)%pos(jj) = ats(ii)%pos(jj) + cell(jj)
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO


       
       DO ii = 1, n_gos ! loop gos
          IF(control%md_algo == velocity_verlet_index)THEN ! velocity verlet algorithm

             ! shift: r+ = v*dt + a*dt**2/2
             shift = control%time_step*gos(ii)%vel + dtdt*g_force(1:3,ii)/gos(ii)%mass

!!$             write(*,*) ii, gos(ii)%vel
!!$             write(*,*) ii, g_force(1:3,ii)
!!$             write(*,*) ii, shift

             gos(ii)%pos = gos(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3 ! loop x, y, z
                IF(pbc(jj))THEN ! periodic boundary
                   IF(gos(ii)%pos(jj) > cell(jj)) gos(ii)%pos(jj) = gos(ii)%pos(jj) - cell(jj)
                   IF(gos(ii)%pos(jj) < 0.d0)     gos(ii)%pos(jj) = gos(ii)%pos(jj) + cell(jj)
                END IF
             END DO
             ! record the max shift in a single coordinate
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))

          ELSE IF(control%md_algo == leapfrog_index)THEN ! leapfrog
             shift = control%time_step*gos(ii)%vel
             gos(ii)%pos = gos(ii)%pos + shift
             ! periodicity
             DO jj = 1, 3
                IF(pbc(jj))THEN
                   IF(gos(ii)%pos(jj) > cell(jj)) gos(ii)%pos(jj) = gos(ii)%pos(jj) - cell(jj)
                   IF(gos(ii)%pos(jj) < 0.d0)     gos(ii)%pos(jj) = gos(ii)%pos(jj) + cell(jj)
                END IF
             END DO
             IF(MAXVAL(shift(:)) > max_shift) max_shift = MAXVAL(shift(:))
          END IF
       END DO


    END IF

    RETURN
  END SUBROUTINE move_particles

  ! Carries out one cg step. The step begins with a linear minimum energy search in
  ! the given direction, followed by determination of the next conjugate direction.  
  ! *mbs MB molecules
  ! *ats atomic particles
  ! *m_dummy dummy MBs
  ! *a_dummy dummy atoms
  ! *n_mbs number of MB molecules
  ! *n_ats number of atomic particles
  ! *n_elems number of different types of atoms
  ! *cell supercell dimensions (x,y,z)
  ! *pbc periodic boundary conditions (x,y,z)
  ! *btype boundary type index
  ! *bval boundary value
  ! *bonds bond counts for MB molecules
  ! *pdone list of MB pairs for bookkeeping
  ! *apdone list of atom pairs for bookkeeping
  ! *is_constr true if constraints have been defined, false otherwise
  ! *m_force forces acting on MB molecules
  ! *m_torque torques acting on MB molecules
  ! *a_force forces acting on atomic particles
  ! *om_force old forces acting on MB molecules
  ! *om_torque old torques acting on MB molecules
  ! *oa_force old forces acting on atomic particles
  ! *mc_force conjugate forces acting on MB molecules
  ! *mc_torque conjugate torques acting on MB molecules
  ! *ac_force conjugate forces acting on atomic particles
  ! *params physical parameters
  ! *control control parameters
  ! *mm_nbors MB-MB neighbors
  ! *ma_nbors MB-atom neighbors
  ! *aa_nbors atom-atom neighbors
  ! *mm_n_ns number of mb-mb neighbors
  ! *ma_n_ns number of mb-atom neighbors
  ! *aa_n_ns number of atom-atom neighbors
  ! *current_minimum to monitor energy convergence, the current minimum should be given, and will be returned
  ! *clock a stepcounter, will be increased by one
  ! *stop_now will be set to true if the convergence criterion is met
  ! *beta the beta parameter for conjugate gradients
  SUBROUTINE cg_step(mbs,ats,m_dummy,a_dummy,n_mbs,n_ats,cell,pbc,btype,bval,&
       params,control,bonds,is_constr,&
       mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,&
       boxed,boxes,box_mbs,box_ats,box_mc,box_ac,&
       m_force,m_torque,a_force,om_force,om_torque,oa_force,mc_force,mc_torque,ac_force,&
       current_minimum,clock,stop_now,beta,n_bond)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:), m_dummy(:)
    TYPE(atom), POINTER :: ats(:), a_dummy(:)
    TYPE(gop), POINTER :: gos(:), g_dummy(:)    
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:), g_force(:,:), &
         om_force(:,:), om_torque(:,:), oa_force(:,:), &
         mc_force(:,:), mc_torque(:,:), ac_force(:,:)
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3)
    REAL(KIND=dp), INTENT(INOUT) :: clock, current_minimum, beta
    INTEGER, INTENT(IN) :: btype(3)
    !LOGICAL, POINTER :: pdone(:,:), apdone(:,:)
    INTEGER, INTENT(IN) :: n_mbs, n_ats
    LOGICAL, INTENT(IN) :: is_constr, pbc(3)
    LOGICAL, INTENT(OUT) :: stop_now
    REAL(KIND=dp), POINTER :: bonds(:),n_bond(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_ns(:), &
         ma_nbors(:,:), ma_n_ns(:), aa_nbors(:,:), aa_n_ns(:), &
         mg_nbors(:,:), mg_n_ns(:), &
         box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), &
         box_mc(:,:,:,:), box_ac(:,:,:,:)
    INTEGER, INTENT(IN) :: boxes(3,2)
    LOGICAL, INTENT(IN) :: boxed
    INTEGER :: ii, jj
    REAL(KIND=dp) :: dot1, dot2, previous_minimum
    
    clock = clock + 1.
    stop_now = .false.

    previous_minimum = current_minimum
    ! line search in the direction of the conjugate forces/gradients
    CALL cg_linesearch(mbs,ats,m_dummy,a_dummy,n_mbs,n_ats,cell,pbc,btype,bval,&
       params,bonds,is_constr,&
       mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,&
       boxed,boxes,box_mbs,box_ats,box_mc,box_ac,&
       mc_force,mc_torque,ac_force,10,0.2d0,current_minimum,drift)

!!$    ! Didn't find anything?
!!$    IF(previous_minimum == current_minimum)THEN
!!$       CALL cg_linesearch(mbs,ats,m_dummy,a_dummy,n_mbs,n_ats,n_elems,cell,pbc,btype,bval,&
!!$            params,control,bonds,pdone,apdone,is_constr,&
!!$            mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,&
!!$            mc_force,mc_torque,ac_force,8,0.01d0,current_minimum,drift)
!!$    END IF

    ! Output
    IF(control%verbose > 0.d0)THEN
       WRITE(*,'(A,I8,A,ES20.4)')      "CG step ",FLOOR(clock+0.5d0),&
            ", change in energy:           ", previous_minimum - current_minimum
       WRITE(*,'(A)') " "
       WRITE(output,'(A,I8,A,ES20.4)') "CG step ",FLOOR(clock+0.5d0),&
           ", change in energy:           ", previous_minimum - current_minimum
       WRITE(output,'(A)') " "
    ELSE
       WRITE(output,'(A,I8,A,ES20.4)') "CG step ",FLOOR(clock+0.5d0),&
           ", change in energy:           ", previous_minimum - current_minimum
       WRITE(output,'(A)') " "
    END IF

    ! If we are close enough to the minimum, stop
    IF(control%time_step < 0.d0)THEN ! check for forces

       stop_now = .true.
       IF(previous_minimum > current_minimum)THEN

          checkmf: DO ii = 1, n_mbs
             DO jj = 1, 3
                IF(mbs(ii)%constrained(jj) /= all_frozen_index .AND. &
                     mbs(ii)%constrained(jj) /= frozen_pos_index )THEN
                   IF( m_force(jj,ii) > -control%time_step )THEN
                      stop_now = .false.
                      EXIT checkmf
                   END IF
                END IF
                IF(mbs(ii)%constrained(jj) /= all_frozen_index )THEN
                   IF( m_torque(jj,ii) > -control%time_step )THEN
                      stop_now = .false.
                      EXIT checkmf
                   END IF
                END IF
             END DO
          END DO checkmf
          IF(stop_now)THEN
             checkaf: DO ii = 1, n_ats
                DO jj = 1, 3
                   IF(mbs(ii)%constrained(jj) /= all_frozen_index .AND. &
                        mbs(ii)%constrained(jj) /= frozen_pos_index )THEN
                      IF( a_force(jj,ii) > -control%time_step**2 )THEN
                         stop_now = .false.
                         EXIT checkaf
                      END IF
                   END IF
                END DO
             END DO checkaf
          END IF

       END IF

    ELSE
       
       IF(previous_minimum - current_minimum < control%time_step)THEN ! check for energy difference

          stop_now = .true.

       END IF
    END IF

    ! end CG if we are close to a minimum
    IF(stop_now)THEN

       IF(control%verbose > 0.d0)THEN
          WRITE(*,'(A)') "CG search reached required accuracy "
          WRITE(*,*) ""
          WRITE(output,'(A)') "CG search reached required accuracy "
          WRITE(output,*) ""
       ELSE
          WRITE(output,'(A)') "CG search reached required accuracy "
          WRITE(output,*) ""
       END IF

       RETURN
    END IF

    om_force = m_force
    om_torque = m_torque
    oa_force = a_force

    ! calculate new forces/gradients
    CALL calc_forces(mbs,ats,gos,&
         mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
         bonds,cell,pbc,btype,bval,control,params,go,is_constr,.false.,&
         m_force,m_torque,a_force,g_force,virial,n_bond)

    ! "mixing parameter"
    dot1 = 0.d0
    dot2 = 0.d0
    beta = 0.d0
    DO ii = 1, n_mbs
       DO jj = 1, 3
          SELECT CASE(mbs(ii)%constrained(jj))
          CASE(all_frozen_index)
          CASE(frozen_pos_index)
             dot1 = dot1 + (m_torque(jj,ii)*(m_torque(jj,ii)-om_torque(jj,ii)))
             dot2 = dot2 + (om_torque(jj,ii)**2)
          CASE DEFAULT
             dot1 = dot1 + (m_force(jj,ii)*(m_force(jj,ii)-om_force(jj,ii)))
             dot2 = dot2 + (om_force(jj,ii)**2)
             dot1 = dot1 + (m_torque(jj,ii)*(m_torque(jj,ii)-om_torque(jj,ii)))
             dot2 = dot2 + (om_torque(jj,ii)**2)
          END SELECT
       END DO
    END DO
    DO ii = 1, n_ats
       DO jj = 1, 3
          SELECT CASE(ats(ii)%constrained(jj))
          CASE(all_frozen_index)
          CASE(frozen_pos_index)
          CASE DEFAULT
             dot1 = dot1 + (a_force(jj,ii)*(a_force(jj,ii)-oa_force(jj,ii)))
             dot2 = dot2 + (oa_force(jj,ii)**2)
          END SELECT
       END DO
    END DO
    IF(ABS(dot2) > denom_tolerance)THEN
       beta = MAX(0.d0,dot1/dot2)
       IF(beta < norm_tolerance) beta = 0.d0
    ELSE
       beta = 0.d0
    END IF
    
    ! new conjugate force/gradient
    mc_force = m_force + beta*mc_force
    mc_torque = m_torque + beta*mc_torque
    ac_force = a_force + beta*ac_force

    CALL checkpoint(time_force,timeA,timeB)

    RETURN
  END SUBROUTINE cg_step

  ! Conjugate gradient linesearch: finds the minimum energy in the given direction.  
  ! *mbs MB molecules
  ! *ats atomic particles
  ! *m_dummy dummy MBs
  ! *a_dummy dummy atoms
  ! *n_mbs number of MB molecules
  ! *n_ats number of atomic particles
  ! *n_elems number of different types of atoms
  ! *cell supercell dimensions (x,y,z)
  ! *pbc periodic boundary conditions (x,y,z)
  ! *btype boundary type index
  ! *bval boundary value
  ! *bonds bond counts for MB molecules
  ! *pdone list of MB pairs for bookkeeping
  ! *apdone list of atom pairs for bookkeeping
  ! *is_constr true if constraints have been defined, false otherwise
  ! *m_force conjugate forces acting on MB molecules
  ! *m_torque conjugate torques acting on MB molecules
  ! *a_force conjugate forces acting on atomic particles
  ! *params physical parameters
  ! *control control parameters
  ! *mm_nbors MB-MB neighbors
  ! *ma_nbors MB-atom neighbors
  ! *aa_nbors atom-atom neighbors
  ! *mm_n_ns number of mb-mb neighbors
  ! *ma_n_ns number of mb-atom neighbors
  ! *aa_n_ns number of atom-atom neighbors
  ! *n_scales the search reduces the steplength once a minimum is found and repeats the search. this  controls the number of levels
  ! *ini_step the step length of the first search
  ! *best_energy the best energy found
  ! *tot_drift measures the maximum length by which the particles have moved, for neighbors list recalculation
  SUBROUTINE cg_linesearch(mbs,ats,m_dummy,a_dummy,n_mbs,n_ats,cell,pbc,btype,bval,&
       params,bonds,is_constr,&
       mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,boxed,boxes,box_mbs,box_ats,box_mc,box_ac,&
       m_force,m_torque,a_force,n_scales,ini_step,best_energy,tot_drift)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:), m_dummy(:)
    TYPE(atom), POINTER :: ats(:), a_dummy(:)
    TYPE(gop), POINTER :: gos(:), g_dummy(:)
    TYPE(mbps), INTENT(IN) :: params    
    TYPE(gops) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:), g_force(:,:)
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3), ini_step
    INTEGER, INTENT(IN) :: btype(3)
    !LOGICAL, POINTER :: pdone(:,:), apdone(:,:)
    INTEGER, INTENT(IN) :: n_mbs, n_ats, n_scales
    LOGICAL, INTENT(IN) :: is_constr, pbc(3)
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_ns(:), &
         ma_nbors(:,:), ma_n_ns(:), aa_nbors(:,:), aa_n_ns(:), &
         mg_nbors(:,:), mg_n_ns(:), &
         box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), box_gos(:,:,:,:,:), &
         box_mc(:,:,:,:), box_ac(:,:,:,:), box_gc(:,:,:,:)
    INTEGER, INTENT(IN) :: boxes(3,2)
    LOGICAL, INTENT(IN) :: boxed
    INTEGER :: ii, jj, i_scale, kk, steps
    !REAL(KIND=dp), INTENT(INOUT) :: best_energy
    REAL(KIND=dp), INTENT(INOUT) :: tot_drift, best_energy
    REAL(KIND=dp) :: dt, spin(3), shift(3), last_energy, energy, &
         ene1,ene2,ene3,ene4,ene5,max_grad,tot_shift

!    CALL pot_energy(mbs,ats,mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,&
!         bonds,cell,pbc,btype,bval,params,pdone,apdone,is_constr,&
!         energy,ene1,ene2,ene3,ene4,ene5)    
    energy = best_energy
    last_energy = energy
    dt = ini_step
    i_scale = 1
    kk = 0
    tot_shift = 0.d0

    CALL checkpoint(time_stat,timeA,timeB)
    
    max_grad = 0.d0
    DO ii = 1, n_mbs
       IF( (m_force(1:3,ii).o.m_force(1:3,ii)) > max_grad**2 )THEN
          max_grad = sqrt(m_force(1:3,ii).o.m_force(1:3,ii))
       END IF
    END DO
    DO ii = 1, n_ats
       IF( (a_force(1:3,ii).o.a_force(1:3,ii)) > max_grad**2 )THEN
          max_grad = sqrt(a_force(1:3,ii).o.a_force(1:3,ii))
       END IF
    END DO

    DO WHILE(dt*max_grad > params%cut_ver*0.5d0)
       dt = dt*0.5d0
    END DO

    CALL checkpoint(time_stat,timeA,timeB)
    steps = 0

    DO WHILE(i_scale <= n_scales)

       kk = kk + 1
       IF(kk > 10000)THEN
          CALL abort("conjugate gradient linesearch stuck")
       END IF

       DO ii = 1, n_mbs
          DO jj = 1, 3
             SELECT CASE (mbs(ii)%constrained(jj)) 
             CASE(frozen_pos_index)
                shift(jj) = 0.d0
                spin(jj) = dt*m_torque(jj,ii)
             CASE(all_frozen_index)
                shift(jj) = 0.d0
                spin(jj) = 0.d0
             CASE DEFAULT
                shift(jj) = dt*m_force(jj,ii)
                spin(jj) = dt*m_torque(jj,ii)
             END SELECT
             m_dummy(ii)%pos(jj) = m_dummy(ii)%pos(jj) + shift(jj)
             ! periodicity
             IF(pbc(jj))THEN
                IF(m_dummy(ii)%pos(jj) > cell(jj))THEN
                   m_dummy(ii)%pos(jj) = m_dummy(ii)%pos(jj) - cell(jj)
                   m_dummy(ii)%inipos(jj) = m_dummy(ii)%inipos(jj) - cell(jj)
                END IF
                IF(m_dummy(ii)%pos(jj) < 0.d0)THEN
                   m_dummy(ii)%pos(jj) = m_dummy(ii)%pos(jj) + cell(jj)
                   m_dummy(ii)%inipos(jj) = m_dummy(ii)%inipos(jj) + cell(jj)
                END IF
             END IF
          END DO
          CALL rotate_molecule(m_dummy(ii),vec2q(spin))
          CALL norm_quaternion(m_dummy(ii)%orientation)
       END DO
       DO ii = 1, n_ats
          DO jj = 1, 3
             SELECT CASE (ats(ii)%constrained(jj)) 
             CASE(frozen_pos_index)
                shift(jj) = 0.d0
             CASE(all_frozen_index)
                shift(jj) = 0.d0
             CASE DEFAULT
                shift(jj) = dt*a_force(ii,jj)
             END SELECT
             a_dummy(ii)%pos(jj) = a_dummy(ii)%pos(jj) + shift(jj)
             ! periodicity
             IF(pbc(jj))THEN
                IF(a_dummy(ii)%pos(jj) > cell(jj))THEN
                   a_dummy(ii)%pos(jj) = a_dummy(ii)%pos(jj) - cell(jj)
                   a_dummy(ii)%inipos(jj) = a_dummy(ii)%inipos(jj) - cell(jj)
                END IF
                IF(a_dummy(ii)%pos(jj) < 0.d0)THEN
                   a_dummy(ii)%pos(jj) = a_dummy(ii)%pos(jj) + cell(jj)
                   a_dummy(ii)%inipos(jj) = a_dummy(ii)%inipos(jj) + cell(jj)
                END IF
             END IF
          END DO
       END DO
       steps = steps + 1

       CALL checkpoint(time_move,timeA,timeB)

       ! check for neighbor updating need
       tot_drift = tot_drift + dt*max_grad
       tot_shift = tot_shift + dt
       IF(tot_drift > 0.5d0*params%cut_ver)THEN
          tot_drift = 0.d0
          CALL update_neighbors(m_dummy,a_dummy,g_dummy,cell,pbc,bonds,&
               boxed,boxes,box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
               mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
               params,go)
       ELSE
          CALL update_bond_numbers(m_dummy,cell,pbc,mm_nbors,mm_n_ns,bonds,params%R_b,params%D_b,params%inv_Db)
       END IF

       CALL checkpoint(time_nbor,timeA,timeB)

       ! check energy
       last_energy = energy

       IF(tot_shift == 0.d0)THEN
          ! we are on the spot where the best energy should be,
          ! but due to numerical error, the positions may be slightly off.
          ! -> correct this
          energy = best_energy
          DO jj = 1, n_mbs
             m_dummy(jj) = mbs(jj)
          END DO
          DO jj = 1, n_ats
             a_dummy(jj) = ats(jj)
          END DO

       ELSE

!!$          CALL pot_energy(m_dummy,a_dummy,g_dummy,&
!!$               mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
!!$               bonds,cell,pbc,btype,bval,params,is_constr,&
!!$               energy,ene1,ene2,ene3,ene4,ene5,ene6,ene7,ene8,ene9,ene10,ene11,ene12)

       END IF

       IF(energy < best_energy)THEN
          ! record the configuration with the best energy
          best_energy = energy
          DO jj = 1, n_mbs
             mbs(jj) = m_dummy(jj)
          END DO
          DO jj = 1, n_ats
             ats(jj) = a_dummy(jj)
          END DO
          tot_shift = 0.d0
       END IF
       
       ! debug
       !WRITE(*,*) kk, i_scale, tot_shift, steps, energy, last_energy
       
       ! if we start to climb again, turn around and take smaller steps
       IF(energy > last_energy .AND. (steps > 1 .OR. i_scale == 1))THEN
          steps = 0
          IF(i_scale == n_scales)THEN ! this was the final level, just end the loop
             i_scale = i_scale + 1
             ! copy the best configuration to dummies
             DO jj = 1, n_mbs
                m_dummy(jj) = mbs(jj)
             END DO
             DO jj = 1, n_ats
                a_dummy(jj) = ats(jj)
             END DO
             ! update neighbors according to the real particles
             CALL update_neighbors(mbs,ats,gos,cell,pbc,bonds,&
                  boxed,boxes,box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
                  mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
                  params,go)
          ELSE
             dt = -0.5d0*dt
             i_scale = i_scale + 1
          END IF
       END IF
       
       CALL checkpoint(time_stat,timeA,timeB)
       
    END DO

    RETURN
  END SUBROUTINE cg_linesearch

  ! Scales the simulation box and coordinates linearly
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *ratio ratio of new and old volumes
  ! *scalable logic list, marks if x-, y- and z-directions can be scaled
  SUBROUTINE scale_box(mbs,ats,cell,ratio,scalable)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(INOUT) :: cell(3)
    REAL(KIND=dp), INTENT(IN) :: ratio
    LOGICAL, INTENT(IN) :: scalable(3)
    REAL(KIND=dp) :: scaler
    INTEGER :: ii, jj, n_mbs, n_ats, scaledim

    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    ! check how many dimensions are scalable
    scaledim = 0
    DO ii = 1, 3
       IF(scalable(ii)) scaledim = scaledim + 1
    END DO
    IF(scaledim == 0) RETURN

    ! calculate the linear scaling ratio
    scaler = ratio**(1.d0/REAL(scaledim,KIND=dp))
    
    DO jj = 1, 3 ! loop x, y, z
       ! scale dimensions of the supercell
       IF(scalable(jj))THEN
          cell(jj) = cell(jj)*scaler
          ! shift particles according the scaling
          DO ii = 1, n_mbs
             mbs(ii)%pos(jj) = mbs(ii)%pos(jj)*scaler
          END DO
          DO ii = 1, n_ats
             ats(ii)%pos(jj) = ats(ii)%pos(jj)*scaler
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE scale_box

  ! Applies a barostat, or calculates the average pressure for future barostat invoke.
  ! This piece of code used to be in the main MD loop, but was moved here to make the
  ! code easier to read. It still uses only global variables.
  SUBROUTINE apply_barostat()
    IMPLICIT NONE

    ! Berendsen barostat
    IF(control_params%md_baro == berendsen_index)THEN
       ! calculate average pressure over baro_interval steps
       CALL kin_energy(molecules,atoms,gos,n_freedom,e_kin,temps,e_lin,e_rot)
       press = press + (0.66666666667d0*e_lin + virial)/(supercell(1)*supercell(2)*supercell(3)) / &
            barosteps
       
       ! apply the barostat only after barotimer has been reached
       IF(simulation_time >= barotimer)THEN

          ! cpu time checkpoint
          CALL checkpoint(time_stat,timeA,timeB)

          ! supercell volume scaling ratio
          scaleratio = (1.d0 + control_params%baro_interval*control_params%time_step/control_params%baro_value*&
               (press - control_params%ini_press))
          ! limit the scaling to a region
          IF(1.d0-scaleratio < -scale_max)THEN
             scaleratio = 1.d0 + scale_max 
          ELSE IF(1.d0-scaleratio > scale_max)THEN
             scaleratio = 1.d0 - scale_max 
          END IF
          ! scale the supercell
          CALL scale_box(molecules,atoms,supercell, &
               scaleratio, &
               control_params%scalable)

          ! the scaling may shift particles
          drift = drift + ABS(1.d0-scaleratio)*MAXVAL(supercell)
          IF(drift > drift_range)THEN ! update neighbor lists if the particles have moved enough
             drift = 0.d0
             CALL update_neighbors(molecules,atoms,gos,supercell,periodic_boundary,mb_bonds,&
                  boxed,boxes,box_mbs,box_ats,box_gos,box_mb_count,box_at_count,box_go_count,&
                  mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,&
                  atat_neighbors,atat_n_nbors,mbgo_neighbors,mbgo_n_nbors,&
                  physical_params,go_params)
          ELSE
             ! bonds must be updated always
             CALL update_bond_numbers(molecules,supercell,periodic_boundary,&
                  mbmb_neighbors,mbmb_n_nbors,mb_bonds,&
                  physical_params%R_b,physical_params%D_b,physical_params%inv_Db)
          END IF
          
          CALL checkpoint(time_move,timeA,timeB)
          
          ! calculate forces according to the new configuration
          CALL calc_forces(molecules,atoms,gos,&
               mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,&
               atat_neighbors,atat_n_nbors,mbgo_neighbors,mbgo_n_nbors,&
               mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
               control_params,physical_params,go_params,&
               is_constrained,.false.,&
               mb_forces,mb_torques,atom_forces,go_forces,virial,n_bond)
          
          CALL checkpoint(time_force,timeA,timeB)
          
          ! reset counters
          barotimer = barotimer + control_params%baro_interval*control_params%time_step
          barosteps = FLOOR((barotimer-simulation_time)/control_params%time_step)+1
          press = 0.d0
       END IF
    END IF ! berendsen barostat
    
  END SUBROUTINE apply_barostat

  ! Checks if it is time to write something out mid-run.
  ! This piece of code used to be in the main MD loop, but was moved here to make the
  ! code easier to read. It still uses only global variables.
  SUBROUTINE interval_outputs()
    IMPLICIT NONE

     ! write rdf after given intervals
     IF(control_params%rdf_writer == ixyz_index)THEN
        IF(simulation_time > rdftimer .OR. finish)THEN
           rdftimer = rdftimer+control_params%rdf_interval
           CALL write_rdf(system_name(1:namelength-1)//RDF_MID//".out",&
                molecules,atoms,supercell,periodic_boundary,&
                mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
                control_params,.false.)

           CALL checkpoint(time_stat,timeA,timeB)

        END IF
     ELSE IF(control_params%rdf_writer == average_index)THEN
        IF(simulation_time > rdftimer .OR. finish)THEN
           rdftimer = rdftimer+control_params%rdf_interval
           CALL checkpoint(time_io,timeA,timeB)
           CALL radial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,&
                mbmb_neighbors,mbmb_n_nbors,mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,physical_params,&
                control_params%rdf_particles,control_params%rdf_range,rdf_temp)
           CALL checkpoint(time_stat,timeA,timeB)
           rdf_ave = rdf_ave + rdf_temp
           rdf_sq = rdf_sq + rdf_temp*rdf_temp
           rdf_count = rdf_count+1           
        END IF
     END IF

     ! write adf after given intervals
     IF(control_params%adf_writer == ixyz_index)THEN
        IF(simulation_time > adftimer .OR. finish)THEN
           adftimer = adftimer+control_params%adf_interval
           CALL write_adf(system_name(1:namelength-1)//ADF_MID//".out",&
                molecules,atoms,supercell,periodic_boundary,physical_params,&
                control_params,.false.)

           CALL checkpoint(time_stat,timeA,timeB)

        END IF
     ELSE IF(control_params%adf_writer == average_index)THEN
        IF(simulation_time > adftimer .OR. finish)THEN
           adftimer = adftimer+control_params%adf_interval
           CALL checkpoint(time_io,timeA,timeB)
           CALL axial_distribution_function(molecules,atoms,n_mols,n_ats,supercell,periodic_boundary,physical_params,&
                control_params,adf_temp)
           CALL checkpoint(time_stat,timeA,timeB)
           adf_ave = adf_ave + adf_temp
           adf_sq = adf_sq + adf_temp*adf_temp
           adf_count = adf_count+1           
        END IF
     END IF

     ! write xyz after given intervals
     IF(control_params%xyz_writer == ixyz_index)THEN
        IF(simulation_time > xyztimer .OR. finish)THEN
           xyztimer = xyztimer+control_params%xyz_interval
           CALL write_xyz(system_name(1:namelength-1)//".xyz",&
                molecules,atoms,gos,physical_params%l_oh,.false.)        
        END IF
     END IF

     ! write bonds after given intervals
     IF(control_params%bond_writer == ixyz_index)THEN
        IF(simulation_time > bondtimer .OR. finish)THEN
           bondtimer = bondtimer+control_params%bond_interval
           CALL write_bonds(system_name(1:namelength-1)//MB_MID//".xyz",system_name(1:namelength-1)//BOND_MID//".out",n_bond,&
                molecules,mbmb_neighbors,mbmb_n_nbors,supercell,periodic_boundary,control_params,physical_params,.false.)
        END IF
     END IF

     ! write continuation file
     IF(control_params%inp_writer == ixyz_index)THEN
        IF(simulation_time > inptimer .OR. finish)THEN
           inptimer = inptimer+control_params%inp_interval
           timelabel = real_to_string(simulation_time)
           CALL write_contfile(system_name(1:namelength-1)//CONT_MID//&
                timelabel(11-MIN(10,MAX(3,1+INT(log(control_params%time_max)/log(10.d0)+0.01))):10 )//"."//MB_IN,&
                molecules,atoms,gos,n_elements,&
                physical_params,control_params,go_params,&
                supercell,boundary_type,boundary_value,is_constrained)
        END IF
     END IF

     CALL checkpoint(time_io,timeA,timeB)

     updated_stats = .false.

     ! Write current situation
     IF(simulation_time > speaktimer .OR. finish)THEN

        ! avoids double calculation of energy etc. if stats are also recorded
        updated_stats = .true.

        ! wall clock
        CALL DATE_AND_TIME(values=newt)
        CALL time_difference(tvector,newt,timediff)
        
        ! Calculate energy
        CALL kin_energy(molecules,atoms,gos,n_freedom,e_kin,temps,e_lin,e_rot)
        CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
             mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
             mbgo_neighbors,mbgo_n_nbors, &
             mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
             physical_params,go_params,is_constrained,&
             e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
             e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)       

        CALL checkpoint(time_stat,timeA,timeB)


        IF(cpu_id == master_cpu)THEN
           IF(control_params%verbose > 0.d0)THEN
              ! virtual time
              IF(control_params%md_algo == cg_index)THEN
                 WRITE(*,'(A,I25)')   "CG steps:              ", FLOOR(simulation_time+0.5d0)
              ELSE
                 WRITE(*,'(A,F25.4)') "Elapsed virtual time:  ", simulation_time
              END IF
              ! time_difference doesn't handle months... maybe fix it?
              IF(timediff(6) /= 0)THEN
                 clocklabel = " - "
                 WRITE(output,'(A)') &
                      "Wall clock time:              "//clocklabel
              ELSE
                 CALL time_string(timediff(1),timediff(2),timediff(3),&
                      timediff(4),timediff(5),clocklabel)
                 WRITE(*,'(A)') &
                      "Wall clock time:              "//clocklabel
              END IF
              WRITE(*,'(A,F20.4)') "Temperature:                ", temps
              WRITE(*,'(A,F20.4)') "Total energy:               ", e_kin+e_pot
              WRITE(*,'(A,F20.4)') "  kinetic:                  ", e_kin
              WRITE(*,'(A,F20.4)') "  potential:                ", e_pot
              WRITE(*,'(A,F20.4)') "     MB-MB Lennard-Jones:   ", e_lenjon
              WRITE(*,'(A,F20.4)') "     MB-MB hydrogen bond:   ", e_bonds
              IF(n_ats > 0)THEN
                 WRITE(*,'(A,F20.4)') "     MB-atom potential:     ", e_mbat
                 WRITE(*,'(A,F20.4)') "     atom-atom potential:   ", e_atat
              END IF
              IF(n_gos > 0)THEN
                 WRITE(*,'(A,F20.4)') "     MB-GO potential:       ", e_mbgo
                 WRITE(*,'(A,F20.4)') "     GO-GO potential:       ", e_go
                 WRITE(*,'(A,F20.4)') "        GO stretching:      ", e_stretch
                 WRITE(*,'(A,F20.4)') "        GO bending:         ", e_bend
                 WRITE(*,'(A,F20.4)') "        GO torsion:         ", e_torsion
                 WRITE(*,'(A,F20.4)') "        GO native:          ", e_native
                 WRITE(*,'(A,F20.4)') "        GO non-native:      ", e_nonnat
              END IF
              IF(is_constrained)THEN
                 WRITE(*,'(A,F20.4)') "     constraint potential:  ", e_constr
              END IF
              WRITE(*,*) ""
           END IF
           ! To file
           
           IF(control_params%md_algo == cg_index)THEN
              WRITE(output,'(A,I25)')   "CG steps:              ", FLOOR(simulation_time+0.5d0)
           ELSE
              WRITE(output,'(A,F25.4)') "Elapsed virtual time:  ", simulation_time
           END IF
           IF(timediff(6) /= 0)THEN
              clocklabel = " - "
              WRITE(output,'(A)') &
                   "Wall clock time:              "//clocklabel
           ELSE
              CALL time_string(timediff(1),timediff(2),timediff(3),&
                   timediff(4),timediff(5),clocklabel)
              WRITE(output,'(A)') &
                   "Wall clock time:              "//clocklabel
           END IF
           WRITE(output,'(A,F20.4)') "Temperature:                ", temps
           WRITE(output,'(A,F20.4)') "Total energy:               ", e_kin+e_pot
           WRITE(output,'(A,F20.4)') "  kinetic:                  ", e_kin
           WRITE(output,'(A,F20.4)') "  potential:                ", e_pot
           WRITE(output,'(A,F20.4)') "     MB-MB Lennard-Jones:   ", e_lenjon
           WRITE(output,'(A,F20.4)') "     MB-MB hydrogen bond:   ", e_bonds
           IF(n_ats > 0)THEN
              WRITE(output,'(A,F20.4)') "     MB-atom potential:     ", e_mbat
              WRITE(output,'(A,F20.4)') "     atom-atom potential:   ", e_atat
           END IF
           IF(is_constrained)THEN
              WRITE(output,'(A,F20.4)') "     constraint potential:  ", e_constr
           END IF
           
           WRITE(output,*) ""
           !CALL flush(output)
        END IF

        speaktimer = speaktimer + ABS(control_params%verbose)

        CALL checkpoint(time_io,timeA,timeB)

     END IF

     ! Write stats
     DO k = 1, control_params%n_statfiles
        IF(simulation_time > stattimer(k) .OR. finish)THEN

           IF(.NOT.updated_stats)THEN
              
              ! wall clock
              CALL DATE_AND_TIME(values=newt)
              CALL time_difference(tvector,newt,timediff)
              
              ! Calculate energy
              CALL kin_energy(molecules,atoms,gos,n_freedom,e_kin,temps,e_lin,e_rot)
              CALL pot_energy(molecules,atoms,gos,mbmb_neighbors,mbmb_n_nbors,&
                   mbat_neighbors,mbat_n_nbors,atat_neighbors,atat_n_nbors,&
                   mbgo_neighbors,mbgo_n_nbors, &
                   mb_bonds,supercell,periodic_boundary,boundary_type,boundary_value,&
                   physical_params,go_params,is_constrained,&
                   e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,&
                   e_mbgo, e_go, e_stretch, e_bend, e_torsion, e_native, e_nonnat)       
              updated_stats = .true.      
              
              CALL checkpoint(time_stat,timeA,timeB)
              
              IF(timediff(6) /= 0)THEN
                 clocklabel = " - "
              ELSE
                 CALL time_string(timediff(1),timediff(2),timediff(3),&
                      timediff(4),timediff(5),clocklabel)
              END IF
              
           END IF
           
           CALL write_stats(statout,molecules,atoms,temps,e_kin,e_pot,&
                e_lin,e_rot,e_lenjon,e_bonds,e_mbat,e_atat,virial,supercell,simulation_time,clocklabel,&
                mb_forces,mb_torques,atom_forces,control_params,k)
           
           stattimer(k) = stattimer(k) + control_params%stat_interval(k)
        END IF
     END DO

     CALL checkpoint(time_io,timeA,timeB)

   END SUBROUTINE interval_outputs

END PROGRAM cashew
