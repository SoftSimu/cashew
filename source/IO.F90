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
! file_handler is a module for handling input and output files. 
! It includes some generally useful file reading   
! subroutines as well as parsing subroutines specific to cashew.
! all the keywords used in the input file are defined here as
! parameters so that one can easily change the keywords.
! <br /><br />Back to <a href="cashew.html">cashew</a>

MODULE file_handler

  USE mb_model
  USE quaternions
  USE parameters
  USE mpi_mod

  ! File extensions
  CHARACTER(LEN=3), PARAMETER :: &
       MAIN_OUT = "out", &
       MB_IN    = "mb"
  CHARACTER(LEN=6), PARAMETER :: CONT_START = "_start"
  CHARACTER(LEN=4), PARAMETER :: CONT_END = "_end"
  CHARACTER(LEN=4), PARAMETER :: CONT_MID = "_mid"
  CHARACTER(LEN=6), PARAMETER :: STATS = "_stats"
  CHARACTER(LEN=4), PARAMETER :: RDF_MID = "_rdf"
  CHARACTER(LEN=4), PARAMETER :: ADF_MID = "_adf"
  CHARACTER(LEN=6), PARAMETER :: BOND_MID = "_bonds"
  CHARACTER(LEN=3), PARAMETER :: MB_MID = "_mb"

  CHARACTER(LEN=70), PARAMETER :: &
       HORILINE = "----------------------------------------------------------------------"

  ! Comment
  CHARACTER, PARAMETER :: comment_char = "!"

  ! Define input tags
  CHARACTER(LEN=labelw), PARAMETER :: &
       mbname = "mb", &
       armname = "hb", &
       grouptag = "to"
       
  CHARACTER(LEN=2), PARAMETER :: &
       lenjon_pot = "lj", &
       cgtag = "cg"

  CHARACTER(LEN=3), PARAMETER :: &
       E_LJtag = "elj", &
       E_HBtag = "ehb", &
       E_PHtag = "efi", &
       R_HBtag = "rhb", &
       L_ohtag = "loh", &
       M_moletag = "mmb", &
       I_moletag = "imb", &
       ENDtag = "end", &
       exp_pot = "exp", &
       pow_pot = "pow", &
       spring_pot = "spr", &
       alltag = "all", &
       postag = "pos", &
       veltag = "vel", &
       runtag = "run", &
       xyzread = "xyz"

  CHARACTER(LEN=4), PARAMETER :: &
!       Temptag = "temp", & ! changed to "temperature"
       Freetag = "free", &
       welltag = "well", &
       walltag = "wall", &
       algotag = "algo", &
       nonetag = "none", &
       testtag = "test", &
       fulltag = "full", &
       atomtag = "atom", &
       timetag = "time", &
       barotag = "baro", &
       armstag = "arms", &
       truetag = "true", &
       hard_pot = "hard", &
       fene_pot = "fene"

  CHARACTER(LEN=5), PARAMETER :: &
       R_btag = "rbond", &
       D_btag = "dbond", &
       Cut_HBtag = "cuthb", &
       Cut_LJtag = "cutlj", &
       Tsimutag = "tsimu", &
       Tsteptag = "tstep", &
       forcetag = "force", &
       starttag = "start", &
       hbondtag = "hbond", &
       cgtoltag = "cgtol", &
       unitstag = "units", &
       fixedtag = "fixed", &
       xaxistag = "xaxis", &
       yaxistag = "yaxis", &
       zaxistag = "zaxis", &
       angletag = "angle", &
       shell_pot = "shell"

  CHARACTER(LEN=6), PARAMETER :: &
       frozentag = "frozen", &
       thermotag = "thermo", &
       verlettag = "verlet", &
       coolertag = "cooler", &
       energytag = "energy", &
       virialtag = "virial", &
       lineartag = "linear", &
       lenjontag = "lenjon", &
       mbatomtag = "mbatom", &
       volumetag = "volume", &
       torquetag = "torque", &
       followtag = "follow"

  CHARACTER(LEN=7), PARAMETER :: &
       Sigma_LJtag = "sigmalj", &
       Sigma_THtag = "sigmath", &
       Sigma_HBtag = "sigmahb", &
       Exp_btag = "expbond", &
       Cut_Vtag = "cutverl", &
       Cut_Atag = "cutatom", &
       seedtag = "rndseed", &
       verbosetag = "verbose", &
       kinetictag = "kinetic", &
       angulartag = "angular", &
       averagetag = "average", &
       cgstepstag = "cgsteps", &
       reducedtag = "reduced", &
       noangletag = "noangle", &
       hardrep_pot = "hardrep"

  CHARACTER(LEN=8), PARAMETER :: &
       periodictag = "periodic", &
       leapfrogtag = "leapfrog", &
       langevintag = "langevin", &
       sendtag     = "startend", &
       intervaltag = "interval", &
       xyztag      = "writexyz", &
       bondtag     = "writebnd", &
       conttag     = "writeinp", &
       rdftag      = "writerdf", &
       adftag      = "writeadf", &
       velotag     = "velocity", &
       coordtag    = "position", &
       atomatomtag = "atomatom", &
       pressuretag = "pressure", &
       forcesumtag = "forcesum"

  CHARACTER(LEN=9), PARAMETER :: &
       clocktag = "wallclock", &
       berendsentag = "berendsen", &
       potentialtag = "potential", &
       writestattag = "writestat"

  CHARACTER(LEN=10), PARAMETER :: &
       rottag = "rotational"

  CHARACTER(LEN=11), PARAMETER :: &
       temptag = "temperature", &
       orientag = "orientation"

  ! input block tagnames
  CHARACTER(LEN=4),  PARAMETER :: cell_block   = "cell"
  CHARACTER(LEN=8),  PARAMETER :: ele_block    = "elements"
  CHARACTER(LEN=10), PARAMETER :: pot_block    = "potentials"
  CHARACTER(LEN=9),  PARAMETER :: pos_block    = "particles"
  CHARACTER(LEN=10), PARAMETER :: vel_block    = "velocities"
  CHARACTER(LEN=11), PARAMETER :: constr_block = "constraints"
  CHARACTER(LEN=7),  PARAMETER :: main_block   = "control"
  CHARACTER(LEN=8),  PARAMETER :: mb_block     = "mb-model"
  CHARACTER(LEN=10), PARAMETER :: stat_block   = "statistics"
  INTEGER, PARAMETER :: total_blocks = 9
  
  ! a weird real used for denoting missing input parameters
  REAL(KIND=dp), PARAMETER :: missing_r = -93678.39285d0

CONTAINS

  ! Reads and parses the main input file.
  ! This subroutine is meant to act as an interface between
  ! the main program and the more specific reading and parsing routines.
  ! Basically, the routine first reads and trims the whole input file
  ! as a character array using "read_lines" and "trim". 
  ! Then, it searches for input blocks and calls the routines that 
  ! parse the information in each block.
  ! Some validity checks are done in this routine,
  ! some in the other parsing routines called from this one.
  ! *filename the name of the input file - the string must be of the right length
  ! *control control parameters
  ! *params physical parameters
  ! *cell supercell dimensions
  ! *pbc true for periodic boundaries
  ! *btype boundary type index
  ! *bval boundary parameter
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mbs number of molecules
  ! *n_ats number of atoms
  ! *n_elems number of different types of atoms
  ! *n_freedom number of degrees of freedom
  ! *found_types the number of atoms of each type
  ! *is_constr true if there are constraints in place
  SUBROUTINE parse_input(filename,control,params,cell,pbc,btype,bval,&
       mbs,ats,n_mbs,n_ats,n_elems,n_freedom,found_types,is_constr)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(cps), INTENT(OUT) :: control
    TYPE(mbps), INTENT(OUT) :: params
    REAL(KIND=dp), INTENT(OUT) :: cell(3), bval(3)
    LOGICAL, INTENT(OUT) :: pbc(3), is_constr
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    INTEGER :: ii, jj, n_lines
    INTEGER, INTENT(OUT) :: n_elems, n_mbs, n_ats, n_freedom, btype(3)
    CHARACTER, POINTER :: data(:,:)
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, POINTER :: lwidths(:), twidths(:), found_types(:)
    CHARACTER(LEN=100) :: readtag
    CHARACTER(LEN=15) :: openname
    CHARACTER(LEN=labelw) :: label1, label2
    INTEGER :: blocklines(total_blocks,2), blnr, blockopen
    INTEGER :: multi(3)
    LOGICAL, POINTER :: random_mb_vel(:,:), random_at_vel(:)
    REAL(KIND=dp) :: adfmax

    ! read the entire input file as a character array
    CALL read_lines(filename,data,lwidths)
    ! trim the data
    CALL trim(data,lwidths,comment_char)
    ! lowercase everything so that the input is not case-sensitive
    CALL lowercase(data)

    n_lines = SIZE(data(1,:))
    blocklines = 0
    blockopen = 0
    openname = ""

    ! search for blocks
    DO ii = 1, n_lines ! loop over data lines
       IF(data(1,ii) == "<")THEN ! a tag?
          CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
          readtag = tokens(1)

          ! check for known tags and store the line numbers for starting and ending the
          ! input blocks in the array "blocklines"
          ! "openname" is the name of the block currently open, and "blockopen" is the corresponding index
          SELECT CASE ( readtag(1:twidths(1)) )
          CASE ("<"//main_block//">") ! start main
             blocklines(1,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 1
             openname = "<"//main_block//">"
          CASE ("</"//main_block//">") ! end main
             blocklines(1,2) = ii
             IF(blockopen /= 1) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//mb_block//">") ! start mb
             blocklines(2,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 2
             openname = "<"//mb_block//">"
          CASE ("</"//mb_block//">") ! end mb 
             blocklines(2,2) = ii
             IF(blockopen /= 2) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//cell_block//">") ! start cell
             blocklines(3,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 3
             openname = "<"//cell_block//">"
          CASE ("</"//cell_block//">") ! end cell
             blocklines(3,2) = ii
             IF(blockopen /= 3) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//ele_block//">") ! start elements
             blocklines(4,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 4
             openname = "<"//ele_block//">"
          CASE ("</"//ele_block//">") ! end elements
             blocklines(4,2) = ii
             IF(blockopen /= 4) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//pos_block//">") ! start positions
             blocklines(5,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 5
             openname = "<"//pos_block//">"
          CASE ("</"//pos_block//">") ! end positions
             blocklines(5,2) = ii
             IF(blockopen /= 5) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//vel_block//">") ! start velocities
             blocklines(6,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 6
             openname = "<"//vel_block//">"
          CASE ("</"//vel_block//">") ! end velocities
             blocklines(6,2) = ii
             IF(blockopen /= 6) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//constr_block//">") ! start constraints
             blocklines(7,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 7
             openname = "<"//constr_block//">"
          CASE ("</"//constr_block//">") ! end constraints
             blocklines(7,2) = ii
             IF(blockopen /= 7) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//pot_block//">") ! start potentials
             blocklines(8,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 8
             openname = "<"//pot_block//">"
          CASE ("</"//pot_block//">") ! end potentials
             blocklines(8,2) = ii
             IF(blockopen /= 8) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          CASE ("<"//stat_block//">") ! start statistics
             blocklines(9,1) = ii
             IF(blockopen /= 0) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 9
             openname = "<"//stat_block//">"
          CASE ("</"//stat_block//">") ! end statistics
             blocklines(9,2) = ii
             IF(blockopen /= 9) CALL abort("found "//readtag(1:twidths(1))//" while in "//openname)
             blockopen = 0
          END SELECT
       END IF       
    END DO
    
    ! Read the blocks in right order
    ! (if they have been defined)
    
    ! main
    blnr = 1
    IF(blocklines(blnr,1) /= 0)THEN ! the block exists
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN ! the block size is positive
          ! parse the data in the block
          CALL parse_main(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               control)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          CALL abort("The block "//main_block//" is empty")
       ELSE         
          CALL abort("The block "//main_block//" is unclosed")
       END IF
    ELSE
       CALL abort("The block "//main_block//" is missing")
    END IF
    
    ! make sure everyone has the same random seed
#ifdef MPI
    CALL MPI_BCAST(control%rnd_seed,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
#endif
    CALL genrand_init(control%rnd_seed)

    ! mb-model
    blnr = 2
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_mb(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               params)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          CALL abort("The block <"//mb_block//"> is empty")
       ELSE         
          CALL abort("The block <"//mb_block//"> is unclosed")
       END IF
    ELSE
       CALL abort("The block <"//mb_block//"> is missing")
    END IF

    ! cell
    blnr = 3
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_cell(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               cell,multi,pbc,btype,bval)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          ! this is not fatal, so only warn
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//cell_block//"> is empty"
          cell = cell_def
          pbc = pbc_def
          btype = free_bound_index
          multi = 1
       ELSE         
          CALL abort("The block <"//cell_block//"> is unclosed")
       END IF
    ELSE ! block missing completely
       cell = cell_def
       pbc = pbc_def
       btype = free_bound_index
       multi = 1
    END IF
 
    ! set scalability to false for free bounds and true for others
    DO ii = 1, 3 ! loop x, y, z
       IF(btype(ii) == free_bound_index)THEN
          control%scalable(ii) = .false.
       ELSE
          control%scalable(ii) = .true.
       END IF
    END DO

    ! elements
    blnr = 4
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_elements(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               params,n_elems)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//ele_block//"> is empty, only MB molecules allowed"
          ! defs
          NULLIFY(params%m_atoms)
          NULLIFY(params%atomic_labels)
          NULLIFY(params%pot_types)
          NULLIFY(params%pot_cut)
          NULLIFY(params%n_pots)
          NULLIFY(params%pot_params)
          n_elems = 0
       ELSE         
          CALL abort("The block <"//ele_block//"> is unclosed")
       END IF
    ELSE ! block is missing
       ! defs
       NULLIFY(params%m_atoms)
       NULLIFY(params%atomic_labels)
       NULLIFY(params%pot_cut)
       NULLIFY(params%pot_types)
       NULLIFY(params%n_pots)
       NULLIFY(params%pot_params)
       n_elems = 0
    END IF

    ! positions
    blnr = 5
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_particles(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               n_elems,mbs,ats,params,found_types)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          CALL abort("The block <"//pos_block//"> is empty")
       ELSE         
          CALL abort("The block <"//pos_block//"> is unclosed")
       END IF
    ELSE
       CALL abort("The block <"//pos_block//"> is missing")
    END IF

    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    ! check for periodicity
    DO ii = 1, n_mbs ! loop over mbs
       DO jj = 1, 3 ! loop over x, y, z
          IF(pbc(jj))THEN 
             ! move the particle by multiples of the cell size to bring it inside the cell, if it is not
             CALL move_in_cell(mbs(ii)%pos(jj),cell(jj))
          END IF
       END DO
       ! record the initial position
       mbs(ii)%inipos = mbs(ii)%pos
    END DO
    DO ii = 1, n_ats ! loop over atoms
       DO jj = 1, 3 ! loop over x, y, z
          IF(pbc(jj))THEN
             ! move the particle by multiples of the cell size to bring it inside the cell, if it is not
             CALL move_in_cell(ats(ii)%pos(jj),cell(jj))
          END IF
       END DO
       ! record the initial position
       ats(ii)%inipos = ats(ii)%pos
    END DO

    ! velocities
    blnr = 6
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_velocities(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               mbs,ats,random_mb_vel,random_at_vel)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//vel_block//"> is empty"
          ! random velocities everywhere
          ALLOCATE(random_mb_vel(SIZE(mbs(:)),2))
          ALLOCATE(random_at_vel(SIZE(ats(:))))
          random_mb_vel = .true.
          random_at_vel = .true.
       ELSE         
          CALL abort("The block <"//pos_block//"> is unclosed")
       END IF
    ELSE
       ! random velocities everywhere
       ALLOCATE(random_mb_vel(SIZE(mbs(:)),2))
       ALLOCATE(random_at_vel(SIZE(ats(:))))
       random_mb_vel = .true.
       random_at_vel = .true.
    END IF

    ! constraints
    blnr = 7
    is_constr = .false.
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_constraints(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               mbs,ats)
          is_constr = .true.
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          ! empty constraints block
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//constr_block//"> is empty"
          DO ii = 1, n_mbs
             mbs(ii)%constrained = no_constr_index
             mbs(ii)%well = 0.d0
          END DO
          DO ii = 1, n_ats
             ats(ii)%constrained = no_constr_index
             ats(ii)%well = 0.d0
          END DO
       ELSE         
          CALL abort("The block <"//constr_block//"> is unclosed")
       END IF
    ELSE
       ! no constraints block
       DO ii = 1, n_mbs
          mbs(ii)%constrained = no_constr_index
          mbs(ii)%well = 0.d0
       END DO
       DO ii = 1, n_ats
          ats(ii)%constrained = no_constr_index
          ats(ii)%well = 0.d0
       END DO
    END IF

    ! Check for velocities of constrained particles
    DO ii = 1, n_mbs ! loop over mbs
       DO jj = 1, 3 ! loop over x, y, z
          IF(mbs(ii)%constrained(jj) == frozen_pos_index)THEN ! frozen position
             mbs(ii)%vel(jj) = 0.d0
             !random_mb_vel(ii,1) = .false. ! no random velocity
          ELSE IF(mbs(ii)%constrained(jj) == all_frozen_index)THEN ! all frozen
             mbs(ii)%vel(jj) = 0.d0 
             mbs(ii)%angvel(jj) = 0.d0
             !random_mb_vel(ii,:) = .false. ! no random velocity or angular velocity
          END IF
       END DO
    END DO
    DO ii = 1, n_ats ! loop over atoms
       DO jj = 1, 3 ! loop over x, y, z
          IF(ats(ii)%constrained(jj) == frozen_pos_index)THEN ! frozen position
             ats(ii)%vel(jj) = 0.d0
             !random_at_vel(ii) = .false.
          ELSE IF(ats(ii)%constrained(jj) == all_frozen_index)THEN ! all frozen
             ats(ii)%vel(jj) = 0.d0 
             !random_at_vel(ii) = .false.
          END IF
       END DO
    END DO

    ! potentials
    blnr = 8
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_potentials(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               n_elems,params,found_types)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//pot_block//"> is empty"
          NULLIFY(params%pot_types)
          NULLIFY(params%pot_cut)
          NULLIFY(params%n_pots)
          NULLIFY(params%pot_params)
          IF(n_ats > 0)THEN
             CALL abort("There are atomic particles with no interactions")
          END IF
          IF(n_elems > 0)THEN
             ALLOCATE(params%pot_cut(-1:n_elems,-1:n_elems))
             ALLOCATE(params%pot_types(1,-1:n_elems,-1:n_elems))
             ALLOCATE(params%n_pots(-1:n_elems,-1:n_elems))
             ALLOCATE(params%pot_params(1,1,-1:n_elems,-1:n_elems))
             params%pot_cut = 0.d0
             params%n_pots = 0
             params%pot_types = 0
             params%pot_params = 0.d0
          END IF
       ELSE         
          CALL abort("The block <"//pot_block//"> is unclosed")
       END IF
    ELSE
       ! no potentials
       NULLIFY(params%pot_cut)
       NULLIFY(params%pot_types)
       NULLIFY(params%n_pots)
       NULLIFY(params%pot_params)
       IF(n_ats > 0)THEN
          CALL abort("There are atomic particles with no interactions")
       END IF
       IF(n_elems > 0)THEN
          ALLOCATE(params%pot_types(1,-1:n_elems,-1:n_elems))
          ALLOCATE(params%pot_cut(-1:n_elems,-1:n_elems))
          ALLOCATE(params%n_pots(-1:n_elems,-1:n_elems))
          ALLOCATE(params%pot_params(1,1,-1:n_elems,-1:n_elems))
          params%n_pots = 0
          params%pot_types = 0
          params%pot_cut = 0.d0
          params%pot_params = 0.d0
       END IF
    END IF

    ! statistics
    blnr = 9
    IF(blocklines(blnr,1) /= 0)THEN
       IF(blocklines(blnr,2) > blocklines(blnr,1)+1)THEN          
          CALL parse_statistics(data(:,blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               lwidths(blocklines(blnr,1)+1:blocklines(blnr,2)-1),&
               control,mbs,ats)
       ELSE IF(blocklines(blnr,2) == blocklines(blnr,1)+1)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "The block <"//stat_block//"> is empty"
          ! no stats
          NULLIFY(control%stats)
          control%n_statfiles = 0
          ALLOCATE(control%n_stats(1))
          ALLOCATE(control%stat_start(1))
          ALLOCATE(control%stat_interval(1))
          control%n_stats = 0
          control%stat_start = control%time_max*2.d0
          control%stat_interval = control%time_max*2.d0
          control%adf_writer = rdf_writer_def
          control%rdf_writer = rdf_writer_def
          control%bond_writer = rdf_writer_def
       ELSE         
          CALL abort("The block <"//ele_block//"> is unclosed")
       END IF
    ELSE ! block is missing
       ! no stats
       NULLIFY(control%stats)
       control%n_statfiles = 0
       ALLOCATE(control%n_stats(1))
       ALLOCATE(control%stat_start(1))
       ALLOCATE(control%stat_interval(1))
       control%n_stats = 0
       control%stat_start = control%time_max*2.d0
       control%stat_interval = control%time_max*2.d0
       control%adf_writer = rdf_writer_def
       control%rdf_writer = rdf_writer_def
       control%bond_writer = rdf_writer_def
    END IF
    IF(n_mbs == 0)THEN
       control%bond_writer = noxyz_index
    END IF

    ! if there are unused elements, remove them from lists
    CALL drop_unused(n_elems,found_types,params)

    ! multiply cell if requested
    ! this will expand the cell vectors and copy the particles in the cell
    CALL copy_cell(mbs,ats,cell,multi,random_mb_vel,random_at_vel,n_mbs,n_ats,found_types)

    ! get the maximum interaction ranges
    params%max_mb_cut = MAX(params%cut_hb,params%cut_lj)+params%cut_ver
    params%max_at_cut = 0.d0
    params%max_ma_cut = 0.d0
    IF(n_elems > 0)THEN ! there are atoms
       DO ii = -1, n_elems ! loop over elements (including mb and mb arms)
          DO jj = -1, n_elems ! -"-
             ! check for a maximum
             IF(MIN(ii,jj) <= 0)THEN ! mb-atom
                IF(params%pot_cut(jj,ii)+params%cut_ver+params%l_oh > params%max_ma_cut)THEN ! max cutoff
                   params%max_ma_cut = params%pot_cut(jj,ii)+params%cut_ver+params%l_oh
                END IF
             ELSE ! atom-atom
                IF(params%pot_cut(jj,ii)+params%cut_ver > params%max_at_cut)THEN
                   params%max_at_cut = params%pot_cut(jj,ii)+params%cut_ver
                END IF
             END IF
             ! check for non-positive cutoffs
             IF(params%n_pots(jj,ii) > 0 .AND. params%pot_cut(jj,ii) < norm_tolerance)THEN
                IF(ii == -1)THEN
                   label1 = armname
                ELSE IF(ii == 0)THEN
                   label1 = mbname
                ELSE
                   label1 = params%atomic_labels(ii)
                END IF
                IF(jj == -1)THEN
                   label2 = armname
                ELSE IF(jj == 0)THEN
                   label2 = mbname
                ELSE
                   label2 = params%atomic_labels(jj)
                END IF
                CALL abort("Potential "//label1//" - "&
                     //label2//" has zero range.")
             END IF
          END DO
       END DO
    END IF
    ! pick the maximum
    params%max_cut = MAX(MAX(params%max_mb_cut,params%max_at_cut),params%max_ma_cut)
    ! check for too small cell:
    ! if the interaction range is more than half the cell length in a periodic direction
    ! it is possible that a particle sees another from two directions.
    ! the code cannot handle this and it must be prevented.
    DO ii = 1, 3 ! loop over x, y, z
       IF(pbc(ii) .AND. cell(ii) < 2.d0*(params%max_cut-params%cut_ver))THEN
          CALL abort("interaction range is longer than size of simulation cell, increase system size")
       END IF
    END DO
    ! get default rdf range if needed
    IF(control%rdf_writer /= noxyz_index)THEN
       IF(control%rdf_range <= 0.d0)THEN
          control%rdf_range = MIN(params%max_cut-params%cut_ver,MINVAL(cell(:))/2.1d0)
       END IF
       IF(control%rdf_range > MINVAL(cell(:))/2.d0)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "rdf range is longer than size of simulation cell, correcting"
          control%rdf_range = MINVAL(cell(:))/2.1d0
       END IF
       
       ! correct rdf type
       IF(n_ats == 0)THEN
          IF(control%rdf_particles == all_rdf_index)THEN
             control%rdf_particles = mb_rdf_index
          ELSE IF(control%rdf_particles == atom_rdf_index)THEN
             IF(CPU_ID == MASTER_CPU) WRITE(*,*) "rdf requested for atoms, but no atoms exist - changing to MBs"
             control%rdf_particles = mb_rdf_index          
          END IF
       END IF
       IF(n_mbs == 0)THEN
          IF(control%rdf_particles == all_rdf_index)THEN
             control%rdf_particles = atom_rdf_index
          ELSE IF(control%rdf_particles == mb_rdf_index)THEN
             IF(CPU_ID == MASTER_CPU) WRITE(*,*) "rdf requested for MBs, but no MBs exist - changing to atoms"
             control%rdf_particles = atom_rdf_index          
          END IF
       END IF
    END IF
    ! get default adf range if needed
    IF(control%adf_writer /= noxyz_index)THEN
       IF(control%adf_range <= 0.d0)THEN
          control%adf_range = MIN(params%max_cut-params%cut_ver,MINVAL(cell(:))/2.1d0)
       END IF
       SELECT CASE(control%adf_axis)
       CASE(x_index)
          adfmax = MIN(cell(2),cell(3))*0.5d0
       CASE(y_index)
          adfmax = MIN(cell(1),cell(3))*0.5d0
       CASE(z_index)
          adfmax = MIN(cell(2),cell(1))*0.5d0
       END SELECT
       IF(control%adf_range > adfmax)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "adf range is longer than size of simulation cell, correcting"
          control%adf_range = adfmax*0.48d0
       END IF
       
       ! correct adf type
       IF(n_ats == 0)THEN
          IF(control%adf_particles == all_rdf_index)THEN
             control%adf_particles = mb_rdf_index
          ELSE IF(control%adf_particles == atom_rdf_index)THEN
             IF(CPU_ID == MASTER_CPU) WRITE(*,*) "adf requested for atoms, but no atoms exist - changing to MBs"
             control%adf_particles = mb_rdf_index          
          END IF
       END IF
       IF(n_mbs == 0)THEN
          IF(control%adf_particles == all_rdf_index)THEN
             control%adf_particles = atom_rdf_index
          ELSE IF(control%adf_particles == mb_rdf_index)THEN
             IF(CPU_ID == MASTER_CPU) WRITE(*,*) "adf requested for MBs, but no MBs exist - changing to atoms"
             control%adf_particles = atom_rdf_index          
          END IF
       END IF
    END IF

    ! assign random velocities
    CALL assign_random_velocities(mbs,ats,random_mb_vel,random_at_vel,control%ini_temp)

    ! Check for velocities of constrained particles
    DO ii = 1, n_mbs ! loop over mbs
       DO jj = 1, 3 ! loop over x, y, z
          IF(mbs(ii)%constrained(jj) == frozen_pos_index)THEN ! frozen position
             mbs(ii)%vel(jj) = 0.d0
             !random_mb_vel(ii,1) = .false. ! no random velocity
          ELSE IF(mbs(ii)%constrained(jj) == all_frozen_index)THEN ! all frozen
             mbs(ii)%vel(jj) = 0.d0 
             mbs(ii)%angvel(jj) = 0.d0
             !random_mb_vel(ii,:) = .false. ! no random velocity or angular velocity
          END IF
       END DO
    END DO
    DO ii = 1, n_ats ! loop over atoms
       DO jj = 1, 3 ! loop over x, y, z
          IF(ats(ii)%constrained(jj) == frozen_pos_index)THEN ! frozen position
             ats(ii)%vel(jj) = 0.d0
             !random_at_vel(ii) = .false.
          ELSE IF(ats(ii)%constrained(jj) == all_frozen_index)THEN ! all frozen
             ats(ii)%vel(jj) = 0.d0 
             !random_at_vel(ii) = .false.
          END IF
       END DO
    END DO

    ! count the degrees of freedom
    IF(is_constr)THEN
       n_freedom = 0
       DO ii = 1, n_mbs ! loop over mbs
          DO jj = 1, 3 ! loop over x, y, z
             SELECT CASE(mbs(ii)%constrained(jj))
             CASE(all_frozen_index) ! no freedom
             CASE(frozen_pos_index) ! rotational degree of freedom
                n_freedom = n_freedom+1                
             CASE DEFAULT ! rotational and linear movement degrees of freedom
                n_freedom = n_freedom+2
             END SELECT
          END DO
       END DO
       DO ii = 1, n_ats ! loop over atoms
          DO jj = 1, 3 ! loop over x, y, z
             SELECT CASE(ats(ii)%constrained(jj))
             CASE(all_frozen_index) ! no freedom
             CASE(frozen_pos_index)
             CASE DEFAULT
                n_freedom = n_freedom+1
             END SELECT
          END DO
       END DO
    ELSE ! if no constraints are present, the mbs have 6 and atoms 3 degrees of freedom
       n_freedom = 6*n_mbs + 3*n_ats
    END IF
    
    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)
    IF(associated(data))THEN
       DEALLOCATE(data)
    END IF
    NULLIFY(data)

    RETURN
  END SUBROUTINE parse_input


  ! copies the system in x, y and z directions the number of times
  ! specified in the argument "multi". This means expanding the supercell
  ! and also copying the molecules and atoms. The copied "image" particles
  ! have the same properties as the ones in the original cell, except for the
  ! positions, of course.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true for periodic boundaries
  ! *multi cell multipliers for x, y and z directions
  ! *gen_mb_vel logic list, true for mb molecules that should have their velocity assigned randomly
  ! *gen_at_vel logic list, true for atoms that should have their velocity assigned randomly
  ! *n_mbs number of molecules
  ! *n_ats number of atoms
  ! *found_types numbers of atoms of each type
  SUBROUTINE copy_cell(mbs,ats,cell,multi,gen_mb_vel,gen_at_vel,n_mbs,n_ats,found_types)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: cell(3)
    TYPE(mb), POINTER :: mbs(:)
    TYPE(mb), ALLOCATABLE :: temp_mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(atom), ALLOCATABLE :: temp_ats(:)
    INTEGER, INTENT(IN) :: multi(3)
    INTEGER, INTENT(INOUT) :: n_mbs, n_ats
    INTEGER, POINTER :: found_types(:)
    INTEGER :: ii, ix, iy, iz, plier, ind
    LOGICAL, POINTER :: gen_mb_vel(:,:), gen_at_vel(:)
    LOGICAL, ALLOCATABLE :: temp_gen_mb(:,:), temp_gen_at(:)

    ! volume multiplyer
    plier = multi(1)*multi(2)*multi(3)

    ALLOCATE(temp_mbs(n_mbs*plier))
    ALLOCATE(temp_ats(n_ats*plier))
    ALLOCATE(temp_gen_mb(n_mbs*plier,2))
    ALLOCATE(temp_gen_at(n_ats*plier))
    
    IF(n_ats > 0)THEN
       DO ii = 1, SIZE(found_types(:))
          found_types(ii) = found_types(ii)*plier
       END DO
    END IF

    ! Multiply the cell in x-, y- and z- directions, including the particles
    DO ix = 1, multi(1) ! loop over duplicates in x direction
       DO iy = 1, multi(2) ! loop over duplicates in y direction
          DO iz = 1, multi(3) ! loop over duplicates in z direction
             DO ii = 1, n_mbs ! loop over mbs
                ind = ii + ( iz-1 + (iy-1)*multi(3) + (ix-1)*multi(3)*multi(2) )*n_mbs ! index of the new particle (in mbs list)
                temp_mbs(ind)%pos(1) = mbs(ii)%pos(1) + (ix-1)*cell(1) ! x coordinate of copied particle
                temp_mbs(ind)%pos(2) = mbs(ii)%pos(2) + (iy-1)*cell(2) ! y -"-
                temp_mbs(ind)%pos(3) = mbs(ii)%pos(3) + (iz-1)*cell(3) ! z -"-
                temp_mbs(ind)%orientation = mbs(ii)%orientation ! orientation remains
                temp_mbs(ind)%vel = mbs(ii)%vel ! velocity remains
                temp_mbs(ind)%angvel = mbs(ii)%angvel ! ang. velocity remains
                temp_mbs(ind)%well = mbs(ii)%well ! constraints remain
                temp_mbs(ind)%constrained = mbs(ii)%constrained
                temp_mbs(ind)%inipos = temp_mbs(ind)%pos ! initial position remains
                temp_mbs(ind)%m_tot = mbs(ii)%m_tot ! mass remains
                temp_mbs(ind)%m_inert = mbs(ii)%m_inert ! mom. of inertia remains
                temp_gen_mb(ind,:) = gen_mb_vel(ii,:) ! random velocity tag remains
                temp_mbs(ind)%index = ind ! index is new
             END DO
             DO ii = 1, n_ats ! loop over atoms
                ind = ii + ( iz-1 + (iy-1)*multi(3) + (ix-1)*multi(3)*multi(2) )*n_ats ! index (in atoms list)
                temp_ats(ind)%pos(1) = ats(ii)%pos(1) + (ix-1)*cell(1)
                temp_ats(ind)%pos(2) = ats(ii)%pos(2) + (iy-1)*cell(2)
                temp_ats(ind)%pos(3) = ats(ii)%pos(3) + (iz-1)*cell(3)
                temp_ats(ind)%vel = ats(ii)%vel
                temp_ats(ind)%well = ats(ii)%well
                temp_ats(ind)%constrained = ats(ii)%constrained
                temp_ats(ind)%inipos = temp_ats(ind)%pos
                temp_ats(ind)%type = ats(ii)%type
                temp_ats(ind)%element = ats(ii)%element
                temp_ats(ind)%mass = ats(ii)%mass
                temp_gen_at(ind) = gen_at_vel(ii)
                temp_ats(ind)%index = ind + n_mbs*plier ! index is the atomic index plus the index of the last mb
             END DO
          END DO
       END DO
    END DO

    ! de- and reallocate arrays
    IF(associated(mbs))THEN
       DEALLOCATE(mbs)
    END IF
    NULLIFY(mbs)
    IF(associated(ats))THEN
       DEALLOCATE(ats)
    END IF
    NULLIFY(ats)
    IF(associated(gen_mb_vel))THEN
       DEALLOCATE(gen_mb_vel)
    END IF
    NULLIFY(gen_mb_vel)
    IF(associated(gen_at_vel))THEN
       DEALLOCATE(gen_at_vel)
    END IF
    NULLIFY(gen_at_vel)
    ALLOCATE(mbs(n_mbs*plier))
    ALLOCATE(ats(n_ats*plier))
    ALLOCATE(gen_mb_vel(n_mbs*plier,2))
    ALLOCATE(gen_at_vel(n_ats*plier))

    mbs = temp_mbs
    ats = temp_ats
    gen_mb_vel = temp_gen_mb
    gen_at_vel = temp_gen_at

    DEALLOCATE(temp_mbs)
    DEALLOCATE(temp_ats)
    DEALLOCATE(temp_gen_mb)
    DEALLOCATE(temp_gen_at)
    
    ! expand the cell
    DO ii = 1,3
       IF(multi(ii) > 1)THEN
          cell(ii) = REAL(multi(ii),KIND=dp)*cell(ii)
       END IF
    END DO

    ! expand the counters for numbers of particles
    n_mbs = n_mbs*plier
    n_ats = n_ats*plier

    RETURN
  END SUBROUTINE copy_cell


  ! parses information in the "control" block
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *control control parameters
  ! *params physical parameters
  SUBROUTINE parse_main(data,lwidths,control)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(cps), INTENT(INOUT) :: control
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, iostat, n_lines, n_tokens, t_ind, t_ind2
    CHARACTER(LEN=100) :: readtag

    n_lines = SIZE(data(1,:))

    ! Set all control params to null or default
    control%ini_temp = ini_temp_def
    control%ini_press = ini_press_def
    control%time_max = missing_r
    control%time_step = missing_r
    control%md_algo = md_algo_def
    control%md_thermo = md_thermo_def
    control%md_baro = md_baro_def
    control%thermo_value = 0.d0
    control%baro_value = 0.d0
    control%baro_interval = baro_interval_def
    control%xyz_writer = xyz_writer_def
    control%xyz_interval = missing_r
    control%inp_writer = inp_writer_def
    control%inp_interval = missing_r
    control%runtype = run_def
    control%verbose = verbose_def
    control%verblevel = verblevel_def
    control%unitscale = units_def

    ! rng seed was already chosen according to time
    ! one can overwrite it, but no def value shouldn't
    ! be set

    DO ii = 1, n_lines ! loop over the data lines
       ! split data to words
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:)) ! number of words

       IF(n_tokens > 1)THEN
          t_ind = 1
          readtag = tokens(t_ind)
          t_ind2 = 2

          ! compare the read tag to known keywords
          SELECT CASE(readtag(1:twidths(t_ind)))
             
          CASE(Temptag) ! temperature
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%ini_temp
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             
          CASE(Pressuretag)
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%ini_press
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))

          CASE(Tsimutag) ! simulation time
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%time_max
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind))) 
             
          CASE(Tsteptag) ! time step
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%time_step
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind))) 
             
          CASE(CGstepstag) ! max cg steps
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%time_max
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind))) 
             
          CASE(CGtoltag) ! cg tolerance
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%time_step
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind))) 

          CASE(algotag) ! algorithm
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
             CASE(verlettag)
                control%md_algo = velocity_verlet_index
             CASE(leapfrogtag)
                control%md_algo = leapfrog_index
             CASE(cgtag)
                control%md_algo = cg_index
             CASE DEFAULT
                CALL abort_value(algotag,main_block)
             END SELECT

          CASE(unitstag) ! unit scale
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
             CASE(truetag)
                control%unitscale = true_units
                kb = true_kb
                m_scale = true_m_scale
             CASE(reducedtag)
                control%unitscale = reduced_units
                kb = 1.d0
                m_scale = 1.d0
             CASE DEFAULT
                CALL abort_value(unitstag,main_block)
             END SELECT

          CASE(thermotag) ! thrmostat
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
             CASE(nonetag)
                control%md_thermo = microcanonical_index
                control%thermo_value = 0.d0
             CASE(langevintag)
                control%md_thermo = langevin_index
                IF(n_tokens > t_ind2)THEN
                   READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%thermo_value
                   IF(iostat /= 0) CALL abort_value(thermotag,main_block)
                ELSE
                   CALL abort_value(thermotag,main_block)
                END IF
             CASE(coolertag)
                control%md_thermo = cooler_index
                IF(n_tokens > t_ind2)THEN
                   READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%thermo_value
                   IF(iostat /= 0) CALL abort_value(thermotag,main_block)
                ELSE
                   CALL abort_value(thermotag,main_block)
                END IF
             CASE DEFAULT
                CALL abort_value(thermotag,main_block)
             END SELECT

          CASE(barotag) ! barostat
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
             CASE(nonetag)
                control%md_baro = microcanonical_index
                control%baro_value = 0.d0
             CASE(berendsentag)
                control%md_baro = berendsen_index
                IF(n_tokens > t_ind2)THEN
                   READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%baro_value
                   IF(iostat /= 0) CALL abort_value(barotag,main_block)
                   IF(n_tokens > t_ind2+1)THEN
                      READ(tokens(t_ind2+2),*,IOSTAT=iostat) control%baro_interval
                      IF(iostat /= 0) CALL abort_value(barotag,main_block)
                   END IF
                ELSE
                   CALL abort_value(barotag,main_block)
                END IF
             CASE DEFAULT
                CALL abort_value(barotag,main_block)
             END SELECT

          CASE(seedtag) ! random number seed
             READ(tokens(t_ind2),*,IOSTAT=iostat) control%rnd_seed
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))

          CASE(verbosetag) ! verbosity
             IF(tokens(t_ind2)(1:twidths(t_ind2)) == nonetag)THEN
                control%verbose = missing_r
                control%verblevel = 0
             ELSE
                READ(tokens(t_ind2),*,IOSTAT=iostat) control%verbose
                IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             END IF
             IF(n_tokens > t_ind2)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%verblevel
                IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             END IF

          CASE(xyztag) ! xyz file writer
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
                CASE(nonetag)
                   control%xyz_writer = noxyz_index
                   control%xyz_interval = 0.d0
                CASE(endtag)
                   control%xyz_writer = exyz_index
                   control%xyz_interval = 0.d0
                CASE(starttag)
                   control%xyz_writer = sxyz_index
                   control%xyz_interval = 0.d0
                CASE(sendtag)
                   control%xyz_writer = sexyz_index
                   control%xyz_interval = 0.d0
                CASE(intervaltag)
                   control%xyz_writer = ixyz_index
                   IF(n_tokens > t_ind2)THEN
                      READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%xyz_interval
                      IF(iostat /= 0) CALL abort("reading xyz interval")
                   ELSE
                      CALL abort_value(xyztag,main_block)
                   END IF
                CASE DEFAULT
                   CALL abort_value(xyztag,main_block)
             END SELECT

          CASE(conttag) ! continuation file writer
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
                CASE(nonetag)
                   control%inp_writer = noxyz_index
                   control%inp_interval = 0.d0
                CASE(endtag)
                   control%inp_writer = exyz_index
                   control%inp_interval = 0.d0
                CASE(starttag)
                   control%inp_writer = sxyz_index
                   control%inp_interval = 0.d0
                CASE(sendtag)
                   control%inp_writer = sexyz_index
                   control%inp_interval = 0.d0
                CASE(intervaltag)
                   control%inp_writer = ixyz_index
                   IF(n_tokens > t_ind2)THEN
                      READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%inp_interval
                      IF(iostat /= 0) CALL abort("reading xyz interval")
                   ELSE
                      CALL abort_value(conttag,main_block)
                   END IF
                CASE DEFAULT
                   CALL abort_value(conttag,main_block)
             END SELECT

          CASE(runtag) ! run type tag
             SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
             CASE(testtag)
                control%runtype = test_index
             CASE(fulltag)
                control%runtype = full_index
             CASE DEFAULT
                CALL abort_value(runtag,main_block)
             END SELECT
             
          CASE DEFAULT ! unrecognized tag
             CALL write_typo(readtag(1:twidths(t_ind)),main_block)
             
          END SELECT
       ELSE ! only one word on the line - the tag is missing its value
          CALL write_miss(tokens(1)(1:twidths(1)),main_block)
       END IF

    END DO

    ! check validity
    IF(control%ini_temp < 0.d0) CALL abort_value(Temptag,main_block)
    IF(control%ini_press < 0.d0) CALL abort_value(Pressuretag,main_block)
    IF(control%time_max <= 0.d0 .AND. control%md_algo == cg_index) control%time_max = cg_max_def
    IF(control%time_max <= 0.d0) CALL abort_value(Tsimutag,main_block)
    IF(control%time_step <= 0.d0 .AND. control%md_algo /= cg_index) CALL abort_value(Tsteptag,main_block)
    IF(control%md_algo == cg_index .AND. control%time_step == missing_r)THEN
       control%time_step = cg_tol_def
    END IF
    IF(control%verbose == missing_r)THEN
        control%verbose = 2.d0*control%time_max
    ELSE IF(control%verbose == verbose_def)THEN
       control%verbose = control%time_max/20.d0
       ! negative verbose -> only print to file
    END IF    
    IF(control%verblevel < -2)THEN
       control%verblevel = -1
    ELSE IF(control%verblevel > 10)THEN
       control%verblevel = 10
    END IF
    IF(control%xyz_writer == ixyz_index .AND. control%xyz_interval <= 0.d0) &
         CALL abort_value(xyztag,main_block)
    IF(control%inp_writer == ixyz_index .AND. control%inp_interval <= 0.d0) &
         CALL abort_value(conttag,main_block)
    IF(control%rnd_seed < 1 .AND. control%rnd_seed /= rnd_seed_def) CALL abort_value(seedtag,main_block)
    IF(control%md_thermo < 1 .OR. control%md_thermo > n_thermos) CALL abort_value(thermotag,main_block)
    IF(control%md_algo < 1 .OR. control%md_algo > n_algos) CALL abort_value(algotag,main_block)
    IF(control%thermo_value < 0.d0) CALL abort_value(thermotag,main_block)
    IF(control%baro_value < 0.d0) CALL abort_value(barotag,main_block)

    IF(control%md_thermo == langevin_index .AND. control%md_algo == leapfrog_index)THEN
       IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Langevin is only implemented for Velocity Verlet - will change the algorithm"
       control%md_algo = velocity_verlet_index
    END IF
    IF(control%md_baro == berendsen_index .AND. control%md_algo == leapfrog_index)THEN
       IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Berendsen is only implemented for Velocity Verlet - will change the algorithm"
       control%md_algo = velocity_verlet_index
    END IF

    IF(control%md_algo == cg_index)THEN
       control%md_thermo = microcanonical_index       
    END IF

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_main

  ! parses information in the "mb-model" block
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *control control parameters
  ! *params physical parameters
  SUBROUTINE parse_mb(data,lwidths,params)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mbps), INTENT(INOUT) :: params
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, iostat, n_lines, n_tokens, t_ind, t_ind2
    CHARACTER(LEN=100) :: readtag

    n_lines = SIZE(data(1,:))

    ! Set all physical params to null of default
    params%e_lj = missing_r
    params%e_hb = missing_r
    params%e_phi = missing_r
    params%r_hb = missing_r
    params%s_lj = missing_r
    params%s_r = missing_r
    params%s_th = missing_r
    params%v_b = missing_r
    params%r_b = missing_r
    params%d_b = missing_r
    params%cut_hb = missing_r
    params%cut_lj = missing_r
    params%cut_ver = missing_r
    params%cut_atom = 0.d0
    params%m_mol = missing_r
    params%I_mol = missing_r
    params%L_oh = L_oh_def


    DO ii = 1, n_lines ! loop over the lines of data
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       t_ind = 1
       readtag = tokens(t_ind)
       t_ind2 = 2

       IF(n_tokens > 1)THEN
          
          ! compare the tags to known keywords
          SELECT CASE(readtag(1:twidths(t_ind)))
             
          CASE(E_LJtag) ! Lennard-Jones energy
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%e_lj
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(E_HBtag) ! H-Bond energy
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%e_hb
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(E_PHtag) ! dihedral angle energy
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%e_phi
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(R_HBtag) ! H-Bond equilibrium distance
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%r_hb
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(Sigma_LJtag) ! Lennard-Jones sigma
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%s_lj
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%inv_slj = 1.d0/params%s_lj
          CASE(Sigma_HBtag) ! H-Bond distance sigma
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%s_r
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%inv_sr = 1.d0/params%s_r
          CASE(Sigma_THtag) ! H-Bond angle sigma
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%s_th
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%inv_sth = 1.d0/params%s_th
          CASE(Exp_btag) ! bond order exponent
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%v_b
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(R_btag) ! bond order inner shell distance
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%r_b
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(D_btag) ! bond order outer shell distance
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%d_b
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%inv_db = 0.5d0*pi/params%d_b
          CASE(Cut_HBtag) ! H-bond cutoff
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%cut_hb
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(Cut_LJtag) ! Lennard-Jones cutoff
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%cut_lj
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(Cut_Vtag) ! Verlet neighbor list marginal
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%cut_ver
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(Cut_Atag) ! Atomic potential cutoff
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%cut_atom
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
          CASE(M_moletag) ! MB mass
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%m_mol
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%m_mol = params%m_mol*m_scale
             params%inv_m = 1.d0 / params%m_mol
          CASE(I_moletag) ! MB moment of inertia
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%i_mol
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             params%i_mol = params%i_mol*m_scale
             params%inv_i = 1.d0 / params%i_mol
          CASE(L_ohtag) ! O-H distance
             READ(tokens(t_ind2),*,IOSTAT=iostat) params%l_oh
             IF(iostat /= 0) CALL abort("reading "//readtag(1:twidths(t_ind)))
             
          CASE DEFAULT ! unrecognized tag
             CALL write_typo(readtag(1:twidths(t_ind)),mb_block)
          END SELECT
       ELSE
          CALL write_miss(tokens(1)(1:twidths(1)),mb_block)
       END IF
    END DO

    ! The physical parameters must be defined, complain if they are not
    IF(params%e_lj == missing_r) CALL abort_value(E_LJtag,mb_block)
    IF(params%e_hb == missing_r) CALL abort_value(E_HBtag,mb_block)
    IF(params%e_phi == missing_r) CALL abort_value(E_PHtag,mb_block)
    IF(params%r_hb <= 0.d0) CALL abort_value(R_HBtag,mb_block)
    IF(params%s_lj <= 0.d0) CALL abort_value(Sigma_LJtag,mb_block)
    IF(params%s_r <= 0.d0) CALL abort_value(Sigma_HBtag,mb_block)
    IF(params%s_th <= 0.d0) CALL abort_value(Sigma_THtag,mb_block)
    IF(params%v_b <= 0.d0) CALL abort_value(Exp_btag,mb_block)
    IF(params%r_b <= 0.d0) CALL abort_value(R_btag,mb_block)
    IF(params%d_b <= 0.d0) CALL abort_value(D_btag,mb_block)
    IF(params%m_mol <= 0.d0) CALL abort_value(M_moletag,mb_block)    
    IF(params%I_mol <= 0.d0) CALL abort_value(I_moletag,mb_block)    
    IF(params%L_oh <= 0.d0) CALL abort_value(I_moletag,mb_block)
    IF(params%cut_hb <= 0.d0)THEN
       CALL abort_value(cut_hbtag,mb_block)
    END IF
    IF(params%cut_lj <= 0.d0)THEN
       CALL abort_value(cut_ljtag,mb_block)
    END IF
    IF(params%cut_ver <= 0.d0)THEN
       CALL abort_value(cut_vtag,mb_block)
    END IF
    IF(params%cut_atom < 0.d0)THEN
       CALL abort_value(cut_atag,mb_block)       
    END IF
    
    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_mb

  ! parses information in the "cell" block. 
  ! that is, reads supercell dimensions, possible multiplers
  ! and boundary types.
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *cell supercell dimensions
  ! *multi supercell multipliers
  ! *pbc true for periodic boundaries
  ! *btype boundary type index
  ! *bval boundary parameter
  SUBROUTINE parse_cell(data,lwidths,cell,multi,pbc,btype,bval)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    REAL(KIND=dp), INTENT(OUT) :: cell(3),bval(3)
    LOGICAL, INTENT(OUT) :: pbc(3)
    INTEGER, INTENT(OUT) :: multi(3),btype(3)
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, iostat, n_lines, n_tokens, t_ind
    CHARACTER(LEN=100) :: readtag

    n_lines = SIZE(data(1,:))
    ! the cell block has a definite structure with three lines
    IF(n_lines /= 3) CALL abort("Block <"//cell_block//"> does not contain 3 lines of data.")

    pbc = pbc_def
    btype = free_bound_index
    multi = 1

    DO ii = 1, 3 ! loop over x, y, z
       ! split the line into tokens
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       ! junk line - this shouldn't happen
       IF(n_tokens == 0) CALL abort("Block <"//cell_block//"> does not contain 3 lines of data.")

       t_ind = 1
       ! read cell vector length
       READ(tokens(t_ind),*,IOSTAT=iostat) cell(ii)
       IF(iostat /= 0) CALL abort("Reading cell vectors")

       ! more tokens?
       IF(n_tokens > 1)THEN

          t_ind = 2
          ! 2nd token is a multiplier?
          IF(is_number(tokens(t_ind)))THEN
             READ(tokens(t_ind),*,IOSTAT=iostat) multi(ii)
             IF(iostat /= 0) CALL abort("Reading cell multipliers")
             IF(multi(ii) < 1) CALL abort("non-positive cell multiplier")
             
             ! 3rd is the boundary condition tag?
             IF(n_tokens > 2)THEN
                t_ind = 3
                readtag = tokens(t_ind)
                ! boundary type
                SELECT CASE(readtag(1:twidths(t_ind)))
                CASE(freetag) ! free
                   pbc(ii) = .FALSE.
                   btype(ii) = free_bound_index
                   bval(ii) = 0.d0
                CASE(periodictag) ! periodic
                   pbc(ii) = .TRUE.
                   btype(ii) = periodic_bound_index
                   bval(ii) = 0.d0
                CASE(walltag) ! wall
                   pbc(ii) = .FALSE.
                   btype(ii) = wall_bound_index
                   IF(n_tokens > t_ind)THEN
                      READ(tokens(t_ind+1),*,IOSTAT=iostat) bval(ii) 
                      IF(iostat /= 0) CALL abort_value(walltag,cell_block)
                   ELSE
                      CALL abort_value(walltag,cell_block)
                   END IF
                CASE DEFAULT ! unrecognized
                   CALL write_typo(readtag(1:twidths(t_ind)),cell_block)
                END SELECT
             END IF

          ELSE ! 2nd token is the boundary condition tag?
             t_ind = 2
             readtag = tokens(t_ind)
             ! boundary type
             SELECT CASE(readtag(1:twidths(t_ind)))
             CASE(freetag) ! free
                pbc(ii) = .FALSE.
                btype(ii) = free_bound_index
                bval(ii) = 0.d0
             CASE(periodictag) ! periodic
                pbc(ii) = .TRUE.
                btype(ii) = periodic_bound_index
                bval(ii) = 0.d0
             CASE(walltag) ! wall
                pbc(ii) = .FALSE.
                btype(ii) = wall_bound_index
                IF(n_tokens > t_ind)THEN
                   READ(tokens(t_ind+1),*,IOSTAT=iostat) bval(ii) 
                   IF(iostat /= 0) CALL abort_value(walltag,cell_block)
                ELSE
                   CALL abort_value(walltag,cell_block)
                END IF
             CASE DEFAULT
                CALL write_typo(readtag(1:twidths(t_ind)),cell_block)
             END SELECT
          END IF
       ELSE ! if no boundary condition is specified, assume free bounds
          pbc(ii) = .false.
          btype(ii) = free_bound_index
          bval(ii) = 0.d0
       END IF

    END DO

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

  END SUBROUTINE parse_cell

  ! parses information in the "elements" block.
  ! that is, reads the atomic labels and masses.
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *n_elems number of different atom types
  ! *params physical parameters
  SUBROUTINE parse_elements(data,lwidths,params,n_elems)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mbps), INTENT(INOUT) :: params
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER, INTENT(OUT) :: n_elems
    INTEGER :: ii, jj, iostat, n_tokens, t_ind
    CHARACTER(LEN=100) :: readtag

    n_elems = SIZE(data(1,:))

    NULLIFY(params%m_atoms)
    NULLIFY(params%invm_atoms)
    NULLIFY(params%atomic_labels)
    ALLOCATE(params%m_atoms(n_elems))
    ALLOCATE(params%atomic_labels(n_elems))
    ALLOCATE(params%invm_atoms(n_elems))
    params%m_atoms = 1.d0
    params%invm_atoms = 1.d0
    params%atomic_labels = "  "

    DO ii = 1, n_elems ! loop over element defining lines
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       
       IF(n_tokens > 1)THEN
          
          t_ind = 1
          readtag = tokens(t_ind)
          IF(readtag(1:labelw) == mbname) CALL abort(mbname//" is a reserved name and cannot be used as label")
          IF(readtag(1:labelw) == armname) CALL abort(armname//" is a reserved name and cannot be used as label")
          
          IF(ii > 1)THEN
             ! check for duplicates
             DO jj = 1, ii-1
                IF(readtag(1:labelw) == params%atomic_labels(jj)) &
                     CALL abort("Duplicate labels found in <"//ele_block//">")             
             END DO
          END IF
          ! a valid label - store it
          params%atomic_labels(ii) = readtag(1:labelw)

          t_ind = 2
          ! read the mass
          READ(tokens(t_ind),*,IOSTAT=iostat) params%m_atoms(ii)
          IF(iostat /= 0) CALL abort("Reading atomic masses")
          params%m_atoms(ii) = params%m_atoms(ii) * m_scale
          IF(params%m_atoms(ii) > 0.d0)THEN
             params%invm_atoms(ii) = 1.d0/params%m_atoms(ii)
          ELSE
             CALL abort_value(tokens(t_ind),ele_block)
          END IF
       ELSE
          CALL abort("Mass missing for element "//tokens(1)(1:labelw))
       END IF

    END DO

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_elements

  ! parses information in the "particles" block.
  ! that is, reads the types and initial positions and orientations 
  ! of particles (and deduces the number of particles from that)
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *n_elems number of different types of atoms
  ! *params physical parameters
  ! *mbs list of mb molecules
  ! *ats list of atoms
  ! *found_types the number of atoms of each type
  SUBROUTINE parse_particles(data,lwidths,n_elems,mbs,ats,params,found_types)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(mbps), INTENT(IN) :: params
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, INTENT(INOUT) :: n_elems
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, jj, imbs, iats, iostat, &
         n_lines, n_tokens, n_mbs, n_ats, t_ind, n_found, xyznumber
    CHARACTER(LEN=100) :: readtag
    LOGICAL :: foundlabel
    INTEGER, ALLOCATABLE :: exhausted_indices(:)
    INTEGER, POINTER :: found_types(:)

    n_lines = SIZE(data(1,:))
    n_mbs = 0
    n_ats = 0
    IF(n_elems > 0)THEN
       ALLOCATE(found_types(n_elems))
       found_types = 0
    END IF

    ! check if we should read the geometry from an xyz file
    IF(n_lines == 1)THEN ! only one line
       CALL tokenize(data(1:lwidths(1),1)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       IF(n_tokens > 1)THEN
          IF(tokens(1)(1:twidths(1)) == xyzread)THEN ! command for reading an xyz geometry
             ! read the xyz geometry file name
             IF(n_tokens > 2)THEN
                ! the number of configuration to be read from a file with several configurations
                READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) xyznumber
                IF(iostat /= 0) CALL abort("reading xyz counter in <"//pos_block//">")     
                IF(xyznumber < 1)THEN
                   xyznumber = 1
                   WRITE(*,*) "non-positive xyz counter in <"//pos_block//">, set to 1"
                END IF           
             ELSE
                xyznumber = 1
             END IF
             ! read and parse the xyz file
             ! note: the current routine uses read_lines and trim which are _very_ slow for large files
             !   - the routine should be rewritten more efficiently exploiting the known format of xyz
             CALL read_xyz_particles(tokens(2)(1:twidths(2)),xyznumber,n_elems,mbs,ats,params,found_types)
             RETURN
          END IF
       END IF
    END IF

    ! read the particle data in the block and count how many particles are there
    DO ii = 1, n_lines ! loop over lines (= particles)
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       IF(n_tokens > 1)THEN
          t_ind = 2
          readtag = tokens(t_ind)
          
          foundlabel = .false.
          IF(readtag(1:labelw) == mbname)THEN ! this line defines an MB molecule
             n_mbs = n_mbs+1
             foundlabel = .true.
          ELSE ! this line defines an atomic particle
             IF(n_elems == 0) CALL abort("Only MB molecules defined, yet other particles found")
             label: DO jj = 1, SIZE(params%atomic_labels(:)) ! find the label from the list of elements
                IF(readtag(1:labelw) == params%atomic_labels(jj))THEN
                   n_ats = n_ats+1
                   found_types(jj) = found_types(jj) + 1
                   foundlabel = .true.
                   EXIT label
                END IF                
             END DO label
          END IF
          IF(.NOT.foundlabel) CALL abort("undefined label "//readtag(1:labelw)) ! unrecognized type
       END IF

    END DO

    IF(n_elems > 0)THEN ! if atomic types have been defined, check if they were really found in the list of particles
       foundlabel = .true. 
       n_found = 0
       DO ii = 1, n_elems
          IF(found_types(ii) == 0)THEN
             foundlabel = .false.
          ELSE
             n_found = n_found+1
          END IF
       END DO
       IF(.NOT.foundlabel)THEN
          IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Note: some particle types were not found in <"//pos_block//">"
       END IF       
    END IF

    ! allocate particle arrays
    NULLIFY(mbs)
    NULLIFY(ats)
    ALLOCATE(mbs(n_mbs))
    ALLOCATE(ats(n_ats))
    ALLOCATE(exhausted_indices(n_mbs+n_ats))
    exhausted_indices = 0

    imbs = 0
    iats = 0
    DO ii = 1, n_lines ! loop over the particles once more - now read and assign positions and orientations
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       IF(n_tokens > 1)THEN
          t_ind = 2
          readtag = tokens(t_ind)
          
          ! MB-molecule
          IF(readtag(1:labelw) == mbname)THEN
             imbs = imbs+1
             ! index
             READ(tokens(1),*,IOSTAT=iostat) mbs(imbs)%index
             IF(iostat /= 0) CALL abort("Reading indices in <"//pos_block//">")
             IF(mbs(imbs)%index < 1) CALL abort("negative particle index")
             ! make sure no index is repeated
             icheck: DO jj = 1, n_mbs+n_ats
                IF(mbs(imbs)%index == exhausted_indices(jj)) &
                     CALL abort("Repeated index in <"//pos_block//">")
             END DO icheck
             ! remember that this index was used
             exhausted_indices(ii) = mbs(imbs)%index
             IF(n_tokens > 4)THEN
                ! position
                READ(tokens(3),*,IOSTAT=iostat) mbs(imbs)%pos(1)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")
                READ(tokens(4),*,IOSTAT=iostat) mbs(imbs)%pos(2)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")
                READ(tokens(5),*,IOSTAT=iostat) mbs(imbs)%pos(3)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")
                
                mbs(imbs)%inipos = mbs(imbs)%pos
               
             ELSE
                CALL abort("Found particle with undefined position")
             END IF
             IF(n_tokens > 8)THEN
                ! orientation
                READ(tokens(6),*,IOSTAT=iostat) mbs(imbs)%orientation%w
                IF(iostat /= 0) CALL abort("Reading orientation in <"//pos_block//">")
                READ(tokens(7),*,IOSTAT=iostat) mbs(imbs)%orientation%x
                IF(iostat /= 0) CALL abort("Reading orientation in <"//pos_block//">")
                READ(tokens(8),*,IOSTAT=iostat) mbs(imbs)%orientation%y
                IF(iostat /= 0) CALL abort("Reading orientation in <"//pos_block//">")
                READ(tokens(9),*,IOSTAT=iostat) mbs(imbs)%orientation%z
                IF(iostat /= 0) CALL abort("Reading orientation in <"//pos_block//">")

                IF(ABS(.norm.(mbs(imbs)%orientation)-1.d0) > norm_tolerance)THEN
                   CALL norm_quaternion(mbs(imbs)%orientation)
                END IF
             ELSE
                ! random orientation
                CALL rand_unit_quaternion(mbs(imbs)%orientation)
             END IF

          ELSE ! atom
             
             ! type
             labelb: DO jj = 1, SIZE(params%atomic_labels(:))
                IF(readtag(1:labelw) == params%atomic_labels(jj))THEN
                   iats = iats+1
                   ats(iats)%type = jj
                   ats(iats)%element = params%atomic_labels(jj)
                   EXIT labelb
                END IF                
             END DO labelb

             ! index
             READ(tokens(1),*,IOSTAT=iostat) ats(iats)%index
             IF(iostat /= 0) CALL abort("Reading indices in <"//pos_block//">")
             IF(ats(iats)%index < 1) CALL abort("negative particle index")
             ! make sure no index is repeated
             icheckb: DO jj = 1, n_mbs+n_ats
                IF(ats(iats)%index == exhausted_indices(jj)) &
                     CALL abort("Repeated index in <"//pos_block//">")
             END DO icheckb
             exhausted_indices(ii) = ats(iats)%index
             IF(n_tokens > 4)THEN
                ! position
                READ(tokens(3),*,IOSTAT=iostat) ats(iats)%pos(1)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")
                READ(tokens(4),*,IOSTAT=iostat) ats(iats)%pos(2)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")
                READ(tokens(5),*,IOSTAT=iostat) ats(iats)%pos(3)
                IF(iostat /= 0) CALL abort("Reading positions in <"//pos_block//">")  

                ats(iats)%inipos = ats(iats)%pos
            
             ELSE
                CALL abort("Found particle with undefined position")
             END IF

          END IF
       END IF

    END DO

    ! assign masses - the current version has same masses for all MB's
    ! and so this is redundant data...
    DO ii = 1, n_mbs
       mbs(ii)%m_tot = params%m_mol
       mbs(ii)%m_inert = params%i_mol
    END DO
    DO ii = 1, n_ats
       ats(ii)%mass = params%m_atoms(ats(ii)%type)
    END DO

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_particles

  ! parses information in the "velocities" block.
  ! that is, reads the initial velocities and angular velocities
  ! for the given particles, identified by their indices.
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *n_elems number of different types of atoms
  ! *control control parameters
  ! *mbs list of mb molecules
  ! *ats list of atoms
  ! *rnd_mb logic list, true for mb molecules that should have their velocity assigned randomly
  ! *rnd_at logic list, true for atoms that should have their velocity assigned randomly
  SUBROUTINE parse_velocities(data,lwidths,mbs,ats,rnd_mb,rnd_at)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, jj, index, p_ind, iostat, sind, eind, &
         n_lines, n_tokens, n_mbs, n_ats, t_ind
    LOGICAL, POINTER :: rnd_mb(:,:), rnd_at(:)
    LOGICAL :: is_mb

    n_lines = SIZE(data(1,:))
    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    ! logic arrays for showing which velocities should be assigned randomly
    ALLOCATE(rnd_mb(n_mbs,2))
    rnd_mb = .true.
    ALLOCATE(rnd_at(n_ats))
    rnd_at = .true.

    DO ii = 1, n_lines ! loop over the lines of data
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       IF(n_tokens > 3)THEN ! enough data to make a velocity vector
          READ(tokens(1),*,IOSTAT=iostat) sind ! read the particle index
          IF(iostat /= 0) &
               CALL abort("invalid index "//tokens(1)(1:twidths(1))//" in <"//vel_block//">")

          ! check for a range of indices (is the second word "to"?)
          IF(tokens(2)(1:twidths(2)) == grouptag)THEN
             READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) eind
             IF(iostat /= 0 .OR. eind < sind) CALL abort("invalid index "//tokens(3)(1:twidths(3))//" in <"//vel_block//">")
             t_ind = 4
          ELSE
             eind = sind
             t_ind = 2
          END IF

          ! go over indices in range s(tart)ind(ex) ... e(nd)ind(eX)
          DO index = sind, eind

             ! search the particle with the given index
             is_mb = .false.
             p_ind = 0
             searchmb: DO jj = 1, n_mbs
                IF(mbs(jj)%index == index)THEN
                   p_ind = jj
                   is_mb = .true.
                   EXIT searchmb
                END IF
             END DO searchmb
             
             IF(.NOT.is_mb)THEN
                searchat: DO jj = 1, n_ats
                   IF(ats(jj)%index == index)THEN
                      p_ind = jj
                      EXIT searchat
                   END IF
                END DO searchat
             END IF
             
             IF(p_ind == 0) &
                  CALL abort("invalid index "//tokens(1)(1:twidths(1))//" in <"//vel_block//">")
             
             ! read the velocity
             IF(is_mb)THEN ! for an mb molecule
                rnd_mb(p_ind,1) = .false. ! no need for a random velocity
                
                IF(n_tokens > t_ind+1)THEN
                   ! velocity
                   READ(tokens(t_ind),*,IOSTAT=iostat) mbs(p_ind)%vel(1)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+1),*,IOSTAT=iostat) mbs(p_ind)%vel(2)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+2),*,IOSTAT=iostat) mbs(p_ind)%vel(3)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                ELSE
                   CALL abort("reading initial velocities")
                END IF
                ! angular velocity, if defined
                IF(n_tokens > t_ind+4)THEN
                   rnd_mb(p_ind,2) = .false.
                   READ(tokens(t_ind+3),*,IOSTAT=iostat) mbs(p_ind)%angvel(1)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+4),*,IOSTAT=iostat) mbs(p_ind)%angvel(2)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+5),*,IOSTAT=iostat) mbs(p_ind)%angvel(3)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                ELSE
                   rnd_mb(p_ind,2) = .true. ! random angular velocity needed (though, this should have been true already)
                END IF                

             ELSE ! atom
                rnd_at(p_ind) = .false. ! no need for random velocity
                
                IF(n_tokens > t_ind+1)THEN
                   ! velocity
                   READ(tokens(t_ind),*,IOSTAT=iostat) ats(p_ind)%vel(1)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+1),*,IOSTAT=iostat) ats(p_ind)%vel(2)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                   READ(tokens(t_ind+2),*,IOSTAT=iostat) ats(p_ind)%vel(3)
                   IF(iostat /= 0) CALL abort("reading initial velocities")
                ELSE
                   CALL abort("reading initial velocities")
                END IF
                IF(n_tokens > t_ind+4)THEN ! if an atomic particle has extra data, print a warning
                   IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Extra data found for an atomic particle in <"//vel_block//">."
                   IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Are you sure the data is not meant for a MB molecule?"
                END IF
             END IF
          END DO
       ELSE
          ! insufficient data on the line
          CALL write_miss(tokens(1)(1:twidths(1)),vel_block)
       END IF
    END DO

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_velocities


  ! parses information in the "constraints" block.
  ! that is, reads the types of constraints
  ! applied to the given particles, identified by their indices.
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *n_elems number of different types of atoms
  ! *params physical parameters
  ! *mbs list of mb molecules
  ! *ats list of atoms
  SUBROUTINE parse_constraints(data,lwidths,mbs,ats)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, jj, iostat, index, p_ind, sind, eind,&
         n_lines, n_tokens, n_mbs, n_ats, t_ind, t_ind2, t_ind_0
    LOGICAL :: is_mb

    n_lines = SIZE(data(1,:))
    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    ! initialize first to no constraints
    DO ii = 1, n_mbs ! loop over mbs
       DO jj = 1, 3 ! loop over x, y, z
          mbs(ii)%constrained(jj) = no_constr_index
          mbs(ii)%well(jj) = 0.d0
       END DO
    END DO
    DO ii = 1, n_ats ! loop over atoms
       DO jj = 1, 3 ! loop over x, y, z
          ats(ii)%constrained(jj) = no_constr_index
          ats(ii)%well(jj) = 0.d0
       END DO
    END DO

    ! read the constraints defined in the block
    DO ii = 1, n_lines ! loop over lines of data
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       
       IF(n_tokens > 1)THEN
          ! index of the particle affected
          READ(tokens(1),*,IOSTAT=iostat) sind
          IF(iostat /= 0) &
               CALL abort("invalid index "//tokens(1)(1:twidths(1))//" in <"//constr_block//">")

          ! check for a range of particles
          IF(tokens(2)(1:twidths(2)) == grouptag)THEN
             READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) eind
             IF(iostat /= 0 .OR. eind < sind) CALL abort("invalid index "//tokens(3)(1:twidths(3))//" in <"//vel_block//">")
             t_ind_0 = 4
          ELSE
             eind = sind
             t_ind_0 = 2
          END IF

          ! go over indices in range s(tart)ind(ex) ... e(nd)ind(eX)
          DO index = sind, eind

             ! search the particle with the given index
             is_mb = .false.
             p_ind = 0
             searchmb: DO jj = 1, n_mbs ! loop over mbs
                IF(mbs(jj)%index == index)THEN
                   p_ind = jj
                   is_mb = .true.
                   EXIT searchmb
                END IF
             END DO searchmb

             IF(.NOT.is_mb)THEN ! if no mb with the given index was found, search in atoms
                searchat: DO jj = 1, n_ats ! loop over atoms
                   IF(ats(jj)%index == index)THEN
                      p_ind = jj
                      EXIT searchat
                   END IF
                END DO searchat
             END IF
          
             IF(p_ind == 0) & ! the index was not found
                  CALL abort("invalid index "//tokens(1)(1:twidths(1))//" in <"//constr_block//">")

             ! read constraints
             IF(is_mb)THEN
                t_ind2 = t_ind_0 - 1
                DO jj = 1, 3 ! x, y, z
                   t_ind = t_ind2 + 1
                   t_ind2 = t_ind
                   
                   ! type of constraint
                   SELECT CASE(tokens(t_ind)(1:twidths(t_ind)))
                   CASE(freetag) ! no constraint
                      mbs(p_ind)%constrained(jj) = no_constr_index
                      mbs(p_ind)%well(jj) = 0.d0
                      
                   CASE(frozentag) ! frozen degrees of freedom
                      t_ind2 = t_ind + 1
                      SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
                      CASE(postag) ! position frozen
                         mbs(p_ind)%constrained(jj) = frozen_pos_index
                         mbs(p_ind)%well(jj) = 0.d0
                      CASE(alltag) ! pos + orientation frozen
                         mbs(p_ind)%constrained(jj) = all_frozen_index
                         mbs(p_ind)%well(jj) = 0.d0
                      CASE(veltag) ! velocity frozen
                         mbs(p_ind)%constrained(jj) = frozen_vel_index
                         mbs(p_ind)%well(jj) = 0.d0      
                      CASE DEFAULT ! unrecognized tag
                         CALL abort("unrecognized constraint "//tokens(t_ind2)(1:twidths(t_ind2)))
                      END SELECT
                      
                   CASE(welltag) ! harmonic well
                      t_ind2 = t_ind + 1
                      mbs(p_ind)%constrained(jj) = harmonic_well_index
                      ! read spring constant
                      READ(tokens(t_ind2)(1:twidths(t_ind2)),*,IOSTAT=iostat) mbs(p_ind)%well(jj)
                      IF(iostat /= 0) CALL abort("reading constraints")
                      IF(mbs(p_ind)%well(jj) < 0.d0) &
                           CALL abort("negative spring constant in <"//constr_block//">")
                      
                   CASE(forcetag) ! external force
                      t_ind2 = t_ind + 1
                      mbs(p_ind)%constrained(jj) = ext_force_index
                      ! read force strength
                      READ(tokens(t_ind2)(1:twidths(t_ind2)),*,IOSTAT=iostat) mbs(p_ind)%well(jj)
                      IF(iostat /= 0) CALL abort("reading constraints")
                      
                   CASE DEFAULT ! unrecognized tag
                      CALL abort("unrecognized constraint "//tokens(t_ind)(1:twidths(t_ind)))

                   END SELECT

                END DO

             ELSE ! atom

                t_ind2 = t_ind_0 - 1
                DO jj = 1, 3 ! x, y, z
                   t_ind = t_ind2 + 1
                   t_ind2 = t_ind
                   
                   ! type of constraint
                   SELECT CASE(tokens(t_ind)(1:twidths(t_ind)))
                   CASE(freetag) ! no constraint
                      ats(p_ind)%constrained(jj) = no_constr_index
                      ats(p_ind)%well(jj) = 0.d0
                   
                   CASE(frozentag) ! frozen degrees of freedom
                      t_ind2 = t_ind + 1
                      SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
                      CASE(postag) ! position frozen
                         ats(p_ind)%constrained(jj) = frozen_pos_index
                         ats(p_ind)%well(jj) = 0.d0
                      CASE(alltag) ! pos + orientation frozen
                         ats(p_ind)%constrained(jj) = all_frozen_index
                         ats(p_ind)%well(jj) = 0.d0
                      CASE(veltag) ! velocity frozen
                         ats(p_ind)%constrained(jj) = frozen_vel_index
                         ats(p_ind)%well(jj) = 0.d0       
                      CASE DEFAULT ! unrecognized tag
                         CALL abort("unrecognized constraint "//tokens(t_ind2)(1:twidths(t_ind2)))               
                      END SELECT
                      
                   CASE(welltag) ! harmonic well
                      t_ind2 = t_ind + 1
                      ats(p_ind)%constrained(jj) = harmonic_well_index
                      ! read spring constant
                      READ(tokens(t_ind2)(1:twidths(t_ind2)),*,IOSTAT=iostat) ats(p_ind)%well(jj)
                      IF(iostat /= 0) CALL abort("reading constraints")
                      IF(ats(p_ind)%well(jj) < 0.d0) &
                           CALL abort("Negative spring constant found in <"//constr_block//">.")
                      
                   CASE(forcetag) ! external force
                      t_ind2 = t_ind + 1
                      ats(p_ind)%constrained(jj) = ext_force_index
                      ! read force strength
                      READ(tokens(t_ind2)(1:twidths(t_ind2)),*,IOSTAT=iostat) ats(p_ind)%well(jj)
                      IF(iostat /= 0) CALL abort("reading constraints")
                      
                   CASE DEFAULT ! unrecognized tag
                      CALL abort("unrecognized constraint "//tokens(t_ind)(1:twidths(t_ind)))
                      
                   END SELECT
                   
                END DO

             END IF
          END DO
       ELSE
          CALL write_miss(tokens(1)(1:twidths(1)),constr_block)
       END IF
    END DO

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_constraints

  
  ! parses information in the "potentials" block.
  ! that is, reads the types of potentials
  ! applied to the given pairs of particles, identified by their labels
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *n_elems number of different types of atoms
  ! *params physical parameters
  ! *found_types numbers of atoms of each type
  SUBROUTINE parse_potentials(data,lwidths,n_elems,params,found_types)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(mbps), INTENT(INOUT) :: params
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:), n_elems
    INTEGER, POINTER :: twidths(:)
    INTEGER, POINTER :: found_types(:)
    INTEGER :: ii, jj, iostat, n_lines, n_tokens, n_pots, t_ind, maxterms, termj
    CHARACTER(LEN=labelw) :: firsttype, secondtype
    LOGICAL, ALLOCATABLE :: interacts(:)
    INTEGER, ALLOCATABLE :: startlines(:), starttypes(:,:)    
    INTEGER :: firstind, secondind
    
    IF(n_elems == 0)THEN ! if there are no atomic species defined, there can be no potentials either
       IF(CPU_ID == MASTER_CPU) WRITE(*,*) "No atomic particles defined, <"//pot_block//"> will be ignored." 
       IF(associated(params%pot_types))THEN
          DEALLOCATE(params%pot_types)
       END IF
       NULLIFY(params%pot_types)
       IF(associated(params%n_pots))THEN
          DEALLOCATE(params%n_pots)
       END IF
       NULLIFY(params%n_pots)
       IF(associated(params%pot_cut))THEN
          DEALLOCATE(params%pot_cut)
       END IF
       NULLIFY(params%pot_cut)
       IF(associated(params%pot_params))THEN
          DEALLOCATE(params%pot_params)
       END IF
       NULLIFY(params%pot_params)
       RETURN
    END IF

    n_lines = SIZE(data(1,:))    

    ! allocate arrays for storing and finding the pair-potential regions in the data
    ALLOCATE(interacts(-1:n_elems))
    ALLOCATE(startlines((n_elems+2)**2+1))
    ALLOCATE(starttypes((n_elems+2)**2,2))
    interacts = .false.
    startlines = 0
    starttypes = -2

    n_pots = 0

    ! find the potential pairs
    DO ii = 1, n_lines ! loop over lines of data
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       
       IF(n_tokens > 1)THEN      
          firsttype = tokens(1)(1:labelw) ! name of the first type in the pair
          secondtype = tokens(2)(1:labelw) ! name of the second type in the pair

          firstind = -2 ! index for unrecognized labels
          secondind = -2
          ! indices for MBs or H-bond arms
          IF(firsttype(1:labelw) == mbname) firstind = 0
          IF(firsttype(1:labelw) == armname) firstind = -1
          IF(secondtype(1:labelw) == mbname) secondind = 0
          IF(secondtype(1:labelw) == armname) secondind = -1
          
          IF(firstind < -1 .OR. secondind < -1)THEN ! there was an unrecognized tag, search the list of elements

             label: DO jj = 1, SIZE(params%atomic_labels(:))
                IF(firsttype(1:labelw) == params%atomic_labels(jj))THEN
                   firstind = jj
                END IF
                IF(secondtype(1:labelw) == params%atomic_labels(jj))THEN
                   secondind = jj
                END IF
                IF(firstind > -2 .AND. secondind > -2) EXIT label
             END DO label
       
          END IF

          ! If both tokens were labels, then we have found the beginning of a new
          ! potential definition
          IF(firstind > -2 .AND. secondind > -2)THEN

             n_pots = n_pots + 1
             startlines(n_pots) = ii ! store the beginning line for the definition of this potential
             starttypes(n_pots,1) = firstind
             starttypes(n_pots,2) = secondind
             interacts(firstind) = .true.
             interacts(secondind) = .true.
             
          END IF

       END IF
    END DO

    ! add this so that we automatically treat the end of data correctly:
    ! potential #n definition starts at startlines(n) and ends at startlines(n+1)-1
    startlines(n_pots+1) = n_lines+1

    ! checks
    IF(n_pots == 0) CALL abort("no valid potentials in <"//pot_block//">")
    DO ii = 1, n_elems
       IF(.NOT.interacts(ii) .AND. found_types(ii) > 0) &
            CALL abort("element "//params%atomic_labels(ii)//" has no interactions")
    END DO

    ! get the max number of terms in a single potential
    maxterms = 0
    DO ii = 1, n_pots ! loop over the found potentials
       IF(maxterms < startlines(ii+1)-startlines(ii)-1) & ! record the maximum numer of terms (each line defines one term)
            maxterms = startlines(ii+1)-startlines(ii)-1
       ! check
       IF(starttypes(ii,1) < 1 .AND. starttypes(ii,2) < 1) &
            CALL abort("additional potentials between MB molecules not allowed")
    END DO

    ! allocate arrays for storing the pair-potential parameters for all pairs
    NULLIFY(params%pot_types)
    ALLOCATE(params%pot_types(maxterms,-1:n_elems,-1:n_elems))
    params%pot_types = 0
    NULLIFY(params%pot_params)
    ALLOCATE(params%pot_params(max_pot_params,maxterms,-1:n_elems,-1:n_elems))
    params%pot_params = 0.d0
    NULLIFY(params%pot_cut)
    ALLOCATE(params%pot_cut(-1:n_elems,-1:n_elems))
    params%pot_cut = 0.d0
    NULLIFY(params%n_pots)
    ALLOCATE(params%n_pots(-1:n_elems,-1:n_elems))
    params%n_pots = 0

    ! get potentials
    DO ii = 1, n_pots ! loop over found potentials
       ! empty ?
       IF(startlines(ii+1) <= startlines(ii)+1) & 
            CALL abort("no data for "//params%atomic_labels(firstind)//" "&
            //params%atomic_labels(secondind)//" in <"//pot_block//">")

       !
       ! read the potential
       !

       ! read the cutoff
       CALL tokenize(data(1:lwidths(startlines(ii)),startlines(ii))," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       IF(n_tokens > 2)THEN
          READ(tokens(3),*,IOSTAT=iostat) params%pot_cut(starttypes(ii,1),starttypes(ii,2))
          IF(iostat /= 0) CALL abort("Reading potential cut-off")
          params%pot_cut(starttypes(ii,2),starttypes(ii,1)) = &
               params%pot_cut(starttypes(ii,1),starttypes(ii,2))
       ELSE ! default cutoff
          params%pot_cut(starttypes(ii,1),starttypes(ii,2)) = params%cut_atom
          params%pot_cut(starttypes(ii,2),starttypes(ii,1)) = params%cut_atom
       END IF

       ! potential terms
       termj = 0
       DO jj = startlines(ii)+1,startlines(ii+1)-1

          termj = termj+1
          CALL tokenize(data(1:lwidths(jj),jj)," ",tokens,twidths)
          n_tokens = SIZE(tokens(:))
          
          ! number of terms in the potential
          params%n_pots(starttypes(ii,1),starttypes(ii,2)) = startlines(ii+1)-startlines(ii)-1
          params%n_pots(starttypes(ii,2),starttypes(ii,1)) = startlines(ii+1)-startlines(ii)-1

          SELECT CASE(tokens(1)(1:twidths(1))) ! potential type
             CASE(LenJon_pot) ! lennard-jones
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                ! store the type
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = lenjon_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = lenjon_index
                DO t_ind = 2, 3 ! parameters
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(exp_pot) ! exponential
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = exp_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = exp_index
                DO t_ind = 2, 3
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(pow_pot) ! power-law
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = pow_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = pow_index
                DO t_ind = 2, 3
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(spring_pot) ! spring
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = spring_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = spring_index
                DO t_ind = 2, 3
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(hard_pot) ! hard
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = hard_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = hard_index
                DO t_ind = 2, 4
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(hardrep_pot) ! repulsive hard
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = hardrep_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = hardrep_index
                DO t_ind = 2, 4
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(fene_pot) ! fene
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = fene_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = fene_index
                DO t_ind = 2, 4
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE(shell_pot) ! shell
                IF(n_tokens < 3) &
                     CALL abort("insufficient data for "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
                params%pot_types(termj,starttypes(ii,1),starttypes(ii,2)) = shell_index
                params%pot_types(termj,starttypes(ii,2),starttypes(ii,1)) = shell_index
                DO t_ind = 2, 4
                   READ(tokens(t_ind)(1:twidths(t_ind)),*,IOSTAT=iostat) &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                   IF(iostat /= 0) CALL abort("reading potentials")
                   params%pot_params(t_ind-1,termj,starttypes(ii,2),starttypes(ii,1)) = &
                        params%pot_params(t_ind-1,termj,starttypes(ii,1),starttypes(ii,2))
                END DO

             CASE DEFAULT
                CALL abort("unrecognized tag "//tokens(1)(1:twidths(1))//" in <"//pot_block//">")
          END SELECT
          

       END DO
    END DO

    DEALLOCATE(interacts)
    DEALLOCATE(startlines)
    DEALLOCATE(starttypes)

    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)

    RETURN
  END SUBROUTINE parse_potentials


  ! parses information in the "statistics" block
  ! *data the contents of the block as a 2D character array
  ! *lwidths widths of the rows in "data"
  ! *control control parameters
  ! *mbs list of molecules
  ! *ats list of atoms
  SUBROUTINE parse_statistics(data,lwidths,control,mbs,ats)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: data(:,:)
    TYPE(cps), INTENT(INOUT) :: control
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, INTENT(IN) :: lwidths(:)
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, jj, kk, ll, iostat, n_lines, n_tokens, firstindex, lastindex, &
         i_stat, n_mbs, n_ats, t_ind, t_ind2, t_ind3, tempindex, maxgroup, groupsize, &
         group_ind, max_n_stats

    n_lines = SIZE(data(1,:))

    ! count number of stat files
    control%n_statfiles = 0
    DO ii = 1, n_lines
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       IF(tokens(1) == writestattag)THEN ! writestat begins the definition of a new stat file
          control%n_statfiles = control%n_statfiles + 1
       END IF
    END DO
    IF(control%n_statfiles > 0)THEN ! if stats are recorded, allocate the arrays for storing the indices and write interval
       ALLOCATE(control%n_stats(control%n_statfiles))
       ALLOCATE(control%stat_start(control%n_statfiles))
       ALLOCATE(control%stat_interval(control%n_statfiles))
    ELSE ! if stats are not recorded, allocate lists of size 1
       ALLOCATE(control%n_stats(1))
       ALLOCATE(control%stat_start(1))
       ALLOCATE(control%stat_interval(1))
    END IF
    control%n_stats = 0
    control%stat_start = 0.d0
    control%stat_interval = ABS(control%verbose)

    ! bonds
    control%bond_writer = xyz_writer_def
    control%bond_interval = missing_r

    ! rdf
    control%rdf_writer = rdf_writer_def
    control%rdf_interval = missing_r
    control%rdf_range = missing_r

    ! adf
    control%adf_writer = rdf_writer_def
    control%adf_interval = missing_r
    control%adf_range = missing_r
    control%adf_axis = z_index
    control%adf_center = 0.d0
    control%adf_follow = 0

    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)
    maxgroup = 0

    group_ind = 0
    ! count the number of stats monitored
    DO ii = 1, n_lines
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))
       IF(tokens(1) == writestattag)THEN ! writestat begins a new group (a new statfile)
          group_ind = group_ind + 1
       END IF

       ! a stat defining tag
       IF(tokens(1)(1:twidths(1)) == timetag .OR. &
            tokens(1)(1:twidths(1)) == clocktag .OR. &
            tokens(1)(1:twidths(1)) == temptag .OR. &
            tokens(1)(1:twidths(1)) == kinetictag .OR. &
            tokens(1)(1:twidths(1)) == potentialtag .OR. &
            tokens(1)(1:twidths(1)) == energytag .OR. &
            tokens(1)(1:twidths(1)) == virialtag .OR. &
            tokens(1)(1:twidths(1)) == lineartag .OR. &
            tokens(1)(1:twidths(1)) == rottag .OR. &
            tokens(1)(1:twidths(1)) == lenjontag .OR. &
            tokens(1)(1:twidths(1)) == hbondtag .OR. &
            tokens(1)(1:twidths(1)) == mbatomtag .OR. &
            tokens(1)(1:twidths(1)) == atomatomtag .OR. &
            tokens(1)(1:twidths(1)) == volumetag .OR. &
            tokens(1)(1:twidths(1)) == pressuretag)THEN

          IF(group_ind == 0)THEN ! no stat writing was requested
             CALL abort("uninitialized stat group in <"//stat_block//">: "//tokens(1)(1:twidths(1)))
          END IF
          control%n_stats(group_ind) = control%n_stats(group_ind) + 1

       ELSE IF(tokens(1)(1:twidths(1)) == intervaltag .OR. &
            tokens(1)(1:twidths(1)) == starttag )THEN ! write-time defining tag

          IF(group_ind == 0)THEN
             CALL abort("uninitialized stat group in <"//stat_block//">: "//tokens(1)(1:twidths(1)))
          END IF 

       ELSE IF(tokens(1)(1:twidths(1)) == coordtag .OR. & ! a stat group defining tag
            tokens(1)(1:twidths(1)) == orientag .OR. &
            tokens(1)(1:twidths(1)) == velotag .OR. &
            tokens(1)(1:twidths(1)) == angulartag .OR. &
            tokens(1)(1:twidths(1)) == forcetag .OR. &
            tokens(1)(1:twidths(1)) == torquetag .OR. &
            tokens(1)(1:twidths(1)) == armstag )THEN
       
          IF(n_tokens > 3)THEN ! count how many particles are monitored - each defines one stat
             IF( is_number(tokens(2)) .AND. tokens(3)(1:twidths(3)) == grouptag .AND. is_number(tokens(4)) )THEN
                READ(tokens(2),*,IOSTAT=iostat) firstindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                IF(lastindex < firstindex) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                control%n_stats = control%n_stats + lastindex - firstindex
             END IF
          END IF
          IF(group_ind == 0)THEN
             CALL abort("uninitialized stat group in <"//stat_block//">: "//tokens(1)(1:twidths(1)))
          END IF 
          control%n_stats(group_ind) = control%n_stats(group_ind) + 1

       ELSE IF(tokens(1)(1:twidths(1)) == forcesumtag)THEN ! a stat defining tag

          ll = 1
          groupsize = 0
          firstindex = 0
          DO WHILE(n_tokens > ll) ! count how many particles are summed to get the maximum sum (group) size
             ll = ll + 1            
             IF(is_number(tokens(ll)))THEN
                groupsize = groupsize + 1
                READ(tokens(ll),*,IOSTAT=iostat) firstindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             ELSE IF(tokens(ll)(1:twidths(ll)) == grouptag)THEN
                IF(n_tokens > ll)THEN
                   ll = ll + 1
                   IF(is_number(tokens(ll)))THEN
                      READ(tokens(ll),*,IOSTAT=iostat) lastindex
                      IF(iostat /= 0 .OR. firstindex <= 0) &
                           CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                      groupsize = groupsize + lastindex - firstindex
                   ELSE
                      CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                   END IF
                END IF
             ELSE
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END DO
          IF(group_ind == 0)THEN
             CALL abort("uninitialized stat group in <"//stat_block//">: "//tokens(1)(1:twidths(1)))
          END IF 
          control%n_stats(group_ind) = control%n_stats(group_ind) + 1
          IF(groupsize > maxgroup) maxgroup = groupsize

       END IF

    END DO
   
    ! maximum number of statistics in any one stat file
    max_n_stats = MAXVAL(control%n_stats(:))

    ! allocate arrays for storing the stat indices
    NULLIFY(control%stats)
    ALLOCATE(control%stats(max_n_stats,control%n_statfiles))
    NULLIFY(control%stat_particles)
    ALLOCATE(control%stat_particles(max_n_stats,control%n_statfiles))
    i_stat = 0
    control%stats = 0
    control%stat_particles = 0
    IF(maxgroup > 0)THEN
       NULLIFY(control%stat_groups)
       ALLOCATE(control%stat_groups(maxgroup,max_n_stats,control%n_statfiles))
       control%stat_groups = 0
    END IF

    group_ind = 0
    ! read the stat tags
    DO ii = 1, n_lines
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       SELECT CASE(tokens(1)(1:twidths(1)))
       CASE(writestattag) ! start a new stat file definition
          group_ind = group_ind + 1
          i_stat = 0

       CASE(intervaltag) ! specify output interval
          IF(n_tokens > 1)THEN
             IF(is_number(tokens(2)))THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) control%stat_interval(group_ind)
                IF(iostat /= 0 .OR. control%stat_interval(group_ind) < 0.d0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF

       CASE(starttag) ! specify output start
          IF(n_tokens > 1)THEN
             IF(is_number(tokens(2)))THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) control%stat_start(group_ind)
                IF(iostat /= 0 .OR. control%stat_interval(group_ind) < 0.d0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF

       CASE(rdftag) ! radial distribution function
          t_ind2 = 2 
          ! writing times
          SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
          CASE(nonetag) ! none
             control%rdf_writer = noxyz_index
             control%rdf_interval = 0.d0
             t_ind3 = t_ind2+1
          CASE(endtag) ! at the end
             control%rdf_writer = exyz_index
             control%rdf_interval = 0.d0
             t_ind3 = t_ind2+1
          CASE(starttag) ! at the start
             control%rdf_writer = sxyz_index
             control%rdf_interval = 0.d0
             t_ind3 = t_ind2+1
          CASE(sendtag) ! at the start and the end
             control%rdf_writer = sexyz_index
             control%rdf_interval = 0.d0
             t_ind3 = t_ind2+1
          CASE(intervaltag) ! periodically once in a defined interval
             control%rdf_writer = ixyz_index
             IF(n_tokens > t_ind2)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%rdf_interval
                IF(iostat /= 0) CALL abort("reading rdf interval")
                t_ind3 = t_ind2+2
             ELSE
                CALL abort_value(rdftag,stat_block)
             END IF
          CASE(averagetag) ! average, at the end
             control%rdf_writer = average_index
             IF(n_tokens > t_ind2+1)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%rdf_interval
                IF(iostat /= 0) CALL abort("reading rdf interval")
                READ(tokens(t_ind2+2),*,IOSTAT=iostat) control%rdf_start
                IF(iostat /= 0)THEN
                   control%rdf_start = 0.d0
                   t_ind3 = t_ind2+2
                ELSE
                   t_ind3 = t_ind2+3
                END IF
             ELSE IF(n_tokens == t_ind2+1)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%rdf_interval
                IF(iostat /= 0) CALL abort("reading rdf interval")
                t_ind3 = t_ind+2
             ELSE
                CALL abort_value(rdftag,stat_block)
             END IF
          CASE DEFAULT ! unrecognized
             CALL abort_value(rdftag,stat_block)
          END SELECT
          IF(n_tokens >= t_ind3)THEN
             ! read the types of particles for which the rdf is measured
             SELECT CASE(tokens(t_ind3)(1:twidths(t_ind3)))
             CASE(alltag)
                control%rdf_particles = all_rdf_index
             CASE(mbname)
                control%rdf_particles = mb_rdf_index
             CASE(atomtag)
                control%rdf_particles = atom_rdf_index
             CASE DEFAULT
                CALL abort_value(rdftag,stat_block)
             END SELECT
          ELSE
             control%rdf_particles = rdf_particles_def
          END IF
          IF(n_tokens > t_ind3)THEN ! range of the rdf
             READ(tokens(t_ind3+1),*,IOSTAT=iostat) control%rdf_range
             IF(iostat /= 0) CALL abort("reading rdf range")
          ELSE
          END IF


       CASE(adftag) ! axial distribution function
          t_ind2 = 2
          ! writing times
          SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
          CASE(nonetag)
             control%adf_writer = noxyz_index
             control%adf_interval = 0.d0
             t_ind = t_ind2+1
          CASE(endtag)
             control%adf_writer = exyz_index
             control%adf_interval = 0.d0
             t_ind = t_ind2+1
          CASE(starttag)
             control%adf_writer = sxyz_index
             control%adf_interval = 0.d0
             t_ind = t_ind2+1
          CASE(sendtag)
             control%adf_writer = sexyz_index
             control%adf_interval = 0.d0
             t_ind = t_ind2+1
          CASE(intervaltag)
             control%adf_writer = ixyz_index
             IF(n_tokens > t_ind2)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%adf_interval
                IF(iostat /= 0) CALL abort("reading adf interval")
                t_ind = t_ind2+2
             ELSE
                CALL abort_value(adftag,stat_block)
             END IF
          CASE(averagetag)                   
             control%adf_writer = average_index
             IF(n_tokens > t_ind2+1)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%adf_interval
                IF(iostat /= 0) CALL abort("reading adf interval")
                READ(tokens(t_ind2+2),*,IOSTAT=iostat) control%adf_start
                IF(iostat /= 0)THEN
                   !CALL abort("reading adf start")
                   control%rdf_start = 0.d0
                   t_ind = t_ind2+2
                ELSE
                   t_ind = t_ind2+3
                END IF
             ELSE
                CALL abort_value(adftag,stat_block)
             END IF
          CASE DEFAULT
             CALL abort_value(adftag,stat_block)
          END SELECT
          ! axis
          IF(control%adf_writer /= noxyz_index)THEN
             IF(n_tokens >= t_ind+2)THEN
                SELECT CASE(tokens(t_ind)(1:twidths(t_ind)))
                CASE(xaxistag)
                   control%adf_axis = x_index
                CASE(yaxistag)
                   control%adf_axis = y_index
                CASE(zaxistag)
                   control%adf_axis = z_index
                CASE DEFAULT
                   CALL abort_value(adftag,stat_block)
                END SELECT
                SELECT CASE(tokens(t_ind+1)(1:twidths(t_ind+1)))
                CASE(angletag)
                   control%adf_angle = .true.
                   t_ind = t_ind+1
                CASE(noangletag)
                   control%adf_angle = .false.
                   t_ind = t_ind+1
                CASE DEFAULT
                   ! no tag, then set false
                   control%adf_angle = .false.
                END SELECT
                SELECT CASE(tokens(t_ind+1)(1:twidths(t_ind+1)))
                CASE(fixedtag) ! fixed axis
                   IF(n_tokens >= t_ind+4)THEN
                      READ(tokens(t_ind+2),*,IOSTAT=iostat) control%adf_center(1)
                      READ(tokens(t_ind+3),*,IOSTAT=iostat) control%adf_center(2)
                      READ(tokens(t_ind+4),*,IOSTAT=iostat) control%adf_center(3)                   
                      t_ind3 = t_ind+5
                   ELSE
                      CALL abort_value(adftag,stat_block)
                   END IF
                CASE(followtag) ! axis follows a particle
                   READ(tokens(t_ind+2),*,IOSTAT=iostat) tempindex
                   t_ind3 = t_ind+3
                   
                   ! find the particle with the corresponding index
                   control%adf_follow = 0
                   findmbA: DO kk = 1, n_mbs
                      IF(mbs(kk)%index == tempindex)THEN
                         control%adf_follow = kk
                         EXIT findmbA
                      END IF
                   END DO findmbA
                   IF(control%adf_follow == 0)THEN
                      findatA: DO kk = 1, n_ats
                         IF(ats(kk)%index == tempindex)THEN
                            control%adf_follow = -kk
                            EXIT findatA
                         END IF
                      END DO findatA
                   END IF
                   IF(control%adf_follow == 0)THEN
                      CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                   END IF
                   
                END SELECT
                
             ELSE
                CALL abort_value(adftag,stat_block)
             END IF
             IF(n_tokens >= t_ind3)THEN
                ! types of particles measured
                SELECT CASE(tokens(t_ind3)(1:twidths(t_ind3)))
                CASE(alltag)
                   control%adf_particles = all_rdf_index
                CASE(mbname)
                   control%adf_particles = mb_rdf_index
                CASE(atomtag)
                   control%adf_particles = atom_rdf_index
                CASE DEFAULT
                   CALL abort_value(adftag,stat_block)
                END SELECT
             ELSE
                control%adf_particles = rdf_particles_def
             END IF
             IF(n_tokens > t_ind3)THEN
                ! adf range
                READ(tokens(t_ind3+1),*,IOSTAT=iostat) control%adf_range
                IF(iostat /= 0) CALL abort("reading adf range")
             ELSE
             END IF
          END IF
          
          ! bond count
       CASE(bondtag)
          t_ind2 = 2
          ! writing times
          SELECT CASE(tokens(t_ind2)(1:twidths(t_ind2)))
          CASE(nonetag)
             control%bond_writer = noxyz_index
             control%bond_interval = 0.d0
          CASE(endtag)
             control%bond_writer = exyz_index
             control%bond_interval = 0.d0
          CASE(starttag)
             control%bond_writer = sxyz_index
             control%bond_interval = 0.d0
          CASE(sendtag)
             control%bond_writer = sexyz_index
             control%bond_interval = 0.d0
          CASE(intervaltag)
             control%bond_writer = ixyz_index
             IF(n_tokens > t_ind2)THEN
                READ(tokens(t_ind2+1),*,IOSTAT=iostat) control%bond_interval
                IF(iostat /= 0) CALL abort("reading bond interval")
             ELSE
                CALL abort_value(bondtag,main_block)
             END IF
          CASE DEFAULT
             CALL abort_value(bondtag,main_block)
          END SELECT
          
          ! stats
       CASE(timetag) ! simulation time
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = time_stat_index
       CASE(clocktag) ! wall clock
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = clock_stat_index
       CASE(temptag) ! temperature
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = temp_stat_index
       CASE(kinetictag) ! kinetic energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = kin_stat_index
       CASE(potentialtag) ! potential energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = pot_stat_index
       CASE(energytag) ! total energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = ene_stat_index
       CASE(virialtag) ! virial
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = virial_stat_index
       CASE(lineartag) ! linear kinetic energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = linear_stat_index
       CASE(rottag) ! rotational kinetic energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = rotational_stat_index
       CASE(lenjontag) ! lennard-jones potential energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = ljpot_stat_index
       CASE(hbondtag) ! h-bond potential energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = hbpot_stat_index
       CASE(mbatomtag) ! mb-atom potential energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = mapot_stat_index
       CASE(atomatomtag) ! atom-atom potential energy
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = aapot_stat_index
       CASE(volumetag) ! volume
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = volume_stat_index
       CASE(pressuretag) ! pressure
          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = press_stat_index
       CASE(coordtag) ! coordinates of a particle
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = coord_stat_index

             findmb1: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb1
                END IF
             END DO findmb1
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                findat1: DO kk = 1, n_ats
                   IF(ats(kk)%index == jj)THEN
                      control%stat_particles(i_stat,group_ind) = -kk
                      EXIT findat1
                   END IF
                END DO findat1
             END IF
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO
       CASE(orientag) ! orientation of an mb
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = orient_stat_index

             findmb2: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb2
                END IF
             END DO findmb2
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO
       CASE(velotag) ! velocity of a particle
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = vel_stat_index

             findmb3: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb3
                END IF
             END DO findmb3
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                findat3: DO kk = 1, n_ats
                   IF(ats(kk)%index == jj)THEN
                      control%stat_particles(i_stat,group_ind) = -kk
                      EXIT findat3
                   END IF
                END DO findat3
             END IF
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO
       CASE(angulartag) ! angular velocity of a particle
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = ang_stat_index

             findmb4: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb4
                END IF
             END DO findmb4
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO

       CASE(armstag) ! h-bond arms of an mb
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = arms_stat_index

             findmb5: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb5
                END IF
             END DO findmb5
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO

       CASE(forcetag) ! force acting on a particle
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = force_stat_index

             findmb6: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb6
                END IF
             END DO findmb6
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                findat6: DO kk = 1, n_ats
                   IF(ats(kk)%index == jj)THEN
                      control%stat_particles(i_stat,group_ind) = -kk
                      EXIT findat6
                   END IF
                END DO findat6
             END IF
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO

       CASE(torquetag) ! torque action on an mb
          READ(tokens(2),*,IOSTAT=iostat) firstindex
          IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
          lastindex = firstindex
          IF(n_tokens > 3)THEN
             IF(is_number(tokens(4)(1:twidths(4))) .AND. tokens(3)(1:twidths(3)) == grouptag)THEN                
                READ(tokens(4),*,IOSTAT=iostat) lastindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END IF
          DO jj = firstindex, lastindex
             i_stat = i_stat + 1
             control%stats(i_stat,group_ind) = torque_stat_index

             findmb7: DO kk = 1, n_mbs
                IF(mbs(kk)%index == jj)THEN
                   control%stat_particles(i_stat,group_ind) = kk
                   EXIT findmb7
                END IF
             END DO findmb7
             IF(control%stat_particles(i_stat,group_ind) == 0)THEN
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF

          END DO

       CASE(forcesumtag) ! sum of forces acting on a group of particles

          i_stat = i_stat + 1
          control%stats(i_stat,group_ind) = forcesum_stat_index

          ll = 1
          groupsize = 0
          firstindex = 0
          DO WHILE(n_tokens > ll)
             ll = ll + 1
             IF(is_number(tokens(ll)))THEN
                groupsize = groupsize + 1
                READ(tokens(ll),*,IOSTAT=iostat) firstindex
                IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)

                control%stat_groups(groupsize,i_stat,group_ind) = 0
                findmb8: DO kk = 1, n_mbs
                   IF(mbs(kk)%index == firstindex)THEN
                      control%stat_groups(groupsize,i_stat,group_ind) = kk
                      EXIT findmb8
                   END IF
                END DO findmb8
                IF(control%stat_groups(groupsize,i_stat,group_ind) == 0)THEN
                   findat8: DO kk = 1, n_ats
                      IF(ats(kk)%index == firstindex)THEN
                         control%stat_groups(groupsize,i_stat,group_ind) = -kk
                         EXIT findat8
                      END IF
                   END DO findat8
                END IF
                IF(control%stat_groups(groupsize,i_stat,group_ind) == 0)THEN
                   CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                END IF

             ELSE IF(tokens(ll)(1:twidths(ll)) == grouptag)THEN
                IF(n_tokens > ll)THEN
                   ll = ll + 1
                   IF(is_number(tokens(ll)))THEN
                      READ(tokens(ll),*,IOSTAT=iostat) lastindex
                      IF(iostat /= 0) CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                      
                      DO jj = firstindex + 1, lastindex
                         groupsize = groupsize + 1
                         control%stat_groups(groupsize,i_stat,group_ind) = 0
                         findmb9: DO kk = 1, n_mbs
                            IF(mbs(kk)%index == jj)THEN
                               control%stat_groups(groupsize,i_stat,group_ind) = kk
                               EXIT findmb9
                            END IF
                         END DO findmb9
                         IF(control%stat_groups(groupsize,i_stat,group_ind) == 0)THEN
                            findat9: DO kk = 1, n_ats
                               IF(ats(kk)%index == jj)THEN
                                  control%stat_groups(groupsize,i_stat,group_ind) = -kk
                                  EXIT findat9
                               END IF
                            END DO findat9
                         END IF
                         IF(control%stat_groups(groupsize,i_stat,group_ind) == 0)THEN
                            CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                         END IF
                      END DO

                   ELSE
                      CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
                   END IF
                END IF
             ELSE
                CALL abort_value(tokens(1)(1:twidths(1)),stat_block)
             END IF
          END DO

       CASE DEFAULT
          CALL write_typo(tokens(1)(1:twidths(1)),stat_block)
       END SELECT
    END DO

    ! validity check and def values
    IF(control%rdf_writer == ixyz_index .AND. control%rdf_interval <= 0.d0 .AND. &
          control%rdf_interval /= missing_r ) CALL abort_value(rdftag,main_block)
    IF(control%rdf_writer == average_index .AND. control%rdf_interval <= 0.d0)THEN
       IF(control%rdf_interval /= missing_r)THEN 
          CALL abort_value(rdftag,main_block)
       ELSE
          control%rdf_interval = control%time_max*0.01d0
       END IF
    END IF
    IF(control%rdf_writer == average_index .AND. control%rdf_start < 0.d0)THEN
       IF(control%rdf_interval /= missing_r)THEN 
          CALL abort_value(rdftag,main_block)
       ELSE
          control%rdf_start = control%time_max*0.5d0
       END IF
    END IF
    IF(control%rdf_writer /= noxyz_index .AND. control%rdf_range <= 0.d0 .AND. &
          control%rdf_range /= missing_r ) CALL abort_value(rdftag,main_block)

    IF(control%bond_writer == ixyz_index .AND. control%bond_interval <= 0.d0) &
         CALL abort_value(bondtag,main_block)

    IF(control%adf_writer == ixyz_index .AND. control%adf_interval <= 0.d0 .AND. &
          control%adf_interval /= missing_r ) CALL abort_value(adftag,main_block)
    IF(control%adf_writer == average_index .AND. control%adf_interval <= 0.d0)THEN
       IF(control%adf_interval /= missing_r)THEN 
          CALL abort_value(adftag,main_block)
       ELSE
          control%adf_interval = control%time_max*0.01d0
       END IF
    END IF
    IF(control%adf_writer == average_index .AND. control%adf_start < 0.d0)THEN
       IF(control%adf_interval /= missing_r)THEN 
          CALL abort_value(adftag,main_block)
       ELSE
          control%adf_start = control%time_max*0.5d0
       END IF
    END IF
    IF(control%adf_writer /= noxyz_index .AND. control%adf_range <= 0.d0 .AND. &
          control%adf_range /= missing_r ) CALL abort_value(adftag,main_block)

    ! no stats, only rdf, bonds etc.
    IF(control%n_statfiles == 0)THEN
       control%stat_start = control%time_max*2.d0
       control%stat_interval = control%time_max*2.d0
    END IF

    RETURN
  END SUBROUTINE parse_statistics


  ! removes unused atom types.
  ! if "particles" or "potentials" contain types that have not been
  ! defined in "elements" (besides mb and hb), the program aborts. if there are types 
  ! defined in "elements" not entering the simulation, this is no problem but there will
  ! be unnecessary elements in some of the parameter tables.
  ! also the number of different types will be wrong and it may cause complications
  ! especially if there are no atoms and the "n_elems" variable suggests there are.
  ! this routine drops the unused parts.
  ! *n_elems number of different types of atoms
  ! *found_types numbers of atoms of each type
  ! *params physical parameters  
  SUBROUTINE drop_unused(n_elems,found_types,params)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: n_elems
    TYPE(mbps), INTENT(INOUT) :: params
    REAL(KIND=dp), ALLOCATABLE :: temp_masses(:), temp_potpars(:,:,:,:), temp_potcut(:,:)
    INTEGER, ALLOCATABLE :: temp_pottypes(:,:,:), temp_npots(:,:), temp_found(:)
    CHARACTER(LEN=labelw), ALLOCATABLE :: temp_labels(:)
    INTEGER, POINTER :: found_types(:)
    INTEGER :: ii, jj, kk, ll, n_found
    LOGICAL :: foundlabel

    IF(n_elems > 0)THEN
       foundlabel = .true. 
       n_found = 0
       DO ii = 1, n_elems ! check which elements are really used in the simulation
          IF(found_types(ii) == 0)THEN
             foundlabel = .false.
          ELSE
             n_found = n_found+1
          END IF
       END DO
       IF(.NOT.foundlabel)THEN ! there is an element defined which does not enter the simulation
          ! allocate temporary arrays for the reduced element list
          ALLOCATE(temp_found(n_found))
          ALLOCATE(temp_labels(n_found))
          ALLOCATE(temp_masses(n_found))
          ALLOCATE(temp_npots(-1:n_found,-1:n_found))
          ALLOCATE(temp_potcut(-1:n_found,-1:n_found))
          ALLOCATE(temp_pottypes(SIZE(params%pot_types(:,1,1)),-1:n_found,-1:n_found))
          ALLOCATE(temp_potpars(SIZE(params%pot_params(:,1,1,1)),SIZE(params%pot_params(1,:,1,1)),&
               -1:n_found,-1:n_found))
          temp_npots = 0
          temp_pottypes = 0
          temp_potcut = 0.d0
          temp_potpars = 0.d0

          kk = 0
          DO ii = 1, n_elems ! copy the data for the elements which really are simulated
             IF(found_types(ii) /= 0)THEN
                kk = kk + 1
                temp_found(kk) = found_types(ii)
                temp_labels(kk) = params%atomic_labels(ii)
                temp_masses(kk) = params%m_atoms(ii)
                temp_npots(-1:0,kk) = params%n_pots(-1:0,ii)
                temp_npots(kk,-1:0) = params%n_pots(ii,-1:0)
                temp_potcut(-1:0,kk) = params%pot_cut(-1:0,ii)
                temp_potcut(kk,-1:0) = params%pot_cut(ii,-1:0)
                temp_pottypes(:,-1:0,kk) = params%pot_types(:,-1:0,ii)
                temp_pottypes(:,kk,-1:0) = params%pot_types(:,ii,-1:0)                
                temp_potpars(:,:,-1:0,kk) = params%pot_params(:,:,-1:0,ii)
                temp_potpars(:,:,kk,-1:0) = params%pot_params(:,:,ii,-1:0)
                ll = 0
                DO jj = 1, n_elems
                   IF(found_types(jj) /= 0)THEN
                      ll = ll + 1
                      temp_npots(ll,kk) = params%n_pots(ll,ii)
                      temp_npots(kk,ll) = params%n_pots(ii,ll)
                      temp_potcut(ll,kk) = params%pot_cut(ll,ii)
                      temp_potcut(kk,ll) = params%pot_cut(ii,ll)
                      temp_pottypes(:,ll,kk) = params%pot_types(:,ll,ii)
                      temp_pottypes(:,kk,ll) = params%pot_types(:,ii,ll)                
                      temp_potpars(:,:,ll,kk) = params%pot_params(:,:,ll,ii)
                      temp_potpars(:,:,kk,ll) = params%pot_params(:,:,ii,ll)
                   END IF
                END DO
             END IF
          END DO          
          
          ! transfer the temporary array data to the real arrays
          IF(associated(found_types))THEN
             DEALLOCATE(found_types)
          END IF
          NULLIFY(found_types)
          IF(associated(params%atomic_labels))THEN
             DEALLOCATE(params%atomic_labels)
          END IF
          NULLIFY(params%atomic_labels)
          IF(associated(params%m_atoms))THEN
             DEALLOCATE(params%m_atoms)
          END IF
          NULLIFY(params%m_atoms)
          IF(associated(params%n_pots))THEN
             DEALLOCATE(params%n_pots)
          END IF
          NULLIFY(params%n_pots)
          IF(associated(params%pot_types))THEN
             DEALLOCATE(params%pot_types)
          END IF
          NULLIFY(params%pot_types)
          IF(associated(params%pot_cut))THEN
             DEALLOCATE(params%pot_cut)
          END IF
          NULLIFY(params%pot_cut)
          IF(associated(params%pot_params))THEN
             DEALLOCATE(params%pot_params)
          END IF
          NULLIFY(params%pot_params)
          ALLOCATE(found_types(n_found))
          ALLOCATE(params%atomic_labels(n_found))
          ALLOCATE(params%m_atoms(n_found))
          ALLOCATE(params%n_pots(-1:n_found,-1:n_found))
          ALLOCATE(params%pot_cut(-1:n_found,-1:n_found))
          ALLOCATE(params%pot_types(SIZE(temp_pottypes(:,1,1)),-1:n_found,-1:n_found))
          ALLOCATE(params%pot_params(SIZE(temp_potpars(:,1,1,1)),SIZE(temp_potpars(1,:,1,1)),&
               -1:n_found,-1:n_found))
          found_types(:) = temp_found(:)
          params%atomic_labels(:) = temp_labels(:)
          params%m_atoms(:) = temp_masses(:)
          params%n_pots(:,:) = temp_npots(:,:)
          params%pot_cut(:,:) = temp_potcut(:,:)
          params%pot_types(:,:,:) = temp_pottypes(:,:,:)
          params%pot_params(:,:,:,:) = temp_potpars(:,:,:,:)
          DEALLOCATE(temp_found)
          DEALLOCATE(temp_labels)
          DEALLOCATE(temp_masses)
          DEALLOCATE(temp_npots)
          DEALLOCATE(temp_potcut)
          DEALLOCATE(temp_pottypes)
          DEALLOCATE(temp_potpars)

          n_elems = n_found

       END IF       
    END IF

  END SUBROUTINE drop_unused


  ! Writes a message saying a tag was not recognized.
  ! This is most likely caused by a typo in the input file.
  ! *tag the unrecognized tag
  ! *inblock the input data block where the tag was found
  SUBROUTINE write_typo(tag,inblock)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: tag, inblock
    IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Unrecognized tag "//tag//" in <"//inblock//">"
    IF(CPU_ID == MASTER_CPU) WRITE(*,*) "   ignoring..."

    RETURN
  END SUBROUTINE write_typo


  ! Writes a message saying a tag is missing its value.
  ! *tag the valueless tag
  ! *inblock the input data block where the tag was found
  SUBROUTINE write_miss(tag,inblock)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: tag, inblock

    IF(CPU_ID == MASTER_CPU) WRITE(*,*) "Tag "//tag//" in <"//inblock//"> is missing its value"
    IF(CPU_ID == MASTER_CPU) WRITE(*,*) "   ignoring..."

    RETURN
  END SUBROUTINE write_miss

  ! Aborts program execution due to an invalid or missing parameter for a required tag.
  ! Before quitting, writes a message stating the invalid parameter.
  ! *tag the invalid tag
  ! *inblock the input data block where the tag was found  
  SUBROUTINE abort_value(tag,inblock)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: tag, inblock
    CALL abort("Parameter "//tag//" in <"//inblock//"> invalid or missing")

  END SUBROUTINE abort_value

  ! Tests whether a string is a number or not
  ! *str the string to be tested
  ! *digit true if "str" is a number
  FUNCTION is_number(str) &
       RESULT(digit)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: digit
    INTEGER :: n_char, ii

    n_char = LEN_TRIM(str)
    digit = .TRUE.

    DO ii = 1, n_char
       IF(str(ii:ii) /= "0" .AND. &
            str(ii:ii) /= "1" .AND. &
            str(ii:ii) /= "2" .AND. &
            str(ii:ii) /= "3" .AND. &
            str(ii:ii) /= "4" .AND. &
            str(ii:ii) /= "5" .AND. &
            str(ii:ii) /= "6" .AND. &
            str(ii:ii) /= "7" .AND. &
            str(ii:ii) /= "8" .AND. &
            str(ii:ii) /= "9" .AND. &
            str(ii:ii) /= "." .AND. &
            str(ii:ii) /= " ") THEN
          digit = .FALSE.
          EXIT
       END IF
    END DO

    RETURN
  END FUNCTION is_number

  ! Changes a 1D character array (up to the given length) to a string.
  ! (There must be an intrinsic command to do the same...?)
  ! *charar the array to be read
  ! *length of the produced string
  ! *str string representation of the character array
  FUNCTION to_string(charar,length) &
       RESULT(str)
    INTEGER, INTENT(IN) :: length
    CHARACTER(LEN=length) :: str
    CHARACTER, INTENT(IN) :: charar(:)
    INTEGER :: readlength, ii

    readlength = MIN(length,SIZE(charar(:)))
    str = ""

    DO ii = 1, readlength
       str(ii:ii) = charar(ii)
    END DO

    RETURN
  END FUNCTION to_string

  ! Reads an entire ASCII text file and returns the read characters as a
  ! 2D character array. Each line in the array has the same length,
  ! so the routine also gives a list of integers stating the 
  ! lengths of the read lines. (the characters on a line beyond this value
  ! are thus all spaces).
  ! *filename the name of the file to be read - the string must be of the right length
  ! *lines the read file as a 2D character array (first index = line, second index = position in the line)
  ! *linelength integer list stating the information lengths of the read lines
  SUBROUTINE read_lines(filename,lines,linelength)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER, POINTER :: lines(:,:)
    INTEGER, POINTER :: linelength(:)
    INTEGER :: ioindex, iostat
    INTEGER :: lcount, wcount, maxw
    CHARACTER :: letter
    LOGICAL :: continue

    ioindex = 86075

    IF(cpu_id == master_cpu)THEN ! only the master reads

       ! open the input
       OPEN(ioindex,FILE=filename,IOSTAT=iostat)
       
       lcount = 1 ! counts the number of lines
       wcount = 0 ! counts the number of characters on a line
       maxw = 0   ! the maximum wcount of all lines
       iostat = 0 ! io status variable
       continue = .TRUE.
       ! read the file, a letter at a time, as long as the file continues
       DO WHILE(continue)
          letter = "e" ! the read letter is stored in "letter"
          READ(ioindex,'(A1)',IOSTAT=iostat,advance='no') letter
          IF(letter == "\t")THEN ! a tab-character is changed to a single space
             letter = " "
          END IF
          IF(iostat == 0)THEN ! if the reading was succesful, the length of the line in increased
             wcount = wcount+1
          ELSE
             !IF(letter /= "e")THEN
             IF(iostat /= -1)THEN ! end of a line is met
                lcount = lcount + 1 ! increase the line count by one
                IF(wcount > maxw) maxw = wcount ! if the width of the line in the largest yet, store the new maximum
                wcount = 0 ! set wcount back to zero, for the new line
             ELSE
                continue = .FALSE. ! if iostat is something else, we have most likely found the end of the file
             END IF
          END IF
       END DO
       
       ! close
       CLOSE(ioindex)
       
       ! allocate the character array for the data according to the measured line count and maximum width
       ALLOCATE(lines(maxw,lcount), STAT = iostat)
       IF(iostat /= 0) CALL abort("too big input file?")
       ALLOCATE(linelength(lcount), STAT = iostat)
       IF(iostat /= 0) CALL abort("too big input file?")
       
    END IF

    ! broadcast the array size to all cpus
#ifdef MPI
    CALL MPI_BCAST(maxw,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    CALL MPI_BCAST(lcount,1,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
#endif

    IF(cpu_id == master_cpu)THEN ! read the data in the character array with the master cpu
       
       lines = " "
       linelength = 0
       
       OPEN(ioindex,FILE=filename,IOSTAT=iostat)
       
       iostat = 0
       lcount = 1
       wcount = 1
       continue = .TRUE.
       DO WHILE(continue)
          letter = "e"       
          READ(ioindex,'(A1)',IOSTAT=iostat,advance='no') letter
          IF(letter == "\t")THEN ! tab to space
             letter = " "
          END IF
          IF(iostat == 0)THEN
             lines(wcount,lcount) = letter ! store the letter
             wcount = wcount+1          
          ELSE
             !IF(letter /= "e")THEN
             IF(iostat /= -1)THEN
                linelength(lcount) = wcount-1
                lcount = lcount + 1
                wcount = 1
             ELSE
                continue = .FALSE.
             END IF
          END IF
       END DO
       
       CLOSE(ioindex)

    ELSE ! in other cpus, allocate the data arrays while the master is reading

       ALLOCATE(lines(maxw,lcount), STAT = iostat)
       IF(iostat /= 0) CALL abort("too big input file?")
       ALLOCATE(linelength(lcount), STAT = iostat)
       IF(iostat /= 0) CALL abort("too big input file?")
       
    END IF

    ! broadcast the read data arrays to all cpus
#ifdef MPI
    CALL MPI_BCAST(lines,SIZE(lines),MPI_CHARACTER,master_cpu,MPI_COMM_WORLD,mpistat)
    CALL MPI_BCAST(linelength,SIZE(linelength),MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
#endif    

    RETURN
  END SUBROUTINE read_lines

  ! Given a set of lines as a 2D character array and the assumed 
  ! lengths of the lines, the routine discards all spaces from the
  ! beginning and end of each line and also removes empty lines.
  ! If a comment character is given, the lines will be also cut
  ! if the character is met and the part to the right of the comment
  ! will be ignored.
  ! *lines the 2D character array to be trimmed
  ! *linelength the information lengths of the lines
  ! *comment an optional comment character after which all information on a line will be ignored
  SUBROUTINE trim(lines,linelength,comment)
    IMPLICIT NONE
    CHARACTER, POINTER :: lines(:,:) 
    CHARACTER, ALLOCATABLE ::templines(:,:)
    INTEGER, POINTER :: linelength(:)
    INTEGER, ALLOCATABLE :: templength(:), inispace(:)
    CHARACTER, INTENT(IN), OPTIONAL :: comment
    INTEGER :: ii, jj, n_lines, allostat, emptylines, maxw
    LOGICAL, ALLOCATABLE :: is_empty(:)
    LOGICAL :: skip

    n_lines = SIZE(lines(1,:))
    ! array for the new lengths of the data lines after trimming
    ALLOCATE(templength(n_lines),STAT=allostat)
    IF(allostat /= 0) CALL abort("handling input data")
    ! array for marking lines that are empty
    ALLOCATE(is_empty(n_lines), STAT=allostat)
    IF(allostat /= 0) CALL abort("handling input data")
    ! array for counting the number of spaces at the beginning of each (non-empty) line
    ALLOCATE(inispace(n_lines), STAT=allostat)
    IF(allostat /= 0) CALL abort("handling input data")

    templength = 0
    inispace = 0
    emptylines = 0 ! counter for number of empty lines
    is_empty = .false.
    DO ii = 1, n_lines ! loop over lines of data

       IF(linelength(ii) > 0)THEN ! some characters are on the line
          
          skip = .false.
          IF(present(comment))THEN ! a comment character is defined
             templength(ii-emptylines) = linelength(ii) ! by default, the linelength remains - but the line index may change
             comm: DO jj = 1, linelength(ii) ! search for a comment character and cut everything after it
                IF(lines(jj,ii) == comment)THEN
                   IF(jj > 1)THEN ! there is something besides the comments
                      templength(ii-emptylines) = jj-1 ! the new, cut length
                      skip = .true.
                      EXIT comm
                   ELSE ! only comments, so the line is empty
                      is_empty(ii) = .true.
                      emptylines = emptylines + 1
                      skip = .true.
                      EXIT comm
                   END IF
                END IF
             END DO comm
          ELSE ! no comments defined
             templength(ii-emptylines) = linelength(ii)
          END IF
          IF(.NOT.is_empty(ii))THEN ! not an empty line
             space: DO jj = templength(ii-emptylines), 1, -1 ! run the line backwards
                IF(lines(jj,ii) /= " ")THEN ! if there are spaces at the end of the line, cut the lenght to drop them
                   templength(ii-emptylines) = jj
                   EXIT space ! end the backwards reading
                END IF
                IF(jj == 1)THEN ! we read through the whole line and found only spaces, it's an empty line
                   emptylines = emptylines + 1
                   is_empty(ii) = .true.
                END IF
             END DO space
          END IF
          IF(.NOT.is_empty(ii))THEN ! still not an empty line
             inis: DO jj = 1, linelength(ii) ! read the line from the beginning
                IF(lines(jj,ii) /= "")THEN ! if there are spaces in the beginning of the line, store the position of the first non-space character
                   inispace(ii) = jj-1
                   templength(ii-emptylines) = templength(ii-emptylines) - (jj-1) ! cut the length by the number of initial spaces
                   EXIT inis
                END IF
             END DO inis
          END IF
       ELSE ! an empty line
          is_empty(ii) = .true.
          emptylines = emptylines + 1
       END IF
       
    END DO

    maxw = 0
    DO ii = 1, n_lines-emptylines ! get the new maximum length of all non-empty lines
       IF(templength(ii) > maxw) maxw = templength(ii)
    END DO

    ! array for the new trimmed data
    ALLOCATE(templines(maxw,n_lines-emptylines), STAT=allostat)
    IF(allostat /= 0) CALL abort("handling input data")
    templines = " "

    emptylines = 0
    DO ii = 1, n_lines ! loop over the lines of data
       
       IF(is_empty(ii))THEN ! skip empty lines
          emptylines = emptylines + 1 ! count empty lines
       ELSE ! not an empty line - store it
          templines(1:templength(ii-emptylines),ii-emptylines) = & !-inispace(ii)) = &
               lines(1+inispace(ii):templength(ii-emptylines)+inispace(ii),ii) ! if there were spaces at the start of the line, shift the data
       END IF

    END DO

    ! new arrays for lines and linelengths
    IF(associated(lines))THEN
       DEALLOCATE(lines)
    END IF
    NULLIFY(lines)
    IF(associated(linelength))THEN
       DEALLOCATE(linelength)
    END IF
    NULLIFY(linelength)
    ALLOCATE(lines(maxw,n_lines-emptylines))
    ALLOCATE(linelength(n_lines-emptylines))

    ! copy the temporary arrays to the real ones
    lines(1:maxw,1:n_lines-emptylines) = &
         templines(1:maxw,1:n_lines-emptylines)
    linelength(1:n_lines-emptylines) = &
         templength(1:n_lines-emptylines)

    DEALLOCATE(is_empty)
    DEALLOCATE(templines)
    DEALLOCATE(templength)
    DEALLOCATE(inispace)

    RETURN
  END SUBROUTINE trim

  ! Given a line of text as a 1D char array and a separator character,
  ! the subroutine splits the line into tokes at any found separator.
  ! The separator character itself will not be included in the tokens and
  ! several separators in a row will be treated as one.
  ! The tokens are returned as strings, and a list with the lengths
  ! of the strings is also given. <br />
  ! For example: the array "hello  world" tokenized using separator " "
  ! will return the tokens "hello" and "world".
  ! *line 1D character array representation of a line of text
  ! *separator the character according to which the the line is split
  ! *tokens list of tokens, obtained by splitting "line" in string format
  ! *tokenwidths list of the lengths of the tokens
  SUBROUTINE tokenize(line,separator,tokens,tokenwidths)
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: line(:)
    CHARACTER, INTENT(IN) :: separator
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, POINTER :: tokenwidths(:)
    INTEGER :: ii, width, n_tokens, maxw, ts, te, allostat
    LOGICAL :: gap

    width = SIZE(line(:))
    gap = .true. ! marker for tracking a cluster of several separators, true when the last character was a separator
    n_tokens = 0
    maxw = 0
    ts = 0
    te = 0
    ! count tokens, i.e., count how many times does the separator appear
    DO ii = 1, width
       IF(line(ii) == separator)THEN ! a separator character was found
          IF(gap)THEN
             ! several separators, treat as one
          ELSE
             ! a new separator, the token ends
             n_tokens = n_tokens+1
             te = ii-1 ! token end index
             IF(maxw < te-ts+1) maxw = te-ts+1 ! max width
          END IF
          gap = .true. ! to mark that the last character was a separator
       ELSE
          IF(gap)THEN
             ! start of a new token
             ts = ii ! token start index
          ELSE
             ! token continues
          END IF
          IF(ii == width)THEN
             ! the line ends, so does the token
             n_tokens = n_tokens + 1
             te = ii ! token end index
             IF(maxw < te-ts+1) maxw = te-ts+1 ! max width
          END IF
          gap = .false. ! to mark that the last character was not a separator
       END IF
    END DO

    ! allocate arrays for the tokens
    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(tokenwidths))THEN
       DEALLOCATE(tokenwidths)
    END IF
    NULLIFY(tokenwidths)
    ALLOCATE(tokens(n_tokens),STAT=allostat)
    ALLOCATE(tokenwidths(n_tokens),STAT=allostat)

    tokens = ""
    tokenwidths = 0
    gap = .true.
    n_tokens = 0
    ts = 0
    te = 0
    ! read tokens
    DO ii = 1, width
       IF(line(ii) == separator)THEN
          IF(gap)THEN
             ! several separators, treat as one
          ELSE
             ! a new separator, the token ends
             n_tokens = n_tokens+1
             te = ii-1 ! token end index
             ! store the token and its width
             tokens(n_tokens) = to_string(line(ts:te),te-ts+1)
             tokenwidths(n_tokens) = te-ts+1
          END IF
          gap = .true.
       ELSE
          IF(gap)THEN
             ! start of a new token
             ts = ii ! token start index
          ELSE
             ! token continues
          END IF
          IF(ii == width)THEN
             ! the line ends, so does the token
             n_tokens = n_tokens + 1
             te = ii ! token end index
             ! store the token and its width
             tokens(n_tokens) = to_string(line(ts:te),te-ts+1)
             tokenwidths(n_tokens) = te-ts+1
             END IF
          gap = .false.
       END IF
    END DO

    RETURN
  END SUBROUTINE tokenize

  ! Changes uppercase letters to lowercase (english alphabet only)
  ! in a 2D character array. This is used to change the entire
  ! input file to lowercase so that parsing would not be case sensitive.
  ! *text the text to be lowecased, as a 2D character array
  SUBROUTINE lowercase(text)
    IMPLICIT NONE
    CHARACTER, POINTER :: text(:,:)
    INTEGER :: ii, jj, ni, nj

    ni = SIZE(text(:,1))
    nj = SIZE(text(1,:))

    DO ii = 1,ni
       DO jj = 1,nj
          SELECT CASE(text(ii,jj))
          CASE("A")
             text(ii,jj) = "a"
          CASE("B")
             text(ii,jj) = "b"
          CASE("C")
             text(ii,jj) = "c"
          CASE("D")
             text(ii,jj) = "d"
          CASE("E")
             text(ii,jj) = "e"
          CASE("F")
             text(ii,jj) = "f"
          CASE("G")
             text(ii,jj) = "g"
          CASE("H")
             text(ii,jj) = "h"
          CASE("I")
             text(ii,jj) = "i"
          CASE("J")
             text(ii,jj) = "j"
          CASE("K")
             text(ii,jj) = "k"
          CASE("L")
             text(ii,jj) = "l"
          CASE("M")
             text(ii,jj) = "m"
          CASE("N")
             text(ii,jj) = "n"
          CASE("O")
             text(ii,jj) = "o"
          CASE("P")
             text(ii,jj) = "p"
          CASE("Q")
             text(ii,jj) = "q"
          CASE("R")
             text(ii,jj) = "r"
          CASE("S")
             text(ii,jj) = "s"
          CASE("T")
             text(ii,jj) = "t"
          CASE("U")
             text(ii,jj) = "u"
          CASE("V")
             text(ii,jj) = "v"
          CASE("W")
             text(ii,jj) = "w"
          CASE("X")
             text(ii,jj) = "x"
          CASE("Y")
             text(ii,jj) = "y"
          CASE("Z")
             text(ii,jj) = "z"
          END SELECT
       END DO
    END DO

    RETURN
  END SUBROUTINE lowercase

  ! Writes time usage statistics to the main output file.
  ! Meant to be called at the end of the simulation.
  ! *ioindex index for the output file (must be open already)
  ! *oldt the time vector for the beginning of the simulation
  ! *time_io cpu time counter for input/output
  ! *time_force cpu time counter for force calculations
  ! *time_stat cpu time counter for statistics calculations
  ! *time_move cpu time counter for system updating
  ! *time_nbor cpu time counter for neighbor finding
  ! *time_tot total cpu time counter  
  SUBROUTINE write_end(ioindex,oldt,&
       time_io, time_force, time_stat, time_nbor, time_move, time_tot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ioindex, oldt(8)
    INTEGER :: tvector(8), dt(6)
    CHARACTER(LEN=20) :: thetime
    REAL(KIND=dp), INTENT(IN) :: time_io, time_force, time_stat, time_nbor, time_move, time_tot
        
    ! get the time
    CALL DATE_AND_TIME(values=tvector)
    CALL date_string(tvector(1),tvector(2),tvector(3),tvector(5),tvector(6),tvector(7),thetime)

    IF(cpu_id == master_cpu)THEN
       WRITE(ioindex,'(A)') "\n"//HORILINE//"\n" ! line
       WRITE(ioindex,'(A)') "  Simulation finished   "//thetime//"\n" ! time
    END IF

    ! time difference from start
    thetime = ""
    CALL time_difference(oldt,tvector,dt)
    CALL time_string(dt(1),dt(2),dt(3),dt(4),dt(5),thetime)

    IF(cpu_id == master_cpu)THEN
       ! time consumption breakdown
       WRITE(ioindex,'(A)') "  Wall clock time "//thetime//"\n"
       WRITE(ioindex,'(A,F20.3,A)') "  Used CPU time:         ",time_tot," s"
       WRITE(ioindex,'(A,F20.3,A)') "     input/output:       ",time_io," s"
       WRITE(ioindex,'(A,F20.3,A)') "     analysis:           ",time_stat," s"
       WRITE(ioindex,'(A,F20.3,A)') "     finding neighbors:  ",time_nbor," s"
       WRITE(ioindex,'(A,F20.3,A)') "     calculating forces: ",time_force," s"
       WRITE(ioindex,'(A,F20.3,A)') "     updating system:    ",time_move," s"
       WRITE(ioindex,'(A)') "\n"
    END IF

    RETURN
  END SUBROUTINE write_end

  ! Writes a valid input file according to the current state of the simulated
  ! system. Meant for writing continuation input files. The actual formatting
  ! and writing of the file is done using "write_input", but this routine
  ! opens the file for output.
  ! *filename the name of the continuation file - the string must be of the right length
  ! *control control parameters
  ! *params physical parameters
  ! *cell supercell dimensions
  ! *pbc true for periodic boundaries
  ! *btype boundary type index
  ! *bval boundary parameter
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_elems number of different types of atoms
  ! *found_types the number of atoms of each type
  ! *is_constr true if there are constraints in place
  SUBROUTINE write_contfile(filename,mbs,ats,n_elems,&
       params,control,cell,btype,bval,is_constr)
    IMPLICIT NONE
    INTEGER :: ioindex
    CHARACTER(LEN=*) :: filename
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3)
    LOGICAL, INTENT(IN) :: is_constr
    INTEGER, INTENT(IN) :: n_elems, btype(3)
    !INTEGER, POINTER :: found_types(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control

    IF(cpu_id == master_cpu)THEN
       ! open the output
       ioindex = 3965
       OPEN(ioindex,FILE=filename)
       
       ! write by calling write_input
       CALL write_input(ioindex,mbs,ats,n_elems,&
            params,control,cell,btype,bval,is_constr)
       
       ! close the output
       CLOSE(ioindex)
    END IF

  END SUBROUTINE write_contfile


  ! Writes the header of the main output file.
  ! The header contains start time information
  ! as well as the parameters in use formatted 
  ! as an input file.
  ! *ioindex index for the output file (must be already open)
  ! *system name of the simulation
  ! *control control parameters
  ! *params physical parameters
  ! *cell supercell dimensions
  ! *pbc true for periodic boundaries
  ! *btype boundary type index
  ! *bval boundary parameter
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_elems number of different types of atoms
  ! *found_types the number of atoms of each type
  ! *is_constr true if there are constraints in place
  ! *tvector simulation start time as given by "date_and_time"
  ! *e_pot total potential energy
  ! *e_lenjon MB-MB lennard-jones energy
  ! *e_bonds MB-MB hyrogen bond energy
  ! *e_mbat MB-atom potential energy
  ! *e_atat atom-atom potential energy
  ! *e_constr potential energy due to constraining potentials/forces
  ! *e_kin kinetic energy
  ! *temper instant temperature
  ! *m_force forces acting on MB molecules
  ! *m_torque torques acting on MB molecules
  ! *a_force forces acting on atoms
  SUBROUTINE write_header(ioindex,system,mbs,ats,n_elems,found_types,tvector,&
       params,control,cell,btype,bval,pbc,is_constr,&
       e_kin,e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,temper,&
       m_force,m_torque,a_force,&
       mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,aa_nbors,aa_n_nbor)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ioindex
    CHARACTER(LEN=*) :: system
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(IN) :: cell(3),bval(3),&
         e_kin,e_pot,e_lenjon,e_bonds,e_mbat,e_atat,e_constr,temper,&
         m_force(:,:), m_torque(:,:), a_force(:,:)
    LOGICAL, INTENT(IN) :: is_constr, pbc(3)
    INTEGER, INTENT(IN) :: tvector(8), n_elems, btype(3)
    INTEGER, POINTER :: found_types(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_nbor(:), &
         ma_nbors(:,:), ma_n_nbor(:), aa_nbors(:,:), aa_n_nbor(:)
    CHARACTER(LEN=20) :: thetime
    INTEGER :: ii, jj, nj, kk, n_mbs, n_ats, closest, closelist(5), listnumber
    REAL(KIND=dp) :: closedist, rangedist(5), testdist, vec(3)

    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    IF(cpu_id == master_cpu)THEN
       WRITE(ioindex,'(A)') & ! title
         "*******************************************************************\n"//&
         "*                         CASHEW "//version//"                        *\n"//&
         "* Coarse Approach Simulator for Hydrogen-bonding Effects in Water *\n"//&
         "*                                                                 *\n"//&
         "*          Dias et al., J. Chem. Phys. 131 (2009) 054505          *\n"//&
         "*                                                                 *\n"//&
         "*                      Teemu Hynninen 2009                        *\n"//&
         "*******************************************************************\n"

       ! time
    CALL date_string(tvector(1),tvector(2),tvector(3),tvector(5),tvector(6),tvector(7),thetime)
    WRITE(ioindex,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)') "  Simulation started   "//thetime//"\n\n"

    ! basic info
    IF(control%verblevel >= 1)THEN
       WRITE(ioindex,'(A)') &
            "Name of simulation: "//system//"\n\n"//& ! name
            "Input read from file "//system//"."//MB_IN//"\n" ! input
#ifdef MPI
       WRITE(ioindex,'(A,I4,A)') &
            "Running on ", n_cpus, " processors\n" ! number of cpus
#else
       WRITE(ioindex,'(A)') &
            "Running the serial version\n" ! number of cpus = 1
#endif
       WRITE(ioindex,'(A)') &
            "Simulation parameters in use: \n"//HORILINE//"\n" ! line

       ! write the current set of input parameters
       CALL write_input(ioindex,mbs,ats,n_elems,&
            params,control,cell,btype,bval,is_constr)
    ELSE ! lower verbosity level
       WRITE(ioindex,'(A)') &
            "Name of simulation: "//system//"\n\n"//&
            "Input read from file "//system//"."//MB_IN//"\n"
#ifdef MPI
       WRITE(ioindex,'(A,I4,A)') &
            "Running on ", n_cpus, " processors\n"
#else
       WRITE(ioindex,'(A)') &
            "Running the serial version\n" ! number of cpus = 1
#endif
    END IF

    
    IF(control%verblevel >= 2)THEN ! write the initial close neighbors lists
       WRITE(ioindex,'(A)') "\n"
       WRITE(ioindex,'(A)') HORILINE//"\n"
       WRITE(ioindex,'(A)') "Close neighbors: \n"
       
       ! generate lists of near neighbors
       DO ii = 1, n_mbs
          closest = 0
          closedist = (cell(1)+cell(2)+cell(3))**2
          DO jj = 1, mm_n_nbor(ii)
             nj = mm_nbors(jj,ii)
             vec = vector(mbs(ii),mbs(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < closedist)THEN
                closest = mbs(nj)%index
                closedist = testdist
             END IF
          END DO
          DO jj = 1, ma_n_nbor(ii)
             nj = ma_nbors(jj,ii)
             vec = vector(mbs(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < closedist)THEN
                closest = ats(nj)%index
                closedist = testdist
             END IF
          END DO
          
          ! find other near neighbors
          !       rangedist = (/ 2.d0, 4.d0, 6.d0, 9.d0, 16.d0 /)
          rangedist = (/ 3.d0, 6.d0, 9.d0, 12.d0, 16.d0 /)
          rangedist = rangedist*closedist
          closelist = 0
          counter1: DO jj = 1, mm_n_nbor(ii)
             nj = mm_nbors(jj,ii)
             vec = vector(mbs(ii),mbs(nj),cell,pbc)
             testdist = (vec.o.vec)
             DO kk = 1, 5
                IF(testdist < rangedist(kk))THEN
                   closelist(kk) = closelist(kk) + 1
                END IF
             END DO
          END DO counter1
          counter2: DO jj = 1, ma_n_nbor(ii)
             nj = ma_nbors(jj,ii)
             vec = vector(mbs(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             DO kk = 1, 5
                IF(testdist < rangedist(kk))THEN
                   closelist(kk) = closelist(kk) + 1
                END IF
             END DO
          END DO counter2
          
          ! write the list of neighbors
          listnumber = 1
          DO jj = 2, 5
             IF(closelist(jj) < 21) listnumber = jj
          END DO
          
          WRITE(ioindex,'(I8,A)',advance='no') mbs(ii)%index," : "
          DO jj = 1, mm_n_nbor(ii)
             nj = mm_nbors(jj,ii)
             vec = vector(mbs(ii),mbs(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < rangedist(listnumber))THEN
                WRITE(ioindex,'(I8,A,F6.3)',advance='no') mbs(nj)%index," ",sqrt(testdist)
             END IF
          END DO
          DO jj = 1, ma_n_nbor(ii)
             nj = ma_nbors(jj,ii)
             vec = vector(mbs(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < rangedist(listnumber))THEN
                WRITE(ioindex,'(I8,A,F6.3)',advance='no') ats(nj)%index," ",sqrt(testdist)
             END IF
          END DO
          WRITE(ioindex,*) ""
          
       END DO
       
       DO ii = 1, n_ats
          closest = 0
          closedist = (cell(1)+cell(2)+cell(3))**2
          DO jj = 1, aa_n_nbor(ii)
             nj = aa_nbors(jj,ii)
             vec = vector(ats(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < closedist)THEN
                closest = ats(nj)%index
                closedist = testdist
             END IF
          END DO
          DO nj = 1, n_mbs
             vec = vector(mbs(nj),ats(ii),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < closedist)THEN
                closest = mbs(nj)%index
                closedist = testdist
             END IF
          END DO
          
          ! find other near neighbors
          !rangedist = (/ 2.d0, 5.d0, 8.d0, 10.d0, 16.d0 /)
          rangedist = (/ 3.d0, 6.d0, 9.d0, 12.d0, 16.d0 /)
          rangedist = rangedist*closedist
          closelist = 0
          counter3: DO jj = 1, aa_n_nbor(ii)
             nj = aa_nbors(jj,ii)
             vec = vector(ats(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             DO kk = 1, 5
                IF(testdist < rangedist(kk))THEN
                   closelist(kk) = closelist(kk) + 1
                END IF
             END DO
          END DO counter3
          counter4: DO nj = 1, n_mbs          
             vec = vector(mbs(nj),ats(ii),cell,pbc)
             testdist = (vec.o.vec)
             DO kk = 1, 5
                IF(testdist < rangedist(kk))THEN
                   closelist(kk) = closelist(kk) + 1
                END IF
             END DO
          END DO counter4
          
          ! write the list of neighbors
          listnumber = 1
          DO jj = 2, 5
             IF(closelist(jj) < 21) listnumber = jj
          END DO
          
          ! write the list of neighbors
          WRITE(ioindex,'(I8,A)',advance='no') ats(ii)%index," : "
          DO jj = 1, aa_n_nbor(ii)
             nj = aa_nbors(jj,ii)
             vec = vector(ats(ii),ats(nj),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < rangedist(listnumber))THEN
                WRITE(ioindex,'(I8,A,F6.3)',advance='no') ats(nj)%index," ",sqrt(testdist)
             END IF
          END DO
          DO nj = 1, n_mbs
             vec = vector(mbs(nj),ats(ii),cell,pbc)
             testdist = (vec.o.vec)
             IF(testdist < rangedist(listnumber))THEN
                WRITE(ioindex,'(I8,A,F6.3)',advance='no') mbs(nj)%index," ",sqrt(testdist)
             END IF
          END DO
          WRITE(ioindex,*) ""
          
       END DO
    END IF


    IF(control%verblevel >= 3)THEN ! write the forces acting on the particles
       WRITE(ioindex,'(A)') "\n"
       WRITE(ioindex,'(A)') HORILINE//"\n"
       WRITE(ioindex,'(A)') "Forces acting on particles: \n"
       DO ii = 1, n_mbs
          WRITE(ioindex,'(I8,F18.6,F18.6,F18.6,A)') mbs(ii)%index, m_force(1:3,ii)," ! "//mbname
       END DO
       DO ii = 1, n_ats
          WRITE(ioindex,'(I8,F18.6,F18.6,F18.6,A)') ats(ii)%index, a_force(1:3,ii)," ! "//params%atomic_labels(ats(ii)%type)
       END DO
       WRITE(ioindex,'(A)') "\nTorques acting on MB molecules: \n"
       DO ii = 1, n_mbs
          WRITE(ioindex,'(I8,F18.6,F18.6,F18.6,A)') mbs(ii)%index, m_torque(1:3,ii)," ! "//mbname
       END DO
    END IF

    ! basic info
    IF(control%verblevel >= 0)THEN
       WRITE(ioindex,'(A)') "\n"
       WRITE(ioindex,'(A)') HORILINE//"\n"
       
       WRITE(ioindex,'(A,F20.4)') "Simulation volume:          ", cell(1)*cell(2)*cell(3) ! volume
       WRITE(ioindex,'(A)') ""
       
       ! number of particles
       WRITE(ioindex,'(A,I20)')   "Number of particles:        ", n_mbs+n_ats 
       WRITE(ioindex,'(A,I20)')   "  water molecules:          ", n_mbs 
       WRITE(ioindex,'(A,I20)')   "  atomic particles:         ", n_ats 
       IF(n_elems > 0)THEN
          DO ii = 1, n_elems
             WRITE(ioindex,'(A,I20)') "     "//params%atomic_labels(ii)//" :                ", found_types(ii)           
          END DO
       END IF
       WRITE(ioindex,'(A)') ""
       
       WRITE(ioindex,'(A,F20.4)') "Initial temperature:        ", temper ! temperature
       !WRITE(ioindex,'(A)') ""
       WRITE(ioindex,'(A,F20.4)') "Initial energy:             ", e_kin+e_pot ! energy
       WRITE(ioindex,'(A,F20.4)') "  kinetic:                  ", e_kin ! kinetic energy
       WRITE(ioindex,'(A,F20.4)') "  potential:                ", e_pot ! potential energy
       WRITE(ioindex,'(A,F20.4)') "     MB-MB Lennard-Jones:   ", e_lenjon
       WRITE(ioindex,'(A,F20.4)') "     MB-MB hydrogen bond:   ", e_bonds
       IF(n_ats > 0)THEN
          WRITE(ioindex,'(A,F20.4)') "     MB-atom potential:     ", e_mbat
          WRITE(ioindex,'(A,F20.4)') "     atom-atom potential:   ", e_atat
       END IF
       IF(is_constr)THEN
          WRITE(ioindex,'(A,F20.4)') "     constraint potential:  ", e_constr
       END IF
    
       WRITE(ioindex,'(A)') "\n"
    END IF
 END IF

    !CALL flush(ioindex)

    RETURN
  END SUBROUTINE write_header


  ! Formats the current simulation parameters and system state to
  ! a valid input file and writes it to a file.
  ! *ioindex index for the output file (must be already open)
  ! *control control parameters
  ! *params physical parameters
  ! *cell supercell dimensions
  ! *pbc true for periodic boundaries
  ! *btype boundary type index
  ! *bval boundary parameter
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_elems number of different types of atoms
  ! *found_types the number of atoms of each type
  ! *is_constr true if there are constraints in place  
  SUBROUTINE write_input(ioindex,mbs,ats,n_elems,&
       params,control,cell,btype,bval,is_constr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ioindex
    !CHARACTER(LEN=*) :: system
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3)
    LOGICAL, INTENT(IN) :: is_constr
    INTEGER, INTENT(IN) :: n_elems, btype(3)    
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    CHARACTER(LEN=20) :: writetag
    CHARACTER(LEN=500) :: writeline
    INTEGER :: ii, jj, kk, ll, n_mbs, n_ats, n_pot, contline, group_ind
    LOGICAL :: is_constrained

    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)


    IF(cpu_id == master_cpu)THEN
    ! control
    WRITE(ioindex,'(A)') "<"//main_block//">"
    IF(control%md_algo == cg_index)THEN ! cg
       WRITE(ioindex,'(A15,I25,A)') CGstepstag//" ",FLOOR(control%time_max+0.5d0)," ! max number of CG steps"
       WRITE(ioindex,'(A15,ES25.5,A)') CGtoltag//" ",control%time_step," ! CG convergence tolerance"
    ELSE ! md
       IF(control%unitscale == true_units)THEN ! real units
          WRITE(ioindex,'(A15,F25.5,A)') Temptag//" ",control%ini_temp," ! [K]"
          WRITE(ioindex,'(A15,ES25.5,A)') Pressuretag//" ",control%ini_press," ! [eV/Ang**3]"
          WRITE(ioindex,'(A15,F25.5,A)') Tsimutag//" ",control%time_max," ! simulation time [fs]"
          WRITE(ioindex,'(A15,F25.5,A)') Tsteptag//" ",control%time_step," ! MD timestep [fs]"
       ELSE ! reduced units
          WRITE(ioindex,'(A15,F25.5,A)') Temptag//" ",control%ini_temp," "
          WRITE(ioindex,'(A15,ES25.5,A)') Pressuretag//" ",control%ini_press," "
          WRITE(ioindex,'(A15,F25.5,A)') Tsimutag//" ",control%time_max," ! simulation time"
          WRITE(ioindex,'(A15,F25.5,A)') Tsteptag//" ",control%time_step," ! MD timestep"
       END IF
    END IF
    SELECT CASE(control%md_algo) ! algorithm
    CASE(velocity_verlet_index)
       WRITE(ioindex,'(A15,A25)') algotag//" ",verlettag
    CASE(leapfrog_index)
       WRITE(ioindex,'(A15,A25)') algotag//" ",leapfrogtag
    CASE(cg_index)
       WRITE(ioindex,'(A15,A25)') algotag//" ",cgtag
    END SELECT
    SELECT CASE(control%md_thermo) ! thermostat
    CASE(microcanonical_index)
       WRITE(ioindex,'(A15,A25)') thermotag//" ",nonetag
    CASE(langevin_index)
       WRITE(ioindex,'(A15,A13,F12.5)') thermotag//" ",langevintag//" ",control%thermo_value
    CASE(cooler_index)
       WRITE(ioindex,'(A15,A13,F12.5)') thermotag//" ",coolertag//" ",control%thermo_value
    END SELECT
    SELECT CASE(control%md_baro) ! barostat
    CASE(microcanonical_index)
       WRITE(ioindex,'(A15,A25)') barotag//" ",nonetag
    CASE(berendsen_index)
       WRITE(ioindex,'(A15,A13,F12.5,F13.5)') barotag//" ",berendsentag//" ",control%baro_value,&
            control%baro_interval
    END SELECT
    SELECT CASE(control%unitscale) ! units
    CASE(true_units)
       WRITE(ioindex,'(A15,A25)') unitstag//" ",truetag       
    CASE(reduced_units)
       WRITE(ioindex,'(A15,A25)') unitstag//" ",reducedtag
    END SELECT
    SELECT CASE(control%runtype) ! run type
    CASE(test_index)
       WRITE(ioindex,'(A15,A25)') runtag//" ",testtag
    CASE(full_index)
       WRITE(ioindex,'(A15,A25)') runtag//" ",fulltag
    END SELECT
    SELECT CASE(control%xyz_writer) ! xyz file writer
    CASE(noxyz_index)
       WRITE(ioindex,'(A15,A25)') xyztag//" ", nonetag
    CASE(sxyz_index)
       WRITE(ioindex,'(A15,A25)') xyztag//" ", starttag
    CASE(exyz_index)
       WRITE(ioindex,'(A15,A25)') xyztag//" ", endtag
    CASE(sexyz_index)
       WRITE(ioindex,'(A15,A25)') xyztag//" ", sendtag
    CASE(ixyz_index)
       WRITE(ioindex,'(A15,A13,F12.5)') xyztag//" ", intervaltag//" ", control%xyz_interval
    END SELECT
    SELECT CASE(control%inp_writer) ! continuation file writer
    CASE(noxyz_index)
       WRITE(ioindex,'(A15,A25)') conttag//" ", nonetag
    CASE(sxyz_index)
       WRITE(ioindex,'(A15,A25)') conttag//" ", starttag
    CASE(exyz_index)
       WRITE(ioindex,'(A15,A25)') conttag//" ", endtag
    CASE(sexyz_index)
       WRITE(ioindex,'(A15,A25)') conttag//" ", sendtag
    CASE(ixyz_index)
       WRITE(ioindex,'(A15,A13,F12.5)') conttag//" ", intervaltag//" ", control%inp_interval
    END SELECT
    WRITE(ioindex,'(A15,I25)') seedtag//" ",control%rnd_seed ! rng seed
    IF(control%verbose > control%time_max)THEN
       WRITE(ioindex,'(A15,A25,I5)') verbosetag//" ",nonetag,control%verblevel
    ELSE
       WRITE(ioindex,'(A15,F25.5,I5)') verbosetag//" ",control%verbose,control%verblevel
    END IF
    WRITE(ioindex,'(A)') "</"//main_block//">"
    ! mb-model
    WRITE(ioindex,'(A)') "<"//mb_block//">"
    IF(control%unitscale == true_units)THEN ! real units
       WRITE(ioindex,'(A10,F20.5,A)') E_LJtag//" ",params%e_lj," ! Lennard-Jones energy coeff. [eV]"
       WRITE(ioindex,'(A10,F20.5,A)') E_HBtag//" ",params%e_hb," ! hydrogen-bond energy coeff. [eV]"
       WRITE(ioindex,'(A10,F20.5,A)') E_PHtag//" ",params%e_phi," ! dihedral angle energy coeff. [1]"
       WRITE(ioindex,'(A10,F20.5,A)') R_HBtag//" ",params%r_hb," ! H-bond eq. length [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_LJtag//" ",params%s_lj," ! Lennard-Jones sigma [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_HBtag//" ",params%s_r," ! H-bond length sigma [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_THtag//" ",params%s_th," ! H-bond angle sigma [1]"
       WRITE(ioindex,'(A10,F20.5,A)') Exp_btag//" ",params%v_b," ! Bond order exponent [1]"
       WRITE(ioindex,'(A10,F20.5,A)') R_btag//" ",params%r_b," ! Bond range [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') D_btag//" ",params%d_b," ! Bond cut width [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') cut_ljtag//" ",params%cut_lj," ! Lennard-Jones cutoff [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') cut_hbtag//" ",params%cut_hb," ! H-bond cutoff [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') cut_vtag//" ",params%cut_ver," ! Neighbor list cutoff [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') cut_atag//" ",params%cut_atom," ! Atomic potential cutoff [Ang]"
       WRITE(ioindex,'(A10,F20.5,A)') m_moletag//" ",params%m_mol/m_scale," ! MB mass [g/mol]"
       WRITE(ioindex,'(A10,F20.5,A)') i_moletag//" ",params%i_mol/m_scale," ! MB mom. of inertia [g/mol*Ang**2]"
       WRITE(ioindex,'(A10,F20.5,A)') l_ohtag//" ",params%l_oh," ! O-H length [Ang]"
    ELSE ! reduced units
       WRITE(ioindex,'(A10,F20.5,A)') E_LJtag//" ",params%e_lj," ! Lennard-Jones energy coeff."
       WRITE(ioindex,'(A10,F20.5,A)') E_HBtag//" ",params%e_hb," ! hydrogen-bond energy coeff."
       WRITE(ioindex,'(A10,F20.5,A)') E_PHtag//" ",params%e_phi," ! dihedral angle energy coeff."
       WRITE(ioindex,'(A10,F20.5,A)') R_HBtag//" ",params%r_hb," ! H-bond eq. length"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_LJtag//" ",params%s_lj," ! Lennard-Jones sigma"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_HBtag//" ",params%s_r," ! H-bond length sigma"
       WRITE(ioindex,'(A10,F20.5,A)') Sigma_THtag//" ",params%s_th," ! H-bond angle sigma"
       WRITE(ioindex,'(A10,F20.5,A)') Exp_btag//" ",params%v_b," ! Bond order exponent"
       WRITE(ioindex,'(A10,F20.5,A)') R_btag//" ",params%r_b," ! Bond range"
       WRITE(ioindex,'(A10,F20.5,A)') D_btag//" ",params%d_b," ! Bond cut width"
       WRITE(ioindex,'(A10,F20.5,A)') cut_ljtag//" ",params%cut_lj," ! Lennard-Jones cutoff"
       WRITE(ioindex,'(A10,F20.5,A)') cut_hbtag//" ",params%cut_hb," ! H-bond cutoff"
       WRITE(ioindex,'(A10,F20.5,A)') cut_vtag//" ",params%cut_ver," ! Neighbor list cutoff"
       WRITE(ioindex,'(A10,F20.5,A)') cut_atag//" ",params%cut_atom," ! Atomic potential cutoff"
       WRITE(ioindex,'(A10,F20.5,A)') m_moletag//" ",params%m_mol/m_scale," ! MB mass"
       WRITE(ioindex,'(A10,F20.5,A)') i_moletag//" ",params%i_mol/m_scale," ! MB mom. of inertia"
       WRITE(ioindex,'(A10,F20.5,A)') l_ohtag//" ",params%l_oh," ! O-H length"
    END IF
    WRITE(ioindex,'(A)') "</"//mb_block//">"
    ! cell
    WRITE(ioindex,'(A)') "<"//cell_block//">"
    DO ii = 1,3 ! loop x, y, z
       IF(btype(ii) == periodic_bound_index)THEN ! periodic boundary
          WRITE(ioindex,'(F20.7,A20)') cell(ii), "  "//periodictag
       ELSE IF(btype(ii) == free_bound_index)THEN ! free boundary
          WRITE(ioindex,'(F20.7,A20)') cell(ii), "  "//freetag
       ELSE IF(btype(ii) == wall_bound_index)THEN ! harmonic wall
          WRITE(ioindex,'(F20.7,A8,F12.5)') cell(ii), "  "//walltag//" ",bval(ii)
       END IF
    END DO
    WRITE(ioindex,'(A)') "</"//cell_block//">"
    ! elements
    IF(n_elems > 0)THEN
       WRITE(ioindex,'(A)') "<"//ele_block//">"
       DO ii = 1, n_elems ! loop over elements
          WRITE(ioindex,'(A7,F10.5)') params%atomic_labels(ii), params%m_atoms(ii)/m_scale
       END DO
    WRITE(ioindex,'(A)') "</"//ele_block//">"
    END IF
    ! particles
    WRITE(ioindex,'(A)') "<"//pos_block//">"
    DO ii = 1, n_mbs ! loop over mbs
       WRITE(ioindex,'(I8,A7,F18.7,F18.7,F18.7,A5,F12.7,F12.7,F12.7,F12.7)') &
            mbs(ii)%index, " "//mbname//" ", mbs(ii)%pos, "     ", &
            mbs(ii)%orientation%w,mbs(ii)%orientation%x,mbs(ii)%orientation%y,mbs(ii)%orientation%z
    END DO
    DO ii = 1, n_ats ! loop over atoms
       WRITE(ioindex,'(I8,A7,F18.7,F18.7,F18.7)') &
            ats(ii)%index, " "//ats(ii)%element//" ", ats(ii)%pos
    END DO
    WRITE(ioindex,'(A)') "</"//pos_block//">"
    ! velocities
    WRITE(ioindex,'(A)') "<"//vel_block//">"
    DO ii = 1, n_mbs ! loop over mbs
       WRITE(ioindex,'(I8,F18.7,F18.7,F18.7,F18.7,F18.7,F18.7,A)') &
            mbs(ii)%index, mbs(ii)%vel, mbs(ii)%angvel, " ! "//mbname
    END DO
    DO ii = 1, n_ats ! loop over atoms
       WRITE(ioindex,'(I8,F18.7,F18.7,F18.7,A)') &
            ats(ii)%index, ats(ii)%vel, " ! "//ats(ii)%element
    END DO
    WRITE(ioindex,'(A)') "</"//vel_block//">"
    ! constraints

    IF(is_constr)THEN ! write the block only if there are constraints
       WRITE(ioindex,'(A)') "<"//constr_block//">"

       DO ii = 1, n_mbs ! loop over mbs
          is_constrained = .false.
          ! check if this particle is constrained
          DO jj = 1, 3 ! loop over x,y,z
             IF(mbs(ii)%constrained(jj) /= no_constr_index)THEN
                is_constrained = .true.
             END IF
          END DO
          IF(is_constrained)THEN ! write the type of constraints if the particle is not free
             WRITE(ioindex,'(I7,A)',advance='no') mbs(ii)%index,"  "
             DO jj = 1, 3
                SELECT CASE(mbs(ii)%constrained(jj))

                CASE(no_constr_index)
                   writetag = freetag
                CASE(all_frozen_index)
                   writetag = frozentag//" "//alltag
                CASE(frozen_pos_index)
                   writetag = frozentag//" "//postag
                CASE(harmonic_well_index)
                   WRITE(writetag,'(A,F10.5)') welltag//" ", mbs(ii)%well(jj)
                CASE(frozen_vel_index)
                   writetag = frozentag//" "//veltag
                CASE(ext_force_index)
                   WRITE(writetag,'(A,F10.5)') forcetag//" ", mbs(ii)%well(jj)

                END SELECT
                
                WRITE(ioindex,'(A17)',advance='no') writetag//"  "
             END DO
             WRITE(ioindex,'(A)') " ! "//mbname
          END IF
       END DO

       DO ii = 1, n_ats ! loop over atoms
          is_constrained = .false.
          ! check if this atom is constrained
          DO jj = 1, 3 ! loop over x, y, z
             IF(ats(ii)%constrained(jj) /= no_constr_index)THEN
                is_constrained = .true.
             END IF
          END DO
          ! write the constraints only if the atom is not free
          IF(is_constrained)THEN
             WRITE(ioindex,'(I7,A)',advance='no') ats(ii)%index,"  "
             DO jj = 1, 3
                SELECT CASE(ats(ii)%constrained(jj))

                CASE(no_constr_index)
                   writetag = freetag
                CASE(all_frozen_index)
                   writetag = frozentag//" "//alltag
                CASE(frozen_pos_index)
                   writetag = frozentag//" "//postag
                CASE(harmonic_well_index)
                   WRITE(writetag,'(A,F10.5)') welltag//" ", ats(ii)%well(jj)
                CASE(frozen_vel_index)
                   writetag = frozentag//" "//veltag
                CASE(ext_force_index)
                   WRITE(writetag,'(A,F10.5)') forcetag//" ", ats(ii)%well(jj)

                END SELECT
                
                WRITE(ioindex,'(A17)',advance='no') writetag//" "
             END DO
             WRITE(ioindex,'(A)') " ! "//ats(ii)%element
          END IF
       END DO

       WRITE(ioindex,'(A)') "</"//constr_block//">"
    END IF

    ! potentials
    IF(n_elems > 0)THEN ! write if there are atoms
       WRITE(ioindex,'(A)') "<"//pot_block//">"
       DO ii = -1, n_elems ! loop over mb, hbonds, and atom types
          DO jj = ii, n_elems ! loop over other types
             IF(params%n_pots(ii,jj) > 0)THEN

                ! the name of the first type
                IF(ii == -1)THEN
                   writetag = " "//armname
                ELSE IF(ii == 0)THEN
                   writetag = " "//mbname
                ELSE
                   writetag = " "//params%atomic_labels(ii)
                END IF

                ! this shouldn't happen
                IF(jj < 1) CALL abort("MB - MB potential not allowed, how did this happen?")

                ! write the line defining the types for which the potential works
                WRITE(ioindex,'(A7,A7,F20.3)') writetag,params%atomic_labels(jj),params%pot_cut(ii,jj)
                
                DO kk = 1, params%n_pots(ii,jj) ! loop over number of potential terms

                   SELECT CASE(params%pot_types(kk,ii,jj)) ! type
                   CASE(lenjon_index) ! lennard-jones
                      n_pot = 2
                      WRITE(ioindex,'(A10)',advance='no') "  "//lenjon_pot//" "
                   CASE(exp_index) ! exponential
                      n_pot = 2
                      WRITE(ioindex,'(A10)',advance='no') "  "//exp_pot//" "
                   CASE(pow_index) ! power law
                      n_pot = 2
                      WRITE(ioindex,'(A10)',advance='no') "  "//pow_pot//" "
                   CASE(spring_index) ! spring
                      n_pot = 2
                      WRITE(ioindex,'(A10)',advance='no') "  "//spring_pot//" "
                   CASE(hard_index) ! hard sphere
                      n_pot = 3
                      WRITE(ioindex,'(A10)',advance='no') "  "//hard_pot//" "
                   CASE(hardrep_index) ! repulsive hard sphere
                      n_pot = 3
                      WRITE(ioindex,'(A10)',advance='no') "  "//hardrep_pot//" "
                   CASE(fene_index) ! fene
                      n_pot = 3
                      WRITE(ioindex,'(A10)',advance='no') "  "//fene_pot//" "
                   END SELECT

                   DO ll = 1, n_pot ! write the potential parameters
                      WRITE(ioindex,'(F10.5,A)',advance='no') params%pot_params(ll,kk,ii,jj)," "
                   END DO
                   WRITE(ioindex,'(A)') ""
                END DO

             END IF
          END DO
       END DO
       WRITE(ioindex,'(A)') "</"//pot_block//">"
    END IF

    ! statistics
    IF(control%n_statfiles > 0 .OR. & ! only write the block if some stats are to be recorded
         control%rdf_writer /= noxyz_index .OR. &
         control%adf_writer /= noxyz_index .OR. &
         control%bond_writer /= noxyz_index )THEN

       WRITE(ioindex,'(A)') "<"//stat_block//">"

       ! bonds
       SELECT CASE(control%bond_writer)
       CASE(noxyz_index)
          WRITE(ioindex,'(A15,A25)') bondtag//" ", nonetag
       CASE(sxyz_index)
          WRITE(ioindex,'(A15,A25)') bondtag//" ", starttag
       CASE(exyz_index)
          WRITE(ioindex,'(A15,A25)') bondtag//" ", endtag
       CASE(sexyz_index)
          WRITE(ioindex,'(A15,A25)') bondtag//" ", sendtag
       CASE(ixyz_index)
          WRITE(ioindex,'(A15,A13,F12.5)') bondtag//" ", intervaltag//" ", control%bond_interval
       END SELECT

       ! rdf
       SELECT CASE(control%rdf_writer)
       CASE(noxyz_index)
          WRITE(ioindex,'(A15,A25)') rdftag//" ", nonetag
       CASE(sxyz_index)
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", starttag, alltag//" ", control%rdf_range
          CASE(mb_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", starttag, mbname//" ", control%rdf_range
          CASE(atom_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", starttag, atomtag//" ", control%rdf_range
          END SELECT
       CASE(exyz_index)
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", endtag, alltag//" ", control%rdf_range
          CASE(mb_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", endtag, mbname//" ", control%rdf_range
          CASE(atom_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", endtag, atomtag//" ", control%rdf_range
          END SELECT
       CASE(sexyz_index)
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", sendtag, alltag//" ", control%rdf_range
          CASE(mb_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", sendtag, mbname//" ", control%rdf_range
          CASE(atom_rdf_index)
             WRITE(ioindex,'(A15,A25,A15,F15.5)') rdftag//" ", sendtag, atomtag//" ", control%rdf_range
          END SELECT
       CASE(ixyz_index)
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,A15,F15.5)') &
                  rdftag//" ", intervaltag//" ", control%rdf_interval, alltag//" ", control%rdf_range
          CASE(mb_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,A15,F15.5)') &
                  rdftag//" ", intervaltag//" ", control%rdf_interval, mbname//" ", control%rdf_range
          CASE(atom_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,A15,F15.5)') &
                  rdftag//" ", intervaltag//" ", control%rdf_interval, atomtag//" ", control%rdf_range
          END SELECT
       CASE(average_index)
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,F12.5,A15,F15.5)') &
                  rdftag//" ", averagetag//" ", control%rdf_interval, control%rdf_start, alltag//" ", control%rdf_range
          CASE(mb_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,F12.5,A15,F15.5)') &
                  rdftag//" ", averagetag//" ", control%rdf_interval, control%rdf_start, mbname//" ", control%rdf_range
          CASE(atom_rdf_index)
             WRITE(ioindex,'(A15,A13,F12.5,F12.5,A15,F15.5)') &
                  rdftag//" ", averagetag//" ", control%rdf_interval, control%rdf_start, atomtag//" ", control%rdf_range
          END SELECT
       END SELECT

       ! adf
       writeline = " "
       SELECT CASE(control%adf_writer)
       CASE(noxyz_index)
          WRITE(writeline(1:40),'(A15,A25)') adftag//" ", nonetag
          contline = 0
       CASE(sxyz_index)
          WRITE(writeline(1:40),'(A15,A25)') adftag//" ", starttag  
          contline = 42    
       CASE(exyz_index)
          WRITE(writeline(1:40),'(A15,A25)') adftag//" ", endtag   
          contline = 42
       CASE(sexyz_index)
          WRITE(writeline(1:40),'(A15,A25)') adftag//" ", sendtag
          contline = 42
       CASE(ixyz_index)
          WRITE(writeline(1:40),'(A15,A13,F12.5)') &
                  adftag//" ", intervaltag//" ", control%adf_interval
          contline = 42
       CASE(average_index)
          WRITE(writeline(1:52),'(A15,A13,F12.5,F12.5)') &
               adftag//" ", averagetag//" ", control%adf_interval, control%adf_start
          contline = 54
       END SELECT
       IF(contline > 0)THEN ! only write the details if not "none"
          SELECT CASE(control%adf_axis)
          CASE(x_index)
             writeline(contline:contline+12) = "    "//xaxistag//"   "
             contline = contline+12
          CASE(y_index)
             writeline(contline:contline+12) = "    "//yaxistag//"   "
             contline = contline+12
          CASE(z_index)
             writeline(contline:contline+12) = "    "//zaxistag//"   "
             contline = contline+12
          END SELECT
          SELECT CASE(control%adf_angle)
          CASE(.true.)
             writeline(contline:contline+9) = "   "//angletag//" "
             contline = contline+9             
          CASE(.false.)
             writeline(contline:contline+9) = " "//noangletag//" "
             contline = contline+9    
          END SELECT
          IF(control%adf_follow /= 0)THEN
             WRITE(writeline(contline:contline+20),'(A10,I10)') followtag//" ", mbs(control%adf_follow)%index
             contline = contline + 20
          ELSE
             WRITE(writeline(contline:contline+55),'(A10,F15.5,F15.5,F15.5)') fixedtag//" ", control%adf_center(1:3)
             contline = contline + 55
          END IF
          SELECT CASE(control%adf_particles)
          CASE(all_rdf_index)
             WRITE(writeline(contline:contline+30),'(A15,F15.5)') alltag//" ", control%adf_range
             contline = contline+30
          CASE(mb_rdf_index)
             WRITE(writeline(contline:contline+30),'(A15,F15.5)') mbname//" ", control%adf_range
             contline = contline+30
          CASE(atom_rdf_index)
             WRITE(writeline(contline:contline+30),'(A15,F15.5)') atomtag//" ", control%adf_range
             contline = contline+30
          END SELECT
       ELSE ! for none, don't write the rest
          contline = 40
       END IF
       WRITE(ioindex,'(A)') writeline(1:contline)

       IF(control%n_statfiles > 0)THEN ! write statfile definitions

          DO group_ind = 1, control%n_statfiles ! loop over statfiles

             WRITE(ioindex,'(A15)') writestattag
             WRITE(ioindex,'(A15,F25.8)') starttag, control%stat_start(group_ind)
             WRITE(ioindex,'(A15,F25.8)') intervaltag, control%stat_interval(group_ind)

             DO ii = 1, control%n_stats(group_ind) ! loop over stats
                SELECT CASE(control%stats(ii,group_ind))
                CASE(time_stat_index)
                   WRITE(ioindex,'(A15)') timetag
                CASE(clock_stat_index)
                   WRITE(ioindex,'(A15)') clocktag
                CASE(temp_stat_index)
                   WRITE(ioindex,'(A15)') temptag
                CASE(kin_stat_index)
                   WRITE(ioindex,'(A15)') kinetictag
                CASE(pot_stat_index)
                   WRITE(ioindex,'(A15)') potentialtag
                CASE(ene_stat_index)
                   WRITE(ioindex,'(A15)') energytag
                CASE(press_stat_index)
                   WRITE(ioindex,'(A15)') pressuretag
                CASE(linear_stat_index)
                   WRITE(ioindex,'(A15)') lineartag
                CASE(rotational_stat_index)
                   WRITE(ioindex,'(A15)') rottag
                CASE(ljpot_stat_index)
                   WRITE(ioindex,'(A15)') lenjontag
                CASE(hbpot_stat_index)
                   WRITE(ioindex,'(A15)') hbondtag
                CASE(mapot_stat_index)
                   WRITE(ioindex,'(A15)') mbatomtag
                CASE(aapot_stat_index)
                   WRITE(ioindex,'(A15)') atomatomtag
                CASE(volume_stat_index)
                   WRITE(ioindex,'(A15)') volumetag
                CASE(virial_stat_index)
                   WRITE(ioindex,'(A15)') virialtag
                CASE(coord_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') coordtag, mbs(control%stat_particles(ii,group_ind))%index
                   ELSE
                      WRITE(ioindex,'(A15,I25)') coordtag, ats(-control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(vel_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') velotag, mbs(control%stat_particles(ii,group_ind))%index
                   ELSE
                      WRITE(ioindex,'(A15,I25)') velotag, ats(-control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(ang_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') angulartag, mbs(control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(orient_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') orientag, mbs(control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(force_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') forcetag, mbs(control%stat_particles(ii,group_ind))%index
                   ELSE
                      WRITE(ioindex,'(A15,I25)') forcetag, ats(-control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(torque_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') torquetag, mbs(control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(arms_stat_index)
                   IF(control%stat_particles(ii,group_ind) > 0)THEN
                      WRITE(ioindex,'(A15,I25)') armstag, mbs(control%stat_particles(ii,group_ind))%index
                   END IF
                CASE(forcesum_stat_index)
                   WRITE(ioindex,'(A15)',advance='no') forcesumtag
                   listing: DO jj = 1, SIZE(control%stat_groups(:,ii,group_ind))
                      IF(control%stat_groups(jj,ii,group_ind) == 0)THEN
                         EXIT listing
                      END IF
                      IF(control%stat_groups(jj,ii,group_ind) > 0)THEN
                         WRITE(ioindex,'(A,I10)',advance='no') " ",mbs(control%stat_groups(jj,ii,group_ind))%index
                      ELSE
                         WRITE(ioindex,'(A,I10)',advance='no') " ",ats(-control%stat_groups(jj,ii,group_ind))%index
                      END IF
                   END DO listing
                   WRITE(ioindex,'(A)',advance='yes') " "             
                   
                END SELECT
             END DO
          END DO
       END IF

       WRITE(ioindex,'(A)') "</"//stat_block//">"
    END IF
    END IF

    RETURN
  END SUBROUTINE write_input

  ! Writes the current system geometry to file
  ! in xyz format. MB molecules are written as one
  ! O atom in the center and 4 H atoms representing
  ! the H-bond arms. Of course, only two of these arms
  ! are really hydrogen atoms - the other two are just
  ! dangling H-bonds - but they are treated equivalently
  ! in the model. The O-H distance must be given.
  ! *filename name of the xyz file to be written - the string must be of the right length
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *hands the length of the h-bond arms
  ! *newfile if true, a new file is written. otherwise, an existing file is appended
  SUBROUTINE write_xyz(filename,mbs,ats,hands,newfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(IN) :: hands
    INTEGER :: n_mbs, n_ats
    INTEGER :: ii, ioindex, iostat
    REAL(KIND=dp) :: h1(3), h2(3), h3(3), h4(3)
    LOGICAL, INTENT(IN) :: newfile


    IF(cpu_id == master_cpu)THEN ! master writes

       ! open output
       ioindex = 261
       IF(newfile)THEN ! start a new file
          OPEN(ioindex,FILE=filename,IOSTAT=iostat)
       ELSE ! continue in an existing file
          OPEN(ioindex,FILE=filename,POSITION="APPEND",IOSTAT=iostat)
       END IF
       IF(iostat /= 0) CALL abort("writing xyz")
       
       n_mbs = mbs_size(mbs)
       n_ats = ats_size(ats)
       WRITE(ioindex,'(A)') ""
       WRITE(ioindex,'(I10)') n_mbs*5+n_ats

       DO ii = 1, n_mbs ! loop over mbs
          CALL get_HB_vectors(mbs(ii),h1,h2,h3,h4)
          ! write the center as O and the tips of the arms as H
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               "O ", mbs(ii)%pos(1), mbs(ii)%pos(2), mbs(ii)%pos(3)  
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               "  H ", mbs(ii)%pos(1)+hands*h1(1), mbs(ii)%pos(2)+hands*h1(2), mbs(ii)%pos(3)+hands*h1(3)  
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               "  H ", mbs(ii)%pos(1)+hands*h2(1), mbs(ii)%pos(2)+hands*h2(2), mbs(ii)%pos(3)+hands*h2(3)  
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               "  H ", mbs(ii)%pos(1)+hands*h3(1), mbs(ii)%pos(2)+hands*h3(2), mbs(ii)%pos(3)+hands*h3(3)  
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               "  H ", mbs(ii)%pos(1)+hands*h4(1), mbs(ii)%pos(2)+hands*h4(2), mbs(ii)%pos(3)+hands*h4(3)       
       END DO
       DO ii = 1, n_ats ! loop over atoms
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') &
               ats(ii)%element, ats(ii)%pos(1), ats(ii)%pos(2), ats(ii)%pos(3)       
       END DO
       
       ! close output
       CLOSE(ioindex)
    END IF

  END SUBROUTINE write_xyz

  ! Writes the current system geometry to file
  ! in xyz format and the bond numbers in another file. Only MB centers are written.
  ! *filename1 name of the xyz file to be written
  ! *filename2 name of the bond count file
  ! *mbs list of molecules
  ! *n_bond an array listing the bond counts for all MBs
  ! *newfile if true, a new file is written. otherwise, an existing file is appended
  SUBROUTINE write_bonds(filename1,filename2,n_bond,mbs,mm_nbors,mm_n_ns,cell,pbc,control,params,newfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename1, filename2
    TYPE(mb), POINTER :: mbs(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_ns(:)
    REAL(KIND=dp), POINTER :: n_bond(:)
    INTEGER :: n_mbs
    INTEGER :: ii, ioindex, iostat
    LOGICAL, INTENT(IN) :: newfile

    IF(cpu_id == master_cpu)THEN
       n_mbs = mbs_size(mbs)
       IF(n_mbs == 0) RETURN ! if there are no mbs, don't do anything
       
       ! write xyz file of just the mb centers
       ioindex = 261
       IF(newfile)THEN ! new file
          OPEN(ioindex,FILE=filename1,IOSTAT=iostat)
       ELSE ! append an existing file
          OPEN(ioindex,FILE=filename1,POSITION="APPEND",IOSTAT=iostat)
       END IF
       IF(iostat /= 0) CALL abort("writing mb xyz")
       
       WRITE(ioindex,'(A)') ""
       WRITE(ioindex,'(I10)') n_mbs
       
       DO ii = 1, n_mbs ! loop over mbs
          WRITE(ioindex,'(A,F20.8,F20.8,F20.8)') & ! write the center coordinates
               "MB ", mbs(ii)%pos(1), mbs(ii)%pos(2), mbs(ii)%pos(3)  
       END DO
       
       ! close output for mb.xyz
       CLOSE(ioindex)
              
       ! write the bond counts to a separate file
       ioindex = 5963
       IF(newfile)THEN ! new file
          OPEN(ioindex,FILE=filename2,IOSTAT=iostat)
       ELSE ! append an existing file
          OPEN(ioindex,FILE=filename2,POSITION="APPEND",IOSTAT=iostat)
       END IF
       IF(iostat /= 0) CALL abort("writing bonds")
       
       DO ii = 1, n_mbs ! loop over mbs
          ! write the bonding data, all on one line
          WRITE(ioindex,'(F8.4)',advance='no') n_bond(ii)
       END DO
       WRITE(ioindex,*) ""
       
       ! close output for bonds.out
       CLOSE(ioindex)
    END IF

  END SUBROUTINE write_bonds

  ! Calculates the radial distribution function and writes it to a given file.
  ! The parameters for RDF generating are stored in the control parameters.
  ! *filename name of the rdf file to be written - the string must be of the right length
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true for periodicc boundaries
  ! *control control parameters
  ! *newfile if true, a new file is written. otherwise, an existing file is appended
  ! *therdf if given, this array will be written as the rdf
  ! *rdf_error if given, this array will be written along therdf and it should be the error estimate for therdf
  SUBROUTINE write_rdf(filename,mbs,ats,cell,pbc,&
       mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,aa_nbors,aa_n_nbor,params,&
       control,newfile,therdf,rdf_error)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(cps), INTENT(IN) :: control
    TYPE(mbps), INTENT(IN) :: params
    INTEGER :: n_mbs, n_ats, ioindex, ii, iostat
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3), newfile
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_nbor(:), &
         ma_nbors(:,:), ma_n_nbor(:), aa_nbors(:,:), aa_n_nbor(:)
    REAL(KIND=dp) :: rdf(rdfbins,4), binwidth
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: therdf(rdfbins,4), rdf_error(rdfbins,4)
    
    
    IF(cpu_id == master_cpu)THEN
       ioindex = 434 ! open output
       IF(newfile)THEN ! new file
          OPEN(ioindex,FILE=filename,IOSTAT=iostat)
       ELSE ! append an existing file
          OPEN(ioindex,FILE=filename,POSITION="APPEND",IOSTAT=iostat)
       END IF
       IF(iostat /= 0) CALL abort("writing rdf")
       
       n_mbs = mbs_size(mbs)
       n_ats = ats_size(ats)

       ! separation between the data points
       binwidth = control%rdf_range/REAL(rdfbins,KIND=dp)
       
       IF(present(therdf))THEN ! rdf is given
          rdf = therdf
          
          IF(present(rdf_error))THEN ! rdf errorbars are given
             
             DO ii = 1, rdfbins ! loop over data points
                SELECT CASE(control%rdf_particles) ! which particles are tracked
                CASE(all_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5)') &
                        binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,1:4), rdf_error(ii,1:4)
                CASE(mb_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                        binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,1), rdf_error(ii,1)
                CASE(atom_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                        binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,2), rdf_error(ii,2)
                END SELECT
             END DO
             
             WRITE(ioindex,'(A)') ""
             
             ! close output
             CLOSE(ioindex)
             
             RETURN ! end here
          
          END IF
       ELSE ! the rdf must be calculated
          CALL radial_distribution_function(mbs(:),ats(:),n_mbs,n_ats,cell,pbc,&
               mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,aa_nbors,aa_n_nbor,params,&
               control%rdf_particles,control%rdf_range,rdf)
       END IF
       
       ! write the rdf
       DO ii = 1, rdfbins ! loop over data points
          SELECT CASE(control%rdf_particles)
          CASE(all_rdf_index)
             WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5,F30.5)') &
                  binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,1:4)
          CASE(mb_rdf_index)
             WRITE(ioindex,'(F30.5,F30.5)') &
                  binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,1)
          CASE(atom_rdf_index)
             WRITE(ioindex,'(F30.5,F30.5)') &
                  binwidth*(REAL(ii,KIND=dp)-0.5d0), rdf(ii,2)
          END SELECT
       END DO
       
       WRITE(ioindex,'(A)') ""
       
       ! close output
       CLOSE(ioindex)
    END IF

    RETURN
  END SUBROUTINE write_rdf

  ! Calculates the axial distribution function and writes it to a given file.
  ! The parameters for ADF generating are stored in the control parameters.
  ! *filename name of the rdf file to be written - the string must be of the right length
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true for periodicc boundaries
  ! *params physical parameters
  ! *control control parameters
  ! *newfile if true, a new file is written. otherwise, an existing file is appended
  ! *theadf if given, this array will be written as the adf
  ! *adf_error if given, this array will be written along theadf and it should be the error estimate for theadf
  SUBROUTINE write_adf(filename,mbs,ats,cell,pbc,params,&
       control,newfile,theadf,adf_error)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(cps), INTENT(IN) :: control
    TYPE(mbps), INTENT(IN) :: params
    INTEGER :: n_mbs, n_ats, ioindex, ii, jj, iostat
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3), newfile
    REAL(KIND=dp) :: adf(rdfbins,rdfbins,3), binwidth(2)
    REAL(KIND=dp), OPTIONAL, INTENT(IN) :: theadf(rdfbins,rdfbins,3), adf_error(rdfbins,rdfbins,3)
    

    IF(cpu_id == master_cpu)THEN
       ! open output
       ioindex = 434
       IF(newfile)THEN
          OPEN(ioindex,FILE=filename,IOSTAT=iostat)
       ELSE
          OPEN(ioindex,FILE=filename,POSITION="APPEND",IOSTAT=iostat)
       END IF
       IF(iostat /= 0) CALL abort("writing adf")
       
       n_mbs = mbs_size(mbs)
       n_ats = ats_size(ats)
       ! data point separation in distance and angle
       binwidth(1) = control%adf_range/REAL(rdfbins,KIND=dp)
       binwidth(2) = twopi/REAL(rdfbins,KIND=dp)
       
       IF(present(theadf))THEN ! adf is given
          adf = theadf
          
          IF(present(adf_error))THEN ! adf error is given
             
             DO ii = 1, rdfbins ! loop over data in distance
                IF(control%adf_angle)THEN ! angle resolved
                   DO jj = 1, rdfbins ! loop over data in angle
                      SELECT CASE(control%adf_particles) ! what particles are measured
                      CASE(all_rdf_index)
                         WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5)') &
                              binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                              adf(jj,ii,1:3), adf_error(jj,ii,1:3)
                      CASE(mb_rdf_index)
                         WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5)') &
                              binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                              adf(jj,ii,1), adf_error(jj,ii,1)
                      CASE(atom_rdf_index)
                         WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5)') &
                              binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                              adf(jj,ii,2), adf_error(jj,ii,2)
                      END SELECT
                   END DO
                   WRITE(ioindex,'(A)') " "
                ELSE ! angle averaged
                   SELECT CASE(control%adf_particles)
                   CASE(all_rdf_index)
                      WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5,F30.5,F30.5,F30.5)') &
                           binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                           adf(1,ii,1:3), adf_error(1,ii,1:3)
                   CASE(mb_rdf_index)
                      WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                           binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                           adf(1,ii,1), adf_error(1,ii,1)
                   CASE(atom_rdf_index)
                      WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                           binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                           adf(1,ii,2), adf_error(1,ii,2)
                   END SELECT
                END IF
             END DO
             
             WRITE(ioindex,'(A)') ""
             
             ! close output
             CLOSE(ioindex)
             
             RETURN
             
          END IF
       ELSE ! the adf must be measured
          CALL axial_distribution_function(mbs(:),ats(:),n_mbs,n_ats,cell,pbc,params,control,adf)
       END IF
       
       DO ii = 1, rdfbins ! loop over the data in distanve
          IF(control%adf_angle)THEN ! angle resolved
             DO jj = 1, rdfbins ! loop over the data in angle
                SELECT CASE(control%adf_particles) ! what particles are measured
                CASE(all_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5,F30.5)') &
                        binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                        adf(jj,ii,1:3)
                CASE(mb_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                        binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                        adf(jj,ii,1)
                CASE(atom_rdf_index)
                   WRITE(ioindex,'(F30.5,F30.5,F30.5)') &
                        binwidth(1)*(REAL(ii,KIND=dp)-0.5d0), binwidth(2)*(REAL(jj,KIND=dp)-0.5d0),&
                        adf(jj,ii,2)
                END SELECT
             END DO
          ELSE ! angle averaged
             SELECT CASE(control%adf_particles)
             CASE(all_rdf_index)
                WRITE(ioindex,'(F30.5,F30.5,F30.5,F30.5)') &
                     binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                     adf(1,ii,1:3)
             CASE(mb_rdf_index)
                WRITE(ioindex,'(F30.5,F30.5)') &
                     binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                     adf(1,ii,1)
             CASE(atom_rdf_index)
                WRITE(ioindex,'(F30.5,F30.5)') &
                     binwidth(1)*(REAL(ii,KIND=dp)-0.5d0),&
                     adf(1,ii,2)
             END SELECT
          END IF
       END DO
       
       WRITE(ioindex,'(A)') ""
       
       ! close output
       CLOSE(ioindex)
    END IF

    RETURN
  END SUBROUTINE write_adf

  ! Writes statistics to a file in the given format. The file must already be open.
  ! *ioindex index for the output (output stream will be ioindex+group_ind)
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *temper temperature
  ! *kinetic kinetic energy
  ! *potential potential energy
  ! *linear linear component of kinetic energy
  ! *rotational rotational component of kinetic energy
  ! *ljpot Lennard-Jones potential energy
  ! *hbpot hydrogen bond potential energy
  ! *mapot mb-atom potential energy
  ! *aapot atom-atom potential energy
  ! *virial the virial
  ! *cell supercell dimensions
  ! *time virtual time
  ! *clock wallclock time
  ! *m_force forces on molecules
  ! *m_torque torques on molecules
  ! *a_force forces on atoms
  ! *control control parameters
  ! *group_ind index of the statistics group (the number of statfile to be written)
  SUBROUTINE write_stats(ioindex,mbs,ats,temper,kinetic,potential,linear,rotational,&
       ljpot,hbpot,mapot,aapot,virial,cell,time,clock,m_force,m_torque,a_force,control,group_ind)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    REAL(KIND=dp), INTENT(IN) :: temper, kinetic, potential, linear, virial, &
         rotational, ljpot, hbpot, mapot, aapot, time, cell(3)
    REAL(KIND=dp), POINTER :: m_force(:,:), m_torque(:,:), a_force(:,:)
    REAL(KIND=dp) :: temparms(3,4), summer(3)
    CHARACTER(LEN=*), INTENT(IN) :: clock
    INTEGER, INTENT(IN) :: ioindex, group_ind
    TYPE(cps), INTENT(IN) :: control
    INTEGER :: ii, jj, kk, n_mbs, n_ats
    
    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)

    ! only the master writes
    IF(cpu_id == master_cpu)THEN
       DO ii = 1, control%n_stats(group_ind) ! loop over stats
          
          SELECT CASE(control%stats(ii,group_ind)) ! select the stat to be written, write on the same line
          CASE(time_stat_index) ! time
             WRITE(ioindex+group_ind,'(A1,F15.4,A1)',advance='no') " ",time," "
          CASE(clock_stat_index) ! wall clock
             WRITE(ioindex+group_ind,'(A)',advance='no') " "//clock//" "
          CASE(temp_stat_index) ! temperature
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",temper," "
          CASE(kin_stat_index) ! kinetic energy
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",kinetic," "
          CASE(pot_stat_index) ! potential energy
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",potential," "
          CASE(ene_stat_index) ! tota energy
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",potential+kinetic," "
          CASE(ljpot_stat_index) ! LJ potential
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",ljpot," "
          CASE(hbpot_stat_index) ! H-bond potential
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",hbpot," "
          CASE(mapot_stat_index) ! Mb-atom potential
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",mapot," "
          CASE(aapot_stat_index) ! atom-atom potential
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",aapot," "
          CASE(virial_stat_index) ! virial
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",virial," "
          CASE(rotational_stat_index) ! rotational kinetic energy
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",rotational," "
          CASE(linear_stat_index) ! linear kinetic energy
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",linear," "
          CASE(volume_stat_index) ! volume
             WRITE(ioindex+group_ind,'(A1,F15.8,A1)',advance='no') " ",cell(1)*cell(2)*cell(3)," "
          CASE(press_stat_index) ! pressure
             WRITE(ioindex+group_ind,'(A1,ES15.8,A1)',advance='no') &
                  " ", (0.66666666667d0*linear + virial)/(cell(1)*cell(2)*cell(3)) ," "
          CASE(coord_stat_index) ! coordinates
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",mbs(control%stat_particles(ii,group_ind))%pos(1),&
                     " ",mbs(control%stat_particles(ii,group_ind))%pos(2),&
                     " ",mbs(control%stat_particles(ii,group_ind))%pos(3)," "
             ELSE IF(control%stat_particles(ii,group_ind) < 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",ats(-control%stat_particles(ii,group_ind))%pos(1),&
                     " ",ats(-control%stat_particles(ii,group_ind))%pos(2),&
                     " ",ats(-control%stat_particles(ii,group_ind))%pos(3)," "
             END IF
          CASE(vel_stat_index) ! velocity
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",mbs(control%stat_particles(ii,group_ind))%vel(1),&
                     " ",mbs(control%stat_particles(ii,group_ind))%vel(2),&
                     " ",mbs(control%stat_particles(ii,group_ind))%vel(3)," "
             ELSE IF(control%stat_particles(ii,group_ind) < 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",ats(-control%stat_particles(ii,group_ind))%vel(1),&
                     " ",ats(-control%stat_particles(ii,group_ind))%vel(2),&
                     " ",ats(-control%stat_particles(ii,group_ind))%vel(3)," "
             END IF
          CASE(orient_stat_index) ! orientation
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F15.8,A1,F15.8,A1,F15.8,A1,F15.8,A1)',advance='no') &
                     " ",mbs(control%stat_particles(ii,group_ind))%orientation%w,&
                     " ",mbs(control%stat_particles(ii,group_ind))%orientation%x,&
                     " ",mbs(control%stat_particles(ii,group_ind))%orientation%y,&
                     " ",mbs(control%stat_particles(ii,group_ind))%orientation%z," "
             END IF
          CASE(ang_stat_index) ! angular velocity
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",mbs(control%stat_particles(ii,group_ind))%angvel(1),&
                     " ",mbs(control%stat_particles(ii,group_ind))%angvel(2),&
                     " ",mbs(control%stat_particles(ii,group_ind))%angvel(3)," "
             END IF
          CASE(arms_stat_index) ! h-bond arms
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                CALL get_arms(mbs(control%stat_particles(ii,group_ind)),temparms)
                DO jj = 1, 4
                   WRITE(ioindex+group_ind,'(A1,F15.8,A1,F15.8,A1,F15.8,A1)',advance='no') &
                        " ",temparms(1,jj),&
                        " ",temparms(2,jj),&
                        " ",temparms(3,jj)," "
                END DO
             END IF
          CASE(force_stat_index) ! force
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",m_force(1,control%stat_particles(ii,group_ind)),&
                     " ",m_force(2,control%stat_particles(ii,group_ind)),&
                     " ",m_force(3,control%stat_particles(ii,group_ind))," "
             ELSE IF(control%stat_particles(ii,group_ind) < 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",a_force(1,-control%stat_particles(ii,group_ind)),&
                     " ",a_force(2,-control%stat_particles(ii,group_ind)),&
                     " ",a_force(3,-control%stat_particles(ii,group_ind))," "
             END IF
          CASE(torque_stat_index) ! torque
             IF(control%stat_particles(ii,group_ind) > 0)THEN
                WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                     " ",m_torque(1,control%stat_particles(ii,group_ind)),&
                     " ",m_torque(2,control%stat_particles(ii,group_ind)),&
                     " ",m_torque(3,control%stat_particles(ii,group_ind))," "
             END IF
          CASE(forcesum_stat_index) ! force sum
             summer = 0.d0
             summing: DO jj = 1, SIZE(control%stat_groups(:,ii,group_ind))
                IF(control%stat_groups(jj,ii,group_ind) == 0)THEN
                   EXIT summing
                END IF
                IF(control%stat_groups(jj,ii,group_ind) > 0)THEN
                   summer = summer + m_force(1:3,control%stat_groups(jj,ii,group_ind))
                ELSE
                   summer = summer + a_force(1:3,-control%stat_groups(jj,ii,group_ind))
                END IF
             END DO summing
             WRITE(ioindex+group_ind,'(A1,F20.8,A1,F20.8,A1,F20.8,A1)',advance='no') &
                  " ",summer(1), &
                  " ",summer(2), &
                  " ",summer(3)," "
             
             
          END SELECT
          
       END DO
       
       ! newline
       WRITE(ioindex+group_ind,'(A)') " "
    END IF
    
    RETURN
  END SUBROUTINE write_stats

  ! Read system geometry from an xyz file. A group of one O follwed by four H will be 
  ! interpreted as an MB molecule, and the orientation is deduced from the relative positions
  ! of the first two Hs with respect to the O. Other particles are interpreted as atomic particles.
  SUBROUTINE read_xyz_particles(filename,geocount,n_elems,mbs,ats,params,found_types)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: geocount
    CHARACTER, POINTER :: data(:,:)
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(mbps), INTENT(IN) :: params
    CHARACTER(LEN=100), POINTER :: tokens(:)
    INTEGER, POINTER :: lwidths(:)
    INTEGER, INTENT(INOUT) :: n_elems
    INTEGER, POINTER :: twidths(:)
    INTEGER :: ii, jj, imbs, iats, iostat, i_geo, i_start, i_end, &
         n_lines, n_tokens, n_mbs, n_ats, t_ind, n_found, xyznumber, &
         running_index
    REAL(KIND=dp) :: c(3), a1(3), a2(3)
    CHARACTER(LEN=100) :: readtag
    LOGICAL :: foundlabel, is_mb
    INTEGER, POINTER :: found_types(:)

    ! it is very inefficient to use these reading routines for a large file.
    ! this should be replaced later
    CALL read_lines(filename,data,lwidths)
    CALL trim(data,lwidths,comment_char)
    CALL lowercase(data)
    n_lines = SIZE(data(1,:))
    n_mbs = 0
    n_ats = 0

    ! find the geometry block to read (if the file contains several)
    i_geo = 0
    i_start = 1
    i_end = 1
    foundlabel = .false.
    DO ii = 1, n_lines
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       IF( is_number(tokens(1)(1:twidths(1))) )THEN
          i_geo = i_geo+1
          IF(i_geo == geocount)THEN
             foundlabel = .true.
             i_start = ii+1
          ELSE IF(i_geo == geocount+1)THEN
             i_end = ii-1
             EXIT
          ELSE IF(i_geo == 1)THEN
             i_start = ii+1
          ELSE IF(i_geo == 2)THEN
             i_end = ii-1
          END IF
       END IF
    END DO
    IF(i_end <= i_start)THEN
       i_end = n_lines
    END IF
    IF(.not.foundlabel)THEN
       WRITE(*,*) "found less geometries in "//filename//" than specified, will read the first one"
    END IF

    ! read the lines and count particles
    ii = i_start
    DO WHILE(ii <= i_end)
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       foundlabel = .false.

       ! an oxygen is found
       IF(tokens(1)(1:twidths(1)) == "o")THEN
          ! search for four following H's to label this as MB
          IF(ii <= i_end-4)THEN
             is_mb = .true.
             DO jj = ii+1, ii+4
                IF(data(1,jj) /= "h" .OR. data(2,jj) /= " ")THEN
                   is_mb = .false.
                END IF
             END DO
          END IF
          ! this is an MB
          IF(is_mb)THEN
             n_mbs = n_mbs + 1
             ii = ii+5
             foundlabel = .true.
          ELSE ! this is just an "O"
             IF(n_elems == 0) CALL abort("Only MB molecules defined, yet other particles found")
             label: DO jj = 1, SIZE(params%atomic_labels(:))
                IF("o" == params%atomic_labels(jj))THEN
                   n_ats = n_ats+1
                   found_types(jj) = found_types(jj) + 1
                   foundlabel = .true.
                   EXIT label
                END IF                
             END DO label             
             ii = ii+1
          END IF
       ELSE ! an arbitrary label is found
          IF(n_elems == 0) CALL abort("Only MB molecules defined, yet other particles found")
          label2: DO jj = 1, SIZE(params%atomic_labels(:))
             IF(tokens(1)(1:twidths(1)) == params%atomic_labels(jj))THEN
                n_ats = n_ats+1
                found_types(jj) = found_types(jj) + 1
                foundlabel = .true.
                EXIT label2
             END IF
          END DO label2
          ii = ii+1
       END IF

       IF(.NOT.foundlabel) CALL abort("undefined label "//tokens(1)(1:twidths(1)))

    END DO

    IF(n_elems > 0)THEN ! check if we found all the elements defined, warn if not
       foundlabel = .true. 
       n_found = 0
       DO ii = 1, n_elems
          IF(found_types(ii) == 0)THEN             
             foundlabel = .false.
          ELSE
             n_found = n_found+1
          END IF
       END DO
       IF(.NOT.foundlabel)THEN
          WRITE(*,*) "Note: some particle types were not found in <"//pos_block//">"
       END IF       
    END IF

    NULLIFY(mbs)
    NULLIFY(ats)
    ALLOCATE(mbs(n_mbs))
    ALLOCATE(ats(n_ats))

    running_index = 1
    imbs = 0
    iats = 0

    ! read the lines and parse particle data
    ii = i_start
    DO WHILE(ii <= i_end)
       CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
       n_tokens = SIZE(tokens(:))

       ! an oxygen is found
       IF(tokens(1)(1:twidths(1)) == "o")THEN
          ! search for four following H's to label this as MB
          IF(ii <= i_end-4)THEN
             is_mb = .true.
             DO jj = ii+1, ii+4
                IF(data(1,jj) /= "h" .OR. data(2,jj) /= " ")THEN
                   is_mb = .false.
                END IF
             END DO
          END IF
          ! this is an MB
          IF(is_mb)THEN

             imbs = imbs + 1
             ! index
             mbs(imbs)%index = running_index
             running_index = running_index + 1

             ! the center
             IF(n_tokens > 3)THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) c(1)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) c(2)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(4)(1:twidths(4)),*,IOSTAT=iostat) c(3)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             ELSE
                CALL abort("Reading xyz geometry in "//filename)
             END IF

             mbs(imbs)%pos = c
             mbs(imbs)%inipos = mbs(imbs)%pos

             ! the first arm
             ii = ii+1
             CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
             n_tokens = SIZE(tokens(:))
             IF(n_tokens > 3)THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) a1(1)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) a1(2)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(4)(1:twidths(4)),*,IOSTAT=iostat) a1(3)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             ELSE
                CALL abort("Reading xyz geometry in "//filename)
             END IF

             ! the second arm
             ii = ii+1
             CALL tokenize(data(1:lwidths(ii),ii)," ",tokens,twidths)
             n_tokens = SIZE(tokens(:))
             IF(n_tokens > 3)THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) a2(1)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) a2(2)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(4)(1:twidths(4)),*,IOSTAT=iostat) a2(3)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             ELSE
                CALL abort("Reading xyz geometry in "//filename)
             END IF

             ! get the quaternion representation for the orientation of the MB
             a1 = a1 - c ! the arms
             a2 = a2 - c
             mbs(imbs)%orientation = get_orientation(a1,a2)

             ! skip the next two lines with the other H's
             ii = ii+3
          ELSE ! this is just an "O"
             
             ! type
             labelb: DO jj = 1, SIZE(params%atomic_labels(:))
                IF(tokens(1)(1:twidths(1)) == params%atomic_labels(jj))THEN
                   iats = iats+1
                   ats(iats)%type = jj
                   ats(iats)%element = params%atomic_labels(jj)
                   EXIT labelb
                END IF                
             END DO labelb

             ! index
             ats(iats)%index = running_index
             running_index = running_index + 1

             ! the center
             IF(n_tokens > 3)THEN
                READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) c(1)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) c(2)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
                READ(tokens(4)(1:twidths(4)),*,IOSTAT=iostat) c(3)
                IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             ELSE
                CALL abort("Reading xyz geometry in "//filename)
             END IF

             ats(iats)%pos = c
             ats(iats)%inipos = ats(iats)%pos
             
             ii = ii+1
          END IF

       ELSE ! an arbitrary label is found

          ! type
          label2b: DO jj = 1, SIZE(params%atomic_labels(:))
             IF(tokens(1)(1:twidths(1)) == params%atomic_labels(jj))THEN
                iats = iats+1
                ats(iats)%type = jj
                ats(iats)%element = params%atomic_labels(jj)
                EXIT label2b
             END IF
          END DO label2b
          
          ! index
          ats(iats)%index = running_index
          running_index = running_index + 1
          
          ! the center
          IF(n_tokens > 3)THEN
             READ(tokens(2)(1:twidths(2)),*,IOSTAT=iostat) c(1)
             IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             READ(tokens(3)(1:twidths(3)),*,IOSTAT=iostat) c(2)
             IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
             READ(tokens(4)(1:twidths(4)),*,IOSTAT=iostat) c(3)
             IF(iostat /= 0) CALL abort("Reading xyz geometry in "//filename)
          ELSE
             CALL abort("Reading xyz geometry in "//filename)
          END IF
          
          ats(iats)%pos = c
          ats(iats)%inipos = ats(iats)%pos
          
          ii = ii+1
       
       END IF
    END DO

    ! assign masses - the current version has same masses for all MB's
    DO ii = 1, n_mbs
       mbs(ii)%m_tot = params%m_mol
       mbs(ii)%m_inert = params%i_mol
    END DO
    DO ii = 1, n_ats
       ats(ii)%mass = params%m_atoms(ats(ii)%type)
    END DO
    
    IF(associated(tokens))THEN
       DEALLOCATE(tokens)
    END IF
    NULLIFY(tokens)
    IF(associated(twidths))THEN
       DEALLOCATE(twidths)
    END IF
    NULLIFY(twidths)
    DEALLOCATE(data)
    NULLIFY(data)
    DEALLOCATE(lwidths)
    NULLIFY(lwidths)
    
  END SUBROUTINE read_xyz_particles


END MODULE file_handler
