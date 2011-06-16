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
! parameters is a module for storing constant parameters.
! it defines physical parameters (e.g. "kb", "pi"), program parameters (e.g. "dp", "norm_tolerance")
! integer indices for various options (e.g. "no_constr_index", "all_frozen_index" etc. are more robust
! names for the possible constraints than 1, 2, etc.), and also two custom types for storing control and model
! parameters. These custom types are used as wrappers for the groups of parameters to allow easy
! passing to subroutines.
! <br /><br />Back to <a href="cashew.html">cashew</a>
! 
! *dp=SELECTED_REAL_KIND double precision
! *pi=3.141592653589793d0 pi
! *twopi=2.d0*pi 2 x pi, 6.283...
! *true_kB=8.61734231E-5 Boltzmann constant in eV/K
! *n_avo Avogadro number
! *sqrt2=sqrt square root of 2
! *sqrt3=sqrt square root of 3
! *true_m_scale=103.64268549109252d0 multiplier for changing from g/mol to eV*fs**2/Ang**2 or g/mol*Ang**2 to eV*fs**2
! *norm_tolerance=1.0E-6 tolerance for treating a norm as zero
! *denom_tolerance=1.0E-3 tolerance for treating a denominator as zero
! *gaussian_limit=8.d0 truncation range for gaussian e**(-x**2/2)
! *scale_max=0.05d0 scaling factor maximum
! *kB=true_kB the value of kB used in the program
! *m_scale=true_m_scale the mass scaling factor used in the program

! *no_constr_index=0 indices for constraint types
! *lenjon_index=1 indices for potential types
! *velocity_verlet_index=1 indices for algorithm types
! *noxyz_index=0 indices for output frequency types
! *microcanonical_index=1 indices for thermostat types
! *test_index=1 indices for run types
! *periodic_bound_index=1 indices for boundary types
! *labelw=5 maximum width for atomic labels
! *ini_temp_def=273.15d0 default values for some of the parameters
! *rdfbins=100 number of bins in rdf histogram
! *mb_rdf_index=1 indices for rdf function
! *time_stat_index=1 indices for statisticss

MODULE parameters
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,100)
!  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,100)
!  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(9,100)  
  REAL(KIND=dp), PARAMETER :: pi = 3.141592653589793d0, &
       twopi = 2.d0*pi, &
       true_kB = 8.61734231E-5, & ! (eV/K)
       n_avo = 6.022E23, &       
       !sqrt2 = sqrt(2.d0), sqrt3 = sqrt(3.d0), &
       sqrt2 = 1.414213562373095d0, sqrt3 = 1.732050807568877d0, &
       true_m_scale = 103.64268549109252d0 ! (g/mol -> eV*fs**2/Å**2 or g/mol*Å**2 -> eV*fs**2)
  REAL(KIND=dp), PARAMETER :: norm_tolerance = 1.0E-6, &
       denom_tolerance = 1.0E-5, &
       gaussian_limit = 8.d0, &
       scale_max = 0.05d0
#ifdef MPI
  CHARACTER(LEN=9), PARAMETER :: version = "   0.32 p"
#else
  CHARACTER(LEN=9), PARAMETER :: version = "     0.32"
#endif

  REAL(KIND=dp) :: kB=true_kB, m_scale=true_m_scale

  ! Constraint indices
  INTEGER, PARAMETER :: &
       no_constr_index = 0, &
       all_frozen_index = 1, &
       frozen_pos_index = 2, &
       harmonic_well_index = 3, &
       frozen_vel_index = 4, &
       ext_force_index = 5

  ! Potential indices
  INTEGER, PARAMETER :: &
       lenjon_index = 1, &
       exp_index = 2, &
       pow_index = 3, &
       spring_index = 4, &
       hard_index = 5, &
       fene_index = 6, &
       hardrep_index = 7, &
       shell_index = 8, &
       max_pot_params = 3

  ! Algorithm indices
  INTEGER, PARAMETER :: &
       velocity_verlet_index = 1, &
       leapfrog_index = 2, &
       cg_index = 3, &
       n_algos = 3 ! number of algorithms
  
  ! output interval indices
  INTEGER, PARAMETER :: &
       noxyz_index = 0, &
       sxyz_index = 1, &
       exyz_index = 2, &
       sexyz_index = 3, &
       ixyz_index = 4, &
       average_index = 5

  ! Thermostat indices
  INTEGER, PARAMETER :: &
       microcanonical_index = 1, &
       langevin_index = 2, &
       cooler_index = 3, &
       n_thermos = 3

  ! Barostat indices
  INTEGER, PARAMETER :: &
       !microcanonical_index = 1, &
       berendsen_index = 2, &
       n_baros = 2

  ! runtypes
  INTEGER, PARAMETER :: &
       test_index = 1, &
       full_index = 2

  ! boundary types
  INTEGER, PARAMETER :: &
       periodic_bound_index = 1, &
       free_bound_index = 2, &
       wall_bound_index = 3

  ! rdf types
  INTEGER, PARAMETER :: rdfbins = 100, &
       mb_rdf_index = 1, &
       atom_rdf_index = 2, &
       all_rdf_index = 3

  ! unit types
  INTEGER, PARAMETER :: true_units = 1, &
       reduced_units = 2

  ! stattypes
  INTEGER, PARAMETER :: &
       time_stat_index = 1, &
       clock_stat_index = 2, &
       temp_stat_index = 3, &
       pot_stat_index = 4, &
       kin_stat_index = 5, &
       press_stat_index = 6, &
       coord_stat_index = 7, &
       orient_stat_index = 8, &
       vel_stat_index = 9, &
       ang_stat_index = 10, &
       arms_stat_index = 11, &
       ene_stat_index = 12, &
       force_stat_index = 13, &
       torque_stat_index = 14, &
       linear_stat_index = 15, &
       rotational_stat_index = 16, &
       ljpot_stat_index = 17, &
       hbpot_stat_index = 18, &
       mapot_stat_index = 19, &
       aapot_stat_index = 20, &
       virial_stat_index = 21, &
       pressure_stat_index = 22, &
       volume_stat_index = 23, &
       forcesum_stat_index = 24, &
       stat_types = 24, &
       stat_digits = 2

  ! axis indices
  INTEGER, PARAMETER :: &
       x_index = 1, &
       y_index = 2, &
       z_index = 3

  ! Max width for atomic labels
  INTEGER, PARAMETER :: labelw = 5
   
  ! default values
  REAL(kind=dp), PARAMETER :: &
       ini_temp_def = 273.15d0, &
       ini_press_def = 1.E-6, &
       L_oh_def = 1.d0, &
       cell_def = 1.d0, &
       baro_interval_def = 10.d0, &
       cg_tol_def = 0.0001d0

  LOGICAL, PARAMETER :: &
       pbc_def = .false.
  
  INTEGER, PARAMETER :: &
       md_algo_def = velocity_verlet_index, &
       md_thermo_def = microcanonical_index, &
       md_baro_def = microcanonical_index, &
       rnd_seed_def = -123456789, &
       xyz_writer_def = noxyz_index, &
       inp_writer_def = exyz_index, &
       rdf_writer_def = noxyz_index, &
       rdf_particles_def = all_rdf_index, &
       verbose_def = -11.d0, &
       verblevel_def = 2, &
       run_def = full_index, &
       units_def = true_units, &
       cg_max_def = 100
  

  ! a custom type storing model parameters
  ! *e_hb hydrogen bond energy
  ! *e_phi dihedral angle energy
  ! *e_lj lennard-jones energy
  ! *r_hb hydrogen bond equilibrium length
  ! *s_lj lennard-jones sigma
  ! *s_r hydrogen bond length sigma
  ! *s_th hydrogen bond angle sigma
  ! *v_b bond order exponent
  ! *r_b bond order range
  ! *d_b bond order transition width
  ! *cut_hb hydrogen bond cutoff
  ! *cut_lj lennard-jones cutoff
  ! *cut_ver neighbor list cutoff (on top of interaction cutoff)
  ! *cut_atom atomic potential cutoff
  ! *m_mol mass of mb molecules
  ! *i_mol moment of inertia of mn molecules
  ! *inv_slj 1/s_lj
  ! *inv_sr 1/s_r
  ! *inv_sth 1/s_th
  ! *inv_db 1/d_b
  ! *inv_m 1/m_mol
  ! *inv_i 1/i_mol
  ! *l_oh O-H distance
  ! *max_cut the maximum cutoff of all interactions
  ! *max_mb_cut the maximum cutoff of mb-mb interactions
  ! *max_ma_cut the maximum cutoff of mb-atom interactions
  ! *max_at_cut the maximum cutoff of atom-atom interactions
  ! *pot_types types of potentials for pairs of particles (#potential,type1,type2)
  ! *n_pots numbers of potentials for pairs of particles (type1,type2)
  ! *pot_params parameters for potentials for pairs of particles (#parameter,#potential,type1,type2)
  ! *m_atoms masses of the different types of atoms
  ! *invm_atoms 1/m_atoms
  ! *pot_cut cutoffs of potentials for pairs of particles (type1,type2)
  ! *atomic_labels names for the atomic types
  TYPE mbps
     REAL(KIND=dp) :: e_hb, e_phi, e_lj, R_hb, &
          s_lj, s_r, s_th, v_b, r_b, d_b, &
          cut_hb, cut_lj, cut_ver, cut_atom, m_mol, I_mol, &          
          inv_slj,inv_sr,inv_sth,inv_db,inv_m,inv_i, &
          L_oh, max_cut, max_mb_cut, max_at_cut, max_ma_cut
     INTEGER, POINTER :: pot_types(:,:,:), n_pots(:,:)
     REAL(KIND=dp), POINTER :: pot_params(:,:,:,:), m_atoms(:), invm_atoms(:), pot_cut(:,:)
     CHARACTER(LEN=labelw), POINTER :: atomic_labels(:)
          
  END TYPE mbps

  ! a custom type storing control parameters
  ! *ini_temp initial (target) temperature
  ! *ini_press target pressure
  ! *time_max simulation length (virtual time) in picoseconds
  ! *time_step length of one (full) simulation timestep in ps
  ! *xyz_interval interval for xyz output (in virtual time) in fs
  ! *inp_interval interval for continuation file output in fs
  ! *rdf_interval interval for continuation file output in fs
  ! *thermo_value parameter for the thermostat (friction coefficient)
  ! *verbose interval for general output in ps
  ! *md_algo index for the molecular dynamics algorithm
  ! *md_thermo index for the thermostat
  ! *rnd_seed seed of the random number generator
  ! *xyz_writer index for xyz output
  ! *runtype index for type of run
  ! *inp_writer index for continuation file output
  ! *rdf_writer index for rdf file output
  ! *rdf_particles index for particle types considered in rdf calculation
  ! *baro_value parameter for barostat
  ! *rdf_range the maximum range for RDF calculation
  ! *stat_interval interval for statistics output in fs
  ! *stat_start starting time for statistics output in fs
  ! *rdf_start starting time for rdf averaging
  ! *stat_start the first point in virtual time in fs for statistics output
  ! *baro_interval interval for applying barostat
  ! *md_baro index for barometer
  ! *n_stats number of stats to be tracked
  ! *n_statfiles number of stats files to be written
  ! *stats indices for stats to be tracked
  ! *stat_particles indices for particles whose stats (e.g. position) are to be tracked
  ! *stat_groups groups of particle indices whose stats are tracked
  ! *scalable logic list, marks if x-, y- and z-directions can be scaled
  ! *unitscale type of units (real/reduced)
  ! *verblevel level of verbosity
  ! *bond_writer index for bond file output
  TYPE cps
     REAL(KIND=dp) :: ini_temp, ini_press, time_max, time_step, xyz_interval, inp_interval, &
          thermo_value, baro_value, verbose, rdf_interval, rdf_start, rdf_range, &
          baro_interval, adf_interval, adf_start, adf_range, adf_center(3), bond_interval
     REAL(KIND=dp), POINTER :: stat_interval(:), stat_start(:)
     INTEGER :: md_algo, md_thermo, md_baro, rnd_seed, xyz_writer, runtype, &
          inp_writer, rdf_writer, rdf_particles, unitscale, n_statfiles, &
          adf_writer, adf_particles, adf_follow, adf_axis, verblevel, bond_writer
     INTEGER, POINTER :: n_stats(:), stats(:,:), stat_particles(:,:), stat_groups(:,:,:)
     LOGICAL :: scalable(3), adf_angle
  END TYPE cps

END MODULE parameters
