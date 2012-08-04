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
! mb_model contains the main physics in cashew.
! it defines the types mb and atom, which denote MB molecules
! and atomic particles, and routines for handling them.
! most importantly, the routines for calculating the energy of
! the simulated system and forces acting on the particles are
! in mb_model. As far as cpu time consumption is concerned,
! the force calculation routine "calc_forces" is the most heavily
! used routine in the program.
! <br /><br />Back to <a href="cashew.html">cashew</a>

! *i1 defines the "no rotation" orientation of one H-bond unit arm of an MB molecule
! *i2 defines the "no rotation" orientation of one H-bond unit arm of an MB molecule
! *i3 defines the "no rotation" orientation of one H-bond unit arm of an MB molecule
! *i4 defines the "no rotation" orientation of one H-bond unit arm of an MB molecule

MODULE mb_model

  USE quaternions
  USE parameters
  USE functions
  USE mpi_mod
  IMPLICIT NONE

  ! unit arms of the mb-molecule at zero rotation
  REAL(KIND=dp), PARAMETER :: &
       i1(3) = (/            0.d0,         0.d0,       1.d0 /), &
       i2(3) = (/ 2.d0*sqrt2/3.d0,         0.d0, -1.d0/3.d0 /), &
       i3(3) = (/     -sqrt2/3.d0,  sqrt2/sqrt3, -1.d0/3.d0 /), &
       i4(3) = (/     -sqrt2/3.d0, -sqrt2/sqrt3, -1.d0/3.d0 /)

  ! the type representing mb molecules
  ! *orientation quaternion representation of the orientation for the molecule
  ! *arms hydrogen bond unit arms (components in a matrix)
  ! *pos position
  ! *vel velocity
  ! *angvel angular velocity
  ! *well parameters for constraining potentials, e.g., spring constant of a harmonic well potential
  ! *inipos initial position of the molecule
  ! *m_tot total mass (currently each particle has its own mass, yet they are the same. this is a waste of memory, but does not affect the speed)
  ! *m_inert total moment of inertia (note that for a tetrahedral molecule, it's a scalar)
  ! *constrained indices defining the types of constraints applied to the molecule in x, y and z directions
  ! *index unique index of the particle
  TYPE mb
     TYPE(qtrn) :: orientation
!     REAL(KIND=dp) :: arms(3,4)
     REAL(KIND=dp) :: pos(3), vel(3), angvel(3), well(3), inipos(3)
     REAL(KIND=dp) :: m_tot, m_inert
     INTEGER :: constrained(3), index
  END TYPE mb

  ! the type representing atomic particles.
  ! since these are point-like particles, they cannot be rotated and therefore
  ! lack the orientation and angular velocity properties of mb molecules
  ! *pos position
  ! *vel velocity
  ! *well parameters for constraining potentials, e.g., spring constant of a harmonic well potential
  ! *inipos initial position of the molecule
  ! *mass total mass
  ! *constrained indices defining the types of constraints applied to the molecule in x, y and z directions
  ! *index unique index of the particle
  ! *element the name of the atom type, e.g., a chemical label
  ! *type index determining the type of the atom (the properties of atoms are defined by types in physical parameters)
  TYPE atom
     REAL(KIND=dp) :: pos(3), vel(3), well(3), inipos(3), mass
     CHARACTER(LEN=labelw) :: element
     INTEGER :: constrained(3), type, index
  END TYPE atom

  TYPE gop
     REAL(KIND=dp) :: pos(3), vel(3), well(3), inipos(3), mass
     INTEGER :: constrained(3), index
  END type gop

  ! interface overloading the function "vector" used for
  ! determining the vector from one particle to another
  ! (taking possible periodic boundaries into account).
  INTERFACE vector
     MODULE PROCEDURE vector_mm, vector_ma, vector_aa, vector_mg, vector_gg
  END INTERFACE

  ! interface overloading the function "displacement"
  ! used for determining the vector from the initial
  ! position of a particle to its current position.
  INTERFACE displacement
     MODULE PROCEDURE displacement_m, displacement_a, displacement_g
  END INTERFACE

CONTAINS

  ! Calculates the hydrogen bond unit arms of the mb-molecule
  ! by rotating the arms of a "non-rotated" molecule
  ! according to the orientation of the given molecule.
  ! *molecule the molecule whose arms are calculated
  ! *j1 one of the arms as a unit vector
  ! *j2 one of the arms as a unit vector
  ! *j3 one of the arms as a unit vector
  ! *j4 one of the arms as a unit vector
  SUBROUTINE get_HB_vectors(molecule,j1,j2,j3,j4)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: molecule
    REAL(KIND=dp), INTENT(OUT) :: j1(3), j2(3), j3(3), j4(3)

    j1 = rotate(i1,molecule%orientation)
    j2 = rotate(i2,molecule%orientation)
    j3 = rotate(i3,molecule%orientation)
    j4 = rotate(i4,molecule%orientation)

  END SUBROUTINE get_HB_vectors

  ! Calculates the hydrogen bond unit arms of the mb-molecule
  ! by rotating the arms of a "non-rotated" molecule
  ! according to the orientation of the given molecule.
  ! *molecule the molecule whose arms are calculated
  ! *arms the arms as a matrix
  SUBROUTINE get_arms(molecule,arms)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: molecule
    REAL(KIND=dp), INTENT(OUT) :: arms(3,4)

    CALL get_HB_vectors(molecule, &
         arms(1:3,1), &
         arms(1:3,2), &
         arms(1:3,3), &
         arms(1:3,4) )

    RETURN
  END SUBROUTINE get_arms

  ! Generates the projections of the hydrogen bond unit arms of a molecule
  ! on the plane perpendicular to the vector rij.
  ! *molecule the molecule whose projected arms are calculated
  ! *j1 one of the arm projection - usually not a unit vector
  ! *j2 one of the arm projection - usually not a unit vector
  ! *j3 one of the arm projection - usually not a unit vector
  ! *j4 one of the arm projection - usually not a unit vector
  ! *rij normal vector of the plane of projection
  SUBROUTINE get_projected_HB_vectors(molecule,rij,j1,j2,j3,j4)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: molecule
    REAL(KIND=dp), INTENT(IN) :: rij(3)
    REAL(KIND=dp), INTENT(OUT) :: j1(3), j2(3), j3(3), j4(3)
    REAL(KIND=dp) :: unit(3), normer

    j1 = rotate(i1,molecule%orientation)
    j2 = rotate(i2,molecule%orientation)
    j3 = rotate(i3,molecule%orientation)
    j4 = rotate(i4,molecule%orientation)
    normer = rij.o.rij
    IF(ABS(normer-1.d0) > norm_tolerance)THEN
       unit = rij/sqrt(normer)
       j1 = j1 - (j1 .o. unit)*unit
       j2 = j2 - (j2 .o. unit)*unit
       j3 = j3 - (j3 .o. unit)*unit
       j4 = j4 - (j4 .o. unit)*unit
    ELSE
       j1 = j1 - (j1 .o. rij)*rij
       j2 = j2 - (j2 .o. rij)*rij
       j3 = j3 - (j3 .o. rij)*rij
       j4 = j4 - (j4 .o. rij)*rij
    END IF

  END SUBROUTINE get_projected_HB_vectors

  ! Rotates the given molecule according to the given
  ! quaternion. The routine updates both the
  ! orientation and the stored unit arms.
  ! *molecule the rotated molecule
  ! *qq representation of the rotation
  SUBROUTINE rotate_molecule(molecule,qq)
    IMPLICIT NONE
    TYPE(mb), INTENT(INOUT) :: molecule
    TYPE(qtrn), INTENT(IN) :: qq

    ! note the order of multiplication
    molecule%orientation = qq*molecule%orientation 
!!$    CALL get_HB_vectors(molecule, &
!!$         molecule%arms(1:3,1), &
!!$         molecule%arms(1:3,2), &
!!$         molecule%arms(1:3,3), &
!!$         molecule%arms(1:3,4) )

  END SUBROUTINE rotate_molecule

  ! Returns the displacement vector of a molecule
  ! with respect to its initial position
  ! *themb the molecule in question
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec the displacement vector from the initial to current position
  FUNCTION displacement_m(themb,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: themb
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = themb%pos(1) - themb%inipos(1)
    vec(2) = themb%pos(2) - themb%inipos(2)
    vec(3) = themb%pos(3) - themb%inipos(3)
!!$    IF(pbc(1))THEN
!!$       IF(2.d0*vec(1) > cell(1))THEN
!!$          vec(1) = vec(1)-cell(1)
!!$       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
!!$          vec(1) = vec(1)+cell(1)
!!$       END IF
!!$    END IF
!!$    IF(pbc(2))THEN
!!$       IF(2.d0*vec(2) > cell(2))THEN
!!$          vec(2) = vec(2)-cell(2)
!!$       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
!!$          vec(2) = vec(2)+cell(2)
!!$       END IF
!!$    END IF
!!$    IF(pbc(3))THEN
!!$       IF(2.d0*vec(3) > cell(3))THEN
!!$          vec(3) = vec(3)-cell(3)
!!$       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
!!$          vec(3) = vec(3)+cell(3)
!!$       END IF
!!$    END IF

  END FUNCTION displacement_m


  ! Returns the displacement vector of a molecule
  ! with respect to its initial position
  ! *thego the molecule in question
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec the displacement vector from the initial to current position
  FUNCTION displacement_g(thego,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(gop), INTENT(IN) :: thego
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = thego%pos(1) - thego%inipos(1)
    vec(2) = thego%pos(2) - thego%inipos(2)
    vec(3) = thego%pos(3) - thego%inipos(3)
!!$    IF(pbc(1))THEN
!!$       IF(2.d0*vec(1) > cell(1))THEN
!!$          vec(1) = vec(1)-cell(1)
!!$       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
!!$          vec(1) = vec(1)+cell(1)
!!$       END IF
!!$    END IF
!!$    IF(pbc(2))THEN
!!$       IF(2.d0*vec(2) > cell(2))THEN
!!$          vec(2) = vec(2)-cell(2)
!!$       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
!!$          vec(2) = vec(2)+cell(2)
!!$       END IF
!!$    END IF
!!$    IF(pbc(3))THEN
!!$       IF(2.d0*vec(3) > cell(3))THEN
!!$          vec(3) = vec(3)-cell(3)
!!$       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
!!$          vec(3) = vec(3)+cell(3)
!!$       END IF
!!$    END IF

  END FUNCTION displacement_g


  ! Returns the displacement vector of an atom
  ! with respect to its initial position
  ! *at the atom in question
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec the displacement vector from the initial to current position
  FUNCTION displacement_a(at,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(atom), INTENT(IN) :: at
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = at%pos(1) - at%inipos(1)
    vec(2) = at%pos(2) - at%inipos(2)
    vec(3) = at%pos(3) - at%inipos(3)
!!$    IF(pbc(1))THEN
!!$       IF(2.d0*vec(1) > cell(1))THEN
!!$          vec(1) = vec(1)-cell(1)
!!$       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
!!$          vec(1) = vec(1)+cell(1)
!!$       END IF
!!$    END IF
!!$    IF(pbc(2))THEN
!!$       IF(2.d0*vec(2) > cell(2))THEN
!!$          vec(2) = vec(2)-cell(2)
!!$       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
!!$          vec(2) = vec(2)+cell(2)
!!$       END IF
!!$    END IF
!!$    IF(pbc(3))THEN
!!$       IF(2.d0*vec(3) > cell(3))THEN
!!$          vec(3) = vec(3)-cell(3)
!!$       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
!!$          vec(3) = vec(3)+cell(3)
!!$       END IF
!!$    END IF

  END FUNCTION displacement_a

  ! Returns the vector to a molecule from the center of mass
  ! of two other molecules, taking into account
  ! the periodic boundaries.
  ! *mb1 position (end) molecule
  ! *mb2 first cm molecules
  ! *mb3 second cm molecule
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from mb2-mb3 cm to mb1
  FUNCTION cm_vector(mb1,mb2,mb3,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1, mb2, mb3
    REAL(KIND=dp) :: vec(3), cm(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    cm = mb2%pos + 0.5d0*vector(mb2,mb3,cell,pbc)
    
    vec(1) = mb1%pos(1) - cm(1)
    vec(2) = mb1%pos(2) - cm(2)
    vec(3) = mb1%pos(3) - cm(3)

    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION cm_vector

  ! Returns the vector from one molecule to another, taking into account
  ! the periodic boundaries.
  ! *mb1 first (start) molecule
  ! *mb2 second (end) molecule
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from mb1 to mb2
  FUNCTION vector_mm(mb1,mb2,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1, mb2
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = mb2%pos(1) - mb1%pos(1)
    vec(2) = mb2%pos(2) - mb1%pos(2)
    vec(3) = mb2%pos(3) - mb1%pos(3)
    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION vector_mm

  ! Returns the vector from a molecule to an atom, taking into account
  ! the periodic boundaries.
  ! *mb1 (start) molecule
  ! *at2 (end) atom
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from mb1 to at2
  FUNCTION vector_ma(mb1,at2,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1
    TYPE(atom), INTENT(IN) :: at2
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = at2%pos(1) - mb1%pos(1)
    vec(2) = at2%pos(2) - mb1%pos(2)
    vec(3) = at2%pos(3) - mb1%pos(3)
    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION vector_ma


  ! Returns the vector from a molecule to a go, taking into account
  ! the periodic boundaries.
  ! *mb1 (start) molecule
  ! *go2 (end) go
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from mb1 to at2
  FUNCTION vector_mg(mb1,go2,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1
    TYPE(gop), INTENT(IN) :: go2
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = go2%pos(1) - mb1%pos(1)
    vec(2) = go2%pos(2) - mb1%pos(2)
    vec(3) = go2%pos(3) - mb1%pos(3)
    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION vector_mg


  ! Returns the vector from a go to a go, taking into account
  ! the periodic boundaries.
  ! *go1 (start) go
  ! *go2 (end) go
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from mb1 to at2
  FUNCTION vector_gg(go1,go2,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(gop), INTENT(IN) :: go1
    TYPE(gop), INTENT(IN) :: go2
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = go2%pos(1) - go1%pos(1)
    vec(2) = go2%pos(2) - go1%pos(2)
    vec(3) = go2%pos(3) - go1%pos(3)
    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION vector_gg




  ! Returns the vector from one atom to another, taking into account
  ! the periodic boundaries.
  ! *at1 first (start) atom
  ! *at2 second (end) atom
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *vec vector from at1 to at2
  FUNCTION vector_aa(at1,at2,cell,pbc) &
       RESULT(vec)
    IMPLICIT NONE
    TYPE(atom), INTENT(IN) :: at1, at2
    REAL(KIND=dp) :: vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec(1) = at2%pos(1) - at1%pos(1)
    vec(2) = at2%pos(2) - at1%pos(2)
    vec(3) = at2%pos(3) - at1%pos(3)
    IF(pbc(1))THEN
       IF(2.d0*vec(1) > cell(1))THEN
          vec(1) = vec(1)-cell(1)
       ELSE IF(2.d0*vec(1) <= -cell(1))THEN
          vec(1) = vec(1)+cell(1)
       END IF
    END IF
    IF(pbc(2))THEN
       IF(2.d0*vec(2) > cell(2))THEN
          vec(2) = vec(2)-cell(2)
       ELSE IF(2.d0*vec(2) <= -cell(2))THEN
          vec(2) = vec(2)+cell(2)
       END IF
    END IF
    IF(pbc(3))THEN
       IF(2.d0*vec(3) > cell(3))THEN
          vec(3) = vec(3)-cell(3)
       ELSE IF(2.d0*vec(3) <= -cell(3))THEN
          vec(3) = vec(3)+cell(3)
       END IF
    END IF

  END FUNCTION vector_aa

  ! Return the distance squared between two molecules 
  ! *mb1 first molecule
  ! *mb2 second molecule
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *dist2 distance squared between mb1 and mb2
  FUNCTION dist_sq(mb1,mb2,cell,pbc) &
       RESULT(dist2)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1, mb2
    REAL(KIND=dp) :: dist2, vec(3)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    vec = vector(mb1,mb2,cell,pbc)
    dist2 = vec(1)**2+vec(2)**2+vec(3)**2

  END FUNCTION dist_sq

  ! Return the distance between two molecules 
  ! *mb1 first molecule
  ! *mb2 second molecule
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *distance distance between mb1 and mb2
  FUNCTION dist(mb1,mb2,cell,pbc) &
       RESULT(distance)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mb1, mb2
    REAL(KIND=dp) :: distance
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)

    distance = sqrt(dist_sq(mb1,mb2,cell,pbc))

  END FUNCTION dist

  ! Initializes arrays storing lists of neighbors and bond numbers
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *mm_nbors MB-MB neighbors list
  ! *ma_nbors MB-atom neighbors list
  ! *aa_nbors atom-atom neighbors list
  ! *mm_n_nbor MB-MB numbers of neighbors
  ! *ma_n_nbor MB-atom numbers of neighbors
  ! *aa_n_nbor atom-atom numbers of neighbors
  ! *params physical parameters
  SUBROUTINE init_data_arrays(mbs,ats,gos,cell,pbc,bonds,boxed,boxes,&
      box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
       mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,&
       mg_nbors,mg_n_nbor,aa_nbors,aa_n_nbor,&
       params,control,go)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_nbor(:), &
         mg_nbors(:,:), mg_n_nbor(:), &
         ma_nbors(:,:), ma_n_nbor(:), aa_nbors(:,:), aa_n_nbor(:), &
         box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), box_gos(:,:,:,:,:),&
         box_mc(:,:,:,:), box_ac(:,:,:,:), box_gc(:,:,:,:)
    INTEGER, INTENT(OUT) :: boxes(3,2)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    LOGICAL, INTENT(OUT) :: boxed
    INTEGER :: ii, jj, nj, n_mb, n_at, n_go, allostat, room, xb, yb, zb, ix, iy, iz, bc(3), nborsize, &
         mm_size, ma_size, mg_size, aa_size, box_size, boxsize
    REAL(KIND=dp) :: distance, vect(3), this_cut, xl(2), yl(2), zl(2), minrange, atomrange
    LOGICAL :: skip

    n_mb = SIZE(mbs(:))
    n_at = SIZE(ats(:))
    n_go = SIZE(gos(:))

    ! allocate arrays
    nborsize = 200
    ALLOCATE(mm_nbors(nborsize,n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(mm_n_nbor(n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(ma_nbors(nborsize,n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(ma_n_nbor(n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(mg_nbors(nborsize,n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(mg_n_nbor(n_mb), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(aa_nbors(nborsize,n_at), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(aa_n_nbor(n_at), STAT=allostat)
    IF(allostat /= 0) CALL abort("allocating data structures")
    ALLOCATE(bonds(n_mb))

    mm_nbors = 0
    mm_n_nbor = 0
    ma_nbors = 0
    ma_n_nbor = 0
    mg_nbors = 0
    mg_n_nbor = 0
    aa_nbors = 0
    aa_n_nbor = 0
    mm_size = nborsize
    ma_size = nborsize
    mg_size = nborsize
    aa_size = nborsize

    ! subcell division:
    ! The supercell may be divided into subcells to find neighbors in O(n) time.
    ! One can go through the list of particles and find the subcell where they
    ! are in O(n) time. Then, one can loop over the particles again and search the nearby
    ! subcells for neighbors in O(n) time. To be efficient, the set of subcells that must be checked
    ! in the neighbor search needs to be smaller than the whole supercell.
    ! Also, if a barostat is used, the cell may become too small.
    ! The barostat problem could easily be circumvented by putting the cell size monitor
    ! in the neighbor update routine - fix?
    minrange = MAXVAL(cell(:))
    IF(params%max_mb_cut > 0.1) minrange = MIN(params%max_mb_cut,minrange)
    IF(params%max_at_cut > 0.1) minrange = MIN(params%max_at_cut,minrange)
    IF(params%max_ma_cut > 0.1) minrange = MIN(params%max_ma_cut,minrange)
    IF(MAXVAL(cell(:)) > 4.d0*minrange .AND. control%md_baro == microcanonical_index)THEN
       boxed = .true.
       IF(params%max_at_cut + params%max_ma_cut > 0.1)THEN
          atomrange = MAX(params%max_at_cut,params%max_ma_cut)
       ELSE
          atomrange = MAXVAL(cell(:))
       END IF
       ! number of subcells in the direction of x, y, z
       boxes(1,1) = MAX(1,FLOOR(cell(1)/params%max_mb_cut))
       boxes(2,1) = MAX(1,FLOOR(cell(2)/params%max_mb_cut))
       boxes(3,1) = MAX(1,FLOOR(cell(3)/params%max_mb_cut))
       boxes(1,2) = MAX(1,FLOOR(cell(1)/atomrange))
       boxes(2,2) = MAX(1,FLOOR(cell(2)/atomrange))
       boxes(3,2) = MAX(1,FLOOR(cell(3)/atomrange))
       ! subcell length in the direction of x, y, z
       xl(1) = cell(1)/REAL(boxes(1,1),KIND=dp)
       yl(1) = cell(2)/REAL(boxes(2,1),KIND=dp)
       zl(1) = cell(3)/REAL(boxes(3,1),KIND=dp)
       xl(2) = cell(1)/REAL(boxes(1,2),KIND=dp)
       yl(2) = cell(2)/REAL(boxes(2,2),KIND=dp)
       zl(2) = cell(3)/REAL(boxes(3,2),KIND=dp)

       boxsize = 200
       IF(n_mb > 0)THEN ! mbs
          ! the lists of molecules in subcells (x,y,z)
          ALLOCATE(box_mbs(boxsize,MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_mbs = 0
          ! the numbers of molecules in subcells (x,y,z)
          ALLOCATE(box_mc(MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_mc = 0
          box_size = boxsize
          DO ii = 1, n_mb ! loop over mbs
             ! find the subcell indices (xb,yb,zb) for the molecule
             CALL box_coordinates(mbs(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_mc(xb,yb,zb,1) = box_mc(xb,yb,zb,1) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_mc(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_mbs,box_size,1)
             END IF
             box_mbs(box_mc(xb,yb,zb,1),xb,yb,zb,1) = ii

             ! find the subcell indices (xb,yb,zb) for the molecule
             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_mc(xb,yb,zb,2) = box_mc(xb,yb,zb,2) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_mc(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_mbs,box_size,1)
             END IF
             box_mbs(box_mc(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF
       box_size = boxsize
       IF(n_go > 0)THEN ! gos
          ! the lists of gos in subcells (x,y,z)
          ALLOCATE(box_gos(boxsize,MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_gos = 0
          ! the numbers of gos in subcells (x,y,z)
          ALLOCATE(box_gc(MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_gc = 0
          box_size = boxsize
          DO ii = 1, n_go ! loop over gos
             ! find the subcell indices (xb,yb,zb) for the go
             CALL box_coordinates(gos(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_gc(xb,yb,zb,1) = box_gc(xb,yb,zb,1) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_gc(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_gos,box_size,1)
             END IF
             box_gos(box_gc(xb,yb,zb,1),xb,yb,zb,1) = ii

             ! find the subcell indices (xb,yb,zb) for the molecule
             CALL box_coordinates(gos(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_gc(xb,yb,zb,2) = box_gc(xb,yb,zb,2) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_gc(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_gos,box_size,1)
             END IF
             box_gos(box_gc(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF
       box_size = boxsize
       IF(n_at > 0)THEN ! atoms
          ! the lists of atoms in subcells (x,y,z)
          ALLOCATE(box_ats(n_at,MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_ats = 0
          ! the numbers of atoms in subcells (x,y,z)
          ALLOCATE(box_ac(MAXVAL(boxes(1,:)),MAXVAL(boxes(2,:)),MAXVAL(boxes(3,:)),2), STAT=allostat)
          IF(allostat /= 0) CALL abort("allocating data structures")
          box_ac = 0
          DO ii = 1, n_at ! loop over atoms
             ! find the subcell indices (xb,yb,zb) for the atom
             CALL box_coordinates(ats(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_ac(xb,yb,zb,1) = box_ac(xb,yb,zb,1) + 1
             IF(box_ac(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_ats,box_size,1)
             END IF
             box_ats(box_ac(xb,yb,zb,1),xb,yb,zb,1) = ii

             ! find the subcell indices (xb,yb,zb) for the atom
             CALL box_coordinates(ats(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_ac(xb,yb,zb,2) = box_ac(xb,yb,zb,2) + 1
             IF(box_ac(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_ats,box_size,2)
             END IF
             box_ats(box_ac(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF

    ELSE
       boxed = .false.
       boxes = 0
    END IF
    

    IF(boxed)THEN ! a subcell division exists

       IF(n_mb > 1)THEN ! mb-mb pairs
          this_cut = params%max_mb_cut ! current cutoff
          DO ii = 1, n_mb-1 ! loop over mbs
             ! find the subcell indices (xb,yb,zb) for the molecule
             CALL box_coordinates(mbs(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             DO ix = -1,1 ! check only the neighboring subcells
                DO iy = -1,1
                   DO iz = -1,1    
                      ! find the neighboring subcell (bc(1:3)) - if the boundaries lead to
                      ! finding the same subcel twice, this function
                      ! notices it and assigns the value .true. to variable
                      ! "skip" only once
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,1),pbc,bc,skip)
                      IF(.not.skip)THEN ! check this subcell
                         DO jj = 1, box_mc(bc(1),bc(2),bc(3),1) ! loop over mbs in the subcell
                            
                            nj = box_mbs(jj,bc(1),bc(2),bc(3),1) ! index of the mb
                            IF(ii < nj)THEN ! prevent double counting
                               ! get the distance between mbs ii and nj
                               vect = vector(mbs(ii),mbs(nj),cell,pbc)
                               distance = .norm.vect
                               IF(distance < this_cut)THEN ! the pair is close to each other
                                  ! store the partners to the lists of neighbors
                                  mm_n_nbor(ii) = mm_n_nbor(ii)+1 ! numbers of neighbors
                                  mm_n_nbor(nj) = mm_n_nbor(nj)+1
                                  IF(mm_n_nbor(ii) >= mm_size-2 .OR. &
                                       mm_n_nbor(jj) >= mm_size-2 )THEN ! if the neighborlist is small, add more space to it
                                     CALL expand_nborlist(mm_nbors,mm_size,1)
                                  END IF
                                  mm_nbors(mm_n_nbor(ii),ii) = nj ! lists of neighbors
                                  mm_nbors(mm_n_nbor(nj),nj) = ii             
                               END IF
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF

       IF(n_go > 0 .AND. n_mb > 0)THEN ! mb-go pairs
          this_cut = params%max_mb_cut ! current cutoff (use mb cutoff)
          DO ii = 1, n_mb
             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip)
                      IF(.not.skip)THEN
                         DO jj = 1, box_gc(bc(1),bc(2),bc(3),2)
                            
                            nj = box_gos(jj,bc(1),bc(2),bc(3),2)
                            vect = vector(mbs(ii),gos(nj),cell,pbc)
                            distance = .norm.vect
                            IF(distance < this_cut)THEN
                               mg_n_nbor(ii) = mg_n_nbor(ii)+1
                               IF(mg_n_nbor(ii) >= mg_size-2)THEN ! if the neighborlist is small, add more space to it
                                  CALL expand_nborlist(mg_nbors,mg_size,1)
                               END IF
                               mg_nbors(mg_n_nbor(ii),ii) = nj
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF


       IF(n_at > 0 .AND. n_mb > 0)THEN ! mb-atom pairs
          !this_cut = params%max_ma_cut ! current cutoff
          DO ii = 1, n_mb
             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip)
                      IF(.not.skip)THEN
                         DO jj = 1, box_ac(bc(1),bc(2),bc(3),2)
                            
                            nj = box_ats(jj,bc(1),bc(2),bc(3),2)
                            vect = vector(mbs(ii),ats(nj),cell,pbc)
                            distance = .norm.vect
                            ! accurate cutoff
                            this_cut = params%cut_ver + params%l_oh + &
                                 MAX(params%pot_cut(0,ats(nj)%type),params%pot_cut(-1,ats(nj)%type))
                            IF(distance < this_cut)THEN
                               ma_n_nbor(ii) = ma_n_nbor(ii)+1
                               IF(ma_n_nbor(ii) >= ma_size-2)THEN ! if the neighborlist is small, add more space to it
                                  CALL expand_nborlist(ma_nbors,ma_size,1)
                               END IF
                               ma_nbors(ma_n_nbor(ii),ii) = nj
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF

       IF(n_at > 1)THEN ! atom-atom pairs
          !this_cut = params%max_at_cut ! current cutoff
          DO ii = 1, n_at             
             CALL box_coordinates(ats(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip)
                      IF(.not.skip)THEN
                         DO jj = 1, box_ac(bc(1),bc(2),bc(3),2)
                            
                            nj = box_ats(jj,bc(1),bc(2),bc(3),2)
                            IF(ii < nj)THEN
                               vect = vector(ats(ii),ats(nj),cell,pbc)
                               distance = .norm.vect
                               ! accurate cutoff
                               this_cut = params%cut_ver + &
                                    params%pot_cut(ats(ii)%type,ats(nj)%type)
                               IF(distance < this_cut)THEN
                                  aa_n_nbor(ii) = aa_n_nbor(ii)+1
                                  aa_n_nbor(nj) = aa_n_nbor(nj)+1
                                  IF(aa_n_nbor(ii) >= aa_size-2 .OR. &
                                       aa_n_nbor(jj) >= aa_size-2 )THEN
                                     CALL expand_nborlist(aa_nbors,aa_size,1)
                                  END IF
                                  aa_nbors(aa_n_nbor(ii),ii) = nj
                                  aa_nbors(aa_n_nbor(nj),nj) = ii             
                               END IF
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO

       END IF

    ELSE ! small system / long interactions with no subcell division

       IF(n_mb > 1)THEN ! mb-mb pairs
          this_cut = params%max_mb_cut ! current cutoff
          DO ii = 1, n_mb-1 ! loop over mbs
             DO jj = ii+1, n_mb ! loop over other mbs ii < jj
                vect = vector(mbs(ii),mbs(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   mm_n_nbor(ii) = mm_n_nbor(ii)+1
                   mm_n_nbor(jj) = mm_n_nbor(jj)+1
                   IF(mm_n_nbor(ii) >= mm_size-2 .OR. &
                        mm_n_nbor(jj) >= mm_size-2 )THEN
                      CALL expand_nborlist(mm_nbors,mm_size,1)
                   END IF
                   mm_nbors(mm_n_nbor(ii),ii) = jj
                   mm_nbors(mm_n_nbor(jj),jj) = ii             
                END IF
             END DO
          END DO
       END IF
       IF(n_go > 0 .AND. n_mb > 0)THEN ! mb-go pairs
          DO ii = 1, n_mb ! loop over mbs
             DO jj = 1, n_go ! loop over gos
                ! accurate cutoff
                this_cut = params%max_mb_cut
                vect = vector(mbs(ii),gos(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   mg_n_nbor(ii) = mg_n_nbor(ii)+1
                   IF(mg_n_nbor(ii) >= mg_size-2)THEN
                      CALL expand_nborlist(mg_nbors,mg_size,1)
                   END IF
                   mg_nbors(mg_n_nbor(ii),ii) = jj
                END IF
             END DO
          END DO
       END IF
       IF(n_at > 0 .AND. n_mb > 0)THEN ! mb-atom pairs
          DO ii = 1, n_mb ! loop over mbs
             DO jj = 1, n_at ! loop over atoms
                ! accurate cutoff
                this_cut = MAX(params%pot_cut(-1,ats(jj)%type),&
                     params%pot_cut(0,ats(jj)%type))+params%cut_ver + params%l_oh
                vect = vector(mbs(ii),ats(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   ma_n_nbor(ii) = ma_n_nbor(ii)+1
                   IF(ma_n_nbor(ii) >= ma_size-2)THEN
                      CALL expand_nborlist(ma_nbors,ma_size,1)
                   END IF
                   ma_nbors(ma_n_nbor(ii),ii) = jj
                END IF
             END DO
          END DO
       END IF
       IF(n_at > 1)THEN ! atom-atom pairs
          !this_cut = params%max_at_cut ! current cutoff
          DO ii = 1, n_at-1 ! loop over atoms
             DO jj = ii+1, n_at ! loop over other atoms (ii < jj)
                ! accurate cutoff
                this_cut = params%pot_cut(ats(ii)%type,ats(jj)%type)+params%cut_ver
                vect = vector(ats(ii),ats(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   aa_n_nbor(ii) = aa_n_nbor(ii)+1
                   aa_n_nbor(jj) = aa_n_nbor(jj)+1
                   IF(aa_n_nbor(ii) >= aa_size-2 .OR. &
                        aa_n_nbor(jj) >= aa_size-2 )THEN
                      CALL expand_nborlist(aa_nbors,aa_size,1)
                   END IF
                   aa_nbors(aa_n_nbor(ii),ii) = jj
                   aa_nbors(aa_n_nbor(jj),jj) = ii
                END IF
             END DO
          END DO
       END IF
    END IF

    ! count the number of neighbors for the mb molecules
    CALL update_bond_numbers(mbs,cell,pbc,mm_nbors,mm_n_nbor,bonds,params%R_b,params%D_b,params%inv_Db)

    RETURN

  END SUBROUTINE init_data_arrays

  ! Finds the coordinates of a point in the box decomposition
  ! *posi position
  ! *xl box length in x direction
  ! *yl box length in y direction
  ! *zl box length in z direction
  ! *xx box index in x direction
  ! *yy box index in y direction
  ! *zz box index in z direction
  SUBROUTINE box_coordinates(posi,xl,yl,zl,xx,yy,zz,pbc,boxes)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: posi(3),xl,yl,zl 
    INTEGER, INTENT(IN) :: boxes(3)
    INTEGER, INTENT(OUT) :: xx, yy, zz
    LOGICAL, INTENT(IN) :: pbc(3)
    
    IF(pbc(1))THEN
       xx = FLOOR(posi(1)/xl)+1
    ELSE
       xx = MIN(boxes(1),MAX(1,FLOOR(posi(1)/xl)+1))
    END IF
    IF(pbc(2))THEN
       yy = FLOOR(posi(2)/yl)+1
    ELSE
       yy = MIN(boxes(2),MAX(1,FLOOR(posi(2)/yl)+1))
    END IF
    IF(pbc(3))THEN
       zz = FLOOR(posi(3)/zl)+1
    ELSE
       zz = MIN(boxes(3),MAX(1,FLOOR(posi(3)/zl)+1))
    END IF

    RETURN
  END SUBROUTINE box_coordinates

  ! Gives the subcell division coordinates of the cell found when moving (ix,iy,iz) cells 
  ! starting from the cell (xb,yb,zb). For finding neighbors, one usually
  ! loops ix = -1, 1; iy = -1, 1; iz = -1, 1 - but this can result in
  ! finding the same cells several times. The logical variable "skip" keeps track of this:
  ! It will be given the value true only once if such a loop finds a cell many times.
  ! *xb x-coordinate of the origin
  ! *yb y-coordinate of the origin
  ! *zb z-coordinate of the origin
  ! *ix shift in x-direction
  ! *iy shift in y-direction
  ! *iz shift in z-direction
  ! *boxes number of cells in x-, y- and z-directions
  ! *pbc periodic bounds
  ! *bc shifted coordinates
  ! *skip logical tag that tracks if neighbor looping from the origin finds the same cells multiple times (is true for only one of them)
  SUBROUTINE box_shift(xb,yb,zb,ix,iy,iz,boxes,pbc,bc,skip)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: xb,yb,zb,ix,iy,iz,boxes(3)
    INTEGER, INTENT(OUT) :: bc(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    LOGICAL, INTENT(INOUT) :: skip
    INTEGER :: jj, orig(3), shift(3)

    ! shift box coordinates by the given amount
    bc(1) = xb+ix
    bc(2) = yb+iy
    bc(3) = zb+iz
    orig = (/xb,yb,zb/)
    shift = (/ix,iy,iz/)
    skip = .false.
    DO jj = 1,3 
       IF(pbc(jj))THEN ! shift over periodic boundary
          DO WHILE(bc(jj) < 1)
             bc(jj) = bc(jj) + boxes(jj)
          END DO
          ! make sure we don't double count over the periodic boundary
          IF( (bc(jj) == orig(jj) .OR. bc(jj) == orig(jj)+1) .AND. shift(jj) < 0)THEN
             skip = .true.
             RETURN
          END IF
          DO WHILE(bc(jj) > boxes(jj))
             bc(jj) = bc(jj) - boxes(jj)
          END DO
          IF( (bc(jj) == orig(jj) .OR. bc(jj) == orig(jj)-1) .AND. shift(jj) > 0 )THEN
             skip = .true.
             RETURN
          END IF
       ELSE ! if the boundary is aperiodic, we cannot move further than to the last box          
          IF(bc(jj) < 1)THEN
             bc(jj) = 1
             skip = .true.  ! this cell will already be found when the shift is zero      
          ELSE IF(bc(jj) > boxes(jj))THEN
             bc(jj) = boxes(jj)
             skip = .true.             
          END IF
       END IF
    END DO
       
    RETURN
  END SUBROUTINE box_shift

  ! Increases the size of a list of neighbors (in case the particles accumulate)
  ! *thelist a 2D array to be expanded: the first index is the one where space is added
  ! *thesize should contain the size of the list, will be assigned the new length
  ! *tenths the factor of expansion in tens of percents
  SUBROUTINE expand_nborlist(thelist,thesize,tenths)
    IMPLICIT NONE
    INTEGER, POINTER :: thelist(:,:)
    INTEGER, ALLOCATABLE :: templist(:,:)
    INTEGER, INTENT(IN) :: tenths
    INTEGER, INTENT(INOUT) :: thesize
    INTEGER :: length
        
    length = SIZE(thelist(1,:))
    ALLOCATE(templist(thesize,length))
    templist(:,:) = thelist(:,:)
    IF(associated(thelist))THEN
       DEALLOCATE(thelist)
    END IF
    NULLIFY(thelist)
    ALLOCATE(thelist(((10+tenths)*thesize)/10,length))
    thelist = 0
    thelist(1:thesize,1:length) = templist(1:thesize,1:length)
    thesize = ((10+tenths)*thesize)/10
    DEALLOCATE(templist)

    RETURN
  END SUBROUTINE expand_nborlist

  ! Increases the size of a list of neighbors to a fixed length (in case the particles accumulate)
  ! *thelist a 2D array to be expanded: the first index is the one where space is added
  ! *thesize should contain the new size of the list
  SUBROUTINE fix_nborlist(thelist,thesize)
    IMPLICIT NONE
    INTEGER, POINTER :: thelist(:,:)
    INTEGER, ALLOCATABLE :: templist(:,:)
    INTEGER, INTENT(IN) :: thesize
    INTEGER :: length, width
        
    length = SIZE(thelist(1,:))
    width = SIZE(thelist(:,1))
    IF(thesize <= width) RETURN
    ALLOCATE(templist(width,length))
    templist(:,:) = thelist(:,:)
    IF(associated(thelist))THEN
       DEALLOCATE(thelist)
    END IF
    NULLIFY(thelist)
    ALLOCATE(thelist(thesize,length))
    thelist = 0
    thelist(1:width,1:length) = templist(1:width,1:length)
    DEALLOCATE(templist)

    RETURN
  END SUBROUTINE fix_nborlist

  ! Increases the size of a list of particles in domain boxes (in case the particles accumulate)
  ! *thelist a 4D array to be expanded: the first index is the one where space is added
  ! *thesize should contain the size of the list, will be assigned the new length
  ! *tenths the factor of expansion in tens of percents
  SUBROUTINE expand_boxlist(thelist,thesize,tenths)
    IMPLICIT NONE
    INTEGER, POINTER :: thelist(:,:,:,:,:)
    INTEGER, ALLOCATABLE :: templist(:,:,:,:,:)
    INTEGER, INTENT(IN) :: tenths
    INTEGER, INTENT(INOUT) :: thesize
    INTEGER :: xl,yl,zl,nlist
        
    xl = SIZE(thelist(1,:,1,1,1))
    yl = SIZE(thelist(1,1,:,1,1))
    zl = SIZE(thelist(1,1,1,:,1))
    nlist = SIZE(thelist(1,1,1,1,:))
    ALLOCATE(templist(thesize,xl,yl,zl,nlist))
    templist(:,:,:,:,:) = thelist(:,:,:,:,:)
    IF(associated(thelist))THEN
       DEALLOCATE(thelist)
    END IF
    NULLIFY(thelist)
    ALLOCATE(thelist(((10+tenths)*thesize)/10,xl,yl,zl,nlist))
    thelist = 0
    thelist(1:thesize,1:xl,1:yl,1:zl,1:nlist) = templist(1:thesize,1:xl,1:yl,1:zl,1:nlist)
    thesize = ((10+tenths)*thesize)/10
    DEALLOCATE(templist)

    RETURN
  END SUBROUTINE expand_boxlist

  ! Updates neighbor lists and counts, but does not initialize the lists
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *mm_nbors MB-MB neighbors list
  ! *ma_nbors MB-atom neighbors list
  ! *aa_nbors atom-atom neighbors list
  ! *mm_n_nbor MB-MB numbers of neighbors
  ! *ma_n_nbor MB-atom numbers of neighbors
  ! *aa_n_nbor atom-atom numbers of neighbors
  ! *params physical parameters
  SUBROUTINE update_neighbors(mbs,ats,gos,cell,pbc,bonds,boxed,boxes,&
       box_mbs,box_ats,box_gos,box_mc,box_ac,box_gc,&
       mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,aa_nbors,aa_n_nbor,mg_nbors,mg_n_nbor,&
       params,go)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_nbor(:), &
         ma_nbors(:,:), ma_n_nbor(:), &
         aa_nbors(:,:), aa_n_nbor(:), &
         mg_nbors(:,:), mg_n_nbor(:), &
         box_mbs(:,:,:,:,:), box_ats(:,:,:,:,:), box_gos(:,:,:,:,:), &
         box_mc(:,:,:,:), box_ac(:,:,:,:), box_gc(:,:,:,:)
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    INTEGER, INTENT(IN) :: boxes(3,2)
    LOGICAL, INTENT(IN) :: pbc(3), boxed
    INTEGER :: ii, jj, nj, n_mb, n_at, n_go, &
         mm_size, ma_size, aa_size, mg_size, box_size, &
         xb, yb, zb, ix, iy, iz, bc(3), &
         mpi_mm_size, mpi_ma_size, mpi_aa_size, mpi_mg_size
    REAL(KIND=dp) :: distance, vect(3), this_cut, xl(2), yl(2), zl(2)
    LOGICAL :: skip

    n_mb = SIZE(mbs(:))
    n_at = SIZE(ats(:))
    n_go = SIZE(gos(:))
    mm_size = SIZE(mm_nbors(:,1))
    ma_size = SIZE(ma_nbors(:,1))
    aa_size = SIZE(aa_nbors(:,1))
    mg_size = SIZE(mg_nbors(:,1))
    mm_nbors = 0
    mm_n_nbor = 0
    ma_nbors = 0
    ma_n_nbor = 0
    aa_nbors = 0
    aa_n_nbor = 0   
    mg_nbors = 0
    mg_n_nbor = 0

    mpi_mm_size = mm_size
    mpi_ma_size = ma_size
    mpi_aa_size = aa_size
    mpi_mg_size = mg_size
    mpi_m_n_nbor = 0
    mpi_a_n_nbor = 0

    ! the subcell decomposition is initialized in init_data_arrays(...), and
    ! the way the division works is explained there
    IF(boxed)THEN ! subcell decomposition
       ! update subcell occupations - this has not been parallellized

       ! cell lengths
       xl(1) = cell(1)/REAL(boxes(1,1),KIND=dp)
       yl(1) = cell(2)/REAL(boxes(2,1),KIND=dp)
       zl(1) = cell(3)/REAL(boxes(3,1),KIND=dp)
       xl(2) = cell(1)/REAL(boxes(1,2),KIND=dp)
       yl(2) = cell(2)/REAL(boxes(2,2),KIND=dp)
       zl(2) = cell(3)/REAL(boxes(3,2),KIND=dp)

       IF(n_mb > 0)THEN ! mbs
          box_mbs = 0
          box_mc = 0
          box_size = SIZE(box_mbs(:,1,1,1,1))
          DO ii = 1, n_mb ! loop over mbs
             ! get the subcell indices
             CALL box_coordinates(mbs(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_mc(xb,yb,zb,1) = box_mc(xb,yb,zb,1) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_mc(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_mbs,box_size,1)
             END IF
             box_mbs(box_mc(xb,yb,zb,1),xb,yb,zb,1) = ii

             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_mc(xb,yb,zb,2) = box_mc(xb,yb,zb,2) + 1
             ! if the allocated space for the lists of particles is filling, expand the available space
             IF(box_mc(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_mbs,box_size,1)
             END IF
             box_mbs(box_mc(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF
       IF(n_at > 0)THEN ! atoms
          box_ats = 0
          box_ac = 0
          box_size = SIZE(box_ats(:,1,1,1,1))
          DO ii = 1, n_at
             CALL box_coordinates(ats(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_ac(xb,yb,zb,1) = box_ac(xb,yb,zb,1) + 1
             IF(box_ac(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_ats,box_size,1)
             END IF
             box_ats(box_ac(xb,yb,zb,1),xb,yb,zb,1) = ii

             CALL box_coordinates(ats(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_ac(xb,yb,zb,2) = box_ac(xb,yb,zb,2) + 1
             IF(box_ac(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_ats,box_size,1)
             END IF
             box_ats(box_ac(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF
       IF(n_at > 0)THEN ! gos
          box_gos = 0
          box_gc = 0
          box_size = SIZE(box_gos(:,1,1,1,1))
          DO ii = 1, n_at
             CALL box_coordinates(gos(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             box_gc(xb,yb,zb,1) = box_gc(xb,yb,zb,1) + 1
             IF(box_gc(xb,yb,zb,1) >= box_size-2)THEN
                CALL expand_boxlist(box_gos,box_size,1)
             END IF
             box_gos(box_gc(xb,yb,zb,1),xb,yb,zb,1) = ii

             CALL box_coordinates(gos(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             box_gc(xb,yb,zb,2) = box_gc(xb,yb,zb,2) + 1
             IF(box_gc(xb,yb,zb,2) >= box_size-2)THEN
                CALL expand_boxlist(box_gos,box_size,1)
             END IF
             box_gos(box_gc(xb,yb,zb,2),xb,yb,zb,2) = ii
          END DO
       END IF

       ! update neighbors
       
       ! mb-mb list
       this_cut = params%max_mb_cut ! current cutoff, start with mb-mb
       IF(first_mpi_mb > 0)THEN
          !DO ii = 1, n_mb-1
          DO ii = first_mpi_mb, last_mpi_mb ! loop over the molecules dedicated for this processor
             CALL box_coordinates(mbs(ii)%pos,xl(1),yl(1),zl(1),xb,yb,zb,pbc,boxes(:,1))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,1),pbc,bc,skip) ! get the shifted box coordinates
                      IF(.not.skip)THEN ! the neighbor loop may find the same box several times at boundaries
                                        ! but each cell must be inspected only once
                         DO jj = 1, box_mc(bc(1),bc(2),bc(3),1) ! loop over the mbs in the box
                         
                            nj = box_mbs(jj,bc(1),bc(2),bc(3),1) ! nj is the index of the mb
                            IF(ii < nj)THEN
                               vect = vector(mbs(ii),mbs(nj),cell,pbc)
                               distance = .norm.vect
                               IF(distance < this_cut)THEN ! if close enough, add to the neighbor list
                                  mm_n_nbor(ii) = mm_n_nbor(ii)+1
                                  mm_n_nbor(nj) = mm_n_nbor(nj)+1
                                  IF(mm_n_nbor(ii) >= mpi_mm_size-2 .OR. &
                                       mm_n_nbor(jj) >= mpi_mm_size-2 )THEN ! if the neighborlist is small, add more space to it
                                     CALL expand_nborlist(mm_nbors,mpi_mm_size,1)
                                  END IF
                                  mm_nbors(mm_n_nbor(ii),ii) = nj
                                  mm_nbors(mm_n_nbor(nj),nj) = ii             
                               END IF
                            END IF

                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF
       
       ! mb-atom list
       IF(first_mpi_mb > 0 .AND. n_as > 0)THEN
          !this_cut = params%max_ma_cut
          DO ii = first_mpi_mb, last_mpi_mb ! loop over the molecules dedicated for this processor
             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip) ! get the shifted box coordinates
                      IF(.not.skip)THEN ! the neighbor loop may find the same box several times at boundaries
                                        ! but each cell must be inspected only once
                         DO jj = 1, box_ac(bc(1),bc(2),bc(3),2) ! loop over the atoms in the box
                            ! accurate cutoff
                            this_cut = MAX(params%pot_cut(-1,ats(jj)%type),&
                                 params%pot_cut(0,ats(jj)%type))+params%cut_ver+params%l_oh
                            nj = box_ats(jj,bc(1),bc(2),bc(3),2)
                            vect = vector(mbs(ii),ats(nj),cell,pbc)
                            distance = .norm.vect
                            IF(distance < this_cut)THEN ! if close enough, add to the neighbor list
                               ma_n_nbor(ii) = ma_n_nbor(ii)+1
                               IF(ma_n_nbor(ii) >= mpi_ma_size-2)THEN ! if the neighborlist is small, add more space to it
                                  CALL expand_nborlist(ma_nbors,mpi_ma_size,1)
                               END IF
                               ma_nbors(ma_n_nbor(ii),ii) = nj
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF



       ! mb-go list
       IF(first_mpi_mb > 0 .AND. n_gs > 0)THEN
          !this_cut = params%max_ma_cut
          DO ii = first_mpi_mb, last_mpi_mb ! loop over the molecules dedicated for this processor
             CALL box_coordinates(mbs(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip) ! get the shifted box coordinates
                      IF(.not.skip)THEN ! the neighbor loop may find the same box several times at boundaries
                                        ! but each cell must be inspected only once
                         DO jj = 1, box_gc(bc(1),bc(2),bc(3),2) ! loop over the gos in the box
                            ! accurate cutoff
                            this_cut = params%cut_lj + params%cut_ver
                            nj = box_gos(jj,bc(1),bc(2),bc(3),2)
                            vect = vector(mbs(ii),gos(nj),cell,pbc)
                            distance = .norm.vect
                            IF(distance < this_cut)THEN ! if close enough, add to the neighbor list
                               mg_n_nbor(ii) = mg_n_nbor(ii)+1
                               IF(mg_n_nbor(ii) >= mpi_mg_size-2)THEN ! if the neighborlist is small, add more space to it
                                  CALL expand_nborlist(mg_nbors,mpi_mg_size,1)
                               END IF
                               mg_nbors(mg_n_nbor(ii),ii) = nj
                            END IF
                            
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END IF



       ! mb-atom list
       IF(first_mpi_atom > 0)THEN
          !this_cut = params%max_at_cut
          DO ii = first_mpi_atom, last_mpi_atom
             CALL box_coordinates(ats(ii)%pos,xl(2),yl(2),zl(2),xb,yb,zb,pbc,boxes(:,2))
             DO ix = -1,1 ! check only the neighboring boxes
                DO iy = -1,1
                   DO iz = -1,1
                      CALL box_shift(xb,yb,zb,ix,iy,iz,boxes(:,2),pbc,bc,skip)
                      IF(.not.skip)THEN
                         DO jj = 1, box_ac(bc(1),bc(2),bc(3),2)
                            nj = box_ats(jj,bc(1),bc(2),bc(3),2)
                            IF(ii < nj)THEN
                               vect = vector(ats(ii),ats(nj),cell,pbc)
                               distance = .norm.vect
                               ! accurate cutoff
                               this_cut = params%pot_cut(ats(ii)%type,ats(nj)%type)+&
                                    params%cut_ver
                               IF(distance < this_cut)THEN
                                  aa_n_nbor(ii) = aa_n_nbor(ii)+1
                                  aa_n_nbor(nj) = aa_n_nbor(nj)+1
                                  IF(aa_n_nbor(ii) >= mpi_aa_size-2 .OR. &
                                       aa_n_nbor(jj) >= mpi_aa_size-2 )THEN
                                     CALL expand_nborlist(aa_nbors,mpi_aa_size,1)
                                  END IF
                                  aa_nbors(aa_n_nbor(ii),ii) = nj
                                  aa_nbors(aa_n_nbor(nj),nj) = ii
                               END IF
                            END IF
                            
                         END DO
                      END IF
                   END DO
                END DO
             END DO
          END DO
          
       END IF

    ELSE ! no subcell decomposition

       IF(first_mpi_mb > 0)THEN ! this cpu should work on mbs
          this_cut = params%max_mb_cut ! current cutoff, start with mb-mb
          DO ii = first_mpi_mb, last_mpi_mb  ! loop over the mbs of this cpu
             DO jj = ii+1, n_mb ! loop over other mbs ii < jj
                vect = vector(mbs(ii),mbs(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   mm_n_nbor(ii) = mm_n_nbor(ii)+1
                   mm_n_nbor(jj) = mm_n_nbor(jj)+1
                   IF(mm_n_nbor(ii) >= mpi_mm_size-2 .OR. &
                        mm_n_nbor(jj) >= mpi_mm_size-2 )THEN
                      CALL expand_nborlist(mm_nbors,mpi_mm_size,1)
                   END IF
                   mm_nbors(mm_n_nbor(ii),ii) = jj
                   mm_nbors(mm_n_nbor(jj),jj) = ii             
                END IF
             END DO
          END DO
       !END IF

       !IF(first_mpi_mb > 0)THEN ! this cpu should work on mbs (mb-atom pairs)
          !this_cut = params%max_ma_cut ! current cutoff
          DO ii = first_mpi_mb, last_mpi_mb
             DO jj = 1, n_at
                ! accurate cutoff
                this_cut = MAX(params%pot_cut(-1,ats(jj)%type),&
                     params%pot_cut(0,ats(jj)%type))+params%cut_ver+params%l_oh
                vect = vector(mbs(ii),ats(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   ma_n_nbor(ii) = ma_n_nbor(ii)+1
                   IF(ma_n_nbor(ii) >= mpi_ma_size-2)THEN
                      CALL expand_nborlist(ma_nbors,mpi_ma_size,1)
                   END IF
                   ma_nbors(ma_n_nbor(ii),ii) = jj
                END IF
             END DO
          END DO

          DO ii = first_mpi_mb, last_mpi_mb
             DO jj = 1, n_go
                ! accurate cutoff
                this_cut = params%cut_lj + params%cut_ver
                vect = vector(mbs(ii),gos(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   mg_n_nbor(ii) = mg_n_nbor(ii)+1
                   IF(mg_n_nbor(ii) >= mpi_mg_size-2)THEN
                      CALL expand_nborlist(mg_nbors,mpi_mg_size,1)
                   END IF
                   mg_nbors(mg_n_nbor(ii),ii) = jj
                END IF
             END DO
          END DO
       END IF



       IF(first_mpi_atom > 0)THEN ! this cpu should work on atoms
          !this_cut = params%max_at_cut ! current cutoff
          DO ii = first_mpi_atom, last_mpi_atom       
             DO jj = ii+1, n_at
                ! accurate cutoff
                this_cut = params%pot_cut(ats(ii)%type,ats(jj)%type)+params%cut_ver
                vect = vector(ats(ii),ats(jj),cell,pbc)
                distance = .norm.vect
                IF(distance < this_cut)THEN
                   aa_n_nbor(ii) = aa_n_nbor(ii)+1
                   aa_n_nbor(jj) = aa_n_nbor(jj)+1
                   IF(aa_n_nbor(ii) >= mpi_aa_size-2 .OR. &
                        aa_n_nbor(jj) >= mpi_aa_size-2 )THEN
                      CALL expand_nborlist(aa_nbors,mpi_aa_size,1)
                   END IF
                   aa_nbors(aa_n_nbor(ii),ii) = jj
                   aa_nbors(aa_n_nbor(jj),jj) = ii             
                END IF
             END DO
          END DO
       END IF
    END IF



#ifdef MPI

    ! gather lists and expand if necessary - mb-mb
    IF(n_ms > 1)THEN

       ! sum the numbers of neighbors found by each cpu to get the totals
       CALL MPI_ALLREDUCE(mm_n_nbor,mpi_m_n_nbor,n_ms,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpistat)
       ! the size of the neighbor list must be at least the current size, and big enough to fit all the found neighbors
       mpi_mm_size = MAX(mpi_mm_size,5*MAXVAL(mpi_m_n_nbor(:))/4)
       ! find and send the maximum of all mm_sizes (if some cpu had its lists expanded during search and therefore
       ! has a bigger mm_size than others, we must expand all lists to match)
       CALL MPI_ALLREDUCE(mpi_mm_size,mm_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpistat)
       ! if some process expanded the lists, sync the list lengths
       CALL fix_nborlist(mm_nbors,mm_size)
       
       ! stack the lists and broadcast them to all cpus
       CALL mpi_stack(mm_nbors,mm_n_nbor,n_ms,mm_size) ! a custom routine for stacking 2d arrays from all cpus, defined in this module
       CALL MPI_BCAST(mm_nbors,n_ms*mm_size,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
       CALL MPI_BCAST(mm_n_nbor,n_ms,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF

    ! mb-atom
    IF(n_ms > 0 .AND. n_as > 0)THEN
       CALL MPI_ALLREDUCE(ma_n_nbor,mpi_m_n_nbor,n_ms,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpistat)
       mpi_ma_size = MAX(mpi_ma_size,5*MAXVAL(mpi_m_n_nbor(:))/4)
       CALL MPI_ALLREDUCE(mpi_ma_size,ma_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpistat)
       CALL fix_nborlist(ma_nbors,ma_size)
       
       CALL mpi_stack(ma_nbors,ma_n_nbor,n_ms,ma_size)
       CALL MPI_BCAST(ma_nbors,n_ms*ma_size,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
       CALL MPI_BCAST(ma_n_nbor,n_ms,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF

    ! mb-go
    IF(n_ms > 0 .AND. n_gs > 0)THEN
       CALL MPI_ALLREDUCE(mg_n_nbor,mpi_m_n_nbor,n_ms,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpistat)
       mpi_mg_size = MAX(mpi_mg_size,5*MAXVAL(mpi_m_n_nbor(:))/4)
       CALL MPI_ALLREDUCE(mpi_mg_size,mg_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpistat)
       CALL fix_nborlist(mg_nbors,mg_size)
       
       CALL mpi_stack(mg_nbors,mg_n_nbor,n_ms,mg_size)
       CALL MPI_BCAST(mg_nbors,n_ms*mg_size,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
       CALL MPI_BCAST(mg_n_nbor,n_ms,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF
    
    ! atom-atom
    IF(n_as > 1)THEN
       CALL MPI_ALLREDUCE(aa_n_nbor,mpi_a_n_nbor,n_as,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpistat)
       mpi_aa_size = MAX(mpi_aa_size,5*MAXVAL(mpi_a_n_nbor(:))/4)
       CALL MPI_ALLREDUCE(mpi_aa_size,aa_size,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpistat)
       CALL fix_nborlist(aa_nbors,aa_size)
       
       CALL mpi_stack(aa_nbors,aa_n_nbor,n_as,aa_size)
       CALL MPI_BCAST(aa_nbors,n_as*aa_size,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
       CALL MPI_BCAST(aa_n_nbor,n_as,MPI_INTEGER,master_cpu,MPI_COMM_WORLD,mpistat)
    END IF

#endif
    
    ! recalculate bond numbers for mb tersoff coefficients
    CALL update_bond_numbers_p(mbs,cell,pbc,mm_nbors,mm_n_nbor,bonds,params%R_b,params%D_b,params%inv_Db)

    RETURN

  END SUBROUTINE update_neighbors

  ! Updates the bond counts (for bond-order term) of MB molecules
  ! *mbs list of molecules
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *nbors MB-MB neighbors list
  ! *n_ns MB-MB numbers of neighbors
  ! *Rb the range for bonding (molecules at distance Rb have 0.5 of a bond)
  ! *Db (half of) the range for changing the bonding value from 1.0 to 0.0
  ! *inv_Db 1/Db
  SUBROUTINE update_bond_numbers(mbs,cell,pbc,nbors,n_ns,bonds,Rb,Db,inv_Db)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: nbors(:,:), n_ns(:)
    REAL(KIND=dp), INTENT(IN) :: Rb,Db,inv_Db
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    INTEGER :: ii, jj, kk, nn, n_mbs, ll
    REAL(KIND=dp) :: distance, vect(3)

    n_mbs = SIZE(n_ns(:))
    DO ii = 1, n_mbs ! loop over mbs
       nn = n_ns(ii)
       bonds(ii) = 0.d0
       DO jj = 1, nn ! loop over neighbors
          kk = nbors(jj,ii) ! kk is the index of the jj:th neighbor of ii
          vect = vector(mbs(ii),mbs(kk),cell,pbc)
          distance = .norm.vect
          IF(distance < Rb+Db)THEN ! kk is close enouhg to ii to contribute
             bonds(ii) = bonds(ii) + proximity(distance,Rb,Db,inv_Db) ! add the proximity function
          END IF    
       END DO       
    END DO

    RETURN
  END SUBROUTINE update_bond_numbers

  ! Updates the bond counts (for bond-order term) of MB molecules using a parallel scheme
  ! *mbs list of molecules
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *nbors MB-MB neighbors list
  ! *n_ns MB-MB numbers of neighbors
  ! *Rb the range for bonding (molecules at distance Rb have 0.5 of a bond)
  ! *Db (half of) the range for changing the bonding value from 1.0 to 0.0
  ! *inv_Db 1/Db
  SUBROUTINE update_bond_numbers_p(mbs,cell,pbc,nbors,n_ns,bonds,Rb,Db,inv_Db)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: nbors(:,:), n_ns(:)
    REAL(KIND=dp), INTENT(IN) :: Rb,Db,inv_Db
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    INTEGER :: ii, jj, kk, nn, n_mbs, ll, index
    REAL(KIND=dp) :: distance, vect(3)

    !n_mbs = SIZE(n_ns(:))
    bonds = 0.d0
    mpi_bonds = 0.d0
    ! equal split
    DO index = 1, bondcount(cpu_id+1) ! handle #bondcount mbs
       ii = bondpos(cpu_id+1)+index ! start from #bondpos + 1
       nn = n_ns(ii)
       DO jj = 1, nn
          kk = nbors(jj,ii)
          vect = vector(mbs(ii),mbs(kk),cell,pbc)
          distance = .norm.vect
          IF(distance < Rb+Db)THEN
             mpi_bonds(index) = mpi_bonds(index) + proximity(distance,Rb,Db,inv_Db)
          END IF    
       END DO       
    END DO

#ifdef MPI
    ! gather the calculated bond numbers (in "mpi_bonds") and broadcast them to all cpus (to "bonds")
    CALL MPI_ALLGATHERV(mpi_bonds,bondcount(cpu_id+1),MPI_DOUBLE_PRECISION,&
         bonds,bondcount,bondpos,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,mpistat)
#else
    bonds = mpi_bonds
#endif

    RETURN
  END SUBROUTINE update_bond_numbers_p

  ! Calculates the kinetic energy and temperature of the system.
  ! (These are directly proportional to each other.)
  ! *mbs molecules
  ! *ats atoms
  ! *ene kinetic energy
  ! *temper temperature
  ! *n_free number of degrees of freedom
  ! *ene_lin the linear part of kinetic energy
  ! *ene_rot the rotational part of kinetic energy
  SUBROUTINE kin_energy(mbs,ats,gos,n_free,ene,temper,ene_lin,ene_rot)       
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    REAL(KIND=dp), INTENT(OUT) :: ene, temper, ene_rot, ene_lin
    REAL(KIND=dp) :: mpi_ene, mpi_ene_rot, mpi_ene_lin
    REAL(KIND=dp) :: v(3), w(3)
    INTEGER, INTENT(IN) :: n_free
    INTEGER :: n_mbs, n_ats, n_gos, ii
    
    n_mbs = mbs_size(mbs)
    n_ats = ats_size(ats)
    n_gos = gos_size(gos)
    ene = 0.d0
    ene_lin = 0.d0
    ene_rot = 0.d0
    mpi_ene = 0.d0
    mpi_ene_lin = 0.d0
    mpi_ene_rot = 0.d0
    v = 0.d0
    w = 0.d0

    IF(first_mpi_mb_split > 0)THEN ! split the task equally
       DO ii = first_mpi_mb_split, last_mpi_mb_split ! mbs
          v = mbs(ii)%vel
          w = mbs(ii)%angvel
          ! Ekin = 1/2 m v^2
          mpi_ene_lin = mpi_ene_lin + mbs(ii)%m_tot*(v.o.v)   ! linear movement 
          mpi_ene_rot = mpi_ene_rot + mbs(ii)%m_inert*(w.o.w) ! rotational movement
       END DO
    END IF
    IF(first_mpi_atom_split > 0)THEN
       DO ii = first_mpi_atom_split, last_mpi_atom_split ! atoms
          v = ats(ii)%vel
          mpi_ene_lin = mpi_ene_lin + ats(ii)%mass*(v.o.v)
       END DO
    END IF
    IF(first_mpi_go_split > 0)THEN
       DO ii = first_mpi_go_split, last_mpi_go_split ! gos
          v = gos(ii)%vel
          mpi_ene_lin = mpi_ene_lin + gos(ii)%mass*(v.o.v)
       END DO
    END IF

    mpi_ene_lin = 0.5d0*mpi_ene_lin
    mpi_ene_rot = 0.5d0*mpi_ene_rot

#ifdef MPI
    ! sum the energy components from all cpus
    CALL MPI_ALLREDUCE(mpi_ene_lin,ene_lin,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_ALLREDUCE(mpi_ene_rot,ene_rot,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
         MPI_COMM_WORLD,mpistat)
#else
    ene_lin = mpi_ene_lin
    ene_rot = mpi_ene_rot
#endif

    ! the total kinetic energy
    ene = ene_lin + ene_rot

    ! temperature = 2*E_kin / N kB, N = 6*N_mb + 3*N_atom (degrees of freedom)
    temper = 2.d0*ene / (REAL(n_free,KIND=dp)*kB)    

    RETURN
  END SUBROUTINE kin_energy

  ! Calculates the potential energy of the system.
  ! The contributions to the potential energy from
  ! MB-MB Lennard-Jones interaction, MB-MB hydrogen bonds,
  ! MB-atom potential, atom-atom potential and constraints
  ! are also given separately.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *btype index for the types of boundaries
  ! *bval parametes for the boundary conditions
  ! *mm_nbors MB-MB neighbors list
  ! *ma_nbors MB-atom neighbors list
  ! *aa_nbors atom-atom neighbors list
  ! *mm_n_ns MB-MB numbers of neighbors
  ! *ma_n_ns MB-atom numbers of neighbors
  ! *aa_n_ns atom-atom numbers of neighbors
  ! *params physical parameters
  ! *pdone logical list of Mb-Mb pairs, will be used for preventing double counting
  ! *apdone logical list for atom-atom pairs
  ! *is_constr true if there are constraints in place
  ! *ene total potential energy
  ! *ene_lj MB-MB Lennard-Jones contribution
  ! *ene_hb MB-MB hydrogen bond contribution
  ! *ene_ma MB-atom potential contribution
  ! *ene_aa atom-atom potential contribution
  ! *ene_con constraint potential contribution
  SUBROUTINE pot_energy(mbs,ats,gos,&
       mm_nbors,mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
       bonds,cell,pbc,btype,bval,params,go,is_constr,&
       ene,ene_lj,ene_hb,ene_ma,ene_aa,ene_con,&
       ene_mbgo, ene_go, ene_stretch, ene_bend, ene_torsion, ene_native, ene_nonnat)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3)
    INTEGER, INTENT(IN) :: btype(3)
    LOGICAL, INTENT(IN) :: pbc(3), is_constr
    REAL(KIND=dp), INTENT(OUT) :: ene, ene_lj, ene_hb, ene_ma, ene_aa, ene_con,&
         ene_mbgo, ene_go, ene_stretch, ene_bend, ene_torsion, ene_native, ene_nonnat
    REAL(KIND=dp) :: mpi_ene, mpi_ene_lj, mpi_ene_hb, mpi_ene_ma, mpi_ene_aa, mpi_ene_con,&
         mpi_ene_mbgo, mpi_ene_go, mpi_ene_stretch, mpi_ene_bend, &
         mpi_ene_torsion, mpi_ene_native, mpi_ene_nonnat
    REAL(KIND=dp) :: denom, phi, rij, inv_rij, hb_ij, hb_i, boi, boj, &
         armsi(3,4), armsj(3,4), hij(3), ratio, &
         vij(3), u_dot_r(2,4), gauss_rij, gauss_th(2,4), cos_phi(4,4), displ(3), &
         vij2(3), rij2, inv_rij2, vij3(3), rij3, inv_rij3, dot, cosine, &
         hij2(3), hij3(3), projection1(3), projection2(3)
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_ns(:), &
         ma_nbors(:,:), ma_n_ns(:), aa_nbors(:,:), aa_n_ns(:), &
         mg_nbors(:,:), mg_n_ns(:)
    INTEGER :: n_mb, n_at, n_go, ii, jj, nj, kk, ll, mm, nn, ns, type, type2
    LOGICAL :: visited, got_displ

    ene = 0.d0
    ene_lj = 0.d0
    ene_hb = 0.d0
    ene_ma = 0.d0
    ene_aa = 0.d0
    ene_mbgo = 0.d0
    ene_go = 0.d0
    ene_stretch = 0.d0
    ene_bend = 0.d0
    ene_torsion = 0.d0
    ene_native = 0.d0
    ene_nonnat = 0.d0
    mpi_ene = 0.d0
    mpi_ene_lj = 0.d0
    mpi_ene_hb = 0.d0
    mpi_ene_ma = 0.d0
    mpi_ene_aa = 0.d0
    mpi_ene_mbgo = 0.d0
    mpi_ene_go = 0.d0
    mpi_ene_stretch = 0.d0
    mpi_ene_bend = 0.d0
    mpi_ene_torsion = 0.d0
    mpi_ene_native = 0.d0
    mpi_ene_nonnat = 0.d0

    n_mb = SIZE(mbs)
    n_at = SIZE(ats)
    n_go = SIZE(gos)

    !!!!!!!!!!!!!!!!!!
    ! MB - MB energy !
    !!!!!!!!!!!!!!!!!!

    IF(first_mpi_mb > 0)THEN
       DO ii = first_mpi_mb, last_mpi_mb ! loop over the particles dedicated for the cpu
          ns = mm_n_ns(ii) ! number of neighbors
          
          CALL get_arms(mbs(ii),armsi) ! shorthand for the hb unit arms (dangling bond vectors)
          boi = bond_order(bonds(ii),params%v_b) ! shorthand for the bond order factor b_i (b(ii))
          hb_i = 0.d0
          
          DO jj = 1,ns ! loop over neighbors of ii
             nj = mm_nbors(jj,ii) ! nj is the index of the jj:th neighbor of ii
             
             IF(ii < nj)THEN ! to prevent double counting, only check ii < nj
                
                boj = bond_order(bonds(nj),params%v_b) ! shorthand for the bond order factor b_j (b(nj))
                
                ! vector and distance between ii ans nj
                vij(1:3) = vector(mbs(ii),mbs(nj),cell,pbc) ! vector from ii to nj
                rij = sqrt(vij.o.vij) ! distance from ii to nj
                inv_rij = 1.d0/rij    ! inverse distance from ii to nj
                vij = vij*inv_rij     ! unit vector from ii to nj
                
                ! Lennard-Jones energy
                IF(rij < params%cut_lj)THEN ! cutoff
                   mpi_ene_lj = mpi_ene_lj + &
                        ( (params%s_lj*inv_rij)**12 - (params%s_lj*inv_rij)**6 )
                END IF
                
                ! hydrogen bond
                IF(rij < params%cut_hb)THEN ! cutoff
                   
                   CALL get_arms(mbs(nj),armsj)
                   ! calculate dot products between the dangling bonds and unit vector vij
                   DO kk = 1, 4
                      u_dot_r(1,kk) = (armsi(1:3,kk).o.vij) ! u_i^k . r_ij
                      u_dot_r(2,kk) = (armsj(1:3,kk).o.vij) ! u_j^l . r_ij
                   END DO
                   ! calculate gaussians
                   !gauss_rij = exp( -0.5d0* (params%inv_sr*(rij-params%R_hb))**2 )             
                   gauss_rij = gauss(params%inv_sr*(rij-params%R_hb))  ! gaussian of distance
                   DO kk = 1, 4
                      gauss_th(1,kk) = gauss( params%inv_sth*(u_dot_r(1,kk)-1.d0) ) ! gaussian of angle for u_i^k
                      gauss_th(2,kk) = gauss( params%inv_sth*(u_dot_r(2,kk)+1.d0) ) ! gaussian of angle for u_j^l
                      !gauss_th(1,kk) = exp( -0.5d0* (params%inv_sth*(u_dot_r(1,kk)-1.d0))**2 )
                      !gauss_th(2,kk) = exp( -0.5d0* (params%inv_sth*(u_dot_r(2,kk)+1.d0))**2 )             
                   END DO
                   ! calculate cos phi's for all (m,n)
                   DO mm = 1, 4
                      DO nn = 1, 4
                         cos_phi(nn,mm) = ( (armsi(1:3,mm).o.armsj(1:3,nn)) - u_dot_r(1,mm)*u_dot_r(2,nn) ) ! eta_mn
                         IF(ABS(cos_phi(nn,mm)) > norm_tolerance)THEN ! check for being zero
                            denom = sqrt( ( 1.d0 - u_dot_r(1,mm)**2 )*( 1.d0 - u_dot_r(2,nn)**2 ) ) ! xi_mn
                            IF(ABS(denom) > denom_tolerance)THEN ! check for being zero
                               cos_phi(nn,mm) = cos_phi(nn,mm) / denom  ! cos phi_mn = eta_mn/xi_mn
                               IF(cos_phi(nn,mm) > 1.d0) cos_phi(nn,mm) = 1.d0 ! correct for numeric error
                               IF(cos_phi(nn,mm) < -1.d0) cos_phi(nn,mm) = -1.d0
                            ELSE ! xi = 0, so should be eta. this shouldn't happen
                               cos_phi(nn,mm) = 0.d0
                            END IF
                         ELSE ! eta = 0, so cos phi = 0 automatically
                            cos_phi(nn,mm) = 0.d0
                         END IF
                      END DO
                   END DO
                   hb_ij = 0.d0
                   
                   ! sum over arms (k,l), i.e. calculate all terms in the h-bond potential
                   DO kk = 1, 4 ! arms of i, u_i^k
                      DO ll = 1, 4 ! arms of j, u_j^l
                         
                         ! calculate sum (1+cos 3phi)
                         phi = 9.d0
                         ! sum over other arms u_i^m, u_j^n (m /= k, n /= l)
                         DO mm = 1, 4
                            IF(mm /= kk)THEN
                               DO nn = 1, 4
                                  IF(nn /= ll)THEN    
                                     ! triple cos: cos 3x = 4cos**3 x - 3cos x
                                     phi = phi + 4*cos_phi(nn,mm)**3 - 3*cos_phi(nn,mm)
                                  END IF
                               END DO
                            END IF
                         END DO
                         phi = phi*params%e_phi*0.5d0 + 1.d0
                         
                         ! multiply by gaussians and sum
                         hb_ij = hb_ij + gauss_th(1,kk)*gauss_th(2,ll)*phi  
                      END DO
                   END DO
                   
                   IF(hb_ij < 1.0E-12)THEN
                      hb_ij = 0.d0
                   ELSE
                      ! multiply by bond order terms and the distance gaussian
                      hb_i = hb_i + hb_ij*gauss_rij*(boi+boj)*0.50                   
                   END IF
                END IF
                
             END IF
          
          END DO

          mpi_ene_hb = mpi_ene_hb + hb_i
          
       END DO
       ! scale the energies according to given parameters
       mpi_ene_lj = 4.d0*params%e_lj*mpi_ene_lj
       mpi_ene_hb = params%e_hb*mpi_ene_hb
       mpi_ene = mpi_ene_lj + mpi_ene_hb
    END IF

    IF(n_go > 0)THEN

       !!!!!!!!!!!!!!!!!!
       ! Mb - go energy !
       !!!!!!!!!!!!!!!!!!

       IF(first_mpi_mb > 0)THEN
          DO ii = first_mpi_mb, last_mpi_mb ! mbs of this cpu
             ns = mg_n_ns(ii) ! number of neighbors
             DO jj = 1, ns ! neighbors of ii
                nj = mg_nbors(jj,ii) ! jj:th neighbor
                type = nj

                ! no double counting
                
                ! vector and distance between ii ans nj
                vij(1:3) = vector(mbs(ii),gos(nj),cell,pbc)
                rij = sqrt(vij.o.vij)
                !inv_rij = 1.d0/rij
                !vij = vij*inv_rij
                   
                IF(rij < params%cut_lj)THEN ! cutoff

                   ratio = go%sigma / rij
                   ratio = ratio*ratio*ratio*ratio*ratio*ratio
                   ! add the potential energy term for this pair
                   mpi_ene_mbgo = mpi_ene_mbgo + &
                        4.d0*go%epsilon_water(type)*( ratio*ratio - ratio )
                      
                END IF
                                   
             END DO

          END DO
       END IF
              
       !!!!!!!!!!!!!!!!!!
       ! go - go energy !
       !!!!!!!!!!!!!!!!!!

       
       IF(first_mpi_go_split > 0)THEN
          DO ii = first_mpi_go_split, last_mpi_go_split ! loop over gos of this cpu
          
             IF(ii < go%length)THEN
                ! vector and distance between ii and ii+1
                vij(1:3) = vector(gos(ii),gos(ii+1),cell,pbc) ! vector from ii to ii+1
                hij = vij
                rij = sqrt(vij.o.vij) ! distance from ii to ii+1
                inv_rij = 1.d0/rij    ! inverse distance from ii to ii+1
                vij = vij*inv_rij     ! unit vector from ii to ii+1
                
                ! stretching energy
                mpi_ene_stretch = mpi_ene_stretch + &
                     go%kr * ( rij-go%r_chain(ii) )**2
                
                IF(ii < go%length-1)THEN
                   ! vector and distance between ii+1 and ii+2
                   vij2(1:3) = vector(gos(ii+1),gos(ii+2),cell,pbc) ! vector from ii+1 to ii+2
                   hij2 = vij2
                   rij2 = sqrt(vij2.o.vij2) ! distance from ii+1 to ii+2
                   inv_rij2 = 1.d0/rij2    ! inverse distance from ii+1 to ii+2
                   vij2 = vij2*inv_rij2     ! unit vector from ii+1 to ii+2

                   ! bending energy
                   cosine = -(vij .o. vij2)
                   mpi_ene_bend = mpi_ene_bend + go%ktheta * (acos(cosine) - go%theta_chain(ii))**2
                   
                   IF(ii < go%length-2)THEN
                      ! vector and distance between ii+2 and ii+3
                      vij3(1:3) = vector(gos(ii+2),gos(ii+3),cell,pbc) ! vector from ii+2 to ii+3
                      hij3 = vij3
                      rij3 = sqrt(vij3.o.vij3) ! distance from ii+2 to ii+3
                      inv_rij3 = 1.d0/rij3    ! inverse distance from ii+2 to ii+3
                      vij3 = vij3*inv_rij3     ! unit vector from ii+2 to ii+3

                      ! torsion energy                      
                      projection1 = -(hij - (hij .o. hij2) * inv_rij2 * inv_rij2 * hij2)
                      projection2 = (hij3 - (hij3 .o. hij2) * inv_rij2 * inv_rij2 * hij2)
                      
                      dot = projection1 .o. projection2 
                      ratio = 1.d0 / ( (.norm.projection1) * (.norm.projection2) )
                      
                      ! cos phi = (p_21 . p_34) / ( |p_21| |p_34| ) = dot*ratio
                      cosine = dot*ratio
                      phi = acos(cosine)
                      mpi_ene_torsion = mpi_ene_torsion + &
                           go%kphi1 * cos( phi - go%phi_chain(ii) ) + &
                           go%kphi2 * cos( 3.d0*( phi - go%phi_chain(ii) ))

                   END IF

                END IF

             END IF

             ! native and non-native bonds / repulsion
             DO jj = ii+4, go%length
                vij(1:3) = vector(gos(ii),gos(jj),cell,pbc) ! vector from ii+2 to ii+3
                rij = sqrt(vij.o.vij) ! distance from ii+2 to ii+3
                inv_rij = 1.d0/rij    ! inverse distance from ii+2 to ii+3
                !vij = vij*inv_rij     ! unit vector from ii+2 to ii+3
                IF(go%r_native(ii,jj) > 0.d0)THEN
                   ratio = go%r_native(ii,jj) * inv_rij
                   mpi_ene_native = mpi_ene_native + &
                        go%epsilon * ( 5.d0 * ratio**12 - 6.d0 * ratio**10 )
                ELSE
                   ratio = go%r_rep * inv_rij
                   ratio = ratio*ratio*ratio*ratio*ratio*ratio
                   mpi_ene_nonnat = mpi_ene_nonnat + &
                        go%epsilon * ( ratio*ratio )                   
                END IF
             END DO
                
          END DO

       END IF

       mpi_ene_go = mpi_ene_stretch + mpi_ene_bend + mpi_ene_torsion + &
            mpi_ene_native + mpi_ene_nonnat

       mpi_ene = mpi_ene + mpi_ene_go + mpi_ene_mbgo

    END IF

    IF(n_at > 0)THEN

       !!!!!!!!!!!!!!!!!!!!
       ! Mb - atom energy !
       !!!!!!!!!!!!!!!!!!!!

       IF(first_mpi_mb > 0)THEN
          DO ii = first_mpi_mb, last_mpi_mb ! mbs of this cpu
             ns = ma_n_ns(ii) ! number of neighbors
             DO jj = 1, ns ! neighbors of ii
                nj = ma_nbors(jj,ii) ! jj:th neighbor
                type = ats(nj)%type
                
                ! no double counting
                
                visited = .false.
                ! do something only if MB interacts with type
                IF(params%n_pots(type,0) > 0)THEN 
                   
                   visited = .true.
                   ! vector and distance between ii ans nj
                   vij(1:3) = vector(mbs(ii),ats(nj),cell,pbc)
                   rij = sqrt(vij.o.vij)
                   !inv_rij = 1.d0/rij
                   !vij = vij*inv_rij
                   
                   IF(rij < params%pot_cut(type,0))THEN ! cutoff

                      ! add the potential energy term for this pair
                      mpi_ene_ma = mpi_ene_ma + atomic_pair_potential(rij,params%n_pots(type,0),&
                           params%pot_types(:,type,0),params%pot_params(:,:,type,0))
                      
                   END IF
                   
                END IF
                
                ! do something only if the "H" interact with type
                IF(params%n_pots(type,-1) > 0)THEN
                   
                   !armsi = mbs(ii)%arms
                   CALL get_arms(mbs(ii),armsi)
                   IF(.NOT.visited)THEN ! if vij was already calculated, no need to do it again
                      vij(1:3) = vector(mbs(ii),ats(nj),cell,pbc)
                      !rij = sqrt(vij.o.vij)                
                   END IF
                   
                   ! sum over the arms
                   DO kk = 1, 4
                      ! vector from the tip of the arm to the atom
                      hij = vij - params%l_oh*armsi(1:3,kk)
                      rij = sqrt(hij.o.hij)
                      
                      IF(rij < params%pot_cut(-1,type))THEN ! cutoff
                         
                         ! add the potential energy
                         mpi_ene_ma = mpi_ene_ma + atomic_pair_potential(rij,params%n_pots(type,-1),&
                              params%pot_types(:,-1,type),params%pot_params(:,:,type,-1))
                         
                      END IF
                      
                   END DO
                   
                END IF
                
             END DO

          END DO
       END IF

       !!!!!!!!!!!!!!!!!!!!!!
       ! atom - atom energy !
       !!!!!!!!!!!!!!!!!!!!!!
       IF(first_mpi_atom > 0)THEN
          DO ii = first_mpi_atom, last_mpi_atom ! atoms of this cpu
             type = ats(ii)%type
             ns = aa_n_ns(ii) ! number of neighbors
             DO jj = 1, ns ! loop over neighbors
                nj = aa_nbors(jj,ii) ! jj:th neighbor
                
                IF(ii < nj)THEN ! no double counting
                   type2 = ats(nj)%type ! atomic type
                   
                   ! do something only if type interacts with type2
                   IF(params%n_pots(type2,type) > 0)THEN 
                      
                      ! vector and distance between ii ans nj
                      vij(1:3) = vector(ats(ii),ats(nj),cell,pbc)
                      rij = sqrt(vij.o.vij)
                      !inv_rij = 1.d0/rij
                      !vij = vij*inv_rij
                      
                      IF(rij < params%pot_cut(type2,type))THEN ! cutoff

                         ! add the potential energy
                         mpi_ene_aa = mpi_ene_aa + atomic_pair_potential(rij,params%n_pots(type2,type),&
                              params%pot_types(:,type2,type),params%pot_params(:,:,type2,type))
                         
                      END IF
                      
                   END IF
                   
                END IF
             END DO
             
          END DO
             
       END IF
       mpi_ene = mpi_ene + mpi_ene_ma + mpi_ene_aa

    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Energy due to boundaries !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mpi_ene_con = 0.d0
    ene_con = 0.d0
    IF(first_mpi_mb_split > 0)THEN
       DO jj = 1, 3 ! loop over x,y,z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall
             DO ii = first_mpi_mb_split, last_mpi_mb_split ! loop over mbs of this cpu
                ! add the potential due to the wall
                IF(mbs(ii)%pos(jj) < 0.d0)THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*mbs(ii)%pos(jj)**2
                ELSE IF(mbs(ii)%pos(jj) > cell(jj))THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*(mbs(ii)%pos(jj)-cell(jj))**2
                END IF
             END DO
          END IF
       END DO
    END IF
    IF(first_mpi_atom_split > 0)THEN    
       DO jj = 1, 3 ! loop over x,y,z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall
             DO ii = first_mpi_atom_split, last_mpi_atom_split ! loop over atoms of this cpu
                ! add the potential doe to the wall
                IF(ats(ii)%pos(jj) < 0.d0)THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*ats(ii)%pos(jj)**2
                ELSE IF(ats(ii)%pos(jj) > cell(jj))THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*(ats(ii)%pos(jj)-cell(jj))**2
                END IF
             END DO
          END IF
       END DO
    END IF
    IF(first_mpi_go_split > 0)THEN    
       DO jj = 1, 3 ! loop over x,y,z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall
             DO ii = first_mpi_go_split, last_mpi_go_split ! loop over atoms of this cpu
                ! add the potential doe to the wall
                IF(gos(ii)%pos(jj) < 0.d0)THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*gos(ii)%pos(jj)**2
                ELSE IF(gos(ii)%pos(jj) > cell(jj))THEN
                   mpi_ene_con = mpi_ene_con + 0.5d0*bval(jj)*(gos(ii)%pos(jj)-cell(jj))**2
                END IF
             END DO
          END IF
       END DO
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Energy due to constraints !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(is_constr)THEN
       IF(first_mpi_mb_split > 0)THEN
          DO ii = first_mpi_mb_split, last_mpi_mb_split ! loop over mbs of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x,y,z
                SELECT CASE(mbs(ii)%constrained(jj)) ! type of constraint
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN ! calculate the displacement if needed
                      displ = displacement(mbs(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   ! add the harmonic potential energy
                   mpi_ene_con = mpi_ene_con+0.5d0*mbs(ii)%well(jj)*(displ(jj))**2
                CASE(ext_force_index) ! external force
                   IF(.NOT.got_displ)THEN
                      displ = displacement(mbs(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   ! add the potential energy of the constant force
                   mpi_ene_con = mpi_ene_con-mbs(ii)%well(jj)*(displ(jj))                
                END SELECT
             END DO
          END DO
       END IF

       IF(first_mpi_atom_split > 0)THEN
          DO ii = first_mpi_atom_split, last_mpi_atom_split ! loop over atoms of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x,y,z
                SELECT CASE(ats(ii)%constrained(jj)) ! type of constraint
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN
                      displ = displacement(ats(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_ene_con = mpi_ene_con+0.5d0*ats(ii)%well(jj)*(displ(jj))**2
                CASE(ext_force_index) ! external force
                   IF(.NOT.got_displ)THEN
                      displ = displacement(ats(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_ene_con = mpi_ene_con-ats(ii)%well(jj)*(displ(jj))                
                END SELECT
             END DO
          END DO
       END IF
       mpi_ene = mpi_ene + mpi_ene_con

       IF(first_mpi_go_split > 0)THEN
          DO ii = first_mpi_go_split, last_mpi_go_split ! loop over atoms of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x,y,z
                SELECT CASE(gos(ii)%constrained(jj)) ! type of constraint
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN
                      displ = displacement(gos(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_ene_con = mpi_ene_con+0.5d0*gos(ii)%well(jj)*(displ(jj))**2
                CASE(ext_force_index) ! external force
                   IF(.NOT.got_displ)THEN
                      displ = displacement(gos(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_ene_con = mpi_ene_con-gos(ii)%well(jj)*(displ(jj))                
                END SELECT
             END DO
          END DO
       END IF
       mpi_ene = mpi_ene + mpi_ene_con
    END IF

#ifdef MPI
    ! sum energies from all cpus and save them at the master cpu
    CALL MPI_REDUCE(mpi_ene,ene,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_lj,ene_lj,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_hb,ene_hb,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_ma,ene_ma,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_aa,ene_aa,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_con,ene_con,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_mbgo,ene_mbgo,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_go,ene_go,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_stretch,ene_stretch,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_bend,ene_bend,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_torsion,ene_torsion,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_native,ene_native,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)
    CALL MPI_REDUCE(mpi_ene_nonnat,ene_nonnat,1,MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
         MPI_COMM_WORLD,mpistat)


#else
    ene = mpi_ene
    ene_lj = mpi_ene_lj
    ene_hb = mpi_ene_hb
    ene_ma = mpi_ene_ma
    ene_aa = mpi_ene_aa
    ene_con = mpi_ene_con
    ene_mbgo = mpi_ene_mbgo
    ene_go = mpi_ene_go
    ene_stretch = mpi_ene_stretch
    ene_bend = mpi_ene_bend
    ene_torsion = mpi_ene_torsion
    ene_native = mpi_ene_native
    ene_nonnat = mpi_ene_nonnat
#endif

    RETURN
  END SUBROUTINE pot_energy

  ! Gives the molecules random velocities according to
  ! Boltzmann distribution (also angular velocities for MB's).
  ! The final velocities and angular velocities will be 
  ! shifted to give SUM P = 0 and SUM L = 0.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *random_mb logical list, true for molecules that should get a new random velocity
  ! *random_at logical list, true for atoms that should get a new random velocity
  ! *temperature the temperature for velocity generating
  SUBROUTINE assign_random_velocities(mbs,ats,gos,random_mb,random_at,random_go,temperature)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    LOGICAL, POINTER :: random_mb(:,:), random_at(:), random_go(:)
    REAL(KIND=dp), INTENT(IN) :: temperature
    REAL(KIND=dp) :: sigma, &
         v(3), w(3), sv(3), sw(3)
    INTEGER :: ii, n_mbs, n_ats, n_gos, n_rand, n_angrand

    sigma = sqrt(kb*temperature)
    IF(associated(mbs))THEN
       n_mbs = SIZE(mbs(:))
    ELSE
       n_mbs = 0
    END IF
    IF(associated(ats))THEN
       n_ats = SIZE(ats(:))
    ELSE
       n_ats = 0
    END IF
    n_gos = gos_size(gos)
    n_rand = 0
    n_angrand = 0
    sv = 0.d0
    sw = 0.d0   

    DO ii = 1, n_mbs ! loop over mbs
       IF(random_mb(ii,1) .OR. random_mb(ii,2))THEN ! the first test is for random velocity, the second for random angular velocity
          CALL rand_2boltzmann_vectors(v,w) ! 2 random velocities (unnormalized)
          IF(random_mb(ii,1))THEN ! velocity
             mbs(ii)%vel = sigma/sqrt(mbs(ii)%m_tot)*v
             n_rand = n_rand+1
             sv = sv+mbs(ii)%m_tot*mbs(ii)%vel ! sum of momenta, needed for assigning sum p = 0
          END IF
          IF(random_mb(ii,2))THEN ! angular velocity
             mbs(ii)%angvel = sigma/sqrt(mbs(ii)%m_inert)*w
             n_angrand = n_angrand+1
             sw = sw+mbs(ii)%m_inert*mbs(ii)%angvel
          END IF

       END IF
    END DO

    DO ii = 1, n_gos ! loop over gos
       IF(random_go(ii))THEN ! random velocity required
          CALL rand_2boltzmann_vectors(v,w)
          gos(ii)%vel = sigma/sqrt(gos(ii)%mass)*v

          n_rand = n_rand+1
          sv = sv+gos(ii)%mass*gos(ii)%vel
       END IF
    END DO

    DO ii = 1, n_ats ! loop over atoms
       IF(random_at(ii))THEN ! random velocity required
          CALL rand_2boltzmann_vectors(v,w)
          ats(ii)%vel = sigma/sqrt(ats(ii)%mass)*v

          n_rand = n_rand+1
          sv = sv+ats(ii)%mass*ats(ii)%vel
       END IF
    END DO

    ! Average momentum and angular momentum
    IF(n_rand > 0)THEN
       sv = sv/REAL(n_rand,KIND=dp)
    END IF
    IF(n_angrand > 0)THEN
       sw = sw/REAL(n_angrand,KIND=dp)
    END IF

    ! Shift the assigned velocities so that the sum is zero.
    ! Predefined velocities stay as they were.
    DO ii = 1, n_mbs
       IF(random_mb(ii,1))THEN
          mbs(ii)%vel = mbs(ii)%vel - sv/mbs(ii)%m_tot
       END IF
       IF(random_mb(ii,2))THEN
          mbs(ii)%angvel = mbs(ii)%angvel - sw/mbs(ii)%m_inert
       END IF
    END DO
    DO ii = 1, n_gos
       IF(random_go(ii))THEN
          gos(ii)%vel = gos(ii)%vel - sv/gos(ii)%mass
       END IF
    END DO
    DO ii = 1, n_ats
       IF(random_at(ii))THEN
          ats(ii)%vel = ats(ii)%vel - sv/ats(ii)%mass
       END IF
    END DO

  END SUBROUTINE assign_random_velocities


  ! Calculates the forces and torques acting on molecules and the
  ! forces acting on atoms. The program spends most of its computational
  ! time executing this routine.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *cell supercell dimensions
  ! *pbc true if periodic boundaries
  ! *bonds list for the bond count of MB molecules (z_i)
  ! *btype index for the types of boundaries
  ! *bval parametes for the boundary conditions
  ! *mm_nbors MB-MB neighbors list
  ! *ma_nbors MB-atom neighbors list
  ! *aa_nbors atom-atom neighbors list
  ! *mm_n_ns MB-MB numbers of neighbors
  ! *ma_n_ns MB-atom numbers of neighbors
  ! *aa_n_ns atom-atom numbers of neighbors
  ! *params physical parameters
  ! *control control parameters
  ! *pair_done logical list of Mb-Mb pairs, will be used for preventing double counting (obsolete)
  ! *atom_pair_done logical list for atom-atom pairs (obsolete)
  ! *is_constr true if there are constraints in place
  ! *m_force forces acting on molecules
  ! *torque torques acting on molecules
  ! *a_force forces acting on atoms
  ! *apply_thermostat if false, forces due to thermostat are not applied
  SUBROUTINE calc_forces(mbs,ats,gos,mm_nbors,&
       mm_n_ns,ma_nbors,ma_n_ns,aa_nbors,aa_n_ns,mg_nbors,mg_n_ns,&
       bonds,cell,pbc,btype,bval,control,params,go,is_constr,apply_thermostat,&
       m_force,torque,a_force,g_force,virial,n_bond)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), INTENT(IN) :: cell(3), bval(3)
    INTEGER, INTENT(IN) :: btype(3)
    LOGICAL, INTENT(IN) :: pbc(3), is_constr
    LOGICAL, INTENT(IN) :: apply_thermostat
    REAL(KIND=dp), POINTER :: m_force(:,:), torque(:,:), a_force(:,:), g_force(:,:), n_bond(:)
    REAL(KIND=dp), INTENT(OUT) :: virial
    REAL(KIND=dp) :: g_ijkl(4,4), phi_ijkl(4,4), phi, rij, inv_rij, &
         nablaphi(3,4,4), nablaphi_term(3,4,4), sum_nphi(3), dtricos, boi, boj, &
         armsi(3,4), armsj(3,4), eta(4,4), xi(4,4), proji(3,4), projj(3,4), &
         vij(3), hij(3), u_dot_r(2,4), gauss_rij, gauss_th(2,4), cos_phi(4,4), triplecos_phi(4,4), &
         Fpair_ij(3), Tpair_i(3), Tpair_j(3), &
         torquephi(3,4,4), torquephi_term(3,4,4), sum_tphi(3), &
         torquephij(3,4,4), torquephi_termj(3,4,4), sum_tphij(3), &
         Many_coeff, sum_g_times_phi, many_coeff_i, many_coeff_j, displ(3), thevirial, limit1, limit2, &
         ratio, vij2(3), hij2(3), rij2, inv_rij2, vij3(3), hij3(3), rij3, inv_rij3, factor, &
         dot, cosine, tmp1(3), tmp2(3), tmp3(3), &
         Fcomponent1(3), Fcomponent2(3), Fcomponent3(3), Fcomponent4(3), projection1(3), projection2(3), &
         v12(3), v23(3), v34(3), inv_p1, inv_p2, inv_r2, ddot(3,4), dp1(3,4), dp2(3,4)
    REAL(KIND=dp), POINTER :: bonds(:)
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_ns(:), &
         ma_nbors(:,:), ma_n_ns(:), &
         aa_nbors(:,:), aa_n_ns(:), &
         mg_nbors(:,:), mg_n_ns(:)
    INTEGER :: n_mb, n_at, n_go, ii, jj, nj, kk, ll, mm, nn, ns, &
         qq, tt, third, nsj, qj, type, type2, ia
    LOGICAL :: visited, got_displ

    ! set sums to zero
    m_force = 0.d0
    torque = 0.d0
    a_force = 0.d0
    g_force = 0.d0
    small_m_force = 0.d0
    small_torque = 0.d0
    small_a_force = 0.d0
    small_g_force = 0.d0
    tiny_m_force = 0.d0
    tiny_torque = 0.d0
    tiny_a_force = 0.d0
    tiny_g_force = 0.d0
    mpi_m_force = 0.d0
    mpi_torque = 0.d0
    mpi_a_force = 0.d0
    mpi_g_force = 0.d0
    virial = 0.d0
    n_mb = SIZE(mbs)
    n_at = SIZE(ats)
    n_go = SIZE(gos)
    
    limit1 = 1.0E-8
    limit2 = 1.0E-16

    ! a coefficient needed for the tersoff-like many body terms
    Many_coeff = 0.25d0*4**(params%v_b)*params%v_b*params%e_hb*params%inv_db         

    ! if bonding statistics are measured, initialize
    IF(control%bond_writer /= noxyz_index)THEN
       mpi_n_bond = 0.0
    END IF

    ! mpi workload check
    CALL start_timer()

    !!!!!!!!!!!!!!!!!!
    ! MB - MB forces !
    !!!!!!!!!!!!!!!!!!

    IF(first_mpi_mb > 0)THEN ! this cpu works on mbs
       DO ii = first_mpi_mb, last_mpi_mb ! sum over molecules dedicated to this cpu
          ns = mm_n_ns(ii) ! number of neighbors
          
          CALL get_arms(mbs(ii),armsi) ! shorthand for hb arm vectors of ii
          boi = bond_order(bonds(ii),params%v_b) ! bond-order term for ii, b_i
          many_coeff_i = many_coeff*bonds(ii)**(-params%v_b-1) ! coefficient for many-body forces (involving b_i)
          
          DO jj = 1,ns ! sum over neighbors of ii
             nj = mm_nbors(jj,ii) ! jj_th neighbor of ii
             
             IF(ii < nj)THEN ! prevent double counting by only checking pairs ii < nj
                many_coeff_j = many_coeff*bonds(nj)**(-params%v_b-1) ! coefficient for many-body forces (involving b_j)

                ! unit vector and distance between ii and nj
                vij(1:3) = vector(mbs(ii),mbs(nj),cell,pbc) ! vector from ii to nj
                rij = sqrt(vij.o.vij) ! distance from ii to nj
                inv_rij = 1.d0/rij ! inverse distance from ii to nj
                vij = vij*inv_rij ! unit vector from ii to nj
                
                Fpair_ij = 0.d0
                Tpair_i = 0.d0
                Tpair_j = 0.d0
                
                ! Lennard-Jones force
                IF(rij < params%cut_lj)THEN ! cutoff
                   Fpair_ij = 24.d0*params%e_lj*inv_rij* &
                        ( (params%s_lj*inv_rij)**6 - 2.d0*(params%s_lj*inv_rij)**12 )*vij(1:3)

                   DO ia = 1, 3
                      IF(ABS(Fpair_ij(ia)) > limit1)THEN
                         ! sum the force
                         mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         mpi_m_force(ia,nj) = mpi_m_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                         ! sum the force
                         small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         small_m_force(ia,nj) = small_m_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                         ! sum the force
                         tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         tiny_m_force(ia,nj) = tiny_m_force(ia,nj) - Fpair_ij(ia)                       
                      END IF
                   END DO

!!$                      ! sum the force
!!$                      mpi_m_force(1:3,ii) = mpi_m_force(1:3,ii) + Fpair_ij
!!$                      ! store also the reaction force
!!$                      mpi_m_force(1:3,nj) = mpi_m_force(1:3,nj) - Fpair_ij

                   ! get the virial
                   virial = virial - (Fpair_ij.o.hij)
                   
                   Fpair_ij = 0.d0
                END IF
                
                ! Pairwise part of the H-bond force
                IF(rij < params%cut_hb)THEN ! cutoff
                   
                   boj = bond_order(bonds(nj),params%v_b) ! bond-order term for nj, b_j
                   CALL get_arms(mbs(nj),armsj) ! shorthand for hb arm vectors of ii
                   
                   ! calculate hb-arm dot products and projections
                   DO kk = 1, 4 ! loop over dangling bond arms
                      u_dot_r(1,kk) = (armsi(1:3,kk).o.vij) ! dot product u_i^k . r_ij
                      u_dot_r(2,kk) = (armsj(1:3,kk).o.vij) ! dot product u_j^l . r_ij
                      proji(1:3,kk) = armsi(1:3,kk) - u_dot_r(1,kk)*vij ! projection ^{r_ij}u_i^k
                      projj(1:3,kk) = armsj(1:3,kk) - u_dot_r(2,kk)*vij ! projection ^{r_ij}u_j^l
                   END DO
                   
                   ! calculate gaussians
                   gauss_rij = gauss(params%inv_sr*(rij-params%R_hb)) ! distance
                   DO kk = 1, 4 ! loop over arms
                      gauss_th(1,kk) = gauss( params%inv_sth*(u_dot_r(1,kk)-1.d0) ) ! angle of u_i^k
                      gauss_th(2,kk) = gauss( params%inv_sth*(u_dot_r(2,kk)+1.d0) ) ! angle of u_j^l
                   END DO
                   
                   ! if bonding stats are recorded, store them now
                   IF(control%bond_writer /= noxyz_index)THEN
                      ! sum the bonds
                      DO kk = 1, 4 ! loop over arms of i
                         DO ll = 1, 4 ! loop over arms of j
                            ! bonding is defined as a sum of the product of gaussians
                            mpi_n_bond(ii) = mpi_n_bond(ii) + gauss_rij*gauss_th(1,kk)*gauss_th(2,ll)
                            mpi_n_bond(nj) = mpi_n_bond(nj) + gauss_rij*gauss_th(1,kk)*gauss_th(2,ll)
                         END DO
                      END DO
                   END IF

                   ! calculate cos phi's (cos phi = eta/xi)
                   DO mm = 1, 4 ! loop over arms of i
                      DO nn = 1, 4 ! loop over arms of j
                         !eta(nn,mm) = ( (armsi(1:3,mm).o.armsj(1:3,nn)) - u_dot_r(1,mm)*u_dot_r(2,nn) )
                         eta(nn,mm) = proji(1:3,mm) .o. projj(1:3,nn)
                      END DO
                   END DO
                   
                   DO mm = 1, 4 ! loop over arms of i
                      DO nn = 1, 4 ! loop over arms of j
                         IF(ABS(eta(nn,mm)) > norm_tolerance)THEN ! if eta = 0, cos phi = 0
                            xi(nn,mm) = sqrt( ( 1.d0 - u_dot_r(1,mm)**2 )*( 1.d0 - u_dot_r(2,nn)**2 ) )
                         ELSE
                            xi(nn,mm) = 0.d0 ! set xi to 0, this will set cos phi to 0 as well
                         END IF
                      END DO
                   END DO
                   
                   WHERE(xi > denom_tolerance) ! if xi = 0, set cos phi to 0, otherwise cos phi = eta/xi
                      cos_phi = eta / xi
                   ELSEWHERE
                      cos_phi = 0.d0
                   END WHERE
                   ! correct possible numeric errors
                   WHERE(cos_phi > 1.d0) cos_phi = 1.d0 
                   WHERE(cos_phi < -1.d0) cos_phi = -1.d0
                   
                   ! triple cos: cos 3x = 4cos**3 x - 3cos x           
                   !DO mm = 1, 4
                   !   DO nn = 1, 4
                   !      triplecos_phi(nn,mm) = 4*cos_phi(nn,mm)**3 - 3*cos_phi(nn,mm)
                   !   END DO
                   !END DO
                   triplecos_phi = 4*cos_phi**3 - 3*cos_phi
                   
                   ! calculate the "derivative of phi" and "rotational derivative of phi" component terms
                   DO mm = 1, 4 ! loop over arms of i
                      DO nn = 1, 4 ! loop over arms of j
                         
                         ! debug: nabla eta = inv_rij*( u_dot_r(2,nn)*proji(mm,1:3) + u_dot_r(1,mm)*projj(nn,1:3) )
                         ! debug: nabla xi = inv_rij/xi(mm,nn)*&
                         !     ( u_dot_r(1,mm)*(1.d0-u_dot_r(2,nn)**2)*proji(mm,1:3) + &
                         !       u_dot_r(2,nn)*(1.d0-u_dot_r(1,mm)**2)*projj(nn,1:3) )                      
                         
                         IF(ABS(xi(nn,mm)) > denom_tolerance)THEN ! cos phi isn't 0
                            dtricos = (12.d0*(cos_phi(nn,mm))**2 - 3.d0) ! D cos 3phi
                            
                            nablaphi_term(1:3,nn,mm) = & ! D Phi_ij^nm
                                 dtricos * inv_rij / xi(nn,mm)**3 * ( &
                                 ( u_dot_r(2,nn)*xi(nn,mm)**2 - u_dot_r(1,mm)*(1.d0-u_dot_r(2,nn)**2)*eta(nn,mm) )*proji(1:3,mm) + &
                                 ( u_dot_r(1,mm)*xi(nn,mm)**2 - u_dot_r(2,nn)*(1.d0-u_dot_r(1,mm)**2)*eta(nn,mm) )*projj(1:3,nn) )
                            
                            ! exact torques - note that the torques of i and j are not symmetric
                            torquephi_term(1:3,nn,mm) = & ! D_phi(i) Phi_ij^nm (derivative with respect to rotation of i)
                                 -dtricos/xi(nn,mm) * &
                              ( armsi(1:3,mm) .x. ( armsj(1:3,nn) + &
                              ( u_dot_r(1,mm)/xi(nn,mm)*( 1.d0-u_dot_r(2,nn)**2 )*cos_phi(nn,mm) - u_dot_r(2,nn) )*vij ) )

                            torquephi_termj(1:3,nn,mm) = & ! D_phi(j) Phi_ij^nm (derivative with respect to rotation of j)
                                 -dtricos/xi(nn,mm) * &
                              ( armsj(1:3,nn) .x. ( armsi(1:3,mm) + &
                              ( u_dot_r(2,nn)/xi(nn,mm)*( 1.d0-u_dot_r(1,mm)**2 )*cos_phi(nn,mm) - u_dot_r(1,mm) )*vij ) )
                            
                         ELSE ! cos phi is zero
                            nablaphi_term(1:3,nn,mm) = 0.d0
                            torquephi_term(1:3,nn,mm) = 0.d0
                            torquephi_termj(1:3,nn,mm) = 0.d0
                         END IF
                      END DO
                   END DO
                   
                   ! calculate G_ij^kl and Phi_ij^kl by summing the needed pre-calculated terms
                   sum_g_times_phi = 0.d0
                   DO kk = 1, 4 ! loop over arms of i
                      DO ll = 1, 4 ! loop over arms of j
                         g_ijkl(ll,kk) = gauss_rij*gauss_th(1,kk)*gauss_th(2,ll) ! G_ij^kl
                         phi = 9.d0
                         sum_nphi = 0.d0
                         sum_tphi = 0.d0
                         sum_tphij = 0.d0
                         DO mm = 1, 4 ! loop over other arms of i (m /= k)
                            IF(mm /= kk)THEN
                               DO nn = 1, 4 ! loop over other arms of j (n /= l)
                                  IF(nn /= ll)THEN 
                                     phi = phi + triplecos_phi(nn,mm) ! sum of cos 3phi
                                     sum_nphi = sum_nphi + nablaphi_term(1:3,nn,mm) ! sum of D Phi
                                     sum_tphi = sum_tphi + torquephi_term(1:3,nn,mm) ! sum of D_phi(i) Phi
                                     sum_tphij = sum_tphij + torquephi_termj(1:3,nn,mm) ! sum of D_phi(j) Phi
                                  END IF
                               END DO
                            END IF
                         END DO
                         ! multiply with the required parameters
                         phi_ijkl(ll,kk) = phi*params%e_phi*0.5d0 + 1.d0
                         nablaphi(1:3,ll,kk) = sum_nphi*params%e_phi*0.5d0
                         torquephi(1:3,ll,kk) = sum_tphi*params%e_phi*0.5d0 
                         torquephij(1:3,ll,kk) = sum_tphij*params%e_phi*0.5d0
                         sum_g_times_phi = sum_g_times_phi + phi_ijkl(ll,kk)*g_ijkl(ll,kk)
                      END DO
                   END DO
                   
                   ! calculate the pairwise force and torque between ii and nj
                   DO kk = 1, 4 ! loop over arms of i
                      DO ll = 1, 4 ! loop over arms of j

                         ! Pairwise force between i and j due to the H-bond potential
                         Fpair_ij = Fpair_ij + g_ijkl(ll,kk)*( &
                              phi_ijkl(ll,kk)*(params%inv_sr**2*(rij-params%R_hb)*vij + &
                              params%inv_sth**2*inv_rij*( (u_dot_r(1,kk)-1)*proji(1:3,kk)+(u_dot_r(2,ll)+1)*projj(1:3,ll) ) ) + &
                              nablaphi(1:3,ll,kk) )

                         ! Torques on i and j due to the H-bond potential
                         Tpair_i = Tpair_i + g_ijkl(ll,kk)*( &
                              phi_ijkl(ll,kk)*params%inv_sth**2*(u_dot_r(1,kk) - 1.d0)*(vij.x.armsi(1:3,kk)) - &
                              torquephi(1:3,ll,kk) )
                         Tpair_j = Tpair_j + g_ijkl(ll,kk)*( &
                              phi_ijkl(ll,kk)*params%inv_sth**2*(u_dot_r(2,ll) + 1.d0)*(vij.x.armsj(1:3,ll)) - &
                              torquephij(1:3,ll,kk) )
                         
                      END DO
                   END DO
                   ! Multiply by energy scale and bond order terms
                   Fpair_ij = -0.5d0*params%e_hb*(boi+boj)*Fpair_ij
                   Tpair_i  = -0.5d0*params%e_hb*(boi+boj)*Tpair_i
                   Tpair_j  = -0.5d0*params%e_hb*(boi+boj)*Tpair_j         
                   

                   DO ia = 1,3
                      IF(ABS(Fpair_ij(ia)) > limit1)THEN
                         ! sum the force
                         mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         mpi_m_force(ia,nj) = mpi_m_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                         ! sum the force
                         small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         small_m_force(ia,nj) = small_m_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                         ! sum the force
                         tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                         ! store also the reaction force
                         tiny_m_force(ia,nj) = tiny_m_force(ia,nj) - Fpair_ij(ia)                     
                      END IF
                      ! torque on i
                      IF(ABS(Tpair_i(ia)) > limit1)THEN
                         mpi_torque(ia,ii) = mpi_torque(ia,ii) + Tpair_i(ia)
                      ELSE IF(ABS(Tpair_i(ia)) > limit2 .AND. ABS(Tpair_i(ia)) <= limit1)THEN
                         small_torque(ia,ii) = small_torque(ia,ii) + Tpair_i(ia)
                      ELSE IF(ABS(Tpair_i(ia)) <= limit2)THEN
                         tiny_torque(ia,ii) = tiny_torque(ia,ii) + Tpair_i(ia)
                      END IF
                      ! torque on j
                      IF(ABS(Tpair_j(ia)) > limit1)THEN
                         mpi_torque(ia,nj) = mpi_torque(ia,nj) + Tpair_j(ia)
                      ELSE IF(ABS(Tpair_j(ia)) > limit2 .AND. ABS(Tpair_j(ia)) <= limit1)THEN
                         small_torque(ia,nj) = small_torque(ia,nj) + Tpair_j(ia)
                      ELSE IF(ABS(Tpair_j(ia)) <= limit2)THEN
                         tiny_torque(ia,nj) = tiny_torque(ia,nj) + Tpair_j(ia)
                      END IF
                   END DO


!!$                      ! sum the force
!!$                      mpi_m_force(1:3,ii) = mpi_m_force(1:3,ii) + Fpair_ij
!!$                      ! store also the reaction force
!!$                      mpi_m_force(1:3,nj) = mpi_m_force(1:3,nj) - Fpair_ij
!!$                      ! torque on i
!!$                      mpi_torque(1:3,ii) = mpi_torque(1:3,ii) + Tpair_i
!!$                      mpi_torque(1:3,nj) = mpi_torque(1:3,nj) + Tpair_j
                   
                   ! get the virial
                   virial = virial - (Fpair_ij.o.hij)
                   
                   ! calculate the many-body terms incorporating the pair
                   ! of molecules ii and nj, i.e. the terms where the sum of g_ijkl*phi_ijkl
                   ! appears - this complicated sum was just calculated
                   IF(bonds(ii) > 4.d0)THEN ! many body force appears only if the neighbor count is more than 4 - first check ii
                      
                      DO tt = 1, ns ! go over neighbors of ii and add the many body contribution
                         third = mm_nbors(tt,ii) ! "third" is the target molecule on which the force acts.
                                                 ! as this is a three-body force, the target may be other than ii or nj
                         
                         ! here, target is not ii, but the tt:th neighbor (it may be nj)
                         ! get the vector and distance from ii to third
                         vij = vector(mbs(third),mbs(ii),cell,pbc)
                         rij = sqrt(vij.o.vij)
                         vij = vij/rij    
                         
                         IF(ABS(rij-params%r_b) < params%d_b)THEN ! force only if the ii-third distance is in the proximity function "shell area"
                            Fpair_ij = many_coeff_i*sum_g_times_phi* & ! this is not a pair force, but the same variable is used for storage
                                 cos(params%inv_db*(rij-params%r_b))*vij  
                            DO ia = 1, 3
                               iF(ABS(Fpair_ij(ia)) > limit1)THEN
                                  mpi_m_force(ia,third) = mpi_m_force(ia,third) + Fpair_ij(ia) ! force acts on "third"
                               ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                                  small_m_force(ia,third) = small_m_force(ia,third) + Fpair_ij(ia)
                               ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                                  tiny_m_force(ia,third) = tiny_m_force(ia,third) + Fpair_ij(ia)
                               END iF
                            END DO

!!$                               mpi_m_force(1:3,third) = mpi_m_force(1:3,third) + Fpair_ij ! force acts on "third"
                            
                            ! get the virial
                            ! the position is measured from the ii,nj center of mass
                            ! to get periodic bounds correctly. the same point of origin is
                            ! used for all many-body forces incorporating the sum of g_ijkl*phi_ijkl
                            vij = cm_vector(mbs(third),mbs(ii),mbs(nj),cell,pbc)
                            virial = virial + (Fpair_ij.o.vij)
                            
                         END IF
                      END DO
                      
                      ! target is ii itself
                      DO qq = 1, ns ! loop over neighbors of ii
                         qj = mm_nbors(qq,ii) ! qq:th neighbor is qj

                         ! vector and distance ii-qj
                         vij = vector(mbs(ii),mbs(qj),cell,pbc)
                         rij = sqrt(vij.o.vij)
                         vij = vij/rij 
                         
                         IF(ABS(rij-params%r_b) < params%d_b)THEN ! test for being in the "shell area"
                            Fpair_ij = many_coeff_i*sum_g_times_phi* &
                                 cos(params%inv_db*(rij-params%r_b))*vij
                            DO ia = 1, 3
                               iF(ABS(Fpair_ij(ia)) > limit1)THEN
                                  mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia) ! force acts on ii
                               ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                                  small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                               ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                                  tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                               END iF
                            END DO

!!$                               mpi_m_force(1:3,ii) = mpi_m_force(1:3,ii) + Fpair_ij ! force acts on ii

                            ! get the virial
                            vij = 0.5d0*vector(mbs(nj),mbs(ii),cell,pbc)
                            virial = virial + (Fpair_ij.o.vij)
                         END IF
                         
                      END DO
                      
                   END IF
                   
                   IF(bonds(nj) > 4)THEN ! contribution also if nj has more than 4 neighbors
                      nsj = mm_n_ns(nj)
                      DO tt = 1, nsj ! go over neighbors of nj and add the many body contribution
                         third = mm_nbors(tt,nj) ! target molecule
                         
                         ! target is not nj, but the tt:th neighbor (it may be ii)
                         ! get the vector third-nj
                         vij = vector(mbs(third),mbs(nj),cell,pbc)
                         rij = sqrt(vij.o.vij)
                         vij = vij/rij    
                         
                         IF(ABS(rij-params%r_b) < params%d_b)THEN ! check for being in the "shell area"
                            Fpair_ij = many_coeff_j*sum_g_times_phi* &
                                 cos(params%inv_db*(rij-params%r_b))*vij
                            DO ia = 1, 3
                               iF(ABS(Fpair_ij(ia)) > limit1)THEN
                                  mpi_m_force(ia,third) = mpi_m_force(ia,third) + Fpair_ij(ia) ! force acts on "third"
                               ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                                  small_m_force(ia,third) = small_m_force(ia,third) + Fpair_ij(ia)
                               ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                                  tiny_m_force(ia,third) = tiny_m_force(ia,third) + Fpair_ij(ia)
                               END iF
                            END DO

!!$                               mpi_m_force(1:3,third) = mpi_m_force(1:3,third) + Fpair_ij ! force acts on third

                            ! get the virial
                            vij = cm_vector(mbs(third),mbs(ii),mbs(nj),cell,pbc)
                            virial = virial + (Fpair_ij.o.vij)
                         END IF
                      END DO
                      
                      ! target is nj itself                      
                      DO qq = 1, nsj ! loop over neighbors
                         qj = mm_nbors(qq,nj)
                         
                         vij = vector(mbs(nj),mbs(qj),cell,pbc)
                         rij = sqrt(vij.o.vij)
                         vij = vij/rij    
                         
                         IF(ABS(rij-params%r_b) < params%d_b)THEN ! check for being in the "shell area"
                            Fpair_ij = many_coeff_j*sum_g_times_phi* &
                                 cos(params%inv_db*(rij-params%r_b))*vij
                            DO ia = 1, 3
                               iF(ABS(Fpair_ij(ia)) > limit1)THEN
                                  mpi_m_force(ia,nj) = mpi_m_force(ia,nj) + Fpair_ij(ia) ! force acts on nj
                               ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                                  small_m_force(ia,nj) = small_m_force(ia,nj) + Fpair_ij(ia)
                               ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                                  tiny_m_force(ia,nj) = tiny_m_force(ia,nj) + Fpair_ij(ia)
                               END iF
                            END DO

!!$                               mpi_m_force(1:3,nj) = mpi_m_force(1:3,nj) + Fpair_ij ! force acts on nj

                            ! get the virial
                            vij = 0.5d0*vector(mbs(ii),mbs(nj),cell,pbc)
                            virial = virial + (Fpair_ij.o.vij)
                         END IF
                         
                      END DO
                      
                   END IF
                   
                END IF
                                
             END IF

          END DO
       END DO

       !!!!!!!!!!!!!!!!!!!!
       ! MB - atom forces !
       !!!!!!!!!!!!!!!!!!!!

       IF(n_at > 0)THEN
          DO ii = first_mpi_mb, last_mpi_mb ! loop over mbs of this cpu
             ns = ma_n_ns(ii)
             DO jj = 1, ns ! loop over neighbors
                nj = ma_nbors(jj,ii)
                Fpair_ij = 0.d0
                type = ats(nj)%type
                
                ! we only loop over MB-atom pairs (not atom-MB), so no double counting occurs
                
                visited = .false.
                IF(params%n_pots(type,0) > 0)THEN  ! MB intercts with type
                   
                   visited = .true.
                   ! vector and distance between ii ans nj
                   vij(1:3) = vector(mbs(ii),ats(nj),cell,pbc)
                   hij = vij
                   rij = sqrt(vij.o.vij)
                   inv_rij = 1.d0/rij
                   vij = vij*inv_rij

                   IF(rij < params%pot_cut(type,0))THEN ! cutoff
                      
                      ! get the pairwise force
                      Fpair_ij = atomic_pair_force(vij,rij,inv_rij,params%n_pots(type,0),&
                           params%pot_types(:,type,0),params%pot_params(:,:,type,0))

                      ! add to the total force
                      !WHERE(Fpair_ij < 1.0E-36) Fpair_ij = 0.d0
                      DO ia = 1, 3
                         IF(ABS(Fpair_ij(ia)) > limit1)THEN
                            mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia)
                            mpi_a_force(ia,nj) = mpi_a_force(ia,nj) - Fpair_ij(ia)
                         ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                            small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                            small_a_force(ia,nj) = small_a_force(ia,nj) - Fpair_ij(ia)
                         ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                            tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                            tiny_a_force(ia,nj) = tiny_a_force(ia,nj) - Fpair_ij(ia)
                         END IF
                      END DO

!!$                         mpi_m_force(1:3,ii) = mpi_m_force(1:3,ii) + Fpair_ij
!!$                         mpi_a_force(1:3,nj) = mpi_a_force(1:3,nj) - Fpair_ij
                      
                      ! get the virial
                      virial = virial - (Fpair_ij.o.hij)
                   END IF
                   
                END IF
                
                ! interactions between "H" and type
                IF(params%n_pots(type,-1) > 0)THEN
                   
                   ! get the dangling bond arm vectors
                   CALL get_arms(mbs(ii),armsi)
                   IF(.NOT.visited)THEN
                      ! vector from ii to nj
                      vij(1:3) = vector(mbs(ii),ats(nj),cell,pbc)
                      !rij = sqrt(vij.o.vij)  
                   ELSE
                      vij = hij
                   END IF
                   
                   ! sum over the arms
                   DO kk = 1, 4
                      ! vector from the arm to nj
                      hij = vij - params%l_oh*armsi(1:3,kk)
                      rij = sqrt(hij.o.hij)
                      inv_rij = 1.d0/rij
                      hij = hij*inv_rij
                      
                      IF(rij < params%pot_cut(type,-1))THEN ! cutoff
                         
                         ! get the pairwise forcw
                         Fpair_ij = atomic_pair_force(hij,rij,inv_rij,params%n_pots(type,-1),&
                              params%pot_types(:,type,-1),params%pot_params(:,:,type,-1))
                         
                         DO ia = 1, 3
                            ! add to the total force
                            IF(ABS(Fpair_ij(ia)) > limit1)THEN
                               mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia)
                               mpi_a_force(ia,nj) = mpi_a_force(ia,nj) - Fpair_ij(ia)
                            ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                               small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                               small_a_force(ia,nj) = small_a_force(ia,nj) - Fpair_ij(ia)
                            ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                               tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                               tiny_a_force(ia,nj) = tiny_a_force(ia,nj) - Fpair_ij(ia)
                            END IF
                         END DO

!!$                            mpi_m_force(1:3,ii) = mpi_m_force(1:3,ii) + Fpair_ij
!!$                            mpi_a_force(1:3,nj) = mpi_a_force(1:3,nj) - Fpair_ij

                         ! since the force acts on the tip of the arm, it excerts a torque on the mb as well
                         mpi_torque(1:3,ii) = mpi_torque(1:3,ii) + params%l_oh*(armsi(1:3,kk).x.Fpair_ij)

                         ! get the virial
                         virial = virial - (Fpair_ij.o.hij)
                         
                      END IF
                      
                   END DO
                   
                END IF
                
             END DO
          END DO          
       END IF




       !!!!!!!!!!!!!!!!!!
       ! MB - go forces !
       !!!!!!!!!!!!!!!!!!

       IF(n_go > 0)THEN
          DO ii = first_mpi_mb, last_mpi_mb ! loop over mbs of this cpu
             ns = mg_n_ns(ii)
             DO jj = 1, ns ! loop over neighbors
                nj = mg_nbors(jj,ii)
                Fpair_ij = 0.d0
                type = nj
                
                ! we only loop over MB-atom pairs (not atom-MB), so no double counting occurs
                
                ! vector and distance between ii ans nj
                vij(1:3) = vector(mbs(ii),gos(nj),cell,pbc)
                hij = vij
                rij = sqrt(vij.o.vij)
                inv_rij = 1.d0/rij
                vij = vij*inv_rij

                IF(rij < params%cut_lj)THEN ! cutoff

                   ratio = go%sigma / rij
                   ratio = ratio*ratio*ratio*ratio*ratio*ratio
                      
                   ! get the pairwise force
                   Fpair_ij = 24.d0*go%epsilon_water(type)*inv_rij* &
                        ( ratio - 2.d0*ratio*ratio )*vij(1:3)

                   ! add to the total force
                   !WHERE(Fpair_ij < 1.0E-36) Fpair_ij = 0.d0
                   DO ia = 1, 3
                      IF(ABS(Fpair_ij(ia)) > limit1)THEN
                         mpi_m_force(ia,ii) = mpi_m_force(ia,ii) + Fpair_ij(ia)
                         mpi_g_force(ia,nj) = mpi_g_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                         small_m_force(ia,ii) = small_m_force(ia,ii) + Fpair_ij(ia)
                         small_g_force(ia,nj) = small_g_force(ia,nj) - Fpair_ij(ia)
                      ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                         tiny_m_force(ia,ii) = tiny_m_force(ia,ii) + Fpair_ij(ia)
                         tiny_g_force(ia,nj) = tiny_g_force(ia,nj) - Fpair_ij(ia)
                      END IF
                   END DO
                      
                   ! get the virial
                   virial = virial - (Fpair_ij.o.hij)
                END IF
                   
                
             END DO
          END DO
       END IF

    END IF


    !!!!!!!!!!!!!!!!!!
    ! go - go forces !
    !!!!!!!!!!!!!!!!!!

    IF(n_go > 0)THEN

       IF(first_mpi_go_split > 0)THEN
          DO ii = first_mpi_go_split, last_mpi_go_split ! loop over gos of this cpu

             IF(ii < go%length)THEN ! skip the last index so that ii+1 stays in bounds
                ! vector and distance between ii ans ii+1
                vij(1:3) = vector(gos(ii),gos(ii+1),cell,pbc)
                hij = vij
                rij = sqrt(vij.o.vij)
                inv_rij = 1.d0/rij
                vij = vij*inv_rij
                
                ! get the pairwise force
                Fpair_ij = go%kr * 2.d0* (rij-go%r_chain(ii))*vij(1:3)

                ! add to the total force
                !WHERE(Fpair_ij < 1.0E-36) Fpair_ij = 0.d0
                DO ia = 1, 3
                   IF(ABS(Fpair_ij(ia)) > limit1)THEN
                      mpi_g_force(ia,ii) = mpi_g_force(ia,ii) + Fpair_ij(ia)
                      mpi_g_force(ia,ii+1) = mpi_g_force(ia,ii+1) - Fpair_ij(ia)
                   ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                      small_g_force(ia,ii) = small_g_force(ia,ii) + Fpair_ij(ia)
                      small_g_force(ia,ii+1) = small_g_force(ia,ii+1) - Fpair_ij(ia)
                   ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                      tiny_g_force(ia,ii) = tiny_g_force(ia,ii) + Fpair_ij(ia)
                      tiny_g_force(ia,ii+1) = tiny_g_force(ia,ii+1) - Fpair_ij(ia)
                   END IF
                END DO
                
                ! get the virial
                virial = virial - (Fpair_ij.o.hij)
                
                IF(ii < go%length-1)THEN
                   ! vector and distance between ii+1 ans ii+2
                   vij2(1:3) = vector(gos(ii+1),gos(ii+2),cell,pbc)
                   hij2 = vij2
                   rij2 = sqrt(vij2.o.vij2)
                   inv_rij2 = 1.d0/rij2
                   vij2 = vij2*inv_rij2
                   
                   dot = -(hij .o. hij2)
                   cosine = dot * inv_rij * inv_rij2

                   tmp1 = hij2 + dot * inv_rij * inv_rij * hij
                   tmp2 = -hij - dot * inv_rij2 * inv_rij2 * hij2

                   ! nabla cos theta:
                   ! Fcomponent1 = -tmp1 * inv_rij * inv_rij2
                   ! Fcomponent2 = -tmp2 * inv_rij * inv_rij2

                   ! nabla theta = - nabla cos theta / sqrt(1 - cos^2 theta)
                   ratio = 1.d0 / sqrt(1.d0 - cosine*cosine)
                   !Fcomponent1 = tmp1 * inv_rij * inv_rij2 * ratio
                   !Fcomponent2 = tmp2 * inv_rij * inv_rij2 * ratio
                   
                   ! nabla K (theta - theta0)^2 = 2K (theta-theta0) * nabla theta
                   factor = 2.d0*go%ktheta * (acos(cosine) - go%theta_chain(ii))
                   Fcomponent1 = tmp1 * inv_rij * inv_rij2 * ratio * factor
                   Fcomponent2 = tmp2 * inv_rij * inv_rij2 * ratio * factor

                   DO ia = 1, 3
                      IF(ABS(Fcomponent1(ia)) > limit1)THEN
                         mpi_g_force(ia,ii) = mpi_g_force(ia,ii) + Fcomponent1(ia)
                         mpi_g_force(ia,ii+1) = mpi_g_force(ia,ii+1) - Fcomponent1(ia)
                      ELSE IF(ABS(Fcomponent1(ia)) > limit2 .AND. ABS(Fcomponent1(ia)) <= limit1)THEN
                         small_g_force(ia,ii) = small_g_force(ia,ii) + Fcomponent1(ia)
                         small_g_force(ia,ii+1) = small_g_force(ia,ii+1) - Fcomponent1(ia)
                      ELSE IF(ABS(Fcomponent1(ia)) <= limit2)THEN
                         tiny_g_force(ia,ii) = tiny_g_force(ia,ii) + Fcomponent1(ia)
                         tiny_g_force(ia,ii+1) = tiny_g_force(ia,ii+1) - Fcomponent1(ia)
                      END IF
                   END DO
                   DO ia = 1, 3
                      IF(ABS(Fcomponent2(ia)) > limit1)THEN
                         mpi_g_force(ia,ii+2) = mpi_g_force(ia,ii+2) + Fcomponent2(ia)
                         mpi_g_force(ia,ii+1) = mpi_g_force(ia,ii+1) - Fcomponent2(ia)
                      ELSE IF(ABS(Fcomponent2(ia)) > limit2 .AND. ABS(Fcomponent2(ia)) <= limit1)THEN
                         small_g_force(ia,ii+2) = small_g_force(ia,ii+2) + Fcomponent2(ia)
                         small_g_force(ia,ii+1) = small_g_force(ia,ii+1) - Fcomponent2(ia)
                      ELSE IF(ABS(Fcomponent2(ia)) <= limit2)THEN
                         tiny_g_force(ia,ii+2) = tiny_g_force(ia,ii+2) + Fcomponent2(ia)
                         tiny_g_force(ia,ii+1) = tiny_g_force(ia,ii+1) - Fcomponent2(ia)
                      END IF
                   END DO

                   ! get the virial
                   virial = virial - (Fcomponent1.o.hij) + (Fcomponent2.o.hij2)

                   IF(ii < go%length-2)THEN

                      ! vector and distance between ii+2 ans ii+3
                      vij3(1:3) = vector(gos(ii+2),gos(ii+3),cell,pbc)
                      hij3 = vij3
                      rij3 = sqrt(vij3.o.vij3)
                      inv_rij3 = 1.d0/rij3
                      vij3 = vij2*inv_rij3

                      ! projections
                      projection1 = -(hij - (hij .o. hij2) * inv_rij2 * inv_rij2 * hij2)
                      projection2 = (hij3 - (hij3 .o. hij2) * inv_rij2 * inv_rij2 *hij2)
                            
                      inv_r2 = inv_rij2 * inv_rij2
                      v12 = hij
                      v23 = hij2
                      v34 = hij3

                      inv_p1 = 1.d0 / (.norm.projection1)
                      inv_p2 = 1.d0 / (.norm.projection2)
                      ratio = inv_p1 * inv_p2
                      dot = projection1 .o. projection2 
                      
                      ! cos phi = (p_21 . p_34) / ( |p_21| |p_34| ) = dot*ratio
                      cosine = dot*ratio

                      ! D p.p' / (|p||p'|) =
                      ! D (p.p') / (|p||p'|) +
                      ! p.p' / (|p||p'|^2) D |p'| +
                      ! p.p' / (|p|^2|p'|) D |p|
                      
                      ! D (p.p') / (|p||p'|) :       
                      
                      ddot = 0.d0
                      ddot(1:3,1) =  -( v34 - (v34 .o. v23) * inv_r2 * v23 )
                      ddot(1:3,2) =  ( v34 - &
                           (v34 .o. v23) * inv_r2 * (v23 - v12) + &
                           (v12 .o. v23) * inv_r2 * v34 - &
                           2.d0 * (v12 .o. v23) * (v34 .o. v23) * inv_r2 * inv_r2 * v23 )
                      ddot(1:3,3) =  -( v12 + &
                           (v34 .o. v23) * inv_r2 * v12 + &
                           (v12 .o. v23) * inv_r2 * (v34 - v23) - &
                           2.d0 * (v12 .o. v23) * (v34 .o. v23) * inv_r2 * inv_r2 * v23 )
                      ddot(1:3,4) =  ( v12 - (v12 .o. v23) *inv_r2 * v23 )
                      ddot = ddot
                      
                      ! p.p' / (|p||p'|^2) D |p'| = p.p' / (2|p||p'|^3) D (p'.p')
                      dp2 = 0.d0
                      dp2(1:3,2) = (v34 .o. v23) * inv_r2 * v34 - &
                           (v34 .o. v23)**2 * inv_r2 * inv_r2 * v23 
                      dp2(1:3,3) = -v34 + &
                           (v34 .o. v23) * inv_r2 * (v23 - v34) + &
                           (v34 .o. v23)**2 * inv_r2 * inv_r2 * v23
                      dp2(1:3,4) = v34 - (v34 .o. v23) * inv_r2 * v23
                      dp2 = dp2 * dot * inv_p2 * inv_p2
    
                      ! p.p' / (|p|^2|p'|) D |p| = p.p' / (2|p|^3|p'|) D (p.p)
                      dp1 = 0.d0
                      dp1(1:3,1) = -v12 + (v12 .o. v23) * inv_r2 * v23
                      dp1(1:3,2) = v12 - &
                           (v12 .o. v23) * inv_r2 * (v23 - v12) - &
                           (v12 .o. v23)**2 * inv_r2 * inv_r2 * v23
                      dp1(1:3,3) = &
                           -(v12 .o. v23) * inv_r2 * v12 + &
                           (v12 .o. v23)**2 * inv_r2 * inv_r2 * v23
                      dp1 = dp1 * dot * inv_p1 * inv_p1
    
                      phi = acos(cosine)
                      ratio = ratio / sqrt(1.d0 - cosine*cosine) * &
                           ( go%kphi1 * sin(phi - go%phi_chain(ii)) + &
                             go%kphi2 * 3.d0*sin( 3.d0*(phi - go%phi_chain(ii)) ) )
                      Fcomponent1(1:3) = ratio * ( ddot(1:3,1) + dp1(1:3,1) + dp2(1:3,1) )
                      Fcomponent2(1:3) = ratio * ( ddot(1:3,2) + dp1(1:3,2) + dp2(1:3,2) )
                      Fcomponent3(1:3) = ratio * ( ddot(1:3,3) + dp1(1:3,3) + dp2(1:3,3) )
                      Fcomponent4(1:3) = ratio * ( ddot(1:3,4) + dp1(1:3,4) + dp2(1:3,4) )
                                            
                      DO ia = 1, 3
                         IF(ABS(Fcomponent1(ia)) > limit1)THEN
                            mpi_g_force(ia,ii) = mpi_g_force(ia,ii) + Fcomponent1(ia)
                         ELSE IF(ABS(Fcomponent1(ia)) > limit2 .AND. ABS(Fcomponent1(ia)) <= limit1)THEN
                            small_g_force(ia,ii) = small_g_force(ia,ii) + Fcomponent1(ia)
                         ELSE IF(ABS(Fcomponent1(ia)) <= limit2)THEN
                            tiny_g_force(ia,ii) = tiny_g_force(ia,ii) + Fcomponent1(ia)
                         END IF
                         IF(ABS(Fcomponent2(ia)) > limit1)THEN
                            mpi_g_force(ia,ii+1) = mpi_g_force(ia,ii+1) + Fcomponent2(ia)
                         ELSE IF(ABS(Fcomponent2(ia)) > limit2 .AND. ABS(Fcomponent2(ia)) <= limit1)THEN
                            small_g_force(ia,ii+1) = small_g_force(ia,ii+1) + Fcomponent2(ia)
                         ELSE IF(ABS(Fcomponent2(ia)) <= limit2)THEN
                            tiny_g_force(ia,ii+1) = tiny_g_force(ia,ii+1) + Fcomponent2(ia)
                         END IF
                         IF(ABS(Fcomponent3(ia)) > limit1)THEN
                            mpi_g_force(ia,ii+2) = mpi_g_force(ia,ii+2) + Fcomponent3(ia)
                         ELSE IF(ABS(Fcomponent3(ia)) > limit2 .AND. ABS(Fcomponent3(ia)) <= limit1)THEN
                            small_g_force(ia,ii+2) = small_g_force(ia,ii+2) + Fcomponent3(ia)
                         ELSE IF(ABS(Fcomponent3(ia)) <= limit2)THEN
                            tiny_g_force(ia,ii+2) = tiny_g_force(ia,ii+2) + Fcomponent3(ia)
                         END IF
                         IF(ABS(Fcomponent4(ia)) > limit1)THEN
                            mpi_g_force(ia,ii+3) = mpi_g_force(ia,ii+3) + Fcomponent4(ia)
                         ELSE IF(ABS(Fcomponent4(ia)) > limit2 .AND. ABS(Fcomponent4(ia)) <= limit1)THEN
                            small_g_force(ia,ii+3) = small_g_force(ia,ii+3) + Fcomponent4(ia)
                         ELSE IF(ABS(Fcomponent4(ia)) <= limit2)THEN
                            tiny_g_force(ia,ii+3) = tiny_g_force(ia,ii+3) + Fcomponent4(ia)
                         END IF
                      END DO

                   END IF

                END IF

             END IF

             ! native and non-native bonds / repulsion
             DO jj = ii+4, go%length
                vij(1:3) = vector(gos(ii),gos(jj),cell,pbc) ! vector from ii+2 to ii+3
                hij = vij
                rij = sqrt(vij.o.vij) ! distance from ii+2 to ii+3
                inv_rij = 1.d0/rij    ! inverse distance from ii+2 to ii+3
                vij = vij*inv_rij     ! unit vector from ii+2 to ii+3
                IF(go%r_native(ii,jj) > 0.d0)THEN
                   ratio = go%r_native(ii,jj) * inv_rij
                   Fpair_ij = go%epsilon * inv_rij * ( 60.d0*ratio**10 - 60.d0 * ratio**12 )*vij
                ELSE
                   ratio = go%r_rep * inv_rij
                   ratio = ratio*ratio*ratio*ratio*ratio*ratio
                   Fpair_ij = 12.d0*go%epsilon * inv_rij * ( -ratio*ratio )*vij
                END IF

                ! add to the total force
                !WHERE(Fpair_ij < 1.0E-36) Fpair_ij = 0.d0
                DO ia = 1, 3
                   IF(ABS(Fpair_ij(ia)) > limit1)THEN
                      mpi_g_force(ia,ii) = mpi_g_force(ia,ii) + Fpair_ij(ia)
                      mpi_g_force(ia,jj) = mpi_g_force(ia,jj) - Fpair_ij(ia)
                   ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                      small_g_force(ia,ii) = small_g_force(ia,ii) + Fpair_ij(ia)
                      small_g_force(ia,jj) = small_g_force(ia,jj) - Fpair_ij(ia)
                   ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                      tiny_g_force(ia,ii) = tiny_g_force(ia,ii) + Fpair_ij(ia)
                      tiny_g_force(ia,jj) = tiny_g_force(ia,jj) - Fpair_ij(ia)
                   END IF
                END DO

                ! get the virial
                virial = virial - (Fpair_ij.o.hij)

             END DO

          END DO
       END IF

    END IF


    !!!!!!!!!!!!!!!!!!!!!!
    ! atom - atom forces !
    !!!!!!!!!!!!!!!!!!!!!!
    IF(first_mpi_atom > 0)THEN
       DO ii = first_mpi_atom, last_mpi_atom ! loop over the atoms of this cpu
          ns = aa_n_ns(ii) ! number of neighbors
          DO jj = 1, ns ! loop over neighbors
             nj = aa_nbors(jj,ii) ! jj:th neighbor of ii
             type = ats(ii)%type ! atomic types
             type2 = ats(nj)%type
             
             ! prevent double counting
             IF(ii < nj)THEN

                IF(params%n_pots(type2,type) > 0)THEN  ! type intercts with type2
                   
                   ! vector and distance between ii ans nj
                   vij(1:3) = vector(ats(ii),ats(nj),cell,pbc)
                   rij = sqrt(vij.o.vij)
                   inv_rij = 1.d0/rij
                   vij = vij*inv_rij
                   
                   IF(rij < params%pot_cut(type2,type))THEN ! cutoff

                      ! get the pairwise force
                      Fpair_ij = atomic_pair_force(vij,rij,inv_rij,params%n_pots(type2,type),&
                           params%pot_types(:,type2,type),params%pot_params(:,:,type2,type))
                      
                      DO ia = 1, 3
                         ! add to the total force
                         IF(ABS(Fpair_ij(ia)) > limit1)THEN
                            mpi_a_force(ia,ii) = mpi_a_force(ia,ii) + Fpair_ij(ia)
                            mpi_a_force(ia,nj) = mpi_a_force(ia,nj) - Fpair_ij(ia)
                         ELSE IF(ABS(Fpair_ij(ia)) > limit2 .AND. ABS(Fpair_ij(ia)) <= limit1)THEN
                            small_a_force(ia,ii) = small_a_force(ia,ii) + Fpair_ij(ia)
                            small_a_force(ia,nj) = small_a_force(ia,nj) - Fpair_ij(ia)
                         ELSE IF(ABS(Fpair_ij(ia)) <= limit2)THEN
                            tiny_a_force(ia,ii) = tiny_a_force(ia,ii) + Fpair_ij(ia)
                            tiny_a_force(ia,nj) = tiny_a_force(ia,nj) - Fpair_ij(ia)
                         END IF
                      END DO

!!$                      ! add to the total force
!!$                         mpi_a_force(1:3,ii) = mpi_a_force(1:3,ii) + Fpair_ij
!!$                         mpi_a_force(1:3,nj) = mpi_a_force(1:3,nj) - Fpair_ij

                      ! get the virial
                      virial = virial - (Fpair_ij.o.hij)

                   END IF

                END IF

             END IF
          END DO
       END DO

    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forces due to boundaries !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(first_mpi_mb_split > 0)THEN
       DO jj = 1, 3 ! loop over x, y, z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall boundary
             DO ii = first_mpi_mb_split, last_mpi_mb_split ! loop over mbs of this cpu
                IF(mbs(ii)%pos(jj) < 0.d0)THEN
                   mpi_m_force(jj,ii) = mpi_m_force(jj,ii) - bval(jj)*mbs(ii)%pos(jj)
                ELSE IF(mbs(ii)%pos(jj) > cell(jj))THEN
                   mpi_m_force(jj,ii) = mpi_m_force(jj,ii) - bval(jj)*(mbs(ii)%pos(jj)-cell(jj))
                END IF
             END DO
          END IF
       END DO
    END IF
    IF(first_mpi_atom_split > 0)THEN
       DO jj = 1, 3 ! loop over x, y, z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall boundary
             DO ii = first_mpi_atom_split, last_mpi_atom_split ! loop over the atoms of this cpu
                IF(ats(ii)%pos(jj) < 0.d0)THEN
                   mpi_a_force(jj,ii) = mpi_a_force(jj,ii) - bval(jj)*ats(ii)%pos(jj)
                ELSE IF(ats(ii)%pos(jj) > cell(jj))THEN
                   mpi_a_force(jj,ii) = mpi_a_force(jj,ii) - bval(jj)*(ats(ii)%pos(jj)-cell(jj))
                END IF
             END DO
          END IF
       END DO
    END IF
    IF(first_mpi_go_split > 0)THEN
       DO jj = 1, 3 ! loop over x, y, z
          IF(btype(jj) == wall_bound_index)THEN ! harmonic wall boundary
             DO ii = first_mpi_go_split, last_mpi_go_split ! loop over the atoms of this cpu
                IF(gos(ii)%pos(jj) < 0.d0)THEN
                   mpi_g_force(jj,ii) = mpi_g_force(jj,ii) - bval(jj)*gos(ii)%pos(jj)
                ELSE IF(gos(ii)%pos(jj) > cell(jj))THEN
                   mpi_g_force(jj,ii) = mpi_g_force(jj,ii) - bval(jj)*(gos(ii)%pos(jj)-cell(jj))
                END IF
             END DO
          END IF
       END DO
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forces due to constraints !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(is_constr)THEN ! ignore if there are no constraints
       IF(first_mpi_mb_split > 0)THEN
          DO ii = first_mpi_mb_split, last_mpi_mb_split ! loop over the mbs of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x, y, z
                SELECT CASE(mbs(ii)%constrained(jj))
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN ! get the displacement if needed
                      displ = displacement(mbs(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_m_force(jj,ii) = mpi_m_force(jj,ii) - mbs(ii)%well(jj)*(displ(jj))
                CASE(ext_force_index) ! external force
                   mpi_m_force(jj,ii) = mpi_m_force(jj,ii) + mbs(ii)%well(jj)
                END SELECT
             END DO
          END DO
       END IF

       IF(first_mpi_atom_split > 0)THEN
          DO ii = first_mpi_atom_split, last_mpi_atom_split ! loop over the atoms of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x, y, z
                SELECT CASE(ats(ii)%constrained(jj))
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN
                      displ = displacement(ats(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_a_force(jj,ii) = mpi_a_force(jj,ii) - ats(ii)%well(jj)*(displ(jj))
                CASE(ext_force_index) ! external force
                   mpi_a_force(jj,ii) = mpi_a_force(jj,ii) + ats(ii)%well(jj)
                END SELECT
             END DO
          END DO
       END IF

       IF(first_mpi_go_split > 0)THEN
          DO ii = first_mpi_go_split, last_mpi_go_split ! loop over the atoms of this cpu
             got_displ = .false.
             DO jj = 1, 3 ! loop over x, y, z
                SELECT CASE(gos(ii)%constrained(jj))
                CASE(harmonic_well_index) ! harmonic well
                   IF(.NOT.got_displ)THEN
                      displ = displacement(gos(ii),cell,pbc)
                      got_displ = .true.
                   END IF
                   mpi_g_force(jj,ii) = mpi_g_force(jj,ii) - gos(ii)%well(jj)*(displ(jj))
                CASE(ext_force_index) ! external force
                   mpi_g_force(jj,ii) = mpi_g_force(jj,ii) + gos(ii)%well(jj)
                END SELECT
             END DO
          END DO
       END IF
    END IF

#ifdef MPI
    ! mpi workload check and balancing
    CALL record_load(timer())
    CALL adjust_load()
#endif

    ! sum forces from all cpus
    IF(n_mb > 0)THEN
       ! add the big and small terms together
       small_m_force = small_m_force + tiny_m_force
       small_torque = small_torque + tiny_torque
       mpi_m_force = mpi_m_force + small_m_force
       mpi_torque = mpi_torque + small_torque
#ifdef MPI
       CALL MPI_ALLREDUCE(mpi_m_force,m_force,SIZE(m_force),MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,mpistat)
       CALL MPI_ALLREDUCE(mpi_torque,torque,SIZE(torque),MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,mpistat)
#else
       m_force = mpi_m_force
       torque = mpi_torque
#endif
    END IF

    IF(n_at > 0)THEN
       small_a_force = small_a_force + tiny_a_force
       mpi_a_force = mpi_a_force + small_a_force
#ifdef MPI
       CALL MPI_ALLREDUCE(mpi_a_force,a_force,SIZE(a_force),MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,mpistat)
#else
       a_force = mpi_a_force
#endif
    END IF


    IF(n_go > 0)THEN
       small_g_force = small_g_force + tiny_g_force
       mpi_g_force = mpi_g_force + small_g_force
#ifdef MPI
       CALL MPI_ALLREDUCE(mpi_g_force,g_force,SIZE(g_force),MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,mpistat)
#else
       g_force = mpi_g_force
#endif
    END IF

#ifdef MPI
    CALL MPI_ALLREDUCE(virial,thevirial,1,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,mpistat)
#else
    thevirial = virial
#endif


!!$    ! sum the small terms
!!$    CALL MPI_ALLREDUCE(small_m_force,mpi_m_force,SIZE(m_force),MPI_DOUBLE_PRECISION,&
!!$         MPI_SUM,MPI_COMM_WORLD,mpistat)
!!$    CALL MPI_ALLREDUCE(small_torque,mpi_torque,SIZE(torque),MPI_DOUBLE_PRECISION,&
!!$         MPI_SUM,MPI_COMM_WORLD,mpistat)
!!$    CALL MPI_ALLREDUCE(small_a_force,mpi_a_force,SIZE(a_force),MPI_DOUBLE_PRECISION,&
!!$         MPI_SUM,MPI_COMM_WORLD,mpistat)
!!$
!!$    ! add the big and small terms together
!!$    m_force = mpi_m_force + m_force
!!$    torque = mpi_torque + torque
!!$    a_force = mpi_a_force + a_force

    ! if bonding stats are recorded, sum them as well
    IF(control%bond_writer /= noxyz_index)THEN
#ifdef MPI
       CALL MPI_REDUCE(mpi_n_bond,n_bond,SIZE(mpi_n_bond),MPI_DOUBLE_PRECISION,MPI_SUM,master_cpu,&
            MPI_COMM_WORLD,mpistat) 
#else
       n_bond = mpi_n_bond
#endif
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forces due to thermostat !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(apply_thermostat)THEN
       
       ! each cpu should add the thermostat by itself, since langevin uses random numbers
       ! and the cpus have identical rnd sequences.
       ! if the random force calculations were split, the forces would be correlated
       CALL add_thermostat_to_forces(mbs,ats,gos,n_mb,n_at,n_go,&
            control,params,go,m_force,torque,a_force,g_force)

    END IF

    virial = thevirial/3.d0

    RETURN
  END SUBROUTINE calc_forces

  ! Adds the effect of a thermostat to the forces (so that one can, e.g., get real forces
  ! first and then later add the thermostat component).
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mb number of molecules
  ! *n_at number of atoms
  ! *control control parameters
  ! *params physical parameters
  ! *m_force forces on molecules
  ! *torque torques on molecules
  ! *a_force forces on atoms
  SUBROUTINE add_thermostat_to_forces(mbs,ats,gos,n_mb,n_at,n_go,&
       control,params,go,m_force,torque,a_force,g_force)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    INTEGER, INTENT(IN) :: n_mb, n_at, n_go
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), torque(:,:), a_force(:,:), g_force(:,:)

    IF(control%md_thermo == langevin_index)THEN ! Langevin thermostat
       CALL add_random_force(ats,n_mb,n_at,n_go,control,params,go,m_force,torque,a_force,g_force)
       CALL add_friction_force(mbs,ats,gos,n_mb,n_at,n_go,control,go,m_force,torque,a_force,g_force)
    ELSE IF(control%md_thermo == cooler_index)THEN ! only friction to cool down the system
       CALL add_friction_force(mbs,ats,gos,n_mb,n_at,n_go,control,go,m_force,torque,a_force,g_force)
    END IF
    
    RETURN
  END SUBROUTINE add_thermostat_to_forces

  ! Adds a friction component to forces - for langevin and "cooler" thermostats.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mb number of molecules
  ! *n_at number of atoms
  ! *control control parameters
  ! *params physical parameters
  ! *m_force forces on molecules
  ! *torque torques on molecules
  ! *a_force forces on atoms
  SUBROUTINE add_friction_force(mbs,ats,gos,n_mb,n_at,n_go,control,go,m_force,torque,a_force,g_force)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    TYPE(gop), POINTER :: gos(:)
    INTEGER, INTENT(IN) :: n_mb, n_at, n_go
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), torque(:,:), a_force(:,:), g_force(:,:)
    INTEGER :: ii

    DO ii = 1, n_mb ! loop over mbs
       m_force(1:3,ii) = m_force(1:3,ii) - mbs(ii)%m_tot*control%thermo_value*mbs(ii)%vel(:) 
       torque(1:3,ii) = torque(1:3,ii) - mbs(ii)%m_inert*control%thermo_value*mbs(ii)%angvel(:)
    END DO
    DO ii = 1, n_at ! loop over atoms
       a_force(1:3,ii) = a_force(1:3,ii) - ats(ii)%mass*control%thermo_value*ats(ii)%vel(:)
    END DO
    DO ii = 1, n_go ! loop over gos
       g_force(1:3,ii) = g_force(1:3,ii) - gos(ii)%mass*control%thermo_value*gos(ii)%vel(:)
    END DO
    
    RETURN
  END SUBROUTINE add_friction_force

  ! Adds the effect of a random force - for langevin thermostat.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mb number of molecules
  ! *n_at number of atoms
  ! *control control parameters
  ! *params physical parameters
  ! *m_force forces on molecules
  ! *torque torques on molecules
  ! *a_force forces on atoms
  SUBROUTINE add_random_force(ats,n_mb,n_at,n_go,control,params,go,m_force,torque,a_force,g_force)
    IMPLICIT NONE
    !TYPE(mb), POINTER :: mbs(:)
    TYPE(atom), POINTER :: ats(:)
    INTEGER, INTENT(IN) :: n_mb, n_at, n_go
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    TYPE(gops), INTENT(IN) :: go
    REAL(KIND=dp), POINTER :: m_force(:,:), torque(:,:), a_force(:,:), g_force(:,:)
    REAL(KIND=dp) :: langevin1, langevin2, random, rnd_vector(3), spare_vector(3)
    INTEGER :: ii, jj
    LOGICAL :: got_spare

    IF(.true.)THEN ! uniformly distributed components for the random forces - faster
       langevin1 = sqrt(24.d0*params%m_mol*control%thermo_value*kb*control%ini_temp/control%time_step) ! range of random forces
       langevin2 = sqrt(24.d0*params%i_mol*control%thermo_value*kb*control%ini_temp/control%time_step) ! range of random torques
       DO ii = 1, n_mb ! loop over mbs
          DO jj = 1, 3 ! loop over x, y, z
             CALL genrand_real1(random)
             m_force(jj,ii) = m_force(jj,ii) + langevin1*(random-0.5d0) ! add the random component for force
             CALL genrand_real1(random)
             torque(jj,ii) = torque(jj,ii) + langevin2*(random-0.5d0) ! add the random component for torque
          END DO
       END DO
       langevin1 = sqrt(24.d0*control%thermo_value*kb*control%ini_temp/control%time_step) ! range of random forces
       DO ii = 1, n_at ! loop over atoms
          DO jj = 1, 3 ! loop over x, y, z
             CALL genrand_real1(random)
             a_force(jj,ii) = a_force(jj,ii) + langevin1*sqrt(params%m_atoms(ats(ii)%type))*(random-0.5) ! add the random component for force
          END DO
       END DO
       DO ii = 1, n_go ! loop over atoms
          DO jj = 1, 3 ! loop over x, y, z
             CALL genrand_real1(random)
             g_force(jj,ii) = g_force(jj,ii) + langevin1*sqrt(go%mass_chain(ii))*(random-0.5) ! add the random component for force
          END DO
       END DO
    ELSE ! gaussian random forces - it's probably nicer than uniformly distributed, but also slower - in principle both work
       rnd_vector = 0.d0
       spare_vector = 0.d0
       langevin1 = sqrt(2.d0*params%m_mol*control%thermo_value*kb*control%ini_temp/control%time_step)
       langevin2 = sqrt(2.d0*params%i_mol*control%thermo_value*kb*control%ini_temp/control%time_step)
       DO ii = 1, n_mb
          CALL rand_2boltzmann_vectors(rnd_vector,spare_vector)
          m_force(1:3,ii) = m_force(1:3,ii) + langevin1*rnd_vector(:)
          torque(1:3,ii) = torque(1:3,ii) + langevin2*spare_vector(:)
       END DO
       langevin1 = 2.d0*control%thermo_value*kb*control%ini_temp/control%time_step
       DO ii = 1, n_at
          got_spare = .false.
          IF(got_spare)THEN
             a_force(1:3,ii) = a_force(1:3,ii) + &
                  sqrt(ats(ii)%mass*langevin1)*spare_vector(:)
             got_spare = .false.
          ELSE
             CALL rand_2boltzmann_vectors(rnd_vector,spare_vector)
             a_force(1:3,ii) = a_force(1:3,ii) + &
                  sqrt(ats(ii)%mass*langevin1)*rnd_vector(:)
             got_spare = .true.
          END IF
       END DO
    END IF
    
    RETURN
  END SUBROUTINE add_random_force


  ! Returns the pairwise potential energy for Mb-atom or atom-atom
  ! interactions defined by the distance between the particles 
  ! and the parameters for the potential.
  ! *rij distance between the particles
  ! *n_pots number of potentials defined
  ! *pot_types list for the types of potentials, denoted as integer indices
  ! *parameters list for the parameters of the potentials
  ! *potene the pair potential energy
  FUNCTION atomic_pair_potential(rij,n_pots,pot_types,parameters) &
       RESULT(potene)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: rij
    REAL(KIND=dp), INTENT(IN) :: parameters(:,:)
    INTEGER, INTENT(IN) :: n_pots
    INTEGER, INTENT(IN) :: pot_types(:)
    REAL(KIND=dp) :: potene
    REAL(KIND=dp) :: range, epsil
    INTEGER :: kk, pot_type

    potene = 0.d0
    ! sum over the defined potentials
    DO kk = 1, n_pots
       pot_type = pot_types(kk) ! type of potential
       IF(pot_type == lenjon_index)THEN ! lennard-jones
          epsil = parameters(1,kk) ! e_lj
          range = parameters(2,kk)/rij ! sigma/rij
          potene = potene + 4.d0*epsil*( range**12 - range**6 )

       ELSE IF(pot_type == exp_index)THEN ! exponential
          epsil = parameters(1,kk) ! e_exp
          range = rij/parameters(2,kk) ! rij/rho
          potene = potene + epsil*exp(-range)

       ELSE IF(pot_type == pow_index)THEN ! power
          epsil = parameters(1,kk) ! e_pow
          range = rij**(-parameters(2,kk)) ! rij**(-n_pow)
          potene = potene + epsil*range

       ELSE IF(pot_type == spring_index)THEN ! spring
          epsil = parameters(1,kk) ! k_spring
          range = (rij-parameters(2,kk)) ! r-r0
          potene = potene + 0.5d0*epsil*range*range

       ELSE IF(pot_type == hard_index)THEN ! hard-sphere + LJ
          epsil = parameters(1,kk) ! e_lj
          IF(rij <= parameters(3,kk)) CALL abort("particle too close to a hard-sphere")
          range = parameters(2,kk)/(rij-parameters(3,kk)) ! sigma/(r-r0)
          potene = potene + 4.d0*epsil*( range**12 - range**6 )

       ELSE IF(pot_type == hardrep_index)THEN ! repulsive hard-sphere
          IF(rij < parameters(3,kk) + 1.22462058d0*parameters(2,kk))THEN
             epsil = parameters(1,kk) ! e_lj
             IF(rij <= parameters(3,kk)) CALL abort("particle too close to a hard-sphere")
             range = parameters(2,kk)/(rij-parameters(3,kk)) ! sigma/(r-r0)
             potene = potene + 4.d0*epsil*( range**12 - range**6 )+epsil
          END IF

       ELSE IF(pot_type == fene_index)THEN ! fene
          epsil = parameters(1,kk)*parameters(2,kk)*parameters(2,kk) ! k*R0**2
          IF(ABS(rij-parameters(3,kk)) >= parameters(2,kk)) CALL abort("particles outside FENE-potential range")
          range = (rij-parameters(3,kk))/parameters(2,kk) ! (r-r0)/R0
          potene = potene - 0.5d0*epsil*LOG(1-range*range)

       ELSE IF(pot_type == shell_index)THEN ! shell
          
          IF(rij < parameters(2,kk))THEN ! close -> fene spring
             epsil = parameters(1,kk)*parameters(2,kk)*parameters(2,kk)
             range = rij/parameters(2,kk)
             potene = potene - 0.5d0*epsil*LOG(1-range*range)
          ELSE ! far -> electrostatic
             potene = potene + parameters(3,kk)/rij
          END IF

       ELSE
          CALL abort("undefined potential encountered while calculating energy")
       END IF
    END DO

    RETURN
  END FUNCTION  atomic_pair_potential


  ! Returns the pairwise force for Mb-atom or atom-atom interactions defined by
  ! distance and potential parameteres.
  ! *rij distance between the particles
  ! *vij unit vector pointing from particle one to particle two - the calculated force will be the force acting on particle one
  ! *inv_rij 1/rij
  ! *n_pots number of potentials defined
  ! *pot_types list for the types of potentials, denoted as integer indices
  ! *parameters list for the parameters of the potentials
  ! *pforce the force acting on the particle from which vij points away
  FUNCTION atomic_pair_force(vij,rij,inv_rij,n_pots,pot_types,parameters) &
       RESULT(pforce)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: vij(3), rij, inv_rij
    REAL(KIND=dp), INTENT(IN) :: parameters(:,:)
    INTEGER, INTENT(IN) :: n_pots
    INTEGER, INTENT(IN) :: pot_types(:)
    REAL(KIND=dp) :: pforce(3)
    REAL(KIND=dp) :: range, epsil, inver
    INTEGER :: kk, pot_type

    pforce = 0.d0
    ! sum over the defined potentials
    DO kk = 1, n_pots
       pot_type = pot_types(kk) ! potential type
       IF(pot_type == lenjon_index)THEN ! lennard-jones
          epsil = parameters(1,kk) ! e_lj
          range = parameters(2,kk)*inv_rij ! sigma/rij
          ! E: 4.d0*epsil*( range**12 - range**6 )
          pforce = pforce + 24.d0*epsil*inv_rij* &
               ( range**6 - 2.d0*range**12 )*vij(1:3)

       ELSE IF(pot_type == exp_index)THEN ! exponential
          epsil = parameters(1,kk) ! e_exp
          range = 1.d0/parameters(2,kk) ! 1/rho
          ! E: epsil*exp(-range*rij)
          pforce = pforce - epsil*exp(-rij*range)*range*vij(1:3)

       ELSE IF(pot_type == pow_index)THEN ! power
          epsil = parameters(1,kk) ! e_pow
          range = parameters(2,kk) ! n_pow
          ! E: + epsil*inv_rij**range
          pforce = pforce - epsil*range*inv_rij**(range+1.d0)*vij(1:3) 

       ELSE IF(pot_type == spring_index)THEN ! spring
          epsil = parameters(1,kk) ! k_spring
          range = rij-parameters(2,kk) ! r-r0          
          ! E: + epsil/2*range**2
          pforce = pforce + epsil*range*vij(1:3) 

       ELSE IF(pot_type == hard_index)THEN ! hard-sphere + LJ
          epsil = parameters(1,kk) ! e_lj
          IF(rij <= parameters(3,kk)) CALL abort("particle too close to a hard-sphere")
          inver = 1/(rij-parameters(3,kk))
          range = parameters(2,kk)*inver ! sigma/(r-r0)
          pforce = pforce + 24.d0*epsil*inver* &
               ( range**6 - 2.d0*range**12 )*vij(1:3)

       ELSE IF(pot_type == hardrep_index)THEN ! repulsive hard-sphere
          IF(rij < parameters(3,kk) + 1.22462048*parameters(2,kk))THEN
             epsil = parameters(1,kk) ! e_lj
             IF(rij <= parameters(3,kk))THEN
                CALL abort("particle too close to a hard-sphere")
             END IF
             inver = 1/(rij-parameters(3,kk))
             range = parameters(2,kk)*inver ! sigma/(r-r0)
             pforce = pforce + 24.d0*epsil*inver* &
                  ( range**6 - 2.d0*range**12 )*vij(1:3)
          END IF

       ELSE IF(pot_type == fene_index)THEN ! fene
          epsil = parameters(1,kk) ! k
          IF(ABS(rij-parameters(3,kk)) >= parameters(2,kk)) CALL abort("particles outside FENE-potential range")
          range = (rij-parameters(3,kk)) ! (r-r0)
          pforce = pforce + epsil*range/(1.d0-(range/parameters(2,kk))**2)*vij(1:3)

       ELSE IF(pot_type == shell_index)THEN ! shell

          IF(rij < parameters(2,kk))THEN ! close -> fene spring
             epsil = parameters(1,kk) ! k
             pforce = pforce + epsil*rij/(1.d0-(rij/parameters(2,kk))**2)*vij(1:3)
          ELSE ! far -> electrostatic
             ! E: + epsil*inv_rij**2
             pforce = pforce - parameters(3,kk)*inv_rij*inv_rij*vij(1:3) 
          END IF

       ELSE
          CALL abort("undefined potential encountered while calculating energy")
       END IF
    END DO

    RETURN
  END FUNCTION  atomic_pair_force

  ! a more secure way to get the size of the pointer array for mb molecules
  ! *thesize the size of the array
  ! *mbs the pointer list of MB molecules
  FUNCTION mbs_size(mbs) &
       RESULT(thesize)
    IMPLICIT NONE
    TYPE(mb), POINTER :: mbs(:)
    INTEGER :: thesize

    IF(associated(mbs))THEN
       thesize = SIZE(mbs)
    ELSE
       thesize = 0
    END IF

    RETURN
  END FUNCTION mbs_size

  ! a more secure way to get the size of the pointer array for atoms
  ! *thesize the size of the array
  ! *ats the pointer list of atoms
  FUNCTION ats_size(ats) &
       RESULT(thesize)
    IMPLICIT NONE
    TYPE(atom), POINTER :: ats(:)
    INTEGER :: thesize

    IF(associated(ats))THEN
       thesize = SIZE(ats)
    ELSE
       thesize = 0
    END IF

    RETURN
  END FUNCTION ats_size


  ! a more secure way to get the size of the pointer array for gos
  ! *thesize the size of the array
  ! *gos the pointer list of gos
  FUNCTION gos_size(gos) &
       RESULT(thesize)
    IMPLICIT NONE
    TYPE(gop), POINTER :: gos(:)
    INTEGER :: thesize

    IF(associated(gos))THEN
       thesize = SIZE(gos)
    ELSE
       thesize = 0
    END IF

    RETURN
  END FUNCTION gos_size



  ! Calculates the rate of proximity [0,1] for two molecules.
  ! That is, this is the bond number function f(r) for
  ! the distance r=riq.
  ! *riq distance between two molecules
  ! *Rb range for bonding
  ! *Db (half of) the range in which the function drops from 1 to 0
  ! *inv_db 1/Db
  ! *fiq rate of proximity, f(riq)
  FUNCTION proximity(riq,Rb,Db,inv_Db) &
       RESULT(fiq)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: riq, Rb, Db, inv_Db
    REAL(KIND=dp) :: fiq

    IF(riq >= Rb+Db)THEN ! f = 0 if r is big
       fiq = 0.d0
    ELSE IF(riq <= Rb-Db)THEN ! f = 1 if r is small
       fiq = 1.d0
    ELSE ! in the shell area, f changes continuously from 1 to 0
       fiq = 0.5d0*(1.d0 - sin( inv_Db*(riq-Rb) ))
    END IF

  END FUNCTION proximity

  ! Returns the bond order factor for z_i bonds according to exponent v_exp.
  ! That is, this is the tersoff coefficient function b_i(z) for z=z_i.
  ! *zi number of bonds (molecules in proximity)
  ! *v_exp scaling exponent for the bond order factor
  ! *bi the bond order (tersoff) coefficient, b_i
  FUNCTION bond_order(zi,v_exp) &
       RESULT(bi)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: zi, v_exp
    REAL(KIND=dp) :: bi

    IF(zi <= 4.d0)THEN ! for z less than 4, b_i is 1
       bi = 1.d0
       RETURN
    ELSE ! for z more than 4, b_i is less than 1
       bi = (4.d0/zi)**v_exp
       RETURN
    END IF

  END FUNCTION bond_order

  ! Calculates the radial distribution function.
  ! The parameters for RDF generating are:
  ! (i) cutrange; determines the length for which the RDF is measured.
  ! (ii) particle_index; determines the particles considered - the value "all_rdf_index"
  ! means all particles are taken into account while for the values
  ! "mb_rdf_index" and "atom_rdf_index" only MB molecules and atoms are considered,
  ! respectively.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mbs number of molecules
  ! *n_ats number of atoms
  ! *cell supercell dimensions
  ! *pbc true for periodicc boundaries
  ! *control control parameters
  ! *particle_index determines which particles are taken into account: all, mb molecules, or atoms
  ! *cutrange the range for measuring the RDF
  ! *rdf the calculated rdfs: (1) mb-mb rdf, (2) atom-atom rdf, (3) mb-atom rdf, (4) total rdf
  SUBROUTINE radial_distribution_function(mbs,ats,n_mbs,n_ats,cell,pbc,&
       mm_nbors,mm_n_nbor,ma_nbors,ma_n_nbor,aa_nbors,aa_n_nbor,params,&
       particle_index,cutrange,rdf)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mbs(:)
    TYPE(atom), INTENT(IN) :: ats(:)
    TYPE(mbps), INTENT(IN) :: params
    INTEGER, INTENT(IN) :: particle_index, n_mbs, n_ats
    INTEGER, POINTER :: mm_nbors(:,:), mm_n_nbor(:), &
         ma_nbors(:,:), ma_n_nbor(:), aa_nbors(:,:), aa_n_nbor(:)
    REAL(KIND=dp), INTENT(IN) :: cutrange, cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    REAL(KIND=dp), INTENT(OUT) :: rdf(rdfbins,4)
    INTEGER :: ii, jj, nj, n_dists(4), binindex
    REAL(KIND=dp) :: binwidth, dd, rr(3), normer

    ! this routine could be parallellized but hasn't been!

    binwidth = cutrange/REAL(rdfbins,KIND=dp) ! width of the bins used for statistics gathering

    rdf(:,:) = 0.d0
    n_dists = 0

    IF(particle_index == mb_rdf_index .OR. & ! measure rdf for mbs
         particle_index == all_rdf_index )THEN

       IF(cutrange <= params%max_mb_cut-params%cut_ver)THEN ! checking only particles in the neighbor list

          DO ii = 1, n_mbs-1 ! loop over mbs
             DO jj = 1, mm_n_nbor(ii) ! loop over neighbors
                nj = mm_nbors(jj,ii) ! jj:th neighbor

                IF(ii < nj)THEN ! prevent double counting
                   rr = vector(mbs(ii),mbs(nj),cell,pbc)
                   dd = .norm.rr
                   
                   binindex = FLOOR(dd/binwidth)+1 ! the index of the bin corresponding to the separation dd
                   
                   IF(binindex <= rdfbins)THEN ! check for being outside the plot range
                      ! record the pair
                      n_dists(1) = n_dists(1) + 1
                      n_dists(4) = n_dists(4) + 1
                      rdf(binindex,1) = rdf(binindex,1) + 1.d0
                      rdf(binindex,4) = rdf(binindex,4) + 1.d0
                   END IF
                END IF
             END DO
          END DO

       ELSE ! check the whole cell
          DO ii = 1, n_mbs-1 ! loop over mbs
             DO jj = ii+1, n_mbs ! loop over other mbs (ii < jj)
                
                rr = vector(mbs(ii),mbs(jj),cell,pbc)
                dd = sqrt(rr .o. rr)
                
                binindex = FLOOR(dd/binwidth)+1
                
                IF(binindex <= rdfbins)THEN
                   n_dists(1) = n_dists(1) + 1
                   n_dists(4) = n_dists(4) + 1
                   rdf(binindex,1) = rdf(binindex,1) + 1.d0
                   rdf(binindex,4) = rdf(binindex,4) + 1.d0
                END IF
                
             END DO
          END DO
       END IF

    END IF


    IF(particle_index == atom_rdf_index .OR. & ! measure rdf for atoms
         particle_index == all_rdf_index )THEN

       IF(cutrange <= params%max_at_cut-params%cut_ver)THEN ! check only neighbors

          DO ii = 1, n_ats-1 ! loop over atoms
             DO jj = 1, aa_n_nbor(ii) ! loop over neighbors
                nj = aa_nbors(jj,ii) ! jj:th neighbor

                IF(ii < nj)THEN ! prevent double counting
                   rr = vector(ats(ii),ats(nj),cell,pbc)
                   dd = .norm.rr
                   
                   binindex = FLOOR(dd/binwidth)+1 ! index of the bin for dd
                   
                   IF(binindex <= rdfbins)THEN
                      n_dists(2) = n_dists(2) + 1
                      rdf(binindex,2) = rdf(binindex,2) + 1.d0
                      n_dists(4) = n_dists(4) + 1
                      rdf(binindex,4) = rdf(binindex,4) + 1.d0
                   END IF
                END IF
             END DO
          END DO

       ELSE ! check the whole cell

          DO ii = 1, n_ats-1 ! loop over atoms
             DO jj = ii+1, n_ats ! loop over other atoms
                
                rr = vector(ats(ii),ats(jj),cell,pbc)
                dd = sqrt(rr .o. rr)
                
                binindex = FLOOR(dd/binwidth)+1
                IF(binindex <= rdfbins)THEN
                   n_dists(2) = n_dists(2) + 1
                   rdf(binindex,2) = rdf(binindex,2) + 1.d0
                   n_dists(4) = n_dists(4) + 1
                   rdf(binindex,4) = rdf(binindex,4) + 1.d0
                END IF
                
             END DO
          END DO
       END IF
    END IF


    IF(particle_index == all_rdf_index )THEN ! measure rdf for mb-atom pairs

       IF(cutrange <= params%max_at_cut-params%cut_ver)THEN ! check only neighbors

          DO ii = 1, n_mbs ! loop over mbs
             DO jj = 1, ma_n_nbor(ii) ! loop over atom neighbors
                nj = ma_nbors(jj,ii) ! jj:th neighbor

                rr = vector(mbs(ii),ats(nj),cell,pbc)
                dd = .norm.rr
                   
                binindex = FLOOR(dd/binwidth)+1
                   
                IF(binindex <= rdfbins)THEN
                   n_dists(3) = n_dists(3) + 1
                   rdf(binindex,3) = rdf(binindex,3) + 1.d0
                   n_dists(4) = n_dists(4) + 1
                   rdf(binindex,4) = rdf(binindex,4) + 1.d0
                END IF
             
             END DO
          END DO
       
       ELSE

          DO ii = 1, n_mbs ! loop over mbs
             DO jj = 1, n_ats ! loop over atoms
                
                rr = vector(mbs(ii),ats(jj),cell,pbc)
                dd = sqrt(rr .o. rr)
                
                binindex = FLOOR(dd/binwidth)+1
                IF(binindex <= rdfbins)THEN
                   n_dists(3) = n_dists(3) + 1
                   rdf(binindex,3) = rdf(binindex,3) + 1.d0
                   n_dists(4) = n_dists(4) + 1
                   rdf(binindex,4) = rdf(binindex,4) + 1.d0
                END IF
                
             END DO
          END DO
       END IF
    END IF

    DO ii = 1, 4 ! loop over the four cases: mb-mb, at-at, mb-at, all-all
       IF(n_dists(ii) > 0)THEN
          ! get the normalizing coefficients
          SELECT CASE(ii)
          CASE(1) !mb-mb
             normer = cell(1)*cell(2)*cell(3)/REAL((n_mbs)*(n_mbs-1),KIND=dp)/binwidth/(2.d0*pi)
          CASE(2) !at-at
             normer = cell(1)*cell(2)*cell(3)/REAL((n_ats)*(n_ats-1),KIND=dp)/binwidth/(2.d0*pi)
          CASE(3) !mb-at
             normer = cell(1)*cell(2)*cell(3)/REAL(n_mbs*n_ats,KIND=dp)/binwidth/(2.d0*pi)
          CASE(4) !any-any
             normer = cell(1)*cell(2)*cell(3)/REAL((n_mbs+n_ats)*(n_mbs+n_ats-1),KIND=dp)/binwidth/(4.d0*pi)
          END SELECT
          DO jj = 1, rdfbins ! resord the normalized rdf
             rdf(jj,ii) = rdf(jj,ii)*normer/((REAL(jj,KIND=dp)-0.5d0)*binwidth)**2
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE radial_distribution_function

  ! Calculates the axial distribution function. This is similar to the radial distribution function
  ! except it records the statistics in cylindric symmetry around a specified axis, not in spherical symmetry
  ! around all particles.
  ! *mbs list of molecules
  ! *ats list of atoms
  ! *n_mbs number of molecules
  ! *n_ats number of atoms
  ! *cell supercell dimensions
  ! *pbc true for periodicc boundaries
  ! *params physical parameters
  ! *control control parameters
  ! *particle_index determines which particles are taken into account: all, mb molecules, or atoms
  ! *adf the calculated adfs: (1) mb-mb adf, (2) atom-atom adf, (3) mb-atom adf, (4) total adf
  SUBROUTINE axial_distribution_function(mbs,ats,n_mbs,n_ats,cell,pbc,&
       params,control,adf)
    IMPLICIT NONE
    TYPE(mb), INTENT(IN) :: mbs(:)
    TYPE(atom), INTENT(IN) :: ats(:)
    TYPE(mbps), INTENT(IN) :: params
    TYPE(cps), INTENT(IN) :: control
    INTEGER, INTENT(IN) :: n_mbs, n_ats
    REAL(KIND=dp), INTENT(IN) :: cell(3)
    LOGICAL, INTENT(IN) :: pbc(3)
    REAL(KIND=dp), INTENT(OUT) :: adf(rdfbins,rdfbins,3)
    INTEGER :: ii, jj, n_dists(3), binindex(2), particle_index
    REAL(KIND=dp) :: binwidth(2), dd, angle, rr(3), normer, center(3), &
         axis(3), angleaxis(3), posaxis(3)

    binwidth(1) = control%adf_range/REAL(rdfbins,KIND=dp) ! bin width in the radial direction
    binwidth(2) = twopi/REAL(rdfbins,KIND=dp) ! bin width in the angular direction
    adf(:,:,:) = 0.d0
    n_dists = 0
    particle_index = control%adf_particles
    ! determine a point through which the symmetry axis should go
    IF(control%adf_follow == 0)THEN
       center = control%adf_center
    ELSE IF(control%adf_follow < 0)THEN
       center = ats(-control%adf_follow)%pos
    ELSE
       center = mbs(-control%adf_follow)%pos
    END IF
    ! pick the direction of the symmetry axis
    SELECT CASE(control%adf_axis)
    CASE(x_index)
       axis = (/ 1.d0, 0.d0, 0.d0 /) ! symmetry axis
       angleaxis = (/ 0.d0, 0.d0, 1.d0 /) ! direction of zero angle
    CASE(y_index)
       axis = (/ 0.d0, 1.d0, 0.d0 /)    
       angleaxis = (/ 0.d0, 0.d0, 1.d0 /)   
    CASE(z_index)
       axis = (/ 0.d0, 0.d0, 1.d0 /)
       angleaxis = (/ 1.d0, 0.d0, 0.d0 /)
    CASE DEFAULT
       CALL abort("axis missing for adf")
    END SELECT
    ! direction of growing angle
    posaxis = (axis .x. angleaxis)

    IF(particle_index == mb_rdf_index .OR. & ! mbs
         particle_index == all_rdf_index )THEN

       DO ii = 1, n_mbs ! loop over mbs
          
          ! displacement vector
          rr = mbs(ii)%pos - center               
          ! perpendicular part
          rr = rr - (rr .o. axis)*axis
          ! periodicity
          DO jj = 1, 3
             IF(pbc(jj))THEN
                IF(rr(jj) > cell(jj)*0.5d0)THEN
                   rr(jj) = rr(jj) - cell(jj)
                ELSE IF(rr(jj) < -cell(jj)*0.5d0)THEN
                   rr(jj) = rr(jj) + cell(jj)
                END IF
             END IF
          END DO
          
          dd = sqrt(rr .o. rr)             
          binindex(1) = FLOOR(dd/binwidth(1))+1 ! index of the bin corresponding to dd
          
          IF(binindex(1) <= rdfbins)THEN
             
             ! track angle
             IF(control%adf_angle)THEN
                ! get the angle
                IF(dd < norm_tolerance)THEN ! at the origin, set angle to zero
                   angle = 0.d0
                ELSE IF((rr.o.posaxis) > 0.d0)THEN ! the angle is positive
                   angle = (rr.o.angleaxis)/dd
                   IF(angle >= 1.d0)THEN ! check for numeric error
                      angle = 0.d0
                   ELSE IF(angle <= -1.d0)THEN
                      angle = pi
                   ELSE
                      angle = acos(angle)
                   END IF
                ELSE ! the angle is negative
                   angle = (rr.o.angleaxis)/dd
                   IF(angle >= 1.d0)THEN
                      angle = twopi
                   ELSE IF(angle <= -1.d0)THEN
                      angle = pi
                   ELSE
                      angle = twopi - acos(angle)
                   END IF
                END IF
                
                binindex(2) = FLOOR(angle/binwidth(2))+1 ! index of the bin corresponding to angle
                ! If angle = 2pi, binindex is rdfbins+1, which is not good.
                ! However, 2pi = 0, so let's make the index 1.
                IF(binindex(2) > rdfbins)THEN
                   binindex(2) = 1
                END IF
             ELSE ! don't track angle at all
                binindex(2) = 1
             END IF

             ! record the particle
             n_dists(1) = n_dists(1) + 1
             n_dists(3) = n_dists(3) + 1
             adf(binindex(2),binindex(1),1) = adf(binindex(2),binindex(1),1) + 1.d0
             adf(binindex(2),binindex(1),3) = adf(binindex(2),binindex(1),3) + 1.d0
          END IF
          
       END DO
    END IF

    IF(particle_index == atom_rdf_index .OR. & ! atoms
         particle_index == all_rdf_index )THEN

       DO ii = 1, n_ats
                
          ! displacement vector
          rr = ats(ii)%pos - center               
          ! perpendicular part
          rr = rr - (rr .o. axis)*axis
          ! periodicity
          DO jj = 1, 3
             IF(pbc(jj))THEN
                IF(rr(jj) > cell(jj)*0.5d0)THEN
                   rr(jj) = rr(jj) - cell(jj)
                ELSE IF(rr(jj) < -cell(jj)*0.5d0)THEN
                   rr(jj) = rr(jj) + cell(jj)
                END IF
             END IF
          END DO

          dd = sqrt(rr .o. rr)             
          binindex(1) = FLOOR(dd/binwidth(1))+1
          
          IF(binindex(1) <= rdfbins)THEN
             
             ! track angle
             IF(control%adf_angle)THEN

                ! get the angle
                IF(dd < norm_tolerance)THEN
                   angle = 0.d0
                ELSE IF((rr.o.posaxis) > 0.d0)THEN
                   angle = (rr.o.angleaxis)/dd
                   IF(angle >= 1.d0)THEN
                      angle = 0.d0
                   ELSE IF(angle <= -1.d0)THEN
                      angle = pi
                   ELSE
                      angle = acos(angle)
                   END IF
                ELSE
                   angle = (rr.o.angleaxis)/dd
                   IF(angle >= 1.d0)THEN
                      angle = twopi
                   ELSE IF(angle <= -1.d0)THEN
                      angle = pi
                   ELSE
                      angle = twopi - acos(angle)
                   END IF
                END IF
                
                binindex(2) = FLOOR(angle/binwidth(2))+1
                ! If angle = 2pi, binindex is rdfbins+1, which is not good.
                ! However, 2pi = 0, so let's make the index 1.
                IF(binindex(2) > rdfbins)THEN
                   binindex(2) = 1
                END IF
             ELSE
                binindex(2) = 1
             END IF

             n_dists(2) = n_dists(2) + 1
             n_dists(3) = n_dists(3) + 1
             adf(binindex(2),binindex(1),2) = adf(binindex(2),binindex(1),2) + 1.d0
             adf(binindex(2),binindex(1),3) = adf(binindex(2),binindex(1),3) + 1.d0
          END IF
          
       END DO
    END IF

    DO ii = 1, 3 ! loop over x,y,z _and_ the 3 cases: mb, atom and all
       IF(n_dists(ii) > 0)THEN          
          SELECT CASE(control%adf_axis) ! get total area
          CASE(x_index)
             normer = cell(2)*cell(3)
          CASE(y_index)
             normer = cell(1)*cell(3)
          CASE(z_index)
             normer = cell(1)*cell(2)
          END SELECT
          SELECT CASE(ii) ! get the normalizing coefficients
          CASE(1) !mb
             normer = normer/( REAL(n_mbs,KIND=dp)*binwidth(1)*binwidth(2) )
          CASE(2) !at
             normer = normer/( REAL(n_ats,KIND=dp)*binwidth(1)*binwidth(2) )
          CASE(3) !any
             normer = normer/&
                  ( REAL(n_mbs+n_ats,KIND=dp)*binwidth(1)*binwidth(2) )
          END SELECT
          IF(control%adf_angle)THEN ! angle resolved adf
             ! record the adf
             DO jj = 1, rdfbins
                adf(:,jj,ii) = adf(:,jj,ii)*normer/((REAL(jj,KIND=dp)-0.5d0)*binwidth(1))
             END DO
          ELSE ! angle averaged adf
             DO jj = 1, rdfbins
                adf(1,jj,ii) = adf(1,jj,ii)*normer/(rdfbins*(REAL(jj,KIND=dp)-0.5d0)*binwidth(1))
             END DO
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE axial_distribution_function


!!$  ! Creates a custom MPI type for mb molecules (not used, but you could)
!!$  SUBROUTINE create_mb_type(new_mpi_type)
!!$    IMPLICIT NONE
!!$    TYPE(mb) :: typer
!!$    INTEGER :: err, blocks(10), displs(10), address(11), types(10), ii
!!$    INTEGER, INTENT(OUT) :: new_mpi_type
!!$
!!$    ! types
!!$    types(1) = MPI_qtrn
!!$    DO ii = 2, 8
!!$       types(ii) = MPI_DOUBLE_PRECISION
!!$    END DO
!!$    types(9) = MPI_INTEGER
!!$    types(10) = MPI_INTEGER
!!$    
!!$    ! block lengths
!!$    blocks(1) = 1
!!$    blocks(2) = 3
!!$    blocks(3) = 3
!!$    blocks(4) = 3
!!$    blocks(5) = 3
!!$    blocks(6) = 3
!!$    blocks(7) = 1
!!$    blocks(8) = 1
!!$    blocks(9) = 3
!!$    blocks(10) = 1
!!$    
!!$    ! addresses
!!$    CALL MPI_Address(typer,address(1),err)
!!$    CALL MPI_Address(typer%orientation,address(2),err)    
!!$    CALL MPI_Address(typer%pos,address(3),err)
!!$    CALL MPI_Address(typer%vel,address(4),err)
!!$    CALL MPI_Address(typer%angvel,address(5),err)
!!$    CALL MPI_Address(typer%well,address(6),err)
!!$    CALL MPI_Address(typer%inipos,address(7),err)
!!$    CALL MPI_Address(typer%m_tot,address(8),err)
!!$    CALL MPI_Address(typer%m_inert,address(9),err)
!!$    CALL MPI_Address(typer%constrained,address(10),err)
!!$    CALL MPI_Address(typer%index,address(11),err)
!!$    DO ii = 1, 10
!!$       displs(ii) = address(ii+1) - address(1)
!!$    END DO
!!$
!!$    ! build
!!$    CALL MPI_TYPE_STRUCT(11,blocks,displs,types,new_mpi_type,err)
!!$    CALL MPI_TYPE_COMMIT(new_mpi_type,err)
!!$
!!$  END SUBROUTINE create_mb_type
!!$
!!$
!!$  ! Creates a custom MPI type for atoms (not used, but you could)
!!$  SUBROUTINE create_atom_type(new_mpi_type)
!!$    IMPLICIT NONE
!!$    TYPE(atom) :: typer
!!$    INTEGER :: err, blocks(9), displs(9), address(10), types(9), ii
!!$    INTEGER, INTENT(OUT) :: new_mpi_type
!!$
!!$    ! types    
!!$    DO ii = 1, 5
!!$       types(ii) = MPI_DOUBLE_PRECISION
!!$    END DO
!!$    types(6) = MPI_INTEGER
!!$    types(7) = MPI_INTEGER
!!$    types(8) = MPI_INTEGER
!!$    types(9) = MPI_CHARACTER
!!$    
!!$    ! block lengths
!!$    blocks(1) = 3
!!$    blocks(2) = 3
!!$    blocks(3) = 3
!!$    blocks(4) = 3
!!$    blocks(5) = 1
!!$    blocks(6) = 3
!!$    blocks(7) = 1
!!$    blocks(8) = 1
!!$    blocks(9) = labelw
!!$    
!!$    ! addresses
!!$    CALL MPI_Address(typer,address(1),err)
!!$    CALL MPI_Address(typer%pos,address(2),err)
!!$    CALL MPI_Address(typer%vel,address(3),err)
!!$    CALL MPI_Address(typer%well,address(4),err)
!!$    CALL MPI_Address(typer%inipos,address(5),err)
!!$    CALL MPI_Address(typer%mass,address(6),err)
!!$    CALL MPI_Address(typer%constrained,address(7),err)
!!$    CALL MPI_Address(typer%type,address(8),err)
!!$    CALL MPI_Address(typer%index,address(9),err)
!!$    CALL MPI_Address(typer%element,address(10),err)
!!$    DO ii = 1, 9
!!$       displs(ii) = address(ii+1) - address(1)
!!$    END DO
!!$
!!$    ! build
!!$    CALL MPI_TYPE_STRUCT(9,blocks,displs,types,new_mpi_type,err)
!!$    CALL MPI_TYPE_COMMIT(new_mpi_type,err)
!!$
!!$  END SUBROUTINE create_atom_type

#ifdef MPI
  ! stacks the "lists" from all cpus together according to the lengths given in "items"
  ! and gathers the complete list to cpu 0
  ! For example:
  ! 
  ! cpu 0          cpu 1          cpu 0
  ! abc....        12.....        abc12..
  ! de.....        3456...   ->   de3456.
  ! fghij..        78.....        fghij78
  !
  ! The stacking is done for the first array index: list(:,1).
  ! The stacking works so that first every cpu 2n+1 sends its data to cpu 2n,
  ! then 2*(2n+1) sends data to 2*2n, and so on, until the final cpu 2^m sends its data to cpu 0:
  !
  ! cpu
  ! 0 1 2 3 4 5 6 7 8 9 10
  ! |-/ |-/ |-/ |-/ |-/ |
  ! |---/   |---/   |---/
  ! |-------/       |
  ! |---------------/
  ! x
  !
  ! *list 2d arrays containing lists to be stacked 
  ! *items the numbers of items to be stacked in each list
  ! *length the number of lists (size of list(1,:))
  ! *max size of lists (size of list(:,1))
  SUBROUTINE mpi_stack(list,items,length,width)
    IMPLICIT NONE
    INTEGER, POINTER :: list(:,:), items(:)
    INTEGER, INTENT(IN) :: length,width
    INTEGER :: remainder, templist(width,length),tempitems(length), ii,level
    LOGICAL :: fine

    fine = .false.
    level = 1 ! level of communications: cpu level*(2n+1) sends to cpu level*n

    DO WHILE(.not.fine) ! "fine" is a marker for ending the receiving-sending loop for this cpu

       level = 2*level ! advance to next level
       remainder = cpu_id - (cpu_id/level)*level ! the remainder for cpu_id/level
       
       ! level gets values 2^k, so the remainders develop as follows:
       ! level  cpu/remainder
       !        0 1 2 3 4 5 6 7 8 9
       ! 2      0 1 0 1 0 1 0 1 0 1
       ! 4      0 1 2 3 0 1 2 3 0 1
       ! 8      0 1 2 3 4 5 6 7 0 1
       ! 16     0 1 2 3 4 5 6 7 8 9
       !  ...
       ! So, the cou should receive or wait when the remainder is zero
       ! and send and finish once the remainder becomes positive

       IF(remainder == 0)THEN
          IF(n_cpus > cpu_id + level/2 )THEN ! is there a cpu that should be sending?
             ! receive the lists and the numbers of elements in the lists
             CALL MPI_RECV(templist,length*width,MPI_INTEGER,cpu_id+level/2,&
                  5501+level,MPI_COMM_WORLD,mpistatus,mpistat)
             CALL MPI_RECV(tempitems,length,MPI_INTEGER,cpu_id+level/2,&
                  5502+level,MPI_COMM_WORLD,mpistatus,mpistat)

             DO ii = 1, length ! stack the lists
                IF(tempitems(ii) > 0)THEN
                   list(items(ii)+1:items(ii)+tempitems(ii),ii) = templist(1:tempitems(ii),ii)
                   items(ii) = items(ii) + tempitems(ii)
                END IF
             END DO             
          ELSE ! there are no more cpus that would send data to this one
             IF(cpu_id /= 0)THEN
                ! since this is not the master cpu, it must wait until it is its turn to send its data forward
                ! so, skip this level and repeat the loop (will advance the level)
             ELSE
                fine = .true. ! the master cpu has received all the data, finish now
             END IF
          END IF
       ELSE ! the first time the remainder is more than zero, send the data and finish for this cpu
          CALL MPI_SEND(list,length*width,MPI_INTEGER,cpu_id-level/2,&
               5501+level,MPI_COMM_WORLD,mpistat)
          CALL MPI_SEND(items,length,MPI_INTEGER,cpu_id-level/2,&
               5502+level,MPI_COMM_WORLD,mpistat)
          fine = .true.
       END IF
    END DO

    RETURN

  END SUBROUTINE mpi_stack
#endif

  ! For two dangling bond arms, returns the quaternion representation of the orientation
  ! for which the first two arms are the given vectors.
  ! *hb1 a dangling bond vector
  ! *hb2 a dangling bond vector
  FUNCTION get_orientation(hb1,hb2) &
       RESULT(qq)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: hb1(3), hb2(3)
    TYPE(qtrn) :: qq
    REAL(KIND=dp) :: ub1(3), ub2(3), axis(3), angle, ri(3), p1(3), p2(3)
    TYPE(qtrn) :: q1, q2

    ! norm hbs (ub = unit bond)
    ub1 = hb1/(.norm.hb1)
    ub2 = hb2/(.norm.hb2)

    ! first rotate i1 to ub1
    axis = (i1 .x. ub1)
    angle = acos( (i1.o.ub1) )
    q1 = rot2q(angle,axis)

    ! i2 after the rotation (ri)
    ri = rotate(i2,q1)

    ! get projections of ri and ub2 perpendicular to ub1 (p1 and p2)
    p2 = ub2 - (ub2 .o. ub1)*ub1
    p1 = ri  - (ri  .o. ub1)*ub1

    ! rotate p1 to p2
    axis = (p1 .x. p2)
    angle = acos( (p1.o.p2)/((.norm.p1)*(.norm.p2)) )
    q2 = rot2q(angle,axis)

    ! get the combined rotation
    qq = q2*q1

    RETURN
  END FUNCTION get_orientation

END MODULE mb_model
