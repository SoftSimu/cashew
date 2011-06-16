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
! functions contains some general functions and subroutines.
! mostly, these are more elaborate random generators, such as
! the random unit vector generator "rand_unit_vector", 
! based on Mersenne Twister. They are
! here since mt95 is not mine and so I do not want to alter it.
! some functions like "gauss" are here just so one can easily try more efficient
! methods for interpolating or evaluating the functions simply by
! changin their definitions here.
! <br /><br />Back to <a href="cashew.html">cashew</a>
!
MODULE functions

  USE parameters
  USE quaternions
  USE mt95
  USE mpi_mod

  IMPLICIT NONE

CONTAINS

  ! Evaluates the non-normalized gaussian exp(-xval**2/2)
  ! *xval x value
  ! *gaussian exp(-xval**2/2)
  FUNCTION gauss(xval) &
       RESULT(gf)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: xval
    REAL(KIND=dp) :: gf

    IF(ABS(xval) > gaussian_limit)THEN
       gf = 0.d0
    ELSE
       gf = exp(-0.5*xval*xval)
    END IF
    RETURN
  END FUNCTION gauss

  ! Generates a random unit quaternion
  ! *random_q random unit quaternion
  SUBROUTINE rand_unit_quaternion(random_q)
    IMPLICIT NONE
    REAL(KIND=dp) :: random_vector(3), rangle
    TYPE(qtrn), INTENT(OUT) :: random_q

    CALL rand_unit_vector(random_vector)
    CALL genrand_real2(rangle)
    rangle = rangle*2.d0*pi
    random_q = rot2q(rangle,random_vector)

    RETURN
  END SUBROUTINE rand_unit_quaternion

  ! Returns a unit vector randomly chosen on a sphere
  ! *random_vector a random unit vector
  SUBROUTINE rand_unit_vector(random_vector)
    IMPLICIT NONE
    REAL(KIND=dp) :: z, phi, theta
    REAL(KIND=dp), INTENT(OUT) :: random_vector(3)

    CALL genrand_real2(phi)
    phi = phi*2*pi
    CALL genrand_real1(z)
    z = 2.d0*(z-0.5d0)

    theta = asin(z)

    random_vector(1) = cos(theta)*cos(phi)
    random_vector(2) = cos(theta)*sin(phi)
    random_vector(3) = z

  END SUBROUTINE rand_unit_vector
  
  ! Generates two normally distributed random numbers
  ! with sigma = 1 using the polar algorithm
  ! *r1 N(0,1) distributed random number
  ! *r2 N(0,1) distributed random number
  SUBROUTINE rand_2gaussians(r1,r2)
    IMPLICIT NONE
    REAL(KIND=dp) :: x1, x2, r1, r2, rr

    rr = 2.d0
    DO WHILE( rr >= 1.d0)
       CALL genrand_real1(x1)
       CALL genrand_real1(x2)
       x1 = 2.d0*x1 - 1.d0
       x2 = 2.d0*x2 - 1.d0
       rr = x1*x1 + x2*x2
    END DO

    rr = SQRT( (-2.d0*LOG(rr))/rr )
    r1 = x1*rr
    r2 = x2*rr

    RETURN

  END SUBROUTINE rand_2gaussians

  ! Generate two vectors according to the maxwell-boltzmann
  ! distribution, i.e., the components are
  ! normally distributed with sigma = 1.
  ! *bv1 random vector with N(0,1) components
  ! *bv2 random vector with N(0,1) components
  SUBROUTINE rand_2boltzmann_vectors(bv1,bv2)
    IMPLICIT NONE
    REAL(KIND=dp) :: g1, g2, g3, g4, g5, g6
    REAL(KIND=dp), INTENT(OUT) :: bv1(3), bv2(3)

    CALL rand_2gaussians(g1,g2)
    CALL rand_2gaussians(g3,g4)
    CALL rand_2gaussians(g5,g6)

    bv1(1) = g1
    bv1(2) = g2
    bv1(3) = g3
    bv2(1) = g4
    bv2(2) = g5
    bv2(3) = g6

    RETURN

  END SUBROUTINE rand_2boltzmann_vectors



  ! Returns the component of vector u projected on the plane
  ! perpendicular to unit vector r
  ! *uu vector to be projected
  ! *rr normal vector of the plane of projection
  ! *proji the projection
  FUNCTION unit_projection(uu,rr) &
       RESULT(proji)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: uu(3), rr(3)
    REAL(KIND=dp) :: proji(3)

    proji = uu - (uu.o.rr)*rr

  END FUNCTION unit_projection

  ! transfer a coordinate inside the range [0,celllength]
  ! by adding or subtracting a suitable multiple of celllength.
  ! the purpose of this is to get the equivalent coordinates
  ! of an arbitrary point inside a periodic supercell.
  ! *coordinate the coordinate to be manipulated
  ! *celllength the length of the cell (range of the value given to "coordinate")
  SUBROUTINE move_in_cell(coordinate,celllength)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: coordinate
    REAL(KIND=dp), INTENT(IN) :: celllength

    DO WHILE(coordinate > celllength)
       coordinate = coordinate-celllength
    END DO
    DO WHILE(coordinate < 0.d0)
       coordinate = coordinate+celllength
    END DO
    RETURN

  END SUBROUTINE move_in_cell


  ! quit execution and print a message
  ! *message (part of) the message to be printed, should tell the user where the error occurred
  SUBROUTINE abort(message)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: message

    !IF(cpu_id == master_cpu)THEN
       WRITE(*,*) "Error encountered: "//message
       WRITE(*,*) "  aborting..."
    !END IF
    STOP 1

  END SUBROUTINE abort

  ! formats a string to display days, hours, minutes, seconds and milliseconds
  ! *days days
  ! *hour hours
  ! *min minutes
  ! *sec seconds
  ! *millisec milliseconds
  ! *thetime a string containing the nicely formatted time
  SUBROUTINE time_string(days,hour,min,sec,millisec,thetime)
    IMPLICIT NONE
    CHARACTER(LEN=18), INTENT(OUT) :: thetime
    INTEGER, INTENT(IN) :: days,hour,min,sec,millisec
    CHARACTER(LEN=2) :: adding
    CHARACTER(LEN=3) :: mstime
    CHARACTER(LEN=4) :: daytime

    thetime = ""
    adding = ""
    daytime = ""

    IF(days >= 1000)THEN
       thetime = "a long time..."
       RETURN
    END IF

    IF(days > 0)THEN
       WRITE(daytime,'(I3)') days
       IF(days < 100)THEN
          daytime(1:1) = "0"
       END IF
       IF(days < 10)THEN
          daytime(2:2) = "0"
       END IF
       thetime(1:6) = daytime//" d "
    ELSE
       thetime(1:6) = "      " 
    END IF

    WRITE(adding,'(I2)') hour
    IF(hour < 10)THEN
       adding(1:1) = "0"
    END IF
    thetime(7:9) = adding//":"

    WRITE(adding,'(I2)') min
    IF(min < 10)THEN
       adding(1:1) = "0"
    END IF
    thetime(10:12) = adding//":"

    WRITE(adding,'(I2)') sec
    IF(sec < 10)THEN
       adding(1:1) = "0"
    END IF
    thetime(13:15) = adding//"."

    WRITE(mstime,'(I3)') millisec
    IF(millisec < 100)THEN
       mstime(1:1) = "0"
    END IF
    IF(millisec < 10)THEN
       mstime(2:2) = "0"
    END IF
    thetime(16:18) = mstime

    RETURN
  END SUBROUTINE time_string

  ! Formats the date and time as dd.mm.yyyy  hh:mm:ss
  ! *year year
  ! *month month
  ! *date day of month
  ! *hour hour
  ! *min minute
  ! *sec second
  ! *thedate a string containing the nicely formatted time
  SUBROUTINE date_string(year,month,date,hour,min,sec,thedate)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(OUT) :: thedate
    INTEGER, INTENT(IN) :: year,month,date,hour,min,sec
    CHARACTER(LEN=2) :: adding
    CHARACTER(LEN=4) :: yearstr
    
    thedate = ""
    adding = ""
    WRITE(adding,'(I2)') date
    IF(date < 10)THEN
       adding(1:1) = "0"
    END IF
    thedate(1:3) = adding//"."
    WRITE(adding,'(I2)') month
    IF(month < 10)THEN
       adding(1:1) = "0"
    END IF
    thedate(4:6) = adding//"."
    WRITE(yearstr,'(I4)') year
    thedate(7:12) = yearstr//"  "
    WRITE(adding,'(I2)') hour
    IF(hour < 10)THEN
       adding(1:1) = "0"
    END IF
    thedate(13:15) = adding//":"
    WRITE(adding,'(I2)') min
    IF(min < 10)THEN
       adding(1:1) = "0"
    END IF
    thedate(16:18) = adding//":"
    WRITE(adding,'(I2)') sec
    IF(sec < 10)THEN
       adding(1:1) = "0"
    END IF
    thedate(19:20) = adding

    RETURN
  END SUBROUTINE date_string

  ! Takes the (rounded) integer part of a real
  ! and formats it to a string
  ! *realnmb the real to be formatted
  ! *str the resulting string
  FUNCTION real_to_string(realnmb) &
       RESULT(str)
    IMPLICIT NONE
    CHARACTER(LEN=10) :: str
    REAL(KIND=dp), INTENT(IN) :: realnmb
    INTEGER :: ii

    WRITE(str,'(I10)') INT(realnmb+0.45)

    DO ii = 1, 10
       IF(str(ii:ii) == " ") str(ii:ii) = "0"
    END DO

    RETURN
  END FUNCTION real_to_string

  ! time difference calculator for wall clock:
  ! calculates the difference "time2"-"time1", where the
  ! times are as given by "date_and_time", and returns it
  ! in the vector "timeout" as days, hours, minutes, seconds and
  ! milliseconds. However, as it is now, the routine doesn't handle months.
  ! The last integer in timeout is a check for success: if it is
  ! 0, all is fine - if 1, the result is wrong.
  ! This should be corrected at some point...
  ! *time1 former time as given by "date_and_time"
  ! *time2 latter time as given by "date_and_time"
  ! *timeout the time difference in days, hours, mins, secs, msecs, check
  SUBROUTINE time_difference(time1,time2,timeout)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: time1(8), time2(8)
    INTEGER, INTENT(OUT) :: timeout(6)
    INTEGER :: memo
    
    ! ms
    IF(time2(8) >= time1(8))THEN
       timeout(5) = time2(8)-time1(8)
       memo = 0
    ELSE
       timeout(5) = time2(8)-time1(8)+1000
       memo = -1
    END IF
    ! s
    IF(time2(7)+memo >= time1(7))THEN
       timeout(4) = time2(7)-time1(7)+memo
       memo = 0
    ELSE
       timeout(4) = time2(7)-time1(7)+60+memo
       memo = -1
    END IF
    ! min
    IF(time2(6)+memo >= time1(6))THEN
       timeout(3) = time2(6)-time1(6)+memo
       memo = 0
    ELSE
       timeout(3) = time2(6)-time1(6)+60+memo
       memo = -1
    END IF
    ! h
    IF(time2(5)+memo >= time1(5))THEN
       timeout(2) = time2(5)-time1(5)+memo
       memo = 0
    ELSE
       timeout(2) = time2(5)-time1(5)+24+memo
       memo = -1
    END IF
    ! d
    IF(time2(3)+memo >= time1(3))THEN
       timeout(1) = time2(3)-time1(3)+memo
       memo = 0
    ELSE
       timeout(1) = time2(3)-time1(3)+24+memo
       memo = -1
    END IF
    timeout(6) = memo

    RETURN
  END SUBROUTINE time_difference

  ! takes the cpu_time reading and adds to the given timer:
  ! tt2 = tt1, tt1 = current time, stopwatch += tt1-tt2
  ! *stopwatch the timer for total time, tt1-tt2 will be added to it
  ! *tt1 should contain the previous time, the current time will be stored here
  ! *tt2 the previous time will be moved here
  SUBROUTINE checkpoint(stopwatch,tt1,tt2)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: stopwatch,tt1,tt2

    tt2 = tt1
    CALL cpu_time(tt1)
    stopwatch = stopwatch + tt1-tt2

  END SUBROUTINE checkpoint

END MODULE functions
