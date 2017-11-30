! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Implementation of the 32 bit Mersenne Twister Pseudo Random Number Generator
! based on the Mersenne Prime 2^19937-1

MODULE mt19937_mod

  ! Desciption:
  !   Implementation of the Mersenne Twister Pseudo Random Number Generator
  !   to compare Model Output of different compilers and computers.
  !
  ! Method:
  !   Provides two overloaded routines:
  !   * mt19937_random_seed   (same interface as intrinsic RANDOM_SEED) 1)
  !   * mt19937_random_number (same interface as intrinsic RANDOM_NUMBER)
  !   1) routine doesn't know the IBM Fortran "generator" keyword
  !
  !   The private routines mt19937_seed, mt19937_extract, and twist are 
  !   based on the Pseudo-Code of
  !   https://en.wikipedia.org/wiki/Mersenne_Twister#Pseudocode
  !
  !   The public subroutine mt19937_random_seed is a wrapper for mt19937_seed
  !   to ensure the same behaviour as intrinsic RANDOM_SEED
  !   The public subroutine mt19937_random_number is in fact an overloaded 
  !   interface, pointing to private subroutines for different possible harvest
  !   parameters (scalar up to 3D array of both real32 and real64)
  !
  ! Code Owner: 
  !   08/11/2017: Holger Wolff holger.wolff@monash.edu
  !
  ! Code Description:
  !   Language: Fortran 2003
  !   This code is close to UMDP3, except that more than one routine is needed
  !   for the random number generator.
  !
  ! Todo:
  !   UMDP3, Section 3.18: Prevention of Code Duplication

  ! iso_fortran_env int32 and real64 are Fortran 2008,
  ! but iso_c_binding is Fortran 2003 and allowed in UMDP3
  USE iso_c_binding, ONLY : int32 => c_int32_t, real32 => c_float, real64 => c_double
  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MT19937_MOD'

  ! The constants for mt19937_32
  INTEGER(int32), PARAMETER :: f = 1812433253
  INTEGER(int32), PARAMETER :: w = 32, n = 624, m = 397, r = 31
  INTEGER(int32), PARAMETER :: a = -1727483681
  INTEGER(int32), PARAMETER :: u = 11, d = -1
  INTEGER(int32), PARAMETER :: s = 7,  b = -1658038656
  INTEGER(int32), PARAMETER :: t = 15, c = -272236544
  INTEGER(int32), PARAMETER :: l = 18
  PRIVATE :: f, w, n, m, r, a, u, d, s, b, t, c, l

  INTEGER(int32) :: mt(0:n-1)
  INTEGER(int32) :: index = n + 1

  INTEGER(int32), PARAMETER :: lower_mask = huge(1_int32)   ! shiftl(1_int32, r) - 1
  INTEGER(int32), PARAMETER :: upper_mask = not(lower_mask)

  PRIVATE :: mt, index, lower_mask, upper_mask

  PRIVATE :: mt19937_seed, mt19937_extract
  PRIVATE :: twist

  INTERFACE mt19937_random_number
    MODULE PROCEDURE random_number_r64, random_number_r32, &
      random_number_1d_r64, random_number_1d_r32, &
      random_number_2d_r64, random_number_2d_r32, &
      random_number_3d_r64, random_number_3d_r32
  END INTERFACE mt19937_random_number
  PRIVATE :: random_number_r64, random_number_1d_r64
  PRIVATE :: random_number_2d_r64, random_number_3d_r64
  PUBLIC :: mt19937_random_number
  PUBLIC :: mt19937_random_seed

  INTEGER(kind=int32), PRIVATE :: my_seed(1)

  INTEGER, PARAMETER, PRIVATE :: maxint = huge(1_int32)
  INTEGER, PARAMETER :: i32_mask = int(huge(1_int32), kind=kind(1)) ! transfer(-1_int32, 1)

CONTAINS

  SUBROUTINE mt19937_seed(seed)
    ! 
    ! Description:
    !   Takes an int32 integer and uses it to seed 
    !   the Mersenne Twister
    !
    !   Currently PRIVATE, but may be exposed if needed.
    IMPLICIT NONE
    INTEGER(int32), INTENT(IN) :: seed
    INTEGER :: i
    index = n
    mt(0) = seed
    DO i = 1, n-1
      mt(i) = f * ieor(mt(i-1), (shiftr( mt(i-1), w-2 ))) + i
    END DO
  END SUBROUTINE mt19937_seed

  FUNCTION mt19937_extract() RESULT(y)
    !
    ! Description:
    !   Returns a random int32 number
    !   If the PRNG has not been seeded yet, it will seed it with
    !   the number 4395.
    !
    !   Currently PRIVATE, but may be exposed if needed.
    IMPLICIT NONE
    INTEGER(int32) :: y
    IF (index > n) THEN
      CALL mt19937_seed(4395_int32)
    END IF
    IF (index == n) THEN
      CALL twist()
    END IF

    y = mt(index)
    y = ieor(y, iand(d, shiftr(y, u)))
    y = ieor(y, iand(b, shiftl(y, s)))
    y = ieor(y, iand(c, shiftl(y, t)))
    y = ieor(y, shiftr(y, l))

    index = index + 1
  END FUNCTION mt19937_extract

  SUBROUTINE twist()
    !
    ! Description:
    !   Creates a new array of numbers on which to base
    !   the next random numbers on.
    !
    !   Currently PRIVATE, no need to expose.
    IMPLICIT NONE
    INTEGER :: i
    INTEGER(int32)  :: x, xA
    DO i = 0, n-1
      x = iand(upper_mask, mt(i)) + &
      iand(lower_mask, mt(mod(i+1, n)))
      xA = shiftr(x, 1)
      IF (modulo(x, 2) /= 0) THEN
        xA = ieor(xA, a)
      END IF
      mt(i) = ieor(xA, mt(modulo(i+m, n)))
    END DO
    index = 0
  END SUBROUTINE twist

  SUBROUTINE mt19937_random_seed(put, get, size)
    !
    ! Description:
    !   Wrapper for currently private mt19937_seed to
    !   give it the same interface as intrinsic RANDOM_SEED
    !   (Only Fortran Standard, no IBM-specific 'generator'
    !   argument.)
    IMPLICIT NONE
    INTEGER, INTENT(IN), OPTIONAL :: put(:)
    INTEGER, INTENT(OUT), OPTIONAL :: get(:)
    INTEGER, INTENT(OUT), OPTIONAL :: size
    LOGICAL :: anything_selected
    REAL :: r
    anything_selected = .FALSE.
    IF (present(size)) THEN
      size = 1
      anything_selected = .TRUE.
    END IF
    IF (present(put)) THEN
      my_seed(1) = int(iand(put(1), i32_mask), kind=int32)
      CALL mt19937_seed(my_seed(1))
      anything_selected = .TRUE.
    END IF
    IF (present(get)) THEN
      get = 0
      get(1) = my_seed(1)
      anything_selected = .TRUE.
    END IF
    IF (.not. anything_selected) THEN
      ! Default behaviour of RANDOM_SEED is to seed
      ! the PRNG with the time.
      CALL RANDOM_SEED
      CALL RANDOM_NUMBER(r)
      my_seed(1) = iand(transfer(r, 1), -1)
      CALL mt19937_seed(my_seed(1))
    END IF
  END SUBROUTINE mt19937_random_seed

    SUBROUTINE random_number_r32(harvest)
      !
      ! Description:
      !   Wrapper to turn random integer number
      !   from mt19937_extract into a real32 
      !   with 0 <= harvest <= 1
      IMPLICIT NONE
      REAL(real32), INTENT(OUT) :: harvest
      INTEGER(int32) :: r, m
      r = mt19937_extract()
      m = iand(r, maxint)
      harvest = real(m, kind=real32) / maxint
    END SUBROUTINE random_number_r32

    SUBROUTINE random_number_1d_r32(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 1d array of real32
      REAL(real32), INTENT(OUT) :: harvest(:)
      INTEGER :: i
      INTEGER(int32) :: r, m
      do i = lbound(harvest, 1), ubound(harvest, 1)
        call random_number_r32(harvest(i))
      end do
    END SUBROUTINE random_number_1d_r32

    SUBROUTINE random_number_2d_r32(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 2d array of real32
      REAL(real32), INTENT(OUT) :: harvest(:,:)
      INTEGER :: i
      DO i = lbound(harvest, 2), ubound(harvest, 2)
        CALL random_number_1d_r32(harvest(:, i))
      END DO
    END SUBROUTINE random_number_2d_r32

    SUBROUTINE random_number_3d_r32(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 3d array of real32
      REAL(real32), INTENT(OUT) :: harvest(:,:,:)
      INTEGER :: i
      DO i = lbound(harvest, 3), ubound(harvest, 3)
        CALL random_number_2d_r32(harvest(:, :, i))
      END DO
    END SUBROUTINE random_number_3d_r32

    SUBROUTINE random_number_r64(harvest)
      !
      ! Description:
      !   Wrapper to turn random integer number
      !   from mt19937_extract into a real64 
      !   with 0 <= harvest <= 1
      IMPLICIT NONE
      REAL(real64), INTENT(OUT) :: harvest
      INTEGER(int32) :: r, m
      r = mt19937_extract()
      m = iand(r, maxint)
      harvest = real(m, kind=real64) / maxint
    END SUBROUTINE random_number_r64

    SUBROUTINE random_number_1d_r64(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 1d array of real64
      REAL(real64), INTENT(OUT) :: harvest(:)
      INTEGER :: i
      INTEGER(int32) :: r, m
      do i = lbound(harvest, 1), ubound(harvest, 1)
        call random_number_r64(harvest(i))
      end do
    END SUBROUTINE random_number_1d_r64

    SUBROUTINE random_number_2d_r64(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 2d array of real64
      REAL(real64), INTENT(OUT) :: harvest(:,:)
      INTEGER :: i
      DO i = lbound(harvest, 2), ubound(harvest, 2)
        CALL random_number_1d_r64(harvest(:, i))
      END DO
    END SUBROUTINE random_number_2d_r64

    SUBROUTINE random_number_3d_r64(harvest)
      !
      ! Description:
      !   Wrapper if harvest is a 3d array of real64
      REAL(real64), INTENT(OUT) :: harvest(:,:,:)
      INTEGER :: i
      DO i = lbound(harvest, 3), ubound(harvest, 3)
        CALL random_number_2d_r64(harvest(:, :, i))
      END DO
    END SUBROUTINE random_number_3d_r64

END MODULE mt19937_mod
