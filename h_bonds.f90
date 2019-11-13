!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM h_bonds
USE OMP_LIB
USE INPUT

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), OHvec_mat(:,:,:), is_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
INTEGER                         :: nb_line_is, nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: OpH_disp_vec(3), OpH_disp_norm
REAL(dp)                        :: oHpO_disp_vec(3), oHpO_disp_norm
REAL(dp)                        :: XpOh_disp_vec(3), XpOh_disp_norm
REAL(dp)                        :: XpO_disp_vec(3), XpO_disp_norm
REAL(dp)                        :: SpOh_disp_vec(3), SpOh_disp_norm
REAL(dp)                        :: SpO_disp_vec(3), SpO_disp_norm
CHARACTER(LEN=2)                :: dummy_char

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, Udonnor_count, Udonnor_count2
INTEGER, ALLOCATABLE            :: nb_max_OHvec(:)

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o, l
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching H-Bonds'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*,"No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM(input_file)
CALL READINPUTSUB(input_file)
file_pos=TRIM(file_pos)
file_is=TRIM(file_is)

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, H-Bonds, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(23,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_el(nb_atm))

OPEN(UNIT=20,FILE=file_pos,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i=1,nb_atm
        atm_mat(2,i,s) = -1
        READ(20, *) atm_el(i), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
        atm_mat(1,i,s) = i
        IF (atm_el(i) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
            atm_mat(3,i,s) = 1
        ELSE IF (atm_el(i) .EQ. "Si") THEN
            atm_mat(2,i,s) = 28
            atm_mat(3,i,s) = 2
        ELSE IF (atm_el(i) .EQ. "SiF") THEN
            atm_mat(2,i,s) = 28
            atm_mat(3,i,s) = 1
        ELSE IF (atm_el(i) .EQ. "OE") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 10
        ELSE IF (atm_el(i) .EQ. "OH") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 11
        ELSE IF (atm_el(i) .EQ. "OA") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 12
        ELSE IF (atm_el(i) .EQ. "OW") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 13
        ELSE IF (atm_el(i) .EQ. "OM") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 14
        ELSE IF (atm_el(i) .EQ. "OP") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 15
        ELSE IF (atm_el(i) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = -1
        ELSE IF (atm_el(i) .EQ. "HW") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 23
        ELSE IF (atm_el(i) .EQ. "HO") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 21
        ELSE IF (atm_el(i) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = -1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN
    ! A ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
    start = OMP_get_wtime()

    OPEN(UNIT=21, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
    nb_line_is = 0
    DO
        READ(21, *, IOSTAT=iostatus)
        IF (iostatus .NE. 0) THEN
            EXIT
        ELSE
            nb_line_is = nb_line_is + 1
        END IF
    END DO
    REWIND(21)
    nb_max_is = CEILING(1.0 * nb_line_is / nb_step) * 2

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

    ! A ----------------------------------------------- Read IS
    start = OMP_get_wtime()

    ALLOCATE(is_mat(5,nb_max_is,nb_step))
    ALLOCATE(nb_is(nb_step))
    is_mat(:,:,:) = 0.0_dp
    nb_is(:) = 0

    OPEN(UNIT=21, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
    DO s = 1, nb_step
        READ(21, *) nb_is(s)
        READ(21, *)
        j = 0
        DO i=1,nb_is(s)
            READ(21, *) dummy_char, is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
            j = j + 1
            is_mat(5,i,s) = j
            IF (is_mat(3,i,s) .LT. 10.0) THEN
                is_mat(4,i,s) = 1
            ELSE
                is_mat(4,i,s) = 2
            END IF
            DO k = 1, 3
                is_mat(k,i,s) = is_mat(k,i,s) - box(k) * ANINT(is_mat(k,i,s)/box(k))
            END DO
        END DO
    END DO
    CLOSE(UNIT=21)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
END IF

! B ----------------------------------------------- OH groups and corresponding values
start = OMP_get_wtime()
ALLOCATE(OHvec_mat(33,nb_o*3,nb_step))
ALLOCATE(nb_max_OHvec(nb_step))

nb_max_OHvec(:) = 0.0
OHvec_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_max_OHvec, nb_step, nb_atm)&
!$OMP SHARED(hb_OpH_rcut)&
!$OMP PRIVATE(s, i, j, k, o)&
!$OMP PRIVATE(OpH_disp_vec, OpH_disp_norm)
DO s = 1, nb_step
    o = 0
    DO i=1,nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            DO j=1,nb_atm
                IF (atm_mat(2,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        OpH_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        OpH_disp_vec(k) = OpH_disp_vec(k) - box(k) * ANINT(OpH_disp_vec(k)/box(k))
                    END DO
                    OpH_disp_norm = NORM2(OpH_disp_vec)
                    IF(OpH_disp_norm .LT. hb_OpH_rcut) THEN
                        o = o + 1
                        OHvec_mat(1,o,s) = atm_mat(1,i,s)
                        OHvec_mat(2,o,s) = atm_mat(1,j,s)
                        OHvec_mat(3,o,s) = atm_mat(3,i,s)
                        OHvec_mat(4,o,s) = atm_mat(3,j,s)
                        DO k = 1, 3
                            OHvec_mat(k+4,o,s) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                            OHvec_mat(k+4,o,s) = OHvec_mat(k+4,o,s) - box(k) * ANINT(OHvec_mat(k+4,o,s)/box(k)) ! Disp
                            OHvec_mat(k+7,o,s) = atm_mat(k+3,i,s) ! O pos
                            OHvec_mat(k+10,o,s) = atm_mat(k+3,j,s) ! H pos
                        END DO
                    END IF
                END IF
            END DO
        END IF
    END DO
    nb_max_OHvec(s) = COUNT(OHvec_mat(1,:,s) .NE. 0, DIM=1)
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH groups:", finish-start, "seconds elapsed"

! C ----------------------------------------------- OH/O Hbonds
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_atm, nb_step)&
!$OMP SHARED(hb_oHpO_rcut)&
!$OMP PRIVATE(s, i, j, k, l)&
!$OMP PRIVATE(oHpO_disp_vec, oHpO_disp_norm, Udonnor_count, Udonnor_count2)
DO s = 1, nb_step
 C1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE C1
        END IF
        Udonnor_count = 0
        Udonnor_count2 = 0
        DO j = 1, nb_atm
            IF (atm_mat(2,j,s) .EQ. 16) THEN
                DO k = 1, 3
                    oHpO_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+10,i,s)
                    oHpO_disp_vec(k) = oHpO_disp_vec(k) - box(k) * ANINT(oHpO_disp_vec(k)/box(k))
                END DO
                oHpO_disp_norm = NORM2(oHpO_disp_vec)
                IF( (oHpO_disp_norm .LE. hb_oHpO_rcut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                    atm_mat(7,j,s) = atm_mat(7,j,s) + 1 ! Acceptor count (O)
                    OHvec_mat(15,i,s) = OHvec_mat(15,i,s) + 1 ! Donnor count (OH)
                    IF (Udonnor_count .EQ. 0) THEN
                        OHvec_mat(16,i,s) = OHvec_mat(16,i,s) + 1 ! Unique donnor count (OH)
                        Udonnor_count = 1
                    END IF
                C21:DO l=1, nb_atm
                        IF (atm_mat(1,l,s) .EQ. OHvec_mat(1,i,s)) THEN
                            atm_mat(8,l,s) = atm_mat(8,l,s) + 1 ! Donnor count (O)
                            IF (Udonnor_count2 .EQ. 0) THEN
                                atm_mat(9,l,s) = atm_mat(9,l,s) + 1 ! Unique donnor count (O)
                                Udonnor_count2 = 1
                            END IF
                            EXIT C21
                        END IF
                    END DO C21
                C22:DO l=1,nb_o*3
                        IF (OHvec_mat(1,l,s) .EQ. 0) THEN
                            CYCLE C22
                        ELSE IF (OHvec_mat(1,l,s) .EQ. atm_mat(1,j,s)) THEN
                            OHvec_mat(14,l,s) = OHvec_mat(14,l,s) + 1 ! Acceptor count (OH)
                        END IF
                    END DO C22
                END IF
            END IF
        END DO
    END DO C1
  C3:DO i = 1, (nb_o*3)
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE C3
        END IF
        DO j = 1, nb_atm
            IF (OHvec_mat(1,i,s) .EQ. atm_mat(1,j,s)) THEN
                OHvec_mat(17,i,s) = atm_mat(7,j,s)
                OHvec_mat(18,i,s) = atm_mat(8,j,s)
                OHvec_mat(19,i,s) = atm_mat(9,j,s)
            END IF
        END DO
    END DO C3
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH/O Hbonds:", finish-start, "seconds elapsed"

! D ----------------------------------------------- Proximity between functionnal groups and any OH group
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_step, nb_atm)&
!$OMP SHARED(hb_XpOh_rcut, hb_CpOh_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XpOh_disp_vec, XpOh_disp_norm)
DO s = 1, nb_step
 D1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE D1
        END IF
     D2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE D2
            END IF
            DO k = 1, 3
                XpOh_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+7,i,s)
                XpOh_disp_vec(k) = XpOh_disp_vec(k) - box(k) * ANINT(XpOh_disp_vec(k)/box(k))
            END DO
            XpOh_disp_norm = NORM2(XpOh_disp_vec)
            IF( (XpOh_disp_norm .LE. hb_XpOh_rcut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(20,i,s) = 1
                    OHvec_mat(24,i,s) = 1
                    IF ( (XpOh_disp_norm .LE. OHvec_mat(32,i,s)) .AND.&
                    (OHvec_mat(32,i,s) .NE. 0) ) THEN
                        atm_mat(1,j,s) = OHvec_mat(31,i,s)
                        XpOh_disp_norm = OHvec_mat(32,i,s)
                        atm_mat(6,j,s) = OHvec_mat(33,i,s)
                    END IF
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    OHvec_mat(21,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    OHvec_mat(22,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    OHvec_mat(23,i,s) = 1
                END IF
            ELSE IF( (XpOh_disp_norm .LE. hb_CpOh_rcut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(24,i,s) = 1
                END IF
            END IF
        END DO D2
    END DO D1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and OH group:", finish-start, "seconds elapsed"

! E ----------------------------------------------- Proximity between functionnal groups and any O atom
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(hb_XpO_rcut, hb_CpO_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XpO_disp_vec, XpO_disp_norm)
DO s = 1, nb_step
    E1:DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .NE. 16) THEN
            CYCLE E1
        END IF
        E2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE E2
            END IF
            DO k = 1, 3
                XpO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                XpO_disp_vec(k) = XpO_disp_vec(k) - box(k) * ANINT(XpO_disp_vec(k)/box(k))
            END DO
            XpO_disp_norm = NORM2(XpO_disp_vec)
            IF ( (XpO_disp_norm .LE. hb_XpO_rcut) .AND. (atm_mat(1,j,s) .NE. atm_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1) THEN ! C
                    atm_mat(10,i,s) = 1
                    atm_mat(14,i,s) = 1
                    IF ( (XpO_disp_norm .LE. atm_mat(22,i,s)) .AND.&
                    (atm_mat(22,i,s) .NE. 0) ) THEN
                        atm_mat(1,j,s) = atm_mat(21,i,s)
                        XpO_disp_norm = atm_mat(22,i,s)
                        atm_mat(6,j,s) = atm_mat(23,i,s)
                    END IF
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    atm_mat(11,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    atm_mat(12,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    atm_mat(13,i,s) = 1
                END IF
            ELSE IF ( (XpO_disp_norm .LE. hb_CpO_rcut) .AND. (atm_mat(1,j,s) .NE. atm_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1) THEN ! C
                    atm_mat(14,i,s) = 1
                END IF
            END IF
        END DO E2
    END DO E1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and O atoms:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN

    ! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, is_mat, nb_is, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpOh_disp_vec, SpOh_disp_norm)
    DO s = 1, nb_step
        F2:DO i = 1, nb_o*3
            IF (OHvec_mat(1,i,s) .EQ. 0) THEN
                CYCLE F2
            END IF
            DO j = 1, nb_is(s)
                IF (is_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        SpOh_disp_vec(k) = is_mat(k,j,s) - OHvec_mat(k+7,i,s)
                        SpOh_disp_vec(k) = SpOh_disp_vec(k) - box(k) * ANINT(SpOh_disp_vec(k)/box(k))
                    END DO
                    SpOh_disp_norm = NORM2(SpOh_disp_vec)
                    IF ( (SpOh_disp_norm .LT. OHvec_mat(25,i,s)) .OR. (OHvec_mat(25,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(25,i,s) = SpOh_disp_norm
                        IF (OHvec_mat(10,i,s) .LT. is_mat(3,j,s)) THEN
                            OHvec_mat(26,i,s) = -1
                        ELSE
                            OHvec_mat(26,i,s) = 1
                        END IF
                        OHvec_mat(27,i,s) = is_mat(5,j,s)

                    END IF
                ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        SpOh_disp_vec(k) = is_mat(k,j,s) - OHvec_mat(k+7,i,s)
                        SpOh_disp_vec(k) = SpOh_disp_vec(k) - box(k) * ANINT(SpOh_disp_vec(k)/box(k))
                    END DO
                    SpOh_disp_norm = NORM2(SpOh_disp_vec)
                    IF ( (SpOh_disp_norm .LT. OHvec_mat(28,i,s)) .OR. (OHvec_mat(28,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(28,i,s) = SpOh_disp_norm
                        IF (OHvec_mat(10,i,s) .GT. is_mat(3,j,s)) THEN
                            OHvec_mat(29,i,s) = -1
                        ELSE
                            OHvec_mat(29,i,s) = 1
                        END IF
                        OHvec_mat(30,i,s) = is_mat(5,j,s)
                    END IF
                END IF
            END DO
        END DO F2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity IS and OH groups:", finish-start, "seconds elapsed"

    ! X ----------------------------------------------- Calculate closest distance between IS and any O atom
    start = OMP_get_wtime()

    !nb_step, nb_atm, always shared.
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_o, is_mat, nb_is, nb_atm, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpO_disp_vec, SpO_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF (atm_mat(2,i,s) .EQ. 16) THEN
                DO j = 1, nb_is(s)
                    IF (is_mat(4,j,s) .EQ. 1) THEN
                        DO k = 1, 3
                            SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                            SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT(SpO_disp_vec(k)/box(k))
                        END DO
                        SpO_disp_norm = NORM2(SpO_disp_vec)
                        IF ( (SpO_disp_norm .LT. atm_mat(15,i,s)) .OR. atm_mat(15,i,s) .EQ. 0.0 ) THEN
                            atm_mat(15,i,s) = SpO_disp_norm
                            IF (atm_mat(6,i,s) .LT. is_mat(3,j,s)) THEN
                                atm_mat(16,i,s) = -1
                            ELSE
                                atm_mat(16,i,s) = 1
                            END IF
                            atm_mat(17,i,s) = is_mat(5,j,s)
                        END IF
                    ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                        DO k = 1, 3
                            SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                            SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT(SpO_disp_vec(k)/box(k))
                        END DO
                        SpO_disp_norm = NORM2(SpO_disp_vec)
                        IF ( (SpO_disp_norm .LT. atm_mat(18,i,s)) .OR. atm_mat(18,i,s) .EQ. 0.0 ) THEN
                            atm_mat(18,i,s) = SpO_disp_norm
                            IF (atm_mat(6,i,s) .GT. is_mat(3,j,s)) THEN
                                atm_mat(19,i,s) = -1
                            ELSE
                                atm_mat(19,i,s) = 1
                            END IF
                            atm_mat(20,i,s) = is_mat(5,j,s)
                        END IF
                    END IF
                END DO
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS and O:", finish-start, "seconds elapsed"

END IF

! X ----------------------------------------------- Write OH BOND
start = OMP_get_wtime()

OPEN(UNIT=31, FILE=suffix//"_O_hbonds.txt")
WRITE(31, '(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A24,A24)')&
 "O_id", "Step","O_type", "Acc", "Don", "UDon"&
, "C", "OE", "OH", "OA", "C9", "dist_IS_down", "dist_IS_up"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            WRITE(31,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E24.14,E24.14)')&
            INT(atm_mat(1,i,s)), s, INT(atm_mat(3,i,s)), INT(atm_mat(7,i,s))&
            , INT(atm_mat(8,i,s)), INT(atm_mat(9,i,s)), INT(atm_mat(10,i,s)), INT(atm_mat(11,i,s))&
            , INT(atm_mat(12,i,s)), INT(atm_mat(13,i,s)), INT(atm_mat(14,i,s))&
            , (atm_mat(15,i,s)*atm_mat(16,i,s)),(atm_mat(18,i,s)*atm_mat(19,i,s))
    END IF
    END DO
END DO
CLOSE(UNIT=31)

OPEN(UNIT=32, FILE = suffix//"_OH_hbonds.txt")
WRITE(32, '(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A24,A24)')&
    "O_id", "O_type", "H_id", "step"&
    , "TAcc", "Don", "UDon"&
    , "C", "OE", "OH", "OA", "C9"&
    , "TAcc", "TDon", "TUDon", "dist_IS_down", "dist_IS_up"
DO s = 1, nb_step
    DO i = 1, nb_max_OHvec(s)
        WRITE(32,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E24.14,E24.14)')&
        INT(OHvec_mat(1,i,s)), INT(OHvec_mat(3,i,s)), INT(OHvec_mat(2,i,s)), s&
        , INT(OHvec_mat(14,i,s)) , INT(OHvec_mat(15,i,s)), INT(OHvec_mat(16,i,s))&
        , INT(OHvec_mat(20,i,s)), INT(OHvec_mat(21,i,s)), INT(OHvec_mat(22,i,s)), INT(OHvec_mat(23,i,s)) , INT(OHvec_mat(24,i,s))&
        , INT(OHvec_mat(17,i,s)), INT(OHvec_mat(18,i,s)), INT(OHvec_mat(19,i,s))&
        , (OHvec_mat(25,i,s)*OHvec_mat(26,i,s)), (OHvec_mat(28,i,s)*OHvec_mat(29,i,s))
    END DO
END DO
CLOSE(UNIT=32)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "O/OH Hbonds output:" ,finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat, nb_is)

DEALLOCATE(OHvec_mat, atm_mat, atm_el, nb_max_OHvec)
END PROGRAM h_bonds