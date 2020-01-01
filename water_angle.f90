!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM water_angle
USE OMP_LIB
USE INPUT

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision ----------------------------------------------- !
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), WAT_mat(:,:,:), is_mat(:,:,:), as_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
INTEGER                         :: nb_line_is, nb_max_is, nb_max_as
INTEGER, ALLOCATABLE            :: nb_is(:), nb_as(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: OpH_disp_vec(3), OpH_disp_norm
REAL(dp)                        :: XpOwat_disp_vec(3), XpOwat_disp_norm
REAL(dp)                        :: SpOwat_disp_vec(3), SpOwat_disp_norm
REAL(dp)                        :: SpS1_disp_vec(3), SpS1_disp_norm, SpS2_disp_vec(3), SpS2_disp_norm
REAL(dp)                        :: tPSvec_down(3), tPSvec_up(3)
REAL(dp)                        :: tPSuvec_down(3), tPSuvec_up(3)
REAL(dp)                        :: WD_vec(3), HpH_disp_vec(3),WD_uvec(3), HpH_disp_uvec(3)
REAL(dp)                        :: tO_pos_iPS_down(3), tO_pos_iPS_up(3), tS_pos_PS_down(3), tS_pos_PS_up(3)
REAL(dp)                        :: tOS_disp_oPS_down(3), tOS_disp_oPS_up(3)
CHARACTER(LEN=2)                :: dummy_char

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o
INTEGER,ALLOCATABLE             :: nb_max_WAT(:)

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o
INTEGER                         :: CAC

!   ----------------------------------------------- Parameters

!   -----------------------------------------------
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'Launching Water Angle'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*, "No input files"
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
PRINT'(A100)', 'Run, Water Angle, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(6,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_el(nb_atm))

OPEN(UNIT=20, FILE=file_pos, STATUS='old', FORM='formatted', ACTION='READ')
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
            nb_line_is = nb_line_is+1
        END IF
    END DO
    REWIND(21)
    nb_max_is = CEILING(1.0 * nb_line_is / nb_step) * 2

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

    ! A ----------------------------------------------- Read IS
    start = OMP_get_wtime()

    ALLOCATE(is_mat(19,nb_max_is,nb_step))
    ALLOCATE(nb_is(nb_step))
    is_mat(:,:,:) = 0.0_dp
    nb_is(:) = 0

    OPEN(UNIT=21, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
    DO s = 1, nb_step
        READ(21, *) nb_is(s)
        READ(21, *)
        j = 0
        DO i = 1, nb_is(s)
            READ(21, *) dummy_char, is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
            j = j + 1
            is_mat(5,i,s) = j
            IF (is_mat(3,i,s) .LT. 10.0) THEN
                is_mat(4,i,s) = 1
            ELSE
                is_mat(4,i,s) = 2
            END IF
            DO k = 1, 3
                is_mat(k,i,s) = is_mat(k,i,s) - box(k) * ANINT(is_mat(k,i,s) / box(k))
            END DO
        END DO
    END DO
    CLOSE(UNIT=21)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
END IF

IF (AS_c .EQ. 'Y') THEN
    ! A ----------------------------------------------- Read AS
        start = OMP_get_wtime()
        ALLOCATE(nb_as(nb_step))
        nb_as(:) = 0
        DO s = 1, nb_step
            nb_as(s) = COUNT(atm_mat(3,:,s) .EQ. 1, DIM=1)
        END DO
        nb_max_as = MAXVAL(nb_as)

        ALLOCATE(as_mat(19,nb_max_as,nb_step))
        as_mat(:,:,:) = 0.0_dp
        j = 0

        DO s = 1, nb_step
            j = 0
            DO i = 1, nb_atm
                IF (atm_mat(3,i,s) .EQ. wa_AS_center_atmnb) THEN
                    j = j + 1
                    as_mat(5,j,s) = atm_mat(1,i,s)
                    DO k = 1, 3
                        as_mat(k,j,s) = atm_mat(k+3,i,s)
                    END DO
                END IF
            END DO
        END DO

        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "AS:", finish-start, "seconds elapsed"
END IF

! B ----------------------------------------------- Water
start = OMP_get_wtime()
ALLOCATE(WAT_mat(46,nb_o*3,nb_step))
ALLOCATE(nb_max_WAT(nb_step))

nb_max_WAT(:) = 0.0
WAT_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_max_WAT, nb_o, nb_step, nb_atm)&
!$OMP PRIVATE(s, i, j, k, o)&
!$OMP PRIVATE(OpH_disp_vec, OpH_disp_norm)
DO s = 1, nb_step
    o = 0
    DO i = 1, nb_atm
        IF (atm_mat(3,i,s) .EQ. 13) THEN
            o = o + 1
            WAT_mat(1,o,s) = atm_mat(1,i,s)
            DO k = 1, 3
                WAT_mat(k+1,o,s) = atm_mat(k+3,i,s)
            END DO
            DO j = 1, nb_atm
                IF (atm_mat(3,j,s) .EQ. 23) THEN
                    DO k = 1, 3
                        OpH_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        OpH_disp_vec(k) = OpH_disp_vec(k) - box(k) * ANINT(OpH_disp_vec(k) / box(k))
                    END DO

                    OpH_disp_norm = NORM2(OpH_disp_vec(:))
                    IF ( ( (OpH_disp_norm .LT. WAT_mat(6,o,s)) .OR.&
                    (WAT_mat(6,o,s) .EQ. 0.0) ) .AND.&
                    (atm_mat(1,j,s) .NE. WAT_mat(5,o,s)) ) THEN
                        WAT_mat(10,o,s) = WAT_mat(5,o,s)
                        WAT_mat(11,o,s) = WAT_mat(6,o,s)
                        WAT_mat(5,o,s) = atm_mat(1,j,s)
                        WAT_mat(6,o,s) = OpH_disp_norm
                        DO k = 1, 3
                            WAT_mat(k+11,o,s) = WAT_mat(k+6,o,s)
                            WAT_mat(k+6,o,s) =  atm_mat(k+3,j,s)
                        END DO

                    ELSE IF ( ( (OpH_disp_norm .LT. WAT_mat(11,o,s)) .OR.&
                        (WAT_mat(6,o,s) .NE. 0.0) .AND.&
                        (WAT_mat(11,o,s) .EQ. 0.0 )) .AND.&
                        (atm_mat(1,j,s) .NE. WAT_mat(10,o,s) ) ) THEN
                        WAT_mat(10,o,s) = atm_mat(1,j,s)
                        WAT_mat(11,o,s) = OpH_disp_norm
                        DO k = 1, 3
                            WAT_mat(k+11,o,s) = atm_mat(k+3,j,s)
                        END DO

                    END IF
                END IF
            END DO
        END IF
    END DO
    DO i=1,nb_o*3
        DO k = 1, 3
            WAT_mat(k+20,i,s) = WAT_mat(k+11,i,s) - WAT_mat(k+6,i,s)
            WAT_mat(k+20,i,s) = WAT_mat(k+20,i,s) - box(k) * ANINT(WAT_mat(k+20,i,s) / box(k)) ! HH Vector
            WAT_mat(k+14,i,s) = WAT_mat(k+6,i,s) + WAT_mat(k+20,i,s)/2.0 ! Center of HH
            WAT_mat(k+17,i,s) = WAT_mat(k+14,i,s) - WAT_mat(k+1,i,s)
            WAT_mat(k+17,i,s) = WAT_mat(k+17,i,s) - box(k) * ANINT(WAT_mat(k+17,i,s) / box(k)) ! (O to center HH) vector (Inverse of of the Dipole Moment)
        END DO
    END DO
    nb_max_WAT(s) = COUNT(WAT_mat(1,:,s) .NE. 0, DIM=1)
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "WAT groups:", finish-start, "seconds elapsed"


! C ----------------------------------------------- Proximity between functionnal groups and any WAT groups
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, nb_step, nb_atm)&
!$OMP SHARED(wa_XpOwat_rcut, wa_CpOwat_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XpOwat_disp_vec, XpOwat_disp_norm)
DO s = 1, nb_step
    C1:DO i = 1, nb_o*3
        IF (WAT_mat(1,i,s) .EQ. 0) THEN
            CYCLE C1
        END IF
        C2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE C2
            END IF
            IF (atm_mat(1,j,s) .EQ. WAT_mat(1,i,s)) THEN
                CYCLE C2
            END IF
            DO k = 1, 3
                XpOwat_disp_vec(k) = atm_mat(k+3,j,s) - WAT_mat(k+1,i,s)
                XpOwat_disp_vec(k) = XpOwat_disp_vec(k) - box(k) * ANINT(XpOwat_disp_vec(k) / box(k))
            END DO
            XpOwat_disp_norm = NORM2(XpOwat_disp_vec)

            IF ( ( (XpOwat_disp_norm .LE. WAT_mat(45,i,s)) .OR.&
            (WAT_mat(45,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                WAT_mat(44,i,s) = atm_mat(1,j,s)
                WAT_mat(45,i,s) = XpOwat_disp_norm
                WAT_mat(46,i,s) = atm_mat(6,j,s)
            END IF
            IF (XpOwat_disp_norm .LE. wa_XpOwat_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C or SiF
                    WAT_mat(27,i,s) = 1
                    WAT_mat(28,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    WAT_mat(24,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    WAT_mat(25,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    WAT_mat(26,i,s) = 1
                END IF
            ELSE IF (XpOwat_disp_norm .LE. wa_CpOwat_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C or SiF
                    WAT_mat(28,i,s) = 1
                END IF
            END IF
        END DO C2
    END DO C1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity WAT groups and FG groups:", finish-start, "seconds elapsed"

!   ----------------------------------------------- IS
IF (IS_c .EQ. 'Y' ) THEN
    ! D ----------------------------------------------- Calculate WAT angle IS
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, is_mat, nb_is, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpOwat_disp_vec, SpOwat_disp_norm, SpS1_disp_vec, SpS1_disp_norm, SpS2_disp_vec, SpS2_disp_norm)&
    !$OMP PRIVATE(tPSuvec_down, tPSuvec_up, tPSvec_down, tPSvec_up, WD_vec, HpH_disp_vec, WD_uvec, HpH_disp_uvec)&
    !$OMP PRIVATE(tO_pos_iPS_down, tO_pos_iPS_up, tS_pos_PS_down, tS_pos_PS_up)&
    !$OMP PRIVATE(tOS_disp_oPS_down, tOS_disp_oPS_up)
    DO s = 1, nb_step
        D2:DO i = 1, nb_o*3
            IF (WAT_mat(1,i,s) .EQ. 0) THEN
                CYCLE D2
            END IF
            DO j = 1, nb_is(s)
                IF (is_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        SpOwat_disp_vec(k) = is_mat(k,j,s) - WAT_mat(k+1,i,s)
                        SpOwat_disp_vec(k) = SpOwat_disp_vec(k) - box(k) * ANINT(SpOwat_disp_vec(k) / box(k))
                    END DO
                    SpOwat_disp_norm = NORM2(SpOwat_disp_vec)
                    IF ( (SpOwat_disp_norm .LT. WAT_mat(29,i,s)) .OR. (WAT_mat(29,i,s) .EQ. 0.0) ) THEN
                        WAT_mat(29,i,s) = SpOwat_disp_norm
                        IF (WAT_mat(4,i,s) .LT. is_mat(3,j,s)) THEN
                            WAT_mat(30,i,s) = -1
                        ELSE
                            WAT_mat(30,i,s) = 1
                        END IF
                        WAT_mat(31,i,s) = is_mat(5,j,s)
                    END IF
                ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        SpOwat_disp_vec(k) = is_mat(k,j,s) - WAT_mat(k+1,i,s)
                        SpOwat_disp_vec(k) = SpOwat_disp_vec(k) - box(k) * ANINT(SpOwat_disp_vec(k) / box(k))
                    END DO
                    SpOwat_disp_norm = NORM2(SpOwat_disp_vec)
                    IF ( (SpOwat_disp_norm .LT. WAT_mat(32,i,s)) .OR. (WAT_mat(32,i,s) .EQ. 0.0) ) THEN
                        WAT_mat(32,i,s) = SpOwat_disp_norm
                        IF (WAT_mat(4,i,s) .GT. is_mat(3,j,s)) THEN
                            WAT_mat(33,i,s) = -1
                        ELSE
                            WAT_mat(33,i,s) = 1
                        END IF
                        WAT_mat(34,i,s) = is_mat(5,j,s)
                    END IF
                END IF
            END DO

            IF ((is_mat(19,INT(WAT_mat(31,i,s)),s) .NE. 1) .OR.&
                (is_mat(19,INT(WAT_mat(34,i,s)),s) .NE. 1)) THEN

                D3:DO j = 1, nb_is(s) ! First one
                    IF ( (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(31,i,s)),s)) .OR.&
                    (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(34,i,s)),s)) ) THEN
                        CYCLE D3
                    END IF
                    IF (is_mat(4,j,s) .EQ. 1) THEN

                        DO k = 1, 3
                            SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(31,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                        END DO
                        SpS1_disp_norm = NORM2(SpS1_disp_vec)

                        IF ( ( (SpS1_disp_norm .LT. is_mat(7,INT(WAT_mat(31,i,s)),s)) .OR.&
                        (is_mat(7,INT(WAT_mat(31,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (is_mat(5,j,s) .NE. is_mat(6,INT(WAT_mat(31,i,s)),s)) ) THEN
                            is_mat(6,INT(WAT_mat(31,i,s)),s) = is_mat(5,j,s)
                            is_mat(7,INT(WAT_mat(31,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                is_mat(k+7,INT(WAT_mat(31,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                            END DO
                        END IF

                    ELSE IF (is_mat(4,j,s) .EQ. 2) THEN

                        DO k = 1, 3
                            SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(34,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                        END DO
                        SpS1_disp_norm = NORM2(SpS1_disp_vec)

                        IF ( ( (SpS1_disp_norm .LT. is_mat(7,INT(WAT_mat(34,i,s)),s)) .OR.&
                        (is_mat(7,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (is_mat(5,j,s) .NE. is_mat(6,INT(WAT_mat(34,i,s)),s)) ) THEN
                            is_mat(6,INT(WAT_mat(34,i,s)),s) = is_mat(5,j,s)
                            is_mat(7,INT(WAT_mat(34,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                is_mat(k+7,INT(WAT_mat(34,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !18,19,20
                            END DO
                        END IF

                    END IF
                END DO D3

            D4:DO j = 1, nb_is(s)

                    IF ( (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(31,i,s)),s)) .OR.&
                    (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(34,i,s)),s)) .OR. &
                    (is_mat(5,j,s) .EQ. is_mat(6,INT(WAT_mat(31,i,s)),s)) .OR. &
                    (is_mat(5,j,s) .EQ. is_mat(6,INT(WAT_mat(34,i,s)),s)) ) THEN
                        CYCLE D4
                    END IF

                    IF (is_mat(4,j,s) .EQ. 1) THEN
                        DO k = 1, 3
                            SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(31,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                            SpS2_disp_vec(k) = is_mat(k+7,INT(WAT_mat(31,i,s)),s)
                            SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                        END DO

                        SpS1_disp_norm = NORM2(SpS1_disp_vec)
                        SpS2_disp_norm = NORM2(SpS2_disp_vec)

                        IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                        (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) ) THEN
                            CYCLE D4
                        END IF

                        IF ( ( (SpS1_disp_norm .LT. is_mat(12,INT(WAT_mat(31,i,s)),s)) .OR.&
                        (is_mat(12,INT(WAT_mat(31,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (is_mat(5,j,s) .NE. is_mat(11,INT(WAT_mat(31,i,s)),s)) ) THEN

                            is_mat(11,INT(WAT_mat(31,i,s)),s) = is_mat(5,j,s)
                            is_mat(12,INT(WAT_mat(31,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                is_mat(k+12,INT(WAT_mat(31,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                            END DO

                        END IF

                    ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                        DO k = 1, 3
                            SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(34,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                            SpS2_disp_vec(k) = is_mat(k+7,INT(WAT_mat(34,i,s)),s)
                            SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                        END DO

                        SpS1_disp_norm = NORM2(SpS1_disp_vec)
                        SpS2_disp_norm = NORM2(SpS2_disp_vec)

                        IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                        (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) ) THEN
                            CYCLE D4
                        END IF

                        IF ( ( (SpS1_disp_norm .LT. is_mat(12,INT(WAT_mat(34,i,s)),s)) .OR.&
                        (is_mat(12,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (is_mat(5,j,s) .NE. is_mat(11,INT(WAT_mat(34,i,s)),s)) ) THEN

                            is_mat(11,INT(WAT_mat(34,i,s)),s) = is_mat(5,j,s)
                            is_mat(12,INT(WAT_mat(34,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                is_mat(k+12,INT(WAT_mat(34,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !18,19,20
                            END DO

                        END IF
                END IF

                END DO D4

                is_mat(16,INT(WAT_mat(31,i,s)),s) =&
                    is_mat(9,INT(WAT_mat(31,i,s)),s) * is_mat(15,INT(WAT_mat(31,i,s)),s) -&
                    is_mat(10,INT(WAT_mat(31,i,s)),s) * is_mat(14,INT(WAT_mat(31,i,s)),s)
                is_mat(17,INT(WAT_mat(31,i,s)),s) =&
                    is_mat(10,INT(WAT_mat(31,i,s)),s) * is_mat(13,INT(WAT_mat(31,i,s)),s) -&
                    is_mat(8,INT(WAT_mat(31,i,s)),s) * is_mat(15,INT(WAT_mat(31,i,s)),s)
                is_mat(18,INT(WAT_mat(31,i,s)),s) =&
                    is_mat(8,INT(WAT_mat(31,i,s)),s) * is_mat(14,INT(WAT_mat(31,i,s)),s) -&
                    is_mat(9,INT(WAT_mat(31,i,s)),s) * is_mat(13,INT(WAT_mat(31,i,s)),s)

                is_mat(16,INT(WAT_mat(34,i,s)),s) =&
                    is_mat(9,INT(WAT_mat(34,i,s)),s) * is_mat(15,INT(WAT_mat(34,i,s)),s) -&
                    is_mat(10,INT(WAT_mat(34,i,s)),s) * is_mat(14,INT(WAT_mat(34,i,s)),s)
                is_mat(17,INT(WAT_mat(34,i,s)),s) =&
                    is_mat(10,INT(WAT_mat(34,i,s)),s) * is_mat(13,INT(WAT_mat(34,i,s)),s) -&
                    is_mat(8,INT(WAT_mat(34,i,s)),s) * is_mat(15,INT(WAT_mat(34,i,s)),s)
                is_mat(18,INT(WAT_mat(34,i,s)),s) =&
                    is_mat(8,INT(WAT_mat(34,i,s)),s) * is_mat(14,INT(WAT_mat(34,i,s)),s) -&
                    is_mat(9,INT(WAT_mat(34,i,s)),s) * is_mat(13,INT(WAT_mat(34,i,s)),s)

                is_mat(19,INT(WAT_mat(31,i,s)),s) = 1
                is_mat(19,INT(WAT_mat(34,i,s)),s) = 1

            END IF
            DO k = 1, 3
                tPSvec_down(k) = is_mat(k+15,INT(WAT_mat(31,i,s)),s)
                tPSvec_up(k) = is_mat(k+15,INT(WAT_mat(34,i,s)),s)
                WD_vec(k) = WAT_mat(k+17,i,s)
                HpH_disp_vec(k) = WAT_mat(k+20,i,s)
            END DO

            tPSuvec_down(:) = tPSvec_down(:) / NORM2(tPSvec_down(:))
            tPSuvec_up(:) = tPSvec_up(:) / NORM2(tPSvec_up(:))
            WD_uvec(:) = WD_vec(:) / NORM2(WD_vec(:))
            HpH_disp_uvec(:) = HpH_disp_vec(:) / NORM2(HpH_disp_vec(:))
            DO k = 1, 3
                ! To check orientation of the normal isace vector
                tO_pos_iPS_down(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_down(k)
                tO_pos_iPS_down(k) = tO_pos_iPS_down(k) - box(k) * ANINT(tO_pos_iPS_down(k) / box(k))
                tS_pos_PS_down(k) = is_mat(k,INT(WAT_mat(31,i,s)),s) + 0.01*tPSuvec_down(k)
                tS_pos_PS_down(k) = tS_pos_PS_down(k) - box(k) * ANINT(tS_pos_PS_down(k) / box(k))
                tOS_disp_oPS_down(k) = tS_pos_PS_down(k) - tO_pos_iPS_down(k)
                tOS_disp_oPS_down(k) = tOS_disp_oPS_down(k) - box(k) * ANINT(tOS_disp_oPS_down(k) / box(k))
                tO_pos_iPS_up(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_up(k)
                tO_pos_iPS_up(k) = tO_pos_iPS_up(k) - box(k) * ANINT(tO_pos_iPS_up(k) / box(k))
                tS_pos_PS_up(k) = is_mat(k,INT(WAT_mat(34,i,s)),s) + 0.01*tPSuvec_up(k)
                tS_pos_PS_up(k) = tS_pos_PS_up(k) - box(k) * ANINT(tS_pos_PS_up(k) / box(k))
                tOS_disp_oPS_up(k) = tS_pos_PS_up(k) - tO_pos_iPS_up(k)
                tOS_disp_oPS_up(k) = tOS_disp_oPS_up(k) - box(k) * ANINT(tOS_disp_oPS_up(k) / box(k))
            END DO

            IF (NORM2(tOS_disp_oPS_down(:)) .GT. WAT_mat(29,i,s)) THEN
                tPSuvec_down(:) = -1.0 * tPSuvec_down(:)
            END IF
            IF (NORM2(tOS_disp_oPS_up(:)) .GT. WAT_mat(32,i,s)) THEN
                tPSuvec_up(:) = -1.0 * tPSuvec_up(:)
            END IF

            IF (ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
            END IF

            WAT_mat(35,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), WD_uvec(:)))
            WAT_mat(36,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:)))

            IF (ACOS(DOT_PRODUCT(tPSuvec_up(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
            END IF

            WAT_mat(37,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), WD_uvec(:)))
            WAT_mat(38,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), HpH_disp_uvec(:)))

        END DO D2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity WAT groups and IS:", finish-start, "seconds elapsed"

    ! E ----------------------------------------------- Write WAT angle IS
    start = OMP_get_wtime()

    OPEN(UNIT=32, FILE = suffix//"_IS_water_angle.txt")
    WRITE(32, '(A10,A10,A10,A10,A10,A10,A20,A20,A20,A20,A20,A20,A10)')&
        "O_id", "cOE", "cOH", "cOA", "cC", "cC9"&
        , "dist_IS_down", "dist_IS_up"&
        , "Angle OH/NIS_down", "Angle OH/NIS_up"&
        , "Angle HH/NIS_down", "Angle HH/NIS_up"&
        , "step"
    DO s = 1, nb_step
        DO i = 1, nb_max_WAT(s)
            WRITE(32, '(I10,I10,I10,I10,I10,I10,E20.7,E20.7,E20.7,E20.7,E20.7,E20.7,I10)')&
            INT(WAT_mat(1,i,s)), INT(WAT_mat(24,i,s)), INT(WAT_mat(25,i,s))&
            , INT(WAT_mat(26,i,s)), INT(WAT_mat(27,i,s)), INT(WAT_mat(28,i,s))&
            , (WAT_mat(29,i,s)*WAT_mat(30,i,s)), (WAT_mat(32,i,s)*WAT_mat(33,i,s))&
            , WAT_mat(35,i,s), WAT_mat(37,i,s)&
            , WAT_mat(36,i,s), WAT_mat(38,i,s)&
            , s
        END DO
    END DO
    CLOSE(UNIT=32)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "WAT IS angles output:", finish-start, "seconds elapsed"

END IF

!   ----------------------------------------------- AS
IF (AS_c .EQ. 'Y' ) THEN
    ! E ----------------------------------------------- Calculate WAT angle AS
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, as_mat, nb_as, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpOwat_disp_vec, SpOwat_disp_norm, SpS1_disp_vec, SpS1_disp_norm, SpS2_disp_vec, SpS2_disp_norm)&
    !$OMP PRIVATE(tPSuvec_down, tPSvec_down, WD_vec, HpH_disp_vec ,WD_uvec, HpH_disp_uvec)&
    !$OMP PRIVATE(tO_pos_iPS_down, tS_pos_PS_down)&
    !$OMP PRIVATE(tOS_disp_oPS_down)
    DO s = 1, nb_step
        E2:DO i = 1, nb_o*3
            IF (WAT_mat(1,i,s) .EQ. 0) THEN
                CYCLE E2
            END IF
            DO j = 1, nb_as(s)
                DO k = 1, 3
                    SpOwat_disp_vec(k) = as_mat(k,j,s) - WAT_mat(k+1,i,s)
                    SpOwat_disp_vec(k) = SpOwat_disp_vec(k) - box(k) * ANINT(SpOwat_disp_vec(k) / box(k))
                END DO
                SpOwat_disp_norm = NORM2(SpOwat_disp_vec)
                IF ( (SpOwat_disp_norm .LT. WAT_mat(39,i,s)) .OR. (WAT_mat(39,i,s) .EQ. 0.0) ) THEN
                    WAT_mat(39,i,s) = SpOwat_disp_norm
                    IF (WAT_mat(4,i,s) .LT. as_mat(3,j,s)) THEN
                        WAT_mat(40,i,s) = -1
                    ELSE
                        WAT_mat(40,i,s) = 1
                    END IF
                    WAT_mat(41,i,s) = as_mat(5,j,s)
                END IF
            END DO

            IF ((as_mat(19,INT(WAT_mat(41,i,s)),s) .NE. 1)) THEN

                E3:DO j = 1, nb_as(s) ! First one
                    IF ( (as_mat(5,j,s) .EQ. as_mat(5,INT(WAT_mat(41,i,s)),s))) THEN
                        CYCLE E3
                    END IF

                    DO k = 1, 3
                        SpS1_disp_vec(k) = as_mat(k,j,s) - as_mat(k,INT(WAT_mat(41,i,s)),s)
                        SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                    END DO
                    SpS1_disp_norm = NORM2(SpS1_disp_vec)

                    IF ( ( (SpS1_disp_norm .LT. as_mat(7,INT(WAT_mat(41,i,s)),s)) .OR.&
                    (as_mat(7,INT(WAT_mat(41,i,s)),s) .EQ. 0.0 ) ) .AND.&
                    (as_mat(5,j,s) .NE. as_mat(6,INT(WAT_mat(41,i,s)),s)) ) THEN
                        as_mat(6,INT(WAT_mat(41,i,s)),s) = as_mat(5,j,s)
                        as_mat(7,INT(WAT_mat(41,i,s)),s) = SpS1_disp_norm
                        DO k = 1, 3
                            as_mat(k+7,INT(WAT_mat(41,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                        END DO
                    END IF

                END DO E3

                E4:DO j = 1, nb_as(s)

                    IF ( (as_mat(5,j,s) .EQ. as_mat(5,INT(WAT_mat(41,i,s)),s)) .OR.&
                    (as_mat(5,j,s) .EQ. as_mat(6,INT(WAT_mat(41,i,s)),s)) ) THEN
                        CYCLE E4
                    END IF

                        DO k = 1, 3
                            SpS1_disp_vec(k) = as_mat(k,j,s) - as_mat(k,INT(WAT_mat(41,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                            SpS2_disp_vec(k) = as_mat(k+7,INT(WAT_mat(41,i,s)),s)
                            SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                        END DO

                        SpS1_disp_norm = NORM2(SpS1_disp_vec)
                        SpS2_disp_norm = NORM2(SpS2_disp_vec)

                        IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                        (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) ) THEN
                            CYCLE E4
                        END IF

                        IF ( ( (SpS1_disp_norm .LT. as_mat(12,INT(WAT_mat(41,i,s)),s)) .OR.&
                        (as_mat(12,INT(WAT_mat(41,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (as_mat(5,j,s) .NE. as_mat(11,INT(WAT_mat(41,i,s)),s)) ) THEN

                            as_mat(11,INT(WAT_mat(41,i,s)),s) = as_mat(5,j,s)
                            as_mat(12,INT(WAT_mat(41,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                as_mat(k+12,INT(WAT_mat(41,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                            END DO

                        END IF

                END DO E4

                as_mat(16,INT(WAT_mat(41,i,s)),s) =&
                    as_mat(9,INT(WAT_mat(41,i,s)),s) * as_mat(15,INT(WAT_mat(41,i,s)),s) -&
                    as_mat(10,INT(WAT_mat(41,i,s)),s) * as_mat(14,INT(WAT_mat(41,i,s)),s)
                as_mat(17,INT(WAT_mat(41,i,s)),s) =&
                    as_mat(10,INT(WAT_mat(41,i,s)),s) * as_mat(13,INT(WAT_mat(41,i,s)),s) -&
                    as_mat(8,INT(WAT_mat(41,i,s)),s) * as_mat(15,INT(WAT_mat(41,i,s)),s)
                as_mat(18,INT(WAT_mat(41,i,s)),s) =&
                    as_mat(8,INT(WAT_mat(41,i,s)),s) * as_mat(14,INT(WAT_mat(41,i,s)),s) -&
                    as_mat(9,INT(WAT_mat(41,i,s)),s) * as_mat(13,INT(WAT_mat(41,i,s)),s)

                as_mat(19,INT(WAT_mat(41,i,s)),s) = 1

            END IF

            DO k = 1, 3
                tPSvec_down(k) = as_mat(k+15,INT(WAT_mat(41,i,s)),s)
                WD_vec(k) = WAT_mat(k+17,i,s)
                HpH_disp_vec(k) = WAT_mat(k+20,i,s)
            END DO
            tPSuvec_down(:) = tPSvec_down(:) / NORM2(tPSvec_down(:))
            WD_uvec(:) = WD_vec(:) / NORM2(WD_vec(:))
            HpH_disp_uvec(:) = HpH_disp_vec(:) / NORM2(HpH_disp_vec(:))

            DO k = 1, 3
                ! To check orientation of the normal isace vector
                tO_pos_iPS_down(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_down(k)
                tO_pos_iPS_down(k) = tO_pos_iPS_down(k) - box(k) * ANINT(tO_pos_iPS_down(k) / box(k))
                tS_pos_PS_down(k) = as_mat(k,INT(WAT_mat(41,i,s)),s) + 0.01*tPSuvec_down(k)
                tS_pos_PS_down(k) = tS_pos_PS_down(k) - box(k) * ANINT(tS_pos_PS_down(k) / box(k))
                tOS_disp_oPS_down(k) = tS_pos_PS_down(k) - tO_pos_iPS_down(k)
                tOS_disp_oPS_down(k) = tOS_disp_oPS_down(k) - box(k) * ANINT(tOS_disp_oPS_down(k) / box(k))
            END DO

            IF (NORM2(tOS_disp_oPS_down(:)) .GT. WAT_mat(39,i,s)) THEN
                tPSuvec_down(:) = -1.0 * tPSuvec_down(:)
            END IF

            IF (ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
            END IF

            WAT_mat(42,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), WD_uvec(:)))
            WAT_mat(43,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:)))

        END DO E2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity WAT groups and AS:", finish-start, "seconds elapsed"

    ! E ----------------------------------------------- Write WAT angle AS
    start = OMP_get_wtime()

    OPEN(UNIT=33, FILE = suffix//"_AS_water_angle.txt")
    WRITE(33, '(A10,A10,A10,A10,A10,A10,A20,A20,A20,A10)')&
        "O_id", "cOE", "cOH", "cOA", "cC", "cC9"&
        , "dist_AS_down", "Angle OH/NAS", "Angle HH/NAS", "step"
    DO s = 1, nb_step
        DO i = 1, nb_max_WAT(s)
            WRITE(33, '(I10,I10,I10,I10,I10,I10,E20.7,E20.7,E20.7, I10)')&
            INT(WAT_mat(1,i,s)), INT(WAT_mat(24,i,s)), INT(WAT_mat(25,i,s))&
            , INT(WAT_mat(26,i,s)), INT(WAT_mat(27,i,s)), INT(WAT_mat(28,i,s))&
            , (WAT_mat(39,i,s)*WAT_mat(40,i,s)), WAT_mat(42,i,s), WAT_mat(43,i,s)&
            , s
        END DO
    END DO
    CLOSE(UNIT=33)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "WAT AS angles output:", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat,nb_is)
IF (AS_c .EQ. 'Y') DEALLOCATE(as_mat,nb_as)

DEALLOCATE(WAT_mat,atm_mat,atm_el,nb_max_WAT)
END PROGRAM water_angle