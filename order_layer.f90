!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM order_layer
USE OMP_LIB
USE INPUT

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file
INTEGER                         :: CAC

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), is_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:)
INTEGER                         :: nb_line_is, nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: SpO_disp_vec(3), SpO_disp_norm
CHARACTER(LEN=3)                :: dummy_char, type

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_h

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j
INTEGER                         :: count_L0, count_L1, count_L2, count_L3

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Order Layer'
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
PRINT'(A100)', 'Run, Order Layer, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(25,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_name(nb_atm))

OPEN(UNIT=20,FILE=file_pos,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i=1,nb_atm
        atm_mat(2,i,s) = -1
        READ(20, *) atm_name(i), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
        atm_mat(1,i,s) = i
        IF (atm_name(i) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
            atm_mat(3,i,s) = 1
        ELSE IF (atm_name(i) .EQ. "Si") THEN
            atm_mat(2,i,s) = 28
            atm_mat(3,i,s) = 2
        ELSE IF (atm_name(i) .EQ. "SiF") THEN
            atm_mat(2,i,s) = 28
            atm_mat(3,i,s) = 1
        ELSE IF (atm_name(i) .EQ. "OE") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 10
        ELSE IF (atm_name(i) .EQ. "OH") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 11
        ELSE IF (atm_name(i) .EQ. "OA") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 12
        ELSE IF (atm_name(i) .EQ. "OW") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 13
        ELSE IF (atm_name(i) .EQ. "OM") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 14
        ELSE IF (atm_name(i) .EQ. "OP") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 15
        ELSE IF (atm_name(i) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = -1
        ELSE IF (atm_name(i) .EQ. "HW") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 23
        ELSE IF (atm_name(i) .EQ. "HO") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 21
        ELSE IF (atm_name(i) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = -1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

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

IF (IS_c .EQ. 'Y' ) THEN

    ! X ----------------------------------------------- Calculate closest distance between IS and any O atom
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_atm, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpO_disp_vec, SpO_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF (atm_mat(3,i,s) .EQ. 13) THEN
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

!   ----------------------------------------------- Print the Order Layer xyz
start = OMP_get_wtime()

OPEN(UNIT=40, FILE = suffix//"_order.xyz")
OPEN(UNIT=41, FILE = suffix//"_cell.cell")
OPEN(UNIT=42, FILE = suffix//"_count.txt")
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_atm-nb_h
    WRITE(40,'(F5.2,A1,F5.2,A1,F5.2)') box(1)," ", box(2)," ", box(3)
    WRITE(41,'(F5.2,A1,F5.2,A1,F5.2)') box(1)," ", box(2)," ", box(3)
    count_L0 = 0
    count_L1 = 0
    count_L2 = 0
    count_L3 = 0
    Z1:DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 1) CYCLE Z1
        IF ( (atm_mat(2,i,s) .EQ. 12) .AND. (atm_mat(3,i,s) .EQ. 1) ) THEN
            type = "C"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 28) .AND. (atm_mat(3,i,s) .EQ. 2) ) THEN
            type ="Si"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 28) .AND. (atm_mat(3,i,s) .EQ. 1) ) THEN
            type ="SiF"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 10) ) THEN
            type = "OE"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 11) ) THEN
            type = "OH"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 12) ) THEN
            type = "OA"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 13) ) THEN
            type = "OW"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 14) ) THEN
            type = "OM"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 15) ) THEN
            type = "OP"
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. -1) ) THEN
            type = "OP"
        END IF
        IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 13) &
            .AND. (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L0_down) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L0_up) ) THEN
            type = "O0"
            count_L0 = count_L0 + 1
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 13) &
            .AND. (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L1_down) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L1_up) ) THEN
            type = "O1"
            count_L1 = count_L1 + 1
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 13) &
            .AND. (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L2_down) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L2_up) ) THEN
            type = "O2"
            count_L2 = count_L2 + 1
        ELSE IF ( (atm_mat(2,i,s) .EQ. 16) .AND. (atm_mat(3,i,s) .EQ. 13) &
            .AND. (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L3_down) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L3_up) ) THEN
            type = "O3"
            count_L3 = count_L3 + 1
        END IF
        WRITE(40,'(A10,E20.7,E20.7,E20.7)') ADJUSTL(type), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
    END DO Z1
    WRITE(42,'(I10,I10,I10,I10,I10)') s, count_L0, count_L1, count_L2, count_L3
END DO
CLOSE(UNIT=40)
CLOSE(UNIT=41)
CLOSE(UNIT=42)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions order output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat, nb_is)

DEALLOCATE(atm_mat, atm_name)
END PROGRAM order_layer