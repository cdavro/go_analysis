!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM order_layer
USE OMP_LIB
USE INPUT
USE SB_GO

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
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
INTEGER                         :: nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: SpO_disp_vec(3), SpO_disp_norm
CHARACTER(LEN=3)                :: type

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_h

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j
INTEGER                         :: count_L0, count_L1, count_L2, count_L3, count_all

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Order Layer'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF ( CAC .EQ. 0 ) THEN
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
ALLOCATE(atm_name(nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

IF ( IS_c .EQ. 'Y' ) THEN
! A ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
    start = OMP_get_wtime()

    CALL sb_count_is(file_is,nb_step,nb_max_is)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

    ! A ----------------------------------------------- Read IS
    start = OMP_get_wtime()

    ALLOCATE(is_mat(5,nb_max_is,nb_step))
    ALLOCATE(nb_is(nb_step))
    is_mat(:,:,:) = 0.0_dp
    nb_is(:) = 0

    CALL sb_read_is(file_is,nb_step,box(:),is_mat(:,:,:),nb_is(:))

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
END IF

IF ( IS_c .EQ. 'Y' ) THEN

    ! X ----------------------------------------------- Calculate closest distance between IS and Water-type Oxygen
    start = OMP_get_wtime()
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_atm, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpO_disp_vec, SpO_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF ( ( atm_mat(3,i,s) .EQ. 34 ) .OR. &
            ( atm_mat(3,i,s) .EQ. 33 ) .OR. &
            ( atm_mat(3,i,s) .EQ. 35 ) ) THEN ! Water-type Oxygen
                DO j = 1, nb_is(s)
                    IF ( is_mat(4,j,s) .EQ. 1 ) THEN
                        DO k = 1, 3
                            SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                            SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT( SpO_disp_vec(k)/box(k) )
                        END DO
                        SpO_disp_norm = NORM2( SpO_disp_vec )
                        IF ( ( SpO_disp_norm .LT. atm_mat(15,i,s) ) .OR. ( atm_mat(15,i,s) .EQ. 0.0 ) ) THEN
                            atm_mat(15,i,s) = SpO_disp_norm
                            IF ( atm_mat(6,i,s) .LT. is_mat(3,j,s) ) THEN
                                atm_mat(16,i,s) = -1
                            ELSE
                                atm_mat(16,i,s) = 1
                            END IF
                            atm_mat(17,i,s) = is_mat(5,j,s)
                        END IF
                    ELSE IF ( is_mat(4,j,s) .EQ. 2 ) THEN
                        DO k = 1, 3
                            SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                            SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT( SpO_disp_vec(k)/box(k) )
                        END DO
                        SpO_disp_norm = NORM2( SpO_disp_vec )
                        IF ( ( SpO_disp_norm .LT. atm_mat(18,i,s) ) .OR. (atm_mat(18,i,s) .EQ. 0.0 ) ) THEN
                            atm_mat(18,i,s) = SpO_disp_norm
                            IF ( atm_mat(6,i,s) .GT. is_mat(3,j,s) ) THEN
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
    PRINT'(A40,F14.2,A20)', "IS and OW-type:", finish-start, "seconds elapsed"

END IF

!   ----------------------------------------------- Print the Order Layer xyz
start = OMP_get_wtime()

OPEN(UNIT=40, FILE = suffix//"_order.xyz")
OPEN(UNIT=41, FILE = suffix//"_order.cell")
OPEN(UNIT=42, FILE = suffix//"_order_count.txt")
WRITE(42,'(A4,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10)') "Traj", "Step", "L0 Count"&
, "L1 Count", "L2 Count", "L3 Count", "LA Count", "TW Count"
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_atm-nb_h
    WRITE(40,'(F7.4,1X,F7.4,1X,F7.4)') box(1), box(2), box(3)
    WRITE(41,'(F7.4,1X,F7.4,1X,F7.4)') box(1), box(2), box(3)
    count_L0 = 0
    count_L1 = 0
    count_L2 = 0
    count_L3 = 0
    count_all = 0
    Z1:DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 1 ) CYCLE Z1
        type = atm_name(i,s)
        IF ( atm_name(i,s) .EQ. "OW" ) THEN
            type = "OW"
            count_all = count_all + 1
        END IF
        IF ( ( atm_name(i,s) .EQ. "OW" ) .AND. &
        (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L0_down ) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L0_up ) ) THEN
            type = "OL0"
            count_L0 = count_L0 + 1
        ELSE IF ( ( atm_name(i,s) .EQ. "OW" ) .AND. &
        (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L1_down ) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L1_up ) ) THEN
            type = "OL1"
            count_L1 = count_L1 + 1
        ELSE IF ( ( atm_name(i,s) .EQ. "OW" ) .AND. &
        (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L2_down ) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L2_up ) ) THEN
            type = "OL2"
            count_L2 = count_L2 + 1
        ELSE IF ( ( atm_name(i,s) .EQ. "OW" ) .AND. &
        (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. L3_down ) .AND.  (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. L3_up ) ) THEN
            type = "OL3"
            count_L3 = count_L3 + 1
        END IF
        WRITE(40,'(A4,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL( type ), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
    END DO Z1
    WRITE(42,'(A4,1X,I10,1X,I10,1X,I10,1X,I10,1X,I10,1X,I10,1X,I10)') suffix, s, count_L0, count_L1, count_L2, count_L3&
    , count_L0+count_L1+count_L2+count_L3, count_all
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
IF ( IS_c .EQ. 'Y' ) DEALLOCATE(is_mat,nb_is)

DEALLOCATE(atm_mat)

END PROGRAM order_layer