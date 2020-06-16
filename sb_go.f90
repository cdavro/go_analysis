MODULE SB_GO
    IMPLICIT NONE
    ! -----------------------------------------------Set Double precision
    INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)

    CONTAINS

    SUBROUTINE sb_read_pos_xyz(sb_file_pos,sb_nb_atm,sb_nb_step,sb_atm_mat)
        CHARACTER(LEN=64), INTENT(IN)           :: sb_file_pos
        INTEGER, INTENT(IN)                     :: sb_nb_atm, sb_nb_step
        REAL(dp), INTENT(INOUT)                 :: sb_atm_mat(:,:,:)
        CHARACTER(LEN=3), ALLOCATABLE           :: sb_atm_name(:)
        INTEGER                                 :: s, i

        ALLOCATE(sb_atm_name(sb_nb_atm))

        OPEN(UNIT=20,FILE=sb_file_pos,STATUS='old',FORM='formatted',ACTION='READ')
        DO s = 1, sb_nb_step
            READ(20, *)
            READ(20, *)
            DO i=1,sb_nb_atm
                sb_atm_mat(2,i,s) = -1
                READ(20, *) sb_atm_name(i), sb_atm_mat(4,i,s), sb_atm_mat(5,i,s), sb_atm_mat(6,i,s)
                sb_atm_mat(1,i,s) = i
                IF (sb_atm_name(i) .EQ. "CC") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 1
                ELSE IF (sb_atm_name(i) .EQ. "C2O") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 3
                ELSE IF (sb_atm_name(i) .EQ. "C3O") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 4
                ELSE IF (sb_atm_name(i) .EQ. "C2A") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 5
                ELSE IF (sb_atm_name(i) .EQ. "C3A") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 6
                ELSE IF (sb_atm_name(i) .EQ. "C2E") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 7
                ELSE IF (sb_atm_name(i) .EQ. "C3E") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 8
                ELSE IF (sb_atm_name(i) .EQ. "OEP") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 10
                ELSE IF (sb_atm_name(i) .EQ. "OET") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 11
                ELSE IF ((sb_atm_name(i) .EQ. "OH2") .OR. &
                    (sb_atm_name(i) .EQ. "OH3")) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 12
                ELSE IF ((sb_atm_name(i) .EQ. "OA2") .OR. &
                    (sb_atm_name(i) .EQ. "OA3")) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 14
                ELSE IF (sb_atm_name(i) .EQ. "OW") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 19
                ELSE IF (sb_atm_name(i) .EQ. "OM") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 18
                ELSE IF (sb_atm_name(i) .EQ. "OP") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 17
                ELSE IF (sb_atm_name(i) .EQ. "HO") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 23
                ELSE IF (sb_atm_name(i) .EQ. "HW") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 29
                ELSE IF (sb_atm_name(i) .EQ. "HM") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 28
                ELSE IF (sb_atm_name(i) .EQ. "HP") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 27
                ELSE IF ((sb_atm_name(i) .EQ. "C").OR.&
                    (sb_atm_name(i) .EQ. "C2").OR.&
                    (sb_atm_name(i) .EQ. "C3").OR.&
                    (sb_atm_name(i) .EQ. "CX")) THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ((sb_atm_name(i) .EQ. "O") .OR.&
                    (sb_atm_name(i) .EQ. "OH") .OR. &
                    (sb_atm_name(i) .EQ. "OE") .OR. &
                    (sb_atm_name(i) .EQ. "OA") .OR. &
                    (sb_atm_name(i) .EQ. "OX")) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ((sb_atm_name(i) .EQ. "H") .OR.&
                    (sb_atm_name(i) .EQ. "H1") .OR.&
                    (sb_atm_name(i) .EQ. "H2") .OR.&
                    (sb_atm_name(i) .EQ. "HX")) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = -1
                END IF
            END DO
        END DO
        CLOSE(UNIT=20)
    DEALLOCATE(sb_atm_name)
    END SUBROUTINE sb_read_pos_xyz


    SUBROUTINE sb_count_is(sb_file_is,sb_nb_step,sb_nb_max_is)
        CHARACTER(LEN=64), INTENT(IN)           :: sb_file_is
        INTEGER, INTENT(IN)                     :: sb_nb_step
        INTEGER, INTENT(OUT)                    :: sb_nb_max_is
        !-- Procedure only
        INTEGER                                 :: sb_nb_line_is=0
        INTEGER                                 :: sb_iostatus=0

        OPEN(UNIT=21, FILE=sb_file_is, STATUS='old', FORM='formatted', ACTION='READ')
        DO
            READ(21, *, IOSTAT=sb_iostatus)
            IF (sb_iostatus .NE. 0) THEN
                EXIT
            ELSE
                sb_nb_line_is = sb_nb_line_is + 1
            END IF
        END DO
        REWIND(21)
        sb_nb_max_is = CEILING(1.0 * sb_nb_line_is / sb_nb_step) * 2
    END SUBROUTINE sb_count_is


    SUBROUTINE sb_read_is(sb_file_is,sb_nb_step,sb_box,sb_is_mat,sb_nb_is)
        CHARACTER(LEN=64), INTENT(IN)           :: sb_file_is
        INTEGER, INTENT(IN)                     :: sb_nb_step
        REAL(dp), INTENT(IN)                    :: sb_box(:)
        REAL(dp), INTENT(INOUT)                 :: sb_is_mat(:,:,:)
        INTEGER, INTENT(INOUT)                  :: sb_nb_is(:)
        !-- Procedure only
        INTEGER                                 :: s, i, j, k
        CHARACTER(LEN=2)                        :: sb_dummy_char

        OPEN(UNIT=21, FILE=sb_file_is, STATUS='old', FORM='formatted', ACTION='READ')
        DO s = 1, sb_nb_step
            READ(21, *) sb_nb_is(s)
            READ(21, *)
            j = 0
            DO i=1,sb_nb_is(s)
                READ(21, *) sb_dummy_char, sb_is_mat(1,i,s), sb_is_mat(2,i,s), sb_is_mat(3,i,s)
                j = j + 1
                sb_is_mat(5,i,s) = j
                IF (sb_dummy_char .EQ. "XD") THEN
                    sb_is_mat(4,i,s) = 1
                ELSE IF (sb_dummy_char .EQ. "XU") THEN
                    sb_is_mat(4,i,s) = 2
                ELSE
                    sb_is_mat(4,i,s) = 0
                END IF
                DO k = 1, 3
                    sb_is_mat(k,i,s) = sb_is_mat(k,i,s) - sb_box(k) * ANINT(sb_is_mat(k,i,s) / sb_box(k))
                END DO
            END DO
        END DO
        CLOSE(UNIT=21)

    END SUBROUTINE sb_read_is

END MODULE SB_GO