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
                IF (sb_atm_name(i) .EQ. "C") THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 1
                ELSE IF (sb_atm_name(i) .EQ. "OE") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 10
                ELSE IF (sb_atm_name(i) .EQ. "OH") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 11
                ELSE IF (sb_atm_name(i) .EQ. "OA") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 12
                ELSE IF (sb_atm_name(i) .EQ. "OW") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 13
                ELSE IF (sb_atm_name(i) .EQ. "OM") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 14
                ELSE IF (sb_atm_name(i) .EQ. "OP") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 15
                ELSE IF (sb_atm_name(i) .EQ. "O") THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = -1
                ELSE IF (sb_atm_name(i) .EQ. "HW") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 23
                ELSE IF (sb_atm_name(i) .EQ. "HO") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 21
                ELSE IF (sb_atm_name(i) .EQ. "H") THEN
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
        CHARACTER(LEN=1)                        :: sb_dummy_char
        
        OPEN(UNIT=21, FILE=sb_file_is, STATUS='old', FORM='formatted', ACTION='READ')
        DO s = 1, sb_nb_step
            READ(21, *) sb_nb_is(s)
            READ(21, *)
            j = 0
            DO i=1,sb_nb_is(s)
                READ(21, *) sb_dummy_char, sb_is_mat(1,i,s), sb_is_mat(2,i,s), sb_is_mat(3,i,s)
                j = j + 1
                sb_is_mat(5,i,s) = j
                IF (sb_is_mat(3,i,s) .LT. 10.0) THEN
                    sb_is_mat(4,i,s) = 1
                ELSE
                    sb_is_mat(4,i,s) = 2
                END IF
                DO k = 1, 3
                    sb_is_mat(k,i,s) = sb_is_mat(k,i,s) - sb_box(k) * ANINT(sb_is_mat(k,i,s) / sb_box(k))
                END DO
            END DO
        END DO
        CLOSE(UNIT=21)

    END SUBROUTINE sb_read_is


END MODULE SB_GO