MODULE SB_GO
    IMPLICIT NONE
    ! -----------------------------------------------Set Double precision
    INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)

    CONTAINS

    SUBROUTINE sb_read_pos_xyz(sb_file_pos,sb_nb_atm,sb_nb_step,sb_atm_mat,sb_atm_name)
        CHARACTER(LEN=64), INTENT(IN)                :: sb_file_pos
        INTEGER, INTENT(IN)                          :: sb_nb_atm, sb_nb_step
        REAL(dp), INTENT(INOUT)                      :: sb_atm_mat(:,:,:)
        CHARACTER(LEN=3), ALLOCATABLE, INTENT(INOUT) :: sb_atm_name(:,:)
        INTEGER                                      :: s, i

        OPEN(UNIT=20,FILE=sb_file_pos,STATUS='old',FORM='formatted',ACTION='READ')
        DO s = 1, sb_nb_step
            READ(20, *)
            READ(20, *)
            DO i=1,sb_nb_atm
                sb_atm_mat(2,i,s) = -1
                READ(20, *) sb_atm_name(i,s), sb_atm_mat(4,i,s), sb_atm_mat(5,i,s), sb_atm_mat(6,i,s)
                sb_atm_mat(1,i,s) = i
                IF (sb_atm_name(i,s) .EQ. "CC") THEN !  Carbon (3C)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 1
                !ELSE IF (sb_atm_name(i,s) .EQ. "CP") THEN ! Carbon (2C)
                !    sb_atm_mat(2,i,s) = 12
                !    sb_atm_mat(3,i,s) = 2
                ELSE IF (sb_atm_name(i,s) .EQ. "C3A") THEN ! Alkoxy Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 3
                ELSE IF (sb_atm_name(i,s) .EQ. "C3O") THEN ! PAlcohol/Alcohol Carbon (3C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 4
                !ELSE IF (sb_atm_name(i,s) .EQ. "C3O") THEN ! P Alcohol Carbon
                !    sb_atm_mat(2,i,s) = 12
                !    sb_atm_mat(3,i,s) = 5
                ELSE IF (sb_atm_name(i,s) .EQ. "C3E") THEN ! Epoxy Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 6
                ELSE IF (sb_atm_name(i,s) .EQ. "C3P") THEN ! Protonated Epoxy Carbon (3C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 7
                ELSE IF (sb_atm_name(i,s) .EQ. "C3F") THEN ! Ether Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 8
                ELSE IF (sb_atm_name(i,s) .EQ. "C3Q") THEN ! Protonated Ether Carbon (3C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 9
                ELSE IF (sb_atm_name(i,s) .EQ. "C2A") THEN ! Ketone/Alkoxy Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 10
                ELSE IF (sb_atm_name(i,s) .EQ. "C2O") THEN ! PKetone/Alcohol/PAlcohol Carbon (2C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 11
                ELSE IF (sb_atm_name(i,s) .EQ. "C2E") THEN ! Epoxy Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 12
                ELSE IF (sb_atm_name(i,s) .EQ. "C2P") THEN ! Protonated Epoxy Carbon (2C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 13
                ELSE IF (sb_atm_name(i,s) .EQ. "C2F") THEN ! Ether Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 14
                ELSE IF (sb_atm_name(i,s) .EQ. "C2Q") THEN ! Protonated Ether Carbon (2C+1O(Hx))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 15
                !!!---
                ELSE IF ((sb_atm_name(i,s) .EQ. "OH2") .OR.& ! PAlcohol/Alcohol Oxygen (1C+xH)
                    (sb_atm_name(i,s) .EQ. "OH3")) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 16 
                ELSE IF (sb_atm_name(i,s) .EQ. "OEN") THEN ! Epoxy Oxygen (2C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 17
                ELSE IF (sb_atm_name(i,s) .EQ. "OFN") THEN ! Ether Oxygen (2C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 18
                ELSE IF (sb_atm_name(i,s) .EQ. "OEP") THEN ! Protonated Epoxy Oxygen (2C+xH)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 19
                ELSE IF (sb_atm_name(i,s) .EQ. "OFP") THEN ! Protonated Ether Oxygen (2C+xH)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 20
                ELSE IF (sb_atm_name(i,s) .EQ. "OA3") THEN ! Alkoxy Oxygen (1C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 21
                ELSE IF (sb_atm_name(i,s) .EQ. "OA2") THEN ! Ketone/Alkoxy Oxygen (1C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 22
                !!!---
                ELSE IF (sb_atm_name(i,s) .EQ. "OW") THEN ! H2O Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 23
                ELSE IF (sb_atm_name(i,s) .EQ. "OM") THEN ! OH Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 24
                ELSE IF (sb_atm_name(i,s) .EQ. "OP") THEN ! H3O Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 25
                !!!---
                ELSE IF (sb_atm_name(i,s) .EQ. "HW") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 26
                ELSE IF (sb_atm_name(i,s) .EQ. "HM") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 27
                ELSE IF (sb_atm_name(i,s) .EQ. "HP") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 28
                ELSE IF (sb_atm_name(i,s) .EQ. "HO") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 29
                ELSE IF (sb_atm_name(i,s) .EQ. "HEP") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 30
                ELSE IF (sb_atm_name(i,s) .EQ. "HET") THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 31
                ELSE IF (sb_atm_name(i,s) .EQ. "NA") THEN
                    sb_atm_mat(2,i,s) = 11
                    sb_atm_mat(3,i,s) = 32
                ELSE IF (sb_atm_name(i,s) .EQ. "CLM") THEN
                    sb_atm_mat(2,i,s) = 17
                    sb_atm_mat(3,i,s) = 33
                ! Unknown carbon types
                ELSE IF ((sb_atm_name(i,s) .EQ. "C").OR.&
                    (sb_atm_name(i,s) .EQ. "C2").OR.&
                    (sb_atm_name(i,s) .EQ. "C3").OR.&
                    (sb_atm_name(i,s) .EQ. "CX")) THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ((sb_atm_name(i,s) .EQ. "O") .OR.&
                    (sb_atm_name(i,s) .EQ. "OH") .OR. &
                    (sb_atm_name(i,s) .EQ. "OE") .OR. &
                    (sb_atm_name(i,s) .EQ. "OA") .OR. &
                    (sb_atm_name(i,s) .EQ. "OX")) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ((sb_atm_name(i,s) .EQ. "H") .OR.&
                    (sb_atm_name(i,s) .EQ. "H1") .OR.&
                    (sb_atm_name(i,s) .EQ. "H2") .OR.&
                    (sb_atm_name(i,s) .EQ. "H3") .OR.&
                    (sb_atm_name(i,s) .EQ. "HX")) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = -1
                END IF
            END DO
        END DO
        CLOSE(UNIT=20)
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