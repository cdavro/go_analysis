!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003

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
        !-- Procedure only
        INTEGER                                      :: s, i

        OPEN(UNIT=20,FILE=sb_file_pos,STATUS='old',FORM='formatted',ACTION='READ')
        DO s = 1, sb_nb_step
            READ(20, *)
            READ(20, *)
            DO i=1,sb_nb_atm
                sb_atm_mat(2,i,s) = -1
                READ(20, *) sb_atm_name(i,s), sb_atm_mat(4,i,s), sb_atm_mat(5,i,s), sb_atm_mat(6,i,s)
                sb_atm_mat(1,i,s) = i
                IF ( sb_atm_name(i,s) .EQ. "C3" ) THEN !  Carbon (3C)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 1
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2" ) THEN ! Carbon (2C)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 2
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1" ) THEN ! Carbon (1C)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 3
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3A" ) THEN ! Alkoxy Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 4
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3B" ) THEN ! Alcohol Carbon (3C+1O(1H))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 5
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3C" ) THEN ! P Alcohol Carbon (3C+1O(2H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 6
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3D" ) THEN ! Ether Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 7
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3E" ) THEN ! Protonated Ether Carbon (3C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 8
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3F" ) THEN ! Epoxy Carbon (3C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 9
                ELSE IF ( sb_atm_name(i,s) .EQ. "C3G" ) THEN ! Protonated Epoxy Carbon (3C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 10
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2A" ) THEN ! Ketone/Alkoxy Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 11
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2B" ) THEN ! PKetone/Alcohol Carbon (2C+1O(1H))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 12
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2C" ) THEN ! PAlcohol Carbon (2C+1O(2H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 13
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2D" ) THEN ! Ether Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 14
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2E" ) THEN ! Protonated Ether Carbon (2C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 15
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2F" ) THEN ! Epoxy Carbon (2C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 16
                ELSE IF ( sb_atm_name(i,s) .EQ. "C2G" ) THEN ! Protonated Epoxy Carbon (2C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 17
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1A" ) THEN ! Ketone/Alkoxy Carbon (1C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 18
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1B" ) THEN ! PKetone/Alcohol Carbon (1C+1O(1H))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 19
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1C" ) THEN ! PAlcohol Carbon (1C+1O(2H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 20
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1D" ) THEN ! Ether Carbon (1C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 21
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1E" ) THEN ! Protonated Ether Carbon (1C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 22
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1F" ) THEN ! Epoxy Carbon (1C+1O)
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 23
                ELSE IF ( sb_atm_name(i,s) .EQ. "C1G" ) THEN ! Protonated Epoxy Carbon (1C+1O(1H+))
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = 24
                !!!----
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OB3" ) .OR. & ! Alcohol Oxygen (1C+1H)
                (sb_atm_name(i,s) .EQ. "OB2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OB1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 25
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OC3" ) .OR. & ! PAlcohol Oxygen (1C+2H+)
                (sb_atm_name(i,s) .EQ. "OC2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OC1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 26
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OD3" ) .OR.& ! Ether Oxygen (1C+0H)
                (sb_atm_name(i,s) .EQ. "OD2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OD1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 27
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OE3" ) .OR.& ! PEther Oxygen (1C+1H+)
                (sb_atm_name(i,s) .EQ. "OE2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OE1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 28
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OF3" ) .OR.& ! Epoxy Oxygen (1C+0H)
                (sb_atm_name(i,s) .EQ. "OF2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OF1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 29
                ELSE IF ( (sb_atm_name(i,s) .EQ. "OG3" ) .OR.& ! PEpoxy Oxygen (1C+1H+)
                (sb_atm_name(i,s) .EQ. "OG2" ) .OR. &
                (sb_atm_name(i,s) .EQ. "OG1" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 30
                ELSE IF ( sb_atm_name(i,s) .EQ. "OA3" ) THEN ! Alkoxy Oxygen (1C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 31
                ELSE IF ( ( sb_atm_name(i,s) .EQ. "OA2" ) .OR. &
                ( sb_atm_name(i,s) .EQ. "OA1" ) ) THEN ! Ketone/Alkoxy Oxygen (1C+0H)
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 32
                !!!---
                ELSE IF ( sb_atm_name(i,s) .EQ. "OM" ) THEN ! OH Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 33
                ELSE IF ( sb_atm_name(i,s) .EQ. "OW" ) THEN ! H2O Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 34
                ELSE IF ( sb_atm_name(i,s) .EQ. "OP" ) THEN ! H3O Oxygen
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = 35
                !!!---
                ELSE IF ( sb_atm_name(i,s) .EQ. "HM" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 36
                ELSE IF ( sb_atm_name(i,s) .EQ. "HW" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 37
                ELSE IF ( sb_atm_name(i,s) .EQ. "HP" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 38
                ELSE IF ( sb_atm_name(i,s) .EQ. "HB" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 39
                ELSE IF ( sb_atm_name(i,s) .EQ. "HC" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 40
                ELSE IF ( sb_atm_name(i,s) .EQ. "HE" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 41
                ELSE IF ( sb_atm_name(i,s) .EQ. "HG" ) THEN
                    sb_atm_mat(2,i,s) = 1
                    sb_atm_mat(3,i,s) = 42
                !!!---
                ELSE IF ( sb_atm_name(i,s) .EQ. "NA" ) THEN
                    sb_atm_mat(2,i,s) = 11
                    sb_atm_mat(3,i,s) = 60
                ELSE IF ( sb_atm_name(i,s) .EQ. "CLM" ) THEN
                    sb_atm_mat(2,i,s) = 17
                    sb_atm_mat(3,i,s) = 61
                ! Unknown carbon types
                ELSE IF ( (sb_atm_name(i,s) .EQ. "C" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "CX" ) ) THEN
                    sb_atm_mat(2,i,s) = 12
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ( (sb_atm_name(i,s) .EQ. "O" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "OA" ) .OR. &
                    (sb_atm_name(i,s) .EQ. "OB" ) .OR. &
                    (sb_atm_name(i,s) .EQ. "OC" ) .OR. &
                    (sb_atm_name(i,s) .EQ. "OD" ) .OR. &
                    (sb_atm_name(i,s) .EQ. "OE" ) .OR. &
                    (sb_atm_name(i,s) .EQ. "OX" ) ) THEN
                    sb_atm_mat(2,i,s) = 16
                    sb_atm_mat(3,i,s) = -1
                ELSE IF ( (sb_atm_name(i,s) .EQ. "H" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "H1" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "H2" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "H3" ) .OR.&
                    (sb_atm_name(i,s) .EQ. "HX" ) ) THEN
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
            IF ( sb_iostatus .NE. 0 ) THEN
                EXIT
            ELSE
                sb_nb_line_is = sb_nb_line_is + 1
            END IF
        END DO
        REWIND(21)
        sb_nb_max_is = CEILING( 1.0 * sb_nb_line_is / sb_nb_step ) * 2
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
                IF ( sb_dummy_char .EQ. "XD" ) THEN
                    sb_is_mat(4,i,s) = 1
                ELSE IF ( sb_dummy_char .EQ. "XU" ) THEN
                    sb_is_mat(4,i,s) = 2
                ELSE
                    sb_is_mat(4,i,s) = 0
                END IF
                DO k = 1, 3
                    sb_is_mat(k,i,s) = sb_is_mat(k,i,s) - sb_box(k) * ANINT( sb_is_mat(k,i,s) / sb_box(k) )
                END DO
            END DO
        END DO
        CLOSE(UNIT=21)

    END SUBROUTINE sb_read_is

    SUBROUTINE sb_dist(xyz_i,xyz_j,box_dim,norm_ij,vec_ij)
        REAL(dp), INTENT(IN) :: xyz_i(:), xyz_j(:), box_dim(:)
        REAL(dp) :: sbo_norm_ij
        REAL(dp) :: sbo_vec_ij(3)
        REAL(dp), INTENT(OUT), OPTIONAL :: norm_ij
        REAL(dp), INTENT(OUT), OPTIONAL :: vec_ij(:)
        sbo_vec_ij = xyz_j - xyz_i
        sbo_vec_ij = sbo_vec_ij - box_dim * ANINT( sbo_vec_ij / box_dim )
        sbo_norm_ij = NORM2 ( sbo_vec_ij )
        IF ( PRESENT( norm_ij)  ) norm_ij = sbo_norm_ij
        IF ( PRESENT( vec_ij ) ) vec_ij = sbo_vec_ij
    END SUBROUTINE sb_dist

    PURE FUNCTION fc_itoc(i) result(res)
        INTEGER, INTENT(IN) :: i
        CHARACTER(:), ALLOCATABLE :: res
        CHARACTER( RANGE( (i)+2 ) ) :: tmp
        WRITE(tmp,'(i0)') i
        res = TRIM( tmp )
    END FUNCTION

    PURE FUNCTION fc_trim_ext(string_in) result(string_out)
        CHARACTER(*), INTENT(IN) :: string_in
        CHARACTER(:), ALLOCATABLE :: string_out
        string_out = TRIM( string_in(1:SCAN( TRIM( string_in ),".", .TRUE. )-1) )
    END FUNCTION

END MODULE SB_GO