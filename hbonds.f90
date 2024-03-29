!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM hbonds
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

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), OHvec_mat(:,:,:), is_mat(:,:,:), OHvec_hbond_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
INTEGER                         :: nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: OpH_disp_vec(3), OpH_disp_norm, OH_disp_vec(3), OH_disp_norm, alpha
REAL(dp)                        :: oHpO_disp_vec(3), oHpO_disp_norm
REAL(dp)                        :: XOh_disp_norm, XHo_disp_norm, XO_disp_norm
REAL(dp)                        :: SpOh_disp_norm, SpO_disp_norm

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, Udonnor_count, Udonnor_count2
INTEGER, ALLOCATABLE            :: nb_max_OHvec(:)

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, j, o, l, n
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching H-Bonds'
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
PRINT'(A100)', 'Run, H-Bonds, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(30,nb_atm,nb_step))
ALLOCATE(atm_name(nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)
DEALLOCATE(atm_name) ! Not Used

nb_o = COUNT( atm_mat(2,:,1) .EQ. 16, DIM=1 )

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

! B ----------------------------------------------- OH groups and corresponding values
start = OMP_get_wtime()
ALLOCATE(OHvec_mat(44,nb_o*3,nb_step))
ALLOCATE(nb_max_OHvec(nb_step))

nb_max_OHvec(:) = 0.0
OHvec_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_max_OHvec, nb_step, nb_atm)&
!$OMP SHARED(hb_OpH_rcut)&
!$OMP PRIVATE(s, i, j, o)&
!$OMP PRIVATE(OpH_disp_norm,OpH_disp_vec)
DO s = 1, nb_step
    o = 0
 OC:DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
            IF ( ( atm_mat(3,i,s) .NE. 34 ) .AND. ( atm_mat(3,i,s) .EQ. 33 ) .AND. &
            ( atm_mat(3,i,s) .NE. 35 ) .AND. ( atm_mat(3,i,s) .EQ. 25 ) .AND. &
            ( atm_mat(3,i,s) .NE. 26 ) .AND. ( atm_mat(3,i,s) .EQ. 28 ) .AND. &
            ( atm_mat(3,i,s) .NE. 30 ) ) THEN
                CYCLE OC
            END IF
         HC:DO j = 1, nb_atm
                IF ( atm_mat(2,j,s) .EQ. 1 ) THEN
                    IF ( ( atm_mat(3,i,s) .EQ. 34 ) .AND. ( atm_mat(3,j,s) .EQ. 37 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 33 ) .AND. ( atm_mat(3,j,s) .EQ. 36 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 35 ) .AND. ( atm_mat(3,j,s) .EQ. 38 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 25 ) .AND. ( atm_mat(3,j,s) .EQ. 39 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 26 ) .AND. ( atm_mat(3,j,s) .EQ. 40 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 28 ) .AND. ( atm_mat(3,j,s) .EQ. 41 ) ) THEN
                        CONTINUE
                    ELSE IF  ( ( atm_mat(3,i,s) .EQ. 30 ) .AND. ( atm_mat(3,j,s) .EQ. 42 ) ) THEN
                        CONTINUE
                    ELSE
                        CYCLE HC
                    END IF

                    CALL sb_dist(atm_mat(4:6,i,s),atm_mat(4:6,j,s),box,norm_ij=OpH_disp_norm,vec_ij=OpH_disp_vec)

                    IF ( OpH_disp_norm .LT. hb_OpH_rcut ) THEN
                        o = o + 1
                        OHvec_mat(1,o,s) = atm_mat(1,i,s)
                        OHvec_mat(2,o,s) = atm_mat(1,j,s)
                        OHvec_mat(3,o,s) = atm_mat(3,i,s)
                        OHvec_mat(4,o,s) = atm_mat(3,j,s)
                        OHvec_mat(5:7,o,s) = OpH_disp_vec  ! Disp
                        OHvec_mat(8:10,o,s) = atm_mat(4:6,i,s) ! O pos
                        OHvec_mat(11:13,o,s) = atm_mat(4:6,j,s) ! H pos
                    END IF
                END IF
            END DO HC
        END IF
    END DO OC
    nb_max_OHvec(s) = COUNT( OHvec_mat(1,:,s) .NE. 0, DIM=1 )
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH groups:", finish-start, "seconds elapsed"

! C ----------------------------------------------- OH/O hbonds
start = OMP_get_wtime()
ALLOCATE(OHvec_hbond_mat(20,nb_o*3,nb_step))
OHvec_hbond_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_atm, nb_step, OHvec_hbond_mat)&
!$OMP SHARED(hb_oHpO_rcut)&
!$OMP PRIVATE(s, i, j, l, n)&
!$OMP PRIVATE(oHpO_disp_vec, oHpO_disp_norm, Udonnor_count, Udonnor_count2, OH_disp_norm, OH_disp_vec)&
!$OMP PRIVATE(alpha)
DO s = 1, nb_step
 C1:DO i = 1, nb_o*3
        IF ( OHvec_mat(1,i,s) .EQ. 0 ) THEN
            CYCLE C1
        END IF
        Udonnor_count = 0
        Udonnor_count2 = 0
        n = 4
        OHvec_hbond_mat(1,i,s) = OHvec_mat(1,i,s)
        OHvec_hbond_mat(2,i,s) = OHvec_mat(3,i,s)
        OHvec_hbond_mat(3,i,s) = OHvec_mat(2,i,s)
        DO j = 1, nb_atm
            IF ( atm_mat(2,j,s) .EQ. 16 ) THEN

                CALL sb_dist(OHvec_mat(11:13,i,s),atm_mat(4:6,j,s),box,norm_ij=oHpO_disp_norm,vec_ij=oHpO_disp_vec) ! HO hbonds vector
                CALL sb_dist(OHvec_mat(8:10,i,s),OHvec_mat(11:13,i,s),box,norm_ij=OH_disp_norm,vec_ij=OH_disp_vec) ! OH vector

                IF( ( oHpO_disp_norm .LE. hb_oHpO_rcut ) .AND. ( atm_mat(1,j,s) .NE. OHvec_mat(1,i,s) ) ) THEN
                    atm_mat(7,j,s) = atm_mat(7,j,s) + 1 ! Acceptor count (O)
                    OHvec_mat(15,i,s) = OHvec_mat(15,i,s) + 1 ! Donnor count (OH)
                    IF ( n .GT. 23 ) THEN
                        PRINT*, "HBONDS OVERFLOW", OHvec_mat(1,i,s)
                    END IF

                    alpha = ACOS(DOT_PRODUCT(oHpO_disp_vec(:), OH_disp_vec(:)) / (oHpO_disp_norm * OH_disp_norm) ) ! alpha
                    OHvec_hbond_mat(n,i,s) = atm_mat(1,j,s)
                    OHvec_hbond_mat(n+1,i,s) = atm_mat(3,j,s)
                    OHvec_hbond_mat(n+2,i,s) = oHpO_disp_norm
                    OHvec_hbond_mat(n+3,i,s) = alpha
                    n = n + 4
                    IF ( Udonnor_count .EQ. 0 ) THEN
                        OHvec_mat(16,i,s) = OHvec_mat(16,i,s) + 1 ! Unique donnor count (OH)
                        Udonnor_count = 1
                    END IF
                C21:DO l=1, nb_atm
                        IF ( atm_mat(1,l,s) .EQ. OHvec_mat(1,i,s) ) THEN
                            atm_mat(8,l,s) = atm_mat(8,l,s) + 1 ! Donnor count (O)
                            IF ( Udonnor_count2 .EQ. 0 ) THEN
                                atm_mat(9,l,s) = atm_mat(9,l,s) + 1 ! Unique donnor count (O)
                                Udonnor_count2 = 1
                            END IF
                            EXIT C21
                        END IF
                    END DO C21
                C22:DO l=1,nb_o*3
                        IF ( OHvec_mat(1,l,s) .EQ. 0 ) THEN
                            CYCLE C22
                        ELSE IF ( OHvec_mat(1,l,s) .EQ. atm_mat(1,j,s) ) THEN
                            OHvec_mat(14,l,s) = OHvec_mat(14,l,s) + 1 ! Acceptor count (OH)
                        END IF
                    END DO C22
                END IF
            END IF
        END DO
    END DO C1
  C3:DO i = 1, (nb_o*3)
        IF ( OHvec_mat(1,i,s) .EQ. 0 ) THEN
            CYCLE C3
        END IF
        DO j = 1, nb_atm
            IF ( OHvec_mat(1,i,s) .EQ. atm_mat(1,j,s) ) THEN
                OHvec_mat(17,i,s) = atm_mat(7,j,s)
                OHvec_mat(18,i,s) = atm_mat(8,j,s)
                OHvec_mat(19,i,s) = atm_mat(9,j,s)
            END IF
        END DO
    END DO C3
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH/O hbonds:", finish-start, "seconds elapsed"

! D ----------------------------------------------- Proximity between functionnal groups and any OH group
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_step, nb_atm)&
!$OMP SHARED(hb_X1Oh_rcut, hb_X2Oh_rcut)&
!$OMP PRIVATE(s, i, j)&
!$OMP PRIVATE(XOh_disp_norm, XHo_disp_norm)
DO s = 1, nb_step
 D1:DO i = 1, nb_o*3
        IF ( OHvec_mat(1,i,s) .EQ. 0 ) THEN
            CYCLE D1
        END IF
     D2:DO j = 1, nb_atm
            IF ( atm_mat(3,j,s) .EQ. -1 ) THEN
                CYCLE D2
            END IF
            IF ( atm_mat(1,j,s) .EQ. OHvec_mat(1,i,s) ) THEN
                CYCLE D2
            END IF

            CALL sb_dist(OHvec_mat(8:10,i,s),atm_mat(4:6,j,s),box,norm_ij=XOh_disp_norm)
            CALL sb_dist(OHvec_mat(11:13,i,s),atm_mat(4:6,j,s),box,norm_ij=XHo_disp_norm)

            IF ( ( (XOh_disp_norm .LE. OHvec_mat(32,i,s) ) .OR. &
            (OHvec_mat(32,i,s) .EQ. 0.0) ) .AND. &
            (atm_mat(2,j,s) .EQ. 12) ) THEN
                OHvec_mat(31,i,s) = atm_mat(1,j,s)
                OHvec_mat(32,i,s) = XOh_disp_norm
                OHvec_mat(33,i,s) = atm_mat(6,j,s)
            END IF

            IF ( ( (XHo_disp_norm .LE. OHvec_mat(36,i,s) ) .OR. &
            (OHvec_mat(36,i,s) .EQ. 0.0) ) .AND. &
            (atm_mat(2,j,s) .EQ. 12) ) THEN
                OHvec_mat(35,i,s) = atm_mat(1,j,s)
                OHvec_mat(36,i,s) = XHo_disp_norm
                OHvec_mat(37,i,s) = atm_mat(6,j,s)
            END IF

            IF ( XOh_disp_norm .LE. hb_X1Oh_rcut ) THEN
                IF ( atm_mat(2,j,s) .EQ. 12 ) THEN ! Close to any carbon
                    OHvec_mat(20,i,s) = 1
                    OHvec_mat(24,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 27) .OR. & ! Close to an ether oxygen (protonated or not)
                (atm_mat(3,j,s) .EQ. 28) ) THEN ! OEP
                    OHvec_mat(21,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 29) .OR. & ! Close to an epoxy oxygen (protonated or not)
                (atm_mat(3,j,s) .EQ. 30) ) THEN ! OFP
                    OHvec_mat(39,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 25 ) THEN ! Close to an alcohol oxygen
                    OHvec_mat(22,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 26 ) THEN ! Close to a protonated alcohol oxygen
                    OHvec_mat(44,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 31 ) THEN ! OA3 Alkoxy
                    OHvec_mat(23,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 32 ) THEN ! OA2/OA1 Ketone/Alkoxy
                    OHvec_mat(40,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 60) .OR. & ! NA
                (atm_mat(3,j,s) .EQ. 61) ) THEN ! CLM
                    OHvec_mat(41,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 33 ) THEN ! OH
                    OHvec_mat(42,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 35 ) THEN ! H3O
                    OHvec_mat(43,i,s) = 1
                END IF
            ELSE IF ( XOh_disp_norm .LE. hb_X2Oh_rcut ) THEN
                IF ( atm_mat(2,j,s) .EQ. 12 ) THEN ! Close to any carbon
                    OHvec_mat(24,i,s) = 1
                END IF
            END IF
        END DO D2
        IF ( OHvec_mat(33,i,s) .GT. OHvec_mat(10,i,s) ) THEN
            OHvec_mat(34,i,s) = -1.0
        ELSE
            OHvec_mat(34,i,s) = 1.0
        END IF
        IF ( OHvec_mat(37,i,s) .GT. OHvec_mat(13,i,s) ) THEN
            OHvec_mat(38,i,s) = -1.0
        ELSE
            OHvec_mat(38,i,s) = 1.0
        END IF
    END DO D1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and OH group:", finish-start, "seconds elapsed"

! E ----------------------------------------------- Proximity between functionnal groups and any O atom
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(hb_X1O_rcut, hb_X2O_rcut)&
!$OMP PRIVATE(s, i, j)&
!$OMP PRIVATE(XO_disp_norm)
DO s = 1, nb_step
    E1:DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .NE. 16 ) THEN
            CYCLE E1
        END IF
        E2:DO j = 1, nb_atm
            IF ( atm_mat(3,j,s) .EQ. -1 ) THEN
                CYCLE E2
            END IF
            IF ( atm_mat(1,j,s) .EQ. atm_mat(1,i,s) ) THEN
                CYCLE E2
            END IF

            CALL sb_dist(atm_mat(4:6,i,s),atm_mat(4:6,j,s),box,norm_ij=XO_disp_norm)

            IF ( ( (XO_disp_norm .LE. atm_mat(22,i,s) ) .OR. &
            ( atm_mat(22,i,s) .EQ. 0.0) ) .AND.&
            ( atm_mat(2,j,s) .EQ. 12) ) THEN
                atm_mat(21,i,s) = atm_mat(1,j,s)
                atm_mat(22,i,s) = XO_disp_norm
                atm_mat(23,i,s) = atm_mat(6,j,s)
            END IF
            IF ( XO_disp_norm .LE. hb_X1O_rcut ) THEN
                IF ( atm_mat(2,j,s) .EQ. 12 ) THEN ! Close to any carbon
                    atm_mat(10,i,s) = 1
                    atm_mat(14,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 27) .OR. & ! Close to an ether oxygen (protonated or not)
                ( atm_mat(3,j,s) .EQ. 28) ) THEN
                    atm_mat(11,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 29) .OR. & ! Close to an epoxy oxygen (protonated or not)
                ( atm_mat(3,j,s) .EQ. 30) ) THEN
                    atm_mat(25,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 25 ) THEN ! Close to an alcohol oxygen
                    atm_mat(12,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 26 ) THEN ! Close to a protonated alcohol oxygen
                    atm_mat(30,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 31 ) THEN ! OA3 Alkoxy
                    atm_mat(13,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 32 ) THEN ! OA2/OA1 Ketone/Alkoxy
                    atm_mat(26,i,s) = 1
                ELSE IF ( (atm_mat(3,j,s) .EQ. 60) .OR. & ! NA
                ( atm_mat(3,j,s) .EQ. 61) ) THEN ! CLM
                    atm_mat(27,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 33 ) THEN ! OH
                    atm_mat(28,i,s) = 1
                ELSE IF ( atm_mat(3,j,s) .EQ. 35 ) THEN ! H3O
                    atm_mat(29,i,s) = 1
                END IF
            ELSE IF ( XO_disp_norm .LE. hb_X2O_rcut ) THEN
                IF ( atm_mat(2,j,s) .EQ. 12 ) THEN ! Close to any carbon
                    atm_mat(14,i,s) = 1
                END IF
            END IF
        END DO E2
        IF ( atm_mat(23,i,s) .GT. atm_mat(6,i,s) ) THEN
            atm_mat(24,i,s) = -1.0
        ELSE
            atm_mat(24,i,s) = 1.0
        END IF
    END DO E1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and O atoms:", finish-start, "seconds elapsed"

IF ( IS_c .EQ. 'Y' ) THEN

    ! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, is_mat, nb_is, nb_step)&
    !$OMP PRIVATE(s, i, j)&
    !$OMP PRIVATE(SpOh_disp_norm)
    DO s = 1, nb_step
        F2:DO i = 1, nb_o*3
            IF ( OHvec_mat(1,i,s) .EQ. 0 ) THEN
                CYCLE F2
            END IF
            DO j = 1, nb_is(s)
                IF ( is_mat(4,j,s) .EQ. 1 ) THEN

                    CALL sb_dist(OHvec_mat(8:10,i,s),is_mat(1:3,j,s),box,norm_ij=SpOh_disp_norm)

                    IF ( ( SpOh_disp_norm .LT. OHvec_mat(25,i,s) ) .OR. ( OHvec_mat(25,i,s) .EQ. 0.0 ) ) THEN
                        OHvec_mat(25,i,s) = SpOh_disp_norm
                        IF ( OHvec_mat(10,i,s) .LT. is_mat(3,j,s) ) THEN
                            OHvec_mat(26,i,s) = -1
                        ELSE
                            OHvec_mat(26,i,s) = 1
                        END IF
                        OHvec_mat(27,i,s) = is_mat(5,j,s)
                    END IF
                ELSE IF ( is_mat(4,j,s) .EQ. 2 ) THEN

                    CALL sb_dist(OHvec_mat(8:10,i,s),is_mat(1:3,j,s),box,norm_ij=SpOh_disp_norm)

                    IF ( ( SpOh_disp_norm .LT. OHvec_mat(28,i,s) ) .OR. ( OHvec_mat(28,i,s) .EQ. 0.0 ) ) THEN
                        OHvec_mat(28,i,s) = SpOh_disp_norm
                        IF ( OHvec_mat(10,i,s) .GT. is_mat(3,j,s) ) THEN
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

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_o, is_mat, nb_is, nb_atm, nb_step)&
    !$OMP PRIVATE(s, i, j)&
    !$OMP PRIVATE(SpO_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
                DO j = 1, nb_is(s)
                    IF ( is_mat(4,j,s) .EQ. 1 ) THEN

                        CALL sb_dist(atm_mat(4:6,i,s),is_mat(1:3,j,s),box,norm_ij=SpO_disp_norm)

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

                        CALL sb_dist(atm_mat(4:6,i,s),is_mat(1:3,j,s),box,norm_ij=SpO_disp_norm)

                        IF ( (SpO_disp_norm .LT. atm_mat(18,i,s) ) .OR. ( atm_mat(18,i,s) .EQ. 0.0 ) ) THEN
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
    PRINT'(A40,F14.2,A20)', "Proximity IS and O atoms:", finish-start, "seconds elapsed"

END IF

! X ----------------------------------------------- Write OH BOND
start = OMP_get_wtime()

OPEN(UNIT=31, FILE=suffix//"_O_HB.txt")
WRITE(31, '(A4,1X,A10,1X,A10,1X,A6&
&,1X,A6,1X,A6,1X,A6&
&,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4&
&,1X,A14,1X,A14,1X,A14)')&
    "Traj", "Step", "O_ID", "O_Type"&
    , "TA_O", "TD_O", "TDU_O"&
    , "cCC", "cOEP", "cOET", "cOH", "cOHP", "cOA3", "cOA", "cION", "cOM", "cH3O", "cCX"&
    , "dist_ISD", "dist_ISU", "dist_AS"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
            WRITE(31,'(A4,1X,I10,1X,I10,1X,I6&
            &,1X,I6,1X,I6,1X,I6&
            &,1X,I4,1X,I4,1X,I4&
            &,1X,I4,1X,I4&
            &,1X,I4,1X,I4&
            &,1X,I4,1X,I4,1X,I4,1X,I4&
            &,1X,E14.5,1X,E14.5,1X,E14.5)')&
            suffix, s, INT(atm_mat(1,i,s)), INT(atm_mat(3,i,s))&
            , INT(atm_mat(7,i,s)), INT(atm_mat(8,i,s)), INT(atm_mat(9,i,s))&
            , INT(atm_mat(10,i,s)), INT(atm_mat(11,i,s)), INT(atm_mat(25,i,s))& ! cCC cOEP cOET
            , INT(atm_mat(12,i,s)), INT(atm_mat(30,i,s))& ! cOH cOH
            , INT(atm_mat(13,i,s)), INT(atm_mat(26,i,s))& ! cOA3 cOA
            , INT(atm_mat(27,i,s)), INT(atm_mat(28,i,s)), INT(atm_mat(29,i,s))& ! cION, cOM, cH3O
            , INT(atm_mat(14,i,s))& ! cCX
            , (atm_mat(15,i,s)*atm_mat(16,i,s)),(atm_mat(18,i,s)*atm_mat(19,i,s))&
            , (atm_mat(22,i,s)*atm_mat(24,i,s))
    END IF
    END DO
END DO
CLOSE(UNIT=31)

OPEN(UNIT=32, FILE = suffix//"_OH_HB.txt")
WRITE(32, '(A4,1X,A10,1X,A10,1X,A10,1X,A6&
&,1X,A6,1X,A6,1X,A6&
&,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4&
&,1X,A6,1X,A6,1X,A6,1X,A14,1X,A14&
&,1X,A14,1X,A14&
&,1X,A6,1X,A6,1X,A5,1X,A6&
&,1X,A6,1X,A6,1X,A5,1X,A6&
&,1X,A6,1X,A6,1X,A5,1X,A6&
&,1X,A6,1X,A6,1X,A5,1X,A6)')&
    "Traj", "Step", "O_ID", "H_ID", "O_Type"&
    , "TA_OH", "TD_OH", "TDU_OH"&
    , "cCC", "cOEP", "cOET", "cOH", "cOHP", "cOA3", "cOA", "cION", "cOM", "cH3O", "cCX"&
    , "TA_O", "TD_O", "TDU_O", "dist_ISD", "dist_ISU"&
    , "dist_AS", "dist_HAS"&
    , "O_ID", "O_Type", "r (A)", "alpha (rad)"&
    , "O_ID", "O_Type", "r (A)", "alpha (rad)"&
    , "O_ID", "O_Type", "r (A)", "alpha (rad)"&
    , "O_ID", "O_Type", "r (A)", "alpha (rad)"
DO s = 1, nb_step
    DO i = 1, nb_max_OHvec(s)
        WRITE(32,'(A4,1X,I10,1X,I10,1X,I10,1X,I6&
        &,1X,I6,1X,I6,1X,I6&
        &,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4&
        &,1X,I6,1X,I6,1X,I6,1X,E14.5,1X,E14.5&
        &,1X,E14.5,1X,E14.5)', ADVANCE = "no")&
        suffix, s, INT(OHvec_mat(1,i,s)), INT(OHvec_mat(2,i,s)), INT(OHvec_mat(3,i,s)) &
        , INT(OHvec_mat(14,i,s)), INT(OHvec_mat(15,i,s)), INT(OHvec_mat(16,i,s))&
        , INT(OHvec_mat(20,i,s)), INT(OHvec_mat(21,i,s)), INT(OHvec_mat(39,i,s))& ! cCC cOEP cOET
        , INT(OHvec_mat(22,i,s)), INT(OHvec_mat(44,i,s))& ! cOH cOH
        , INT(OHvec_mat(23,i,s)), INT(OHvec_mat(40,i,s))& ! cOA3 cOA
        , INT(OHvec_mat(41,i,s)), INT(OHvec_mat(42,i,s)), INT(OHvec_mat(43,i,s))& ! cION, cOM, cH3O
        , INT(OHvec_mat(24,i,s))& ! cCX
        , INT(OHvec_mat(17,i,s)), INT(OHvec_mat(18,i,s)), INT(OHvec_mat(19,i,s))&
        , (OHvec_mat(25,i,s)*OHvec_mat(26,i,s)), (OHvec_mat(28,i,s)*OHvec_mat(29,i,s))&
        , (OHvec_mat(32,i,s)*OHvec_mat(34,i,s)), (OHvec_mat(36,i,s)*OHvec_mat(38,i,s))
        DO j = 4, 19, 4
            WRITE(32, '(1X,I6,1X,I6,1X,F5.3,1X,F6.4)', ADVANCE = "no")&
            INT(OHvec_hbond_mat(j,i,s)),INT(OHvec_hbond_mat(j+1,i,s)), OHvec_hbond_mat(j+2,i,s), OHvec_hbond_mat(j+3,i,s)
        END DO
        WRITE(32,'()')
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
IF ( IS_c .EQ. 'Y') DEALLOCATE(is_mat,nb_is)

DEALLOCATE(OHvec_mat,atm_mat,nb_max_OHvec)

END PROGRAM hbonds