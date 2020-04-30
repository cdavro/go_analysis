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
INTEGER                         :: nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: OpH_disp_vec(3), OpH_disp_norm, OH_disp_vec(3), OH_disp_norm, alpha
REAL(dp)                        :: oHpO_disp_vec(3), oHpO_disp_norm
REAL(dp)                        :: XOh_disp_vec(3), XOh_disp_norm
REAL(dp)                        :: XHo_disp_vec(3), XHo_disp_norm
REAL(dp)                        :: XO_disp_vec(3), XO_disp_norm
REAL(dp)                        :: SpOh_disp_vec(3), SpOh_disp_norm
REAL(dp)                        :: SpO_disp_vec(3), SpO_disp_norm

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, Udonnor_count, Udonnor_count2
INTEGER, ALLOCATABLE            :: nb_max_OHvec(:)

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o, l, n
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
ALLOCATE(atm_mat(24,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:))

nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN
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
ALLOCATE(OHvec_mat(38,nb_o*3,nb_step))
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

! C ----------------------------------------------- OH/O hbonds
start = OMP_get_wtime()
ALLOCATE(OHvec_hbond_mat(20,nb_o*3,nb_step))
OHvec_hbond_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_atm, nb_step, OHvec_hbond_mat)&
!$OMP SHARED(hb_oHpO_rcut)&
!$OMP PRIVATE(s, i, j, k, l, n)&
!$OMP PRIVATE(oHpO_disp_vec, oHpO_disp_norm, Udonnor_count, Udonnor_count2, OH_disp_norm, OH_disp_vec)&
!$OMP PRIVATE(alpha)
DO s = 1, nb_step
 C1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE C1
        END IF
        Udonnor_count = 0
        Udonnor_count2 = 0
        n = 4
        OHvec_hbond_mat(1,i,s) = OHvec_mat(1,i,s)
        OHvec_hbond_mat(2,i,s) = OHvec_mat(3,i,s)
        OHvec_hbond_mat(3,i,s) = OHvec_mat(2,i,s)
        DO j = 1, nb_atm
            IF (atm_mat(2,j,s) .EQ. 16) THEN
                DO k = 1, 3
                    oHpO_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+10,i,s)
                    oHpO_disp_vec(k) = oHpO_disp_vec(k) - box(k) * ANINT(oHpO_disp_vec(k)/box(k)) ! HO hbonds vector
                    OH_disp_vec(k) = OHvec_mat(k+10,i,s) - OHvec_mat(k+7,i,s)
                    OH_disp_vec(k) = OH_disp_vec(k) - box(k) * ANINT(OH_disp_vec(k)/box(k)) ! OH vector
                END DO

                oHpO_disp_norm = NORM2(oHpO_disp_vec) ! r
                OH_disp_norm = NORM2(OH_disp_vec)

                IF( (oHpO_disp_norm .LE. hb_oHpO_rcut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                    atm_mat(7,j,s) = atm_mat(7,j,s) + 1 ! Acceptor count (O)
                    OHvec_mat(15,i,s) = OHvec_mat(15,i,s) + 1 ! Donnor count (OH)
                    IF (n .GT. 23) THEN
                        PRINT*, "HBONDS OVERFLOW", OHvec_mat(1,i,s)
                    END IF

                    alpha = ACOS(DOT_PRODUCT(oHpO_disp_vec(:), OH_disp_vec(:)) / (oHpO_disp_norm * OH_disp_norm) ) ! alpha
                    OHvec_hbond_mat(n,i,s) = atm_mat(1,j,s)
                    OHvec_hbond_mat(n+1,i,s) = atm_mat(3,j,s)
                    OHvec_hbond_mat(n+2,i,s) = oHpO_disp_norm
                    OHvec_hbond_mat(n+3,i,s) = alpha
                    n = n + 4
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
PRINT'(A40,F14.2,A20)', "OH/O hbonds:", finish-start, "seconds elapsed"

! D ----------------------------------------------- Proximity between functionnal groups and any OH group
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_step, nb_atm)&
!$OMP SHARED(hb_X1Oh_rcut, hb_X2Oh_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XOh_disp_vec, XOh_disp_norm, XHo_disp_vec, XHo_disp_norm)
DO s = 1, nb_step
 D1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE D1
        END IF
     D2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE D2
            END IF
            IF (atm_mat(1,j,s) .EQ. OHvec_mat(1,i,s)) THEN
                CYCLE D2
            END IF
            DO k = 1, 3
                XOh_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+7,i,s)
                XOh_disp_vec(k) = XOh_disp_vec(k) - box(k) * ANINT(XOh_disp_vec(k)/box(k))
                XHo_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+10,i,s)
                XHo_disp_vec(k) = XHo_disp_vec(k) - box(k) * ANINT(XHo_disp_vec(k)/box(k))
            END DO
            XOh_disp_norm = NORM2(XOh_disp_vec)
            XHo_disp_norm = NORM2(XHo_disp_vec)

            IF ( ( (XOh_disp_norm .LE. OHvec_mat(32,i,s)) .OR.&
            (OHvec_mat(32,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                OHvec_mat(31,i,s) = atm_mat(1,j,s)
                OHvec_mat(32,i,s) = XOh_disp_norm
                OHvec_mat(33,i,s) = atm_mat(6,j,s)
            END IF

            IF ( ( (XHo_disp_norm .LE. OHvec_mat(36,i,s)) .OR.&
            (OHvec_mat(36,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                OHvec_mat(35,i,s) = atm_mat(1,j,s)
                OHvec_mat(36,i,s) = XHo_disp_norm
                OHvec_mat(37,i,s) = atm_mat(6,j,s)
            END IF

            IF (XOh_disp_norm .LE. hb_X1Oh_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(20,i,s) = 1
                    OHvec_mat(24,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    OHvec_mat(21,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    OHvec_mat(22,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    OHvec_mat(23,i,s) = 1
                END IF
            ELSE IF (XOh_disp_norm .LE. hb_X2Oh_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(24,i,s) = 1
                END IF
            END IF
        END DO D2
        IF (OHvec_mat(33,i,s) .GT. OHvec_mat(10,i,s) ) THEN
            OHvec_mat(34,i,s) = -1.0
        ELSE
            OHvec_mat(34,i,s) = 1.0
        END IF
        IF (OHvec_mat(37,i,s) .GT. OHvec_mat(13,i,s) ) THEN
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
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XO_disp_vec, XO_disp_norm)
DO s = 1, nb_step
    E1:DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .NE. 16) THEN
            CYCLE E1
        END IF
        E2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE E2
            END IF
            IF (atm_mat(1,j,s) .EQ. atm_mat(1,i,s)) THEN
                CYCLE E2
            END IF
            DO k = 1, 3
                XO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                XO_disp_vec(k) = XO_disp_vec(k) - box(k) * ANINT(XO_disp_vec(k)/box(k))
            END DO
            XO_disp_norm = NORM2(XO_disp_vec)

            IF ( ( (XO_disp_norm .LE. atm_mat(22,i,s)) .OR.&
            (atm_mat(22,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                atm_mat(21,i,s) = atm_mat(1,j,s)
                atm_mat(22,i,s) = XO_disp_norm
                atm_mat(23,i,s) = atm_mat(6,j,s)
            END IF
            IF (XO_disp_norm .LE. hb_X1O_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1) THEN ! C
                    atm_mat(10,i,s) = 1
                    atm_mat(14,i,s) = 1

                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    atm_mat(11,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    atm_mat(12,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    atm_mat(13,i,s) = 1
                END IF
            ELSE IF (XO_disp_norm .LE. hb_X2O_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1) THEN ! C
                    atm_mat(14,i,s) = 1
                END IF
            END IF
        END DO E2
        IF (atm_mat(23,i,s) .GT. atm_mat(6,i,s) ) THEN
            atm_mat(24,i,s) = -1.0
        ELSE
            atm_mat(24,i,s) = 1.0
        END IF
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
WRITE(31, '(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A20,A20,A20)')&
"Step", "Oid", "OType", "TA_O", "TD_O", "TDU_O"&
, "cC", "cOE", "cOH", "cOA", "cCX", "dist_ISD", "dist_ISU", "dist_AS"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            WRITE(31,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E20.7,E20.7,E20.7)')&
            s, INT(atm_mat(1,i,s)), INT(atm_mat(3,i,s)), INT(atm_mat(7,i,s))&
            , INT(atm_mat(8,i,s)), INT(atm_mat(9,i,s)), INT(atm_mat(10,i,s)), INT(atm_mat(11,i,s))&
            , INT(atm_mat(12,i,s)), INT(atm_mat(13,i,s)), INT(atm_mat(14,i,s))&
            , (atm_mat(15,i,s)*atm_mat(16,i,s)),(atm_mat(18,i,s)*atm_mat(19,i,s))&
            , (atm_mat(22,i,s)*atm_mat(24,i,s))
    END IF
    END DO
END DO
CLOSE(UNIT=31)

OPEN(UNIT=32, FILE = suffix//"_OH_hbonds.txt")
WRITE(32, '(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A20,A20,A20,A20)')&
    "Step","Oid", "Hid", "OType"&
    , "TA_OH", "TD_OH", "TDU_OH"&
    , "cC", "cOE", "cOH", "cOA", "cCX"&
    , "TA_O", "TD_O", "TDU_O", "dist_ISD", "dist_ISU"&
    , "dist_AS", "dist_HAS"
DO s = 1, nb_step
    DO i = 1, nb_max_OHvec(s)
        WRITE(32,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E20.7,E20.7,E20.7,E20.7)')&
        s, INT(OHvec_mat(1,i,s)), INT(OHvec_mat(2,i,s)), INT(OHvec_mat(3,i,s)) &
        , INT(OHvec_mat(14,i,s)), INT(OHvec_mat(15,i,s)), INT(OHvec_mat(16,i,s))&
        , INT(OHvec_mat(20,i,s)), INT(OHvec_mat(21,i,s)), INT(OHvec_mat(22,i,s)), INT(OHvec_mat(23,i,s)) , INT(OHvec_mat(24,i,s))&
        , INT(OHvec_mat(17,i,s)), INT(OHvec_mat(18,i,s)), INT(OHvec_mat(19,i,s))&
        , (OHvec_mat(25,i,s)*OHvec_mat(26,i,s)), (OHvec_mat(28,i,s)*OHvec_mat(29,i,s))&
        , (OHvec_mat(32,i,s)*OHvec_mat(34,i,s)), (OHvec_mat(36,i,s)*OHvec_mat(38,i,s))
    END DO
END DO
CLOSE(UNIT=32)

OPEN(UNIT=33, FILE = suffix//"_OH_hbonds_dist.txt")
WRITE(33, '(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10)')&
    "Step","OHid","Hid", "OType"&
    , "Oid", "OType", "r", "alpha"&
    , "Oid", "OType", "r", "alpha"&
    , "Oid", "OType", "r", "alpha"&
    , "Oid", "OType", "r", "alpha"

DO s = 1, nb_step
    DO i = 1, nb_max_OHvec(s)
        WRITE(33, '(I10,I10,I10,I10)', ADVANCE = "no")&
            s, INT(OHvec_hbond_mat(1,i,s)), INT(OHvec_hbond_mat(3,i,s)), INT(OHvec_hbond_mat(2,i,s))
        DO j = 4, 19, 4
            WRITE(33, '(I10,I10,F10.3,F10.3)', ADVANCE = "no")&
            INT(OHvec_hbond_mat(j,i,s)),INT(OHvec_hbond_mat(j+1,i,s)), OHvec_hbond_mat(j+2,i,s), OHvec_hbond_mat(j+3,i,s)
        END DO
        WRITE(33,'()')
    END DO
END DO
CLOSE(UNIT=33)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "O/OH Hbonds output:" ,finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat,nb_is)

DEALLOCATE(OHvec_mat, atm_mat,nb_max_OHvec)
END PROGRAM hbonds