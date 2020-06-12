!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM vvcf
USE OMP_LIB
USE INPUT

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start, finish, start_i, finish_i, avg_timigs, timings_t
REAL(dp)                        :: start_tot, finish_tot
REAL(dp), ALLOCATABLE           :: timings(:)

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), OHvec_mat(:,:,:), is_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:)
INTEGER                         :: nb_line_is, nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: OpH_disp_vec(3), OpH_disp_norm
REAL(dp)                        :: oHpO_disp_vec(3), oHpO_disp_norm
REAL(dp)                        :: XOh_disp_vec(3), XOh_disp_norm
REAL(dp)                        :: SpOh_disp_vec(3), SpOh_disp_norm
CHARACTER(LEN=2)                :: dummy_char

! ----------------------------------------------- VVCF
INTEGER                         :: mcs, mcsb
REAL(dp)                        :: OH_vel_vec(3), OH_disp_vec(3)
REAL(dp), ALLOCATABLE           :: vvcf_xxz(:)
REAL(dp)                        :: tij_vec(3), trij, vvcf_xxz_t
INTEGER, ALLOCATABLE            :: v(:,:)

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, Udonnor_count, Udonnor_count2

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o, l, t
INTEGER                         :: CAC, keep

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching VVCF'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
start_tot = OMP_get_wtime()

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
file_vel=TRIM(file_vel)
file_is=TRIM(file_is)

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, VVCF, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(12,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_name(nb_atm))

OPEN(UNIT=20, FILE=file_pos, STATUS='old', FORM='formatted', ACTION='READ')
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

nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)

DEALLOCATE(atm_name)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Read velocities
start = OMP_get_wtime()

OPEN(UNIT=21, FILE=file_vel, STATUS='old', FORM='formatted', ACTION='READ')
DO s = 1, nb_step
    READ(21, *)
    READ(21, *)
    DO i = 1, nb_atm
        READ(21, *) dummy_char, atm_mat(7,i,s), atm_mat(8,i,s), atm_mat(9,i,s)
    END DO
END DO
CLOSE(UNIT=21)

! Velocities are in bohr/aut
DO k = 1, 3
    atm_mat(k+6,:,:) = atm_mat(k+6,:,:) * bohr_to_angstrom / aut_to_fs
END DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Velocities:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN
    ! A ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
    start = OMP_get_wtime()

    OPEN(UNIT=22, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
    nb_line_is=0
    DO
        READ(22, *, IOSTAT=iostatus)
        IF (iostatus .NE. 0) THEN
            EXIT
        ELSE
            nb_line_is = nb_line_is + 1
        END IF
    END DO
    REWIND(22)
    nb_max_is = CEILING(1.0 * nb_line_is / nb_step) * 2

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

    ! A ----------------------------------------------- Read IS
    start = OMP_get_wtime()

    ALLOCATE(is_mat(5,nb_max_is,nb_step))
    ALLOCATE(nb_is(nb_step))
    is_mat(:,:,:) = 0.0_dp
    nb_is(:) = 0

    OPEN(UNIT=22, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
    DO s = 1, nb_step
        READ(22, *) nb_is(s)
        READ(22, *)
        j = 0
        DO i=1,nb_is(s)
            READ(22, *) dummy_char, is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
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
    CLOSE(UNIT=22)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
END IF

! B ----------------------------------------------- OH groups and corresponding values
start = OMP_get_wtime()
ALLOCATE(OHvec_mat(39,nb_o*3,nb_step))

OHvec_mat(:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_step, nb_atm)&
!$OMP SHARED(vvcf_OpH_rcut)&
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
                    IF(OpH_disp_norm .LT. vvcf_OpH_rcut) THEN
                        o = o + 1
                        OHvec_mat(1,o,s) = atm_mat(1,i,s)
                        OHvec_mat(2,o,s) = atm_mat(1,j,s)
                        OHvec_mat(3,o,s) = atm_mat(3,i,s)
                        OHvec_mat(4,o,s) = atm_mat(3,j,s)
                        DO k = 1, 3
                            OHvec_mat(k+4,o,s) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                            OHvec_mat(k+4,o,s) = OHvec_mat(k+4,o,s) - box(k) * ANINT(OHvec_mat(k+4,o,s)/box(k)) ! Disp
                            OHvec_mat(k+7,o,s) = atm_mat(k+6,j,s) - atm_mat(k+6,i,s)
                            OHvec_mat(k+10,o,s) = atm_mat(k+3,i,s) ! O pos
                            OHvec_mat(k+13,o,s) = atm_mat(k+3,j,s) ! H pos
                            OHvec_mat(k+36,o,s) = OHvec_mat(k+10,o,s) + (1.00784/(15.999+1.00784)) * OHvec_mat(k+4,o,s)
                        END DO
                    END IF
                END IF
            END DO
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH groups:", finish-start, "seconds elapsed"

! C ----------------------------------------------- OH/O Hbonds
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, OHvec_mat, nb_o, nb_atm, nb_step)&
!$OMP SHARED(vvcf_oHpO_rcut)&
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
                    oHpO_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+13,i,s)
                    oHpO_disp_vec(k) = oHpO_disp_vec(k) - box(k) * ANINT(oHpO_disp_vec(k)/box(k))
                END DO
                oHpO_disp_norm = NORM2(oHpO_disp_vec)
                IF( (oHpO_disp_norm .LE. vvcf_oHpO_rcut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                    atm_mat(10,j,s) = atm_mat(10,j,s) + 1 ! Acceptor count (O)
                    OHvec_mat(18,i,s) = OHvec_mat(18,i,s) + 1 ! Donnor count (OH)
                    IF (Udonnor_count .EQ. 0) THEN
                        OHvec_mat(19,i,s) = OHvec_mat(19,i,s) + 1 ! Unique donnor count (OH)
                        Udonnor_count = 1
                    END IF
                C21:DO l=1, nb_atm
                        IF (atm_mat(1,l,s) .EQ. OHvec_mat(1,i,s)) THEN
                            atm_mat(11,l,s) = atm_mat(11,l,s) + 1 ! Donnor count (O)
                            IF (Udonnor_count2 .EQ. 0) THEN
                                atm_mat(12,l,s) = atm_mat(12,l,s) + 1 ! Unique donnor count (O)
                                Udonnor_count2 = 1
                            END IF
                            EXIT C21
                        END IF
                    END DO C21
                C22:DO l=1,nb_o*3
                        IF (OHvec_mat(1,l,s) .EQ. 0) THEN
                            CYCLE C22
                        ELSE IF (OHvec_mat(1,l,s) .EQ. atm_mat(1,j,s)) THEN
                            OHvec_mat(17,l,s) = OHvec_mat(17,l,s) + 1 ! Acceptor count (OH)
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
                OHvec_mat(20,i,s) = atm_mat(10,j,s)
                OHvec_mat(21,i,s) = atm_mat(11,j,s)
                OHvec_mat(22,i,s) = atm_mat(12,j,s)
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
!$OMP SHARED(vvcf_X1Oh_rcut, vvcf_X2Oh_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XOh_disp_vec, XOh_disp_norm)
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
                XOh_disp_vec(k) = atm_mat(k+3,j,s) - OHvec_mat(k+10,i,s)
                XOh_disp_vec(k) = XOh_disp_vec(k) - box(k) * ANINT(XOh_disp_vec(k)/box(k))
            END DO
            XOh_disp_norm = NORM2(XOh_disp_vec)

            IF ( ( (XOh_disp_norm .LE. OHvec_mat(35,i,s)) .OR.&
            (OHvec_mat(35,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                OHvec_mat(34,i,s) = atm_mat(1,j,s)
                OHvec_mat(35,i,s) = XOh_disp_norm
                OHvec_mat(36,i,s) = atm_mat(6,j,s)
            END IF
            IF (XOh_disp_norm .LE. vvcf_X1Oh_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(23,i,s) = 1
                    OHvec_mat(27,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    OHvec_mat(24,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    OHvec_mat(25,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    OHvec_mat(26,i,s) = 1
                END IF
            ELSE IF (XOh_disp_norm .LE. vvcf_X2Oh_rcut) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(27,i,s) = 1
                END IF
            END IF
        END DO D2
    END DO D1
END DO
!$OMP END PARALLEL DO

DEALLOCATE(atm_mat)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and OH group:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN

    ! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(box, OHvec_mat, nb_o, is_mat, nb_is, nb_step)&
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
                        SpOh_disp_vec(k) = is_mat(k,j,s) - OHvec_mat(k+10,i,s)
                        SpOh_disp_vec(k) = SpOh_disp_vec(k) - box(k) * ANINT(SpOh_disp_vec(k)/box(k))
                    END DO
                    SpOh_disp_norm = NORM2(SpOh_disp_vec)
                    IF ( (SpOh_disp_norm .LT. OHvec_mat(28,i,s)) .OR. (OHvec_mat(28,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(28,i,s) = SpOh_disp_norm
                        IF (OHvec_mat(13,i,s) .LT. is_mat(3,j,s)) THEN
                            OHvec_mat(29,i,s) = -1
                        ELSE
                            OHvec_mat(29,i,s) = 1
                        END IF
                        OHvec_mat(30,i,s) = is_mat(5,j,s)

                    END IF
                ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        SpOh_disp_vec(k) = is_mat(k,j,s) - OHvec_mat(k+10,i,s)
                        SpOh_disp_vec(k) = SpOh_disp_vec(k) - box(k) * ANINT(SpOh_disp_vec(k)/box(k))
                    END DO
                    SpOh_disp_norm = NORM2(SpOh_disp_vec)
                    IF ( (SpOh_disp_norm .LT. OHvec_mat(31,i,s)) .OR. (OHvec_mat(31,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(31,i,s) = SpOh_disp_norm
                        IF (OHvec_mat(13,i,s) .GT. is_mat(3,j,s)) THEN
                            OHvec_mat(32,i,s) = -1
                        ELSE
                            OHvec_mat(32,i,s) = 1
                        END IF
                        OHvec_mat(33,i,s) = is_mat(5,j,s)
                    END IF
                END IF
            END DO
        END DO F2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity IS and OH groups:", finish-start, "seconds elapsed"

END IF

IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat, nb_is)

! H ----------------------------------------------- VVAF/VVCF
start = OMP_get_wtime()
mcs = INT(mct / timestep_fs)
mcsb = INT(mctb / timestep_fs) + 1

IF (mcs+1 .GE. nb_step) THEN
    PRINT*, "Error: Max correlation time > trajectory time"
    STOP
END IF
ALLOCATE(vvcf_xxz(mcs+1))
ALLOCATE(timings(mcs+1))
timings(:) = 0.0_dp
vvcf_xxz(:) = 0.0_dp

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(NONE) SHARED(vvcf_xxz, mcsb, mcs, OHvec_mat, nb_o, box, timings, nb_step)&
!$OMP SHARED(layers_s, hbonds_s, intracorr, autocorr, IS_c, layer_up, layer_down, timestep_fs, vvcf_rcut)&
!$OMP SHARED(water_only, close_c_only, up_down_only, close_gl_only, close_gol_only, close_ol_only, hbonds_double)&
!$OMP SHARED(dcle, dcgt, acle, acgt, dcle2, dcgt2, acle2, acgt2, vvcf_X2Oh_rcut)&
!$OMP PRIVATE(t, s, i, j, l, v)&
!$OMP PRIVATE(tij_vec, trij, OH_vel_vec, OH_disp_vec, start_i, finish_i, keep, vvcf_xxz_t, timings_t)
DO t = mcsb, mcs+1
    start_i = OMP_get_wtime()
    ALLOCATE(v(nb_o*3,nb_step))
    v(:,:) = 0.0_dp
    vvcf_xxz_t = 0.0_dp

    DO s = 1, nb_step
        IF (s-1+t-1 .LT. nb_step) THEN
         H1:DO i = 1, nb_o*3

                IF (OHvec_mat(1,i,s) .EQ. 0) THEN ! Skip empty
                    CYCLE H1
                END IF

                IF ( (water_only .EQ. "Y") .AND.& ! Select water only
                (OHvec_mat(3,i,s) .NE. 13) ) THEN
                    CYCLE H1
                ELSE IF ( (water_only .EQ. "E") .AND.& ! Exclude water
                (OHvec_mat(3,i,s) .EQ. 13) ) THEN
                    CYCLE H1
                END IF

                IF ( (close_c_only .EQ. "Y") .AND.& ! Select OHvec close to the C surface (<=value)
                (OHvec_mat(31,i,s) .GT. vvcf_X2Oh_rcut) ) THEN
                    CYCLE H1
                END IF

                IF ( (close_gl_only .EQ. "Y") .AND. (&
                (OHvec_mat(23,i,s) .EQ. 0) .OR.&
                ( (OHvec_mat(24,i,s) .EQ. 1) .OR.&
                (OHvec_mat(25,i,s) .EQ. 1) .OR.&
                (OHvec_mat(26,i,s) .EQ. 1) ) ) )THEN
                    CYCLE H1
                END IF

                IF ( (close_gol_only .EQ. "Y") .AND. (&
                (OHvec_mat(23,i,s) .EQ. 0) .OR.&
                ( (OHvec_mat(24,i,s) .EQ. 0) .OR.&
                (OHvec_mat(25,i,s) .EQ. 0) .OR.&
                (OHvec_mat(26,i,s) .EQ. 0) ) ) )THEN
                    CYCLE H1
                END IF

                IF ( (close_ol_only .EQ. "Y") .AND. (&
                (OHvec_mat(23,i,s) .EQ. 1) .OR.&
                ( (OHvec_mat(24,i,s) .EQ. 0) .OR.&
                (OHvec_mat(25,i,s) .EQ. 0) .OR.&
                (OHvec_mat(26,i,s) .EQ. 0) ) ) )THEN
                    CYCLE H1
                END IF

                IF ( (up_down_only .EQ. "D") .AND.&
                (OHvec_mat(13,i,s) .GT. OHvec_mat(36,i,s)) ) THEN ! Down only
                    CYCLE H1
                ELSE IF ( (up_down_only .EQ. "U") .AND.&
                (OHvec_mat(13,i,s) .LT. OHvec_mat(36,i,s)) ) THEN ! Up only
                    CYCLE H1
                END IF

                IF ((layers_s .EQ. "Y") .AND. (IS_c .EQ. 'Y' )) THEN

                    IF ( (OHvec_mat(31,i,s)*OHvec_mat(32,i,s) .GT. layer_up) .OR.&
                    (OHvec_mat(31,i,s)*OHvec_mat(32,i,s) .LE. layer_down)) THEN ! Surface distance (go)
                        CYCLE H1
                    END IF

                END IF

                IF (hbonds_s .EQ. "Y") THEN

                    IF (hbonds_double .EQ. "N") THEN
                        IF ( ( (OHvec_mat(21,i,s) .GT. DCGT) .OR. (OHvec_mat(21,i,s) .LE. DCLE) ) .OR.&
                        ( (OHvec_mat(20,i,s) .GT. ACGT ) .OR. (OHvec_mat(20,i,s) .LE. ACLE) ) ) THEN
                            CYCLE H1
                        END IF

                    ELSE IF (hbonds_double .EQ. "Y") THEN
                        keep=0
                        IF ( ( (OHvec_mat(21,i,s) .LE. DCGT) .AND. (OHvec_mat(21,i,s) .GT. DCLE) ) .AND.&
                        ( (OHvec_mat(20,i,s) .LE. ACGT ) .AND. (OHvec_mat(20,i,s) .GT. ACLE) ) ) THEN
                            keep=1
                        ELSE IF ( ( (OHvec_mat(21,i,s) .LE. DCGT2) .AND. (OHvec_mat(21,i,s) .GT. DCLE2) ) .AND.&
                        ( (OHvec_mat(20,i,s) .LE. ACGT2 ) .AND. (OHvec_mat(20,i,s) .GT. ACLE2) ) ) THEN
                            keep=1
                        ELSE
                            keep=0
                        END IF

                        IF (keep .EQ. 0) CYCLE H1

                    END IF

                END IF

             H2:DO j = 1, nb_o * 3

                    IF ( (OHvec_mat(1,j,s) .EQ. 0) ) THEN ! Skip empty
                        CYCLE H2
                    END IF

                    ! Trick for intra cross-correlation only
                    IF (intracorr .EQ. "Y") THEN
                        IF (OHvec_mat(1,i,s) .NE. OHvec_mat(1,j,s)) THEN
                            CYCLE H2
                        END IF
                    END IF

                    ! Trick for autocorrelation
                    IF (autocorr .EQ. "Y") THEN
                        IF (OHvec_mat(1,i,s) .NE. OHvec_mat(1,j,s)) THEN
                            CYCLE H2
                        ELSE IF (OHvec_mat(2,i,s) .NE. OHvec_mat(2,j,s)) THEN
                            CYCLE H2
                        END IF
                    END IF

                    DO k = 1, 3
                        tij_vec(k) = OHvec_mat(k+36,j,s) - OHvec_mat(k+36,i,s)
                        tij_vec(k) = tij_vec(k) - box(k) * ANINT(tij_vec(k)/box(k))
                    END DO
                    trij = NORM2(tij_vec(:))

                    IF ( trij .LT. vvcf_rcut ) THEN

                     H3:DO l = 1, nb_o *3

                            IF ( (OHvec_mat(1,j,s) .NE. OHvec_mat(1,l,s+t-1)) .OR.&
                             (OHvec_mat(2,j,s) .NE. OHvec_mat(2,l,s+t-1)) ) THEN
                                CYCLE H3
                            END IF

                            DO k = 1, 3
                                OH_vel_vec(k) =  OHvec_mat(k+7,l,s+t-1)
                                OH_disp_vec(k) = OHvec_mat(k+4,l,s+t-1)
                            END DO

                            vvcf_xxz_t = vvcf_xxz_t + (OHvec_mat(10,i,s) &
                            * &
                            (DOT_PRODUCT(OH_vel_vec(:),OH_disp_vec(:)) &
                            / &
                            NORM2(OH_disp_vec(:))))

                        END DO H3 ! j at t

                    END IF

                END DO H2 ! j at 0

            END DO H1 ! i at 0

        END IF

    END DO ! Time

    DEALLOCATE(v)

    vvcf_xxz_t = vvcf_xxz_t * 1.0 / (nb_step - (t-1))

    finish_i = OMP_get_wtime()
    timings_t = finish_i-start_i

    IF (MODULO(t,25) .EQ. 0) THEN
        PRINT('(I10,A1,I10,E24.14,E24.14,E24.14)'), t, "/", mcs+1, (t-1)*timestep_fs, vvcf_xxz_t, timings_t
    END IF
    vvcf_xxz(t) = vvcf_xxz_t
    timings(t) = timings_t

END DO ! Corr
!$OMP END PARALLEL DO

avg_timigs = (SUM(timings(:)) / (mcs+1-mcsb) )

DEALLOCATE(timings)

finish = OMP_get_wtime()

IF (autocorr .EQ. "Y") THEN
    PRINT'(A40,F14.2,A20,A20,F14.2)', "Done with VVAF_xxz:", finish-start, "seconds elapsed", "avg per step:", avg_timigs
ELSE
    PRINT'(A40,F14.2,A20,A20,F14.2)', "Done with VVCF_xxz:", finish-start, "seconds elapsed", "avg per step:", avg_timigs
END IF

start = OMP_get_wtime()

IF (autocorr .EQ. "Y") THEN
    OPEN(UNIT=51, FILE=suffix//"_vvaf_xxz.txt")
    WRITE(51,'(A20,A20)') "Time (fs)", "VVAF_xxz (Å2.fs2)"
ELSE
    OPEN(UNIT=51, FILE=suffix//"_vvcf_xxz.txt")
    WRITE(51,'(A20,A20)') "Time (fs)", "VVCF_xxz (Å2.fs2)"
END IF
    DO t = mcsb, mcs+1
        WRITE(51,'(E24.14,E24.14)') (t-1)*timestep_fs, vvcf_xxz(t)
    END DO
CLOSE(UNIT=51)

finish = OMP_get_wtime()

IF (autocorr .EQ. "Y") THEN
    PRINT'(A40,F14.2,A20)', "Done with VVAF_xxz output:", finish-start, "seconds elapsed"
ELSE
    PRINT'(A40,F14.2,A20)', "Done with VVCF_xxz output:", finish-start, "seconds elapsed"
END IF

DEALLOCATE(vvcf_xxz)

finish_tot = OMP_get_wtime()
!   ----------------------------------------------- End

PRINT'(A40,F14.2,A20)', "Total time:", finish_tot-start_tot, "seconds elapsed"
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit

END PROGRAM vvcf