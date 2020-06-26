!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM assign
USE OMP_LIB
USE INPUT

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:), atm_type(:,:), temp_atm_type(:,:)
REAL(dp), ALLOCATABLE           :: center_avgpos(:,:)

!   ----------------------------------------------- Temporary variables
REAL(dp)                        :: z_shift, temp_dist2(4)
REAL(dp)                        :: ij_disp_vec(3), ij_disp_norm

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h, nb_c, up, temp_id1, temp_id2(4)
INTEGER, ALLOCATABLE            :: nb_center(:)
INTEGER, ALLOCATABLE            :: nb_epoxide(:), nb_alcohol(:), nb_alkoxide(:)
INTEGER, ALLOCATABLE            :: nb_water(:), nb_hydronium(:), nb_hydroxide(:)
INTEGER, ALLOCATABLE            :: nb_oxygen_group(:), nb_ether(:)
REAL(dp), ALLOCATABLE           :: conn_mat(:,:,:)

!   ----------------------------------------------- Counters
INTEGER                         :: s, i, j, k, l, m
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Assign'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

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
file_vel=TRIM(file_vel)

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, Assign, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(15,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

!   ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_name(nb_atm))
ALLOCATE(atm_type(nb_atm,nb_step))

OPEN(UNIT=20, FILE=file_pos, STATUS='old', FORM='formatted', ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i = 1, nb_atm
        READ(20, *) atm_name(i), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        atm_mat(1,i,s) = i
        IF (atm_name(i)(1:1) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
        ELSE IF (atm_name(i)(1:1) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
        ELSE IF (atm_name(i)(1:1) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
        END IF
        IF (atm_name(i) .EQ. assign_center_name) THEN
            atm_mat(15,i,s) = 1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_c = COUNT(atm_mat(2,:,1) .EQ. 12, DIM=1)
nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

ALLOCATE(nb_center(nb_step))
DO s = 1, nb_step
    nb_center(s) = COUNT(atm_mat(15,:,s) .EQ. 1, DIM=1)
END DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Read velocities
IF (file_vel .NE. '0') THEN
    start = OMP_get_wtime()

    OPEN(UNIT=21, FILE=file_vel, STATUS='old', FORM='formatted', ACTION='READ')
    DO s = 1, nb_step
        READ(21, *)
        READ(21, *)
        DO i = 1, nb_atm
            READ(21, *) atm_name(i), atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
        END DO
    END DO
    CLOSE(UNIT=21)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Velocities:", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- Calculate geom center of go carbons and center so Zcarbon_avg = 0.0
IF ( WRAP_C .EQ. "Y") THEN
    start = OMP_get_wtime()

    ALLOCATE(center_avgpos(3,nb_step))

    DO s = 1, nb_step
        center_avgpos(:,s) = 0.0_dp

        DO i = 1, nb_atm
            IF (atm_mat(15,i,s) .EQ. 1) THEN
                DO k = 1, 3
                    center_avgpos(k,s) = center_avgpos(k,s) + atm_mat(k+2,i,s)
                END DO
            END IF
        END DO

        center_avgpos(:,s) = center_avgpos(:,s) / nb_center(s)
        z_shift = 0.0_dp - center_avgpos(3,s)

        DO i = 1, nb_atm
            atm_mat(5,i,s) = atm_mat(5,i,s) + z_shift
            DO k = 1, 3
                atm_mat(k+2,i,s) = atm_mat(k+2,i,s) - box(k) * ANINT(atm_mat(k+2,i,s) / box(k))
            END DO
        END DO

    END DO

    DEALLOCATE(center_avgpos,nb_center)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Center/Wrap:", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- Search the topology, water and oxygen groups
start = OMP_get_wtime()

ALLOCATE(conn_mat(17,nb_atm,nb_step))
ALLOCATE(nb_epoxide(nb_step))
ALLOCATE(nb_alcohol(nb_step))
ALLOCATE(nb_alkoxide(nb_step))
ALLOCATE(nb_water(nb_step))
ALLOCATE(nb_hydronium(nb_step))
ALLOCATE(nb_hydroxide(nb_step))
ALLOCATE(nb_ether(nb_step))
ALLOCATE(nb_oxygen_group(nb_step))
ALLOCATE(temp_atm_type(nb_atm,2))

!atm_mat(9,:,:) = 0 ! NbH
!atm_mat(10,:,:) = 0 ! NbC
!atm_mat(11,:,:) = 0 ! NbO
!atm_mat(12,:,:) = 0 ! NbC
!atm_mat(13,:,:) = 0 ! NbO
!atm_mat(14,:,:) = 0 ! NbH
conn_mat(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, nb_atm, nb_step, box,conn_mat, atm_type)&
!$OMP SHARED(assign_HO_rcut, assign_OC_rcut, assign_OO_rcut)&
!$OMP SHARED(assign_HC_rcut, assign_CC_rcut, assign_HH_rcut)&
!$OMP SHARED(nb_epoxide, nb_hydroxide, nb_alcohol, nb_water, nb_hydronium, nb_alkoxide, nb_oxygen_group)&
!$OMP SHARED(nb_ether)&
!$OMP PRIVATE(s, i, j, k, l, m)&
!$OMP PRIVATE(ij_disp_norm, ij_disp_vec, up, temp_id1, temp_id2, temp_dist2,temp_atm_type)
DO s = 1, nb_step
    DO i = 1, nb_atm-1
        DO j = i+1, nb_atm
            DO k = 1, 3
                ij_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                ij_disp_vec(k) = ij_disp_vec(k) - box(k) * ANINT(ij_disp_vec(k) / box(k))
            END DO
            ij_disp_norm = NORM2(ij_disp_vec)
            IF ( ij_disp_norm .LT. 2.00 ) THEN
                IF (atm_mat(2,i,s) .EQ. 16) THEN
                    IF (atm_mat(2,j,s) .EQ. 1) THEN
                        IF(ij_disp_norm .LE. assign_HO_rcut) THEN
                            atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                            atm_mat(11,j,s) = atm_mat(11,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 12) THEN
                        IF(ij_disp_norm .LE. assign_OC_rcut) THEN
                            atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                            atm_mat(11,j,s) = atm_mat(11,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 16) THEN
                        IF(ij_disp_norm .LE. assign_OO_rcut) THEN
                            atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                            atm_mat(11,j,s) = atm_mat(11,j,s) + 1
                            up = 1
                        END IF
                    END IF
                ELSE IF (atm_mat(2,i,s) .EQ. 12) THEN
                    IF (atm_mat(2,j,s) .EQ. 1) THEN
                        IF(ij_disp_norm .LE. assign_HC_rcut) THEN
                            atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                            atm_mat(10,j,s) = atm_mat(10,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 12) THEN
                        IF(ij_disp_norm .LE. assign_CC_rcut) THEN
                            atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                            atm_mat(10,j,s) = atm_mat(10,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 16) THEN
                        IF(ij_disp_norm .LE. assign_OC_rcut) THEN
                            atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                            atm_mat(10,j,s) = atm_mat(10,j,s) + 1
                            up = 1
                        END IF
                    END IF
                ELSE IF (atm_mat(2,i,s) .EQ. 1) THEN
                    IF (atm_mat(2,j,s) .EQ. 1) THEN
                        IF(ij_disp_norm .LE. assign_HH_rcut) THEN
                            atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                            atm_mat(9,j,s) = atm_mat(9,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 12) THEN
                        IF(ij_disp_norm .LE. assign_HC_rcut) THEN
                            atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                            atm_mat(9,j,s) = atm_mat(9,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF (atm_mat(2,j,s) .EQ. 16) THEN
                        IF(ij_disp_norm .LE. assign_HO_rcut) THEN
                            atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                            atm_mat(10,j,s) = atm_mat(10,j,s) + 1
                            up = 1
                        END IF
                    END IF
                END IF
            END IF
            IF (up .EQ. 1) THEN
                conn_mat(1,i,s) = atm_mat(1,i,s)
                conn_mat(2,i,s) = atm_mat(2,i,s)
                IF (conn_mat(3,i,s) .EQ. 0) THEN
                    conn_mat(3,i,s) = atm_mat(1,j,s)
                    conn_mat(4,i,s) = atm_mat(2,j,s)
                    conn_mat(5,i,s) = ij_disp_norm
                ELSE IF (conn_mat(6,i,s) .EQ. 0) THEN
                    conn_mat(6,i,s) = atm_mat(1,j,s)
                    conn_mat(7,i,s) = atm_mat(2,j,s)
                    conn_mat(8,i,s) = ij_disp_norm
                ELSE IF (conn_mat(9,i,s) .EQ. 0) THEN
                    conn_mat(9,i,s) = atm_mat(1,j,s)
                    conn_mat(10,i,s) = atm_mat(2,j,s)
                    conn_mat(11,i,s) = ij_disp_norm
                ELSE IF (conn_mat(12,i,s) .EQ. 0) THEN
                    conn_mat(12,i,s) = atm_mat(1,j,s)
                    conn_mat(13,i,s) = atm_mat(2,j,s)
                    conn_mat(14,i,s) = ij_disp_norm
                ELSE IF (conn_mat(15,i,s) .EQ. 0) THEN
                    conn_mat(15,i,s) = atm_mat(1,j,s)
                    conn_mat(16,i,s) = atm_mat(2,j,s)
                    conn_mat(17,i,s) = ij_disp_norm
                END IF
                conn_mat(1,j,s) = atm_mat(1,j,s)
                conn_mat(2,j,s) = atm_mat(2,j,s)
                IF (conn_mat(3,j,s) .EQ. 0) THEN
                    conn_mat(3,j,s) = atm_mat(1,i,s)
                    conn_mat(4,j,s) = atm_mat(2,i,s)
                    conn_mat(5,j,s) = ij_disp_norm
                ELSE IF (conn_mat(6,j,s) .EQ. 0) THEN
                    conn_mat(6,j,s) = atm_mat(1,i,s)
                    conn_mat(7,j,s) = atm_mat(2,i,s)
                    conn_mat(8,j,s) = ij_disp_norm
                ELSE IF (conn_mat(9,j,s) .EQ. 0) THEN
                    conn_mat(9,j,s) = atm_mat(1,i,s)
                    conn_mat(10,j,s) = atm_mat(2,i,s)
                    conn_mat(11,j,s) = ij_disp_norm
                ELSE IF (conn_mat(12,j,s) .EQ. 0) THEN
                    conn_mat(12,j,s) = atm_mat(1,i,s)
                    conn_mat(13,j,s) = atm_mat(2,i,s)
                    conn_mat(14,j,s) = ij_disp_norm
                ELSE IF (conn_mat(15,j,s) .EQ. 0) THEN
                    conn_mat(15,j,s) = atm_mat(1,i,s)
                    conn_mat(16,j,s) = atm_mat(2,i,s)
                    conn_mat(17,j,s) = ij_disp_norm
                END IF
            END IF
            up = 0
        END DO
    END DO

    DO i = 1, nb_atm
        DO j = 4, 18, 3
            IF (conn_mat(j,i,s) .EQ. 12) THEN
                atm_mat(12,i,s) = atm_mat(12,i,s) + 1
            ELSE IF (conn_mat(j,i,s) .EQ. 16) THEN
                atm_mat(13,i,s) = atm_mat(13,i,s) + 1
            ELSE IF (conn_mat(j,i,s) .EQ. 1) THEN
                atm_mat(14,i,s) = atm_mat(14,i,s) + 1
            END IF
        END DO
    END DO
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 2)) THEN
                atm_type(i,s) = "OW"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 3)) THEN
                atm_type(i,s) = "OP"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 1)) THEN
                atm_type(i,s) = "OM"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 1) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .GE. 1)) THEN
                atm_type(i,s) = "OH"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 1) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "OA"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 2) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .GE. 0)) THEN
                atm_type(i,s) = "OE"
            ELSE
                atm_type(i,s) = "OX"
                print*, s, i, atm_type(i,s), atm_mat(12,i,s), atm_mat(13,i,s), atm_mat(14,i,s)
            END IF
        ELSE IF (atm_mat(2,i,s) .EQ. 12) THEN
            IF ( (atm_mat(12,i,s) .EQ. 3) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "CC"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 3) .AND. (atm_mat(13,i,s) .EQ. 1) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "C3"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 2) .AND. (atm_mat(13,i,s) .EQ. 1) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "C2"
            ELSE IF ( (atm_mat(12,i,s) .EQ. 1) .AND. (atm_mat(13,i,s) .EQ. 1) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "C2" ! To Correct
            ELSE
                atm_type(i,s) = "CX"
                print*, s, i, atm_type(i,s), atm_mat(12,i,s), atm_mat(13,i,s), atm_mat(14,i,s)
            END IF
        ELSE IF (atm_mat(2,i,s) .EQ. 1) THEN
            IF ((atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 1) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "H1"
            ELSE IF ((atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 2) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "H2"
            ELSE IF ((atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 3) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                atm_type(i,s) = "H3"
            ELSE
                atm_type(i,s) = "HX"
                print*, s, i, atm_type(i,s), atm_mat(12,i,s), atm_mat(13,i,s), atm_mat(14,i,s)
            END IF
        END IF
    END DO
    temp_atm_type(:,:) = "XXX"
    temp_atm_type(:,1) = atm_type(:,s)
    DO l = 1, 100
        k = 0
        DO i = 1, nb_atm
            IF (temp_atm_type(i,1) .EQ. temp_atm_type(i,2)) THEN
                k = k + 1
            END IF
            temp_atm_type(i,2) = temp_atm_type(i,1)
        END DO
        IF (k .EQ. nb_atm) THEN
            EXIT
        END IF
        IF (l .EQ. 100) THEN
            PRINT'(A40,I3,A20)', "Topology not converged in step:", s
            EXIT
        END IF

     AS:DO i = 1, nb_atm
            ! ASSIGN CARBONS
            IF ( (atm_type(i,s) .EQ. "C3") .OR.&
                (atm_type(i,s) .EQ. "C3O") .OR.&
                (atm_type(i,s) .EQ. "C3A") .OR.&
                (atm_type(i,s) .EQ. "C3E"))THEN
                DO j = 4, 18, 3
                    IF (conn_mat(j,i,s) .EQ. 16) THEN
                        temp_id1 = INT(conn_mat(j-1,i,s))
                        IF ((atm_type(temp_id1,s) .EQ. "OH") .OR.&
                            (atm_type(temp_id1,s) .EQ. "OH3"))THEN
                            atm_type(temp_id1,s) = "OH3"
                            atm_type(i,s) = "C3O"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OA") .OR.&
                            (atm_type(temp_id1,s) .EQ. "OA3")) THEN
                            atm_type(temp_id1,s) = "OA3"
                            atm_type(i,s) = "C3A"
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OE") .AND. atm_mat(14,temp_id1,s) .LT. 1).OR.&
                            (atm_type(temp_id1,s) .EQ. "OEP" )) THEN
                            atm_type(temp_id1,s) = "OEP"
                            atm_type(i,s) = "C3E"
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OE") .AND. atm_mat(14,temp_id1,s) .GE. 1).OR.&
                            (atm_type(temp_id1,s) .EQ. "OEQ" )) THEN
                            atm_type(temp_id1,s) = "OEQ"
                            atm_type(i,s) = "C3E"
                        END IF
                    END IF
                END DO
            ! ASSIGN CARBONS
            ELSE IF ( (atm_type(i,s) .EQ. "C2") .OR.&
                (atm_type(i,s) .EQ. "C2O") .OR.&
                (atm_type(i,s) .EQ. "C2A") .OR.&
                (atm_type(i,s) .EQ. "C2E"))THEN
                DO j = 4, 18, 3
                    IF (conn_mat(j,i,s) .EQ. 16) THEN
                        temp_id1 = INT(conn_mat(j-1,i,s))
                        IF ((atm_type(temp_id1,s) .EQ. "OH") .OR.&
                            (atm_type(temp_id1,s) .EQ. "OH2")) THEN
                            atm_type(temp_id1,s) = "OH2"
                            atm_type(i,s) = "C2O"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OA") .OR.&
                            (atm_type(temp_id1,s) .EQ. "OA2")) THEN
                            atm_type(temp_id1,s) = "OA2"
                            atm_type(i,s) = "C2A"
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OE") .AND. atm_mat(14,temp_id1,s) .LT. 1).OR.&
                            (atm_type(temp_id1,s) .EQ. "OET" )) THEN
                            atm_type(temp_id1,s) = "OET"
                            atm_type(i,s) = "C2E"
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OE") .AND. atm_mat(14,temp_id1,s) .GE. 1).OR.&
                            (atm_type(temp_id1,s) .EQ. "OEU" )) THEN
                            atm_type(temp_id1,s) = "OEU"
                            atm_type(i,s) = "C2E"
                        END IF
                    END IF
                END DO
            ! ASSIGN MONOVALENT H
            ELSE IF ((atm_type(i,s) .EQ. "H1") .OR.&
                (atm_type(i,s) .EQ. "HE") .OR.&
                (atm_type(i,s) .EQ. "HU") .OR.&
                (atm_type(i,s) .EQ. "HO") .OR.&
                (atm_type(i,s) .EQ. "HW") .OR.&
                (atm_type(i,s) .EQ. "HP") .OR.&
                (atm_type(i,s) .EQ. "HM")) THEN
                DO j = 4, 18, 3
                    IF (conn_mat(j,i,s) .EQ. 16) THEN
                        temp_id1 = INT(conn_mat(j-1,i,s))
                        IF ((atm_type(temp_id1,s) .EQ. "OH") .OR.&
                        (atm_type(temp_id1,s) .EQ. "OH2") .OR.&
                        (atm_type(temp_id1,s) .EQ. "OH3")) THEN
                            atm_type(i,s) = "HO"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OM")) THEN
                            atm_type(i,s) = "HM"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OW")) THEN
                            atm_type(i,s) = "HW"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OP")) THEN
                            atm_type(i,s) = "HP"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OEQ")) THEN
                            atm_type(i,s) = "HE"
                        ELSE IF ((atm_type(temp_id1,s) .EQ. "OEU")) THEN
                            atm_type(i,s) = "HU"
                        END IF
                    END IF
                END DO
            ! CHECK DIVALENT H
            ELSE IF ((atm_type(i,s) .EQ. "H2")) THEN
                k = 1
                temp_id2(:) = 0
                temp_dist2(:) = 0.0_dp
                DO j = 4, 18, 3
                    IF (conn_mat(j,i,s) .EQ. 16) THEN
                        temp_id2(k) = INT(conn_mat(j-1,i,s))
                        temp_dist2(k) = conn_mat(j+1,i,s)
                        k = k + 1
                    END IF
                END DO
                IF ( (temp_dist2(1) .EQ. 0) .OR.&
                (temp_dist2(2) .EQ. 0)) THEN
                    PRINT*,"WIERD",s,i
                    CYCLE AS
                END IF
                atm_type(i,s) = "H1"
                atm_mat(14,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = atm_mat(14,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) - 1
                DO j = 4, 18, 3
                    IF ( temp_id2(MAXLOC(temp_dist2, DIM=1)) .EQ. conn_mat(j-1,i,s) ) THEN
                        conn_mat(j,i,s) = 0
                        conn_mat(j+1,i,s) = 0
                        conn_mat(j-1,i,s) = 0
                    END IF
                    IF ( atm_mat(1,i,s) .EQ. conn_mat(j-1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) ) THEN
                        conn_mat(j,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                        conn_mat(j+1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                        conn_mat(j-1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                    END IF
                END DO
                m=1
            ! CHECK TRIVALENT H
            ELSE IF ((atm_type(i,s) .EQ. "H3")) THEN
                k = 1
                temp_id2(:) = 0
                temp_dist2(:) = 0.0_dp
                DO j = 4, 18, 3
                    IF (conn_mat(j,i,s) .EQ. 16) THEN
                        temp_id2(k) = INT(conn_mat(j-1,i,s))
                        temp_dist2(k) = conn_mat(j+1,i,s)
                        k = k + 1
                    END IF
                END DO
                IF ( (temp_dist2(1) .EQ. 0) .OR.&
                    (temp_dist2(2) .EQ. 0) .OR.&
                    (temp_dist2(3) .EQ. 0)) THEN
                    PRINT*,"WIERD",s,i
                    CYCLE AS
                END IF
                atm_type(i,s) = "H2"
                atm_mat(14,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = atm_mat(14,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) - 1
                DO j = 4, 18, 3
                    IF ( temp_id2(MAXLOC(temp_dist2, DIM=1)) .EQ. conn_mat(j-1,i,s) ) THEN
                        conn_mat(j,i,s) = 0
                        conn_mat(j+1,i,s) = 0
                        conn_mat(j-1,i,s) = 0
                    END IF
                    IF ( atm_mat(1,i,s) .EQ. conn_mat(j-1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) ) THEN
                        conn_mat(j,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                        conn_mat(j+1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                        conn_mat(j-1,temp_id2(MAXLOC(temp_dist2, DIM=1)),s) = 0
                    END IF
                END DO
                m=1
            END IF
            temp_atm_type(i,1) = atm_type(i,s)
        END DO AS
        ! Reupdate Oxygen
        IF (m .EQ. 1) THEN
            DO i = 1, nb_atm
                IF ((atm_mat(2,i,s) .EQ. 16) .AND. (m .EQ. 1) )THEN
                    IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 2)) THEN
                        atm_type(i,s) = "OW"
                    ELSE IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 3)) THEN
                        atm_type(i,s) = "OP"
                    ELSE IF ( (atm_mat(12,i,s) .EQ. 0) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 1)) THEN
                        atm_type(i,s) = "OM"
                    ELSE IF ( (atm_mat(12,i,s) .EQ. 1) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .GE. 1)) THEN
                        atm_type(i,s) = "OH"
                    ELSE IF ( (atm_mat(12,i,s) .EQ. 1) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .EQ. 0)) THEN
                        atm_type(i,s) = "OA"
                    ELSE IF ( (atm_mat(12,i,s) .EQ. 2) .AND. (atm_mat(13,i,s) .EQ. 0) .AND. (atm_mat(14,i,s) .GE. 0)) THEN
                        atm_type(i,s) = "OE"
                    ELSE
                        atm_type(i,s) = "OX"
                        PRINT*, s, i, atm_type(i,s), atm_mat(12,i,s), atm_mat(13,i,s), atm_mat(14,i,s)
                    END IF
                END IF
                temp_atm_type(i,1) = atm_type(i,s)
            END DO
        END IF
        m = 0

    END DO
    nb_epoxide(s) = COUNT( (atm_type(:,s) .EQ. "OEP"), DIM=1) +  COUNT( (atm_type(:,s) .EQ. "OEQ"), DIM=1)
    nb_ether(s) = COUNT( (atm_type(:,s) .EQ. "OET"), DIM=1) +  COUNT( (atm_type(:,s) .EQ. "OEU"), DIM=1)
    nb_alcohol(s) = COUNT( (atm_type(:,s) .EQ. "OH2"), DIM=1) + COUNT( (atm_type(:,s) .EQ. "OH3"), DIM=1)
    nb_alkoxide(s) = COUNT( (atm_type(:,s) .EQ. "OA2"), DIM=1) + COUNT( (atm_type(:,s) .EQ. "OA3"), DIM=1)
    nb_water(s) = COUNT( (atm_type(:,s) .EQ. "OW"), DIM=1)
    nb_hydroxide(s) = COUNT( (atm_type(:,s) .EQ. "OM"), DIM=1)
    nb_hydronium(s) = COUNT( (atm_type(:,s) .EQ. "OP"), DIM=1)
    nb_oxygen_group(s) = nb_epoxide(s) + nb_alcohol(s) + nb_ether(s) + &
        nb_alkoxide(s) + nb_hydroxide(s) + nb_hydronium(s) + nb_water(s)

END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Topology:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print counts
start = OMP_get_wtime()

OPEN(UNIT=50, FILE = suffix//"_oxygen_groups_population.txt")
WRITE(50, '(A4,1X,A10,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4)') &
    "Traj", "Step", "EPO", "ETH", "ALC", "ALK", "H2O", "OH-", "H3O+", "TSum", "TO"
DO s = 1, nb_step
    WRITE(50, '(A4,1X,I10,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4)') &
    suffix, s, nb_epoxide(s), nb_ether(s), nb_alcohol(s), nb_alkoxide(s)&
    , nb_water(s), nb_hydroxide(s), nb_hydronium(s), nb_oxygen_group(s), nb_o
END DO
CLOSE(UNIT=50)

DEALLOCATE(nb_epoxide,nb_alcohol,nb_ether,nb_alkoxide,nb_water,nb_hydronium,nb_hydroxide,nb_oxygen_group)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Oxygen groups topology output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the xyz and velocities files
start = OMP_get_wtime()
IF ( WRAP_C .EQ. "Y") THEN
    OPEN(UNIT=40, FILE = suffix//"_wrapped_"//file_pos)
    IF (file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_wrapped_"//file_vel)
ELSE 
    OPEN(UNIT=40, FILE = suffix//"_nonwrapped_"//file_pos)
    IF (file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_nonwrapped_"//file_vel)
END IF
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_atm
    IF (file_vel .NE. '0') WRITE(41,'(I10)') nb_atm
    WRITE(40,'(A10,I10)') "Step nb:", s
    IF (file_vel .NE. '0') WRITE(41,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_atm
        WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL(atm_type(i,s)), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        IF (file_vel .NE. '0') WRITE(41,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL(atm_type(i,s))&
        , atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
    END DO
END DO
CLOSE(UNIT=40)
IF (file_vel .NE. '0') CLOSE(UNIT=41)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions/Velocities output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the waterlist (mask indices for surf)
IF (waterlist .EQ. 1) THEN
    start = OMP_get_wtime()

    OPEN(UNIT=42, FILE = suffix//"_waterlist.txt")
    DO s = 1, nb_step
        WRITE(42,'(A14)', ADVANCE="no")"mask = indices "
        DO i = 1, nb_atm
            IF (atm_type(i,s) .EQ. "OW") THEN
                WRITE(42,'(I5)', ADVANCE="no") INT(atm_mat(1,i,s) - 1)
            END IF
        END DO
        WRITE(42,*)
    END DO
    CLOSE(UNIT=42)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Waterlist (surf) output:", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_name,atm_mat,conn_mat)

END PROGRAM assign