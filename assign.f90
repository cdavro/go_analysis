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
REAL(dp)                        :: z_shift
REAL(dp)                        :: ij_disp_vec(3), ij_disp_norm
REAL(dp)                        :: temp_dist
INTEGER                         :: temp_el, temp_id

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h, nb_c, up, temp_id1
INTEGER, ALLOCATABLE            :: nb_center(:)
INTEGER, ALLOCATABLE            :: nb_epoxide(:), nb_pepoxide(:)
INTEGER, ALLOCATABLE            :: nb_ether(:),  nb_pether(:), nb_palcohols(:)
INTEGER, ALLOCATABLE            :: nb_alcohol(:), nb_alkoxide(:), nb_ketone(:)
INTEGER, ALLOCATABLE            :: nb_water(:), nb_hydronium(:), nb_hydroxide(:)
INTEGER, ALLOCATABLE            :: nb_oxygen_group(:)
REAL(dp), ALLOCATABLE           :: conn_mat(:,:,:)

!   ----------------------------------------------- Counters
INTEGER                         :: s, i, j, k, l, n
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Assign'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF ( CAC .EQ. 0 ) THEN
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
ALLOCATE(atm_mat(17,nb_atm,nb_step))
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
        IF ( atm_name(i)(1:2) .EQ. "Na" ) THEN
            atm_mat(2,i,s) = 11
            atm_type(i,s) = "NA"
        ELSE IF ( atm_name(i)(1:2) .EQ. "Cl" ) THEN
            atm_mat(2,i,s) = 17
            atm_type(i,s) = "CLM"
        ELSE IF ( atm_name(i)(1:1) .EQ. "C" ) THEN
            atm_mat(2,i,s) = 12
        ELSE IF ( atm_name(i)(1:1) .EQ. "O" ) THEN
            atm_mat(2,i,s) = 16
        ELSE IF ( atm_name(i)(1:1) .EQ. "H" ) THEN
            atm_mat(2,i,s) = 1
        END IF
        IF ( atm_name(i) .EQ. assign_center_name ) THEN
            atm_mat(17,i,s) = 1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_c = COUNT( atm_mat(2,:,1) .EQ. 12, DIM=1 )
nb_o = COUNT( atm_mat(2,:,1) .EQ. 16, DIM=1 )
nb_h = COUNT( atm_mat(2,:,1) .EQ. 1, DIM=1 )

ALLOCATE(nb_center(nb_step))
DO s = 1, nb_step
    nb_center(s) = COUNT( atm_mat(17,:,s) .EQ. 1, DIM=1 )
END DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Read velocities
IF ( file_vel .NE. '0' ) THEN
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
IF ( WRAP_C .EQ. "Y" ) THEN
    start = OMP_get_wtime()

    ALLOCATE(center_avgpos(3,nb_step))

    DO s = 1, nb_step
        center_avgpos(:,s) = 0.0_dp

        DO i = 1, nb_atm
            IF ( atm_mat(17,i,s) .EQ. 1 ) THEN
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
                atm_mat(k+2,i,s) = atm_mat(k+2,i,s) - box(k) * ANINT( atm_mat(k+2,i,s) / box(k) )
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
ALLOCATE(nb_pepoxide(nb_step))
ALLOCATE(nb_ether(nb_step))
ALLOCATE(nb_pether(nb_step))
ALLOCATE(nb_alcohol(nb_step))
ALLOCATE(nb_alkoxide(nb_step))
ALLOCATE(nb_ketone(nb_step))
ALLOCATE(nb_water(nb_step))
ALLOCATE(nb_hydronium(nb_step))
ALLOCATE(nb_hydroxide(nb_step))
ALLOCATE(nb_palcohols(nb_step))
ALLOCATE(nb_oxygen_group(nb_step))
ALLOCATE(temp_atm_type(nb_atm,2))

conn_mat(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, nb_atm, nb_step, box, conn_mat, atm_type)&
!$OMP SHARED(assign_OC_rcut, assign_OO_rcut)&
!$OMP SHARED(assign_CC_rcut, assign_HH_rcut)&
!$OMP SHARED(nb_epoxide, nb_pepoxide, nb_ether, nb_pether, nb_palcohols)&
!$OMP SHARED(nb_hydroxide, nb_alcohol, nb_ketone, nb_water, nb_hydronium, nb_alkoxide, nb_oxygen_group)&
!$OMP PRIVATE(s, i, j, k, l, n)&
!$OMP PRIVATE(ij_disp_norm, ij_disp_vec, up, temp_dist, temp_id, temp_el)&
!$OMP PRIVATE(temp_id1, temp_atm_type)
DO s = 1, nb_step
    ! First to assign hydrogen atom (closest atom, one bond only)
    DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 1 ) THEN
            temp_dist = 0.0_dp
            DO j = 1, nb_atm
                IF ( i .EQ. j) CYCLE
                !IF ( atm_mat(2,j,s) .NE. 1 ) THEN
                    DO k = 1, 3
                        ij_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        ij_disp_vec(k) = ij_disp_vec(k) - box(k) * ANINT( ij_disp_vec(k) / box(k) )
                    END DO
                    ij_disp_norm = NORM2( ij_disp_vec )
                    IF ( ( ( ij_disp_norm .LT. 2.00 ) .AND. (temp_dist .NE. 0 ) .AND.&
                     ( ij_disp_norm .LT. temp_dist ) .OR. (temp_dist .EQ. 0 ) ) ) THEN
                        temp_dist = ij_disp_norm
                        temp_id = INT(atm_mat(1,j,s))
                        temp_el = INT(atm_mat(2,j,s))
                    END IF
                !END IF
            END DO
            !IF (temp_dist .EQ. 0 ) CYCLE
            IF ( temp_el .EQ. 16 ) THEN
                atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                atm_mat(10,temp_id,s) = atm_mat(10,temp_id,s) + 1
            ELSE IF ( temp_el .EQ. 12 ) THEN
                atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                atm_mat(10,temp_id,s) = atm_mat(10,temp_id,s) + 1
            ELSE IF ( temp_el .EQ. 1) THEN
                atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                atm_mat(10,temp_id,s) = atm_mat(10,temp_id,s) + 1
            ELSE
                atm_mat(12,i,s) = atm_mat(12,i,s) + 1
                atm_mat(10,temp_id,s) = atm_mat(10,temp_id,s) + 1
            END IF
            ! Fill the connectivity matrix for i
            conn_mat(1,i,s) = atm_mat(1,i,s)
            conn_mat(2,i,s) = atm_mat(2,i,s)
            k=3
            DO WHILE ( k .LE. 15 )
                IF ( conn_mat(k,i,s) .EQ. 0) THEN
                    conn_mat(k,i,s) = atm_mat(1,temp_id,s)
                    conn_mat(k+1,i,s) = atm_mat(2,temp_id,s)
                    conn_mat(k+2,i,s) = temp_dist
                    EXIT
                END IF
                k=k+3
            END DO
            ! Fill the connectivity matrix for j
            conn_mat(1,temp_id,s) = atm_mat(1,temp_id,s)
            conn_mat(2,temp_id,s) = atm_mat(2,temp_id,s)
            k=3
            DO WHILE ( k .LE. 15 )
                IF ( conn_mat(k,temp_id,s) .EQ. 0) THEN
                    conn_mat(k,temp_id,s) = atm_mat(1,i,s)
                    conn_mat(k+1,temp_id,s) = atm_mat(2,i,s)
                    conn_mat(k+2,temp_id,s) = temp_dist
                    EXIT
                END IF
                k=k+3
            END DO
        END IF
    END DO
    ! Second Loop to assign heavy atom between them (ij pair)
 T2i:DO i = 1, nb_atm - 1
        IF ( atm_mat(2,i,s) .EQ. 1 ) CYCLE T2i ! Escape H
    T2j:DO j = i+1, nb_atm
            IF ( atm_mat(2,j,s) .EQ. 1 ) CYCLE T2j ! Escape H
            DO k = 1, 3
                ij_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                ij_disp_vec(k) = ij_disp_vec(k) - box(k) * ANINT( ij_disp_vec(k) / box(k) )
            END DO
            ij_disp_norm = NORM2( ij_disp_vec )
            IF ( ij_disp_norm .LT. 2.00 ) THEN
                IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
                    IF ( ( atm_mat(2,j,s) .EQ. 12 ) .AND. ( ij_disp_norm .LE. assign_OC_rcut) ) THEN
                        IF ( ij_disp_norm .LE. assign_OC_rcut ) THEN
                            atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                            atm_mat(11,j,s) = atm_mat(11,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF ( ( atm_mat(2,j,s) .EQ. 16 ) .AND. ( ij_disp_norm .LE. assign_OO_rcut) ) THEN
                        IF ( ij_disp_norm .LE. assign_OO_rcut ) THEN
                            atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                            atm_mat(11,j,s) = atm_mat(11,j,s) + 1
                            up = 1
                        END IF
                    END IF
                    ! Check if i is Carbon
                ELSE IF ( atm_mat(2,i,s) .EQ. 12 ) THEN
                    IF ( ( atm_mat(2,j,s) .EQ. 12 ) .AND. ( ij_disp_norm .LE. assign_CC_rcut) ) THEN
                        IF ( ij_disp_norm .LE. assign_CC_rcut ) THEN
                            atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                            atm_mat(9,j,s) = atm_mat(9,j,s) + 1
                            up = 1
                        END IF
                    ELSE IF ( ( atm_mat(2,j,s) .EQ. 16 ) .AND. ( ij_disp_norm .LE. assign_OC_rcut) ) THEN
                        IF ( ij_disp_norm .LE. assign_OC_rcut ) THEN
                            atm_mat(11,i,s) = atm_mat(11,i,s) + 1
                            atm_mat(9,j,s) = atm_mat(9,j,s) + 1
                            up = 1
                        END IF
                    END IF
                END IF
                IF ( up .EQ. 1 ) THEN
                    ! Fill the connectivity matrix for i
                    conn_mat(1,i,s) = atm_mat(1,i,s)
                    conn_mat(2,i,s) = atm_mat(2,i,s)
                    k=3
                    DO WHILE ( k .LE. 15 )
                        IF ( conn_mat(k,i,s) .EQ. 0) THEN
                            conn_mat(k,i,s) = atm_mat(1,j,s)
                            conn_mat(k+1,i,s) = atm_mat(2,j,s)
                            conn_mat(k+2,i,s) = ij_disp_norm
                            EXIT
                        END IF
                        k=k+3
                    END DO
                    ! Fill the connectivity matrix for j
                    conn_mat(1,j,s) = atm_mat(1,j,s)
                    conn_mat(2,j,s) = atm_mat(2,j,s)
                    k=3
                    DO WHILE ( k .LE. 15 )
                        IF ( conn_mat(k,j,s) .EQ. 0) THEN
                            conn_mat(k,j,s) = atm_mat(1,i,s)
                            conn_mat(k+1,j,s) = atm_mat(2,i,s)
                            conn_mat(k+2,j,s) = ij_disp_norm
                            EXIT
                        END IF
                        k=k+3
                    END DO
                    up = 0
                END IF
            END IF
        END DO T2j
    END DO T2i
    DO i = 1, nb_atm
        DO j = 4, 18, 3
            IF ( conn_mat(j,i,s) .EQ. 12 ) THEN
                atm_mat(13,i,s) = atm_mat(13,i,s) + 1
            ELSE IF ( conn_mat(j,i,s) .EQ. 1 ) THEN
                atm_mat(14,i,s) = atm_mat(14,i,s) + 1
            ELSE IF ( conn_mat(j,i,s) .EQ. 16 ) THEN
                atm_mat(15,i,s) = atm_mat(15,i,s) + 1
            ELSE IF ( conn_mat(j,i,s) .NE. 0 ) THEN
                atm_mat(16,i,s) = atm_mat(16,i,s) + 1
            END IF
        END DO
        ! Initial Check
        IF ( ( atm_mat(9,i,s) .NE. atm_mat(13,i,s) ) .OR.&
        ( atm_mat(10,i,s) .NE. atm_mat(14,i,s) ) .OR.&
        ( atm_mat(11,i,s) .NE. atm_mat(15,i,s) ) .OR.&
        ( atm_mat(12,i,s) .NE. atm_mat(16,i,s) ) ) THEN
            PRINT'(A28,1X,I4,1X,A8,1X,I10)', "Initial count error on atom:", i, "in step:", s
            print*,conn_mat(:,i,s)
            PRINT*, atm_mat(9,i,s), atm_mat(13,i,s)
            PRINT*, atm_mat(10,i,s), atm_mat(14,i,s)
            PRINT*, atm_mat(11,i,s), atm_mat(15,i,s)
            PRINT*, atm_mat(12,i,s), atm_mat(16,i,s)
        END IF
    END DO
    ! Now conectivity.
    DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
            IF ( ( atm_mat(9,i,s) .EQ. 0 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 1 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OM" ! HO- Oxygen: 0C, 0O, 1H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 0 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 2 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OW" ! Water oxygen: 0C, 0O, 2H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 0 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 3 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OP" ! H3O+ Oxygen: 0C, 0O, 3H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 1 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OA" ! Alkoxy/Ketone Oxygen: 1C, 0O, 0H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 1 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 1 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OB" ! Alcohol Oxygen: 1C, 0O, 1H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 1 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .GT. 1 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OC" ! Protonated Alcohol Oxygen: 1C, 0O, =>1H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 2 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OD" ! Epoxy/Ether Oxygen: 2C, 0O, 0H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 2 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .GT. 1 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "OE" ! Protonated Epoxy/Ether Oxygen: 2C, 0O, =>1H
            ELSE
                atm_type(i,s) = "OX" ! Unknown
            END IF
        ELSE IF ( atm_mat(2,i,s) .EQ. 12 ) THEN
            IF ( ( atm_mat(9,i,s) .EQ. 3 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "C3" ! Carbon: 3C, 0O, 0H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 2 ) .AND. ( atm_mat(11,i,s) .EQ. 0 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "C2" ! Carbon: 2C, 0O, 0H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 3 ) .AND. ( atm_mat(11,i,s) .EQ. 1 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "C3O" ! Carbon: 3C, 1O, 0H
            ELSE IF ( ( atm_mat(9,i,s) .EQ. 2 ) .AND. ( atm_mat(11,i,s) .EQ. 1 ) .AND. ( atm_mat(10,i,s) .EQ. 0 )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "C2O" ! Carbon: 2C, 1O, 0H
            ELSE
                atm_type(i,s) = "CX" ! Unknown
            END IF
        ELSE IF ( atm_mat(2,i,s) .EQ. 1 ) THEN
            IF ( ( ( atm_mat(9,i,s) .EQ. 1 ) .OR. ( atm_mat(11,i,s) .EQ. 1 ) .OR. ( atm_mat(10,i,s) .EQ. 1 ) )&
            .AND. ( atm_mat(12,i,s) .EQ. 0 ) ) THEN
                atm_type(i,s) = "H1"
            ELSE IF ( ( atm_mat(9,i,s) .GT. 1 ) .OR. ( atm_mat(11,i,s) .GT. 1 ) .OR. ( atm_mat(10,i,s) .GT. 1 )&
            .OR. ( atm_mat(12,i,s) .GT. 1 ) ) THEN
                atm_type(i,s) = "H2"
            ELSE IF ( atm_mat(12,i,s) .EQ. 1 ) THEN
                atm_type(i,s) = "HX"
            ELSE
            END IF
        END IF
    END DO

    temp_atm_type(:,:) = "XXX"
    temp_atm_type(:,1) = atm_type(:,s)

    DO l = 1, 100
        k = 0
        DO i = 1, nb_atm
            IF ( temp_atm_type(i,1) .EQ. temp_atm_type(i,2 ) ) THEN
                k = k + 1
            END IF
            temp_atm_type(i,2) = temp_atm_type(i,1)
        END DO
        IF ( k .EQ. nb_atm ) THEN
            DO i = 1, nb_atm
                IF ( (atm_type(i,s) .EQ. "HX" ) .OR.&
                (atm_type(i,s) .EQ. "H2" ) .OR.&
                (atm_type(i,s) .EQ. "CX" ) .OR.&
                (atm_type(i,s) .EQ. "OX" ) ) THEN
                    PRINT*, s, i, atm_type(i,s), atm_mat(13,i,s), atm_mat(14,i,s), atm_mat(15,i,s), atm_mat(16,i,s)
                    PRINT*, conn_mat(:,i,s)
                END IF
            END DO
            EXIT
        END IF
        IF ( l .EQ. 100 ) THEN
            PRINT'(A31,1X,I10)', "Topology not converged in step:", s
            EXIT
        END IF
     AS:DO i = 1, nb_atm
            ! ASSIGN CARBONS LINKED WITH 3 CARBONS
            IF ( (atm_type(i,s) .EQ. "C3" ) .OR.& ! Init
                (atm_type(i,s) .EQ. "C3O" ) .OR.&
                (atm_type(i,s) .EQ. "C3A" ) .OR.& ! Alkoxy
                (atm_type(i,s) .EQ. "C3B" ) .OR.& ! Alcohol
                (atm_type(i,s) .EQ. "C3C" ) .OR.& ! PAlcohol
                (atm_type(i,s) .EQ. "C3D" ) .OR.& ! Ether
                (atm_type(i,s) .EQ. "C3E" ) .OR.& ! PEther
                (atm_type(i,s) .EQ. "C3F" ) .OR.& ! Epoxy 
                (atm_type(i,s) .EQ. "C3G" ) ) THEN ! PEpoxy
                DO j = 4, 18, 3
                    IF ( conn_mat(j,i,s) .EQ. 16 ) THEN
                        temp_id1 = INT( conn_mat(j-1,i,s) )
                        IF ( (atm_type(temp_id1,s) .EQ. "OA" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OA3" ) ) THEN
                            atm_type(temp_id1,s) = "OA3" ! Alkoxy Oxygen (1C+0H)
                            atm_type(i,s) = "C3A" ! Alkoxy Carbon (3C+1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OB" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OB3" ) ) THEN
                            atm_type(temp_id1,s) = "OB3" ! Alcohol Oxygen (1C+1H)
                            atm_type(i,s) = "C3B" ! Alcohol Carbon (3C+1O(H))
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OC" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OC3" ) ) THEN
                            atm_type(temp_id1,s) = "OC3" ! Protonated Alcohol Oxygen (1C+2H)
                            atm_type(i,s) = "C3C" ! Protonated Alcohol Carbon (3C+1O(H2))
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OD" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OD3" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OF3" ) ) THEN
                            atm_type(temp_id1,s) = "OD3" ! Default to Ether Oxygen (2C+0H)
                            atm_type(i,s) = "C3D" ! Default to Ether Carbon (3C+1O)
                            DO k = 4, 18, 3
                                IF ( conn_mat(k,temp_id1,s) .EQ. 12 ) THEN
                                    DO n = 4, 18, 3
                                        IF ( conn_mat(n,i,s) .EQ. 12 ) THEN
                                            IF ( (conn_mat(k-1,temp_id1,s) .EQ. conn_mat(n-1,i,s) ) .AND.&
                                            (conn_mat(k-1,temp_id1,s) .NE. 0 ) .AND.&
                                            (conn_mat(n-1,i,s) .NE. 0 ) ) THEN
                                                atm_type(temp_id1,s) = "OF3" ! Epoxy Oxygen (2C+0H)
                                                atm_type(i,s) = "C3F" ! Epoxy Carbon (3C+1O)
                                            END IF
                                        END IF
                                    END DO
                                END IF
                            END DO
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OE" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OE3" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OG3" ) ) THEN
                            atm_type(temp_id1,s) = "OE3" ! Default to Protonated Ether Oxygen (2C+xH)
                            atm_type(i,s) = "C3E" ! Default to Protonated Ether Carbon (3C+1O(Hx))
                            DO k = 4, 18, 3
                                IF ( conn_mat(k,temp_id1,s) .EQ. 12 ) THEN
                                    DO n = 4, 18, 3
                                        IF ( conn_mat(n,i,s) .EQ. 12 ) THEN
                                            IF ( (conn_mat(k-1,temp_id1,s) .EQ. conn_mat(n-1,i,s) ) .AND.&
                                            (conn_mat(k-1,temp_id1,s) .NE. 0 ) .AND.&
                                            (conn_mat(n-1,i,s) .NE. 0 ) ) THEN
                                                atm_type(temp_id1,s) = "OG3" ! Protonated Epoxy Oxygen (2C+xH)
                                                atm_type(i,s) = "C3G" ! Protonated Epoxy Carbon (3C+1O(Hx))
                                            END IF
                                        END IF
                                    END DO
                                END IF
                            END DO
                        END IF
                    END IF
                END DO
            ! ASSIGN CARBONS
            ELSE IF ( (atm_type(i,s) .EQ. "C2" ) .OR.& ! Init
                (atm_type(i,s) .EQ. "C2O" ) .OR.&
                (atm_type(i,s) .EQ. "C2A" ) .OR.& ! Alkoxy/Ketone
                (atm_type(i,s) .EQ. "C2B" ) .OR.& ! PKetone/Alcohol
                (atm_type(i,s) .EQ. "C2C" ) .OR.& ! PAlcohol 
                (atm_type(i,s) .EQ. "C2D" ) .OR.& ! Ether
                (atm_type(i,s) .EQ. "C2E" ) .OR.& ! PEther
                (atm_type(i,s) .EQ. "C2F" ) .OR.& ! Epoxy
                (atm_type(i,s) .EQ. "C2G" ) ) THEN ! PEpoxy
                DO j = 4, 18, 3
                    IF ( conn_mat(j,i,s) .EQ. 16 ) THEN
                        temp_id1 = INT( conn_mat(j-1,i,s) )
                        IF ( (atm_type(temp_id1,s) .EQ. "OA" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OA2" ) ) THEN
                            atm_type(temp_id1,s) = "OA2" ! Alkoxy/Ketone Oxygen (1C+0H)
                            atm_type(i,s) = "C2A" ! Alkoxy/Ketone Carbon (3C+1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OB" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OB2" ) ) THEN
                            atm_type(temp_id1,s) = "OB2" ! PKetone/Alcohol Oxygen (1C+1H)
                            atm_type(i,s) = "C2B" ! PKetone/Alcohol Carbon (2C+1O(H))
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OC" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OC2" ) ) THEN
                            atm_type(temp_id1,s) = "OC2" ! Protonated Alcohol Oxygen (1C+2H)
                            atm_type(i,s) = "C2C" ! Protonated Alcohol Carbon (3C+1O(H2))
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OD" ) .AND. atm_mat(14,temp_id1,s) .LT. 1 ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OD2" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OF2" ) ) THEN
                            atm_type(temp_id1,s) = "OD2" ! Default to Ether Oxygen (2C+0H)
                            atm_type(i,s) = "C2D" ! Default to Ether Carbon (2C+1O)
                            DO k = 4, 18, 3
                                IF ( conn_mat(k,temp_id1,s) .EQ. 12 ) THEN
                                    DO n = 4, 18, 3
                                        IF ( conn_mat(n,i,s) .EQ. 12 ) THEN
                                            IF ( (conn_mat(k-1,temp_id1,s) .EQ. conn_mat(n-1,i,s) ) .AND.&
                                            (conn_mat(k-1,temp_id1,s) .NE. 0 ) .AND.&
                                            (conn_mat(n-1,i,s) .NE. 0 ) ) THEN
                                                atm_type(temp_id1,s) = "OF2" ! Epoxy Oxygen (2C+0H)
                                                atm_type(i,s) = "C2F" ! Epoxy Carbon (2C+1O)
                                            END IF
                                        END IF
                                    END DO
                                END IF
                            END DO
                        ELSE IF ( ( (atm_type(temp_id1,s) .EQ. "OE" ) .AND. atm_mat(14,temp_id1,s) .GE. 1 ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OE2" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OG2" ) ) THEN
                            atm_type(temp_id1,s) = "OE2" ! Default to Protonated Ether Oxygen (2C+xH)
                            atm_type(i,s) = "C2E" ! Default to Protonated Ether Carbon (2C+1O(Hx) )
                            DO k = 4, 18, 3
                                IF ( conn_mat(k,temp_id1,s) .EQ. 12 ) THEN
                                    DO n = 4, 18, 3
                                        IF ( conn_mat(n,i,s) .EQ. 12 ) THEN
                                            IF ( (conn_mat(k-1,temp_id1,s) .EQ. conn_mat(n-1,i,s) ) .AND.&
                                            (conn_mat(k-1,temp_id1,s) .NE. 0 ) .AND.&
                                            (conn_mat(n-1,i,s) .NE. 0 ) ) THEN
                                                atm_type(temp_id1,s) = "OG2" ! Protonated Epoxy Oxygen (2C+xH)
                                                atm_type(i,s) = "C2G" ! Protonated Epoxy Carbon (2C+1O(Hx) )
                                            END IF
                                        END IF
                                    END DO
                                END IF
                            END DO
                        END IF
                    END IF
                END DO
            ! ASSIGN MONOVALENT H
            ELSE IF ( (atm_type(i,s) .EQ. "H" ) .OR.&
                (atm_type(i,s) .EQ. "H1" ) .OR.&
                (atm_type(i,s) .EQ. "HM" ) .OR.&
                (atm_type(i,s) .EQ. "HW" ) .OR.&
                (atm_type(i,s) .EQ. "HP" ) .OR.&
                (atm_type(i,s) .EQ. "HB" ) .OR.&
                (atm_type(i,s) .EQ. "HC" ) .OR.&
                (atm_type(i,s) .EQ. "HE" ) .OR.&
                (atm_type(i,s) .EQ. "HG" ) ) THEN
                DO j = 4, 18, 3
                    IF ( conn_mat(j,i,s) .EQ. 16 ) THEN
                        temp_id1 = INT( conn_mat(j-1,i,s) )
                        IF ( (atm_type(temp_id1,s) .EQ. "OM" ) ) THEN
                            atm_type(i,s) = "HM" ! OH Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OW" ) ) THEN
                            atm_type(i,s) = "HW" ! H2O Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OP" ) ) THEN
                            atm_type(i,s) = "HP" ! H3O Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OB" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OB2" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OB3" ) ) THEN
                            atm_type(i,s) = "HB" ! Alcohol/PKetone Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OC" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OC2" ) .OR.&
                        (atm_type(temp_id1,s) .EQ. "OC3" ) ) THEN
                            atm_type(i,s) = "HC" ! PAlcohol Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OE2" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OE3" ) ) THEN
                            atm_type(i,s) = "HE" ! Protonated Epoxy Hydrogen (1O)
                        ELSE IF ( (atm_type(temp_id1,s) .EQ. "OG2" ) .OR.&
                            (atm_type(temp_id1,s) .EQ. "OG3" ) ) THEN
                            atm_type(i,s) = "HG" ! Protonated Epoxy Hydrogen (1O)
                        END IF
                    END IF
                END DO
            END IF
            temp_atm_type(i,1) = atm_type(i,s)
        END DO AS
    END DO

    nb_epoxide(s) = COUNT( (atm_type(:,s) .EQ. "OF2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OF3"), DIM=1 )
    nb_pepoxide(s) = COUNT( (atm_type(:,s) .EQ. "OG2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OG3"), DIM=1 )
    nb_ether(s) = COUNT( (atm_type(:,s) .EQ. "OD2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OD3"), DIM=1 )
    nb_pether(s) = COUNT( (atm_type(:,s) .EQ. "OE2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OE3"), DIM=1 )
    nb_alcohol(s) = COUNT( (atm_type(:,s) .EQ. "OB2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OB3"), DIM=1 )
    nb_palcohols(s) = COUNT( (atm_type(:,s) .EQ. "OC2"), DIM=1 ) + COUNT( (atm_type(:,s) .EQ. "OC3"), DIM=1 )
    nb_alkoxide(s) = COUNT( (atm_type(:,s) .EQ. "OA3"), DIM=1 )
    nb_ketone(s) = COUNT( (atm_type(:,s) .EQ. "OA2"), DIM=1 )
    nb_hydroxide(s) = COUNT( (atm_type(:,s) .EQ. "OM"), DIM=1 )
    nb_water(s) = COUNT( (atm_type(:,s) .EQ. "OW"), DIM=1 )
    nb_hydronium(s) = COUNT( (atm_type(:,s) .EQ. "OP"), DIM=1 )
    nb_oxygen_group(s) = nb_epoxide(s) + nb_pepoxide(s) + nb_ether(s) +  nb_pether(s) +&
        nb_alcohol(s) + nb_alkoxide(s) + nb_ketone(s) + nb_hydroxide(s) + nb_hydronium(s) + nb_water(s)

END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Topology:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print counts
start = OMP_get_wtime()

OPEN(UNIT=50, FILE = suffix//"_oxygen_groups_population.txt")
WRITE(50, '(A4,1X,A10,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4)') &
    "Traj", "Step", "EPO", "PEPO", "ETH", "PETH", "PALC", "ALC", "ALK", "AKL2", "H2O", "OH-", "H3O+", "TSum", "TO"
DO s = 1, nb_step
    WRITE(50, '(A4,1X,I10,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4)') &
    suffix, s, nb_epoxide(s), nb_pepoxide(s), nb_ether(s), nb_pether(s), nb_palcohols(s), nb_alcohol(s), nb_alkoxide(s)&
    , nb_ketone(s), nb_water(s), nb_hydroxide(s), nb_hydronium(s), nb_oxygen_group(s), nb_o
END DO
CLOSE(UNIT=50)

DEALLOCATE(nb_epoxide,nb_pepoxide,nb_alcohol,nb_ether,nb_pether,nb_alkoxide,nb_ketone,nb_water&
,nb_hydronium,nb_hydroxide,nb_oxygen_group, nb_palcohols)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Oxygen groups topology output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the xyz and velocities files
start = OMP_get_wtime()

IF ( WRAP_C .EQ. "Y" ) THEN
    OPEN(UNIT=40, FILE = suffix//"_wrapped_"//file_pos)
    IF ( file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_wrapped_"//file_vel)
ELSE
    OPEN(UNIT=40, FILE = suffix//"_nonwrapped_"//file_pos)
    IF ( file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_nonwrapped_"//file_vel)
END IF
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_atm
    IF ( file_vel .NE. '0') WRITE(41,'(I10)') nb_atm
    WRITE(40,'(A10,I10)') "Step nb:", s
    IF ( file_vel .NE. '0') WRITE(41,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_atm
        WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL( atm_type(i,s) ), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        IF ( file_vel .NE. '0') WRITE(41,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL( atm_type(i,s) )&
        , atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
    END DO
END DO
CLOSE(UNIT=40)
IF ( file_vel .NE. '0') CLOSE(UNIT=41)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions/Velocities output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the waterlist (mask indices for surf)
IF ( waterlist .EQ. 1 ) THEN
    start = OMP_get_wtime()

    OPEN(UNIT=42, FILE = suffix//"_waterlist.txt")
    DO s = 1, nb_step
        WRITE(42,'(A14)', ADVANCE="no")"mask = indices "
        DO i = 1, nb_atm
            IF ( atm_type(i,s) .EQ. "OW" ) THEN
                WRITE(42,'(I5)', ADVANCE="no") INT( atm_mat(1,i,s) - 1 )
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