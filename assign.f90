!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM assign
USE OMP_LIB
USE INPUT_MOD

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_el(:)
REAL(dp), ALLOCATABLE           :: center_avgpos(:,:)
INTEGER, ALLOCATABLE            :: nb_max_OHvec(:)

!   ----------------------------------------------- Temporary variables
REAL(dp)                        :: z_shift
REAL(dp)                        :: OaH_disp_vec(3), OaH_disp_norm
REAL(dp)                        :: OaC_disp_vec(3), OaC_disp_norm
REAL(dp)                        :: OaSi_disp_vec(3), OaSi_disp_norm
CHARACTER(LEN=3)                :: type

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h, nb_c
INTEGER, ALLOCATABLE            :: nb_center(:)
INTEGER, ALLOCATABLE            :: nb_epoxide(:), nb_alcohol(:), nb_alkoxide(:)
INTEGER, ALLOCATABLE            :: nb_water(:), nb_hydronium(:), nb_hydroxide(:)
INTEGER, ALLOCATABLE            :: nb_oxygen_group(:)

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j
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
! DEFINE AS: atm_id, atm_nb, atm_x, atm_y, atm_z, vel_x, vel_y, vel_z, nb_H, nb_C
ALLOCATE(atm_mat(13,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

!   ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_el(nb_atm))

OPEN(UNIT=20, FILE=file_pos, STATUS='old', FORM='formatted', ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i = 1, nb_atm
        READ(20, *) atm_el(i), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        atm_mat(1,i,s) = i
        IF (atm_el(i)(1:1) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
        ELSE IF (atm_el(i)(1:1) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
        ELSE IF (atm_el(i)(1:1) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
        ELSE IF (atm_el(i)(1:2) .EQ. "Si") THEN
            atm_mat(2,i,s) = 28
        END IF
        IF (atm_el(i) .EQ. assign_center_el) THEN
            atm_mat(13,i,s) = 1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)
nb_c = COUNT(atm_mat(2,:,1) .EQ. 12, DIM=1)

ALLOCATE(nb_center(nb_step))
DO s = 1, nb_step
    nb_center(s) = COUNT(atm_mat(13,:,s) .EQ. 1, DIM=1)
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
            READ(21, *) atm_el(i), atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
        END DO
    END DO
    CLOSE(UNIT=21)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Velocities:", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- Calculate geom center of go carbons and center so Zcarbon_avg = 0.0
start = OMP_get_wtime()

ALLOCATE(center_avgpos(3,nb_step))

DO s = 1, nb_step
    center_avgpos(:,s) = 0.0_dp

    DO i = 1, nb_atm
        IF (atm_mat(13,i,s) .EQ. 1) THEN
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

!   ----------------------------------------------- Search the topology, water and oxygen groups
start = OMP_get_wtime()
ALLOCATE(nb_max_OHvec(nb_step))

ALLOCATE(nb_epoxide(nb_step))
ALLOCATE(nb_alcohol(nb_step))
ALLOCATE(nb_alkoxide(nb_step))
ALLOCATE(nb_water(nb_step))
ALLOCATE(nb_hydronium(nb_step))
ALLOCATE(nb_hydroxide(nb_step))
ALLOCATE(nb_oxygen_group(nb_step))

nb_max_OHvec = 0
atm_mat(9,:,:) = -1
atm_mat(10,:,:) = -1
atm_mat(11,:,:) = -1
atm_mat(12,:,:) = -1

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, nb_atm, nb_step, box)&
!$OMP SHARED(assign_OpH_rcut, assign_OpC_rcut, assign_OpSi_rcut)&
!$OMP SHARED(nb_epoxide, nb_hydroxide, nb_alcohol, nb_water, nb_hydronium, nb_alkoxide, nb_oxygen_group)&
!$OMP PRIVATE(s,i,j,k)&
!$OMP PRIVATE(OaC_disp_vec, OaC_disp_norm, OaH_disp_vec, OaH_disp_norm, OaSi_disp_vec, OaSi_disp_norm)
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN ! Select only Oxygens

            atm_mat(9,i,s) = 0
            atm_mat(10,i,s) = 0

            DO j = 1, nb_atm
                IF (atm_mat(2,j,s) .EQ. 1) THEN ! Compare with Hydrogens
                    DO k = 1, 3
                        OaH_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        OaH_disp_vec(k) = OaH_disp_vec(k) - box(k) * ANINT(OaH_disp_vec(k) / box(k))
                    END DO
                    OaH_disp_norm = NORM2(OaH_disp_vec)
                    IF(OaH_disp_norm .LE. assign_OpH_rcut) THEN
                        atm_mat(9,i,s) = atm_mat(9,i,s) + 1
                        atm_mat(11,j,s) = atm_mat(1,i,s)
                    END IF

                ELSE IF (atm_mat(2,j,s) .EQ. 12) THEN ! Compare with Carbons
                    DO k = 1, 3
                        OaC_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        OaC_disp_vec(k) = OaC_disp_vec(k) - box(k) * ANINT(OaC_disp_vec(k) / box(k))
                    END DO
                    OaC_disp_norm = NORM2(OaC_disp_vec)
                    IF(OaC_disp_norm .LE. assign_OpC_rcut) THEN
                        atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                    END IF

                ELSE IF (atm_mat(2,j,s) .EQ. 28) THEN ! Compare with Carbons
                    DO k = 1, 3
                        OaSi_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        OaSi_disp_vec(k) = OaSi_disp_vec(k) - box(k) * ANINT(OaSi_disp_vec(k) / box(k))
                    END DO
                    OaSi_disp_norm = NORM2(OaSi_disp_vec)
                    IF(OaSi_disp_norm .LE. assign_OpSi_rcut) THEN
                        atm_mat(10,i,s) = atm_mat(10,i,s) + 1
                    END IF

                END IF
            END DO
            DO j = 1, nb_atm
                IF (atm_mat(1,i,s) .EQ. atm_mat(11,j,s)) THEN
                    atm_mat(12,j,s) = atm_mat(9,i,s)
                END IF
            END DO
        END IF
    END DO

    nb_epoxide(s) = COUNT((atm_mat(10,:,s) .EQ. 2) .AND. (atm_mat(9,:,s) .EQ. 0), DIM=1)
    nb_alcohol(s) = COUNT((atm_mat(10,:,s) .EQ. 1) .AND. (atm_mat(9,:,s) .EQ. 1), DIM=1)
    nb_alkoxide(s) = COUNT((atm_mat(10,:,s) .EQ. 1) .AND. (atm_mat(9,:,s) .EQ. 0), DIM=1)
    nb_water(s) = COUNT((atm_mat(10,:,s) .EQ. 0) .AND. (atm_mat(9,:,s) .EQ. 2), DIM=1)
    nb_hydroxide(s) = COUNT((atm_mat(10,:,s) .EQ. 0) .AND. (atm_mat(9,:,s) .EQ. 1), DIM=1)
    nb_hydronium(s) = COUNT((atm_mat(10,:,s) .EQ. 0) .AND. (atm_mat(9,:,s) .EQ. 3), DIM=1)
    nb_oxygen_group(s) = nb_epoxide(s) + nb_alcohol(s) + nb_alkoxide(s) + nb_hydroxide(s) + nb_hydronium(s) + nb_water(s)

END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Oxygen groups topologies:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print counts
start = OMP_get_wtime()

OPEN(UNIT=50, FILE = suffix//"_oxygen_groups_population.txt")
WRITE(50, '(A10,A10,A10,A10,A10,A10,A10,A10,A10)') "Step", "Epoxide", "Alcohol", "Alkoxide", "Water"&
, "Hydroxide", "Hydronium", "Total", "Total_O"
DO s = 1, nb_step
    WRITE(50, '(I10,I10,I10,I10,I10,I10,I10,I10,I10)') s, nb_epoxide(s), nb_alcohol(s), nb_alkoxide(s)&
    , nb_water(s),nb_hydroxide(s), nb_hydronium(s), nb_oxygen_group(s), nb_o
END DO
CLOSE(UNIT=50)

DEALLOCATE(nb_max_OHvec,nb_epoxide,nb_alcohol,nb_alkoxide,nb_water,nb_hydronium,nb_hydroxide,nb_oxygen_group)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Oxygen groups topologies output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the xyz and velocities files
start = OMP_get_wtime()

IF (file_vel .NE. '0') OPEN(UNIT=40, FILE = suffix//"_wrapped_"//file_pos)
OPEN(UNIT=41, FILE = suffix//"_wrapped_"//file_vel)
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_atm
    IF (file_vel .NE. '0') WRITE(41,'(I10)') nb_atm
    WRITE(40,'(A10,I10)') "Step nb:", s
    IF (file_vel .NE. '0') WRITE(41,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_atm
        IF ((atm_mat(10,i,s) .EQ. 2) .AND. (atm_mat(9,i,s) .EQ. 0)) THEN
            type = "OE"
        ELSE IF ((atm_mat(10,i,s) .EQ. 1) .AND. (atm_mat(9,i,s) .EQ. 1)) THEN
            type = "OH"
        ELSE IF ((atm_mat(10,i,s) .EQ. 1) .AND. (atm_mat(9,i,s) .EQ. 0)) THEN
            type = "OA"
        ELSE IF ((atm_mat(10,i,s) .EQ. 0) .AND. (atm_mat(9,i,s) .EQ. 2)) THEN
            type = "OW"
        ELSE IF ((atm_mat(10,i,s) .EQ. 0) .AND. (atm_mat(9,i,s) .EQ. 1)) THEN
            type = "OM"
        ELSE IF ((atm_mat(10,i,s) .EQ. 0) .AND. (atm_mat(9,i,s) .EQ. 3)) THEN
            type = "OP"
        ELSE IF ((atm_mat(10,i,s) .EQ. -1) .AND. (atm_mat(12,i,s) .EQ. 2)) THEN
            type = "HW"
        ELSE IF ((atm_mat(10,i,s) .EQ. -1) .AND. (atm_mat(12,i,s) .EQ. 1)) THEN
            type = "HO"
        ELSE IF ((atm_mat(10,i,s) .EQ. -1) .AND. (atm_mat(9,i,s) .EQ. -1)) THEN
            type = atm_el(i)
        ELSE
            type = "O"
        END IF
        WRITE(40,'(A10,E24.14,E24.14,E24.14)') ADJUSTL(type), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        IF (file_vel .NE. '0') WRITE(41,'(A10,E24.14,E24.14,E24.14)') ADJUSTL(type), atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
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
        WRITE(42,'(A14)',advance="no")"mask = indices "
        DO i = 1, nb_atm
            IF (atm_mat(10,i,s) .EQ. 0) THEN
                WRITE(42,'(I5)',advance="no") INT(atm_mat(1,i,s))
            END IF
        END DO
        WRITE(42,*)
    END DO
    CLOSE(UNIT=42)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Waterlist (surf):", finish-start, "seconds elapsed"
END IF

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_el,atm_mat)
END PROGRAM assign