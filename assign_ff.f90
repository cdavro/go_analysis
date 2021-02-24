!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM assign_ff
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
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:), atm_type(:,:)
REAL(dp), ALLOCATABLE           :: center_avgpos(:,:)

!   ----------------------------------------------- Temporary variables
REAL(dp)                        :: z_shift

!   ----------------------------------------------- Count variables
INTEGER                         :: nb_o
INTEGER, ALLOCATABLE            :: nb_center(:)
INTEGER, ALLOCATABLE            :: nb_epoxide(:), nb_alcohol(:), nb_alkoxide(:)
INTEGER, ALLOCATABLE            :: nb_water(:), nb_hydronium(:), nb_hydroxide(:)
INTEGER, ALLOCATABLE            :: nb_oxygen_group(:), nb_ether(:)
REAL(dp), ALLOCATABLE           :: conn_mat(:,:,:)

!   ----------------------------------------------- Counters
INTEGER                         :: s, i, k
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
ALLOCATE(atm_mat(9,nb_atm,nb_step))
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
        IF ( atm_name(i) .EQ. "126" ) THEN
            atm_mat(2,i,s) = 12
            atm_type(i,s) = "C3E"
        ELSE IF ( atm_name(i) .EQ. "108" ) THEN
            atm_mat(2,i,s) = 12
            atm_type(i,s) = "C3O"
        ELSE IF ( atm_name(i) .EQ. "88" ) THEN
            atm_mat(2,i,s) = 12
            atm_type(i,s) = "CC"
        ELSE IF ( atm_name(i) .EQ. "122" ) THEN
            atm_mat(2,i,s) = 16
            atm_type(i,s) = "OEP"
        ELSE IF ( atm_name(i) .EQ. "109" ) THEN
            atm_mat(2,i,s) = 16
            atm_type(i,s) = "OH3"
        ELSE IF ( atm_name(i) .EQ. "110" ) THEN
            atm_mat(2,i,s) = 1
            atm_type(i,s) = "HO"
        ELSE IF ( atm_name(i) .EQ. "76" ) THEN
            atm_mat(2,i,s) = 16
            atm_type(i,s) = "OW"
        ELSE IF ( atm_name(i) .EQ. "77" ) THEN
            atm_mat(2,i,s) = 1
            atm_type(i,s) = "HW"
        END IF
        IF ( atm_mat(2,i,s) .EQ. assign_center_nb ) THEN
            atm_mat(9,i,s) = 1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

nb_o = COUNT( atm_mat(2,:,1) .EQ. 16, DIM=1 )

ALLOCATE(nb_center(nb_step))
DO s = 1, nb_step
    nb_center(s) = COUNT( atm_mat(9,:,s) .EQ. 1, DIM=1 )
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
            IF ( atm_mat(9,i,s) .EQ. 1 ) THEN
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
ALLOCATE(nb_epoxide(nb_step))
ALLOCATE(nb_alcohol(nb_step))
ALLOCATE(nb_alkoxide(nb_step))
ALLOCATE(nb_water(nb_step))
ALLOCATE(nb_hydronium(nb_step))
ALLOCATE(nb_hydroxide(nb_step))
ALLOCATE(nb_ether(nb_step))
ALLOCATE(nb_oxygen_group(nb_step))

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(nb_step, atm_type)&
!$OMP SHARED(nb_epoxide, nb_hydroxide, nb_alcohol, nb_water, nb_hydronium, nb_alkoxide, nb_oxygen_group)&
!$OMP SHARED(nb_ether)&
!$OMP PRIVATE(s)

DO s = 1, nb_step
    nb_epoxide(s) = COUNT( (atm_type(:,s) .EQ. "OEP"), DIM=1) +  COUNT( (atm_type(:,s) .EQ. "OEQ"), DIM=1 )
    nb_ether(s) = COUNT( (atm_type(:,s) .EQ. "OET"), DIM=1) +  COUNT( (atm_type(:,s) .EQ. "OEU"), DIM=1 )
    nb_alcohol(s) = COUNT( (atm_type(:,s) .EQ. "OH2"), DIM=1) + COUNT( (atm_type(:,s) .EQ. "OH3"), DIM=1 )
    nb_alkoxide(s) = COUNT( (atm_type(:,s) .EQ. "OA2"), DIM=1) + COUNT( (atm_type(:,s) .EQ. "OA3"), DIM=1 )
    nb_water(s) = COUNT( (atm_type(:,s) .EQ. "OW"), DIM=1 )
    nb_hydroxide(s) = COUNT( (atm_type(:,s) .EQ. "OM"), DIM=1 )
    nb_hydronium(s) = COUNT( (atm_type(:,s) .EQ. "OP"), DIM=1 )
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
IF ( WRAP_C .EQ. "Y" ) THEN
    OPEN(UNIT=40, FILE = suffix//"_wrapped_"//fc_trim_ext(file_pos)//".xyz")
    IF ( file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_wrapped_"//fc_trim_ext(file_vel)//".vel")
ELSE
    OPEN(UNIT=40, FILE = suffix//"_nonwrapped_"//fc_trim_ext(file_pos)//".xyz")
    IF ( file_vel .NE. '0') OPEN(UNIT=41, FILE = suffix//"_nonwrapped_"//fc_trim_ext(file_vel)//".vel")
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

END PROGRAM assign_ff