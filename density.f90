!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM density
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
INTEGER                         :: CAC

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), is_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
INTEGER, ALLOCATABLE            :: nb_is(:)
INTEGER                         :: nb_max_is
REAL(dp), ALLOCATABLE           :: avg_z(:)
REAL(dp), ALLOCATABLE           :: dens_down(:,:), dens_up(:,:), avg_dens_down(:), avg_dens_up(:)
REAL(dp), ALLOCATABLE           :: dens_down_avgz_c(:,:), avg_dens_down_avgz_c(:)

REAL(dp)                        :: ISaO_disp_vec(3), ISaO_disp_norm
REAL(dp)                        :: r
INTEGER                         :: count_dens_down, count_dens_up, count_dens_down_avgz_c

!   ----------------------------------------------- Counters
INTEGER                         :: s, i, j, k, o

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Density'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*, "No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM( input_file )
CALL READINPUTSUB(input_file)
file_pos=TRIM( file_pos )
file_is=TRIM( file_is )

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, Density, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Allocate the atm_mat array
ALLOCATE(atm_mat(12,nb_atm,nb_step))
ALLOCATE(atm_name(nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)
DEALLOCATE(atm_name) ! Not Used

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

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

CALL sb_read_is(file_is,nb_step,box(:),is_mat(:,:,:),nb_is(:) )

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"

! B ----------------------------------------------- Calculate closest distance between IS and any O atom
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_step, nb_atm)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(ISaO_disp_vec, ISaO_disp_norm)
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. 16 ) THEN
            DO j = 1, nb_is(s)
                IF ( is_mat(4,j,s) .EQ. 1 ) THEN
                    DO k = 1, 3
                        ISaO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                        ISaO_disp_vec(k) = ISaO_disp_vec(k) - box(k) * ANINT( ISaO_disp_vec(k) / box(k) )
                    END DO
                    ISaO_disp_norm = NORM2( ISaO_disp_vec)
                    IF ( (ISaO_disp_norm .LT. atm_mat(7,i,s) ) .OR. atm_mat(7,i,s) .EQ. 0.0 ) THEN
                        atm_mat(7,i,s) = ISaO_disp_norm
                        IF ( atm_mat(6,i,s) .LT. is_mat(3,j,s) ) THEN
                            atm_mat(8,i,s) = -1
                        ELSE
                            atm_mat(8,i,s) = 1
                        END IF
                        atm_mat(9,i,s) = is_mat(5,j,s)
                    END IF
                ELSE IF ( is_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        ISaO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                        ISaO_disp_vec(k) = ISaO_disp_vec(k) - box(k) * ANINT( ISaO_disp_vec(k) / box(k) )
                    END DO
                    ISaO_disp_norm = NORM2( ISaO_disp_vec)
                    IF ( (ISaO_disp_norm .LT. atm_mat(10,i,s) ) .OR. atm_mat(10,i,s) .EQ. 0.0 ) THEN
                        atm_mat(10,i,s) = ISaO_disp_norm
                        IF ( atm_mat(6,i,s) .GT. is_mat(3,j,s) ) THEN
                            atm_mat(11,i,s) = -1
                        ELSE
                            atm_mat(11,i,s) = 1
                        END IF
                        atm_mat(12,i,s) = is_mat(5,j,s)
                    END IF
                END IF
            END DO
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity IS and O atoms:", finish-start, "seconds elapsed"

! C ----------------------------------------------- Density profiles
start = OMP_get_wtime()

! GO-WATER
ALLOCATE(dens_down(dens_step,nb_step))
dens_down(:,:) = 0.0_dp

! AIR WATER
ALLOCATE(dens_up(dens_step,nb_step))
dens_up(:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(dens_rstart, dens_dr, dens_step)&
!$OMP SHARED(dens_down, dens_up)&
!$OMP PRIVATE(s, r, j, i)&
!$OMP PRIVATE(count_dens_down, count_dens_up)
DO s = 1, nb_step
    r = dens_rstart
    DO j = 1, dens_step
        count_dens_down = 0
        count_dens_up = 0
    C1:DO i = 1, nb_atm
            IF ( atm_mat(3,i,s) .NE. 23 ) THEN !Water Oxygen
                CYCLE C1
            END IF
            IF ( (atm_mat(7,i,s)*atm_mat(8,i,s) .GE. r ) .AND. (atm_mat(7,i,s)*atm_mat(8,i,s) .LT. r+dens_dr) ) THEN
                count_dens_down = count_dens_down + 1
            END IF
            IF ( (atm_mat(10,i,s)*atm_mat(11,i,s) .GE. r ) .AND. (atm_mat(10,i,s)*atm_mat(11,i,s) .LT. r+dens_dr) ) THEN
                count_dens_up = count_dens_up + 1
            END IF
        END DO C1
        dens_down(j,s) = ((15.999+ (2*1.00784) ) * count_dens_down ) / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        dens_up(j,s) = ((15.999+ (2*1.00784) ) * count_dens_up ) / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        r = r + dens_dr
    END DO
END DO
!$OMP END PARALLEL DO

ALLOCATE(avg_dens_down(dens_step))
avg_dens_down(:) = 0.0_dp
ALLOCATE(avg_dens_up(dens_step))

DO j = 1, dens_step
    avg_dens_down(j) = SUM(dens_down(j,:) ) / nb_step
    avg_dens_up(j) = SUM(dens_up(j,:) ) / nb_step
END DO

OPEN(UNIT=41, FILE = suffix//"_density_profile_IS_D_step.txt")
OPEN(UNIT=42, FILE = suffix//"_density_profile_IS_U_step.txt")
WRITE(41, '(A4,1X,A10,1X,A10,1X,A10,1X,A14)') "Traj", "Step", ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
WRITE(42, '(A4,1X,A10,1X,A10,1X,A10,1X,A14)') "Traj", "Step", ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
DO s = 1, nb_step
    DO j = 1, dens_step
        WRITE(41,'(A4,1X,I10,1X,F10.3,1X,F10.3,1X,E14.5)') suffix, s, (dens_rstart + (j-1) * dens_dr)&
        , (dens_rstart + j * dens_dr), dens_down(j,s)
        WRITE(42,'(A4,1X,I10,1X,F10.3,1X,F10.3,1X,E14.5)') suffix, s, (dens_rstart + (j-1) * dens_dr)&
        , (dens_rstart + j * dens_dr), dens_up(j,s)
    END DO
END DO
CLOSE(UNIT=41)
CLOSE(UNIT=42)

DEALLOCATE(dens_down,dens_up)

OPEN(UNIT=43, FILE = suffix//"_density_profile_IS_D_avg.txt")
OPEN(UNIT=44, FILE = suffix//"_density_profile_IS_U_avg.txt")


WRITE(43, '(A4,1X,A10,1X,A10,1X,A14)') 'Traj', ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
WRITE(44, '(A4,1X,A10,1X,A10,1X,A14)') 'Traj', ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
DO j = 1, dens_step
    WRITE(43, '(A4,F10.3,1X,F10.3,1X,E14.5)') suffix, (dens_rstart + (j-1) * dens_dr)&
    , (dens_rstart + j * dens_dr), avg_dens_down(j)
    WRITE(44, '(A4,F10.3,1X,F10.3,1X,E14.5)') suffix, (dens_rstart + (j-1) * dens_dr)&
    , (dens_rstart + j * dens_dr), avg_dens_up(j)
END DO
CLOSE(UNIT=43)
CLOSE(UNIT=44)

DEALLOCATE(avg_dens_down,avg_dens_up)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Density profiles (IS):", finish-start, "seconds elapsed"

! D -----------------------------------------------
start = OMP_get_wtime()

! AVGz of sheet
ALLOCATE(dens_down_avgz_c(dens_step,nb_step))
ALLOCATE(avg_z(nb_step))

dens_down_avgz_c(:,:) = 0.0_dp
avg_z(:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(dens_rstart, dens_dr, dens_step, dens_center_atmnb)&
!$OMP SHARED(avg_z, dens_down_avgz_c)&
!$OMP PRIVATE(s, r, j, i, o)&
!$OMP PRIVATE(count_dens_down_avgz_c)
DO s = 1, nb_step
    o = 0
    DO i = 1, nb_atm
        IF ( atm_mat(2,i,s) .EQ. dens_center_atmnb ) THEN
            avg_z(s) = avg_z(s) + atm_mat(6,i,s)
            o = o + 1
        END IF
    END DO
    avg_z(s) = avg_z(s) / o
    r = dens_rstart
    DO j = 1, dens_step
        count_dens_down_avgz_c = 0
    D1:DO i = 1, nb_atm
            IF ( atm_mat(3,i,s) .NE. 26 ) THEN !Water Oxygen
                CYCLE D1
            END IF
            IF ( ((atm_mat(6,i,s) - avg_z(s) ) .GE. r ) .AND. ( (atm_mat(6,i,s) - avg_z(s) ) .LT. r+dens_dr) ) THEN
                count_dens_down_avgz_c = count_dens_down_avgz_c + 1
            END IF
        END DO D1
        dens_down_avgz_c(j,s) = ((15.999+ (2*1.00784) ) * count_dens_down_avgz_c ) &
        / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        r = r + dens_dr
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(avg_z)

ALLOCATE(avg_dens_down_avgz_c(dens_step))
avg_dens_down_avgz_c(:) = 0.0_dp

DO j = 1, dens_step
    avg_dens_down_avgz_c(j) = SUM(dens_down_avgz_c(j,:) ) / nb_step
END DO

OPEN(UNIT=41, FILE = suffix//"_density_profile_avgZC_step.txt")
WRITE(41, '(A4,1X,A10,1X,A10,1X,A10,1X,A14)') "Traj", "Step", ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
DO s = 1, nb_step
    DO j = 1, dens_step
        WRITE(41,'(A4,1X,I10,1X,F10.3,1X,F10.3,1X,E14.5)') suffix, s, (dens_rstart + (j-1) * dens_dr)&
        , (dens_rstart + j * dens_dr), dens_down_avgz_c(j,s)
    END DO
END DO
CLOSE(UNIT=41)

DEALLOCATE(dens_down_avgz_c)

OPEN(UNIT=43, FILE = suffix//"_density_profile_avgZC_avg.txt")
WRITE(43, '(A4,1X,A10,1X,A10,1X,A14)') 'Traj', ">= r (A)", "< r (A)", "ρ/ρ(bulk)"
DO j = 1, dens_step
    WRITE(43, '(A4,1X,F10.3,1X,F10.3,1X,E14.5)') suffix, (dens_rstart + (j-1) * dens_dr)&
    , (dens_rstart + j * dens_dr), avg_dens_down_avgz_c(j)
END DO
CLOSE(UNIT=43)

DEALLOCATE(avg_dens_down_avgz_c)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Density profiles (avgZC):", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_mat,is_mat,nb_is)

END PROGRAM density