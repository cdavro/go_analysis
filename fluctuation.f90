!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM fluctuation
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
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:)
REAL(dp)                        :: avg_z, CO_disp_vec(3), CO_disp_norm

!   ----------------------------------------------- Counters
INTEGER                         :: i, s, o ,j, k
INTEGER, ALLOCATABLE            :: error_c(:)
INTEGER                         :: error_t

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Fluctuation'
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
file_is=TRIM(file_is)

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, Fluctuation, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Allocate the atm_mat array
ALLOCATE(atm_mat(8,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:))

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

! D -----------------------------------------------
start = OMP_get_wtime()
ALLOCATE(error_c(nb_step))
error_t = 0
error_c(:) = 0

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(fluct_center_atmnb, error_c, fluct_OpC_rcut)&
!$OMP PRIVATE(CO_disp_norm, CO_disp_vec, avg_z)&
!$OMP PRIVATE(s, i, o, j)
DO s = 1, nb_step
    o = 0
    avg_z = 0.0_dp
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            avg_z = avg_z + atm_mat(6,i,s)
            o = o + 1
        END IF
    END DO
    avg_z = avg_z / o
    D1:DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            atm_mat(7,i,s) = atm_mat(6,i,s) - avg_z
         D2:DO j = 1, nb_atm
                IF (atm_mat(3,j,s) .EQ. 10) THEN
                    DO k=1,3
                        CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                    END DO
                    CO_disp_norm = NORM2(CO_disp_vec)
                    IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                        IF (atm_mat(8,i,s) .NE. 0.0) THEN
                            error_c(s) = error_c(s) + 1
                        END IF
                        atm_mat(8,i,s) = 10
                    END IF
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN
                    DO k=1,3
                        CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                    END DO
                    CO_disp_norm = NORM2(CO_disp_vec)
                    IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                        IF (atm_mat(8,i,s) .NE. 0.0) THEN
                            error_c(s) = error_c(s) + 1
                        END IF
                        atm_mat(8,i,s) = 11
                    END IF
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN
                    DO k=1,3
                        CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                    END DO
                    CO_disp_norm = NORM2(CO_disp_vec)
                    IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                        IF (atm_mat(8,i,s) .NE. 0.0) THEN
                            error_c(s) = error_c(s) + 1
                        END IF
                        atm_mat(8,i,s) = 12
                    END IF
                END IF
            END DO D2
        END IF
    END DO D1
END DO
!$OMP END PARALLEL DO

error_t = SUM(error_c(:))

OPEN(UNIT=41, FILE = suffix//"_fluct_avgz_c.txt")
WRITE(41, '(A10,A10,A20,A10)') "Step", "id", "z_fluct", "type"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            WRITE(41, '(I10,I10,E20.7,I10)') s, INT(atm_mat(1,i,s))&
            , atm_mat(7,i,s), INT(atm_mat(8,i,s))
        END IF
    END DO
END DO
CLOSE(UNIT=41)

finish = OMP_get_wtime()
PRINT'(A40,I10)', "Errors:", error_t
PRINT'(A40,F14.2,A20)', "Fluctuation profiles:", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_mat)

END PROGRAM fluctuation