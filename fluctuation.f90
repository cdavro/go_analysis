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
REAL(dp)                        :: avg_z

!   ----------------------------------------------- Counters
INTEGER                         :: s, i, o

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

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(fluct_center_atmnb)&
!$OMP PRIVATE(avg_z)&
!$OMP PRIVATE(s, i, o)
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
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            atm_mat(7,i,s) = atm_mat(6,i,s) - avg_z
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

OPEN(UNIT=41, FILE = suffix//"_fluct_avgZC.txt")
WRITE(41, '(A4,1X,A10,1X,A10,1X,A6,1X,A14)') "Traj", "Step", "C_ID", "C_Type", "zFluct"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            WRITE(41, '(A4,1X,I10,1X,I10,1X,I6,1X,E14.5)') suffix, s, INT(atm_mat(1,i,s))&
            , INT(atm_mat(3,i,s)) , atm_mat(7,i,s)
        END IF
    END DO
END DO
CLOSE(UNIT=41)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Fluctuation profiles:", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_mat)

END PROGRAM fluctuation