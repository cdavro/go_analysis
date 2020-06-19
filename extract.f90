!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM extract
USE OMP_LIB
USE INPUT
USE SB_GO

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file, out_file
CHARACTER(LEN=10)               :: stepi_s, stepf_s
INTEGER                         :: CAC

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
INTEGER                         :: dotpos

!   ----------------------------------------------- Counters
INTEGER                         :: s, i

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Extract'
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

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, Extract, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Allocate the atm_mat array
ALLOCATE(atm_mat(12,nb_atm,nb_step))
ALLOCATE(atm_name(nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"


!   ----------------------------------------------- Print the xyz and velocities files
start = OMP_get_wtime()

dotpos = SCAN(TRIM(file_pos),".", BACK= .true.)
WRITE (stepi_s, "(I10)") stepi
WRITE (stepf_s, "(I10)") stepf
stepi_s = TRIM(stepi_s)
stepf_s = TRIM(stepf_s)

IF ( dotpos > 0 ) out_file = file_pos(1:dotpos-1)//"-"//stepi_s//"-"//stepf_s//".xyz"

OPEN(UNIT=40, FILE = out_file)
DO s = stepi, stepf
    WRITE(40,'(I10)') nb_atm
    WRITE(40,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_atm
        WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL(atm_name(i,s)), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
    END DO
END DO
CLOSE(UNIT=40)
IF (file_vel .NE. '0') CLOSE(UNIT=41)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions/Velocities output:", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_mat,atm_name)

END PROGRAM extract