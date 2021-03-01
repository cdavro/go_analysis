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
CHARACTER(LEN=100)              :: input_file
INTEGER                         :: CAC

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), is_mat(:,:,:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
INTEGER                         :: nb_max_is
INTEGER, ALLOCATABLE            :: nb_pt_is(:)

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

IF ( CAC .EQ. 0 ) THEN
    PRINT*, "No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM(input_file)
CALL READINPUTSUB(input_file)
in_file=TRIM(in_file)

!   ----------------------------------------------- Controls
! To Do

!   -----------------------------------------------
PRINT'(A100)', 'Run, Extract, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

IF (is_XYZ .EQ. 'Y') THEN

!   ----------------------------------------------- Allocate the atm_mat array
    ALLOCATE(atm_mat(6,nb_atm,nb_step))
    ALLOCATE(atm_name(nb_atm,nb_step))
    atm_mat(:,:,:) = 0.0_dp

!   ----------------------------------------------- Read positions
    start = OMP_get_wtime()
    print*,in_file
    CALL sb_read_pos_xyz(in_file,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "XYZ File:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Write the xyz file
    start = OMP_get_wtime()
    
    OPEN(UNIT=40, FILE = fc_trim_ext(in_file)//"-"//fc_itoc( stepi )//"-"//fc_itoc( stepf )//".xyz")
    DO s = stepi, stepf
        WRITE(40,'(I10)') nb_atm
        WRITE(40,'(A10,I10)') "Step nb:", s
        DO i = 1, nb_atm
            WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') ADJUSTL( atm_name(i,s) ), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
        END DO
    END DO
    CLOSE(UNIT=40)
    
    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "XYZ output:", finish-start, "seconds elapsed"
    DEALLOCATE(atm_mat,atm_name)

ELSE IF (is_SURF .EQ. 'Y') THEN

!   ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
    start = OMP_get_wtime()

    CALL sb_count_is(in_file,nb_step,nb_max_is)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Read IS
    start = OMP_get_wtime()

    ALLOCATE(is_mat(4,nb_max_is,nb_step))
    ALLOCATE(nb_pt_is(nb_step))
    is_mat(:,:,:) = 0.0_dp
    nb_pt_is(:) = 0

    CALL sb_read_is(in_file,nb_step,box,is_mat(:,:,:),nb_pt_is)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the IS surface
    start = OMP_get_wtime()

    OPEN(UNIT=40, FILE = fc_trim_ext(in_file)//"-"//fc_itoc( stepi )//"-"//fc_itoc( stepf )//".surf")
    DO s = stepi, stepf
        WRITE(40,'(I10)')  nb_pt_is(s)
        WRITE(40,'(A10,I10)') "Step nb:", s
        DO i = 1, nb_pt_is(s)
            IF ( is_mat(4,i,s) .EQ. 1.0 ) THEN
                WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') "XD", is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
            ELSE IF ( is_mat(4,i,s) .EQ. 2.0 ) THEN
                WRITE(40,'(A3,1X,E14.5,1X,E14.5,1X,E14.5)') "XU", is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
            END IF
        END DO
    END DO
    CLOSE(UNIT=40)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS output:", finish-start, "seconds elapsed"
    DEALLOCATE(is_mat,nb_pt_is)

END IF

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit

END PROGRAM extract