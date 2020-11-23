!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM surface_wrap
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
REAL(dp), ALLOCATABLE           :: is_mat(:,:,:)
INTEGER, ALLOCATABLE            :: nb_pt_is(:)
INTEGER                         :: nb_line_is, nb_max_is

!   ----------------------------------------------- Counters
INTEGER                         :: s, i
CHARACTER(LEN=64)               :: dummy
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'Launching Surface Wrap'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF ( CAC .EQ. 0 ) THEN
    PRINT*, "No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM(input_file)
CALL READINPUTSUB(input_file)

file_is=TRIM(file_is)

!   ----------------------------------------------- Controls
IF (file_is .EQ. '0') THEN
    PRINT*, "No surface file provided"
    STOP
END IF

!   -----------------------------------------------
PRINT'(A100)', 'Run, Surface Wrap, Run!'
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'

!   ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
start = OMP_get_wtime()

OPEN(UNIT=20, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
nb_line_is=0_dp
DO
    READ(20,*,IOSTAT=iostatus)
    IF ( iostatus .NE. 0 ) THEN
        EXIT
    ELSE
        nb_line_is = nb_line_is + 1
    END IF
END DO
REWIND(20)
nb_max_is = CEILING(1.0 * nb_line_is / nb_step) * 2

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(is_mat(4,nb_max_is,nb_step))
ALLOCATE(nb_pt_is(nb_step))
is_mat(:,:,:) = 0.0_dp
nb_pt_is(:) = 0

!   ----------------------------------------------- Read positions
start = OMP_get_wtime()

OPEN(UNIT=20, FILE=file_is, STATUS='old', FORM='formatted', ACTION='READ')
DO s = 1, nb_step
    READ(20, *) nb_pt_is(s)
    READ(20, *) dummy, dummy
    DO i = 1, nb_pt_is(s)
        READ(20, *) dummy, is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s), is_mat(4,i,s), dummy
    END DO
END DO
CLOSE(UNIT=20)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Wrap and write
start = OMP_get_wtime()
IF ( WRAP_C .EQ. "Y" ) THEN
    DO s = 1, nb_step
        DO i = 1, nb_pt_is(s)
            is_mat(:,i,s) = is_mat(:,i,s) - box(:) * ANINT( is_mat(:,i,s) / box(:) )
        END DO
    END DO
END IF

IF ( WRAP_C .EQ. "Y") THEN
    OPEN(UNIT=40, FILE = suffix//"_wrapped_surf.xyz")
ELSE 
    OPEN(UNIT=40, FILE = suffix//"_nonwrapped_surf.xyz")
END IF
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_pt_is(s)
    WRITE(40,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_pt_is(s)
        IF ( is_mat(4,i,s) .GT. 0.0 ) THEN
            WRITE(40, '(A4,1X,E14.5,1X,E14.5,1X,E14.5)') "XD", is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
        ELSE IF ( is_mat(4,i,s) .LT. 0.0 ) THEN
            WRITE(40, '(A4,1X,E14.5,1X,E14.5,1X,E14.5)') "XU", is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
        ELSE
            WRITE(40, '(A4,1X,E14.5,1X,E14.5,1X,E14.5)') "X0", is_mat(1,i,s), is_mat(2,i,s), is_mat(3,i,s)
        END IF
    END DO
END DO
CLOSE(UNIT=40)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Center/Wrap:", finish-start, "seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(is_mat,nb_pt_is)

END PROGRAM surface_wrap