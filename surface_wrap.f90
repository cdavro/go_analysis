!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM surface_wrap
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
REAL(dp), ALLOCATABLE           :: pt_mat(:,:,:)
INTEGER, ALLOCATABLE            :: nb_pt(:)
INTEGER                         :: nb_line, nb_max_pt

!   ----------------------------------------------- Counters
INTEGER                         :: i, s
CHARACTER(LEN=64)               :: dummy
INTEGER                         :: CAC

!   ----------------------------------------------- 
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching Surface Wrap'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*,"No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM(input_file)
CALL READINPUTSUB(input_file)

file_surf=TRIM(file_surf)

!   ----------------------------------------------- Controls
IF (file_surf .EQ. '0') THEN
    PRINT*, "No surface file provided"
    STOP
END IF

!   ----------------------------------------------- 
PRINT'(A100)', 'Run, Surface Wrap, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
start = OMP_get_wtime()

OPEN(UNIT=20,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
nb_line=0_dp
DO
    READ(20,*,IOSTAT=iostatus)
    IF (iostatus .NE. 0) THEN
        EXIT
    ELSE
        nb_line=nb_line+1
    END IF
END DO
REWIND(20)
nb_max_pt=CEILING(1.0*nb_line/nb_step)*2

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS grid:",finish-start,"seconds elapsed"

!   ----------------------------------------------- Allocate function for reading files
! DEFINED AS: pt_x, pt_y, pt_z
ALLOCATE(pt_mat(3,nb_max_pt,nb_step))
ALLOCATE(nb_pt(nb_step))
pt_mat(:,:,:) = 0.0_dp
nb_pt(:) = 0

!   ----------------------------------------------- Read positions
start = OMP_get_wtime()

OPEN(UNIT=20,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
DO s=1,nb_step
    READ(20, *) nb_pt(s)
    READ(20, *) dummy, dummy
    DO i=1,nb_pt(s)
        READ(20, *) dummy, pt_mat(1,i,s), pt_mat(2,i,s), pt_mat(3,i,s), dummy, dummy
    END DO
END DO
CLOSE(UNIT=20)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS:",finish-start,"seconds elapsed"

!   ----------------------------------------------- Wrap and write
start = OMP_get_wtime()

DO s=1,nb_step
    DO i=1,nb_pt(s)
        pt_mat(:,i,s) = pt_mat(:,i,s) - box(:) * ANINT(pt_mat(:,i,s)/box(:))
    END DO
END DO

OPEN(UNIT=40, FILE = suffix//"_wrapped_surf.xyz")
DO s = 1, nb_step
    WRITE(40,'(I10)') nb_pt(s)
    WRITE(40,'(A10,I10)') "Step nb:", s
    DO i = 1, nb_pt(s)
        WRITE(40,'(A10,E24.14,E24.14,E24.14)') "X", pt_mat(1,i,s), pt_mat(2,i,s), pt_mat(3,i,s)
    END DO
END DO
CLOSE(UNIT=40)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Center/Wrap:",finish-start,"seconds elapsed"

!   ----------------------------------------------- End
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
DEALLOCATE(pt_mat,nb_pt)
END PROGRAM surface_wrap