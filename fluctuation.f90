!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM fluctuation
USE OMP_LIB
USE INPUT

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
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
REAL(dp)                        :: avg_z, CO_disp_vec(3),CO_disp_norm

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

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(9,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_el(nb_atm))

OPEN(UNIT=20,FILE=file_pos,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i=1,nb_atm
        atm_mat(2,i,s) = -1
        READ(20, *) atm_el(i), atm_mat(4,i,s), atm_mat(5,i,s), atm_mat(6,i,s)
        atm_mat(1,i,s) = i
        IF (atm_el(i) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
            atm_mat(3,i,s) = 1
        ELSE IF (atm_el(i) .EQ. "OE") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 10
        ELSE IF (atm_el(i) .EQ. "OH") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 11
        ELSE IF (atm_el(i) .EQ. "OA") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 12
        ELSE IF (atm_el(i) .EQ. "OW") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 13
        ELSE IF (atm_el(i) .EQ. "OM") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 14
        ELSE IF (atm_el(i) .EQ. "OP") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = 15
        ELSE IF (atm_el(i) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
            atm_mat(3,i,s) = -1
        ELSE IF (atm_el(i) .EQ. "HW") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 23
        ELSE IF (atm_el(i) .EQ. "HO") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = 21
        ELSE IF (atm_el(i) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
            atm_mat(3,i,s) = -1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

! D -----------------------------------------------
start = OMP_get_wtime()
ALLOCATE(error_c(nb_step))

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(fluct_center_atmnb, error_c)&
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
        atm_mat(7,i,s) = avg_z
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            atm_mat(8,i,s) = atm_mat(6,i,s) - avg_z
         D2:DO j = 1, nb_atm
            IF (atm_mat(2,j,s) .EQ. 10) THEN
                DO k=1,3
                    CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                    CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                END DO
                CO_disp_norm = NORM2(CO_disp_vec)
                IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                    IF (atm_mat(9,i,s) .NE. 0.0) THEN
                        error_c = error_c + 1
                    END IF
                    atm_mat(9,i,s) = 10
                END IF
            ELSE IF (atm_mat(2,j,s) .EQ. 11) THEN
                DO k=1,3
                    CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                    CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                END DO
                CO_disp_norm = NORM2(CO_disp_vec)
                IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                    IF (atm_mat(9,i,s) .NE. 0.0) THEN
                        error_c = error_c + 1
                    END IF
                    atm_mat(9,i,s) = 11
                END IF
            ELSE IF (atm_mat(2,j,s) .EQ. 12) THEN
                DO k=1,3
                    CO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                    CO_disp_vec(k) = CO_disp_vec(k) - box(k) * ANINT(CO_disp_vec(k)/box(k))
                END DO
                CO_disp_norm = NORM2(CO_disp_vec)
                IF (CO_disp_norm .LT. fluct_OpC_rcut) THEN
                    IF (atm_mat(9,i,s) .NE. 0.0) THEN
                        error_c = error_c + 1
                    END IF
                    atm_mat(9,i,s) = 12
                END IF
            END IF
            END DO D2
        END IF
    END DO D1
END DO
!$OMP END PARALLEL DO

error_t = SUM(error_c(:))

OPEN(UNIT=41, FILE = suffix//"_fluct_avgz_c.txt")
WRITE(41, '(A10,A10,A24,A24,A10)') "Step", "id", "avg_z", "z_fluct", "type"
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. fluct_center_atmnb) THEN
            WRITE(41, '(I10,I10,E20.7,E20.7,I10)') s, INT(atm_mat(1,i,s))&
            , atm_mat(7,i,s), atm_mat(8,i,s), INT(atm_mat(9,i,s))
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