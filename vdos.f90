!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM vdos
USE OMP_LIB
USE INPUT
USE SB_GO

IMPLICIT NONE

!   ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

!   ----------------------------------------------- Timings
REAL(dp)                        :: start, finish, start_i, finish_i, avg_timigs
REAL(dp), ALLOCATABLE           :: timings(:)

!   ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file

!   ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), is_mat(:,:,:)
INTEGER                         :: nb_max_is
INTEGER, ALLOCATABLE            :: nb_is(:)

!   ----------------------------------------------- Temp variables
REAL(dp)                        :: atm_vel0_vec(3), atm_vel_vec(3)
REAL(dp)                        :: XO_disp_vec(3), XO_disp_norm
REAL(dp)                        :: SpO_disp_vec(3), SpO_disp_norm
CHARACTER(LEN=2)                :: dummy_char

! ----------------------------------------------- VDOS
INTEGER                         :: mcs, mcsb
REAL(dp), ALLOCATABLE           :: vdos_mat(:,:), vdos_atm(:), vdos_cnt(:)

!   ----------------------------------------------- Count variables

!   ----------------------------------------------- Counters
INTEGER                         :: i, j, k, s, t
INTEGER                         :: CAC

!   -----------------------------------------------
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'
PRINT'(A100)', 'Launching VDOS'
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
file_pos=TRIM(file_pos)
file_vel=TRIM(file_vel)
file_is=TRIM(file_is)

!   ----------------------------------------------- Controls
mcs = INT(mct / timestep_fs)

IF (mcs+1 .GE. nb_step) THEN
    PRINT*, "Error: Max correlation time > trajectory time"
    STOP
END IF

!   -----------------------------------------------
PRINT'(A100)', 'Run, VDOS, Run!'
PRINT'(A100)','--------------------------------------------------'&
,'--------------------------------------------------'

!   ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(25,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()

CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:))

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Read velocities
start = OMP_get_wtime()

OPEN(UNIT=21, FILE=file_vel, STATUS='old', FORM='formatted', ACTION='READ')
DO s = 1, nb_step
    READ(21, *)
    READ(21, *)
    DO i = 1, nb_atm
        READ(21, *) dummy_char, atm_mat(7,i,s), atm_mat(8,i,s), atm_mat(9,i,s)
    END DO
END DO
CLOSE(UNIT=21)

! Velocities are in bohr/aut
DO k = 1, 3
    atm_mat(k+6,:,:) = atm_mat(k+6,:,:) * bohr_to_angstrom / aut_to_fs
END DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Velocities:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN
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
        
        CALL sb_read_is(file_is,nb_step,box(:),is_mat(:,:,:),nb_is(:))
        
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
    END IF

! E ----------------------------------------------- Proximity between functionnal groups and any atom
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, nb_step, nb_atm)&
!$OMP SHARED(hb_X1O_rcut, hb_X2O_rcut)&
!$OMP PRIVATE(s, i, j, k)&
!$OMP PRIVATE(XO_disp_vec, XO_disp_norm)
DO s = 1, nb_step
    E1:DO i = 1, nb_atm
        E2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE E2
            END IF
            IF (atm_mat(1,j,s) .EQ. atm_mat(1,i,s)) THEN
                CYCLE E2
            END IF
            DO k = 1, 3
                XO_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                XO_disp_vec(k) = XO_disp_vec(k) - box(k) * ANINT(XO_disp_vec(k)/box(k))
            END DO
            XO_disp_norm = NORM2(XO_disp_vec)

            IF ( ( (XO_disp_norm .LE. atm_mat(22,i,s)) .OR.&
            (atm_mat(22,i,s) .EQ. 0.0) ) .AND.&
            (atm_mat(3,j,s) .EQ. 1.) ) THEN
                atm_mat(21,i,s) = atm_mat(1,j,s)
                atm_mat(22,i,s) = XO_disp_norm
                atm_mat(23,i,s) = atm_mat(6,j,s)
            END IF
            IF (XO_disp_norm .LE. hb_X1O_rcut) THEN
                IF (atm_mat(2,j,s) .EQ. 12) THEN ! CC
                    atm_mat(10,i,s) = 1
                    atm_mat(14,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OEP
                    atm_mat(11,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OET
                    atm_mat(25,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OH
                    atm_mat(12,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 14) THEN ! OA
                    atm_mat(13,i,s) = 1
                END IF
            ELSE IF (XO_disp_norm .LE. hb_X2O_rcut) THEN
                IF (atm_mat(2,j,s) .EQ. 12) THEN ! CX
                    atm_mat(14,i,s) = 1
                END IF
            END IF
        END DO E2
    END DO E1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and O atoms:", finish-start, "seconds elapsed"

IF (IS_c .EQ. 'Y' ) THEN

    ! X ----------------------------------------------- Calculate closest distance between IS and any atom
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_atm, nb_step)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(SpO_disp_vec, SpO_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            DO j = 1, nb_is(s)
                IF (is_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                        SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT(SpO_disp_vec(k)/box(k))
                    END DO
                    SpO_disp_norm = NORM2(SpO_disp_vec)
                    IF ( (SpO_disp_norm .LT. atm_mat(15,i,s)) .OR. atm_mat(15,i,s) .EQ. 0.0 ) THEN
                        atm_mat(15,i,s) = SpO_disp_norm
                        IF (atm_mat(6,i,s) .LT. is_mat(3,j,s)) THEN
                            atm_mat(16,i,s) = -1
                        ELSE
                            atm_mat(16,i,s) = 1
                        END IF
                        atm_mat(17,i,s) = is_mat(5,j,s)
                    END IF
                ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        SpO_disp_vec(k) = is_mat(k,j,s) - atm_mat(k+3,i,s)
                        SpO_disp_vec(k) = SpO_disp_vec(k) - box(k) * ANINT(SpO_disp_vec(k)/box(k))
                    END DO
                    SpO_disp_norm = NORM2(SpO_disp_vec)
                    IF ( (SpO_disp_norm .LT. atm_mat(18,i,s)) .OR. atm_mat(18,i,s) .EQ. 0.0 ) THEN
                        atm_mat(18,i,s) = SpO_disp_norm
                        IF (atm_mat(6,i,s) .GT. is_mat(3,j,s)) THEN
                            atm_mat(19,i,s) = -1
                        ELSE
                            atm_mat(19,i,s) = 1
                        END IF
                        atm_mat(20,i,s) = is_mat(5,j,s)
                    END IF
                END IF
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS and O:", finish-start, "seconds elapsed"

END IF


! H ----------------------------------------------- VVCF
start = OMP_get_wtime()
mcsb = INT(mctb / timestep_fs) + 1


ALLOCATE(vdos_mat(nb_atm,mcs+1))
ALLOCATE(vdos_cnt(mcs+1))
ALLOCATE(vdos_atm(mcs+1))
ALLOCATE(timings(mcs+1))
timings(:) = 0.0_dp

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(NONE) SHARED(vdos_mat, mcsb, mcs, atm_mat, box, timings, nb_step, nb_atm)&
!$OMP SHARED(vdos_atm, vdos_cnt, timestep_fs)&
!$OMP PRIVATE(t, s, i, k)&
!$OMP PRIVATE(atm_vel0_vec, atm_vel_vec, start_i, finish_i)
DO t = mcsb, mcs+1
    start_i = OMP_get_wtime()
    vdos_mat(:,t) = 0.0_dp
    vdos_atm(t) = 0.0_dp
    vdos_cnt(t) = 0.0_dp
    DO s = 1, nb_step
        IF (s-1+t-1 .LT. nb_step) THEN
         H1:DO i = 1, nb_atm

                !IF ( ( (water_only .EQ. "Y") .AND.& ! Select water only
                !( (atm_mat(3,i,s) .NE. 13) .OR. atm_mat(3,i,s) .NE. 23) ) ) THEN
                !    CYCLE H1
                !ELSE IF ( ( (water_only .EQ. "E") .AND.& ! NOT WATER
                !    ( (atm_mat(3,i,s) .EQ. 13) .OR. atm_mat(3,i,s) .EQ. 23) ) ) THEN
                !    CYCLE H1
                !END IF

                !IF ( (close_c_only .EQ. "Y") .AND.& ! Select OHvec close to the C surface (<=9A)
                !(atm_mat(14,i,s) .NE. 1) ) THEN
                !    CYCLE H1
                !END IF

                !IF ( (up_down_only .EQ. "D") .AND.&
                !(atm_mat(6,i,s) .GT. atm_mat(22,i,s)) ) THEN ! Down only
                !    CYCLE H1
                !ELSE IF ( (up_down_only .EQ. "U") .AND.&
                !(atm_mat(6,i,s) .LT. atm_mat(22,i,s)) ) THEN ! Up only
                !    CYCLE H1
                !END IF

                !IF ((layers .EQ. "Y") .AND. (IS_c .EQ. 'Y' )) THEN

                !    IF ( (atm_mat(15,i,s)*atm_mat(16,i,s) .GT. layer_up) .OR.&
                !    (atm_mat(15,i,s)*atm_mat(16,i,s) .LE. layer_down)) THEN ! Surface distance (go)
                !        CYCLE H1
                !    END IF

                !END IF

                vdos_cnt(t) = vdos_cnt(t) + 1

                DO k = 1, 3
                    atm_vel0_vec(k) =  atm_mat(k+6,i,s)
                    atm_vel_vec(k) =  atm_mat(k+6,i,s+t-1)
                END DO

                vdos_mat(i,t) = vdos_mat(i,t) + DOT_PRODUCT(atm_vel0_vec,atm_vel_vec)

            END DO H1 ! i at 0

        END IF
    END DO ! Time

    vdos_mat(:,t) = vdos_mat(:,t) * 1.0 / (nb_step - (t-1))

    vdos_atm(t) = SUM(vdos_mat(:,t)) * 1.0 / (vdos_cnt(t))

    finish_i = OMP_get_wtime()
    timings(t)=finish_i-start_i

    IF (MODULO(t,25) .EQ. 0) THEN
        PRINT('(I10,A1,I10,E20.7,E20.7,E20.7)'),t,"/",mcs+1, (t-1)*timestep_fs, vdos_atm(t), timings(t)
    END IF

END DO ! Corr
!$OMP END PARALLEL DO

avg_timigs = (SUM(timings(:)) / (mcs+1-mcsb) )
finish = OMP_get_wtime()

PRINT'(A40,F14.2,A20,A20,F14.2)', "Done with VDOS:",finish-start,"seconds elapsed","avg per step:",avg_timigs

start = OMP_get_wtime()

open(UNIT=51, FILE = suffix//"_vdos_xxz.txt")
WRITE(51,'(A20,A20,A20)') "Time (fs)","VDOS","VDOS_cnt"
DO t = mcsb, mcs+1
    WRITE(51,'(E20.7,E20.7,E20.7)') (t-1)*timestep_fs, vdos_atm(t), vdos_cnt(t)
END DO
CLOSE(UNIT=51)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Done with VDOS output:",finish-start,"seconds elapsed"

DEALLOCATE(vdos_mat, timings, vdos_atm, vdos_cnt)

!   ----------------------------------------------- End
PRINT'(A100)', '--------------------------------------------------'&
, '--------------------------------------------------'
PRINT'(A100)', 'The END'

!   ----------------------------------------------- Deallocate and exit
IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat,nb_is)

DEALLOCATE(atm_mat)
END PROGRAM vdos