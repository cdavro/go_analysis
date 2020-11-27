!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM react_event
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
    CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
    INTEGER, ALLOCATABLE            :: event_t(:,:), event_mat(:,:,:)
    INTEGER                         :: event

    !   ----------------------------------------------- Counters
    INTEGER                         :: s, i, k

    !   -----------------------------------------------
    PRINT'(A100)','--------------------------------------------------'&
    ,'--------------------------------------------------'
    PRINT'(A100)', 'Launching React Event'
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

    !   ----------------------------------------------- Controls
    ! To Do

    !   -----------------------------------------------
    PRINT'(A100)', 'Run,  React Event, Run!'
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'

    !   ----------------------------------------------- Allocate the atm_mat array
    ALLOCATE(atm_mat(6,nb_atm,nb_step))
    ALLOCATE(atm_name(nb_atm,nb_step))

    ! A ----------------------------------------------- Read positions
    start = OMP_get_wtime()

    CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"


    ! B ----------------------------------------------- Count reactive events
    start = OMP_get_wtime()

    ALLOCATE(event_t(11,nb_step))
    ALLOCATE(event_mat(3,nb_atm,nb_step))
    event = 0
    event_t(:,:) = 0
    event_mat(:,:,:) = 0
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, nb_step, nb_atm)&
    !$OMP SHARED(event, event_t, event_mat)&
    !$OMP PRIVATE(s, i, k)
    DO s = 2, nb_step
        DO i = 1, nb_atm
            IF ( ( atm_mat(3,i,s) .NE. atm_mat(3,i,s-1) ) .AND. ( atm_mat(2,i,s) .EQ. 16 ) ) THEN
                event = event + 1
                event_t(1,s) = event_t(1,s) + 1
                DO k = 1, 11
                    IF ( atm_mat(3,i,s-1) .EQ. k+24 ) THEN
                        event_t(k+1,s) = event_t(k+1,s) + 1
                    END IF
                END DO
                event_mat(1,i,s) = INT( atm_mat(1,i,s) )
                event_mat(2,i,s) = INT( atm_mat(3,i,s-1) )
                event_mat(3,i,s) = INT( atm_mat(3,i,s) )
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Event time count:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print the closest distance between IS and any Na/Cl/OH/H3O atom
    start = OMP_get_wtime()

    OPEN(UNIT=40, FILE = suffix//"_event_count.txt")
    WRITE(40,'(A4,1X,A10,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6)') &
        "Traj", "Step", "TEvents", "OH", "OH2", "ETH", "PETH", "EPO", "PEPO", "OA3", "OA", "OM", "OW", "OP"
    DO s = 1, nb_step
        WRITE(40,'(A4,I10)', ADVANCE = "no") &
        suffix, s
        DO k = 1, 12
            WRITE(40, '(1X,I6)', ADVANCE = "no")&
            event_t(k,s)
        END DO
        WRITE(40,'()')
    END DO
    CLOSE(UNIT=40)

    OPEN(UNIT=41, FILE = suffix//"_total_event_count.txt")
    WRITE(41,'(A10,1X,A10)') &
        "Total", "Total 5 ts"
    WRITE(41,'(I10,1X,I10)') &
        event
    CLOSE(UNIT=41)

    OPEN(UNIT=42, FILE = suffix//"_event.txt")
    WRITE(42,'(A4,1X,A10,1X,A6,1X,A3,1X,A3)') &
        "Traj", "Step", "Index", "T-1", "T"
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF (event_mat(1,i,s) .NE. 0) THEN
                WRITE(42,'(A4,1X,I10,1X,I6,1X,I3,1X,I3)') &
                suffix, s, event_mat(1,i,s), event_mat(2,i,s), event_mat(3,i,s)
            END IF
        END DO
    END DO
    CLOSE(UNIT=42)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Writing react events:", finish-start, "seconds elapsed"

        !   ----------------------------------------------- End
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    PRINT'(A100)', 'The END'

    !   ----------------------------------------------- Deallocate and exit
    DEALLOCATE(atm_mat)
END PROGRAM react_event