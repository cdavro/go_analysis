!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM dist
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

    REAL(dp)                        :: ISaX_disp_norm, XaX_disp_norm

    !   ----------------------------------------------- Counters
    INTEGER                         :: s, i, j

    !   -----------------------------------------------
    PRINT'(A100)','--------------------------------------------------'&
    ,'--------------------------------------------------'
    PRINT'(A100)', 'Launching Dist'
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
    file_is=TRIM(file_is)

    !   ----------------------------------------------- Controls
    ! To Do

    !   -----------------------------------------------
    PRINT'(A100)', 'Run, Dist, Run!'
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'

    !   ----------------------------------------------- Allocate the atm_mat array
    ALLOCATE(atm_mat(20,nb_atm,nb_step))
    ALLOCATE(atm_name(nb_atm,nb_step))
    atm_mat(:,:,:) = 0.0_dp

    ! A ----------------------------------------------- Read positions
    start = OMP_get_wtime()

    CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

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

    CALL sb_read_is(file_is,nb_step,box(:),is_mat(:,:,:),nb_is(:))

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"

    ! B ----------------------------------------------- Calculate closest distance between IS and any Na/Cl/OH/H3O atom
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_step, nb_atm)&
    !$OMP PRIVATE(s, i, j)&
    !$OMP PRIVATE(ISaX_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF ( ( atm_mat(3,i,s) .EQ. 60 ) .OR. & ! Na
            (atm_mat(3,i,s) .EQ. 61) .OR. & ! Cl
            (atm_mat(3,i,s) .EQ. 33) .OR. & ! OH
            (atm_mat(3,i,s) .EQ. 35) ) THEN ! H3O
                DO j = 1, nb_is(s)
                    IF ( is_mat(4,j,s) .EQ. 1 ) THEN

                        CALL sb_dist(atm_mat(4:6,i,s),is_mat(1:3,j,s),box,norm_ij=ISaX_disp_norm)

                        IF ( ( ISaX_disp_norm .LT. atm_mat(7,i,s) ) .OR. ( atm_mat(7,i,s) .EQ. 0.0 ) ) THEN
                            atm_mat(7,i,s) = ISaX_disp_norm
                            IF ( atm_mat(6,i,s) .LT. is_mat(3,j,s ) ) THEN
                                atm_mat(8,i,s) = -1
                            ELSE
                                atm_mat(8,i,s) = 1
                            END IF
                            atm_mat(9,i,s) = is_mat(5,j,s)
                        END IF
                    ELSE IF ( is_mat(4,j,s) .EQ. 2 ) THEN
                        
                        CALL sb_dist(atm_mat(4:6,i,s),is_mat(1:3,j,s),box,norm_ij=ISaX_disp_norm)

                        IF ( ( ISaX_disp_norm .LT. atm_mat(10,i,s) ) .OR. ( atm_mat(10,i,s) .EQ. 0.0 ) ) THEN
                            atm_mat(10,i,s) = ISaX_disp_norm
                            IF ( atm_mat(6,i,s) .GT. is_mat(3,j,s ) ) THEN
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
    PRINT'(A40,F14.2,A20)', "Proximity IS and charged species:", finish-start, "seconds elapsed"

    ! B ----------------------------------------------- Calculate closest distance between Ions
    start = OMP_get_wtime()

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_step, nb_atm)&
    !$OMP PRIVATE(s, i, j)&
    !$OMP PRIVATE(XaX_disp_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF ( ( atm_mat(3,i,s) .EQ. 60 ) .OR. & ! Na
            (atm_mat(3,i,s) .EQ. 61) .OR. & ! Cl
            (atm_mat(3,i,s) .EQ. 33) .OR. & ! OH
            (atm_mat(3,i,s) .EQ. 35) ) THEN ! H3O
                DO j = 1, nb_atm
                    IF ( i .EQ. j ) CYCLE
                    IF ( ( atm_mat(3,j,s) .EQ. 60 ) .OR. & ! Na
                    ( atm_mat(3,j,s) .EQ. 61 ) .OR. & ! Cl
                    ( atm_mat(3,j,s) .EQ. 33 ) .OR. & ! OH
                    ( atm_mat(3,j,s) .EQ. 35 ) ) THEN ! H3O
                        CALL sb_dist(atm_mat(4:6,i,s),atm_mat(4:6,j,s),box,norm_ij=XaX_disp_norm)
                        IF ( atm_mat(3,j,s) .EQ. 60 ) THEN ! Ion-Na pair
                            IF ( ( XaX_disp_norm .LT. atm_mat(13,i,s) ) .OR. ( atm_mat(13,i,s) .EQ. 0.0 ) ) THEN
                                atm_mat(13,i,s) = XaX_disp_norm
                                atm_mat(14,i,s) = atm_mat(1,j,s)
                            END IF
                        ELSE IF ( atm_mat(3,j,s) .EQ. 61 ) THEN ! Ion-Cl pair
                            IF ( ( XaX_disp_norm .LT. atm_mat(15,i,s) ) .OR. ( atm_mat(15,i,s) .EQ. 0.0 ) ) THEN
                                atm_mat(15,i,s) = XaX_disp_norm
                                atm_mat(16,i,s) = atm_mat(1,j,s)
                            END IF
                        ELSE IF ( atm_mat(3,j,s) .EQ. 33 ) THEN ! Ion-OH pair
                            IF ( ( XaX_disp_norm .LT. atm_mat(17,i,s) ) .OR. ( atm_mat(17,i,s) .EQ. 0.0 ) ) THEN
                                atm_mat(17,i,s) = XaX_disp_norm
                                atm_mat(18,i,s) = atm_mat(1,j,s)
                            END IF
                        ELSE IF ( atm_mat(3,j,s) .EQ. 35 ) THEN ! Ion-H3O pair
                            IF ( ( XaX_disp_norm .LT. atm_mat(19,i,s) ) .OR. ( atm_mat(19,i,s) .EQ. 0.0 ) ) THEN
                                atm_mat(19,i,s) = XaX_disp_norm
                                atm_mat(20,i,s) = atm_mat(1,j,s)
                            END IF
                        END IF
                    END IF
                END DO
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Distance between charged species:", finish-start, "seconds elapsed"

!   ----------------------------------------------- Print output
    start = OMP_get_wtime()

    OPEN(UNIT=40, FILE = suffix//"_anydist.txt")
    WRITE(40,'(A4,1X,A10,1X,A6,1X,A6,1X,A6,1X,A14,1X,A14&
    &,1X,A6,1X,A14,1X,A6,1X,A14&
    &,1X,A6,1X,A14,1X,A6,1X,A14)') &
    "Traj", "Step", "ID", "EL", "TYPE", "dist_ISD", "dist_ISU"&
    , "ID_Na", "dist_Na", "ID_Cl", "dist_Cl", "ID_OH", "dist_OH", "ID_H3O", "dist_H3O"
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF ( ( atm_mat(3,i,s) .EQ. 60 ) .OR. & ! Na
            (atm_mat(3,i,s) .EQ. 61) .OR. & ! Cl
            (atm_mat(3,i,s) .EQ. 33) .OR. & ! OH
            (atm_mat(3,i,s) .EQ. 35) ) THEN ! H3O
                WRITE(40,'(A4,1X,I10,1X,I6,1X,I6,1X,I6,1X,E14.5,1X,E14.5&
                &,1X,I6,1X,E14.5,1X,I6,1X,E14.5&
                &,1X,I6,1X,E14.5,1X,I6,1X,E14.5)') &
                suffix, s, INT( atm_mat(1,i,s) ), INT( atm_mat(2,i,s) ), INT( atm_mat(3,i,s) )&
                , atm_mat(7,i,s)*atm_mat(8,i,s), atm_mat(10,i,s)*atm_mat(11,i,s)&
                , INT( atm_mat(14,i,s) ), atm_mat(13,i,s), INT( atm_mat(16,i,s) ), atm_mat(15,i,s)&
                , INT( atm_mat(18,i,s) ), atm_mat(17,i,s), INT( atm_mat(20,i,s) ), atm_mat(19,i,s)
            END IF
        END DO
    END DO
    CLOSE(UNIT=40)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Dist output:", finish-start, "seconds elapsed"

    !   ----------------------------------------------- End
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    PRINT'(A100)', 'The END'

    !   ----------------------------------------------- Deallocate and exit
    DEALLOCATE(atm_mat,is_mat,nb_is)
END PROGRAM dist