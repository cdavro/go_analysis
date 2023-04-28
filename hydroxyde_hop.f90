!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM proton_hop
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
    REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), is_mat(:,:,:), OH_mat(:,:,:), HOPP(:,:,:)
    CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
    INTEGER, ALLOCATABLE            :: nb_is(:), nb_max_OH(:)
    INTEGER                         :: nb_max_is, nb_o, nb_oh, hopp_c, nb_max_HOPP

    REAL(dp)                        :: ISaX_disp_norm
    REAL(dp)                        :: prev_new_norm

    !   ----------------------------------------------- Counters
    INTEGER                         :: s, i, j, k, o

    !   -----------------------------------------------
    PRINT'(A100)','--------------------------------------------------'&
    ,'--------------------------------------------------'
    PRINT'(A100)', 'Launching Proton Hop'
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
    PRINT'(A100)', 'Run, Proton Hop, Run!'
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'

    !   ----------------------------------------------- Allocate the atm_mat array
    ALLOCATE(atm_mat(20,nb_atm,nb_step))
    ALLOCATE(atm_name(nb_atm,nb_step))
    atm_mat(:,:,:) = 0.0_dp

    ! A ----------------------------------------------- Read positions
    start = OMP_get_wtime()

    CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)

    nb_o = COUNT( atm_mat(2,:,1) .EQ. 16, DIM=1 )
    nb_oh = COUNT( atm_mat(3,:,:) .EQ. 33 )

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"

    IF ( nb_oh .NE. 0) THEN
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

        CALL sb_read_is(file_is,nb_step,box(:),is_mat(:,:,:),nb_is(:) )

        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"

        ! B ----------------------------------------------- Allocate and calculate closest distance between IS and OH atom
        start = OMP_get_wtime()

        ALLOCATE(OH_mat(12,nb_o,nb_step))
        ALLOCATE(nb_max_OH(nb_step))
        OH_mat(:,:,:) = 0.0_dp
        nb_max_OH(:) = 0
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, is_mat, nb_is, nb_step, nb_atm, OH_mat, nb_max_OH)&
        !$OMP PRIVATE(s, i, j, o)&
        !$OMP PRIVATE(ISaX_disp_norm)
        DO s = 1, nb_step
            o = 0
            DO i = 1, nb_atm
                IF (atm_mat(3,i,s) .EQ. 33) THEN ! OH
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
                    o = o + 1
                    OH_mat(1:12,o,s) = atm_mat(1:12,i,s)
                END IF
            END DO
            nb_max_OH(s) = COUNT( OH_mat(1,:,s) .NE. 0, DIM=1 )
        END DO
        !$OMP END PARALLEL DO

        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "Proximity IS and OH+:", finish-start, "seconds elapsed"

        ! B ----------------------------------------------- Calculate the p-hop
        start = OMP_get_wtime()

        ALLOCATE(HOPP(12,MAXVAL(nb_max_OH),nb_step))
        HOPP(:,:,:) = 0.0_dp
        hopp_c = 0
        nb_max_HOPP = 0

        W:DO s = 1, nb_step
            IF ( nb_max_OH(s) .EQ. 0 ) THEN ! Skip when no proton at all in step s
                CYCLE W
            END IF
            X:DO i = 1, nb_max_OH(s) ! Cycle all protons in step s
                IF (OH_mat(1,i,s) .EQ. 0 ) THEN ! Skip empty
                    CYCLE X
                END IF
                IF ( s .EQ. 1 ) THEN ! Init
                    HOPP(2,i,s) = OH_mat(1,i,s)
                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                    hopp_c = hopp_c + 1
                    HOPP(12,i,s) = hopp_c
                ELSE
                    IF ( nb_max_OH(s-1) .EQ. 0 ) THEN ! Init when no proton at all in step s-1
                        HOPP(2,i,s) = OH_mat(1,i,s)
                        HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                        hopp_c = hopp_c + 1
                        HOPP(12,i,s) = hopp_c
                        CYCLE X
                    END IF
                Y:DO j = 1, nb_max_OH(s-1) ! Cycle all protons in step s-1
                        IF (HOPP(2,j,s-1)  .EQ. 0 ) THEN  ! Skip empty
                            CYCLE Y
                        END IF
                        IF ( OH_mat(1,i,s) .EQ. HOPP(2,j,s-1) ) THEN ! If index ID s is same as s-1, update HOPP(s)
                            DO k = 1,10
                                HOPP(k,i,s) = HOPP(k,j,s-1)
                            END DO
                            HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                            HOPP(12,i,s) = HOPP(12,j,s-1)
                            CYCLE X
                        ELSE

                            CALL sb_dist(OH_mat(4:6,i,s),OH_mat(4:6,j,s-1),box,norm_ij=prev_new_norm) ! calculate distance between index s and index s-1

                            IF ( prev_new_norm .LT. 4.00) THEN ! If less than 3 and not in the memory, jump
                                IF ( ( OH_mat(1,i,s) .NE. HOPP(3,j,s-1) ) .AND.&
                                ( OH_mat(1,i,s)  .NE. HOPP(4,j,s-1) ).AND.&
                                ( OH_mat(1,i,s)  .NE. HOPP(5,j,s-1) ).AND.&
                                ( OH_mat(1,i,s)  .NE. HOPP(6,j,s-1) ) ) THEN
                                    HOPP(1,i,s) = HOPP(1,j,s-1) + 1
                                    HOPP(2,i,s) = OH_mat(1,s,s)
                                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                                    DO k = 3,10
                                        HOPP(k,i,s) = HOPP(k-1,j,s-1)
                                    END DO
                                    HOPP(12,i,s) = HOPP(12,j,s-1)
                                    CYCLE X
                                ELSE IF ( OH_mat(1,i,s) .EQ. HOPP(3,j,s-1) ) THEN
                                    HOPP(1,i,s) = HOPP(1,j,s-1) - 1
                                    HOPP(2,i,s) = OH_mat(1,i,s)
                                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                                    DO k = 3,9
                                        HOPP(k,i,s) = HOPP(k+1,j,s-1)
                                    END DO
                                    HOPP(10,i,s) = 0
                                    HOPP(12,i,s) = HOPP(12,j,s-1)
                                    CYCLE X
                                ELSE IF ( OH_mat(1,i,s) .EQ. HOPP(4,j,s-1) ) THEN
                                    HOPP(1,i,s) = HOPP(1,j,s-1) - 2
                                    HOPP(2,i,s) = OH_mat(1,i,s)
                                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                                    DO k = 3,8
                                        HOPP(k,i,s) = HOPP(k+1,j,s-1)
                                    END DO
                                    DO k = 9,10
                                        HOPP(k,i,s) = 0
                                    END DO
                                    HOPP(12,i,s) = HOPP(12,j,s-1)
                                    CYCLE X
                                ELSE IF ( OH_mat(1,i,s) .EQ. HOPP(5,j,s-1) ) THEN
                                    HOPP(1,i,s) = HOPP(1,j,s-1) - 3
                                    HOPP(2,i,s) = OH_mat(1,i,s)
                                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                                    DO k = 3,7
                                        HOPP(k,i,s) = HOPP(k+1,j,s-1)
                                    END DO
                                    DO k = 8,10
                                        HOPP(k,i,s) = 0
                                    END DO
                                    HOPP(12,i,s) = HOPP(12,j,s-1)
                                    CYCLE X
                                ELSE IF ( OH_mat(1,i,s) .EQ. HOPP(6,j,s-1) ) THEN
                                    HOPP(1,i,s) = HOPP(1,j,s-1) - 4
                                    HOPP(2,i,s) = OH_mat(1,i,s)
                                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                                    DO k = 3,6
                                        HOPP(k,i,s) = HOPP(k+1,j,s-1)
                                    END DO
                                    DO k = 7,10
                                        HOPP(k,i,s) = 0
                                    END DO
                                    HOPP(12,i,s) = HOPP(12,j,s-1)
                                    CYCLE X
                                ELSE
                                    PRINT*, "ERROR",s,i, OH_mat(:,i,s)
                                END IF
                            END IF
                        END IF
                    END DO Y
                    ! This part should only happend if there is a new proton
                    HOPP(1,i,s) = 0
                    HOPP(2,i,s) = OH_mat(1,i,s)
                    HOPP(11,i,s) = OH_mat(7,i,s)*OH_mat(8,i,s)
                    hopp_c = hopp_c + 1
                    HOPP(12,i,s) = hopp_c
                    CYCLE X
                END IF
            END DO X
        END DO W
        nb_max_HOPP = INT( MAXVAL( HOPP(12,:,:) ) )
        IF (hopp_c .NE. nb_max_HOPP) THEN
            PRINT*, "ERROR on max hopp", nb_max_HOPP, hopp_c
        END IF

        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "Hopping function:", finish-start, "seconds elapsed"

    !   ----------------------------------------------- Print output
        start = OMP_get_wtime()

        OPEN(UNIT=40, FILE = suffix//"_hop_all.txt")
        WRITE(40,'(A4,1X,A10,1X,A6,1X,A6,1X,A6,1X,A14)')&
        "Traj", "Step", "H(t)","ID","H_ID","dist_ISD"
        DO s = 1, nb_step
            DO i = 1, MAXVAL(nb_max_OH)
                IF (HOPP(2,i,s) .NE. 0) THEN
                    WRITE(40,'(A4,1X,I10,1X,I6,1X,I6,1X,I6,1X,E14.5)') &
                    suffix, s, INT( HOPP(1,i,s) ), INT( HOPP(2,i,s) ), INT( HOPP(12,i,s) ), HOPP(11,i,s)
                END IF
            END DO
        END DO
        CLOSE(UNIT=40)
        DO j = 1, nb_max_HOPP
            OPEN(UNIT=41, FILE = suffix//"_hop_"//fc_itoc(j)//".txt")
            WRITE(41,'(A4,1X,A10,1X,A6,1X,A6,1X,A6,1X,A14)')&
            "Traj", "Step", "H(t)","ID","H_ID","dist_ISD"
            DO s = 1, nb_step
                DO i = 1, MAXVAL(nb_max_OH)
                    IF (HOPP(2,i,s) .NE. 0) THEN
                        IF (HOPP(12,i,s) .EQ. j) THEN
                            WRITE(41,'(A4,1X,I10,1X,I6,1X,I6,1X,I6,1X,E14.5)') &
                            suffix, s, INT( HOPP(1,i,s) ), INT( HOPP(2,i,s) ), INT( HOPP(12,i,s) ), HOPP(11,i,s)
                        END IF
                    END IF
                END DO
            END DO
            CLOSE(UNIT=41)
        END DO
        DO j = 1, nb_max_HOPP
            OPEN(UNIT=42, FILE = suffix//"_hop_"//fc_itoc(j)//"_debug.txt")
            WRITE(42,'(A4,1X,A10,1X,A6,1X,A6,1X,A6,1X,A14&
            &,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6,1X,A6)')&
            "Traj", "Step", "H(t)","ID","H_ID","dist_ISD"&
            &,"ID-1","ID-2","ID-3","ID-4","ID-5","ID-6","ID-7","ID-8"
            DO s = 1, nb_step
                DO i = 1, MAXVAL(nb_max_OH)
                    IF (HOPP(2,i,s) .NE. 0) THEN
                        IF (HOPP(12,i,s) .EQ. j) THEN
                            WRITE(42,'(A4,1X,I10,1X,I6,1X,I6,1X,I6,1X,E14.5)', ADVANCE = "no")&
                            suffix, s, INT( HOPP(1,i,s) ), INT( HOPP(2,i,s) ), INT( HOPP(12,i,s) ), HOPP(11,i,s)
                            DO k= 3,10
                                WRITE(42,'(1X,I6)', ADVANCE = "no")&
                                INT( HOPP(k,i,s) )
                            END DO
                            WRITE(42,'()')
                        END IF
                    END IF
                END DO
            END DO
            CLOSE(UNIT=42)
        END DO

        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "Dist output:", finish-start, "seconds elapsed"
    ELSE
        PRINT'(A40)', "NO OH- EVER"
    END IF

    !   ----------------------------------------------- End
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    PRINT'(A100)', 'The END'

    !   ----------------------------------------------- Deallocate and exit
    DEALLOCATE(atm_mat,atm_name)
    IF (nb_oh .GT. 0) THEN
        DEALLOCATE(is_mat,nb_is,OH_mat,HOPP,nb_max_OH)
    END IF
END PROGRAM proton_hop