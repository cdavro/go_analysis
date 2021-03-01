!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM water_angle
    USE OMP_LIB
    USE INPUT
    USE SB_GO
    
    IMPLICIT NONE
    
    !   ----------------------------------------------- Set Double precision ----------------------------------------------- !
    INTEGER, PARAMETER              :: dp=KIND(0.0d0)
    
    !   ----------------------------------------------- Timings
    REAL(dp)                        :: start,finish
    
    !   ----------------------------------------------- Input files
    CHARACTER(LEN=100)              :: input_file
    
    !   ----------------------------------------------- Infos/properties
    REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), WAT_mat(:,:,:), is_mat(:,:,:), as_mat(:,:,:)
    CHARACTER(LEN=3), ALLOCATABLE   :: atm_name(:,:)
    INTEGER                         :: nb_max_is, nb_max_as
    INTEGER, ALLOCATABLE            :: nb_is(:), nb_as(:)
    
    !   ----------------------------------------------- Temp variables
    REAL(dp)                        :: OpH_disp_norm
    REAL(dp)                        :: XOwat_disp_norm
    REAL(dp)                        :: SpOwat_disp_vec(3), SpOwat_disp_norm
    REAL(dp)                        :: SpS1_disp_vec(3), SpS1_disp_norm, SpS2_disp_vec(3), SpS2_disp_norm
    REAL(dp)                        :: tPSvec_down(3), tPSvec_up(3)
    REAL(dp)                        :: tPSuvec_down(3), tPSuvec_up(3)
    REAL(dp)                        :: WD_vec(3), HpH_disp_vec(3),WD_uvec(3), HpH_disp_uvec(3)
    REAL(dp)                        :: OH1_disp_vec(3), OH1_disp_uvec(3), OH2_disp_vec(3), OH2_disp_uvec(3)
    REAL(dp)                        :: tO_pos_iPS_down(3), tO_pos_iPS_up(3), tS_pos_PS_down(3), tS_pos_PS_up(3)
    REAL(dp)                        :: tOS_disp_oPS_down(3), tOS_disp_oPS_up(3)
    
    !   ----------------------------------------------- Count variables
    INTEGER                         :: nb_o
    INTEGER,ALLOCATABLE             :: nb_max_WAT(:)
    
    !   ----------------------------------------------- Counters
    INTEGER                         :: i, s, k, j, o
    INTEGER                         :: CAC
    
    !   ----------------------------------------------- Parameters
    
    !   -----------------------------------------------
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    PRINT'(A100)', 'Launching Water Angle'
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    
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
    PRINT'(A100)', 'Run, Water Angle, Run!'
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    
    !   ----------------------------------------------- Allocate function for reading files
    ALLOCATE(atm_mat(6,nb_atm,nb_step))
    ALLOCATE(atm_name(nb_atm,nb_step))
    atm_mat(:,:,:) = 0.0_dp
    
    ! A ----------------------------------------------- Read positions
    start = OMP_get_wtime()
    
    CALL sb_read_pos_xyz(file_pos,nb_atm,nb_step,atm_mat(1:6,:,:),atm_name)
    DEALLOCATE(atm_name) ! Not Used
    
    nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
    
    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Positions:", finish-start, "seconds elapsed"
    
    IF (IS_c .EQ. 'Y' ) THEN
    ! A ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
        start = OMP_get_wtime()
    
        CALL sb_count_is(file_is,nb_step,nb_max_is)
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "IS grid:", finish-start, "seconds elapsed"
    
    ! A ----------------------------------------------- Read IS
        start = OMP_get_wtime()
    
        ALLOCATE(is_mat(19,nb_max_is,nb_step))
        ALLOCATE(nb_is(nb_step))
        is_mat(:,:,:) = 0.0_dp
        nb_is(:) = 0
    
        CALL sb_read_is(file_is,nb_step,box,is_mat(1:5,:,:),nb_is)
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "IS:", finish-start, "seconds elapsed"
    END IF
    
    IF (AS_c .EQ. 'Y') THEN
        ! A ----------------------------------------------- Read AS
            start = OMP_get_wtime()
            ALLOCATE(nb_as(nb_step))
            nb_as(:) = 0
            DO s = 1, nb_step
                nb_as(s) = COUNT(atm_mat(2,:,s) .EQ. wa_AS_center_atmnb, DIM=1)
            END DO
            nb_max_as = MAXVAL(nb_as)
    
            ALLOCATE(as_mat(19,nb_max_as,nb_step))
            as_mat(:,:,:) = 0.0_dp
            j = 0
    
            DO s = 1, nb_step
                j = 0
                DO i = 1, nb_atm
                    IF (atm_mat(2,i,s) .EQ. wa_AS_center_atmnb) THEN
                        j = j + 1
                        as_mat(5,j,s) = atm_mat(1,i,s)
                        DO k = 1, 3
                            as_mat(k,j,s) = atm_mat(k+3,i,s)
                        END DO
                    END IF
                END DO
            END DO
    
            finish = OMP_get_wtime()
            PRINT'(A40,F14.2,A20)', "AS:", finish-start, "seconds elapsed"
    END IF
    
    ! B ----------------------------------------------- Water
    start = OMP_get_wtime()
    ALLOCATE(WAT_mat(58,nb_o*3,nb_step))
    ALLOCATE(nb_max_WAT(nb_step))
    
    nb_max_WAT(:) = 0.0
    WAT_mat(:,:,:) = 0.0_dp
    
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_max_WAT, nb_o, nb_step, nb_atm)&
    !$OMP PRIVATE(s, i, j, k, o)&
    !$OMP PRIVATE(OpH_disp_norm)
    DO s = 1, nb_step
        o = 0
        DO i = 1, nb_atm
            IF (atm_mat(3,i,s) .EQ. 34) THEN
                o = o + 1
                WAT_mat(1,o,s) = atm_mat(1,i,s)
                DO k = 1, 3
                    WAT_mat(k+1,o,s) = atm_mat(k+3,i,s)
                END DO
                DO j = 1, nb_atm
                    IF (atm_mat(3,j,s) .EQ. 37) THEN
        
                        CALL sb_dist(atm_mat(4:6,i,s),atm_mat(4:6,j,s),box,OpH_disp_norm)
    
                        IF ( ( (OpH_disp_norm .LT. WAT_mat(6,o,s)) .OR.&
                        (WAT_mat(6,o,s) .EQ. 0.0) ) .AND.&
                        (atm_mat(1,j,s) .NE. WAT_mat(5,o,s)) ) THEN
                            WAT_mat(10,o,s) = WAT_mat(5,o,s)
                            WAT_mat(11,o,s) = WAT_mat(6,o,s)
                            WAT_mat(5,o,s) = atm_mat(1,j,s)
                            WAT_mat(6,o,s) = OpH_disp_norm
                            DO k = 1, 3
                                WAT_mat(k+11,o,s) = WAT_mat(k+6,o,s)
                                WAT_mat(k+6,o,s) =  atm_mat(k+3,j,s)
                            END DO
    
                        ELSE IF ( ( (OpH_disp_norm .LT. WAT_mat(11,o,s)) .OR.&
                            (WAT_mat(6,o,s) .NE. 0.0) .AND.&
                            (WAT_mat(11,o,s) .EQ. 0.0 )) .AND.&
                            (atm_mat(1,j,s) .NE. WAT_mat(10,o,s) ) ) THEN
                            WAT_mat(10,o,s) = atm_mat(1,j,s)
                            WAT_mat(11,o,s) = OpH_disp_norm
                            DO k = 1, 3
                                WAT_mat(k+11,o,s) = atm_mat(k+3,j,s)
                            END DO
    
                        END IF
                    END IF
                END DO
            END IF
        END DO
        DO i=1,nb_o*3
            call sb_dist(WAT_mat(7:9,i,s),WAT_mat(12:14,i,s),box,vec_ij=WAT_mat(21:23,i,s)) ! HH Vector
            WAT_mat(15:17,i,s) = WAT_mat(7:9,i,s) + WAT_mat(21:23,i,s)/2.0 ! Center of HH
            call sb_dist(WAT_mat(2:4,i,s),WAT_mat(15:17,i,s),box,vec_ij=WAT_mat(18:20,i,s)) ! HOH Bisector
        END DO
        nb_max_WAT(s) = COUNT(WAT_mat(1,:,s) .NE. 0, DIM=1)
    END DO
    !$OMP END PARALLEL DO
    
    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "WAT groups:", finish-start, "seconds elapsed"
    
    
    ! C ----------------------------------------------- Proximity between functionnal groups and any WAT groups
    start = OMP_get_wtime()
    
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, nb_step, nb_atm)&
    !$OMP SHARED(wa_X1Owat_rcut, wa_X2Owat_rcut)&
    !$OMP PRIVATE(s, i, j, k)&
    !$OMP PRIVATE(XOwat_disp_norm)
    DO s = 1, nb_step
        C1:DO i = 1, nb_o*3
            IF (WAT_mat(1,i,s) .EQ. 0) THEN
                CYCLE C1
            END IF
            C2:DO j = 1, nb_atm
                IF (atm_mat(3,j,s) .EQ. -1) THEN
                    CYCLE C2
                END IF
                IF (atm_mat(1,j,s) .EQ. WAT_mat(1,i,s)) THEN
                    CYCLE C2
                END IF
    
                CALL sb_dist(WAT_mat(2:4,i,s),atm_mat(4:6,j,s),box,norm_ij=XOwat_disp_norm)
    
                IF ( ( (XOwat_disp_norm .LE. WAT_mat(51,i,s)) .OR.&
                (WAT_mat(51,i,s) .EQ. 0.0) ) .AND.&
                (atm_mat(2,j,s) .EQ. 12) ) THEN
                    WAT_mat(50,i,s) = atm_mat(1,j,s)
                    WAT_mat(51,i,s) = XOwat_disp_norm
                    WAT_mat(52,i,s) = atm_mat(6,j,s)
                END IF
                IF (XOwat_disp_norm .LE. wa_X1Owat_rcut) THEN
                    IF (atm_mat(2,j,s) .EQ. 12) THEN ! Close to any carbon
                        WAT_mat(27,i,s) = 1
                        WAT_mat(28,i,s) = 1
                    ELSE IF ( (atm_mat(3,j,s) .EQ. 27) .OR.& ! Close to an ether oxygen (protonated or not)
                    (atm_mat(3,j,s) .EQ. 28) )THEN
                        WAT_mat(24,i,s) = 1
                    ELSE IF ( (atm_mat(3,j,s) .EQ. 29) .OR.& ! Close to an epoxy oxygen (protonated or not)
                    (atm_mat(3,j,s) .EQ. 30) )THEN
                        WAT_mat(53,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 25) THEN ! Close to an alcohol oxygen
                        WAT_mat(25,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 26) THEN ! Close to a protonated alcohol oxygen
                        WAT_mat(58,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 31) THEN ! OA3 Alkoxy
                        WAT_mat(26,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 32) THEN ! OA2/OA1 Ketone/Alkoxy
                        WAT_mat(54,i,s) = 1
                    ELSE IF ( (atm_mat(3,j,s) .EQ. 60) .OR.& ! NA
                    (atm_mat(3,j,s) .EQ. 61) )THEN ! CLM
                        WAT_mat(55,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 33) THEN ! OH
                        WAT_mat(56,i,s) = 1
                    ELSE IF (atm_mat(3,j,s) .EQ. 35) THEN ! H3O
                        WAT_mat(57,i,s) = 1
                    END IF
                ELSE IF (XOwat_disp_norm .LE. wa_X2Owat_rcut) THEN
                    IF (atm_mat(2,j,s) .EQ. 12) THEN ! Close to any carbon
                        WAT_mat(28,i,s) = 1
                    END IF
                END IF
            END DO C2
        END DO C1
    END DO
    !$OMP END PARALLEL DO
    
    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity WAT groups and FG groups:", finish-start, "seconds elapsed"
    
    !   ----------------------------------------------- IS
    IF (IS_c .EQ. 'Y' ) THEN
        ! D ----------------------------------------------- Calculate WAT angle IS
        start = OMP_get_wtime()
    
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, is_mat, nb_is, nb_step)&
        !$OMP PRIVATE(s, i, j, k)&
        !$OMP PRIVATE(SpOwat_disp_vec, SpOwat_disp_norm, SpS1_disp_vec, SpS1_disp_norm, SpS2_disp_vec, SpS2_disp_norm)&
        !$OMP PRIVATE(tPSuvec_down, tPSuvec_up, tPSvec_down, tPSvec_up, WD_vec, HpH_disp_vec, WD_uvec, HpH_disp_uvec)&
        !$OMP PRIVATE(OH1_disp_vec, OH2_disp_vec, OH1_disp_uvec, OH2_disp_uvec)&
        !$OMP PRIVATE(tO_pos_iPS_down, tO_pos_iPS_up, tS_pos_PS_down, tS_pos_PS_up)&
        !$OMP PRIVATE(tOS_disp_oPS_down, tOS_disp_oPS_up)
        DO s = 1, nb_step
            D2:DO i = 1, nb_o*3
                IF (WAT_mat(1,i,s) .EQ. 0) THEN
                    CYCLE D2
                END IF
                DO j = 1, nb_is(s)
                    IF (is_mat(4,j,s) .EQ. 1) THEN
    
                        CALL sb_dist(WAT_mat(2:4,i,s),is_mat(1:3,j,s),box,norm_ij=SpOwat_disp_norm)
    
                        IF ( (SpOwat_disp_norm .LT. WAT_mat(29,i,s)) .OR. (WAT_mat(29,i,s) .EQ. 0.0) ) THEN
                            WAT_mat(29,i,s) = SpOwat_disp_norm
                            IF (WAT_mat(4,i,s) .LT. is_mat(3,j,s)) THEN
                                WAT_mat(30,i,s) = -1
                            ELSE
                                WAT_mat(30,i,s) = 1
                            END IF
                            WAT_mat(31,i,s) = is_mat(5,j,s)
                        END IF
                    ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                        
                        CALL sb_dist(WAT_mat(2:4,i,s),is_mat(1:3,j,s),box,norm_ij=SpOwat_disp_norm)
    
                        IF ( (SpOwat_disp_norm .LT. WAT_mat(32,i,s)) .OR. (WAT_mat(32,i,s) .EQ. 0.0) ) THEN
                            WAT_mat(32,i,s) = SpOwat_disp_norm
                            IF (WAT_mat(4,i,s) .GT. is_mat(3,j,s)) THEN
                                WAT_mat(33,i,s) = -1
                            ELSE
                                WAT_mat(33,i,s) = 1
                            END IF
                            WAT_mat(34,i,s) = is_mat(5,j,s)
                        END IF
                    END IF
                END DO
    
                IF ((is_mat(19,INT(WAT_mat(31,i,s)),s) .NE. 1) .OR.&
                    (is_mat(19,INT(WAT_mat(34,i,s)),s) .NE. 1)) THEN
    
                    D3:DO j = 1, nb_is(s) ! First one
                        IF ( (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(31,i,s)),s)) .OR.&
                        (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(34,i,s)),s)) ) THEN
                            CYCLE D3
                        END IF
                        IF (is_mat(4,j,s) .EQ. 1) THEN
    
                            CALL sb_dist(is_mat(1:3,INT(WAT_mat(31,i,s)),s),is_mat(1:3,j,s),box &
                            , norm_ij=SpS1_disp_norm,vec_ij=SpS1_disp_vec)
    
                            IF ( ( (SpS1_disp_norm .LT. is_mat(7,INT(WAT_mat(31,i,s)),s)) .OR.&
                            (is_mat(7,INT(WAT_mat(31,i,s)),s) .EQ. 0.0 ) ) .AND.&
                            (is_mat(5,j,s) .NE. is_mat(6,INT(WAT_mat(31,i,s)),s)) ) THEN
                                is_mat(6,INT(WAT_mat(31,i,s)),s) = is_mat(5,j,s)
                                is_mat(7,INT(WAT_mat(31,i,s)),s) = SpS1_disp_norm
                                is_mat(8:10,INT(WAT_mat(31,i,s)),s) = SpS1_disp_vec(:) / SpS1_disp_norm !8,9,10
                            END IF
    
                        ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
    
                            CALL sb_dist(is_mat(1:3,INT(WAT_mat(34,i,s)),s),is_mat(1:3,j,s),box &
                            , norm_ij=SpS1_disp_norm,vec_ij=SpS1_disp_vec)
    
                            IF ( ( (SpS1_disp_norm .LT. is_mat(7,INT(WAT_mat(34,i,s)),s)) .OR.&
                            (is_mat(7,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                            (is_mat(5,j,s) .NE. is_mat(6,INT(WAT_mat(34,i,s)),s)) ) THEN
                                is_mat(6,INT(WAT_mat(34,i,s)),s) = is_mat(5,j,s)
                                is_mat(7,INT(WAT_mat(34,i,s)),s) = SpS1_disp_norm
                                is_mat(8:10,INT(WAT_mat(34,i,s)),s) = SpS1_disp_vec(:) / SpS1_disp_norm !18,19,20
                            END IF
    
                        END IF
                    END DO D3
    
                D4:DO j = 1, nb_is(s)
    
                        IF ( (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(31,i,s)),s)) .OR.&
                        (is_mat(5,j,s) .EQ. is_mat(5,INT(WAT_mat(34,i,s)),s)) .OR. &
                        (is_mat(5,j,s) .EQ. is_mat(6,INT(WAT_mat(31,i,s)),s)) .OR. &
                        (is_mat(5,j,s) .EQ. is_mat(6,INT(WAT_mat(34,i,s)),s)) ) THEN
                            CYCLE D4
                        END IF
    
                        IF (is_mat(4,j,s) .EQ. 1) THEN
    
                            DO k = 1, 3
                                SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(31,i,s)),s)
                                SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                                SpS2_disp_vec(k) = is_mat(k+7,INT(WAT_mat(31,i,s)),s)
                                SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                            END DO
    
                            SpS1_disp_norm = NORM2(SpS1_disp_vec)
                            SpS2_disp_norm = NORM2(SpS2_disp_vec)
                            !CALL sb_dist(is_mat(1:3,INT(WAT_mat(31,i,s)),s),is_mat(1:3,j,s),box &
                            !, norm_ij=SpS1_disp_norm,vec_ij=SpS1_disp_vec)
                            !print*, "NIGGER"
                            !CALL sb_dist((/0._dp,0._dp,0._dp/),is_mat(8:10,INT(WAT_mat(31,i,s)),s),box &
                            !, norm_ij=SpS2_disp_norm,vec_ij=SpS2_disp_vec)
    
                            IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                            (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) )&
                            THEN
                                CYCLE D4
                            END IF
    
                            IF ( ( (SpS1_disp_norm .LT. is_mat(12,INT(WAT_mat(31,i,s)),s)) .OR.&
                            (is_mat(12,INT(WAT_mat(31,i,s)),s) .EQ. 0.0 ) ) .AND.&
                            (is_mat(5,j,s) .NE. is_mat(11,INT(WAT_mat(31,i,s)),s)) ) THEN
    
                                is_mat(11,INT(WAT_mat(31,i,s)),s) = is_mat(5,j,s)
                                is_mat(12,INT(WAT_mat(31,i,s)),s) = SpS1_disp_norm
                                DO k = 1, 3
                                    is_mat(k+12,INT(WAT_mat(31,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                                END DO
    
                            END IF
    
                        ELSE IF (is_mat(4,j,s) .EQ. 2) THEN
                            DO k = 1, 3
                                SpS1_disp_vec(k) = is_mat(k,j,s) - is_mat(k,INT(WAT_mat(34,i,s)),s)
                                SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                                SpS2_disp_vec(k) = is_mat(k+7,INT(WAT_mat(34,i,s)),s)
                                SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                            END DO
    
                            SpS1_disp_norm = NORM2(SpS1_disp_vec)
                            SpS2_disp_norm = NORM2(SpS2_disp_vec)
    
                            IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                            (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) )&
                            THEN
                                CYCLE D4
                            END IF
    
                            IF ( ( (SpS1_disp_norm .LT. is_mat(12,INT(WAT_mat(34,i,s)),s)) .OR.&
                            (is_mat(12,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                            (is_mat(5,j,s) .NE. is_mat(11,INT(WAT_mat(34,i,s)),s)) ) THEN
    
                                is_mat(11,INT(WAT_mat(34,i,s)),s) = is_mat(5,j,s)
                                is_mat(12,INT(WAT_mat(34,i,s)),s) = SpS1_disp_norm
                                DO k = 1, 3
                                    is_mat(k+12,INT(WAT_mat(34,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !18,19,20
                                END DO
    
                            END IF
                    END IF
    
                    END DO D4
    
                    is_mat(16,INT(WAT_mat(31,i,s)),s) =&
                        is_mat(9,INT(WAT_mat(31,i,s)),s) * is_mat(15,INT(WAT_mat(31,i,s)),s) -&
                        is_mat(10,INT(WAT_mat(31,i,s)),s) * is_mat(14,INT(WAT_mat(31,i,s)),s)
                    is_mat(17,INT(WAT_mat(31,i,s)),s) =&
                        is_mat(10,INT(WAT_mat(31,i,s)),s) * is_mat(13,INT(WAT_mat(31,i,s)),s) -&
                        is_mat(8,INT(WAT_mat(31,i,s)),s) * is_mat(15,INT(WAT_mat(31,i,s)),s)
                    is_mat(18,INT(WAT_mat(31,i,s)),s) =&
                        is_mat(8,INT(WAT_mat(31,i,s)),s) * is_mat(14,INT(WAT_mat(31,i,s)),s) -&
                        is_mat(9,INT(WAT_mat(31,i,s)),s) * is_mat(13,INT(WAT_mat(31,i,s)),s)
    
                    is_mat(16,INT(WAT_mat(34,i,s)),s) =&
                        is_mat(9,INT(WAT_mat(34,i,s)),s) * is_mat(15,INT(WAT_mat(34,i,s)),s) -&
                        is_mat(10,INT(WAT_mat(34,i,s)),s) * is_mat(14,INT(WAT_mat(34,i,s)),s)
                    is_mat(17,INT(WAT_mat(34,i,s)),s) =&
                        is_mat(10,INT(WAT_mat(34,i,s)),s) * is_mat(13,INT(WAT_mat(34,i,s)),s) -&
                        is_mat(8,INT(WAT_mat(34,i,s)),s) * is_mat(15,INT(WAT_mat(34,i,s)),s)
                    is_mat(18,INT(WAT_mat(34,i,s)),s) =&
                        is_mat(8,INT(WAT_mat(34,i,s)),s) * is_mat(14,INT(WAT_mat(34,i,s)),s) -&
                        is_mat(9,INT(WAT_mat(34,i,s)),s) * is_mat(13,INT(WAT_mat(34,i,s)),s)
    
                    is_mat(19,INT(WAT_mat(31,i,s)),s) = 1
                    is_mat(19,INT(WAT_mat(34,i,s)),s) = 1
    
                END IF
    
                DO k = 1, 3
                    tPSvec_down(k) = is_mat(k+15,INT(WAT_mat(31,i,s)),s)
                    tPSvec_up(k) = is_mat(k+15,INT(WAT_mat(34,i,s)),s)
                    WD_vec(k) = WAT_mat(k+17,i,s)
                    HpH_disp_vec(k) = WAT_mat(k+20,i,s)
                    OH1_disp_vec(k) = WAT_mat(k+6,i,s) - WAT_mat(k+1,i,s)
                    OH1_disp_vec(k) = OH1_disp_vec(k) - box(k) * ANINT(OH1_disp_vec(k) / box(k))
                    OH2_disp_vec(k) = WAT_mat(k+11,i,s) - WAT_mat(k+1,i,s)
                    OH2_disp_vec(k) = OH2_disp_vec(k) - box(k) * ANINT(OH2_disp_vec(k) / box(k))
                END DO
    
                tPSuvec_down(:) = tPSvec_down(:) / NORM2(tPSvec_down(:))
                tPSuvec_up(:) = tPSvec_up(:) / NORM2(tPSvec_up(:))
                WD_uvec(:) = WD_vec(:) / NORM2(WD_vec(:))
                HpH_disp_uvec(:) = HpH_disp_vec(:) / NORM2(HpH_disp_vec(:))
                OH1_disp_uvec(:) = OH1_disp_vec(:) / NORM2(OH1_disp_vec(:))
                OH2_disp_uvec(:) = OH2_disp_vec(:) / NORM2(OH2_disp_vec(:))
    
                DO k = 1, 3
                    ! To check orientation of the normal isace vector
                    tO_pos_iPS_down(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_down(k)
                    tO_pos_iPS_down(k) = tO_pos_iPS_down(k) - box(k) * ANINT(tO_pos_iPS_down(k) / box(k))
                    tS_pos_PS_down(k) = is_mat(k,INT(WAT_mat(31,i,s)),s) + 0.01*tPSuvec_down(k)
                    tS_pos_PS_down(k) = tS_pos_PS_down(k) - box(k) * ANINT(tS_pos_PS_down(k) / box(k))
                    tOS_disp_oPS_down(k) = tS_pos_PS_down(k) - tO_pos_iPS_down(k)
                    tOS_disp_oPS_down(k) = tOS_disp_oPS_down(k) - box(k) * ANINT(tOS_disp_oPS_down(k) / box(k))
                    tO_pos_iPS_up(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_up(k)
                    tO_pos_iPS_up(k) = tO_pos_iPS_up(k) - box(k) * ANINT(tO_pos_iPS_up(k) / box(k))
                    tS_pos_PS_up(k) = is_mat(k,INT(WAT_mat(34,i,s)),s) + 0.01*tPSuvec_up(k)
                    tS_pos_PS_up(k) = tS_pos_PS_up(k) - box(k) * ANINT(tS_pos_PS_up(k) / box(k))
                    tOS_disp_oPS_up(k) = tS_pos_PS_up(k) - tO_pos_iPS_up(k)
                    tOS_disp_oPS_up(k) = tOS_disp_oPS_up(k) - box(k) * ANINT(tOS_disp_oPS_up(k) / box(k))
                END DO
    
                IF (NORM2(tOS_disp_oPS_down(:)) .GT. WAT_mat(29,i,s)) THEN
                    tPSuvec_down(:) = -1.0 * tPSuvec_down(:)
                END IF
                IF (NORM2(tOS_disp_oPS_up(:)) .GT. WAT_mat(32,i,s)) THEN
                    tPSuvec_up(:) = -1.0 * tPSuvec_up(:)
                END IF
    
                IF (ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                    HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
                END IF
    
                WAT_mat(35,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), WD_uvec(:)))
                WAT_mat(36,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:)))
                WAT_mat(37,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), OH1_disp_uvec(:) ))
                WAT_mat(38,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), OH2_disp_uvec(:) ))
    
                IF (ACOS(DOT_PRODUCT(tPSuvec_up(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                    HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
                END IF
    
                WAT_mat(39,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), WD_uvec(:)))
                WAT_mat(40,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), HpH_disp_uvec(:)))
                WAT_mat(41,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), OH1_disp_uvec(:)))
                WAT_mat(42,i,s) = ACOS(DOT_PRODUCT(tPSuvec_up(:), OH2_disp_uvec(:)))
    
            END DO D2
        END DO
        !$OMP END PARALLEL DO
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "Proximity WAT groups and IS:", finish-start, "seconds elapsed"
    
        ! E ----------------------------------------------- Write WAT angle IS
        start = OMP_get_wtime()
    
        OPEN(UNIT=32, FILE = suffix//"_IS_WA.txt")
        WRITE(32, '(A4,1X,A10,1X,A10,1X,A10,1X,A10,1X&
        &,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X&
        &,A14,1X,A14,1X,A14,1X,A14,1X,A14,1X,A14,1X,A14,1X,A14,1X,A14,1X,A14)')&
            "Traj", "Step", "O_ID", "H1_ID", "H2_ID"&
            , "cCC", "cOEP", "cOET", "cOH", "cOHP", "cOA3", "cOA2", "cION", "cOM", "cH3O", "cCX"&
            , "dist_ISD", "dist_ISU"&
            , "DW/NISD", "DW/NISU"&
            , "HH/NISD", "HH/NISU"&
            , "OH1/NISD", "OH1/NISU"&
            , "OH2/NISD", "OH2/NISD"
    
        DO s = 1, nb_step
            DO i = 1, nb_max_WAT(s)
                WRITE(32, '(A4,1X,I10,1X,I10,1X,I10,1X,I10,1X&
                &,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X&
                &,E14.5,1X,E14.5,1X,E14.5,1X,E14.5,1X,E14.5,1X&
                &,E14.5,1X,E14.5,1X,E14.5,1X,E14.5,1X,E14.5)')&
                suffix, s, INT(WAT_mat(1,i,s)), INT(WAT_mat(5,i,s)), INT(WAT_mat(10,i,s))&
                , INT(WAT_mat(27,i,s)), INT(WAT_mat(53,i,s)), INT(WAT_mat(24,i,s))&
                , INT(WAT_mat(25,i,s)), INT(WAT_mat(58,i,s)), INT(WAT_mat(26,i,s))&
                , INT(WAT_mat(54,i,s)), INT(WAT_mat(55,i,s)), INT(WAT_mat(56,i,s))&
                , INT(WAT_mat(57,i,s)), INT(WAT_mat(28,i,s))&
                , (WAT_mat(29,i,s)*WAT_mat(30,i,s)), (WAT_mat(32,i,s)*WAT_mat(33,i,s))&
                , WAT_mat(35,i,s), WAT_mat(39,i,s)&
                , WAT_mat(36,i,s), WAT_mat(40,i,s)&
                , WAT_mat(37,i,s), WAT_mat(41,i,s)&
                , WAT_mat(38,i,s), WAT_mat(42,i,s)
            END DO
        END DO
        CLOSE(UNIT=32)
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "WAT IS angles output:", finish-start, "seconds elapsed"
    
    END IF
    
    !   ----------------------------------------------- AS
    IF (AS_c .EQ. 'Y' ) THEN
        ! E ----------------------------------------------- Calculate WAT angle AS
        start = OMP_get_wtime()
    
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat, box, WAT_mat, nb_o, as_mat, nb_as, nb_step)&
        !$OMP PRIVATE(s, i, j, k)&
        !$OMP PRIVATE(SpOwat_disp_vec, SpOwat_disp_norm, SpS1_disp_vec, SpS1_disp_norm, SpS2_disp_vec, SpS2_disp_norm)&
        !$OMP PRIVATE(tPSuvec_down, tPSvec_down, WD_vec, HpH_disp_vec , WD_uvec, HpH_disp_uvec)&
        !$OMP PRIVATE(OH1_disp_vec, OH2_disp_vec, OH1_disp_uvec, OH2_disp_uvec)&
        !$OMP PRIVATE(tO_pos_iPS_down, tS_pos_PS_down)&
        !$OMP PRIVATE(tOS_disp_oPS_down)
        DO s = 1, nb_step
            E2:DO i = 1, nb_o*3
                IF (WAT_mat(1,i,s) .EQ. 0) THEN
                    CYCLE E2
                END IF
                DO j = 1, nb_as(s)
                    DO k = 1, 3
                        SpOwat_disp_vec(k) = as_mat(k,j,s) - WAT_mat(k+1,i,s)
                        SpOwat_disp_vec(k) = SpOwat_disp_vec(k) - box(k) * ANINT(SpOwat_disp_vec(k) / box(k))
                    END DO
                    SpOwat_disp_norm = NORM2(SpOwat_disp_vec)
                    IF ( (SpOwat_disp_norm .LT. WAT_mat(43,i,s)) .OR. (WAT_mat(43,i,s) .EQ. 0.0) ) THEN
                        WAT_mat(43,i,s) = SpOwat_disp_norm
                        IF (WAT_mat(4,i,s) .LT. as_mat(3,j,s)) THEN
                            WAT_mat(44,i,s) = -1
                        ELSE
                            WAT_mat(44,i,s) = 1
                        END IF
                        WAT_mat(45,i,s) = as_mat(5,j,s)
                    END IF
                END DO
    
                IF ((as_mat(19,INT(WAT_mat(45,i,s)),s) .NE. 1)) THEN
    
                    E3:DO j = 1, nb_as(s) ! First one
                        IF ( (as_mat(5,j,s) .EQ. as_mat(5,INT(WAT_mat(45,i,s)),s))) THEN
                            CYCLE E3
                        END IF
    
                        DO k = 1, 3
                            SpS1_disp_vec(k) = as_mat(k,j,s) - as_mat(k,INT(WAT_mat(45,i,s)),s)
                            SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                        END DO
                        SpS1_disp_norm = NORM2(SpS1_disp_vec)
    
                        IF ( ( (SpS1_disp_norm .LT. as_mat(7,INT(WAT_mat(45,i,s)),s)) .OR.&
                        (as_mat(7,INT(WAT_mat(45,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (as_mat(5,j,s) .NE. as_mat(6,INT(WAT_mat(45,i,s)),s)) ) THEN
                            as_mat(6,INT(WAT_mat(45,i,s)),s) = as_mat(5,j,s)
                            as_mat(7,INT(WAT_mat(45,i,s)),s) = SpS1_disp_norm
                            DO k = 1, 3
                                as_mat(k+7,INT(WAT_mat(45,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                            END DO
                        END IF
    
                    END DO E3
    
                    E4:DO j = 1, nb_as(s)
    
                        IF ( (as_mat(5,j,s) .EQ. as_mat(5,INT(WAT_mat(45,i,s)),s)) .OR.&
                        (as_mat(5,j,s) .EQ. as_mat(6,INT(WAT_mat(45,i,s)),s)) ) THEN
                            CYCLE E4
                        END IF
    
                            DO k = 1, 3
                                SpS1_disp_vec(k) = as_mat(k,j,s) - as_mat(k,INT(WAT_mat(45,i,s)),s)
                                SpS1_disp_vec(k) = SpS1_disp_vec(k) - box(k) * ANINT(SpS1_disp_vec(k) / box(k))
                                SpS2_disp_vec(k) = as_mat(k+7,INT(WAT_mat(45,i,s)),s)
                                SpS2_disp_vec(k) = SpS2_disp_vec(k) - box(k) * ANINT(SpS2_disp_vec(k) / box(k))
                            END DO
    
                            SpS1_disp_norm = NORM2(SpS1_disp_vec)
                            SpS2_disp_norm = NORM2(SpS2_disp_vec)
    
                            IF ( (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .LT. 0.50).OR.&
                            (ACOS(DOT_PRODUCT(SpS1_disp_vec(:)/SpS1_disp_norm,SpS2_disp_vec(:)/SpS2_disp_norm)) .GT. c_pi-0.50) )&
                            THEN
                                CYCLE E4
                            END IF
    
                            IF ( ( (SpS1_disp_norm .LT. as_mat(12,INT(WAT_mat(45,i,s)),s)) .OR.&
                            (as_mat(12,INT(WAT_mat(45,i,s)),s) .EQ. 0.0 ) ) .AND.&
                            (as_mat(5,j,s) .NE. as_mat(11,INT(WAT_mat(45,i,s)),s)) ) THEN
    
                                as_mat(11,INT(WAT_mat(45,i,s)),s) = as_mat(5,j,s)
                                as_mat(12,INT(WAT_mat(45,i,s)),s) = SpS1_disp_norm
                                DO k = 1, 3
                                    as_mat(k+12,INT(WAT_mat(45,i,s)),s) = SpS1_disp_vec(k) / SpS1_disp_norm !8,9,10
                                END DO
    
                            END IF
    
                    END DO E4
    
                    as_mat(16,INT(WAT_mat(45,i,s)),s) =&
                        as_mat(9,INT(WAT_mat(45,i,s)),s) * as_mat(15,INT(WAT_mat(45,i,s)),s) -&
                        as_mat(10,INT(WAT_mat(45,i,s)),s) * as_mat(14,INT(WAT_mat(45,i,s)),s)
                    as_mat(17,INT(WAT_mat(45,i,s)),s) =&
                        as_mat(10,INT(WAT_mat(45,i,s)),s) * as_mat(13,INT(WAT_mat(45,i,s)),s) -&
                        as_mat(8,INT(WAT_mat(45,i,s)),s) * as_mat(15,INT(WAT_mat(45,i,s)),s)
                    as_mat(18,INT(WAT_mat(45,i,s)),s) =&
                        as_mat(8,INT(WAT_mat(45,i,s)),s) * as_mat(14,INT(WAT_mat(45,i,s)),s) -&
                        as_mat(9,INT(WAT_mat(45,i,s)),s) * as_mat(13,INT(WAT_mat(45,i,s)),s)
    
                    as_mat(19,INT(WAT_mat(45,i,s)),s) = 1
    
                END IF
    
                DO k = 1, 3
                    tPSvec_down(k) = as_mat(k+15,INT(WAT_mat(45,i,s)),s)
                    WD_vec(k) = WAT_mat(k+17,i,s)
                    HpH_disp_vec(k) = WAT_mat(k+20,i,s)
                    OH1_disp_vec(k) = WAT_mat(k+6,i,s) - WAT_mat(k+1,i,s)
                    OH1_disp_vec(k) = OH1_disp_vec(k) - box(k) * ANINT(OH1_disp_vec(k) / box(k))
                    OH2_disp_vec(k) = WAT_mat(k+11,i,s) - WAT_mat(k+1,i,s)
                    OH2_disp_vec(k) = OH2_disp_vec(k) - box(k) * ANINT(OH2_disp_vec(k) / box(k))
                END DO
    
                tPSuvec_down(:) = tPSvec_down(:) / NORM2(tPSvec_down(:))
                WD_uvec(:) = WD_vec(:) / NORM2(WD_vec(:))
                HpH_disp_uvec(:) = HpH_disp_vec(:) / NORM2(HpH_disp_vec(:))
                HpH_disp_uvec(:) = HpH_disp_vec(:) / NORM2(HpH_disp_vec(:))
                OH1_disp_uvec(:) = OH1_disp_vec(:) / NORM2(OH1_disp_vec(:))
                OH2_disp_uvec(:) = OH2_disp_vec(:) / NORM2(OH2_disp_vec(:))
    
                DO k = 1, 3
                    ! To check orientation of the normal isace vector
                    tO_pos_iPS_down(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_down(k)
                    tO_pos_iPS_down(k) = tO_pos_iPS_down(k) - box(k) * ANINT(tO_pos_iPS_down(k) / box(k))
                    tS_pos_PS_down(k) = as_mat(k,INT(WAT_mat(45,i,s)),s) + 0.01*tPSuvec_down(k)
                    tS_pos_PS_down(k) = tS_pos_PS_down(k) - box(k) * ANINT(tS_pos_PS_down(k) / box(k))
                    tOS_disp_oPS_down(k) = tS_pos_PS_down(k) - tO_pos_iPS_down(k)
                    tOS_disp_oPS_down(k) = tOS_disp_oPS_down(k) - box(k) * ANINT(tOS_disp_oPS_down(k) / box(k))
                END DO
    
                IF (NORM2(tOS_disp_oPS_down(:)) .GT. WAT_mat(43,i,s)) THEN
                    tPSuvec_down(:) = -1.0 * tPSuvec_down(:)
                END IF
    
                IF (ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:))) .LT. c_pi/2.0) THEN
                    HpH_disp_uvec(:) = -1.0 * HpH_disp_uvec(:)
                END IF
    
                WAT_mat(46,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), WD_uvec(:)))
                WAT_mat(47,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), HpH_disp_uvec(:)))
                WAT_mat(48,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), OH1_disp_uvec(:)))
                WAT_mat(49,i,s) = ACOS(DOT_PRODUCT(tPSuvec_down(:), OH2_disp_uvec(:)))
    
            END DO E2
        END DO
        !$OMP END PARALLEL DO
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "Proximity WAT groups and AS:", finish-start, "seconds elapsed"
    
        ! E ----------------------------------------------- Write WAT angle AS
        start = OMP_get_wtime()
    
        OPEN(UNIT=33, FILE = suffix//"_AS_WA.txt")
        WRITE(33, '(A4,1X,A10,1X,A10,1X,A10,1X,A10,1X&
            &,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X&
            &,A14,1X,A14,1X,A14,1X,A14,1X,A14)')&
            "Traj", "Step", "O_ID", "H1_ID", "H2_ID"&
            , "cCC", "cOEP", "cOET", "cOH", "cOHP", "cOA3", "cOA2", "cION", "cOM", "cH3O", "cCX"&
            , "dist_AS", "DW/NAS", "HH/NAS", "OH1/NAS", "OH2/NAS"
        DO s = 1, nb_step
            DO i = 1, nb_max_WAT(s)
                WRITE(33, '(A4,1X,I10,1X,I10,1X,I10,1X,I10,1X&
                &,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X&
                &,E14.5,1X,E14.5,1X,E14.5,1X,E14.5,1X,E14.5)')&
                suffix, s, INT(WAT_mat(1,i,s)), INT(WAT_mat(5,i,s)), INT(WAT_mat(10,i,s))&
                , INT(WAT_mat(27,i,s)), INT(WAT_mat(53,i,s)), INT(WAT_mat(24,i,s))&
                , INT(WAT_mat(25,i,s)), INT(WAT_mat(58,i,s)), INT(WAT_mat(26,i,s))&
                , INT(WAT_mat(54,i,s)), INT(WAT_mat(55,i,s)), INT(WAT_mat(56,i,s))&
                , INT(WAT_mat(57,i,s)), INT(WAT_mat(28,i,s))&
                , (WAT_mat(43,i,s)*WAT_mat(44,i,s)), WAT_mat(46,i,s), WAT_mat(47,i,s)&
                , WAT_mat(48,i,s), WAT_mat(49,i,s)
            END DO
        END DO
        CLOSE(UNIT=33)
    
        finish = OMP_get_wtime()
        PRINT'(A40,F14.2,A20)', "WAT AS angles output:", finish-start, "seconds elapsed"
    END IF
    
    !   ----------------------------------------------- End
    PRINT'(A100)', '--------------------------------------------------'&
    , '--------------------------------------------------'
    PRINT'(A100)', 'The END'
    
    !   ----------------------------------------------- Deallocate and exit
    IF (IS_c .EQ. 'Y') DEALLOCATE(is_mat,nb_is)
    IF (AS_c .EQ. 'Y') DEALLOCATE(as_mat,nb_as)
    
    DEALLOCATE(WAT_mat, atm_mat, nb_max_WAT)
    
END PROGRAM water_angle