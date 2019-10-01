PROGRAM vvcf
USE OMP_LIB

IMPLICIT NONE

! ----------------------------------------------- Set Double precision ----------------------------------------------- !
INTEGER, PARAMETER              :: dp=KIND(0.0d0)
INTEGER                         :: iostatus

! ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

! ----------------------------------------------- Filenames
CHARACTER(LEN=64)               :: file_xyz, file_vel, file_surf

! ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), OHvec_mat(:,:,:), surf_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
INTEGER                         :: nb_line, nb_max_pt
INTEGER, ALLOCATABLE            :: nb_surf(:)

! ----------------------------------------------- Temp variables
REAL(dp)                        :: tOH_disp_vec(3), tOH_norm, tXOHvec_disp_vec(3), tXOHvec_norm
REAL(dp)                        :: tSOHvec_disp_vec(3), tSOHvec_norm
REAL(dp)                        :: tSS_disp_vec(3), tSS_norm
REAL(dp)                        :: tPS_uvec_go(3), tPS_OOHvec_pos_vec_go(3), tPS_HHvec_pos_vec(3)
REAL(dp)                        :: tPS_uvec_air(3), tPS_OOHvec_pos_vec_air(3)
REAL(dp)                        :: tOHvec_disp_vec(3), tOHvec_disp_uvec(3),tHHvec_disp_vec(3)
CHARACTER(LEN=2)                :: dummy_char

! ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h
INTEGER,ALLOCATABLE             :: nb_max_OHvec(:)

! ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o
INTEGER                         :: CAC

! ----------------------------------------------- Parameters
REAL(dp), PARAMETER             :: bohr_to_angstrom=0.529177249_dp
REAL(dp), PARAMETER             :: aut_to_fs=0.0241888432658569977_dp

! ----------------------------------------------- SYSTEM DEPENDANT
INTEGER, PARAMETER              :: nb_atm=1039, nb_step=1000
CHARACTER(LEN=2)                :: suffix="00"
REAL(dp), PARAMETER             :: xlo=0.0_dp,xhi=21.8489966560_dp
REAL(dp), PARAMETER             :: ylo=0.0_dp,yhi=21.2373561788_dp
REAL(dp), PARAMETER             :: zlo=0.0_dp,zhi=70.0_dp
REAL(dp), PARAMETER             :: rOH_cut=1.3, rOHbond_cut=2.5, rXOH_cut=3.5, rXO_cut=3.5
REAL(dp), PARAMETER             :: rCOH_cut=9.0, rCO_cut=9.0 ! Close to Carbon (interface)
REAL(dp), PARAMETER             :: surf_up=12345, surf_down=67890
REAL(dp)                        :: box(3)
INTEGER, PARAMETER              :: hbond_output = 1
INTEGER, PARAMETER              :: layers = 1


! ----------------------------------------------- Allocate function for reading files
! DEFINE AS: atm_id, atm_nb, atm_x, atm_y, atm_z, vel_x, vel_y, vel_z, nb_OH
ALLOCATE(atm_mat(23,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()
file_surf="0"

IF (CAC .EQ. 0) THEN
    PRINT*,"NO ARGUMENT"
    STOP
ELSE IF (CAC .EQ. 1) THEN
    PRINT*, "Need positions and velocties files"
    STOP
ELSE IF (CAC .EQ. 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,file_surf)
    file_surf=TRIM(file_surf)
END IF

CALL GET_COMMAND_ARGUMENT(1,file_xyz)
file_xyz=TRIM(file_xyz)
CALL GET_COMMAND_ARGUMENT(2,file_vel)
file_vel=TRIM(file_vel)
IF ( (file_surf .EQ. "0") .AND. (layers .EQ. 0) ) THEN
    PRINT*, "Need surface to vvcf layers"
END IF

! ----------------------------------------------- Calculate the box size
box(1) = xhi - xlo
box(2) = yhi - ylo
box(3) = zhi - zlo

! A ----------------------------------------------- Read positions
start = OMP_get_wtime()
ALLOCATE(atm_el(nb_atm))

OPEN(UNIT=20,FILE=file_xyz,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(20, *)
    READ(20, *)
    DO i=1,nb_atm
        atm_mat(2,i,s) = -1
        READ(20, *) atm_el(i), atm_mat(3,i,s), atm_mat(4,i,s), atm_mat(5,i,s)
        atm_mat(1,i,s) = i
        IF (atm_el(i) .EQ. "C") THEN
            atm_mat(2,i,s) = 12
            atm_mat(9,i,s) = 1
        ELSE IF (atm_el(i) .EQ. "OE") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 10
        ELSE IF (atm_el(i) .EQ. "OH") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 11
        ELSE IF (atm_el(i) .EQ. "OA") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 12
        ELSE IF (atm_el(i) .EQ. "OW") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 13
        ELSE IF (atm_el(i) .EQ. "OM") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 14
        ELSE IF (atm_el(i) .EQ. "OP") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = 15
        ELSE IF (atm_el(i) .EQ. "O") THEN
            atm_mat(2,i,s) = 16
            atm_mat(9,i,s) = -1
        ELSE IF (atm_el(i) .EQ. "HW") THEN
            atm_mat(2,i,s) = 1
            atm_mat(9,i,s) = 23
        ELSE IF (atm_el(i) .EQ. "HO") THEN
            atm_mat(2,i,s) = 1
            atm_mat(9,i,s) = 21
        ELSE IF (atm_el(i) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
            atm_mat(9,i,s) = -1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:",finish-start,"seconds elapsed"

! A ----------------------------------------------- Read velocities
start = OMP_get_wtime()

OPEN(UNIT=21,FILE=file_vel,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(21, *)
    READ(21, *)
    DO i=1,nb_atm
        READ(21, *) atm_el(i), atm_mat(6,i,s), atm_mat(7,i,s), atm_mat(8,i,s)
    END DO
END DO
CLOSE(UNIT=21)

! Velocities are in bohr/aut
DO k = 1, 3
    atm_mat(k+5,:,:) = atm_mat(k+5,:,:) * bohr_to_angstrom / aut_to_fs
END DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Velocities:",finish-start,"seconds elapsed"

! ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
IF (file_surf .NE. "0") THEN
    start = OMP_get_wtime()

    OPEN(UNIT=22,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
    nb_line=0
    DO
        READ(22,*,IOSTAT=iostatus)
        IF (iostatus .NE. 0) THEN
            EXIT
        ELSE
            nb_line=nb_line+1
        END IF
    END DO
    REWIND(22)
    nb_max_pt=CEILING(1.0*nb_line/nb_step)*2

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS grid:",finish-start,"seconds elapsed"
END IF

! A ----------------------------------------------- Read surface
IF (file_surf .NE. "0") THEN

    ALLOCATE(surf_mat(32,nb_max_pt,nb_step))
    ALLOCATE(nb_surf(nb_step))
    surf_mat(:,:,:) = 0.0_dp
    nb_surf(:) = 0

    start = OMP_get_wtime()

    OPEN(UNIT=22,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
    DO s = 1, nb_step
        READ(22, *) nb_surf(s)
        READ(22, *)
        j = 0
        DO i=1,nb_surf(s)
            READ(22, *) dummy_char, surf_mat(1,i,s), surf_mat(2,i,s), surf_mat(3,i,s)
            j = j + 1
            surf_mat(5,i,s) = j
            IF (surf_mat(3,i,s) .LT. 10.0) THEN
                surf_mat(4,i,s) = 1
            ELSE
                surf_mat(4,i,s) = 2
            END IF
            DO k = 1, 3
                surf_mat(k,i,s) = surf_mat(k,i,s) - box(k) * ANINT(surf_mat(k,i,s)/box(k))
            END DO
        END DO
    END DO
    CLOSE(UNIT=22)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "IS:",finish-start,"seconds elapsed"
END IF

! ----------------------------------------------- Get number of relevent atom types
nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

! B ----------------------------------------------- OH groups and corresponding values
start = OMP_get_wtime()
ALLOCATE(OHvec_mat(45,nb_o*3,nb_step))
ALLOCATE(nb_max_OHvec(nb_step))

nb_max_OHvec(:) = 0.0
OHvec_mat(:,:,:) = 0.0_dp

!nb_step, nb_atm, rOH_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_max_OHvec,nb_o)&
!$OMP PRIVATE(s,i,j,k,o,tOH_disp_vec,tOH_norm)
DO s = 1, nb_step
    o = 0
    DO i=1,nb_atm
        IF (atm_mat(9,i,s) .EQ. 13) THEN
            o = o + 1
            OHvec_mat(1,o,s) = atm_mat(1,i,s)
            DO k=1,3
                OHvec_mat(k+1,o,s) = atm_mat(k+2,i,s)
            END DO
            DO j=1,nb_atm
                IF (atm_mat(9,j,s) .EQ. 23) THEN
                    DO k = 1, 3
                        tOH_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        tOH_disp_vec(k) = tOH_disp_vec(k) - box(k) * ANINT(tOH_disp_vec(k)/box(k))
                    END DO
                    tOH_norm = NORM2(tOH_disp_vec(:))
                    IF ( ( (tOH_norm .LT. OHvec_mat(6,o,s)) .OR.&
                    (OHvec_mat(6,o,s) .EQ. 0.0) ) .AND.&
                    (atm_mat(1,j,s) .NE. OHvec_mat(5,i,s)) ) THEN
                        OHvec_mat(10,o,s) = OHvec_mat(5,o,s)
                        OHvec_mat(11,o,s) = OHvec_mat(6,o,s)
                        OHvec_mat(5,o,s) = atm_mat(1,j,s)
                        OHvec_mat(6,o,s) = tOH_norm
                        DO k=1,3
                            OHvec_mat(k+11,o,s) = OHvec_mat(k+6,o,s)
                            OHvec_mat(k+6,o,s) =  atm_mat(k+2,j,s)
                        END DO

                    ELSE IF ( ( (tOH_norm .LT. OHvec_mat(11,o,s)) .OR.&
                        (OHvec_mat(6,o,s) .NE. 0.0) .AND.&
                        (OHvec_mat(11,o,s) .EQ. 0.0 )) .AND.&
                        (atm_mat(1,j,s) .NE. OHvec_mat(10,o,s) ) ) THEN
                        OHvec_mat(10,o,s) = atm_mat(1,j,s)
                        OHvec_mat(11,o,s) = tOH_norm
                        DO k=1,3
                            OHvec_mat(k+11,o,s) = atm_mat(k+2,j,s)
                        END DO

                    END IF
                END IF
            END DO
        END IF
    END DO
    DO i=1,nb_o*3
        DO k=1,3
            OHvec_mat(k+39,i,s) = OHvec_mat(k+11,i,s) - OHvec_mat(k+6,i,s)
            OHvec_mat(k+39,i,s) = OHvec_mat(k+39,i,s) - box(k) * ANINT(OHvec_mat(k+39,i,s)/box(k))
            OHvec_mat(k+14,i,s) = (OHvec_mat(k+11,i,s) + OHvec_mat(k+6,i,s)) / 2.0
            OHvec_mat(k+17,i,s) = OHvec_mat(k+14,i,s) - OHvec_mat(k+1,i,s)
            OHvec_mat(k+17,i,s) = OHvec_mat(k+17,i,s) - box(k) * ANINT(OHvec_mat(k+17,i,s)/box(k))
        END DO
    END DO
    nb_max_OHvec(s) = COUNT(OHvec_mat(1,:,s) .NE. 0, DIM=1)
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "OH groups:"&
    ,finish-start,"seconds elapsed"


! D ----------------------------------------------- Proximity between functionnal groups and any OH group
start = OMP_get_wtime()

!nb_step, nb_atm, rXOH_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_o)&
!$OMP PRIVATE(s,i,j,k,tXOHvec_disp_vec,tXOHvec_norm)
DO s = 1, nb_step
    D1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE D1
        END IF
        D2:DO j = 1, nb_atm
            IF (atm_mat(9,j,s) .EQ. -1) THEN
                CYCLE D2
            END IF
            DO k = 1, 3
                tXOHvec_disp_vec(k) = atm_mat(k+2,j,s) - OHvec_mat(k+1,i,s)
                tXOHvec_disp_vec(k) = tXOHvec_disp_vec(k) - box(k) * ANINT(tXOHvec_disp_vec(k)/box(k))
            END DO
            tXOHvec_norm = NORM2(tXOHvec_disp_vec)
            IF( (tXOHvec_norm .LE. rXOH_cut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(38,i,s) = 1
                    OHvec_mat(39,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 10) THEN ! OE
                    OHvec_mat(21,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 11) THEN ! OH
                    OHvec_mat(22,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 12) THEN ! OA
                    OHvec_mat(23,i,s) = 1
                END IF
            ELSE IF( (tXOHvec_norm .LE. rCOH_cut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(39,i,s) = 1
                END IF
            END IF
        END DO D2
    END DO D1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity FG and OH group:"&
    ,finish-start,"seconds elapsed"

! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
IF (file_surf .NE. "0") THEN
    start = OMP_get_wtime()

    !nb_step, nb_atm, always shared.
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_o,surf_mat,nb_surf)&
    !$OMP PRIVATE(s,i,j,k,tSOHvec_disp_vec,tSOHvec_norm,tSS_disp_vec,tss_norm)&
    !$OMP PRIVATE(tPS_uvec_go,tPS_uvec_air,tPS_OOHvec_pos_vec_go,tPS_OOHvec_pos_vec_air)&
    !$OMP PRIVATE(tPS_HHvec_pos_vec)&
    !$OMP PRIVATE(tOHvec_disp_vec,tOHvec_disp_uvec,thhvec_disp_vec)
    DO s = 1, nb_step
        F2:DO i = 1, nb_o*3
            IF (OHvec_mat(1,i,s) .EQ. 0) THEN
                CYCLE F2
            END IF
            DO j = 1, nb_surf(s)
                IF (surf_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        tSOHvec_disp_vec(k) = surf_mat(k,j,s) - OHvec_mat(k+1,i,s)
                        tSOHvec_disp_vec(k) = tSOHvec_disp_vec(k) - box(k) * ANINT(tSOHvec_disp_vec(k)/box(k))
                    END DO
                    tSOHvec_norm = NORM2(tSOHvec_disp_vec)
                    IF ( (tSOHvec_norm .LT. OHvec_mat(24,i,s)) .OR. (OHvec_mat(24,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(24,i,s) = tSOHvec_norm
                        IF (OHvec_mat(4,i,s) .LT. surf_mat(3,j,s)) THEN
                            OHvec_mat(25,i,s) = -1
                        ELSE
                            OHvec_mat(25,i,s) = 1
                        END IF
                        OHvec_mat(26,i,s) = surf_mat(3,j,s)
                        OHvec_mat(34,i,s) = surf_mat(5,j,s)
                    END IF
                ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        tSOHvec_disp_vec(k) = surf_mat(k,j,s) - OHvec_mat(k+1,i,s)
                        tSOHvec_disp_vec(k) = tSOHvec_disp_vec(k) - box(k) * ANINT(tSOHvec_disp_vec(k)/box(k))
                    END DO
                    tSOHvec_norm = NORM2(tSOHvec_disp_vec)
                    IF ( (tSOHvec_norm .LT. OHvec_mat(27,i,s)) .OR. (OHvec_mat(27,i,s) .EQ. 0.0) ) THEN
                        OHvec_mat(27,i,s) = tSOHvec_norm
                        IF (OHvec_mat(4,i,s) .GT. surf_mat(3,j,s)) THEN
                            OHvec_mat(28,i,s) = -1
                        ELSE
                            OHvec_mat(28,i,s) = 1
                        END IF
                        OHvec_mat(29,i,s) = surf_mat(3,j,s)
                        OHvec_mat(35,i,s) = surf_mat(5,j,s)
                    END IF
                END IF
            END DO

            IF ((surf_mat(32,INT(OHvec_mat(34,i,s)),s) .NE. 1) .OR.&
             (surf_mat(33,INT(OHvec_mat(35,i,s)),s) .NE. 1)) THEN

             F3:DO j = 1, nb_surf(s)

                    IF ( (surf_mat(5,j,s) .EQ. surf_mat(5,INT(OHvec_mat(34,i,s)),s)) .OR.&
                    (surf_mat(5,j,s) .EQ. surf_mat(5,INT(OHvec_mat(35,i,s)),s)) ) THEN
                        CYCLE F3
                    END IF

                    IF (surf_mat(4,j,s) .EQ. 1) THEN
                        DO k = 1, 3
                            tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(OHvec_mat(34,i,s)),s)
                            tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                        END DO
                        tSS_norm = NORM2(tSS_disp_vec)

                        IF ( ( (tSS_norm .LT. surf_mat(7,INT(OHvec_mat(34,i,s)),s)) .OR.&
                        (surf_mat(7,INT(OHvec_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (surf_mat(5,j,s) .NE. surf_mat(6,INT(OHvec_mat(34,i,s)),s)) ) THEN
                            surf_mat(11,INT(OHvec_mat(34,i,s)),s) = surf_mat(6,INT(OHvec_mat(34,i,s)),s)
                            surf_mat(12,INT(OHvec_mat(34,i,s)),s) = surf_mat(7,INT(OHvec_mat(34,i,s)),s)
                            surf_mat(6,INT(OHvec_mat(34,i,s)),s) = surf_mat(5,j,s)
                            surf_mat(7,INT(OHvec_mat(34,i,s)),s) = tSS_norm
                            DO k=1,3
                                surf_mat(k+12,INT(OHvec_mat(34,i,s)),s) = surf_mat(k+7,INT(OHvec_mat(34,i,s)),s) !13/14/15
                                surf_mat(k+7,INT(OHvec_mat(34,i,s)),s) = tSS_disp_vec(k) / tSS_norm !8,9,10
                            END DO
                        ELSE IF ( ( (tSS_norm .LT. surf_mat(12,INT(OHvec_mat(34,i,s)),s)) .OR.&
                        ( (surf_mat(7,INT(OHvec_mat(34,i,s)),s) .NE. 0.0) .AND.&
                        (surf_mat(12,INT(OHvec_mat(34,i,s)),s) .EQ. 0.0) ) ) .AND.&
                        (surf_mat(5,j,s) .NE. surf_mat(6,INT(OHvec_mat(34,i,s)),s)) ) THEN
                            surf_mat(11,INT(OHvec_mat(34,i,s)),s) = surf_mat(5,j,s)
                            surf_mat(12,INT(OHvec_mat(34,i,s)),s) = tSS_norm
                            DO k=1,3
                                surf_mat(k+12,INT(OHvec_mat(34,i,s)),s) = tSS_disp_vec(k) / tSS_norm !13/14/15
                            END DO
                        END IF

                    ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                        DO k = 1, 3
                            tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(OHvec_mat(35,i,s)),s)
                            tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                        END DO
                        tSS_norm = NORM2(tSS_disp_vec)

                        IF ( ( (tSS_norm .LT. surf_mat(17,INT(OHvec_mat(35,i,s)),s)) .OR.&
                        (surf_mat(17,INT(OHvec_mat(35,i,s)),s) .EQ. 0.0 ) ) .AND.&
                        (surf_mat(5,j,s) .NE. surf_mat(16,INT(OHvec_mat(35,i,s)),s)) ) THEN
                            surf_mat(21,INT(OHvec_mat(35,i,s)),s) = surf_mat(16,INT(OHvec_mat(35,i,s)),s)
                            surf_mat(22,INT(OHvec_mat(35,i,s)),s) = surf_mat(17,INT(OHvec_mat(35,i,s)),s)
                            surf_mat(16,INT(OHvec_mat(35,i,s)),s) = surf_mat(5,j,s)
                            surf_mat(17,INT(OHvec_mat(35,i,s)),s) = tSS_norm
                            DO k=1,3
                                surf_mat(k+22,INT(OHvec_mat(35,i,s)),s) = surf_mat(k+17,INT(OHvec_mat(35,i,s)),s) !23/24/25
                                surf_mat(k+17,INT(OHvec_mat(35,i,s)),s) = tSS_disp_vec(k) / tSS_norm !18,19,20
                            END DO
                        ELSE IF ( ( (tSS_norm .LT. surf_mat(22,INT(OHvec_mat(35,i,s)),s)) .OR.&
                        ( (surf_mat(17,INT(OHvec_mat(35,i,s)),s) .NE. 0.0) .AND.&
                        (surf_mat(22,INT(OHvec_mat(35,i,s)),s) .EQ. 0.0) ) ) .AND.&
                        (surf_mat(5,j,s) .NE. surf_mat(16,INT(OHvec_mat(35,i,s)),s)) ) THEN
                            surf_mat(21,INT(OHvec_mat(35,i,s)),s) = surf_mat(5,j,s)
                            surf_mat(22,INT(OHvec_mat(35,i,s)),s) = tSS_norm
                            DO k=1,3
                                surf_mat(k+22,INT(OHvec_mat(35,i,s)),s) = tSS_disp_vec(k) / tSS_norm !23/24/25
                            END DO
                        END IF

                    END IF
                END DO F3

                surf_mat(26,INT(OHvec_mat(34,i,s)),s) = surf_mat(9,INT(OHvec_mat(34,i,s)),s) *&
                surf_mat(15,INT(OHvec_mat(34,i,s)),s) - surf_mat(10,INT(OHvec_mat(34,i,s)),s) *&
                 surf_mat(14,INT(OHvec_mat(34,i,s)),s)
                surf_mat(27,INT(OHvec_mat(34,i,s)),s) = surf_mat(10,INT(OHvec_mat(34,i,s)),s) *&
                surf_mat(13,INT(OHvec_mat(34,i,s)),s) - surf_mat(8,INT(OHvec_mat(34,i,s)),s) *&
                 surf_mat(15,INT(OHvec_mat(34,i,s)),s)
                surf_mat(28,INT(OHvec_mat(34,i,s)),s) = surf_mat(8,INT(OHvec_mat(34,i,s)),s) *&
                surf_mat(14,INT(OHvec_mat(34,i,s)),s) - surf_mat(9,INT(OHvec_mat(34,i,s)),s) *&
                 surf_mat(13,INT(OHvec_mat(34,i,s)),s)

                surf_mat(29,INT(OHvec_mat(35,i,s)),s) = surf_mat(19,INT(OHvec_mat(35,i,s)),s) *&
                surf_mat(25,INT(OHvec_mat(35,i,s)),s) - surf_mat(20,INT(OHvec_mat(35,i,s)),s) *&
                 surf_mat(24,INT(OHvec_mat(35,i,s)),s)
                surf_mat(30,INT(OHvec_mat(35,i,s)),s) = surf_mat(20,INT(OHvec_mat(35,i,s)),s) *&
                surf_mat(23,INT(OHvec_mat(35,i,s)),s) - surf_mat(28,INT(OHvec_mat(35,i,s)),s) *&
                 surf_mat(25,INT(OHvec_mat(35,i,s)),s)
                surf_mat(31,INT(OHvec_mat(35,i,s)),s) = surf_mat(28,INT(OHvec_mat(35,i,s)),s) *&
                surf_mat(24,INT(OHvec_mat(35,i,s)),s) - surf_mat(29,INT(OHvec_mat(35,i,s)),s) *&
                 surf_mat(23,INT(OHvec_mat(35,i,s)),s)

                surf_mat(32,INT(OHvec_mat(34,i,s)),s) = 1
                surf_mat(33,INT(OHvec_mat(35,i,s)),s) = 1

            END IF

            DO k=1,3
                tPS_uvec_go(k) = surf_mat(k+25,INT(OHvec_mat(34,i,s)),s)
                tPS_uvec_air(k) = surf_mat(k+28,INT(OHvec_mat(35,i,s)),s)
                tPS_OOHvec_pos_vec_go(k) = OHvec_mat(k+1,i,s) + surf_mat(k+25,INT(OHvec_mat(34,i,s)),s)
                tPS_OOHvec_pos_vec_go(k) = tPS_OOHvec_pos_vec_go(k) - box(k) * ANINT(tPS_OOHvec_pos_vec_go(k)/box(k))
                tPS_OOHvec_pos_vec_air(k) = OHvec_mat(k+1,i,s) + surf_mat(k+28,INT(OHvec_mat(35,i,s)),s)
                tPS_OOHvec_pos_vec_air(k) = tPS_OOHvec_pos_vec_air(k) - box(k) * ANINT(tPS_OOHvec_pos_vec_air(k)/box(k))
                tPS_HHvec_pos_vec(k) = OHvec_mat(k+1,i,s) + OHvec_mat(k+39,i,s)
                tPS_HHvec_pos_vec(k) = tPS_HHvec_pos_vec(k) - box(k) * ANINT(tPS_HHvec_pos_vec(k)/box(k))
                tOHvec_disp_vec(k) = OHvec_mat(k+17,i,s)
                tHHvec_disp_vec(k) = OHvec_mat(k+39,i,s)
            END DO
            tOHvec_disp_uvec(:) = tOHvec_disp_vec(:) / NORM2(tOHvec_disp_vec(:))
            IF (NORM2(tPS_OOHvec_pos_vec_go(:)) .LT. OHvec_mat(24,i,s)) THEN
                tPS_uvec_go(:) = -1.0 * tPS_uvec_go(:)
            END IF
            IF (NORM2(tPS_OOHvec_pos_vec_air(:)) .LT. OHvec_mat(27,i,s)) THEN
                tPS_uvec_air(:) = -1.0 * tPS_uvec_air(:)
            END IF

            IF (NORM2(tPS_HHvec_pos_vec(:)) .GT. OHvec_mat(24,i,s)) THEN
                tHHvec_disp_vec(:) = -1.0 *  tHHvec_disp_vec(:)
            END IF

            OHvec_mat(36,i,s) = ACOS(DOT_PRODUCT(tPS_uvec_go(:),  tOHvec_disp_vec(:)))
            OHvec_mat(44,i,s) = ACOS(DOT_PRODUCT(tPS_uvec_go(:),  tHHvec_disp_vec(:)))

            IF (NORM2(tPS_HHvec_pos_vec(:)) .GT. OHvec_mat(27,i,s)) THEN
                tHHvec_disp_vec(:) = -1.0 *  tHHvec_disp_vec(:)
            END IF
        
            OHvec_mat(37,i,s) = ACOS(DOT_PRODUCT(tPS_uvec_air(:), tOHvec_disp_uvec(:)))
            OHvec_mat(45,i,s) = ACOS(DOT_PRODUCT(tPS_uvec_air(:),  tHHvec_disp_vec(:)))

        END DO F2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "Proximity IS and OH groups:"&
        ,finish-start,"seconds elapsed"
END IF

! X ----------------------------------------------- Write OH BOND
IF (hbond_output .EQ. 1) THEN
    start = OMP_get_wtime()

    OPEN(UNIT=32, FILE = suffix//"_HOH_hbonds.txt")
    WRITE(32,'(A10,A10,A24,A24,A24,A24,A24,A24)')&
        "O_id", "C9", "dist_IS_go", "dist_IS_air"&
    , "Angle OH/NIS_go", "Angle OH/NIS_air"&
    , "Angle HH/NIS_go", "Angle HH/NIS_air"
    DO s = 1, nb_step
        DO i = 1, nb_max_OHvec(s)
            WRITE(32,'(I10,I10,E24.14,E24.14,E24.14,E24.14,E24.14,E24.14)')&
                INT(OHvec_mat(1,i,s)),INT(OHvec_mat(39,i,s))&
            , (OHvec_mat(24,i,s)*OHvec_mat(25,i,s)), (OHvec_mat(27,i,s)*OHvec_mat(28,i,s)), OHvec_mat(36,i,s), OHvec_mat(37,i,s)&
            , OHvec_mat(44,i,s), OHvec_mat(45,i,s)
        END DO
    END DO
    CLOSE(UNIT=32)

    finish = OMP_get_wtime()
    PRINT'(A40,F14.2,A20)', "O/OH Hbonds output:"&
        ,finish-start,"seconds elapsed"
END IF
! ----------------------------------------------- Deallocate and exit
DEALLOCATE(OHvec_mat,atm_mat)
END PROGRAM vvcf
