PROGRAM vvcf
USE OMP_LIB

IMPLICIT NONE

! ----------------------------------------------- Set Double precision ----------------------------------------------- !
INTEGER, PARAMETER              :: dp=KIND(0.0d0)
INTEGER                         :: iostatus

! ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

! ----------------------------------------------- Filenames
CHARACTER(LEN=64)               :: file_xyz, file_surf

! ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), WAT_mat(:,:,:), surf_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
INTEGER                         :: nb_line, nb_max_pt
INTEGER, ALLOCATABLE            :: nb_surf(:)

! ----------------------------------------------- Temp variables
REAL(dp)                        :: tOH_disp_vec(3), tOH_norm
REAL(dp)                        :: tXWAT_disp_vec(3), tXWAT_disp_norm
REAL(dp)                        :: tSWAT_disp_vec(3), tSWAT_norm
REAL(dp)                        :: tSS_disp_vec(3), tSS_norm, tSS1_disp_vec(3), tSS1_norm
REAL(dp)                        :: tPSuvec_go(3), tPSuvec_air(3)
REAL(dp)                        :: tWDvec(3), tHHvec(3),tWDuvec(3), tHHuvec(3)
REAL(dp)                        :: tO_pos_iPS_go(3), tO_pos_iPS_air(3), tS_pos_PS_go(3), tS_pos_PS_air(3)
REAL(dp)                        :: tOS_disp_oPS_go(3), tOS_disp_oPS_air(3)
CHARACTER(LEN=2)                :: dummy_char

! ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h
INTEGER,ALLOCATABLE             :: nb_max_WAT(:)

! ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o
INTEGER                         :: CAC

! ----------------------------------------------- Parameters
REAL(dp), PARAMETER             :: bohr_to_angstrom=0.529177249_dp
REAL(dp), PARAMETER             :: aut_to_fs=0.0241888432658569977_dp
REAL(dp), PARAMETER             :: pi=4.0_dp*DATAN(1.0_dp)
REAL(dp)                        :: box(3)

! ----------------------------------------------- SYSTEM DEPENDANT
INTEGER, PARAMETER              :: nb_atm=1039, nb_step=1000
CHARACTER(LEN=2)                :: suffix="00"
REAL(dp), PARAMETER             :: xlo=0.0_dp,xhi=21.8489966560_dp
REAL(dp), PARAMETER             :: ylo=0.0_dp,yhi=21.2373561788_dp
REAL(dp), PARAMETER             :: zlo=0.0_dp,zhi=70.0_dp
REAL(dp), PARAMETER             :: rOH_cut=1.3, rOHbond_cut=2.5, rXWAT_cut=3.5
REAL(dp), PARAMETER             :: rCWAT_cut=9.0 ! Close to Carbon (interface)

! ----------------------------------------------- Allocate function for reading files
! DEFINE AS: atm_id, atm_nb, atm_x, atm_y, atm_z, vel_x, vel_y, vel_z, nb_OH
ALLOCATE(atm_mat(6,nb_atm,nb_step))
atm_mat(:,:,:) = 0.0_dp

! ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*,"NO ARGUMENT"
    STOP
ELSE IF (CAC .EQ. 1) THEN
    PRINT*, "Need positions and IS files"
    STOP
END if

CALL GET_COMMAND_ARGUMENT(1,file_xyz)
file_xyz=TRIM(file_xyz)
CALL GET_COMMAND_ARGUMENT(2,file_surf)
file_surf=TRIM(file_surf)

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

nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Positions:",finish-start,"seconds elapsed"

! ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
start = OMP_get_wtime()

OPEN(UNIT=21,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
nb_line=0
DO
    READ(21,*,IOSTAT=iostatus)
    IF (iostatus .NE. 0) THEN
        EXIT
    ELSE
        nb_line=nb_line+1
    END IF
END DO
REWIND(21)
nb_max_pt=CEILING(1.0*nb_line/nb_step)*2

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS grid:",finish-start,"seconds elapsed"

! A ----------------------------------------------- Read surface
start = OMP_get_wtime()

ALLOCATE(surf_mat(19,nb_max_pt,nb_step))
ALLOCATE(nb_surf(nb_step))
surf_mat(:,:,:) = 0.0_dp
nb_surf(:) = 0

OPEN(UNIT=21,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
DO s = 1, nb_step
    READ(21, *) nb_surf(s)
    READ(21, *)
    j = 0
    DO i=1,nb_surf(s)
        READ(21, *) dummy_char, surf_mat(1,i,s), surf_mat(2,i,s), surf_mat(3,i,s)
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
CLOSE(UNIT=21)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "IS:",finish-start,"seconds elapsed"

! B -------------WAT---------------------------------- Water
start = OMP_get_wtime()
ALLOCATE(WAT_mat(45,nb_o*3,nb_step))
ALLOCATE(nb_max_WAT(nb_step))

nb_max_WAT(:) = 0.0
WAT_mat(:,:,:) = 0.0_dp

!nb_step, nb_atm, rWAT_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,WAT_mat,nb_max_WAT,nb_o)&
!$OMP PRIVATE(s,i,j,k,o,tOH_disp_vec,tOH_norm)
DO s = 1, nb_step
    o = 0
    DO i=1,nb_atm
        IF (atm_mat(3,i,s) .EQ. 13) THEN
            o = o + 1
            WAT_mat(1,o,s) = atm_mat(1,i,s)
            DO k=1,3
                WAT_mat(k+1,o,s) = atm_mat(k+3,i,s)
            END DO
            DO j=1,nb_atm
                IF (atm_mat(3,j,s) .EQ. 23) THEN
                    DO k = 1, 3
                        tOH_disp_vec(k) = atm_mat(k+3,j,s) - atm_mat(k+3,i,s)
                        tOH_disp_vec(k) = tOH_disp_vec(k) - box(k) * ANINT(tOH_disp_vec(k)/box(k))
                    END DO
                    tOH_norm = NORM2(tOH_disp_vec(:))
                    IF ( ( (tOH_norm .LT. WAT_mat(6,o,s)) .OR.&
                    (WAT_mat(6,o,s) .EQ. 0.0) ) .AND.&
                    (atm_mat(1,j,s) .NE. WAT_mat(5,i,s)) ) THEN
                        WAT_mat(10,o,s) = WAT_mat(5,o,s)
                        WAT_mat(11,o,s) = WAT_mat(6,o,s)
                        WAT_mat(5,o,s) = atm_mat(1,j,s)
                        WAT_mat(6,o,s) = tOH_norm
                        DO k=1,3
                            WAT_mat(k+11,o,s) = WAT_mat(k+6,o,s)
                            WAT_mat(k+6,o,s) =  atm_mat(k+3,j,s)
                        END DO

                    ELSE IF ( ( (tOH_norm .LT. WAT_mat(11,o,s)) .OR.&
                        (WAT_mat(6,o,s) .NE. 0.0) .AND.&
                        (WAT_mat(11,o,s) .EQ. 0.0 )) .AND.&
                        (atm_mat(1,j,s) .NE. WAT_mat(10,o,s) ) ) THEN
                        WAT_mat(10,o,s) = atm_mat(1,j,s)
                        WAT_mat(11,o,s) = tOH_norm
                        DO k=1,3
                            WAT_mat(k+11,o,s) = atm_mat(k+3,j,s)
                        END DO

                    END IF
                END IF
            END DO
        END IF
    END DO
    DO i=1,nb_o*3
        DO k=1,3
            WAT_mat(k+39,i,s) = WAT_mat(k+11,i,s) - WAT_mat(k+6,i,s)
            WAT_mat(k+39,i,s) = WAT_mat(k+39,i,s) - box(k) * ANINT(WAT_mat(k+39,i,s)/box(k)) ! HH Vector
            WAT_mat(k+14,i,s) = WAT_mat(k+6,i,s) + WAT_mat(k+39,i,s)/2.0 ! Center of HH
            WAT_mat(k+17,i,s) = WAT_mat(k+14,i,s) - WAT_mat(k+1,i,s)
            WAT_mat(k+17,i,s) = WAT_mat(k+17,i,s) - box(k) * ANINT(WAT_mat(k+17,i,s)/box(k)) ! (O to center HH) vector (Inverse of of the Dipole Moment)
        END DO
    END DO
    nb_max_WAT(s) = COUNT(WAT_mat(1,:,s) .NE. 0, DIM=1)
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "WAT groups:"&
    ,finish-start,"seconds elapsed"


! D ----------------------------------------------- Proximity between functionnal groups and any WAT groups
start = OMP_get_wtime()

!nb_step, nb_atm, rXOH_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,WAT_mat,nb_o)&
!$OMP PRIVATE(s,i,j,k,tXWAT_disp_vec,tXWAT_disp_norm)
DO s = 1, nb_step
    D1:DO i = 1, nb_o*3
        IF (WAT_mat(1,i,s) .EQ. 0) THEN
            CYCLE D1
        END IF
        D2:DO j = 1, nb_atm
            IF (atm_mat(3,j,s) .EQ. -1) THEN
                CYCLE D2
            END IF
            DO k = 1, 3
                tXWAT_disp_vec(k) = atm_mat(k+3,j,s) - WAT_mat(k+1,i,s)
                tXWAT_disp_vec(k) = tXWAT_disp_vec(k) - box(k) * ANINT(tXWAT_disp_vec(k)/box(k))
            END DO
            tXWAT_disp_norm = NORM2(tXWAT_disp_vec)
            IF( (tXWAT_disp_norm .LE. rXWAT_cut) .AND. (atm_mat(1,j,s) .NE. WAT_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    WAT_mat(38,i,s) = 1
                    WAT_mat(39,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 10) THEN ! OE
                    WAT_mat(21,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 11) THEN ! OH
                    WAT_mat(22,i,s) = 1
                ELSE IF (atm_mat(3,j,s) .EQ. 12) THEN ! OA
                    WAT_mat(23,i,s) = 1
                END IF
            ELSE IF( (tXWAT_disp_norm .LE. rCWAT_cut) .AND. (atm_mat(1,j,s) .NE. WAT_mat(1,i,s))) THEN
                IF (atm_mat(3,j,s) .EQ. 1.) THEN ! C
                    WAT_mat(39,i,s) = 1
                END IF
            END IF
        END DO D2
    END DO D1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity WAT groups and FG groups:"&
    ,finish-start,"seconds elapsed"

! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
start = OMP_get_wtime()

!nb_step, nb_atm, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,WAT_mat,nb_o,surf_mat,nb_surf)&
!$OMP PRIVATE(s,i,j,k,tSWAT_disp_vec,tSWAT_norm,tSS_disp_vec,tss_norm,tSS1_disp_vec,tSS1_norm)&
!$OMP PRIVATE(tPSuvec_go,tPSuvec_air,tWDvec,tHHvec,tWDuvec,tHHuvec)&
!$OMP PRIVATE(tO_pos_iPS_go,tO_pos_iPS_air,tS_pos_PS_go,tS_pos_PS_air)&
!$OMP PRIVATE(tOS_disp_oPS_go,tOS_disp_oPS_air)
DO s = 1, nb_step
    F2:DO i = 1, nb_o*3
        IF (WAT_mat(1,i,s) .EQ. 0) THEN
            CYCLE F2
        END IF
        DO j = 1, nb_surf(s)
            IF (surf_mat(4,j,s) .EQ. 1) THEN
                DO k = 1, 3
                    tSWAT_disp_vec(k) = surf_mat(k,j,s) - WAT_mat(k+1,i,s)
                    tSWAT_disp_vec(k) = tSWAT_disp_vec(k) - box(k) * ANINT(tSWAT_disp_vec(k)/box(k))
                END DO
                tSWAT_norm = NORM2(tSWAT_disp_vec)
                IF ( (tSWAT_norm .LT. WAT_mat(24,i,s)) .OR. (WAT_mat(24,i,s) .EQ. 0.0) ) THEN
                    WAT_mat(24,i,s) = tSWAT_norm
                    IF (WAT_mat(4,i,s) .LT. surf_mat(3,j,s)) THEN
                        WAT_mat(25,i,s) = -1
                    ELSE
                        WAT_mat(25,i,s) = 1
                    END IF
                    WAT_mat(26,i,s) = surf_mat(3,j,s)
                    WAT_mat(34,i,s) = surf_mat(5,j,s)
                END IF
            ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                DO k = 1, 3
                    tSWAT_disp_vec(k) = surf_mat(k,j,s) - WAT_mat(k+1,i,s)
                    tSWAT_disp_vec(k) = tSWAT_disp_vec(k) - box(k) * ANINT(tSWAT_disp_vec(k)/box(k))
                END DO
                tSWAT_norm = NORM2(tSWAT_disp_vec)
                IF ( (tSWAT_norm .LT. WAT_mat(27,i,s)) .OR. (WAT_mat(27,i,s) .EQ. 0.0) ) THEN
                    WAT_mat(27,i,s) = tSWAT_norm
                    IF (WAT_mat(4,i,s) .GT. surf_mat(3,j,s)) THEN
                        WAT_mat(28,i,s) = -1
                    ELSE
                        WAT_mat(28,i,s) = 1
                    END IF
                    WAT_mat(29,i,s) = surf_mat(3,j,s)
                    WAT_mat(35,i,s) = surf_mat(5,j,s)
                END IF
            END IF
        END DO

        IF ((surf_mat(19,INT(WAT_mat(34,i,s)),s) .NE. 1) .OR.&
            (surf_mat(19,INT(WAT_mat(35,i,s)),s) .NE. 1)) THEN

            F3:DO j = 1, nb_surf(s) ! First one
                IF ( (surf_mat(5,j,s) .EQ. surf_mat(5,INT(WAT_mat(34,i,s)),s)) .OR.&
                (surf_mat(5,j,s) .EQ. surf_mat(5,INT(WAT_mat(35,i,s)),s)) ) THEN
                    CYCLE F3
                END IF
                IF (surf_mat(4,j,s) .EQ. 1) THEN

                    DO k = 1, 3
                        tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(WAT_mat(34,i,s)),s)
                        tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                    END DO
                    tSS_norm = NORM2(tSS_disp_vec)

                    IF ( ( (tSS_norm .LT. surf_mat(7,INT(WAT_mat(34,i,s)),s)) .OR.&
                    (surf_mat(7,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                    (surf_mat(5,j,s) .NE. surf_mat(6,INT(WAT_mat(34,i,s)),s)) ) THEN
                        surf_mat(6,INT(WAT_mat(34,i,s)),s) = surf_mat(5,j,s)
                        surf_mat(7,INT(WAT_mat(34,i,s)),s) = tSS_norm
                        DO k=1,3
                            surf_mat(k+7,INT(WAT_mat(34,i,s)),s) = tSS_disp_vec(k) / tSS_norm !8,9,10
                        END DO
                    END IF

                ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN

                    DO k = 1, 3
                        tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(WAT_mat(35,i,s)),s)
                        tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                    END DO
                    tSS_norm = NORM2(tSS_disp_vec)

                    IF ( ( (tSS_norm .LT. surf_mat(7,INT(WAT_mat(35,i,s)),s)) .OR.&
                    (surf_mat(7,INT(WAT_mat(35,i,s)),s) .EQ. 0.0 ) ) .AND.&
                    (surf_mat(5,j,s) .NE. surf_mat(6,INT(WAT_mat(35,i,s)),s)) ) THEN
                        surf_mat(6,INT(WAT_mat(35,i,s)),s) = surf_mat(5,j,s)
                        surf_mat(7,INT(WAT_mat(35,i,s)),s) = tSS_norm
                        DO k=1,3
                            surf_mat(k+7,INT(WAT_mat(35,i,s)),s) = tSS_disp_vec(k) / tSS_norm !18,19,20
                        END DO
                    END IF

                END IF
            END DO F3

         F4:DO j = 1, nb_surf(s)

                IF ( (surf_mat(5,j,s) .EQ. surf_mat(5,INT(WAT_mat(34,i,s)),s)) .OR.&
                (surf_mat(5,j,s) .EQ. surf_mat(5,INT(WAT_mat(35,i,s)),s)) .OR. &
                (surf_mat(5,j,s) .EQ. surf_mat(6,INT(WAT_mat(34,i,s)),s)) .OR. &
                (surf_mat(5,j,s) .EQ. surf_mat(6,INT(WAT_mat(35,i,s)),s)) ) THEN
                    CYCLE F4
                END IF

                IF (surf_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(WAT_mat(34,i,s)),s)
                        tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                        tSS1_disp_vec(k) = surf_mat(k+7,INT(WAT_mat(34,i,s)),s)
                        tSS1_disp_vec(k) = tSS1_disp_vec(k) - box(k) * ANINT(tSS1_disp_vec(k)/box(k))
                    END DO

                    tSS_norm = NORM2(tSS_disp_vec)
                    tSS1_norm = NORM2(tSS1_disp_vec)

                    IF ( (ACOS(DOT_PRODUCT(tSS_disp_vec(:)/tSS_norm,tSS1_disp_vec(:)/tSS1_norm)) .LT. 0.50).OR.&
                    (ACOS(DOT_PRODUCT(tSS_disp_vec(:)/tSS_norm,tSS1_disp_vec(:)/tSS1_norm)) .GT. pi-0.50) ) THEN
                        CYCLE F4
                    END IF

                    IF ( ( (tSS_norm .LT. surf_mat(12,INT(WAT_mat(34,i,s)),s)) .OR.&
                    (surf_mat(12,INT(WAT_mat(34,i,s)),s) .EQ. 0.0 ) ) .AND.&
                    (surf_mat(5,j,s) .NE. surf_mat(11,INT(WAT_mat(34,i,s)),s)) ) THEN

                        surf_mat(11,INT(WAT_mat(34,i,s)),s) = surf_mat(5,j,s)
                        surf_mat(12,INT(WAT_mat(34,i,s)),s) = tSS_norm
                        DO k=1,3
                            surf_mat(k+12,INT(WAT_mat(34,i,s)),s) = tSS_disp_vec(k) / tSS_norm !8,9,10
                        END DO

                    END IF

                ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        tSS_disp_vec(k) = surf_mat(k,j,s) - surf_mat(k,INT(WAT_mat(35,i,s)),s)
                        tSS_disp_vec(k) = tSS_disp_vec(k) - box(k) * ANINT(tSS_disp_vec(k)/box(k))
                        tSS1_disp_vec(k) = surf_mat(k+7,INT(WAT_mat(35,i,s)),s)
                        tSS1_disp_vec(k) = tSS1_disp_vec(k) - box(k) * ANINT(tSS1_disp_vec(k)/box(k))
                    END DO

                    tSS_norm = NORM2(tSS_disp_vec)
                    tSS1_norm = NORM2(tSS1_disp_vec)

                    IF ( (ACOS(DOT_PRODUCT(tSS_disp_vec(:)/tSS_norm,tSS1_disp_vec(:)/tSS1_norm)) .LT. 0.50).OR.&
                    (ACOS(DOT_PRODUCT(tSS_disp_vec(:)/tSS_norm,tSS1_disp_vec(:)/tSS1_norm)) .GT. pi-0.50) ) THEN
                        CYCLE F4
                    END IF

                    IF ( ( (tSS_norm .LT. surf_mat(12,INT(WAT_mat(35,i,s)),s)) .OR.&
                    (surf_mat(12,INT(WAT_mat(35,i,s)),s) .EQ. 0.0 ) ) .AND.&
                    (surf_mat(5,j,s) .NE. surf_mat(11,INT(WAT_mat(35,i,s)),s)) ) THEN

                        surf_mat(11,INT(WAT_mat(35,i,s)),s) = surf_mat(5,j,s)
                        surf_mat(12,INT(WAT_mat(35,i,s)),s) = tSS_norm
                        DO k=1,3
                            surf_mat(k+12,INT(WAT_mat(35,i,s)),s) = tSS_disp_vec(k) / tSS_norm !18,19,20
                        END DO

                    END IF
            END IF

            END DO F4

            surf_mat(16,INT(WAT_mat(34,i,s)),s) =&
                surf_mat(9,INT(WAT_mat(34,i,s)),s) * surf_mat(15,INT(WAT_mat(34,i,s)),s) -&
                surf_mat(10,INT(WAT_mat(34,i,s)),s) * surf_mat(14,INT(WAT_mat(34,i,s)),s)
            surf_mat(17,INT(WAT_mat(34,i,s)),s) =&
                surf_mat(10,INT(WAT_mat(34,i,s)),s) * surf_mat(13,INT(WAT_mat(34,i,s)),s) -&
                surf_mat(8,INT(WAT_mat(34,i,s)),s) * surf_mat(15,INT(WAT_mat(34,i,s)),s)
            surf_mat(18,INT(WAT_mat(34,i,s)),s) =&
                surf_mat(8,INT(WAT_mat(34,i,s)),s) * surf_mat(14,INT(WAT_mat(34,i,s)),s) -&
                surf_mat(9,INT(WAT_mat(34,i,s)),s) * surf_mat(13,INT(WAT_mat(34,i,s)),s)

            surf_mat(16,INT(WAT_mat(35,i,s)),s) =&
                surf_mat(9,INT(WAT_mat(35,i,s)),s) * surf_mat(15,INT(WAT_mat(35,i,s)),s) -&
                surf_mat(10,INT(WAT_mat(35,i,s)),s) * surf_mat(14,INT(WAT_mat(35,i,s)),s)
            surf_mat(17,INT(WAT_mat(35,i,s)),s) =&
                surf_mat(10,INT(WAT_mat(35,i,s)),s) * surf_mat(13,INT(WAT_mat(35,i,s)),s) -&
                surf_mat(8,INT(WAT_mat(35,i,s)),s) * surf_mat(15,INT(WAT_mat(35,i,s)),s)
            surf_mat(18,INT(WAT_mat(35,i,s)),s) =&
                surf_mat(8,INT(WAT_mat(35,i,s)),s) * surf_mat(14,INT(WAT_mat(35,i,s)),s) -&
                surf_mat(9,INT(WAT_mat(35,i,s)),s) * surf_mat(13,INT(WAT_mat(35,i,s)),s)

            surf_mat(19,INT(WAT_mat(34,i,s)),s) = 1
            surf_mat(19,INT(WAT_mat(35,i,s)),s) = 1

        END IF

        DO k=1,3
            tPSuvec_go(k) = surf_mat(k+15,INT(WAT_mat(34,i,s)),s)
            tPSuvec_air(k) = surf_mat(k+15,INT(WAT_mat(35,i,s)),s)
            tWDvec(k) = WAT_mat(k+17,i,s)
            tHHvec(k) = WAT_mat(k+39,i,s)
        END DO
        tWDuvec(:) = tWDvec(:) / NORM2(tWDvec(:))
        tHHuvec(:) = tHHvec(:) / NORM2(tHHvec(:))
        DO k=1,3
            ! To check orientation of the normal surface vector
            tO_pos_iPS_go(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_go(k)
            tO_pos_iPS_go(k) = tO_pos_iPS_go(k) - box(k) * ANINT(tO_pos_iPS_go(k)/box(k))
            tS_pos_PS_go(k) = surf_mat(k,INT(WAT_mat(34,i,s)),s) + 0.01*tPSuvec_go(k)
            tS_pos_PS_go(k) = tS_pos_PS_go(k) - box(k) * ANINT(tS_pos_PS_go(k)/box(k))
            tOS_disp_oPS_go(k) = tS_pos_PS_go(k) - tO_pos_iPS_go(k)
            tOS_disp_oPS_go(k) = tOS_disp_oPS_go(k) - box(k) * ANINT(tOS_disp_oPS_go(k)/box(k))
            tO_pos_iPS_air(k) = WAT_mat(k+1,i,s) - 0.01*tPSuvec_air(k) 
            tO_pos_iPS_air(k) = tO_pos_iPS_air(k) - box(k) * ANINT(tO_pos_iPS_air(k)/box(k))
            tS_pos_PS_air(k) = surf_mat(k,INT(WAT_mat(35,i,s)),s) + 0.01*tPSuvec_air(k) 
            tS_pos_PS_air(k) = tS_pos_PS_air(k) - box(k) * ANINT(tS_pos_PS_air(k)/box(k))
            tOS_disp_oPS_air(k) = tS_pos_PS_air(k) - tO_pos_iPS_air(k)
            tOS_disp_oPS_air(k) = tOS_disp_oPS_air(k) - box(k) * ANINT(tOS_disp_oPS_air(k)/box(k))
        END DO



        IF (NORM2(tOS_disp_oPS_go(:)) .GT. WAT_mat(24,i,s)) THEN
            tPSuvec_go(:) = -1.0 * tPSuvec_go(:)
        END IF
        IF (NORM2(tOS_disp_oPS_air(:)) .GT. WAT_mat(27,i,s)) THEN
            tPSuvec_air(:) = -1.0 * tPSuvec_air(:)
        END IF

        IF (ACOS(DOT_PRODUCT(tPSuvec_go(:), tHHuvec(:))) .LT. pi/2.0) THEN
            tHHuvec(:) = -1.0 * tHHuvec(:)
        END IF

        WAT_mat(36,i,s) = ACOS(DOT_PRODUCT(tPSuvec_go(:), tWDuvec(:)))
        WAT_mat(44,i,s) = ACOS(DOT_PRODUCT(tPSuvec_go(:), tHHuvec(:)))

        IF (ACOS(DOT_PRODUCT(tPSuvec_air(:), tHHuvec(:))) .LT. pi/2.0) THEN
            tHHuvec(:) = -1.0 * tHHuvec(:)
        END IF

        WAT_mat(37,i,s) = ACOS(DOT_PRODUCT(tPSuvec_air(:), tWDuvec(:)))
        WAT_mat(45,i,s) = ACOS(DOT_PRODUCT(tPSuvec_air(:), tHHuvec(:)))

    END DO F2
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity WAT groups and IS:"&
    ,finish-start,"seconds elapsed"

! X ----------------------------------------------- Write OH BOND
start = OMP_get_wtime()

OPEN(UNIT=32, FILE = suffix//"_water-angle.txt")
WRITE(32,'(A10,A10,A20,A20,A20,A20,A20,A20)')&
    "O_id", "C9", "dist_IS_go", "dist_IS_air"&
, "Angle OH/NIS_go", "Angle OH/NIS_air"&
, "Angle HH/NIS_go", "Angle HH/NIS_air"
DO s = 1, nb_step
    DO i = 1, nb_max_WAT(s)
        WRITE(32,'(I10,I10,E20.7,E20.7,E20.7,E20.7,E20.7,E20.7)')&
            INT(WAT_mat(1,i,s)),INT(WAT_mat(39,i,s))&
        , (WAT_mat(24,i,s)*WAT_mat(25,i,s)), (WAT_mat(27,i,s)*WAT_mat(28,i,s)), WAT_mat(36,i,s), WAT_mat(37,i,s)&
        , WAT_mat(44,i,s), WAT_mat(45,i,s)
    END DO
END DO
CLOSE(UNIT=32)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "WAT angles output:"&
    ,finish-start,"seconds elapsed"

! ----------------------------------------------- Deallocate and exit
DEALLOCATE(WAT_mat,atm_mat,surf_mat)
END PROGRAM vvcf