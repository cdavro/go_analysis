PROGRAM vvcf
USE OMP_LIB

IMPLICIT NONE

! ----------------------------------------------- Set Double precision ----------------------------------------------- !
INTEGER, PARAMETER              :: dp=KIND(0.0d0)
INTEGER                         :: iostatus

! ----------------------------------------------- Timings
REAL(dp)                        :: start,finish,start_i,finish_i,avg_timigs
REAL(dp),ALLOCATABLE            :: timings(:)

! ----------------------------------------------- Filenames
CHARACTER(LEN=64)               :: file_xyz, file_vel, file_surf

! ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), OHvec_mat(:,:,:), surf_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
REAL(dp), ALLOCATABLE           :: dens_go(:,:), dens_air(:,:), avg_dens_go(:), avg_dens_air(:)
INTEGER                         :: nb_line, nb_max_pt
INTEGER, ALLOCATABLE            :: nb_surf(:)

! ----------------------------------------------- Temp variables
REAL(dp)                        :: tOH_disp_vec(3), tOH_norm, tOOHvec_disp_vec(3), tOOHvec_norm
REAL(dp)                        :: tXOHvec_disp_vec(3), tXOHvec_norm, tXOvec_disp_vec(3), tXOvec_norm
REAL(dp)                        :: tSOHvec_disp_vec(3), tSOHvec_norm, tSOvec_disp_vec(3), tSOvec_norm
INTEGER                         :: Udonnor_count, Udonnor_count2
CHARACTER(LEN=2)                :: dummy_char
REAL(dp)                        :: r
INTEGER                         :: count_dens_go, count_dens_air

! ----------------------------------------------- VVCF
REAL(dp)                        :: mct, mctb
INTEGER                         :: mcs, mcsb
REAL(dp)                        :: OH_vel_vec(3), OH_disp_vec(3)
REAL(dp), ALLOCATABLE           :: vvcf_xxz(:)
REAL(dp)                        :: tij_vec(3), trij
INTEGER, ALLOCATABLE            :: v(:,:)

! ----------------------------------------------- Count variables
INTEGER                         :: nb_o, nb_h
INTEGER,ALLOCATABLE             :: nb_max_OHvec(:)

! ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o, l, t, u
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
REAL(dp), PARAMETER             :: r_vvcf_cut=2.0 ! Cut for cross-correlation (with trick)
REAL(dp), PARAMETER             :: surf_up=12345, surf_down=67890
! Select 0 0 = 0, 0 1 = 1N, 1 0 = 1S, 1 1 = 2S, 1 2 = 3S, 2 0 = 2D, 2 1 = 3D, 2 2 = 4D, 2 3 = 5D
INTEGER, PARAMETER              :: UDC=2, AC=3
REAL(dp)                        :: box(3)
REAL(dp), PARAMETER             :: dr = 0.5_dp
INTEGER, PARAMETER              :: dens_step = 75
REAL(dp), PARAMETER             :: timestep_fs = 0.5_dp
INTEGER, PARAMETER              :: hbond_output = 1, density_output = 1
INTEGER, PARAMETER              :: vvcf_c = 0
INTEGER, PARAMETER              :: hbonds = 0, hbonds_c = 0
INTEGER, PARAMETER              :: layers = 0, layers_c = 0


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
        ELSE IF (atm_el(i) .EQ. "H") THEN
            atm_mat(2,i,s) = 1
            atm_mat(9,i,s) = -1
        END IF
    END DO
END DO
CLOSE(UNIT=20)

finish = OMP_get_wtime()
PRINT'(A40,F15.5,A20)', "Positions:",finish-start,"seconds elapsed"

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
PRINT'(A40,F15.5,A20)', "Velocities:",finish-start,"seconds elapsed"

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
    PRINT'(A40,F15.5,A20)', "IS grid:",finish-start,"seconds elapsed"
END IF

! A ----------------------------------------------- Read surface
IF (file_surf .NE. "0") THEN

    ALLOCATE(surf_mat(4,nb_max_pt,nb_step))
    ALLOCATE(nb_surf(nb_step))
    surf_mat(:,:,:) = 0.0_dp
    nb_surf(:) = 0

    start = OMP_get_wtime()

    OPEN(UNIT=22,FILE=file_surf,STATUS='old',FORM='formatted',ACTION='READ')
    DO s = 1, nb_step
        READ(22, *) nb_surf(s)
        READ(22, *)
        DO i=1,nb_surf(s)
            READ(22, *) dummy_char, surf_mat(1,i,s), surf_mat(2,i,s), surf_mat(3,i,s)
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
    PRINT'(A40,F15.5,A20)', "IS:",finish-start,"seconds elapsed"
END IF

! ----------------------------------------------- Get number of relevent atom types
nb_o = COUNT(atm_mat(2,:,1) .EQ. 16, DIM=1)
nb_h = COUNT(atm_mat(2,:,1) .EQ. 1, DIM=1)

! B ----------------------------------------------- OH groups and corresponding values
start = OMP_get_wtime()
ALLOCATE(OHvec_mat(33,nb_o*3,nb_step))
ALLOCATE(nb_max_OHvec(nb_step))

nb_max_OHvec(:) = 0.0
OHvec_mat(:,:,:) = 0.0_dp

!nb_step, nb_atm, rOH_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_max_OHvec)&
!$OMP PRIVATE(s,i,j,k,o,tOH_disp_vec,tOH_norm)
DO s = 1, nb_step
    o = 0
    DO i=1,nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            DO j=1,nb_atm
                IF (atm_mat(2,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        tOH_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                        tOH_disp_vec(k) = tOH_disp_vec(k) - box(k) * ANINT(tOH_disp_vec(k)/box(k))
                    END DO
                    tOH_norm = NORM2(tOH_disp_vec)
                    IF(tOH_norm .LT. rOH_cut) THEN
                        o = o + 1
                        OHvec_mat(1,o,s) = atm_mat(1,i,s)
                        OHvec_mat(23,o,s) = atm_mat(9,i,s)
                        OHvec_mat(2,o,s) = atm_mat(1,j,s)

                        DO k = 1, 3
                            OHvec_mat(k+2,o,s) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                            OHvec_mat(k+2,o,s) = OHvec_mat(k+2,o,s) - box(k) * ANINT(OHvec_mat(k+2,o,s)/box(k)) ! Disp
                            OHvec_mat(k+5,o,s) = atm_mat(k+5,j,s) - atm_mat(k+5,i,s) ! Vel
                            OHvec_mat(k+8,o,s) = atm_mat(k+2,i,s) ! O pos
                            OHvec_mat(k+11,o,s) = atm_mat(k+2,j,s) ! H pos
                        END DO

                    END IF
                END IF
            END DO
        END IF
    END DO
    nb_max_OHvec(s) = COUNT(OHvec_mat(1,:,s) .NE. 0, DIM=1)
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F15.5,A20)', "OH groups:"&
    ,finish-start,"seconds elapsed"

! C ----------------------------------------------- OH/O Hbonds
start = OMP_get_wtime()

!nb_step, nb_atm, rOH_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_o)&
!$OMP PRIVATE(s,i,j,k,l,tOOHvec_disp_vec,tOOHvec_norm,Udonnor_count,Udonnor_count2)
DO s = 1, nb_step
 C1:DO i = 1, nb_o*3
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE C1
        END IF
        Udonnor_count = 0
        Udonnor_count2 = 0
        DO j = 1, nb_atm
            IF (atm_mat(2,j,s) .EQ. 16) THEN
                DO k = 1, 3
                    tOOHvec_disp_vec(k) = atm_mat(k+2,j,s) - OHvec_mat(k+11,i,s)
                    tOOHvec_disp_vec(k) = tOOHvec_disp_vec(k) - box(k) * ANINT(tOOHvec_disp_vec(k)/box(k))
                END DO
                tOOHvec_norm = NORM2(tOOHvec_disp_vec)
                IF( (tOOHvec_norm .LE. rOHbond_cut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                    atm_mat(10,j,s) = atm_mat(10,j,s) + 1 ! Acceptor count (O)
                    OHvec_mat(17,i,s) = OHvec_mat(17,i,s) + 1 ! Donnor count (OH)
                    IF (Udonnor_count .EQ. 0) THEN
                        OHvec_mat(18,i,s) = OHvec_mat(18,i,s) + 1 ! Unique donnor count (OH)
                        Udonnor_count = 1
                    END IF
                C21:DO l=1, nb_atm
                        IF (atm_mat(1,l,s) .EQ. OHvec_mat(1,i,s)) THEN
                            atm_mat(11,l,s) = atm_mat(11,l,s) + 1 ! Donnor count (O)
                            IF (Udonnor_count2 .EQ. 0) THEN
                                atm_mat(12,l,s) = atm_mat(12,l,s) + 1 ! Unique donnor count (O)
                                Udonnor_count2 = 1
                            END IF
                            EXIT C21
                        END IF
                    END DO C21
                C22:DO l=1,nb_o*3
                        IF (OHvec_mat(1,l,s) .EQ. 0) THEN
                            CYCLE C22
                        ELSE IF (OHvec_mat(1,l,s) .EQ. atm_mat(1,j,s)) THEN
                            OHvec_mat(16,l,s) = OHvec_mat(16,l,s) + 1 ! Acceptor count (OH)
                        END IF
                    END DO C22
                END IF
            END IF
        END DO
    END DO C1
  C3:DO i = 1, (nb_o*3)
        IF (OHvec_mat(1,i,s) .EQ. 0) THEN
            CYCLE C3
        END IF
        DO j = 1, nb_atm
            IF (OHvec_mat(1,i,s) .EQ. atm_mat(1,j,s)) THEN
                OHvec_mat(31,i,s) = atm_mat(11,j,s)
                OHvec_mat(32,i,s) = atm_mat(12,j,s)
                OHvec_mat(33,i,s) = atm_mat(10,j,s)
            END IF
        END DO
    END DO C3
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F15.5,A20)', "OH/O Hbonds:"&
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
                tXOHvec_disp_vec(k) = atm_mat(k+2,j,s) - OHvec_mat(k+8,i,s)
                tXOHvec_disp_vec(k) = tXOHvec_disp_vec(k) - box(k) * ANINT(tXOHvec_disp_vec(k)/box(k))
            END DO
            tXOHvec_norm = NORM2(tXOHvec_disp_vec)
            IF( (tXOHvec_norm .LE. rXOH_cut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(19,i,s) = 1
                    OHvec_mat(30,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 10) THEN ! OE
                    OHvec_mat(20,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 11) THEN ! OH
                    OHvec_mat(21,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 12) THEN ! OA
                    OHvec_mat(22,i,s) = 1
                END IF
            ELSE IF( (tXOHvec_norm .LE. rCOH_cut) .AND. (atm_mat(1,j,s) .NE. OHvec_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1.) THEN ! C
                    OHvec_mat(30,i,s) = 1
                END IF
            END IF
        END DO D2
    END DO D1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F15.5,A20)', "Proximity FG and OH group:"&
    ,finish-start,"seconds elapsed"

! E ----------------------------------------------- Proximity between functionnal groups and any O atom
start = OMP_get_wtime()

!nb_step, nb_atm, rXO_cut, always shared.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box)&
!$OMP PRIVATE(s,i,j,k,tXOvec_disp_vec,tXOvec_norm)
DO s = 1, nb_step
    E1:DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .NE. 16) THEN
            CYCLE E1
        END IF
        E2:DO j = 1, nb_atm
            IF (atm_mat(9,j,s) .EQ. -1) THEN
                CYCLE E2
            END IF
            DO k = 1, 3
                tXOvec_disp_vec(k) = atm_mat(k+2,j,s) - atm_mat(k+2,i,s)
                tXOvec_disp_vec(k) = tXOvec_disp_vec(k) - box(k) * ANINT(tXOvec_disp_vec(k)/box(k))
            END DO
            tXOvec_norm = NORM2(tXOvec_disp_vec)
            IF ( (tXOvec_norm .LE. rXO_cut) .AND. (atm_mat(1,j,s) .NE. atm_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1) THEN ! C
                    atm_mat(13,i,s) = 1
                    atm_mat(17,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 10) THEN ! OE
                    atm_mat(14,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 11) THEN ! OH
                    atm_mat(15,i,s) = 1
                ELSE IF (atm_mat(9,j,s) .EQ. 12) THEN ! OA
                    atm_mat(16,i,s) = 1
                END IF
            ELSE IF ( (tXOvec_norm .LE. rCO_cut) .AND. (atm_mat(1,j,s) .NE. atm_mat(1,i,s))) THEN
                IF (atm_mat(9,j,s) .EQ. 1) THEN ! C
                    atm_mat(17,i,s) = 1
                END IF
            END IF
        END DO E2
    END DO E1
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F15.5,A20)', "Proximity FG and O atoms:"&
    ,finish-start,"seconds elapsed"


! F ----------------------------------------------- Calculate closest distance between IS and any OH groups
IF (file_surf .NE. "0") THEN
    start = OMP_get_wtime()

    !nb_step, nb_atm, always shared.
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,OHvec_mat,nb_o,surf_mat,nb_surf)&
    !$OMP PRIVATE(s,i,j,k,tSOHvec_disp_vec,tSOHvec_norm)
    DO s = 1, nb_step
        F2:DO i = 1, nb_o*3
            IF (OHvec_mat(1,i,s) .EQ. 0) THEN
                CYCLE F2
            END IF
            DO j = 1, nb_surf(s)
                IF (surf_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        tSOHvec_disp_vec(k) = surf_mat(k,j,s) - OHvec_mat(k+8,i,s)
                        tSOHvec_disp_vec(k) = tSOHvec_disp_vec(k) - box(k) * ANINT(tSOHvec_disp_vec(k)/box(k))
                    END DO
                    tSOHvec_norm = NORM2(tSOHvec_disp_vec)
                    IF ( (tSOHvec_norm .LT. OHvec_mat(24,i,s)) .OR. OHvec_mat(24,i,s) .EQ. 0.0 ) THEN
                        OHvec_mat(24,i,s) = tSOHvec_norm
                        IF (OHvec_mat(11,i,s) .LT. surf_mat(3,j,s)) THEN
                            OHvec_mat(25,i,s) = -1
                        ELSE
                            OHvec_mat(25,i,s) = 1
                        END IF
                        OHvec_mat(26,i,s) = surf_mat(3,j,s)
                    END IF
                ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        tSOHvec_disp_vec(k) = surf_mat(k,j,s) - OHvec_mat(k+8,i,s)
                        tSOHvec_disp_vec(k) = tSOHvec_disp_vec(k) - box(k) * ANINT(tSOHvec_disp_vec(k)/box(k))
                    END DO
                    tSOHvec_norm = NORM2(tSOHvec_disp_vec)
                    IF ( (tSOHvec_norm .LT. OHvec_mat(27,i,s)) .OR. OHvec_mat(27,i,s) .EQ. 0.0 ) THEN
                        OHvec_mat(27,i,s) = tSOHvec_norm
                        IF (OHvec_mat(11,i,s) .GT. surf_mat(3,j,s)) THEN
                            OHvec_mat(28,i,s) = -1
                        ELSE
                            OHvec_mat(28,i,s) = 1
                        END IF
                        OHvec_mat(29,i,s) = surf_mat(3,j,s)
                    END IF
                END IF
            END DO
        END DO F2
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F15.5,A20)', "Proximity IS and OH groups:"&
        ,finish-start,"seconds elapsed"
END IF

! X ----------------------------------------------- Calculate closest distance between IS and any O atom
IF (file_surf .NE. "0") THEN
    start = OMP_get_wtime()

    !nb_step, nb_atm, always shared.
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,nb_o,surf_mat,nb_surf)&
    !$OMP PRIVATE(s,i,j,k,tSOvec_disp_vec,tSOvec_norm)
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF (atm_mat(2,i,s) .EQ. 16) THEN
                DO j = 1, nb_surf(s)
                    IF (surf_mat(4,j,s) .EQ. 1) THEN
                        DO k = 1, 3
                            tSOvec_disp_vec(k) = surf_mat(k,j,s) - atm_mat(k+2,i,s)
                            tSOvec_disp_vec(k) = tSOvec_disp_vec(k) - box(k) * ANINT(tSOvec_disp_vec(k)/box(k))
                        END DO
                        tSOvec_norm = NORM2(tSOvec_disp_vec)
                        IF ( (tSOvec_norm .LT. atm_mat(18,i,s)) .OR. atm_mat(18,i,s) .EQ. 0.0 ) THEN
                            atm_mat(18,i,s) = tSOvec_norm
                            IF (atm_mat(5,i,s) .LT. surf_mat(3,j,s)) THEN
                                atm_mat(19,i,s) = -1
                            ELSE
                                atm_mat(19,i,s) = 1
                            END IF
                            atm_mat(20,i,s) = surf_mat(3,j,s)
                        END IF
                    ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                        DO k = 1, 3
                            tSOvec_disp_vec(k) = surf_mat(k,j,s) - atm_mat(k+2,i,s)
                            tSOvec_disp_vec(k) = tSOvec_disp_vec(k) - box(k) * ANINT(tSOvec_disp_vec(k)/box(k))
                        END DO
                        tSOvec_norm = NORM2(tSOvec_disp_vec)
                        IF ( (tSOvec_norm .LT. atm_mat(21,i,s)) .OR. atm_mat(21,i,s) .EQ. 0.0 ) THEN
                            atm_mat(21,i,s) = tSOvec_norm
                            IF (atm_mat(5,i,s) .GT. surf_mat(3,j,s)) THEN
                                atm_mat(22,i,s) = -1
                            ELSE
                                atm_mat(22,i,s) = 1
                            END IF
                            atm_mat(23,i,s) = surf_mat(3,j,s)
                        END IF
                    END IF
                END DO
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO

    finish = OMP_get_wtime()
    PRINT'(A40,F15.5,A20)', "Proximity IS and O atoms:"&
        ,finish-start,"seconds elapsed"
END IF

! X ----------------------------------------------- Write OH BOND
IF (hbond_output .EQ. 1) THEN
    start = OMP_get_wtime()

    OPEN(UNIT=31, FILE = suffix//"_O_hbonds.dat")
    WRITE(31,'(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10)') "O_id","step","Don","UDon","Acc"&
    ,"C","OE","OH","OA","C9","O_type", "DistToISGo", "DistToISAir"
    DO s = 1, nb_step
        DO i = 1, nb_atm
            IF (atm_mat(2,i,s) .EQ. 16) THEN
                WRITE(31,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E24.14,E24.14)')&
                INT(atm_mat(1,i,s)), s, INT(atm_mat(11,i,s))&
            , INT(atm_mat(12,i,s)), INT(atm_mat(10,i,s)), INT(atm_mat(13,i,s)), INT(atm_mat(14,i,s))&
            , INT(atm_mat(15,i,s)), INT(atm_mat(16,i,s)), INT(atm_mat(17,i,s)), INT(atm_mat(9,i,s))&
            , (atm_mat(18,i,s)*atm_mat(19,i,s)),(atm_mat(21,i,s)*atm_mat(22,i,s))
        END IF
        END DO
    END DO
    CLOSE(UNIT=31)

    OPEN(UNIT=32, FILE = suffix//"_OH_hbonds.dat")
    WRITE(32,'(A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A24,A24)')&
        "O_id","O_type","H_id","step","Don","UDon","TAcc","C","OE"&
    ,"OH","OA","C9","TDon", "TUDon", "TAcc", "DistToISGo", "DistToISAir"
    DO s = 1, nb_step
        DO i = 1, nb_max_OHvec(s)
            WRITE(32,'(I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,I10,E24.14,E24.14)')&
                INT(OHvec_mat(1,i,s)), INT(OHvec_mat(23,i,s))&
            , INT(OHvec_mat(2,i,s)), s, INT(OHvec_mat(17,i,s)) , INT(OHvec_mat(18,i,s)), INT(OHvec_mat(16,i,s))&
            , INT(OHvec_mat(19,i,s)), INT(OHvec_mat(20,i,s)), INT(OHvec_mat(21,i,s)), INT(OHvec_mat(22,i,s))&
            , INT(OHvec_mat(30,i,s)), INT(OHvec_mat(31,i,s)), INT(OHvec_mat(32,i,s)), INT(OHvec_mat(33,i,s))&
            , (OHvec_mat(24,i,s)*OHvec_mat(25,i,s)), (OHvec_mat(27,i,s)*OHvec_mat(28,i,s))
        END DO
    END DO
    CLOSE(UNIT=32)

    finish = OMP_get_wtime()
    PRINT'(A40,F15.5,A20)', "O/OH Hbonds output:"&
        ,finish-start,"seconds elapsed"
END IF
    
! G ----------------------------------------------- Density profiles
IF ( (file_surf .NE. "0") .AND. (density_output .EQ. 1) ) THEN
    start = OMP_get_wtime()

    ! GO-WATER
    ALLOCATE(dens_go(dens_step,nb_step))
    dens_go(:,:) = 0.0_dp

    ! AIR WATER
    ALLOCATE(dens_air(dens_step,nb_step))
    dens_air(:,:) = 0.0_dp

    !nb_step, always shared.
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(box,OHvec_mat,nb_o,dens_go,dens_air)&
    !$OMP PRIVATE(s,r,j,i,count_dens_go,count_dens_air)
    DO s = 1, nb_step
        r = -10.0_dp
        DO j = 1, dens_step
            count_dens_go = 0
            count_dens_air = 0
        G1:DO i = 1, nb_o*3
                IF (OHvec_mat(1,i,s) .EQ. 0) THEN
                    CYCLE G1
                END IF
                IF (OHvec_mat(23,i,s) .NE. 13) THEN
                    CYCLE G1
                END IF
                IF ( (OHvec_mat(24,i,s)*OHvec_mat(25,i,s) .GE. r) .AND. (OHvec_mat(24,i,s)*OHvec_mat(25,i,s) .LT. r+dr) ) THEN
                    count_dens_go = count_dens_go + 1
                END IF
                IF ( (OHvec_mat(27,i,s)*OHvec_mat(28,i,s) .GE. r) .AND. (OHvec_mat(27,i,s)*OHvec_mat(28,i,s) .LT. r+dr) ) THEN
                    count_dens_air = count_dens_air + 1
                END IF
            END DO G1
            dens_go(j,s) = (18.0 * count_dens_go ) / (box(1) *  box(2) * dr * 1d-24 * 6.02214086d23)
            dens_air(j,s) = (18.0 * count_dens_air ) / (box(1) *  box(2) * dr * 1d-24 * 6.02214086d23)
            r = r + dr
        END DO
    END DO
    !$OMP END PARALLEL DO

    dens_go = dens_go / 2.0
    dens_air = dens_air / 2.0

    ALLOCATE(avg_dens_go(dens_step))
    avg_dens_go(:) = 0.0_dp
    ALLOCATE(avg_dens_air(dens_step))
    avg_dens_air(:) = 0.0_dp

    DO j = 1, dens_step
        avg_dens_go(j) = SUM(dens_go(j,:)) / nb_step
        avg_dens_air(j) = SUM(dens_air(j,:)) / nb_step
    END DO

    OPEN(UNIT=41, FILE = suffix//"_density_profile_go.dat")
    OPEN(UNIT=42, FILE = suffix//"_density_profile_air.dat")
    WRITE(41, '(A24,A24,A24,A24)') "step","[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
    WRITE(42, '(A24,A24,A24,A24)') "step","[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
    DO s = 1, nb_step
        DO j = 1, dens_step
            WRITE(41, '(I24,E24.14,E24.14,E24.14)') s, (-10.0 + (j-1) * dr), (-10.0 + j * dr), dens_go(j,s)
            WRITE(42, '(I24,E24.14,E24.14,E24.14)') s, (-10.0 + (j-1) * dr), (-10.0 + j * dr), dens_air(j,s)
        END DO
    END DO
    CLOSE(UNIT=41)
    CLOSE(UNIT=42)

    DEALLOCATE(dens_go,dens_air)

    OPEN(UNIT=43, FILE = suffix//"_avg_density_profile_go.dat")
    OPEN(UNIT=44, FILE = suffix//"_avg_density_profile_air.dat")
    WRITE(43, '(A24,A24,A24)') "[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
    WRITE(44, '(A24,A24,A24)') "[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
    DO j = 1, dens_step
        WRITE(43, '(E24.14,E24.14,E24.14)') (-10.0 + (j-1) * dr), (-10.0 + j * dr), avg_dens_go(j)
        WRITE(44, '(E24.14,E24.14,E24.14)') (-10.0 + (j-1) * dr), (-10.0 + j * dr), avg_dens_air(j)
    END DO
    CLOSE(UNIT=43)
    CLOSE(UNIT=44)

    DEALLOCATE(avg_dens_go,avg_dens_air)

    finish = OMP_get_wtime()
    PRINT'(A40,F15.5,A20)', "Density profiles:"&
        ,finish-start,"seconds elapsed"
END IF

! H ----------------------------------------------- VVCF
IF (vvcf_c .EQ. 1 ) THEN
    start = OMP_get_wtime()
    mct = 2500 ! in femtosecond
    mctb = 0 ! in femtosecond
    mcs = INT(mct / timestep_fs)
    mcsb = INT(mctb / timestep_fs) + 1

    IF (mcs+1 .GE. nb_step) THEN
        PRINT*, "Error: Max correlation time > trajectory time"
        STOP
    END IF

    ALLOCATE(vvcf_xxz(mcs+1))
    ALLOCATE(timings(mcs+1))
    timings(:) = 0.0_dp
    !nb_step is a parameter so always shared
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(NONE) SHARED(vvcf_xxz,mcsb,mcs,OHvec_mat,nb_o,box,timings)&
    !$OMP PRIVATE(t,s,i,j,tij_vec,trij, OH_vel_vec, OH_disp_vec,l,v,u,start_i,finish_i)
    DO t = mcsb, mcs+1
        start_i = OMP_get_wtime()
        vvcf_xxz(t) = 0.0_dp
        ALLOCATE(v(nb_o*3,nb_step))
        v(:,:) = 0
        DO s = 1, nb_step
            IF (s-1+t-1 .LT. nb_step) THEN
             H1:DO i = 1, nb_o*3

                    IF (OHvec_mat(1,i,s) .EQ. 0) THEN ! Skip empty
                        CYCLE H1
                    END IF
                    IF (OHvec_mat(23,i,s) .NE. 13) THEN ! Select water only
                        CYCLE H1
                    END IF
                    IF (OHvec_mat(30,i,s) .NE. 1) THEN ! Select water close to the C surface (<=9A)
                        CYCLE H1
                    END IF

                    IF (layers .EQ. 1 ) THEN

                        IF ( (OHvec_mat(24,i,s)*OHvec_mat(25,i,s) .GT. surf_up) .OR.&
                        (OHvec_mat(24,i,s)*OHvec_mat(25,i,s) .LE. surf_down)) THEN ! Surface distance (go)
                            CYCLE H1
                        END IF

                        IF (layers_c .EQ. 1 ) THEN

                            IF (s .EQ. 1) THEN ! init
                                DO u = s, s+t-1
                                    IF ( (OHvec_mat(24,i,u)*OHvec_mat(25,i,u) .GT. surf_up) .OR.&
                                    (OHvec_mat(24,i,u)*OHvec_mat(25,i,u) .LE.  surf_down)) THEN
                                        CYCLE H1
                                    END IF
                                    v(i,s) = v(i,s) + 1
                                END DO
                                IF (v(i,s) .NE. t) THEN
                                    PRINT*, "Very incorrect"
                                    CYCLE H1
                                END IF
                            ELSE IF (v(i,s-1) .LE. 1) THEN
                                DO u = s, s+t-1
                                    IF ( (OHvec_mat(24,i,u)*OHvec_mat(25,i,u) .GT. surf_up) .OR.&
                                    (OHvec_mat(24,i,u)*OHvec_mat(25,i,u) .LE.  surf_down)) THEN
                                        CYCLE H1
                                    END IF
                                    v(i,s) = v(i,s) + 1
                                END DO
                                IF (v(i,s) .NE. t) THEN
                                    PRINT*, "Very incorrect"
                                    CYCLE H1
                                END IF
                            ELSE IF (v(i,s-1) .EQ. t) THEN
                                IF ( (OHvec_mat(24,i,s+t-1)*OHvec_mat(25,i,s+t-1) .GT. surf_up) .OR.&
                                (OHvec_mat(24,i,s+t-1)*OHvec_mat(25,i,s+t-1) .LE.  surf_down)) THEN
                                    CYCLE H1
                                END IF
                                v(i,s) = t
                            ELSE
                                CYCLE H1
                            END IF

                        END IF
                    END IF

                    IF (hbonds .EQ. 1) THEN

                        IF ( (OHvec_mat(32,i,s) .NE. UDC ) .OR.& ! Exclude S (1) or D (2) (Unique hbond, O-H water engaged in any nb of Hbonds)
                        (OHvec_mat(33,i,s) .NE. AC ) ) THEN ! Exclude nb of incoming Hbonds (Acceptor)
                            CYCLE H1
                        END IF

                        IF (hbonds_c .EQ. 1) THEN

                            IF (s .EQ. 1) THEN ! init
                                DO u = s, s+t-1
                                    IF ( (OHvec_mat(32,i,u) .NE. UDC ) .OR.&
                                    (OHvec_mat(33,i,u) .NE. AC ) ) THEN
                                        CYCLE H1
                                    END IF
                                    v(i,s) = v(i,s) + 1
                                END DO
                                IF (v(i,s) .NE. t) THEN
                                    PRINT*, "Very incorrect"
                                    CYCLE H1
                                END IF
                            ELSE IF (v(i,s-1) .LE. 1) THEN
                                DO u = s, s+t-1
                                    IF ( (OHvec_mat(32,i,u) .NE. UDC ) .OR.&
                                    (OHvec_mat(33,i,u) .NE. AC ) ) THEN
                                        CYCLE H1
                                    END IF
                                    v(i,s) = v(i,s) + 1
                                END DO
                                IF (v(i,s) .NE. t) THEN
                                    PRINT*, "Very incorrect"
                                    CYCLE H1
                                END IF
                            ELSE IF (v(i,s-1) .EQ. t) THEN
                                IF ( (OHvec_mat(32,i,s+t-1) .NE. UDC ) .OR.&
                                (OHvec_mat(33,i,s+t-1) .NE. AC ) ) THEN
                                    CYCLE H1
                                END IF
                                v(i,s) = t
                            ELSE
                                CYCLE H1
                            END IF

                        END IF
                    END IF

                 H2:DO j = 1, nb_o * 3

                        IF ( (OHvec_mat(1,j,s) .EQ. 0) ) THEN ! Skip empty
                            CYCLE H2
                        END IF

                        ! Trick for intra cross-correlation only
                        IF (r_vvcf_cut .LE. 2.05) THEN
                            IF (OHvec_mat(1,i,s) .NE. OHvec_mat(1,j,s)) THEN
                                CYCLE H2
                            END IF
                        END IF

                        DO k = 1, 3
                            tij_vec(k) = OHvec_mat(k+8,j,s) - OHvec_mat(k+8,i,s)
                            tij_vec(k) = tij_vec(k) - box(k) * ANINT(tij_vec(k)/box(k))
                        END DO
                        trij = NORM2(tij_vec(:))

                        IF ( trij .LT. r_vvcf_cut ) THEN

                         H3:DO l = 1, nb_o *3
                                IF ( (OHvec_mat(1,j,s) .NE. OHvec_mat(1,l,s+t-1)) .OR.&
                                 (OHvec_mat(2,j,s) .NE. OHvec_mat(2,l,s+t-1)) ) THEN
                                    CYCLE H3
                                END IF

                                DO k = 1, 3
                                    OH_vel_vec(k) =  OHvec_mat(k+5,l,s+t-1)
                                    OH_disp_vec(k) = OHvec_mat(k+2,l,s+t-1)
                                END DO

                                vvcf_xxz(t) = vvcf_xxz(t) + (OHvec_mat(8,i,s) &
                                * &
                                (DOT_PRODUCT(OH_vel_vec(:),OH_disp_vec(:)) &
                                / &
                                NORM2(OH_disp_vec(:))))

                            END DO H3 ! j at t

                        END IF

                    END DO H2 ! j at 0

                END DO H1 ! i at 0

            END IF
        END DO ! Time
        DEALLOCATE(v)
        vvcf_xxz(t) = vvcf_xxz(t) * 1.0 / (nb_step - (t-1))

        finish_i = OMP_get_wtime()
        timings(t)=finish_i-start_i

        IF (MODULO(t,25) .EQ. 0) THEN
            PRINT('(I10,A1,I10,E24.14,E24.14,E24.14)'),t,"/",mcs+1, (t-1)*timestep_fs, vvcf_xxz(t),timings(t)
        END IF

    END DO ! Corr
    !$OMP END PARALLEL DO

    avg_timigs = (SUM(timings(:)) / (mcs+1-mcsb) )
    finish = OMP_get_wtime()
    
    PRINT'(A40,F15.5,A20,A20,F15.5)', "Done with VVCF_xxz:",finish-start,"seconds elapsed","avg per step:",avg_timigs

    start = OMP_get_wtime()

    open(UNIT=51, FILE = suffix//"_vvcf_xxz.dat")
    WRITE(51,'(A20,A20)') "Time (fs)","VVCF_xxz (Å2.fs2)"
    DO t = mcsb, mcs+1
        WRITE(51,'(E24.14,E24.14)') (t-1)*timestep_fs, vvcf_xxz(t)
    END DO
    CLOSE(UNIT=51)

    finish = OMP_get_wtime()
    PRINT'(A40,F15.5,A20)', "Done with VVCF_xxz output:",finish-start,"seconds elapsed"

DEALLOCATE(vvcf_xxz)
END IF

! ----------------------------------------------- Deallocate and exit
DEALLOCATE(OHvec_mat,atm_mat)
END PROGRAM vvcf
