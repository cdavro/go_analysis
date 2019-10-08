!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

PROGRAM density
USE OMP_LIB
USE INPUT_MOD

IMPLICIT NONE

! ----------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)

! ----------------------------------------------- Timings
REAL(dp)                        :: start,finish

! ----------------------------------------------- Input files
CHARACTER(LEN=100)              :: input_file
INTEGER                         :: CAC
CHARACTER(LEN=2)                :: dummy_char

! ----------------------------------------------- Infos/properties
REAL(dp), ALLOCATABLE           :: atm_mat(:,:,:), surf_mat(:,:,:)
CHARACTER(LEN=2), ALLOCATABLE   :: atm_el(:)
INTEGER                         :: nb_line, nb_max_pt
INTEGER, ALLOCATABLE            :: nb_surf(:)
REAL(dp), ALLOCATABLE           :: avg_z(:)
REAL(dp), ALLOCATABLE           :: dens_go(:,:), dens_air(:,:), avg_dens_go(:), avg_dens_air(:)
REAL(dp), ALLOCATABLE           :: dens_go_c(:,:), avg_dens_go_c(:)

REAL(dp)                        :: tSOvec_disp_vec(3), tSOvec_norm
REAL(dp)                        :: r
INTEGER                         :: count_dens_go, count_dens_air, count_dens_go_c

! ----------------------------------------------- Counters
INTEGER                         :: i, s, k, j, o

! ----------------------------------------------- Get arguments (filenames, choices)
CAC = COMMAND_ARGUMENT_COUNT()

IF (CAC .EQ. 0) THEN
    PRINT*,"No input files"
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, input_file)
input_file=TRIM(input_file)
CALL READINPUTSUB(input_file)
file_pos=TRIM(file_pos)
file_surf=TRIM(file_surf)

! ----------------------------------------------- Controls
! To Do

! ----------------------------------------------- Allocate function for reading files
ALLOCATE(atm_mat(12,nb_atm,nb_step))
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
PRINT'(A40,F14.2,A20)', "Positions:",finish-start,"seconds elapsed"

! A ----------------------------------------------- Since the number of points for the IS isn't constant, count it.
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

ALLOCATE(surf_mat(4,nb_max_pt,nb_step))
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

! B ----------------------------------------------- Calculate closest distance between IS and any O atom
start = OMP_get_wtime()

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,box,surf_mat,nb_surf,nb_step,nb_atm)&
!$OMP PRIVATE(s,i,j,k,tSOvec_disp_vec,tSOvec_norm)
DO s = 1, nb_step
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 16) THEN
            DO j = 1, nb_surf(s)
                IF (surf_mat(4,j,s) .EQ. 1) THEN
                    DO k = 1, 3
                        tSOvec_disp_vec(k) = surf_mat(k,j,s) - atm_mat(k+3,i,s)
                        tSOvec_disp_vec(k) = tSOvec_disp_vec(k) - box(k) * ANINT(tSOvec_disp_vec(k)/box(k))
                    END DO
                    tSOvec_norm = NORM2(tSOvec_disp_vec)
                    IF ( (tSOvec_norm .LT. atm_mat(7,i,s)) .OR. atm_mat(7,i,s) .EQ. 0.0 ) THEN
                        atm_mat(7,i,s) = tSOvec_norm
                        IF (atm_mat(6,i,s) .LT. surf_mat(3,j,s)) THEN
                            atm_mat(8,i,s) = -1
                        ELSE
                            atm_mat(8,i,s) = 1
                        END IF
                        atm_mat(9,i,s) = surf_mat(5,j,s)
                    END IF
                ELSE IF (surf_mat(4,j,s) .EQ. 2) THEN
                    DO k = 1, 3
                        tSOvec_disp_vec(k) = surf_mat(k,j,s) - atm_mat(k+3,i,s)
                        tSOvec_disp_vec(k) = tSOvec_disp_vec(k) - box(k) * ANINT(tSOvec_disp_vec(k)/box(k))
                    END DO
                    tSOvec_norm = NORM2(tSOvec_disp_vec)
                    IF ( (tSOvec_norm .LT. atm_mat(10,i,s)) .OR. atm_mat(10,i,s) .EQ. 0.0 ) THEN
                        atm_mat(10,i,s) = tSOvec_norm
                        IF (atm_mat(6,i,s) .GT. surf_mat(3,j,s)) THEN
                            atm_mat(11,i,s) = -1
                        ELSE
                            atm_mat(11,i,s) = 1
                        END IF
                        atm_mat(12,i,s) = surf_mat(5,j,s)
                    END IF
                END IF
            END DO
        END IF
    END DO
END DO
!$OMP END PARALLEL DO

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Proximity IS and O atoms:"&
    ,finish-start,"seconds elapsed"

! C ----------------------------------------------- Density profiles
start = OMP_get_wtime()

! GO-WATER
ALLOCATE(dens_go(dens_step,nb_step))
dens_go(:,:) = 0.0_dp

! AIR WATER
ALLOCATE(dens_air(dens_step,nb_step))
dens_air(:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,nb_atm,nb_step,box,dens_go,dens_air,dens_rstart,dens_dr,dens_step)&
!$OMP PRIVATE(s,r,j,i,count_dens_go,count_dens_air)
DO s = 1, nb_step
    r = dens_rstart
    DO j = 1, dens_step
        count_dens_go = 0
        count_dens_air = 0
    C1:DO i = 1, nb_atm
            IF (atm_mat(3,i,s) .NE. 13) THEN
                CYCLE C1
            END IF
            IF ( (atm_mat(7,i,s)*atm_mat(8,i,s) .GE. r) .AND. (atm_mat(7,i,s)*atm_mat(8,i,s) .LT. r+dens_dr) ) THEN
                count_dens_go = count_dens_go + 1
            END IF
            IF ( (atm_mat(10,i,s)*atm_mat(11,i,s) .GE. r) .AND. (atm_mat(10,i,s)*atm_mat(11,i,s) .LT. r+dens_dr) ) THEN
                count_dens_air = count_dens_air + 1
            END IF
        END DO C1
        dens_go(j,s) = (18.0 * count_dens_go ) / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        dens_air(j,s) = (18.0 * count_dens_air ) / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        r = r + dens_dr
    END DO
END DO
!$OMP END PARALLEL DO

ALLOCATE(avg_dens_go(dens_step))
avg_dens_go(:) = 0.0_dp
ALLOCATE(avg_dens_air(dens_step))
avg_dens_air(:) = 0.0_dp

DO j = 1, dens_step
    avg_dens_go(j) = SUM(dens_go(j,:)) / nb_step
    avg_dens_air(j) = SUM(dens_air(j,:)) / nb_step
END DO

OPEN(UNIT=41, FILE = suffix//"_density_profile_go.txt")
OPEN(UNIT=42, FILE = suffix//"_density_profile_air.txt")
WRITE(41, '(A24,A24,A24,A24)') "step","[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
WRITE(42, '(A24,A24,A24,A24)') "step","[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
DO s = 1, nb_step
    DO j = 1, dens_step
        WRITE(41, '(I24,E24.14,E24.14,E24.14)') s, (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), dens_go(j,s)
        WRITE(42, '(I24,E24.14,E24.14,E24.14)') s, (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), dens_air(j,s)
    END DO
END DO
CLOSE(UNIT=41)
CLOSE(UNIT=42)

DEALLOCATE(dens_go,dens_air)

OPEN(UNIT=43, FILE = suffix//"_avg_density_profile_go.txt")
OPEN(UNIT=44, FILE = suffix//"_avg_density_profile_air.txt")
WRITE(43, '(A24,A24,A24)') "[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
WRITE(44, '(A24,A24,A24)') "[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
DO j = 1, dens_step
    WRITE(43, '(E24.14,E24.14,E24.14)') (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), avg_dens_go(j)
    WRITE(44, '(E24.14,E24.14,E24.14)') (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), avg_dens_air(j)
END DO
CLOSE(UNIT=43)
CLOSE(UNIT=44)

DEALLOCATE(avg_dens_go,avg_dens_air)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Density profiles:"&
    ,finish-start,"seconds elapsed"

! D ----------------------------------------------- 
start = OMP_get_wtime()

! AVGz of sheet
ALLOCATE(dens_go_c(dens_step,nb_step))
ALLOCATE(avg_z(nb_step))

dens_go_c(:,:) = 0.0_dp
avg_z(:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(atm_mat,nb_atm,nb_step,box,dens_go_c,avg_z,dens_rstart,dens_dr,dens_step)&
!$OMP PRIVATE(s,r,j,i,count_dens_go_c,o)
DO s = 1, nb_step
    o = 0
    DO i = 1, nb_atm
        IF (atm_mat(2,i,s) .EQ. 12) THEN
            avg_z(s) = avg_z(s) + atm_mat(6,i,s) 
            o = o + 1
        END IF
    END DO
    avg_z(s) = avg_z(s) / o
    r = dens_rstart
    DO j = 1, dens_step
        count_dens_go_c = 0
    D1:DO i = 1, nb_atm
            IF (atm_mat(3,i,s) .NE. 13) THEN
                CYCLE D1
            END IF
            IF ( ((atm_mat(6,i,s) - avg_z(s)) .GE. r) .AND. ( (atm_mat(6,i,s) - avg_z(s)) .LT. r+dens_dr)) THEN
                count_dens_go_c = count_dens_go_c + 1
            END IF
        END DO D1
        dens_go_c(j,s) = (18.0 * count_dens_go_c ) / (box(1) *  box(2) * dens_dr * 1d-24 * 6.02214086d23)
        r = r + dens_dr
    END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(avg_z)

ALLOCATE(avg_dens_go_c(dens_step))
avg_dens_go_c(:) = 0.0_dp

DO j = 1, dens_step
    avg_dens_go_c(j) = SUM(dens_go_c(j,:)) / nb_step
END DO

OPEN(UNIT=41, FILE = suffix//"_density_profile_go_c.txt")
WRITE(41, '(A24,A24,A24,A24)') "step","[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
DO s = 1, nb_step
    DO j = 1, dens_step
        WRITE(41, '(I24,E24.14,E24.14,E24.14)') s, (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), dens_go_c(j,s)
    END DO
END DO
CLOSE(UNIT=41)

DEALLOCATE(dens_go_c)

OPEN(UNIT=43, FILE = suffix//"_avg_density_profile_go_c.txt")
WRITE(43, '(A24,A24,A24)') "[ r (Å)","r+dr (Å) [", "ρ/ρ(bulk)"
DO j = 1, dens_step
    WRITE(43, '(E24.14,E24.14,E24.14)') (dens_rstart + (j-1) * dens_dr), (dens_rstart + j * dens_dr), avg_dens_go_c(j)
END DO
CLOSE(UNIT=43)

DEALLOCATE(avg_dens_go_c)

finish = OMP_get_wtime()
PRINT'(A40,F14.2,A20)', "Density profiles:"&
    ,finish-start,"seconds elapsed"

! ----------------------------------------------- Deallocate and exit
DEALLOCATE(atm_mat,surf_mat)

END PROGRAM density