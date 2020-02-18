! Author:Rolf David
! License:MIT
! UTF-8, CRLF, Fortran2003, OpenMP

MODULE INPUT

IMPLICIT NONE
! -----------------------------------------------Set Double precision
INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)
! -----------------------------------------------Reading related variables
CHARACTER(LEN=100)          :: label, value
INTEGER                     :: iostatus=0
INTEGER                     :: line_c=0
INTEGER                     :: delim
INTEGER, PARAMETER          :: file_unit=99
! ----------------------------------------------- Parameters
REAL(dp), PARAMETER         :: bohr_to_angstrom=0.529177249_dp
REAL(dp), PARAMETER         :: aut_to_fs=0.0241888432658569977_dp
REAL(dp), PARAMETER         :: c_pi=4.0_dp*DATAN(1.0_dp)

! -----------------------------------------------Variables
CHARACTER(LEN=64)           :: file_pos='0', file_vel='0', file_is='0'
INTEGER                     :: nb_atm, nb_step
CHARACTER(LEN=2)            :: suffix
REAL(dp)                    :: xlo=0.0_dp, xhi=0.0_dp
REAL(dp)                    :: ylo=0.0_dp, yhi=0.0_dp
REAL(dp)                    :: zlo=0.0_dp, zhi=0.0_dp
REAL(dp)                    :: box(3)
! Surface related
CHARACTER(LEN=1)            :: IS_c
CHARACTER(LEN=1)            :: AS_c
! Assign
CHARACTER(LEN=3)            :: assign_center_el
REAL(dp)                    :: assign_OpH_rcut, assign_OpC_rcut, assign_OpSi_rcut
INTEGER                     :: waterlist
! Density
REAL(dp)                    :: dens_dr=0.0_dp, dens_rstart=0.0_dp
INTEGER                     :: dens_step=75
INTEGER                     :: dens_center_atmnb
! Water Angle
REAL(dp)                    :: wa_X1Owat_rcut, wa_X2Owat_rcut
INTEGER                     :: wa_AS_center_atmnb
! Hbonds
REAL(dp)                    :: hb_oHpO_rcut, hb_OpH_rcut
REAL(dp)                    :: hb_X1Oh_rcut, hb_X2Oh_rcut
REAL(dp)                    :: hb_X1O_rcut, hb_X2O_rcut
! VVCF
REAL(dp)                    :: vvcf_oHpO_rcut, vvcf_OpH_rcut
REAL(dp)                    :: vvcf_X1Oh_rcut, vvcf_X2Oh_rcut
REAL(dp)                    :: vvcf_rcut
REAL(dp)                    :: timestep_fs=0.5_dp, mct, mctb
CHARACTER(LEN=1)            :: hbonds_s, hbonds_c, layers_s, layers_c
CHARACTER(LEN=1)            :: hbonds_double="N"
CHARACTER(LEN=1)            :: intra_only, water_only, close_c_only
CHARACTER(LEN=1)            :: close_gl_only, close_gol_only, close_ol_only
CHARACTER(LEN=1)            :: up_down_only="U"
INTEGER                     :: ACLE, ACGT, DCLE, DCGT
INTEGER                     :: ACLE2, ACGT2, DCLE2, DCGT2
REAL(dp)                    :: layer_down, layer_up
! -----------------------------------------------Counters
INTEGER, PRIVATE            :: i

CONTAINS

SUBROUTINE READINPUTSUB(input_file)

    CHARACTER(LEN=100)      :: input_file

    OPEN(file_unit, file=input_file)
    PRINT'(A100)', 'Input parameters'
    PRINT'(A100)', '--------------------------------------------------' &
    , '--------------------------------------------------'
    RL:DO WHILE (iostatus == 0)
        READ(file_unit, '(A)', IOSTAT=iostatus) value
        IF (iostatus == 0)then
            line_c=line_c + 1

            delim=SCAN(value, '    ')
            label=value(1:delim)
            value=value(delim + 1:)
            IF (label(1:1) == '!')THEN
                CYCLE RL
            END IF
            SELECT CASE (label)
                CASE ('file_pos')
                    READ(value, * , IOSTAT=iostatus) file_pos
                    PRINT'(A50,A64)', 'XYZ positions file:', ADJUSTR(file_pos)
                CASE ('file_vel')
                    READ(value, * , IOSTAT=iostatus) file_vel
                    PRINT'(A50,A64)', 'XYZ velocities file:', ADJUSTR(file_vel)
                CASE ('file_is')
                    READ(value, * , IOSTAT=iostatus) file_is
                    PRINT'(A50,A64)', 'XYZ velocities file:', ADJUSTR(file_is)
                CASE ('nb_atm')
                    READ(value, * , IOSTAT=iostatus) nb_atm
                    PRINT'(A50,I64)', 'Number of Atoms: ', nb_atm
                CASE ('nb_step')
                    READ(value, * , IOSTAT=iostatus) nb_step
                    PRINT'(A50,I64)', 'Number of steps: ', nb_step
                CASE('IS_c')
                    READ(value, * , IOSTAT=iostatus) IS_c
                    PRINT'(A50,A64)', 'IS_c:', IS_c
                CASE('AS_c')
                    READ(value, * , IOSTAT=iostatus) AS_c
                    PRINT'(A50,A64)', 'AS_c:', AS_c
                CASE ('suffix')
                    READ(value, * , IOSTAT=iostatus) suffix
                    PRINT'(A50,A64)', 'Suffix: ', suffix
                CASE ('xlo')
                    READ(value, * , IOSTAT=iostatus) xlo
                    PRINT'(A50,E64.10)', 'Box_X_low: ', xlo
                CASE ('xhi')
                    READ(value, * , IOSTAT=iostatus) xhi
                    PRINT'(A50,E64.10)', 'Box_X_high: ', xhi
                CASE ('ylo')
                    READ(value, * , IOSTAT=iostatus) ylo
                    PRINT'(A50,E64.10)', 'Box_Y_low: ', ylo
                CASE ('yhi')
                    READ(value, * , IOSTAT=iostatus) yhi
                    PRINT'(A50,E64.10)', 'Box_Y_high: ', yhi
                CASE ('zlo')
                    READ(value, * , IOSTAT=iostatus) zlo
                    PRINT'(A50,E64.10)', 'Box_Z_low: ', zlo
                CASE ('zhi')
                    READ(value, * , IOSTAT=iostatus) zhi
                    PRINT'(A50,E64.10)', 'Box_Z_high: ', zhi
! Assign
                CASE('assign_center_el')
                    READ(value, * , IOSTAT=iostatus) assign_center_el
                    PRINT'(A50,A64)', 'assign_center_el:', ADJUSTR(assign_center_el)
                CASE ('assign_OpH_rcut')
                    READ(value, * , IOSTAT=iostatus) assign_OpH_rcut
                    PRINT'(A50,E64.2)', 'assign_OpH_rcut: ', assign_OpH_rcut
                CASE ('assign_OpC_rcut')
                    READ(value, * , IOSTAT=iostatus) assign_OpC_rcut
                    PRINT'(A50,E64.2)', 'assign_OpC_rcut: ', assign_OpC_rcut
                CASE ('assign_OpSi_rcut')
                    READ(value, * , IOSTAT=iostatus) assign_OpSi_rcut
                    PRINT'(A50,E64.2)', 'assign_OpSi_rcut: ', assign_OpSi_rcut
                CASE ('waterlist')
                    READ(value, * , IOSTAT=iostatus) waterlist
                    PRINT'(A50,I64)', 'waterlist: ', waterlist
! Density
                CASE ('dens_center_atmnb')
                    READ(value, * , IOSTAT=iostatus) dens_center_atmnb
                    PRINT'(A50,I64)', 'dens_center_atmnb: ', dens_center_atmnb
                CASE ('dens_dr')
                    READ(value, * , IOSTAT=iostatus) dens_dr
                    PRINT'(A50,E64.2)', 'dens_dr: ', dens_dr
                CASE ('dens_step')
                    READ(value, * , IOSTAT=iostatus) dens_step
                    PRINT'(A50,I64)', 'dens_step: ', dens_step
                CASE ('dens_rstart')
                    READ(value, * , IOSTAT=iostatus) dens_rstart
                    PRINT'(A50,E64.2)', 'dens_rstart: ', dens_rstart
! Water Angle
                CASE ('wa_AS_center_atmnb')
                    READ(value, * , IOSTAT=iostatus) wa_AS_center_atmnb
                    PRINT'(A50,I64)', 'wa_AS_center_atmnb: ', wa_AS_center_atmnb
                CASE ('wa_X1Owat_rcut')
                    READ(value, * , IOSTAT=iostatus) wa_X1Owat_rcut
                    PRINT'(A50,E64.2)', 'wa_X1Owat_rcut: ', wa_X1Owat_rcut
                CASE ('wa_X2Owat_rcut')
                    READ(value, * , IOSTAT=iostatus) wa_X2Owat_rcut
                    PRINT'(A50,E64.2)', 'wa_X2Owat_rcut: ', wa_X2Owat_rcut
! HBonds
                CASE ('hb_oHpO_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_oHpO_rcut
                    PRINT'(A50,E64.2)', 'hb_oHpO_rcut: ', hb_oHpO_rcut
                CASE ('hb_OpH_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_OpH_rcut
                    PRINT'(A50,E64.2)', 'hb_OpH_rcut: ', hb_OpH_rcut
                CASE ('hb_X1Oh_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_X1Oh_rcut
                    PRINT'(A50,E64.2)', 'hb_X1Oh_rcut: ', hb_X1Oh_rcut
                CASE ('hb_X2Oh_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_X2Oh_rcut
                    PRINT'(A50,E64.2)', 'hb_X2Oh_rcut: ', hb_X2Oh_rcut
                CASE ('hb_X1O_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_X1O_rcut
                    PRINT'(A50,E64.2)', 'hb_X1O_rcut: ', hb_X1O_rcut
                CASE ('hb_X2O_rcut')
                    READ(value, * , IOSTAT=iostatus) hb_X2O_rcut
                    PRINT'(A50,E64.2)', 'hb_X2O_rcut: ', hb_X2O_rcut
! VVCF
                CASE ('vvcf_oHpO_rcut')
                    READ(value, * , IOSTAT=iostatus) vvcf_oHpO_rcut
                    PRINT'(A50,E64.2)', 'vvcf_oHpO_rcut: ', vvcf_oHpO_rcut
                CASE ('vvcf_OpH_rcut')
                    READ(value, * , IOSTAT=iostatus) vvcf_OpH_rcut
                    PRINT'(A50,E64.2)', 'vvcf_OpH_rcut: ', vvcf_OpH_rcut
                CASE ('vvcf_X1Oh_rcut')
                    READ(value, * , IOSTAT=iostatus) vvcf_X1Oh_rcut
                    PRINT'(A50,E64.2)', 'vvcf_X1Oh_rcut: ', vvcf_X1Oh_rcut
                CASE ('vvcf_X2Oh_rcut')
                    READ(value, * , IOSTAT=iostatus) vvcf_X2Oh_rcut
                    PRINT'(A50,E64.2)', 'vvcf_X2Oh_rcut: ', vvcf_X2Oh_rcut
                CASE ('vvcf_rcut')
                    READ(value, * , IOSTAT=iostatus) vvcf_rcut
                    PRINT'(A50,E64.2)', 'vvcf_rcut: ', vvcf_rcut
                CASE ('mct')
                    READ(value, * , IOSTAT=iostatus) mct
                    PRINT'(A50,E64.2)', 'mct: ', mct
                CASE ('mctb')
                    READ(value, * , IOSTAT=iostatus) mctb
                    PRINT'(A50,E64.2)', 'mctb: ', mctb
                CASE ('timestep_fs')
                    READ(value, * , IOSTAT=iostatus) timestep_fs
                    PRINT'(A50,E64.2)', 'timestep_fs: ', timestep_fs
                CASE ('hbonds')
                    READ(value, * , IOSTAT=iostatus) hbonds_s
                    PRINT'(A50,A64)', 'hbonds: ', hbonds_s
                CASE ('hbonds_c')
                    READ(value, * , IOSTAT=iostatus) hbonds_c
                    PRINT'(A50,A64)', 'hbonds_c: ', hbonds_c
                CASE ('layers')
                    READ(value, * , IOSTAT=iostatus) layers_s
                    PRINT'(A50,A64)', 'layers: ', layers_s
                CASE ('layers_c')
                    READ(value, * , IOSTAT=iostatus) layers_c
                    PRINT'(A50,A64)', 'layers_c: ', layers_c
                CASE ('intra_only')
                    READ(value, * , IOSTAT=iostatus) intra_only
                    PRINT'(A50,A64)', 'intra_only: ', intra_only
                CASE ('water_only')
                    READ(value, * , IOSTAT=iostatus) water_only
                    PRINT'(A50,A64)', 'water_only: ', water_only
                CASE ('close_c_only')
                    READ(value, * , IOSTAT=iostatus) close_c_only
                    PRINT'(A50,A64)', 'close_c_only: ', close_c_only
                CASE ('close_gl_only')
                    READ(value, * , IOSTAT=iostatus) close_gl_only
                    PRINT'(A50,A64)', 'close_gl_only: ', close_gl_only
                CASE ('close_gol_only')
                    READ(value, * , IOSTAT=iostatus) close_gol_only
                    PRINT'(A50,A64)', 'close_gol_only: ', close_gol_only
                CASE ('close_ol_only')
                    READ(value, * , IOSTAT=iostatus) close_ol_only
                    PRINT'(A50,A64)', 'close_ol_only: ', close_ol_only
                CASE ('up_down_only')
                    READ(value, * , IOSTAT=iostatus) up_down_only
                    PRINT'(A50,A64)', 'up_down_only: ', up_down_only
                CASE ('DCLE')
                    READ(value, * , IOSTAT=iostatus) DCLE
                    PRINT'(A50,I64)', 'DCLE: ', DCLE
                CASE ('ACLE')
                    READ(value, * , IOSTAT=iostatus) ACLE
                    PRINT'(A50,I64)', 'ACLE: ', ACLE
                CASE ('DCGT')
                    READ(value, * , IOSTAT=iostatus) DCGT
                    PRINT'(A50,I64)', 'DCGT: ', DCGT
                CASE ('ACGT')
                    READ(value, * , IOSTAT=iostatus) ACGT
                    PRINT'(A50,I64)', 'ACGT: ', ACGT
                CASE ('hbonds_double')
                    READ(value, * , IOSTAT=iostatus) hbonds_double
                    PRINT'(A50,A64)', 'hbonds_double: ', hbonds_double
                CASE ('DCLE2')
                    READ(value, * , IOSTAT=iostatus) DCLE2
                    PRINT'(A50,I64)', 'DCLE2: ', DCLE2
                CASE ('ACLE2')
                    READ(value, * , IOSTAT=iostatus) ACLE2
                    PRINT'(A50,I64)', 'ACLE2: ', ACLE2
                CASE ('DCGT2')
                    READ(value, * , IOSTAT=iostatus) DCGT2
                    PRINT'(A50,I64)', 'DCGT2: ', DCGT2
                CASE ('ACGT2')
                    READ(value, * , IOSTAT=iostatus) ACGT2
                    PRINT'(A50,I64)', 'ACGT2: ', ACGT2
                CASE ('layer_down')
                    READ(value, * , IOSTAT=iostatus) layer_down
                    PRINT'(A50,E64.2)', 'layer_down: ', layer_down
                CASE ('layer_up')
                    READ(value, * , IOSTAT=iostatus) layer_up
                    PRINT'(A50,E64.2)', 'layer_up: ', layer_up
                CASE DEFAULT
                    PRINT'(A50,I64)', 'Invalid label, line:', line_c
            END SELECT
        END IF
    END DO RL
    PRINT'(A100)', '--------------------------------------------------' &
    , '--------------------------------------------------'
    ! -----------------------------------------------Calculate the box size
    box(1)=xhi - xlo
    box(2)=yhi - ylo
    box(3)=zhi - zlo

    DO i=1, 3
        IF (box(i) .LT. 0.1)THEN
            PRINT'(A100)', '--------------------------------------------------' &
            , '--------------------------------------------------'
            PRINT'(A100)', "The box seems very small, please check your input file"
            PRINT'(A100)', '--------------------------------------------------' &
            , '--------------------------------------------------'
            STOP
        END IF
    END DO

END SUBROUTINE READINPUTSUB

END MODULE INPUT