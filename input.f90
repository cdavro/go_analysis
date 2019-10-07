!Author: Rolf David
!License: MIT
!UTF-8, CRLF, Fortran2003, OpenMP

MODULE INPUT_MOD

IMPLICIT NONE

! ----------------------------------------------- Set Double precision
INTEGER, PARAMETER   :: dp=KIND(0.0d0)

! ----------------------------------------------- Reading related variables
CHARACTER(LEN=100)   :: label, value
INTEGER              :: iostatus = 0
INTEGER              :: line_c = 0
INTEGER              :: delim
INTEGER, PARAMETER   :: file_unit = 99

! ----------------------------------------------- Variables
CHARACTER(LEN=64)    :: file_pos='0', file_vel='0'
INTEGER              :: nb_atm, nb_step
CHARACTER(LEN=2)     :: suffix
REAL(dp)             :: xlo=0.0_dp,xhi=0.0_dp
REAL(dp)             :: ylo=0.0_dp,yhi=0.0_dp
REAL(dp)             :: zlo=0.0_dp,zhi=0.0_dp
REAL(dp)             :: rOH_cut, rOC_cut
REAL(dp)             :: box(3)
INTEGER              :: nb_Cgo
INTEGER              :: waterlist

CONTAINS

SUBROUTINE READINPUTSUB(input_file)

    CHARACTER(LEN=100) :: input_file

    OPEN(file_unit, file=input_file)

    RL:DO WHILE (iostatus == 0)
        READ(file_unit, '(A)', IOSTAT=iostatus) value
        IF (iostatus == 0) then
            line_c = line_c + 1

            delim = SCAN(value, '    ')
            label = value(1:delim)
            value = value(delim+1:)
            IF (label(1:1) == '!' ) THEN
                CYCLE RL
            END IF
            SELECT CASE (label)
                CASE ('file_pos')
                    READ(value, *, IOSTAT=iostatus) file_pos
                    PRINT *, 'XYZ positions file:', file_pos
                CASE ('file_vel')
                    READ(value, *, IOSTAT=iostatus) file_vel
                    PRINT *, 'XYZ velocities file:', file_vel
                CASE ('nb_atm')
                    READ(value, *, IOSTAT=iostatus) nb_atm
                    PRINT *, 'Number of Atoms: ', nb_atm
                CASE ('nb_step')
                    READ(value, *, IOSTAT=iostatus) nb_step
                    PRINT *, 'nb_step: ', nb_step
                CASE ('suffix')
                    READ(value, *, IOSTAT=iostatus) suffix
                    PRINT *, 'suffix: ', suffix
                CASE ('nb_Cgo')
                    READ(value, *, IOSTAT=iostatus) nb_cgo
                    PRINT *, 'nb_Cgo: ', nb_cgo
                CASE ('xlo')
                    READ(value, *, IOSTAT=iostatus) xlo
                    PRINT *, 'xlo: ', xlo
                CASE ('xhi')
                    READ(value, *, IOSTAT=iostatus) xhi
                    PRINT *, 'xhi: ', xhi
                CASE ('ylo')
                    READ(value, *, IOSTAT=iostatus) ylo
                    PRINT *, 'ylo: ', ylo
                CASE ('yhi')
                    READ(value, *, IOSTAT=iostatus) yhi
                    PRINT *, 'yhi: ', yhi
                CASE ('zlo')
                    READ(value, *, IOSTAT=iostatus) zlo
                    PRINT *, 'zlo: ', zlo
                CASE ('zhi')
                    READ(value, *, IOSTAT=iostatus) zhi
                    PRINT *, 'zhi: ', zhi
                CASE ('rOH_cut')
                    READ(value, *, IOSTAT=iostatus) rOH_cut
                    PRINT *, 'rOH_cut: ', rOH_cut
                CASE ('rOC_cut')
                    READ(value, *, IOSTAT=iostatus) rOC_cut
                    PRINT *, 'rOC_cut: ', rOC_cut
                CASE ('waterlist')
                    READ(value, *, IOSTAT=iostatus) waterlist
                    PRINT *, 'waterlist: ', waterlist
                CASE DEFAULT
                    PRINT *, 'Invalid label, line:', line_c
            END SELECT
        END IF
    END DO RL

    ! ----------------------------------------------- Calculate the box size
    box(1) = xhi - xlo
    box(2) = yhi - ylo
    box(3) = zhi - zlo

END SUBROUTINE READINPUTSUB

END MODULE INPUT_MOD