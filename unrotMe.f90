PROGRAM unrotMe
  IMPLICIT NONE
  DOUBLE PRECISION :: lon_rot_in, lat_rot_in, lon_pol_in, lat_pol_in
  DOUBLE PRECISION :: phirot2phi, rlarot2rla
  DOUBLE PRECISION :: phirot2phi_out,   rlarot2rla_out
  INTEGER          :: phirot2phi_deg,   rlarot2rla_deg
  DOUBLE PRECISION :: phirot2phi_min,   rlarot2rla_min

  CALL GET_COORDS(lat_rot_in,lon_rot_in,lat_pol_in,lon_pol_in)

  lat_pol_in = 89.208767
  lon_pol_in = 179.050450D0
  phirot2phi_out = phirot2phi(lat_rot_in,lon_rot_in,lat_pol_in,lon_pol_in,0.0D0)
  rlarot2rla_out = rlarot2rla(lat_rot_in,lon_rot_in,lat_pol_in,lon_pol_in,0.0D0)
  phirot2phi_deg = INT(phirot2phi_out)
  rlarot2rla_deg = INT(rlarot2rla_out)
  phirot2phi_min = (phirot2phi_out - DBLE(INT(phirot2phi_out)))  * 60.0D0
  rlarot2rla_min = (rlarot2rla_out - DBLE(INT(rlarot2rla_out)))  * 60.0D0

  IF ( phirot2phi_min >= 59.9995D0 )THEN
     phirot2phi_deg = phirot2phi_deg + 1.0D0
     phirot2phi_min = 0.D0
  END IF
  IF ( rlarot2rla_min >= 59.9995D0 )THEN
     rlarot2rla_deg = rlarot2rla_deg + 1.0D0
     rlarot2rla_min = 0.D0
  END IF

  !  WRITE(*,'(A2,F13.10,A2,F13.10)') &
  !       ' N',phirot2phi_out,  &
  !       ' E',rlarot2rla_out
  WRITE(*,*)
  WRITE(*,*) 'The given point in the rotated system is located at:'
  WRITE(*,"(' N',I2,'°',F6.3,A1,' E',I3.3,'°',F6.3,A1)") &
       phirot2phi_deg, phirot2phi_min, "'", &
       rlarot2rla_deg, rlarot2rla_min, "'"
  WRITE(*,*) 'in the unrotated geographic system!'
  WRITE(*,*) 'Happy hunting!'

END PROGRAM unrotMe

SUBROUTINE GET_COORDS(lat,lon,pollat,pollon)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(OUT) :: lat,lon,pollat,pollon
  DOUBLE PRECISION             :: lat_deg, lat_min, lon_deg, lon_min
  DOUBLE PRECISION             :: pol_lat_deg, pol_lat_min, pol_lon_deg, pol_lon_min
  INTEGER :: io_err

  WRITE(*,*)
  WRITE(*,*) '------------------------'
  WRITE(*,*) 'Input format : DD°MM.MMM'
  WRITE(*,*) '------------------------'
  WRITE(*,*)

  WRITE(*,'(A)')'Enter LATITUDE of the point in ROTATED system'
  DO
     WRITE(*,'(A$)')'Degrees     (DD) : '
     READ (*,*,iostat=io_err) lat_deg
     if ( io_err == 0 ) then
        exit
     else
        cycle
     end if
  END DO
  DO
     WRITE(*,'(A$)')'Minutes (MM.MMM) : '
     READ (*,*,iostat=io_err) lat_min
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  WRITE(*,'(A)')'Enter LONGIUDE of the point in ROTATED system'
  DO
     WRITE(*,'(A$)')'Degrees     (DD) : '
     READ (*,*,iostat=io_err) lon_deg
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  DO
     WRITE(*,'(A$)')'Minutes (MM.MMM) : '
     READ (*,*,iostat=io_err) lon_min
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  lat = lat_deg + (lat_min / 60.0D0 )
  lon = lon_deg + (lon_min / 60.0D0 )

  WRITE(*,'(A)')'Enter LATITUDE of the rotated North Pole in UNROTATED system'
  DO
     WRITE(*,'(A$)')'Degrees     (DD) : '
     READ (*,*,iostat=io_err) pol_lat_deg
     if ( io_err == 0 ) then
        exit
     else
        cycle
     end if
  END DO
  DO
     WRITE(*,'(A$)')'Minutes (MM.MMM) : '
     READ (*,*,iostat=io_err) pol_lat_min
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  WRITE(*,'(A)')'Enter LONGIUDE of the rotated North Pole in UNROTATED system'
  DO
     WRITE(*,'(A$)')'Degrees     (DD) : '
     READ (*,*,iostat=io_err) pol_lon_deg
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  DO
     WRITE(*,'(A$)')'Minutes (MM.MMM) : '
     READ (*,*,iostat=io_err) pol_lon_min
     if ( io_err == 0) then
        exit
     else
        cycle
     end if
  END DO
  pollat = pol_lat_deg + (pol_lat_min / 60.0D0 )
  pollon = pol_lon_deg + (pol_lon_min / 60.0D0 )

  RETURN
END SUBROUTINE GET_COORDS
DOUBLE PRECISION FUNCTION  phirot2phi ( phirot, rlarot, polphi, pollam, polgam )

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This function converts phi from one rotated system to phi in another
  !   system. If the optional argument polgam is present, the other system
  !   can also be a rotated one, where polgam is the angle between the two
  !   north poles.
  !   If polgam is not present, the other system is the real geographical
  !   system.
  !
  ! Method:
  !   Transformation formulas for converting between these two systems.
  !
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !
  ! Declarations:
  !
  !------------------------------------------------------------------------------

  ! Parameter list:
  DOUBLE PRECISION , INTENT (IN)      ::        &
       polphi,   & ! latitude of the rotated north pole
       pollam,   & ! longitude of the rotated north pole
       phirot,   & ! latitude in the rotated system
       rlarot      ! longitude in the rotated system

  DOUBLE PRECISION , INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  DOUBLE PRECISION                    ::        &
       phirot2phi  ! latitude in the geographical system

  DOUBLE PRECISION                    ::        &
       zsinpol, zcospol, zphis, zrlas, zarg, zgam

  DOUBLE PRECISION        ::        &
       zrpi18, zpir18

  ! Begin function phirot2phi

  pi = 2.0D0 * ACOS(0.0D0)
  zrpi18 = 180.0D0 / pi
  zpir18 = pi / 180.0D0

  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)

  zphis       = zpir18 * phirot
  IF (rlarot > 180.0D0) THEN
     zrlas = rlarot - 360.0D0
  ELSE
     zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas

  IF (polgam /= 0.0D0) THEN
     zgam  = zpir18 * polgam
     zarg  = zsinpol*SIN (zphis) +                                           &
          zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
     zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF

  phirot2phi  = zrpi18 * ASIN (zarg)

  RETURN
END FUNCTION phirot2phi

FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This function converts lambda from one rotated system to lambda in another
  !   system. If the optional argument polgam is present, the other system
  !   can also be a rotated one, where polgam is the angle between the two
  !   north poles.
  !   If polgam is not present, the other system is the real geographical
  !   system.
  !
  ! Method:
  !   Transformation formulas for converting between these two systems.
  !
  ! Modules used:    NONE
  !
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !
  ! Declarations:
  !
  !------------------------------------------------------------------------------

  ! Parameter list:
  DOUBLE PRECISION , INTENT (IN)      ::        &
       polphi,   & ! latitude of the rotated north pole
       pollam,   & ! longitude of the rotated north pole
       phirot,   & ! latitude in the rotated system
       rlarot      ! longitude in the rotated system

  DOUBLE PRECISION , INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  DOUBLE PRECISION                    ::        &
       rlarot2rla  ! latitude in the geographical system

  ! Local variables
  DOUBLE PRECISION                    ::        &
       zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

  DOUBLE PRECISION        ::        &
       zrpi18, zpir18

  !------------------------------------------------------------------------------

  ! Begin function rlarot2rla

  pi = 2.0D0 * ACOS(0.0D0)
  zrpi18 = 180.0D0 / pi
  zpir18 = pi / 180.0D0

  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)

  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0D0) THEN
     zrlas = rlarot - 360.0D0
  ELSE
     zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas

  IF (polgam /= 0.0) THEN
     zgam    = zpir18 * polgam
     zarg1   = SIN (zlampol) *                                                &
          (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
          + zcospol * SIN(zphis))                                               &
          - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
     zarg2   = COS (zlampol) *                                                &
          (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
          + zcospol * SIN(zphis))                                               &
          + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
  ELSE
     zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
          zcospol *              SIN(zphis)) -    &
          COS (zlampol) *             SIN(zrlas) * COS(zphis)
     zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
          zcospol *              SIN(zphis)) +   &
          SIN (zlampol) *             SIN(zrlas) * COS(zphis)
  ENDIF

  IF (zarg2 == 0.0) zarg2 = 1.0D-20

  rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)

  RETURN
END FUNCTION rlarot2rla
