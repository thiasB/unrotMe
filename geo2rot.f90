PROGRAM geo2rot
  IMPLICIT NONE
  DOUBLE PRECISION :: lon_geo_in, lat_geo_in,  lon_pol_in,  lat_pol_in
  DOUBLE PRECISION :: phi2phirot, rla2rlarot
  DOUBLE PRECISION :: phi2phirot_out, rla2rlarot_out
  INTEGER          :: phi2phirot_deg, rla2rlarot_deg
  DOUBLE PRECISION :: phi2phirot_min, rla2rlarot_min

  lat_geo_in = 52.16666D0
  lon_geo_in = 13.33333D0
  lat_pol_in  = 89.208767D0
  lon_pol_in  = 179.050450D0

  phi2phirot_out = phi2phirot(lat_geo_in,lon_geo_in,lat_pol_in,lon_pol_in)
  rla2rlarot_out = rla2rlarot(lat_geo_in,lon_geo_in,lat_pol_in,lon_pol_in,0.0D0)
  phi2phirot_deg = INT(phi2phirot_out)
  rla2rlarot_deg = INT(rla2rlarot_out)
  phi2phirot_min = (phi2phirot_out - DBLE(INT(phi2phirot_out)))  * 60.0D0
  rla2rlarot_min = (rla2rlarot_out - DBLE(INT(rla2rlarot_out)))  * 60.0D0
  WRITE(*,"(' N',I2,'°',F6.3,A1,' E',I3.3,'°',F6.3,A1)") &
       phi2phirot_deg, phi2phirot_min, "'", &
       rla2rlarot_deg, rla2rlarot_min, "'"

  WRITE(*,'(A2,f9.5,A2,f9.5)') &
       ' N', phi2phirot(lat_geo_in,lon_geo_in,lat_pol_in,lon_pol_in), &
       ' E', rla2rlarot(lat_geo_in,lon_geo_in,lat_pol_in,lon_pol_in,0.0D0)

END PROGRAM geo2rot

DOUBLE PRECISION FUNCTION  phi2phirot ( phi, rla, polphi, pollam )
  ! Parameter list:
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (IN)      ::        &
       polphi,  & ! latitude of the rotated north pole
       pollam,  & ! longitude of the rotated north pole
       phi,     & ! latitude in the geographical system
       rla        ! longitude in the geographical system

  DOUBLE PRECISION                   ::        &
       phi2phirot ! longitude in the rotated system

  ! Local variables
  DOUBLE PRECISION                       ::    &
       zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1, pi

  DOUBLE PRECISION            ::    &
       zrpi18, zpir18

  ! Begin function phi2phirot

  pi = 2.0D0 * ACOS(0.0D0)
  zrpi18 = 180.0D0 / pi
  zpir18 = pi / 180.0D0
  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0D0) THEN
     zrla1  = rla - 360.0D0
  ELSE
     zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)
  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

  RETURN
END FUNCTION phi2phirot

DOUBLE PRECISION FUNCTION  rla2rlarot ( phi, rla, polphi, pollam, polgam )
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (IN)      ::        &
       polphi,  & ! latitude of the rotated north pole
       pollam,  & ! longitude of the rotated north pole
       phi,     & ! latitude in geographical system
       rla        ! longitude in geographical system

  DOUBLE PRECISION, INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  DOUBLE PRECISION                   ::        &
       rla2rlarot ! latitude in the the rotated system

  ! Local variables
  DOUBLE PRECISION                   ::    &
       zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1, pi

  DOUBLE PRECISION            ::    &
       zrpi18, zpir18

  ! Begin function rla2rlarot

  pi = 2.0D0 * ACOS(0.0D0)
  zrpi18 = 180.0D0 / pi
  zpir18 = pi / 180.0D0

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0D0) THEN
     zrla1  = rla - 360.0D0
  ELSE
     zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0) zarg2 = 1.0E-30

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

  IF (polgam /= 0.0D0 ) THEN
     rla2rlarot = polgam + rla2rlarot
     IF (rla2rlarot > 180.0D0) rla2rlarot = rla2rlarot -360.0D0
  ENDIF
  RETURN
END FUNCTION rla2rlarot
