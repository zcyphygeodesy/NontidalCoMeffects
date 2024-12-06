## Fortran codes for nontidal geocentric variation effects of mass on all-element geodetic variations
https://www.zcyphygeodesy.com/en/h-nd-123.html
## [Algorithm purpose]
    Given the longitude, latitude, ellipsoidal height and time of the calculation point, using the Earth's mass centric variation time series from the measured SLR, compute the Earth's centric variation effects of mass on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (10μE) or horizontal gravity gradient (NW, to the north and to the west, 10μE).
    The variation of the Earth's center of mass is equal to the first-degree term of Earth's loading deformation, which excites the variations of all the geometric and physical geodetic elements in the Earth's space with time, rather than can be simply expressed as the ground site displacement of pure geometric quantity.
    Improve the algorithm of Earth's centric variation effects of mass in the IERS Conventions (2010) and then compute the nontidal geocentric variation effects of mass on all-element geodetic variations in the whole Earth space.
## [Computation Output]
    tdn(14): the nontidal geocentric variation effects of mass on all-element geodetic variations.
    tdn(1:14) stores the nontidal geocentric variation effects of mass on 10 kinds of geodetic variations, which are the solid tidal effects on height anomaly tdn(1) (mm), ground gravity #tdn(2) (μGal), gravity disturbance tdn(3) (μGal), ground tilt #tdn(4:5) (SW, to the south and to the west, mas), vertical deflection tdn(6:7) (SW, to the south and to the west, mas), horizontal displacement #tdn(8:9) (EN, to the east and to the north, mm), ground radial displacement #tdn(10) (mm), ground normal or orthometric height #tdn(11) (mm), radial gravity gradient tdn(12 )(10μE) or horizontal gravity gradient tdn(13:14) (NW, to the north and to the west, 10μE).
    The calculation point can be on the ground, low altitude, satellite, ocean or underwater space. The geodetic variations abvove marked with # are valid only when the site is fixed with the solid Earth.
## [Geophysical models]
    (1) the Earth's mass centric variation time series from the measured SLR.
    The Earth's mass centric variation time series file Monthly_geocenter_MK.txt (ITRF2014) from Center for Space Research in University of Texas in USA (UT/CSR) from LAGEOS-1/2, Stella, Starlette, AJISAI, BEC and LARES Satellite Laser Ranging (SLR) measured.
## [Main program for test entrance]
    NontidalCoMeffects.f90
    The record of the test output file reslt.txt: the long integer time agreed by ETideLoad, difference between the MJD day and starting MJD0, Earth's mass centric variation dx,dy,dz (mm), tdn(1:14)
## (1) Calculation module for normal Earth’s gravity field
    GNormalfd(BLH, NFD, GRS)
    Input parameters: BLH(3) - latitude, longitude (decimal degrees) and ellipsoidal height (m) at the calculation point.
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (2) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (3) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (4) Algorithm library for converting of time system
    CAL2JD (IY0, IM0, ID0, DJM, J); JD2CAL(DJ1, DJ2, IY, IM, ID, FD, J)
## (5) Other auxiliary modules
    PickRecord(str0, kln, rec, nn); tmcnt(tm, iyr, imo, idy, ihr, imn, sec)
    mjdtotm(mjd0, ltm); tmtostr(tm, tmstr)
    spline1d(x,y,xa,ya,n); spline3( x, y, n, sx, f, f1, m)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] ETideLoad4.5 User Reference https://www.zcyphygeodesy.com/en/
    8.5.3 Earth's mass centric variation effects on all-element geodetic variations

DOS executable test file, geophysical models and all input and output data.
