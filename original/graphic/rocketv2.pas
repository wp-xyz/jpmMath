 {*************************************************************************
  *             SIMULATION OF THE TRAJECTORY OF A V2 ROCKET                *
  * ---------------------------------------------------------------------- *
  *   From a DOS program By J-P Dumont, adapted to Windows By J-P Moreau   *
  *      (WINCRT version with type Longarray and printing capability)      *
  *                                                                        *
  *                                         Version 2.1, August 1996       *
  *                                             (www.jpmoreau.fr)          *
  * ---------------------------------------------------------------------- *
  * TUTORIAL:                                                              *
  * Geographic longitudes are Westward positive. They can vary between     *
  * -180 and +180 degrees. Latitudes are positive in northern hemisphere.  *
  * They can vary between -90° and +90°.                                   *
  *                                                                        *
  * Hypothesis: The rotating Earth is considered as a flattened ellipsoid. *
  * The rocket V2 is launched vertically from a platform. After some       *
  * seconds, there is a "swing" of trajectory to reach a maximal range.    *
  * The vehicle incidence constantly equals 0. The atmospheric model is    *
  * the ARDC model, 1962. The propulsion thrust is stopped when the speed  *
  * allows to reach the target.                                            *
  *                                                                        *
  * Version 2.1:     Zoom at impact zone,                                  *
  *                  Graph of dynamic pressure,                            *
  *                  Initial cape corrected from Coriolis effect,          *
  *                  Possibility of direct (slope = 45°) or diving         *
  *                  trajectories (slope <> 45°).                          *
  *************************************************************************}
  PROGRAM RocketV2;
  {=========================================================================}
  USES WinCrtMy,WinDos,WinTypes,WObjects,WinProcs,Strings,Time,LargArra,
  Type_def,Graph_2D,WinPrint;  {units used by program}
  {=========================================================================}
  {Maximal number of numeric values stored in Table of type (object) LongArray}
  CONST TablSize = 8000;     {16000 for DtSave = 1 pt/sec}
  blue = $00FF0000;
  red  = $000000FF;
  line = ' =============================================================';

  TYPE
  LongPtr = ^REAL_AR;             {Cf. unit Type_def}
  LongARRAY = OBJECT( LgARRAY )   {Cf. unit LargArra}
    CONSTRUCTOR Init(iNum : LongInt);
    FUNCTION Index(N : LongInt) : LongPtr;
  END;

  {Methods of object LongArray}
  CONSTRUCTOR LongARRAY.Init(iNum : LongInt);
  BEGIN
    IF NOT LgARRAY.Init(iNum, SizeOF(REAL_AR)) THEN
    Fail;
  END;

  FUNCTION LongARRAY.Index(N : LongInt) : LongPtr;
  BEGIN
    Index := LongPtr(LgARRAY.Index(N));
  END;

  {global variables}
  VAR
  OUT                     : TEXT;              {ouput file (screen or printer}
  IprntSw                 : BOOLEAN;               {TRUE=printer,FALSE=screen}
  Table                   : LongARRAY;              {to store numeric results}
  x_y                     : RV;        {pointer to auxiliary table for graphs}
  CrtPen,OldPen,InvPen    : HPen;                                {Pen handles}
  OldTextCoul             : Longint;                      {default text color}
  Pinfo                   : PPrinterInfo;           {pointer to printer infos}
  {==========================================================================}
  CONST

  {Number of simultaneous equations to integrate}
  NbSMb = 6;
  {Number of general parameters to store}
  nVAR  =  8;
  {menu index} 
  trajV =  9;
  TrajH = 10;
  Zoom  = 11;
  Fini  = 12;
  {Maximal dimension of state vector}
  NVAR_Max = 6;
  {Maximal number of different items in Table}
  Maxc = 20;

  {half main axis of geoid}
  a = 6.378140E6;
  {Flatness of geoid = 1/298.257}
  e = 3.352813E-3;
  {gravitationnal constant by Earth mass}
  GM = 3.98603E14;
  {gravitationnal constant by x Difference between main inertial moments}
  GA_B = 17.536E24;
  {Sideral angular Earth speed}
  Omega = 72.9212E-6;
  {atmospheric model: 1- exponential 2- ARDC}
  Model = 2;
  {Altitude scaling for Earth atmosphere}
  Hauteur = 9870;
  {Diameter of missile V2}
  D = 1.65;
  {Air specific mass at ground and at 20° Celsius}
  Rho = 1.2;
  {Drag coefficient}
  Cx = 0.15;
  {Mass of structures}
  Ms = 3870;  {4170 according to french Professor Lacaze}
  {Mass of ergols (Ethanol and liquid oxygen) }
  Me = 8670;  {JONES data: 8700 = 4900 LOX + 3800 Ethanol}
  {Thrust in Newton}
  Thrust = 265E3;  {245.25E3 JONES}
  {Duration in sec. of propulsion phase}
  Tprop = 65;
  {trajectory slope in degrees}
  Pente = 45;
  {Swinging time of trajectory}
  Tbascul = 4;
  {
  Launching point latitude - Rotterdam: 51ø55'19" 
  Phi1 = 51.9219444;
  Launching point longitude - Rotterdam: -4ø29'23"
  Lambda1 = -4.489722;
  }
  {Launching point latitude - Wassenaar - }
  Phi1 = 52.05669;
  {Launching point Longitude - Wassenar - }
  Lambda1 = -4.4096;
  {Target point latitude London St Paul's Cathedral: 51ø30'54" }
  Phi2 = 51.515;
  {Target point longitude London St Paul's Cathedral: 0ø06'21" }
  Lambda2 = 0.105833333;
  cc = 101300.0/760.0;

  Half = 0.5;

  Variable : ARRAY[1..NVAR] OF STRING[30] =
  ('Distance in longitude [km]',
  'Distance in latitude [km]',
  'Altitude [km]',
  'Speed  [m/s]',
  'Dynamic pressure [MPa]',
  'Distance to Target [km]',
  'Nozzle inclination [°]',
  'Azimut of Target [°]');

  {------------------------------------------------------------------------}
  TYPE Vecteurs_Etat = ARRAY[1..NVAR_Max] OF REAL_AR;
  Vectors     = ARRAY[1..3] OF REAL_AR;
  Pointers    = ARRAY[0..Maxc] OF WORD;

  {------------------------------------------------------------------------}
  VAR Alpha,Altitude,Azimut,beta,
  Cazim,Sazim,
  Deg,Distance,dt,dtmIN,dtsave,
  eps,ergols,gamma,gr,gphi,h1,
  MaitreCouple,Masse,Omega2,Pdyn,Q,Rad,
  Lambda,Phi,Lambdab,Phib,
  Sphi,Cphi,Slambda,Clambda,R,R2,
  TgPente,Temps,Textinct,t2,Tmax,
  V,V2,X,Y,Z                       : REAL_AR;
  i,istep,nbad,nok,nsmb,
  nstep,NstepMax                   : INTEGER;
  EpuisementErgols                 : BOOLEAN;
  Etat                             : Vecteurs_Etat;
  ipt                              : Pointers;
  {==========================================================================}
  FUNCTION  Tan(x: REAL_AR) : REAL_AR;
  VAR Cx : REAL_AR;
  BEGIN
    Cx:=Cos(x);
    IF (Cx<>0) THEN Tan:=sIN(x)/Cx
    ELSE
    BEGIN
      WRITELN(' Argument of Tan odd multiple of Pi/2');
      Halt;
    END
  END;
  {--------------------------------------------------------------------------}
  PROCEDURE ATan2 {Returns an argument (phase) between -PI et +PI}
  (Numerateur,denominateur:REAL_AR;
  VAR Phase:REAL_AR);
  CONST DEMIPI=PI/2;
  BEGIN
    IF (Abs(Denominateur)<Macheps) THEN
    IF (numerateur < 0) THEN Phase:=-DEMIPI
    ELSE Phase:= DEMIPI
    ELSE
    BEGIN
      Phase:=Arctan(numerateur/denominateur);
      IF (denominateur <0 ) THEN Phase := Phase + Pi;
    END;       {else}
  END;         {ATAN2}
  {=========================================================================}
  PROCEDURE AzimutBut
  (Lambda1,Lambda2,Phi1,Phi2: REAL_AR;
  VAR Azimut : REAL_AR);
  {-------------------------------------------------------------------------
  Calculate the bearing of target in local launch pad coordinates.
  The x axis is the local meridian northward oriented.
  The y axis is the local parallel westward oriented.
  The z axis is the local vertical upwards oriented.
  The bearing is calculated from Ox axis in positive direction around Oz axis
  (angular unit: radians).
  Référence: ITT - REFERENCE DATA FOR RADIO ENGINEERS page 26.10
  --------------------------------------------------------------------------}
  CONST Half = 0.5;
  VAR AlphaM,AlphaP,Den,DiffLat,
  DiffLong,Num,Rad,SumLat,TgCs2    : REAL_AR;

  BEGIN
    DiffLong := Lambda2-Lambda1;
    DiffLat  := Phi2-Phi1;
    SumLat   := Phi1+Phi2;
    TgCs2    := Tan(Half*DiffLong);
    Num := sIN(Half*DiffLat);
    Den := cos(Half*SumLat)*TgCs2;
    Atan2(Num,Den,AlphaM);
    Num := cos(Half*DiffLat);
    Den := sIN(Half*SumLat)*TgCs2;
    Atan2(Num,Den,AlphaP);
    Azimut := AlphaP - AlphaM;
  END;
  {=========================================================================}
  PROCEDURE DistanceBut
  (Lambda1,Lambda2,Phi1,Phi2:REAL_AR;
  VAR Distance:REAL_AR);
  {-------------------------------------------------------------------------
  Calculate the distance along a great Earth circle to reach a target with a
  longitude lambda2 and a latitude Phi2.
  The radius of a great earth circle is taken as the half longitudinal axis
  of terrestial ellipsoid.
  Version 1.1 of 07/24/1992: this version corrects a bug in the case where
  the lauch pad and the target are on the same parallel.
  --------------------------------------------------------------------------}
  VAR AlphaM,AlphaP,Den,DiffLat,Result,
  DiffLong,Num,Rad,SumLat,Tgz,z   : REAL_AR;
  Ntour : LONGINT;

  BEGIN
    DiffLong := Lambda2-Lambda1;
    DiffLat  := Phi2-Phi1;
    SumLat   := Phi1+Phi2;
    Num := sIN(0.5*DiffLat);
    Den := cos(0.5*SumLat)*Tan(0.5*DiffLong);
    Atan2(Num,Den,AlphaM);
    Num := cos(0.5*DiffLat);
    Den := sIN(0.5*SumLat)*Tan(0.5*DiffLong);
    Atan2(Num,Den,AlphaP);

    IF(Abs(DiffLat) > Macheps) THEN
    BEGIN
      Num:= tan(0.5*DiffLat)*sIN(AlphaP);
      Den:= sIN(AlphaM);
      ATAN2(Num,Den,Result);
      z := 2*Result;
    END
    ELSE
    BEGIN  {case Phi1 = Phi2}
      IF (Abs(SumLat)> Macheps) THEN
      BEGIN
        Num:= Cos(AlphaP);
        Den:= Tan(Phi1);
        ATAN2(Num,Den,Result);
        z := 2*Result;
      END
      ELSE z := DiffLong;  {2 points on equator}
    END;
    z := Abs(z);

    IF(z > 2*PI) THEN
    BEGIN
      Ntour:= Trunc(z/(2*PI));
      z := z - Ntour*2*PI;
    END;

    Distance := a*z;
  END;

  {=========================================================================}
  PROCEDURE CoordonneesGeocentriquesInitiales
  (Phi,Lambda  : REAL_AR;
  VAR X,Y,Z : REAL_AR);
  {--------------------------------------------------------------------------
  Calculate the coordinates X,Y,Z, knowing the longitude Lambda and the local
  local geocentric latitude.
  Origine O of coordinates: Earth center.
  Ox axis: intersection of Equator and Greenwich meridian.
  OY axis: eastward, such as OXYZ is direct.
  OZ axis: along Earth axis, towards North pole.
  Local variable: Rayon = distance to Earth center.
  Units : Distances in meters.
  --------------------------------------------------------------------------}
  CONST un = 1;
  VAR Rayon,RCoef : REAL_AR;
  BEGIN
    Cphi := Cos(Phi);
    SPhi := SIN(Phi);
    CLambda := Cos(Lambda);
    SLambda := SIN(Lambda);
    RCoef := un + (un/sqr(un-e)-un)*SPhi*Sphi;
    Rayon := a/Sqrt(Rcoef);
    X :=  Rayon*Cphi*CLambda;
    Y := -Rayon*Cphi*SLambda; (* Sens particulier des longitudes *)
    Z :=  Rayon*SPhi;
  END;

  {=========================================================================}
  FUNCTION CorrectionDeCap(Cap:REAL_AR): REAL_AR;
  {--------------------------------------------------------------------------
  The cape correction is calculated with the hypothesis of a flat Earth and a
  speed at end of propulsion phase on a trajectory with a slope of 45 degrees.
  --------------------------------------------------------------------------}
  {pesanteur = 9.81  Vzero = 1600}
  CONST C=1.70247963E8;    { C = 4*Vzero*sqr(Vzero/sqr(pesanteur) }
  VAR   Cazim,DeltaY,Sgam : REAL_AR;
  BEGIN
    Sgam   := sIN(Gamma);
    Cazim  := Cos(Cap);
    DeltaY := C*Omega;
    DeltaY := DeltaY*(Cphi*CAzim*Sgam/3-Sphi*cos(Gamma))*Sqr(Sgam);
    CorrectionDeCap := DeltaY/Distance;
  END;
  {=========================================================================}
  PROCEDURE ConditionsInitiales
  (Phi,Lambda : REAL_AR;
  VAR Etat: Vecteurs_Etat);
  BEGIN
    CoordonneesGeocentriquesInitiales(Phi,Lambda,X,Y,Z);
    Etat[1] := X;
    Etat[2] := Y;
    Etat[3] := Z;
    Etat[4] := 0;
    Etat[5] := 0;
    Etat[6] := 0;
  END;
  {=========================================================================}
  PROCEDURE Atmosphere_ARDC_1962(z:REAL_AR;VAR zm,Tm,pm,zp,Tp,pp: REAL_AR);
  {--------------------------------------------------------------------------
  Seek in table nearest inferior and superior altitudes w.r.t. current
  altitude and associated values of temperature and pressure.
  Table units: Altitude: 1000 ft - Temperature: °F - Pressure: Torr or mm Hg
  --------------------------------------------------------------------------}
  CONST Ardc_Table : ARRAY[1..74,1..3] OF REAL_AR =
  ((0,     59.0,       760.11),
  (2.5,   50.1,       693.81),
  (5.0,   41.2,       632.38),
  (7.5,   32.3,       575.45),
  (10,    23.4,       522.75),
  (12.5,  14.5,       474.04),
  (15,    5.5,        429.08),
  (17.5,  -3.4,       387.66),
  (20,    -12.3,      349.53),
  (22.5,  -21.2,      314.51),
  (25,    -30,        282.40),
  (27.5,  -38.9,      253.00),
  (30,    -47.8,      226.13),
  (32.5,   -56.7,      201.63),
  (35,    -65.6,      179.33),
  (40,    -69.7,      141.18),
  (45,    -69.7,      111.13),
  (50,    -69.7,       87.49),
  (55,    -69.7,       68.88),
  (60,    -69.7,       54.24),
  (70,    -67.4,       33.66),
  (80,    -61.9,       21.01),
  (90,    -56.5,       13.21),
  (100,   -51.1,       8.36),
  (110,   -41.3,       5.33),
  (120,   -26.1,       3.45),
  (130,   -10.9,       2.27),
  (140,     4.3,       1.51),
  (150,     19.4,      1.02),
  (160,     27.5,      0.697),
  (170,     27.5,      0.478),
  (180,     18.9,      0.326),
  (190,      8.1,      0.221),
  (200,     -2.7,      0.148),
  (210,     -21.9,     9.85E-2),
  (220,     -43.5,     6.41E-2),
  (230,     -64.9,     4.08E-2),
  (240,     -86.4,     2.53E-2),
  (250,     -107.8,    1.53E-2),
  (260,     -129.3,    8.92E-3),
  (270,     -134.5,    5.09E-3),
  (280,     -134.5,    2.91E-3),
  (290,     -134.5,    1.66E-3),
  (300,     -126.8,    9.5E-4),
  (325,     -86.5,     2.6E-4),
  (350,     -24.5,     8.5E-5),
  (375,     75.5,      3.3E-5),
  (400,     234,       1.6E-5),
  (450,     734,       6.3E-6),
  (500,    1204,       3.5E-6),
  (550,    1492,       2.2E-6),
  (600,    1647,       1.5E-6),
  (650,    1755,       1.0E-6),
  (700,    1836,       7.4E-7),
  (750,    1913,       5.4E-7),
  (800,    1964,       4.0E-7),
  (850,    2011,       2.9E-7),
  (900,    2053,       2.2E-7),
  (950,    2093,       1.7E-7),
  (1000,   2125,       1.3E-7),
  (1100,   2160,       7.9E-8),
  (1200,   2189,       5.0E-8),
  (1300,   2215,       3.2E-8),
  (1400,   2217,       2.1E-8),
  (1500,   2221,       1.4E-8),
  (1600,   2232,       9.6E-9),
  (1700,   2234,       6.6E-9),
  (1800,   2233,       4.6E-9),
  (1900,   2241,       3.3E-9),
  (2000,   2251,       2.3E-9),
  (2100,   2252,       1.7E-9),
  (2200,   2254,       1.2E-9),
  (2300,   2254,       8.8E-10),
  (2320,   2253,       8.3E-10));

  VAR i,j : INTEGER;
  BEGIN
    i:= 1;
    j:= 1;
    zp := Ardc_Table[i,j];
    REPEAT
      INC(i);
      zm := zp;
      zp := Ardc_Table[i,j];
    UNTIL zp > z;
    Tm := Ardc_Table[i-1,2];
    Tp := Ardc_Table[i,2];
    pm := Ardc_Table[i-1,3];
    pp := Ardc_Table[i,3];
  END;
  {=========================================================================}
  FUNCTION AirDensity(Altitude:REAL_AR):REAL_AR;
  CONST InvFeet = 3.28084;
  TCONST1 = 255.38;
  TCONST2 = 0.5555555;
  AirMol  = 0.0288;  {Molar mass of 22.4 l of air in kg }
  GasCONST= 8.3169;  {Gaz constant in J} 
  VAR interp,z,zm,zp,Tm,Tp,Temperature,pm,pp,pressure: REAL_AR;
  BEGIN
    CASE Model OF
      1: AirDensity := Rho*exp(-Altitude/Hauteur);
      2: BEGIN
        z:= 1.0E-3*Altitude*InvFeet; (* Conversion de mètres en Kft *)
        Atmosphere_ARDC_1962(z,zm,tm,pm,zp,tp,pp);
        interp := (z-zm)/(zp-zm);
        Temperature:=tm+interp*(tp-tm);  (* degrés Fahrenheit *)
        pressure   :=pm+interp*(pp-pm);  (* torr *)
        Temperature:=TCONST1+TCONST2*Temperature; (* øK *)
        Pressure   :=pressure*cc;        (* cc=101300.0/760.0  (Pa) *)
        AirDensity :=Pressure*AirMol/(GasCONST*Temperature);
      END;
    END; {Case}
  END; {AirDensity}

  {=========================================================================}
  PROCEDURE AirDrag
  (Altitude,Vx,Vy,Vz:REAL_AR;
  VAR Dragx,Dragy,Dragz: REAL_AR);
  {-------------------------------------------------------------------------}
  VAR DragCoef, DrV : REAL_AR;
  BEGIN
    DragCoef := 0.5*AirDensity(Altitude)*MaitreCouple*Cx;
    V2       := Vx*Vx+Vy*Vy+Vz*Vz;
    V        := sqrt(V2); 
    DrV      := -DragCoef*V;
    Dragx    := DrV*Vx;
    Dragy    := DrV*Vy;
    Dragz    := DrV*Vz;
  END;

  {=========================================================================}
  PROCEDURE Gravitation(X,Y,Z: REAL_AR; VAR gx,gy,gz : REAL_AR);
  {-------------------------------------------------------------------------
  Calculate the geocentric gravity components with the hypothesis that the
  geoid is a revolution body defined by its mass and the difference (A-B)
  of its main inertial moments. The GM and G(A-B) values are given by G
  measurements at several latitudes.
  G = Universal Gravity Constant
  M = Earth Mass
  -------------------------------------------------------------------------}
  CONST Trois = 3;
  un    = 1;
  Half  = 0.5;
  VAR   GABSR4  : REAL_AR;

  BEGIN
    GABSR4 := Trois*GA_B/(R2*R2);
    gr := -GM/R2 - Half*GABSR4*(un-Trois*SPhi*SPhi);
    gphi := -SPhi*CPhi*GABSR4;
    gx := gr*CPhi-gphi*SPhi;
    gy := -gx*SLambda;
    gx := gx*CLambda;
    gz := gr*SPhi+gPhi*CPhi;
  END;
  {=========================================================================}
  PROCEDURE NewGeographicCoordinates
  (X,Y,Z:REAL_AR;
  VAR Phi,Lambda:REAL_AR);

  VAR Rp,Rp2 : REAL_AR;

  BEGIN
    Rp2:= X*X+Y*Y;
    R2 := Rp2+Z*Z;
    Rp := Sqrt(Rp2);
    Atan2(Z,Rp,Phi);
    SPhi := SIN(Phi);
    CPhi := Cos(Phi);
    Atan2(Y,X,Lambda);
    Lambda := - Lambda;        {sign convention for longitudes}
    SLambda := SIN(Lambda);
    CLambda := Cos(Lambda);
  END;
  {=========================================================================}
  PROCEDURE NewAltitude
  (X,Y,Z:REAL_AR;
  VAR altitude:REAL_AR);
  CONST un = 1;
  VAR Rayon,RCoef : REAL_AR;
  BEGIN
    RCoef := un + (un/sqr(un-e)-un)*SPhi*Sphi;
    Rayon := a/Sqrt(Rcoef);
    R2:=X*X+Y*Y+Z*Z;
    R:=sqrt(R2);
    altitude := R-Rayon
  END;
  {=========================================================================}
  FUNCTION MasseInstantanee(t: REAL_AR):REAL_AR;
  BEGIN
    ergols := Me -Q*t;
    IF (ergols <= 0) THEN
    IF (EpuisementErgols) THEN ergols := 0
    ELSE
    BEGIN
      ergols := 0;
      EpuisementErgols := true;
      Textinct := t;
      GoToXY(1,7);
      WRITELN(line);
      WRITELN(' Ergols exhausted at t = ',Textinct:6:3,' seconds');
      DistanceBut(Lambda,Lambdab,Phi,Phib,Distance);
      WRITELN(' Distance to target    = ',1e-3*distance:6:3,' [km]');
      WRITELN(line);
    END;
    MasseInstantanee := Ms + ergols;
  END;

  {=========================================================================}
  PROCEDURE ThrustVectors(t,X,Y,Z: REAL_AR;
  VAR Px,Py,Pz: REAL_AR);
  {-------------------------------------------------------------------------
  Nozzle orientation is controled in order to have a thrust angle of 45°
  with the local horizontal plane.
  Thrust is in the plane of the great Earth circle joining the launch pad
  to the target.
  The beta calculation takes into account the missile instantaneous mass.
  Thrust is supposed to pass through the misssile inertial center.
  --------------------------------------------------------------------------}
  CONST un = 1;

  VAR den,eta,num,
  Sbeta,Cbeta,
  CazimCbeta,SazimCbeta,
  CPhiSbeta,SPhiCbeta,temp  : REAL_AR;

  BEGIN
    IF (t >= Tbascul) THEN
    BEGIN
      eta := Abs(Masse*gr)/Thrust;
      num := -un+sqrt(un+Sqr(TgPente)-eta*eta);
      den :=  TgPente-eta;
      Atan2(Num,Den,beta);
      beta := 2*beta;      {nozzle inclination on local horizontal} 
      Sbeta := SIN(beta);
      Cbeta := Cos(beta);
      CAzimCbeta := CAzim*Cbeta;
      SAzimCbeta := SAzim*Cbeta;
      CPhiSbeta := CPhi*Sbeta;
      SphiCbeta := Sphi*Cbeta;
      temp := CPhiSbeta-SPhi*CAzimCbeta;
      Px:=Thrust*( CLambda*temp-SLambda*SAzimCbeta);
      Py:=Thrust*(-SLambda*temp-CLambda*SAzimCbeta);
      Pz:=Thrust*( CPhi*CAzimCbeta+SPhi*Sbeta);
    END
    ELSE
    BEGIN
      beta := Rad*90;    {At first, the nozzle is aligned with missile axis}
      Px:= Thrust*CPhi*CLambda;
      Py:=-Thrust*CPhi*SLambda;
      Pz:= Thrust*SPhi;
    END;
  END;   {ThrustVectors}


  {===================================================================}
  PROCEDURE Derivs (t: REAL_AR; XT: vecteurs_etat; VAR XF: vecteurs_etat);
  {==================================================================
  This procedure calculates the second members of the 1st order
  differential system:  dy/dt = F(y,t)
  ------------------------------------------------------------------
  INPUTS:
    t      = value of instantaneous speed at beginning of time step
    XT     = State vectors at beginning of integration step
  -------------------------------------------------------------------
  OUTPUT:
    XF     = second member vectors F(y,t)
  -------------------------------------------------------------------}

  CONST Two  = 2;

  VAR Dragx,Dragy,Dragz,
  gx,gy,gz,Px,Py,Pz,Vx,Vy,Vz          : REAL_AR;
  XV2                                 : Vectors;

  BEGIN
    X     := XT[1];        (* Etat du système *)
    Y     := XT[2];
    Z     := XT[3];
    Vx    := XT[4];
    Vy    := XT[5];
    Vz    := XT[6];

    NewGeographicCoordinates(X,Y,Z,Phi,Lambda);
    NewAltitude(X,Y,Z,Altitude);   (* calcule aussi R et R2 *)
    AirDrag(Altitude,Vx,Vy,Vz,Dragx,Dragy,Dragz);
    Gravitation(X,Y,Z,gx,gy,gz);
    Masse :=  MasseInstantanee(t);
    IF(NOT EpuisementErgols) THEN ThrustVectors(t,X,Y,Z,Px,Py,Pz)
    ELSE
    BEGIN
      Px := 0;
      Py := 0;
      Pz := 0;
    END;
    (* Calcul des seconds membres *)
    XF[1]:= Vx;                                   (* dX/dt  *)
    XF[2]:= Vy;                                   (* dY/dt  *)
    XF[3]:= Vz;                                   (* dZ/dt  *)
    XF[4]:= Px/Masse+gx+Dragx/Masse+Two*Omega*Vy+Sqr(Omega)*X;
    XF[5]:= Py/Masse+gy+Dragy/Masse-Two*Omega*Vx+Sqr(Omega)*Y;
    XF[6]:= Pz/Masse+gz+Dragz/Masse;
  END;

  {==================================================================}
  PROCEDURE rk4(VAR y,dydt : Vecteurs_Etat;
  n: INTEGER;
  t,h: REAL_AR;
  VAR yout: Vecteurs_Etat);
  {==================================================================
  | Adapted from W.H.PRESS,B.P.FLANNERY,S.A.TEUKOLSKY,W.T.VETTERLING:
  | Numerical Recipes - Cambridge University Press - 1988 -
  ------------------------------------------------------------------
  | This procedure integrates the 1st order differential system:
  |      dy/dt = F(y,t)
  | by a Runge-Kutta method of fourth order to advance the solution
  | on h interval of independant variable t
  ------------------------------------------------------------------
  | INPUTS:
  | y     = State vector at beginning of integration step
  | dydt  = its derivate at same point
  | n     = number of equations of the system
  | t     = value of instantaneous speed at beginning of time step
  | h     = time step
  -------------------------------------------------------------------
  | OUTPUT:
  | yout  =  State vector at end of integration step.
  -------------------------------------------------------------------}
  VAR    i   : INTEGER;
  th,hh,h6   : REAL_AR;
  dym,dyt,yt : Vecteurs_etat;

  BEGIN
    hh := Half* h;
    h6 :=h/6;
    th := t+ hh;
    FOR i := 1 TO n DO yt[i] := y[i] + hh * dydt[i];
    derivs(th,yt,dyt);
    FOR i := 1 TO n DO yt[i] := y[i] + hh * dyt[i];
    derivs(th,yt,dym);
    FOR i := 1 TO n DO BEGIN
      yt[i] := y[i]+ h*dym[i];
      dym[i]:= dyt[i] + dym[i];
    END;
    derivs(t+h,yt,dyt);
    FOR i := 1 TO n DO yout[i] := y[i]+h6*(dydt[i]+dyt[i]+2*dym[i]);
  END; { RK4 }

  {======================================================================}
  PROCEDURE rkqc(VAR y,dydt: Vecteurs_etat;
  n: INTEGER;
  VAR t: REAL_AR;
  htry,eps: REAL_AR;
  yscal: Vecteurs_etat;
  VAR hdid,hnext: REAL_AR);
  {=====================================================================
  | Adapted from W.H.PRESS,B.P.FLANNERY,S.A.TEUKOLSKY,W.T.VETTERLING:
  | Numerical Recipes - Cambridge University Press - 1988 -
  ----------------------------------------------------------------------
  | Runge-Kutta integration step with local error control of truncation
  | and adjustment of integration step to respect a required precision
  ----------------------------------------------------------------------
  | INPUTS:
  | y      = State Vector of dimension n
  | dydt   = its derivative at begin value of independant variable, t
  | XJ     = main inertial vectors.
  | XV     = speed vectors of CoG in fixed coordinates.
  | n      = number of equations of the  system.
  | t      = begin value of the independant variable.
  | htry   = time step proposed as a try.
  | eps    = precision requirement:
  |          Max (ycalc[i] - yvrai[i])/yscal[i] < eps
  | yscal  = normalization vector of solution .
  ----------------------------------------------------------------------
  | OUTPUTS:
  |  y     = final state vecto
  |  t     = end value of instantaneous speed
  |  hdid  = actual time step (increase of i.s.)
  |  hnext = advised time step for next step.
  =====================================================================}

  LABEL 1;

  CONST
  pgrow=-0.20;
  pshrnk=-0.25;
  fcOR=0.06666666;   {1/15 }
  un = 1.0;
  safety=0.9;
  errcon=6.0E-4;
  tiny= Macheps;

  VAR
  i,kount              : INTEGER;
  tsav,hh,h,temp,errmax: REAL_AR;
  dysav,ysav,ytemp     : Vecteurs_etat;

  BEGIN
    tsav:= t;          {Save begin time}
    FOR i:=1 TO n DO BEGIN
      ysav[i] := y[i];
      dysav[i]:= dydt[i];
    END;
    h:= htry;         {adjust try value increment}
    1:   hh := 0.5*h; {take 2 half time step}
    rk4(ysav,dysav,n,tsav,hh,ytemp);
    t:= tsav + hh;
    derivs(t,ytemp,dydt);
    rk4(ytemp,dydt,n,t,hh,y);
    t:= tsav + h;
    IF (t = tsav) THEN HALT(1);
    rk4(ysav,dysav,n,tsav,h,ytemp);
    errmax := 0; {Evaluate error}
    temp :=0;
    FOR i := 1 TO n DO BEGIN
      ytemp[i] := y[i] - ytemp[i];   {ytemp = estimated error}
      IF (yscal[i]>tiny) THEN temp := abs(ytemp[i]/yscal[i]);
      IF ( errmax < temp) THEN errmax := temp;
    END; { i }
    errmax:= errmax/eps;   {real error / requirement}
    IF (errmax > un) THEN  {Error too big, reduce h}
    BEGIN
      h := safety*h*exp(pshrnk*ln(errmax));
      GOTO 1;  {start again}
    END
    ELSE BEGIN   {step was a success}
      hdid := h; {Calculate next time step}
      IF (errmax > errcon) THEN
      hnext:=safety*h*exp(pgrow*ln(errmax))
      ELSE hnext:= 4.0*h;
    END;
    FOR i := 1 TO n DO y[i]:=y[i]+ytemp[i]*fcOR;
  END;  {rkqc}

  {==============================================================}
  PROCEDURE odeint(VAR ystart : Vecteurs_etat;
  NVAR  : INTEGER;
  t1,t2,eps,h1,hmIN: REAL_AR;
  VAR nok,nbad: INTEGER);
  {================================================================
  | Adapted from W.H.PRESS,B.P.FLANNERY,S.A.TEUKOLSKY,W.T.VETTERLING:
  | Numerical Recipes - Cambridge University Press - 1988 -
  ----------------------------------------------------------------------
  |  This procedure integrates the 1st order differential system:
  |       dy/dt = F(y,t)
  |  where y and F are  vectors of dimension Nvar, between times
  |  t1 and t2.
  |-----------------------------------------------------------------
  |  INPUTS:
  |  ystart= state vector at beginning of integration step.
  |  XJ    = object main inertial vectors.
  |  XV    = speed vectors of CoG in fixed coordinates.
  |  NVAR  = number of equations
  |  t1    = begin time of integration
  |  t2    = end time of integration
  |          (t2 may be > or < t1)
  |  eps   = absolute precision requirement for solution
  |  h1    = time step proposed as a try
  |  hmIN  = minimum time increment.
  |------------------------------------------------------------------
  |  OUTPUTS:
  |  ystart= vector of coordinates and speeds at end of step
  |  nok   = number of unchanged time steps
  |  nbad  = number of modified time steps
  ===================================================================}

  CONST maxstp=10000;
  two   = 2;
  zero  = 0;

  VAR   nstp,i               : INTEGER;
  tsav,t,hnext,hdid,h  : REAL_AR;
  yscal,y,dydt         : vecteurs_etat;

  BEGIN
    t:= t1;
    IF ( t2 > t1) THEN h :=   abs(h1)
    ELSE h := - abs(h1);
    nok:= 0;
    nbad :=0;
    FOR i := 1 TO NVAR DO y[i] := ystart[i];
    tsav := t - dtsave*two;
    FOR nstp :=1 TO maxstp DO
    BEGIN
      Derivs(t,y,dydt);
      FOR i:=1 TO NVAR DO yscal[i]:=abs(y[i])+abs(dydt[i]*h);
      IF (((t+h-t2)*(t+h-t1)) > zero) THEN h := t2 - t;
      rkqc(y,dydt,NVAR,t,h,eps,yscal,hdid,hnext);
      IF (hdid = h ) THEN nok  := nok +1
      ELSE nbad := nbad +1;
      IF (((t-t2)*(t2-t1)) >= zero) THEN
      BEGIN
        FOR i := 1 TO NVAR DO ystart[i] := y[i];
        EXIT; { c'est fini }
      END;
      IF (abs(hnext) < hmIN) THEN HALT(2);
      h := hnext;
    END; {nstp}
    HALT(3);
  END; { odeint }

  {===========================================================}
  PROCEDURE InitPointers(VAR ipt:Pointers);
  {===========================================================
  | Calculate pointers to store results in Table              }

  VAR i  : INTEGER;

  BEGIN
    ipt[0]:=0;
    FOR i:= 1 TO NVAR DO ipt[i]:=ipt[i-1]+NstepMax;
  END;

  {=========================================================================}
  PROCEDURE StoreResults(Temps: REAL_AR;
  Etat : Vecteurs_etat);
  {================================================================}
  VAR i : INTEGER;
  temp : single;

  BEGIN
    Table.Index(ipt[0]+iStep)^:=Temps;
    CPhi := Cos(Phi);
    Table.Index(ipt[1]+iStep)^:=1e-3*a*CPhi*(Lambda-Lambdab);
    Table.Index(ipt[2]+iStep)^:=1e-3*a*(Phi-Phib);
    Table.Index(ipt[3]+iStep)^:=1e-3*Altitude;
    Table.Index(ipt[4]+iStep)^:=V;
    Table.Index(ipt[5]+iStep)^:=Pdyn;
    DistanceBut(Lambda,lambdab,Phi,Phib,Distance);
    AzimutBut(Lambda,Lambdab,Phi,Phib,Azimut);
    Table.Index(ipt[6]+iStep)^:=1e-3*Distance;
    Table.Index(ipt[7]+iStep)^:=Deg*beta;
    Table.Index(ipt[8]+iStep)^:=Deg*Azimut;
  END;

  {=========================================================================}
  PROCEDURE Zoomer;  {draw end section of trajectory}

  CONST
        Range   = 10;
        DeltaR  = 1;
        Nptzoom = 30;

  {Zoom in the +/- Range [km] zone around target}

  VAR
        i,iVAR,Radius,trajH,trajV,OldMaxX,OldMaxY,Zoom,R : INTEGER;
        Titre,TitreX,TitreY    : ARRAY[0..20] OF CHAR;

  BEGIN

    clrscr;
    {define graph zone in physical coordinates}
    InitFenetre(CrtDC,10,-Range,Range,-Range,Range);
    {select blue pen}
    OldPen:=SelectObject(CrtDC,CrtPen);

    {draw 5 circles around target}
    FOR Radius := 5 DOWNTO 1 DO
    BEGIN
      R:=Radius*DeltaR;
      Circle(CrtDC,0,0,R,TRUE)    {Procedure Circle is defined in Graph_2D}
    END;

    SelectObject(CrtDC,InvPen);   {invisible white pen}
    LineXY(crtDC,-11.0,3.25);
    SelectObject(CrtDC,CrtPen);   {again blue pen}

    For i:=nstep-Nptzoom to nstep do
      x_y^[i]:=  Table.Index(ipt[2]+i)^;

    {prepare graph titles}
    StrPCopy(Titre, 'Zoom on impact Zone');
    StrPCopy(TitreX,'longitude [km]');
    StrPCopy(TitreY,'latitude [km]');

    For i:=nstep-Nptzoom to nstep do
      LineXY(CrtDC,Table.Index(ipt[1]+i)^,x_y^[i]);

    {write graph titles}
    Legendes(CrtDC,titre,titreX,titreY);

    SelectObject(CrtDC,OldPen);          {return to default black pen}
    SortieGraphique;                     {exit graph options}

    IF rep='n' THEN   {exit option}
    BEGIN
      DeleteObject(CrtPen);   {unit Time takes in charge to free Table}
      DeleteObject(InvPen);
      Donewincrt
    END;

    IF rep='i' THEN                                           {send to printer}
    BEGIN
      Xratio:=XrEPS;                        {values for EPSON Stylus Color 760}
      Yratio:=YrEPS;
      New(Pinfo,Init);                     {pointer to printer infos structure}
      WITH Pinfo^ DO                       {see unit WinPrint                 }
      BEGIN
        StartDoc('ROCKET V2');
        OldMaxX:=MaxX; OldMaxY:=MaxY;              {save current screen values}
        MaxX:=3600; MaxY:=2520;                               {360 dpi printer}
        InitFenetre(PrintDC,11,-Range,Range,-Range,Range);  {n° 11 for printer}
        {select blue pen}
        OldPen:=SelectObject(PrintDC,CrtPen);
        FOR Radius := 5 DOWNTO 1 DO
        BEGIN
          R:=Radius*DeltaR;
          Circle(PrintDC,0,0,R,TRUE)    {Procedure Circle is defined in Graph_2D}
        END;
        SelectObject(PrintDC,InvPen);   {invisible white pen}
        LineXY(PrintDC,-11.0,3.25);
        SelectObject(PrintDC,CrtPen);   {again blue pen}
        For i:=nstep-Nptzoom to nstep do
          LineXY(PrintDC,Table.Index(ipt[1]+i)^,x_y^[i]);
        {print graph titles}
        Legendes(PrintDC,titre,titreX,titreY);
        NewFrame;
        EndDoc
      END;
      Dispose(Pinfo,Done);
      {restore screen values}
      MaxX:=OldMaxX; MaxY:=OldMaxY;
      Xratio:=XrIBM; Yratio:=YrIBM; Fen11:=FALSE
    END
  END; {Zoomer}

  {=========================================================================}
  PROCEDURE DrawResults;

  LABEL Menu,EndProcedure;

  CONST Half = 0.5;

  VAR
       err,i, iVAR, OldMaxX, OldMaxY : INTEGER;
       Titre,TitreX,TitreY : ARRAY[0..30] OF CHAR;
       x1,x2,y1,y2 : REAL_AR;  c:STRING[2];

  BEGIN
    Menu: ClrScr;

    Rectangle(CrtDC,10,35,300,235);

    WRITELN;
    WRITELN('         DRAWING OPTIONS');
    WRITELN;
    FOR i:=1 TO NVAR DO
    begin
      gotoxy(3,3+i);
      WRITELN(i:3,' ',variable[i])
    end;
    gotoxy(3,12); WRITELN(TrajV:3,' Vertical trajectory');
    gotoxy(3,13); WRITELN(TrajH:3,' Horizontal trajectory');
    gotoxy(3,14); WRITELN(Zoom:3,' Zoom at impact point');
    gotoxy(3,15); WRITELN(fini:3,' Exit');
    WRITELN;

    IVAR:=0;
    Repeat
      gotoxy(7,17);
      WRITE('Your choice (1 to ',fini:2,') : '); Clreol;
      {$I-} Readln(c); {$I+} Val(c,IVAR,err);
      if (err<>0) or (IVAR<1) or (IVAR>12) then MessageBeep(0)
    Until IVAR IN [1..12];

    OldPen:=SelectObject(CrtDC,CrtPen);  {select blue pen}

    CASE iVAR OF
      1..nVAR : BEGIN
        FOR i:=1 TO nstep DO
          x_y^[i]:=Table.Index(ipt[iVAR]+i)^;
        clrscr;
        MinMax(nstep,x_y,y1,y2);
        InitFenetre(CrtDC,10,Table.Index(ipt[0]+1)^,Table.Index(ipt[0]+nstep)^,y1,y2);
        StrPCopy(Titre,variable[iVAR]);
        StrCopy(TitreX,'Time [s]');
        StrCopy(TitreY,' ');
        OldPen:=SelectObject(CrtDC,CrtPen);                    {select blue pen}
        MoveXY(CrtDC,Table.Index(ipt[0]+1)^,x_y^[1]);
        For i:=2 to nstep do
          LineXY(CrtDC,Table.Index(ipt[0]+i)^,x_y^[i]);
        Legendes(CrtDC,Titre,TitreX,TitreY);
        SelectObject(CrtDC,OldPen);                   {select default black pen}
        SortieGraphique;
        IF rep='n' THEN GOTO EndProcedure;                         {exit option}
        IF rep='i' THEN                                        {send to printer}
        BEGIN
          Xratio:=XrEPS; Yratio:=YrEPS;      {values for EPSON Stylus Color 760}
          New(Pinfo,Init);                  {pointer to printer infos structure}
          WITH Pinfo^ DO                    {see unit WinPrint                 }
          BEGIN
            StartDoc('ROCKET V2');
            OldMaxX:=MaxX; OldMaxY:=MaxY;
            MaxX:=3600; MaxY:=2520;                            {360 dpi printer}
            InitFenetre(PrintDC,11,Table.Index(ipt[0]+1)^,Table.Index(ipt[0]+nstep)^,y1,y2);
            MoveXY(PrintDC,Table.Index(ipt[0]+1)^,x_y^[1]);
            For i:=2 to nstep do
              LineXY(PrintDC,Table.Index(ipt[0]+i)^,x_y^[i]);
            Legendes(PrintDC,Titre,TitreX,TitreY);
            NewFrame;
            EndDoc
          END;
          Dispose(Pinfo,Done);
          {restore screen values}
          MaxX:=OldMaxX; MaxY:=OldMaxY;
          Xratio:=XrIBM; Yratio:=YrIBM; Fen11:=FALSE
        END
      END;
      TrajV   : BEGIN
        FOR i:=1 TO nstep DO
          x_y^[i]:=Table.Index(ipt[6]+i)^;
        MinMax(nstep,x_y,x1,x2);
        FOR i:=1 TO nstep DO
          x_y^[i]:=Table.Index(ipt[3]+i)^;
        MinMax(nstep,x_y,y1,y2);
        clrscr;
        InitFenetre(CrtDC,10,x1,x2,y1,y2);
        StrCopy(Titre,'Vertical trajectory');
        StrCopy(TitreX,'Dist. to target km');
        StrCopy(TitreY,'Z km');
        OldPen:=SelectObject(CrtDC,CrtPen);                    {select blue pen}
        MoveXY(CrtDC,Table.Index(ipt[6]+1)^,x_y^[1]);
        For i:=2 to nstep do
          LineXY(CrtDC,Table.Index(ipt[6]+i)^,x_y^[i]);
        Legendes(CrtDC,Titre,TitreX,TitreY);
        SelectObject(CrtDC,OldPen);                   {select default black pen}
        SortieGraphique;
        IF rep='n' THEN GOTO EndProcedure;                         {exit option}
        IF rep='i' THEN                                        {send to printer}
        BEGIN
          Xratio:=XrEPS; Yratio:=YrEPS;      {values for EPSON Stylus Color 760}
          New(Pinfo,Init);                
          WITH Pinfo^ DO
          BEGIN
            StartDoc('ROCKET V2');
            OldMaxX:=MaxX; OldMaxY:=MaxY;
            MaxX:=3600; MaxY:=2520;                            {360 dpi printer}
            InitFenetre(PrintDC,11,x1,x2,y1,y2);      {window n° 11 for printer}
            MoveXY(PrintDC,Table.Index(ipt[6]+1)^,x_y^[1]);
            For i:=2 to nstep do
              LineXY(PrintDC,Table.Index(ipt[6]+i)^,x_y^[i]);
            Legendes(PrintDC,Titre,TitreX,TitreY);
            NewFrame;
            EndDoc
          END;
          Dispose(Pinfo,Done);
          {restore screen values}
          MaxX:=OldMaxX; MaxY:=OldMaxY;
          Xratio:=XrIBM; Yratio:=YrIBM; Fen11:=FALSE
        END
      END;
      TrajH   : BEGIN
        FOR i:=1 TO nstep DO
          x_y^[i]:=Table.Index(ipt[1]+i)^;
        MinMax(nstep,x_y,x1,x2);
        FOR i:=1 TO nstep DO
          x_y^[i]:=Table.Index(ipt[2]+i)^;
        MinMax(nstep,x_y,y1,y2);
        clrscr;
        InitFenetre(CrtDC,10,x1,x2,y1,y2);
        StrCopy(Titre,'Horizontal trajectory');
        StrCopy(TitreX,'longitude km');
        StrCopy(TitreY,'latitude km');
        OldPen:=SelectObject(CrtDC,CrtPen);                    {select blue pen}
        MoveXY(CrtDC,Table.Index(ipt[1]+1)^,x_y^[1]);
        For i:=2 to nstep do
          LineXY(CrtDC,Table.Index(ipt[1]+i)^,x_y^[i]);
        Legendes(CrtDC,Titre,TitreX,TitreY);
        SelectObject(CrtDC,OldPen);                   {select default black pen}
        SortieGraphique;
        IF rep='n' THEN GOTO EndProcedure;
        IF rep='i' THEN                                        {send to printer}
        BEGIN
          Xratio:=XrEPS; Yratio:=YrEPS;      {values for EPSON Stylus Color 760}
          New(Pinfo,Init);                     
          WITH Pinfo^ DO                       
          BEGIN
            StartDoc('ROCKET V2');
            OldMaxX:=MaxX; OldMaxY:=MaxY;
            MaxX:=3600; MaxY:=2520;                            {360 dpi printer}
            InitFenetre(PrintDC,11,x1,x2,y1,y2);      {window n° 11 for printer}
            MoveXY(PrintDC,Table.Index(ipt[1]+1)^,x_y^[1]);
            For i:=2 to nstep do
              LineXY(PrintDC,Table.Index(ipt[1]+i)^,x_y^[i]);
            Legendes(PrintDC,Titre,TitreX,TitreY);
            NewFrame;
            EndDoc
          END;
          Dispose(Pinfo,Done);
          {restore screen values}
          MaxX:=OldMaxX; MaxY:=OldMaxY;
          Xratio:=XrIBM; Yratio:=YrIBM; Fen11:=FALSE
        END
      END;
      Zoom : Zoomer;
      Fini : GOTO EndProcedure;
    END; { Case }
    GOTO Menu;
    EndProcedure:                      {exit program and free memory}
    DeleteObject(CrtPen);   {Unit Time takes in charge to free Table}
    DeleteObject(InvPen);
    Dispose(x_y);
    DoneWinCrt
  END;   {DrawResults}

  {=========================================================================}
  PROCEDURE GroundZero;  {interpolate last iteration}

  VAR Azm,Azp,Deltat,
  DxPlus,DxMoins,DyPlus,DyMoins,
  EcartX,EcartY,t2,
  TPlus,TMoins,ZPlus,ZMoins,Qm,Qp,Vm,Vp     : REAL_AR;

  BEGIN
    TPlus  := Table.Index(ipt[0]+istep)^;
    Tmoins := Table.Index(ipt[0]+istep-1)^;
    ZPlus  := 1E3*Table.Index(ipt[3]+istep)^;    {back to meters}
    ZMoins := 1E3*Table.Index(ipt[3]+istep-1)^;  {back to meters}
    DxPlus := Table.Index(ipt[2]+istep)^;
    DxMoins:= Table.Index(ipt[2]+istep-1)^;
    DyPlus := Table.Index(ipt[1]+istep)^;
    DyMoins:= Table.Index(ipt[1]+istep-1)^;
    Vm     := Table.Index(ipt[4]+istep-1)^;
    Vp     := Table.Index(ipt[4]+istep)^;
    Qm     := Table.Index(ipt[5]+istep-1)^;
    Qp     := Table.Index(ipt[5]+istep)^;
    Azm    := Table.Index(ipt[8]+istep-1)^;
    Azp    := Table.Index(ipt[8]+istep)^;
    T2     := (ZPlus*Tmoins-TPlus*Zmoins)/(ZPLus-ZMoins);
    Temps  := T2;
    Deltat:= ZMoins/(ZMoins-ZPlus);
    Table.Index(ipt[0]+istep)^:= Temps;
    EcartY := DyMoins+Deltat*(DyPlus-DyMoins);
    EcartX := DxMoins+Deltat*(DxPlus-DxMoins);
    Table.Index(ipt[2]+istep)^ := EcartX;
    Table.Index(ipt[1]+istep)^ := EcartY;
    Distance := Sqrt(EcartX*EcartX+EcartY*EcartY);
    V :=  Vm+Deltat*(Vp-Vm);
    Table.Index(ipt[4]+istep)^:= V;
    Pdyn := Qm+Deltat*(Qp-Qm);
    Table.Index(ipt[5]+istep)^:= Pdyn;
    Table.Index(ipt[8]+istep)^:= Azm+Deltat*(Azp-Azm);
    Altitude := ZMoins+Deltat*(ZPlus-ZMoins);
    Table.Index(ipt[3]+istep)^:= 1e-3*Altitude;
    GoToXY(1,14);
    WRITE(t2:10:2,1e-3*Altitude:10:3,V:10:2,Pdyn:10:3,Distance:14:3);
    ClrEol;
    IF IprntSw THEN
    WRITELN(OUT,t2:10:2,1e-3*Altitude:10:3,V:10:2,Pdyn:10:3,Distance:14:3);
    GoToXY(1,16);
    WRITELN(line);
    WRITELN(' Target error [km]      in Longitude    in Latitude ');
    WRITELN('  ',Distance:8:3,'              ',EcartY:8:3,'        ',EcartX:8:3);
  END;  {GroundZero}


  VAR
  Impact : BOOLEAN;
  S : ARRAY[0..20] OF CHAR;
  s1: STRING[10];


  BEGIN  {of main program}
    WinCrtInit(' ROCKET  V2');
    New(x_y);
    CrtPen:=CreatePen(ps_Solid,1,blue);              {blue line thickness 1}
    InvPen:=CreatePen(ps_Solid,1,RGB(255,255,255));  {invisible white line}
    IprntSw  := false;    {true to send numeric results to printer}
    CheckBreak := true;
    ClrScr;

    IF IprntSw THEN
    Assign(OUT,'PRN')
    ELSE
    AssignCrt(OUT);
    ReWRITE(OUT);

    Table.Init(TablSize);              {allocate memory for Table}
    Impact := False;
    StartTiming;        {begin computing time - See unit Time.pas}
    nsmb := NbSMb;
    NstepMax := (TablSize DIV (4*(NVAR+1)) ) -1;
    eps := 1E-4;                     {required relative precision}
    Omega2 := Sqr(Omega);            {Earth angular speed power 2}
    MaitreCouple := Pi*D*D/4;               {missile main section}
    Q := Me/Tprop;                    {ergols mass output [kg/s] }
    Rad := PI/180;                        {degree to radian ratio}
    Deg := 180/PI;                        {radian to degree ratio}
    gamma := Rad*Pente;              {trajectory slope in radians}
    TgPente := Tan(Gamma);
    EpuisementErgols := False;

    InitPointers(ipt);

    Phi     := Rad*Phi1;
    Lambda  := Rad*Lambda1;
    Phib    := Rad*Phi2;
    Lambdab := Rad*Lambda2;

    AzimutBut(Lambda,Lambdab,Phi,Phib,Azimut);

    WRITELN;
    WRITELN(' Target azimut.................:',Deg*Azimut:10:2,' ø');
    DistanceBut(Lambda,Lambdab,Phi,Phib,Distance);
    WRITELN(' Distance to target............:',1E-3*Distance:10:2,' km');
    Tmax   := 400;                  {end integration time in sec.}
    dtSave := 2;                            {results every 2 sec.}

    Temps  := 0;
    h1     := 0.1;
    Istep  := 1;
    ConditionsInitiales(Phi,Lambda,Etat);
    Alpha := Azimut - CorrectionDeCap(Azimut);
    Sazim := SIN(Alpha);
    Cazim := Cos(Alpha);
    WRITELN(' Azimut corrected from Coriolis:',Deg*Alpha:10:2,' ø');

    Altitude := 0;
    V := 0;
    Beta := Rad*90;               {nozzle aligned on missile axis}
    StoreResults(Temps,Etat);

    dt:=h1;
    dtmIN :=0.001*h1;                 {smallest time step allowed}

    IF IprntSw THEN WRITELN('  Printer option on...');
    GoToXY(1,11);

    IF IprntSw THEN  {printer option}
    BEGIN
      WRITELN(OUT,'============================================================');
      WRITELN(OUT,'    Time      Altitude   Speed    Dyn. P    Distance Target ');
      WRITELN(OUT,'     [s]        [km]     [m/s]     [Mpa]         [km]');
      {screen display}
      WRITELN('    Time      Altitude   Speed     Dyn. P    Distance Target');
      WRITELN('     [s]        [km]     [m/s]      [Mpa]         [km]      ')
    END
    ELSE
    BEGIN    {screen default option}
      WRITELN('    Time      Altitude   Speed     Dyn. P    Distance Target ');
      WRITELN('     [s]        [km]     [m/s]      [Mpa]         [km]       ')
    END;

    {main integration loop}
    OldTextCoul:=SetTextColOR(CrtDC,red);      {write time value in red}
    WHILE (temps <= Tmax) AND (NOT Impact) AND (istep<NstepMax) DO
    BEGIN
      Inc(istep,1);
      t2:=temps+dtsave;
      odeint(Etat,nsmb,temps,t2,eps,dt,dtmIN,nok,nbad);
      Pdyn := 1e-6*Half*V2*AirDensity(Altitude);        {MPa}
      Str(t2:8:0,s1);
      StrPCopy(S,s1);
      TextOut(CrtDC,20,194,S,strlen(S));
      gotoxy(12,14);
      WRITE(1e-3*Altitude:10:3,V:10:2,Pdyn:10:3,(1e-3*Distance):14:3);
      ClrEol;
      IF IprntSw THEN
      WRITELN(OUT,t2:10:2,1e-3*Altitude:10:3,V:10:2,Pdyn:10:3,(1e-3*Distance):14:3);
      IF (nok >=0 ) THEN
      BEGIN
        Temps:=t2;
        StoreResults(temps,Etat);
      END
      ELSE
      BEGIN
        Dec(istep,1);
        Impact := true;                       {a way to stop calculations}
      END;
      IF (Altitude <= 0) THEN
      BEGIN
        Impact := true;
        GroundZero;                   {Interpolate state vector at impact}
      END
    END; (* WHILE *)
    SetTextColOR(CrtDC,OldTextCoul);                  {default text color}
    nstep:=iStep;                            {actual number of time steps}
    {DateHeure(OUT);}
    StopTiming;                                      {stop computing time}
    IF IprntSw THEN Close(OUT);
    WRITELN;
    WRITELN(' Computing time: ',Elapsed,' sec.');
    WRITELN;
    WRITELN(' Press any key to continue...');
    ReadKey;

    DrawResults;                            {takes in charge program exit}

  END.

{end of file rocketv2.pas}