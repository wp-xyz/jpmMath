{*******************************************************
* Calculate relativistic mass and speed of an electron *
* accelerated in an electron gun                       *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* Give electron gun voltage in volts: 1e6              *
*                                                      *
* Relativistic mass (kg) and speed (m/s):              *
*  2.69339461904430E-0030  2.82124894363268E+0008      *
*                                                      *
* ---------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90 By David R.   *  
*        Brooks, Springer-Verlag New York, 1997".      *
*                                                      *
*                  TPW Release By J-P Moreau, Paris.   *
*                          (www.jpmoreau.fr)           *
********************************************************
!Explanations:
!------------
!
! An electron accelerated by a voltage V in an electron gun
!
!                              2     2
! acquires an energy of Ve = mc - m c,  where
!                                  0
!                                        -19
!  charge on an electron   e = 1.602 x 10    coulomb
!
!                                         -31
!  rest mass               m  = 9.109 x 10    kg
!                           0
!                                         8 
!  speed of light          c = 2.9979 x 10  m/s
!
! The speed v of an electron of relativistic mass m (kg) is
! obtained from
!                                          2
!                 m/m  = 1 / sqrt(1 - (v/c) )
!                    0
!
!----------------------------------------------------------}
Program Relativistic_Mass;

Uses WinCrt1;

Var
 
  rest_mass, rela_mass,         { kg }
  voltage,                      { volt }
  speed,                        { m/s }
  e,                            { electron charge in coulomb }
  c: double;                    { speed of light in m/s }

Begin

  e:=1.602e-19; c:=2.9979e8; rest_mass:=9.109e-31;

  writeln;
  write(' Give electron gun voltage in volts: ');

  readln(voltage);

  rela_mass := (voltage*e + rest_mass*c*c) / (c*c);

  speed := c*sqrt(1.0 - Sqr(rest_mass/rela_mass));

  writeln;
  writeln(' Relativistic mass (kg) and speed (m/s):');
  writeln(' ', rela_mass, ' ', speed);
  writeln;

  ReadKey;
  DoneWinCrt

End.

{ End of file rel_mass.pas }