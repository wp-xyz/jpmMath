{**************************************************************
!*  Calculate angle of refraction for an incident ray using   *
!*  Snell's Law                                               *
!* ---------------------------------------------------------- *
!* SAMPLE RUN:                                                *
!*                                                            *
!* Give indices for incident and refracting medium: 1.0 1.5   *
!*                                                            *
!* What is the angle of incidence? 30                         *
!*                                                            *
!* Refracted angle =  19.4712                                 *
!*                                                            *
!* ---------------------------------------------------------- *
!* Ref.: "Problem Solving with Fortran 90 By David R. Brooks, * 
!*        Springer-Verlag New York, 1997".                    *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
*                                 (www.jpmoreau.fr)           *
!**************************************************************
! Explanations:
! ------------               \ i! 
!                    incident \ !
!                      beam    \!
!                   ---------------------------
!                               !*
!                               ! *  refracted
!                               !r *   beam
!
!
!               Snell's Law:  n sin(i) = n sin(r)
!                              i          r 
!                           -1
!               or   r = sin  [(n / n ) sin(i)]
!                                i   r
!
!-------------------------------------------------------------}
Program Refract;

Uses WinCrt;

Var
  ni, nr:double;             {indices of refraction (dimensionless) }
  incident,refracted:double; {angles from perpendicular (deg) }
  DegToRad:double;           {required for trig functions }
  sinus,cosinus: double;     {auxiliary trig values }

Begin

  DegToRad := pi/180.0;

  Writeln;
  Write(' Give indices for incident and refracting medium: ');
  Readln(ni, nr);

  Writeln;
  Write(' What is the angle of incidence? ');
  Readln(incident);

{ Convert refracted angle to degrees before displaying its value }

  sinus := ni/nr*sin(incident*DegToRad);

  cosinus := sqrt(1.0-sinus*sinus);

  refracted := arctan(sinus/cosinus);

  Writeln;
  Writeln(' Refracted angle = ', (refracted/DegToRad):8:4);

  ReadKey;
  DoneWinCrt

End.

{end of file refract.pas}