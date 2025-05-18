{*****************************************************
*    This program caculates the absolute stellar     *
*     based on relative magintude and distance       *
* -------------------------------------------------- *
* SAMPLE RUN:                                        *
*                                                    *
* To calculate absolute stellar magnitude:           *
* Give relative mag. and distance in parsecs: 10 3   *
*                                                    *
* A star with relative magnitude  10.00              *
* at a distance of   3.0 parsecs                     *
* has an absolute magnitude of  12.61                *
*                                                    *
* -------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90 By David R. *
*        Brooks, Springer-Verlag New York, 1997".    *
*                                                    *
*                 TPW Release By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)           *
******************************************************
! Explanations;
! ------------
!   The absolute magnitude M of a star is related to its
! relative magnitude m and and the distance to a star r,
! measured in parsecs, where 1 parsec = 3.26 light years,
! by the equation:
!
!                  M = m + 5 - 5 log10(r)
! 
!   According to this equation, a star with a relative
! magnitude of +1 at a distance of 10 parsecs has an 
! absolute magnitude of +1. The larger the magnitude, the 
! dimmer the star. Sirius is a very bright star with a
! relative magnitude of -1.58. Stars visible to the naked
! eye range mostly from about +1 to +6 in relative magni-
! tude. 
!   The dimmest star that can be seen with the 200-inch
! Hale telescope has a magnitude of about +23.
!--------------------------------------------------------}  
Program Starmag;

Uses WinCrt;

Var
    abs_mag, rel_mag, parsecs: Real;

Begin

  writeln;
  writeln(' To calculate absolute stellar magnitude:');
  write(' Give relative mag. and distance in parsecs: ');
		 
  readln(rel_mag, parsecs);
  
  abs_mag := rel_mag + 5.0 - 5.0 * Ln(parsecs)/2.302585;

  writeln;
  writeln(' A star with relative magnitude ',rel_mag:6:2);
  writeln(' at a distance of ', parsecs:5:1,' parsecs');
  writeln(' has an absolute magnitude of ', abs_mag:6:2);
  writeln;

  ReadKey;
  DoneWinCrt

End.

{end of file starmag.pas}