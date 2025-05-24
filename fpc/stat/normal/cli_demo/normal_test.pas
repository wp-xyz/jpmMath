{*********************************************************************
*         NORMAL AND INVERSE NORMAL PROBABILITY FUNCTIONS            *
* ------------------------------------------------------------------ *
* SAMPLE RUN:                                                        *
*                                                                    *
* Value of variable: 0.2                                             *
*                                                                    *
* Probability =  3.91042693975456E-0001                              *
*                                                                    *
* Verify:                                                            *
*   X =  2.00000286102295E-0001                                      *
*                                                                    *
*                             Pascal Version By J-P Moreau, Paris.   *
*                                      (www.jpmoreau.fr)             *
*********************************************************************}
program normal_test; // named "normal" in the original code

uses
  jpmTypes, jpmStats;
var
  P, X, Y: Float;
  s: String;
  res: Integer;
begin
  repeat
    Writeln;
    Write(' Value of variable (or ENTER to close): ');
    Readln(s);
    val(s, X, res);
    if res <> 0 then
      break;
    Writeln;
    P := NormalDist(X);
    Writeln(' Probability = ', P, ' (', P:0:3, ')');
//    Writeln;
    Y := InvNormalDist(P);
    Writeln(' Verify:   X = ', Y, ' (', Y:0:3, ')');
  until false;
end.

