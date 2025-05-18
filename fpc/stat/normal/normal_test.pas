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
program normal_test;

uses
  uNormal;

var
  P, X, Y: Double;

begin
  Writeln;
  Write(' Value of variable: ');
  Readln(X);
  Writeln;
  P := phi(X);
  Writeln(' Probability = ', P);
  Writeln;
  Writeln(' Verify:');
  Y := Normal(P);
  Writeln('   X = ', Y);
  Writeln;
  Write('Press ENTER to close...');
  ReadLn;
end.

