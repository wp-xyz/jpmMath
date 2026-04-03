{*****************************************************************
*      Purpose: This program computes the beta function          * 
*               B(p,q) for p > 0 and q > 0 using                 *
*               subroutine BETA                                  *
*      Input :  p  --- Parameter  ( p > 0 )                      *
*               q  --- Parameter  ( q > 0 )                      *
*      Output:  BT --- B(p,q)                                    *
*      Examples:                                                 *
*                p       q           B(p,q)                      *
*              ---------------------------------                 *
*               1.5     2.0     .2666666667D+00                  *
*               2.5     2.0     .1142857143D+00                  *
*               1.5     3.0     .1523809524D+00                  *
* -------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special        *
*             Functions jin.ece.uiuc.edu/routines/routines.html" *
*                                                                *
*                              TPW Release By J-P Moreau, Paris. *
*                                      (www.jpmoreau.fr)         *
*****************************************************************}
Program betaLaz;

uses
  jpmTypes, jpmSpecialFunc;
var
  BT, P, Q: Float;
begin
  WriteLn;
  WriteLn('     p       q              B(p,q)         ');
  WriteLn('  -----------------------------------------');
  P := 1.5;
  Q := 2.0;
  BT := Beta(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '    ', BT);
  P := 2.5;
  Q := 2.0;
  BT := Beta(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '    ', BT);
  P := 1.5;
  Q := 3.0;
  BT := Beta(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '    ', BT);
  WriteLn;
  Write('Press ENTER to Close...');
  ReadLn;
end.

{ end of file mbeta.pas }
