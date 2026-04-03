{**********************************************************************
!*       Purpose: This program computes Airy functions and their      *
!*                derivatives using subroutine AIRYA                  *
!*       Input:   xstart, xend --- Arguments of Airy function         * 
!*                xstep --- x increment                               *
!*       Output:  AI --- Ai(x)                                        * 
!*                BI --- Bi(x)                                        *
!*                AD --- Ai'(x)                                       *
!*                BD --- Bi'(x)                                       *
!*       Example:                                                     *
!*                xstart =0 xend = 30 xstep = 10                      *
!*                                                                    *
!*   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)      *
!*  ----------------------------------------------------------------  *
!*   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00   *
!*  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10   *
!*  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26   *
!*  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48   *
!*                                                                    *
!*   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)     *
!*  ----------------------------------------------------------------  *
!*   0       .35502805      .61492663     -.25881940      .44828836   *
!*  10       .04024124     -.31467983      .99626504      .11941411   *
!*  20      -.17640613     -.20013931      .89286286     -.79142903   *
!*  30      -.08796819     -.22444694     1.22862060     -.48369473   *
!* ------------------------------------------------------------------ *
!* REFERENCE: "Fortran Routines for Computation of Special Functions, *
!*             jin.ece.uiuc.edu/routines/routines.html".              *
!*                                                                    *
!*                               Pascal Release By J-P Moreau, Paris. *
!*                                        (www.jpmoreau.fr)           *
!*********************************************************************}        
Program AiryLaz;

uses
  jpmTypes, jpmSpecialFunc;
var
  s: String;
  x, xStart, xEnd, xStep: Double;
  AI, BI, AD, BD: Double;
  res: Integer;

begin
  WriteLn;
  Write(' Please enter x start (0): ');
  ReadLn(s);
  if s = '' then
    xStart := 0
  else
  begin
    Val(s, xStart, res);
    if res <> 0 then Halt;
  end;
  Write('               x end (30): ');
  ReadLn(s);
  if s = '' then
    xEnd := 30
  else
  begin
    Val(s, xEnd, res);
    if res <> 0 then Halt;
  end;
  Write('          and x step (10): ');
  ReadLn(s);
  if s = '' then
    xStep := 10
  else
  begin
    Val(s, xStep, res);
    if res <> 0 then Halt;
  end;

  x := xStart;
  while x <= xEnd do
  begin
    AiryA(x, AI, BI, AD, BD);
    if (x = xStart) then
    begin
      WriteLn;
      WriteLn('    x         Ai(x)             Bi(x)            Ai''(x)            Bi''(x)');
      WriteLn('  -----------------------------------------------------------------------------')
    end;
    WriteLn(x:5:0, '   ' , AI:15, '   ' , BI:15, '   ' , AD:15, '   ', BD:15);
    x := x + xStep;
  end;

  x := xStart;
  while x <= xEnd do
  begin
    AiryA(-x, AI, BI, AD, BD);
    if (x = xStart) then
    begin
      WriteLn;
      WriteLn('    x         Ai(x)             Bi(x)            Ai''(x)            Bi''(x)');
      WriteLn('  -----------------------------------------------------------------------------')
    end;
    WriteLn(  -x:5:0, '   ', AI:13:8, '   ' , BI:15:8, '   ' , AD:14:8, '   ', BD:15:8);
    x := x + xStep;
  end;

  WriteLn;
  Write('Press [ENTER] to close...');
  ReadLn;
end.

{end of file AiryLaz.pas}
