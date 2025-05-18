{******************************************************
* Program to demonstrate the use of multi-dimensional *
*     Steepest Descent Optimization subroutine        *
* --------------------------------------------------- *
*  Reference: BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*                  TPW Version By J-P Moreau, Paris.  *
*                         (www.jpmoreau.fr)           *
* --------------------------------------------------- *
* Example:   Find a local maximum of                  *
*            F(x,y,z) = sin(x)+2*cos(y)-sin(z)        *
*                                                     *
* SAMPLE RUN:                                         *
*                                                     *
* How many dimensions: 3                              *
*                                                     *
* Convergence criterion: .000000001                   *
*                                                     *
* Maximum number of iterations: 50                    *
*                                                     *
* Starting constant: 1                                *
*                                                     *
* Input the starting point:                           *
*                                                     *
*     X(1) = 1                                        *
*     X(2) = 1                                        *
*     X(3) = 1                                        *
*                                                     *
* The results are:                                    *
*                                                     *
*     X(1) = 1.5707973                                *
*     X(2) = -0.0000170                               *
*     X(3) = -4.7123856                               *
*                                                     *
* Local maximum found = 4.0000000                     *
*                                                     *
* The number of steps was: 33                         *
*                                                     *
******************************************************}
PROGRAM DEMO_STEEPDS;
Uses WinCrt;

Var
        D, Y: Array[1..3] of DOUBLE;
        X,X1: Array[1..10] of DOUBLE;
        e,xk : DOUBLE;
        i,l,m,n : INTEGER;


{**********************************************
  Function subroutine                         }
  FUNCTION Eval:DOUBLE;
  Begin
    Eval := SIN(X[1])+2.0*COS(X[2])-SIN(X[3])
  End;
{*********************************************}


{**************************************************
*    Steepest descent optimization subroutine     *
* ----------------------------------------------- *
* This routine finds the local maximum or minimum *
* of an L-dimensional function using the method   *
* of steepest decent, or the gradient.            *
* The function must be available in subroutine    *
* 1000. In this version, finite differences are   *
* used to calculate the L partial derivatives.    *
* ----------------------------------------------- *
* INPUTS:                                         *
*   l - The dimension of function to study        *
*   e - The convergence criteria                  *
*   m - The maximum number of iterations          *
*   xk - A starting constant                      *
*   X(i) - Initial values of variables            *
* OUTPUTS:                                        *
*   X(i) - The locally optimum set                *
*   Eval - The local maximum found                *
*   n - The number of iterations performed,       *
**************************************************}
PROCEDURE Steepda;
Label 50,51,100,200,fin;
Const MACHEPS = 1E-15;
Var dd : DOUBLE;
    i,j  : INTEGER;

  Procedure Utilit1;
  Var i:INTEGER;
  Begin
    {Find the magnitude of the gradient}
    dd:=0.0;
    for i:=1 to l do
      dd:=dd+D[i]*D[i];
    dd:=SQRT(dd)
  End;

  Procedure Utilit2;
  Var i:INTEGER;
  Begin
    {Update the X[i) }
    for i:=1 to l do
    begin
      {Save old values}
      X1[i]:=X[i];
      X[i]:=X[i]+xk*D[i]/dd
    end;
    Y[3]:=Eval
  End;

  { Find approximations of partial derivatives D(i)
    by finite differences }
  Procedure Derivatives;
  Var a,b,yy : DOUBLE;
      i : INTEGER;
  Begin
    for i:=1 to l do
    begin
      {Save X(i) }
      a:=X[i];
      {Find increment}
      b:=D[i]*xk/(2.0*dd);
      {Move increment in X(i) }
      X[i]:=X[i]+b;
      {Obtain yy}
      yy:=Eval;
      {Guard against divide by zero near maximum}
      if b=0 then b:=1e-12;
      {Update D(i) }
      D[i]:=(yy-Y[3])/b;
      {Guard against locked up derivative}
      if D[i]=0 then D[i]:=1e-5;
      {Restore X(i) and yy}
      X[i]:=a; yy:=Y[3]
    end;
    {Obtain dd}
    Utilit1;
  End;

Begin {Steepds}
  n:=0;
  {The routine needs three values of Y to get started
   Generate starting D(i) values
   These are not even good guesses and slow the program a little}
  dd:=1.0;
  D[1]:=1.0/SQRT(l);
  for i:=2 to l do D[i]:=D[i-1];
  {Start initial probe}
  for j:=1 to 3 do
  begin
    {Obtain yy and D(i) }
    Y[j]:=Eval;
    {Update X(i) }
    Utilit1;
    Utilit2;
  end;
  {We now have a history to base the subsequent search on
   Accelerate search if approach is monotonic }
50: if ABS(Y[2]-Y[1])<MACHEPS then goto 51;
  if (Y[3]-Y[2])/(Y[2]-Y[1])>0 then xk:=xk*1.2;
  {Decelerate if heading the wrong way}
51: if Y[3]<Y[2] then xk:=xk/2.0;
  {Update the Y[i) if value has increased}
  if Y[3]>Y[2] then goto 100;
  {Restore the X(i) }
  for i:=1 to l do X[i]:=X1[i];
  goto 200;
100: Y[1]:=Y[2]; Y[2]:=Y[3];
  {Obtain new values}
200: Y[3]:=Eval;
  Derivatives; {Get D(i) }
  {if dd=0 then the precision limit of the computer has been reached}
  if dd=0 then goto fin;
  {Update X(i) }
  Utilit2;
  {Check for maximum iterations and convergence}
  n:=n+1;
  if n>=m then goto fin;
  if ABS(Y[3]-Y[2])<e then goto fin;
  {Try another iteration}
  goto 50;
fin: End; {of subroutine Steepds}


{main program}
BEGIN
  clrscr;
  writeln;
  write(' How many dimensions: '); read(l);
  writeln;
  write(' Convergence criterion: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(m);
  writeln;
  write(' Starting constant: '); read(xk);
  writeln;
  writeln(' Input the starting point: ');
  writeln;
  for i:=1 to l do
  begin
    write('    X(',i,') = '); read(X[i])
  end;

  Steepda;   {Call steepest descent optimization subroutine}

  writeln;; writeln;
  writeln(' The results are:');
  writeln;
  for i:=1 to l do
  begin
    writeln('    X(',i,') = ',X[i]:10:7)
  end;
  writeln;
  Writeln(' Local maximum found = ',Eval:10:7);
  writeln;
  writeln(' The number of iterations was ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file steepda.pas}





