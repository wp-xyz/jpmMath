{****************************************************
*    Program to demonstrate the LIN subroutine      *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
* ------------------------------------------------- *
* Example: Find two complex roots of polynomial:    *
*            f(x) = x^5+5x^4+10x^3+10x^2+5x+1       *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
* Input order of polynomial: 5                      *
*                                                   *
* Input the polynomial coefficients:                *
*                                                   *
*    A(0) = ? 1                                     *
*    A(1) = ? 5                                     *
*    A(2) = ? 10                                    *
*    A(3) = ? 10                                    *
*    A(4) = ? 5                                     *
*    A(5) = ? 1                                     *
*                                                   *
* Convergence factor: 0                             *
* Maximum number of iterations: 1000                *
*                                                   *
* The roots found are:                              *
*                                                   *
*      X1 = -9.37616596879390E-0001                 *
*      Y1 =  2.58117301321899E-0002                 *
*                                                   *
*      X2 = -9.37616596879390E-0001                 *
*      Y2 = -2.58117301321899E-0002                 *
*                                                   *
* The number of iterations was: 1000                *
*                                                   *
*                 TPW Version By J-P Moreau, Paris. *
*                         ww.jpmoreau.fr)           *
****************************************************}
PROGRAM DEMO_LIN;
Uses WinCrt;

Var
        A,B,C : Array[0..10] of DOUBLE;
        i,k,m,n : INTEGER;
        aa,bb,e,x1,x2,y1,y2 : DOUBLE;



{************************************************
*     Polynomial complex roots subroutine       *
* --------------------------------------------- *
*     This routine uses the Lin's method        *
* --------------------------------------------- *
* Reference: A Practical Guide to Computer      *
*            Method for Engineers by Shoup      *
* --------------------------------------------- *
* INPUTS:                                       *
*  Polynomial coefficients      : A(0) to A(m)  *
*  Order of polynomial          : m             *
*  Initial guess                : a and b       *
*  Convergence factor           : e             *
*  Maximum number of iterations : n             *
* OUTPUTS:                                      *
*  Two complex roots            : x,y1  x2,y2   *
*  Number of iterations         : k             *
************************************************}
PROCEDURE LIN;
Label 100,200,300,400,500;
Var i,j:INTEGER;
    a1,b1,cc:DOUBLE;
Begin
  {Normalize the A(i) series}
  for i:=0 to m do
    C[i]:=A[i]/A[m];
  {Start iteration
   Set initial guess for the quadratic coefficients}
  B[0]:=0; B[1]:=0;
100: B[m-1]:=C[m-1]-aa;
  B[m-2]:=C[m-2]-aa*B[m-1]-bb;
  for j:=3 to m do
    B[m-j]:=C[m-j]-aa*B[m+1-j]-bb*B[m+2-j];
  {Guard against divide by zero}
  if B[2]<>0 then goto 200;
  aa:=aa+1e-7; bb:=bb-1e-7;
200: a1:=(C[1]-bb*B[3])/B[2];
  b1:=C[0]/B[2];
  k:=k+1;
  {Test for the number of iterations}
  if k>=n then goto 300;
  {Test for convergence}
  if ABS(aa-a1)+ABS(bb-b1) < e*e then goto 300;
  aa:=a1; bb:=b1;
  {Next iteration}
  goto 100;
300: aa:=a1; bb:=b1;
  cc:=aa*aa-4.0*bb;
  {Is there an imaginary part?}
  if cc>0 then goto 400;
  y1:=SQRT(-cc); y2:=-y1; x1:=-aa; x2:=x1;
  goto 500;
400: y1:=0; y2:=y1; x1:=-aa+SQRT(cc); x2:=-aa-SQRT(cc);
500: x1:=x1/2.0;
  x2:=x2/2.0; y1:=y1/2.0; y2:=y2/2.0
End;

{main program}
BEGIN
  clrscr;
  Writeln;
  Write(' Order of polynomial: '); read(m);
  Writeln;
  Writeln(' Input the polynomial coefficients:');
  Writeln;
  for i:=0 to m do
  begin
    Write('     A(',i,') = '); readln(A[i])
  end;
  Writeln;
  Write(' Convergence factor: '); read(e);
  Write(' Maximum number of iterations: '); read(n);

  aa:=PI; bb:=SQRT(2.0);

  LIN;       {Call LIN subroutione}

  Writeln;
  Writeln(' The roots found are:');
  Writeln;
  Writeln('      X1 = ', x1);
  Writeln('      Y1 = ', y1);
  Writeln;
  Writeln('      X2 = ', x2);
  Writeln('      Y2 = ', y2);
  Writeln;
  Writeln(' The number of iterations was: ', k);
  Writeln;
  ReadKey; DoneWinCrt
END.

{End of file lin.pas}