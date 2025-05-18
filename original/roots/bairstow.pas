{****************************************************
*  Program to demonstrate the BAIRSTOW subroutine   *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*                TPW version by J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)          *
* ------------------------------------------------- *
* Example: Find two complex roots of polynomial:    *
*            f(x) = x^5-10x^4+35x^3-50x^2+24x       *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
* Input order of polynomial: 5                      *
*                                                   *
* Input the polynomial coefficients:                *
*                                                   *
*    A(0) =  0                                      *
*    A(1) =  24                                     *
*    A(2) =  -50                                    *
*    A(3) =  35                                     *
*    A(4) =  -10                                    *
*    A(5) =  1                                      *
*                                                   *
* Convergence factor: 0.000001                      *
* Maximum number of iterations: 50                  *
*                                                   *
* The roots found are:                              *
*                                                   *
*      X1 =  1.00000000000000E+0000                 *
*      Y1 =  0.00000000000000E+0000                 *
*                                                   *
*      X2 =  0.00000000000000E+0000                 *
*      Y2 =  0.00000000000000E+0000                 *
*                                                   *
* The number of iterations was:  11                 *
*                                                   *
****************************************************}
PROGRAM DEMO_BAIRSTOW;
Uses WinCrt;

Var
        A,B,C,D : Array[0..10] of DOUBLE;
        i,k,m,n : INTEGER;
        aa,bb,e,x1,x2,y1,y2 : DOUBLE;


{*************************************************
*       Bairstow complex root subroutine         *
* ---------------------------------------------- *
* This routine finds the complex conjugate roots *
* of a polynomial having real coefficients.      *
* ---------------------------------------------- *
* Reference: Computer Methods for Science and    *
*            Engineering by R.L. Lafara.         *
* ---------------------------------------------- *
* INPUTS:                                        *
*  Polynomial coefficients      : A(0) to A(m)   *
*  Order of polynomial (>=4)    : m              *
*  Initial guess                : a and b        *
*  Convergence factor           : e              *
*  Maximum number of iterations : n              *
* OUTPUTS:                                       *
*  Two conjugate complex roots  : x1,y1  x2,y2   *
*  Number of iterations         : k              *
*************************************************}
PROCEDURE Bairstow;
Label 100,200,300,400;
Var i,j : INTEGER;
    a1,b1,cc,dd,d2 : DOUBLE;
Begin
  {Normalize the A[i] series}
  for i:=0 to m do
    C[i]:=A[i]/A[m];
  {Take initial estimates for aa and bb}
  k:=0; B[m]:=1.0;
  {Start iteration sequence}
100: B[m-1]:=C[m-1]-aa;
  for j:=2 to m-1 do
    B[m-j]:=C[m-j]-aa*B[m+1-j]-bb*B[m+2-j];
  B[0]:=C[0]-bb*B[2];
  D[m-1]:=-1.0; D[m-2]:=-B[m-1]+aa;
  for j:=3 to m-1 do
    D[m-j]:=-B[m+1-j]-aa*D[m+1-j]-bb*D[m+2-j];
  D[0]:=-bb*D[2];
  d2:=-B[2]-bb*D[3];
  dd:=D[1]*d2-D[0]*D[2];
  a1:=-B[1]*d2+B[0]*D[2]; a1:=a1/dd;
  b1:=-D[1]*B[0]+D[0]*B[1]; b1:=b1/dd;
  aa:=aa+a1; bb:=bb+b1; k:=k+1;
  {Test for the number of iterations}
  if k>=n then goto 200;
  {Test for convergence}
  if ABS(a1)+ABS(b1)>e*e then goto 100;
  {Extract roots from quadratic equation}
200: cc:=aa*aa-4.0*bb;
  {Test to see if a complex root}
  if cc>0 then goto 300;
  x1:=-aa; x2:=x1; y1:=SQRT(-cc); y2:=-y1;
  goto 400;
300: x1:=-aa+SQRT(cc);
  x2:=-aa-SQRT(cc); y1:=0; y2:=y1;
400: x1:=x1/2.0; x2:=x2/2.0; y1:=y1/2.0; y2:=y2/2.0
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

  Bairstow;     {Call Bairstow subroutione}

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

{End of file bairstow.pas}