{***************************************************
*  Program to demonstrate Lagrange interpolation   *
*     of Function SIN(X) in double precision       *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  Lagrange interpolation of SIN(X):               *
*                                                  *
*    Input X: 0.5                                  *
*    Order of interpolation: 4                     *
*                                                  *
*    SIN(X) = 0.47942552 (exact value: 0.47942554) *
*                                                  *
***************************************************}
PROGRAM Lagrange;
Uses WinCrt;

Label   100;

Var
        
        X,Y : Array[0..14] of DOUBLE;
        xx,yy : DOUBLE;
        iv,n  : integer;


{*******************************************************
*          Lagrange interpolation subroutine           *
* n is the level of the interpolation ( Ex. n=2 is     *
* quadratic ). v is the total number of table values.  *
* X(i), Y(i) are the coordinate table values, Y(i)     *
* being the dependant variable. The X(i) may be arbi-  *
* trarily spaced.  x is the interpolation point which  *
* is assumed to be in the interval  with at least one  *
* table value to the left, and n to the right. If this *
* is violated, n will be set to zero. It is assumed    *
* that the table values are in ascending X(i) order.   *
*******************************************************}
PROCEDURE Interpol_Lagrange;
Label 100,200,300,400,fin;
Var i,j,k : INTEGER;
    XL    : Array[0..9] of DOUBLE;
Begin
  {Check to see if interpolation point is correct}
  if xx < X[1] then goto 100; 
  if xx <= X[iv-n] then goto 200;
  {An error has been encountered}
100: n:=0 ; goto fin;
  {Find the relevant table interval}
200: i:=0;
300: i:=i+1;
  if xx > X[i] then goto 300;
  i:=i-1;
  {Begin interpolation}
  for j:=0 to n do  XL[j]:=1.0;
  yy:=0.0;
  for k:=0 to n do
  begin
    for j:=0 to n do
    begin
      if j=k then goto 400;
      XL[k]:=XL[k]*(xx-X[j+i])/(X[i+k]-X[j+i]);
400: end;
    yy:=yy+XL[k]*Y[i+k]
  end;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Lagrange interpolation of SIN(X):');
  iv:=14;
  {Input sine table
  -----------------------------------------------------------------
   Sine table values from  Handbook of mathematical functions
   by M. Abramowitz and I.A. Stegun, NBS, june 1964
  ----------------------------------------------------------------}
  X[1]:=0.000; Y[1]:=0.00000000; X[2]:=0.125; Y[2]:=0.12467473;
  X[3]:=0.217; Y[3]:=0.21530095; X[4]:=0.299; Y[4]:=0.29456472;
  X[5]:=0.376; Y[5]:=0.36720285; X[6]:=0.450; Y[6]:=0.43496553;
  X[7]:=0.520; Y[7]:=0.49688014; X[8]:=0.589; Y[8]:=0.55552980;
  X[9]:=0.656; Y[9]:=0.60995199; X[10]:=0.721; Y[10]:=0.66013615;
  X[11]:=0.7853981634; Y[11]:=0.7071067812;
  X[12]:=0.849; Y[12]:=0.75062005; X[13]:=0.911; Y[13]:=0.79011709;
  X[14]:=0.972; Y[14]:=0.82601466;
  {----------------------------------------------------------------

   Input interpolation point }
100: writeln;

  write('   Input X: '); read(xx);
  write('   Order of interpolation: '); read(n);

  Interpol_Lagrange;

  if n<>0 then
  begin
    writeln;
    writeln('   SIN(X) = ',yy:10:8,' (exact value: ',sin(xx):10:8,')');
    goto 100
  end;
  writeln;
  DoneWinCrt
END.

{End of file lagrange.pas}
