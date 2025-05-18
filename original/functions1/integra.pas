{***************************************************
* Program to demonstrate the general integration   *
* subroutine. Example is the integral of SIN(X)    *
* from X1 to X2 in double precision.               *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*             Pascal version by J-P Moreau, Paris  *
*                      (www.jpmoreau.fr)           *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  Integration of SIN(X)                           *
*                                                  *
*  Input start point  end point: 0  0.5            *
*  Integral from 0.000 to 0.500 equals  0.1224231  *
*  Error check: 1                                  *
*                                                  * 
***************************************************}
PROGRAM Integra;
Uses WinCrt;

CONST
        SIZE = 14;

VAR
        X  : Array[0..SIZE] of double;
        Y  : Array[0..SIZE] of double;
        XM : Array[0..SIZE+3] of double;
        Z  : Array[0..SIZE] of double;

        a,b,x1,x2,xx,yy,zz : double;
        i,iv,n,n1,z1 : INTEGER;



{*******************************************************
*          Akima spline fitting subroutine             *
* ---------------------------------------------------- *
* The input table is X(i), Y(i), where Y(i) is the     *
* dependant variable. The interpolation point is xx,   *
* which is assumed to be in the interval of the table  *
* with at least one table value to the left, and three *
* to the right. The interpolated returned value is yy. *
* n is returned as an error check (n=0 implies error). *
* It is also assumed that the X(i) are in ascending    *
* order.                                               *
*******************************************************}
PROCEDURE Interpol_Akima;
Label 100,200,300,fin;
Var   i:integer;
Begin
  n:=1;
  {special case xx:=0}
  if xx=0.0 then
  begin
    yy:=0.0; goto fin
  end;
  {Check to see if interpolation point is correct}
  if (xx<X[1]) or (xx>=X[iv-3]) then
  begin
    n:=0 ; goto fin
  end;
  X[0]:=2.0*X[1]-X[2];
  {Calculate Akima coefficients, a and b}
  for i:=1 to iv-1 do
    {Shift i to i+2}
    XM[i+2]:=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
  XM[iv+2]:=2.0*XM[iv+1]-XM[iv];
  XM[iv+3]:=2.0*XM[iv+2]-XM[iv+1];
  XM[2]:=2.0*Xm[3]-XM[4];
  XM[1]:=2.0*XM[2]-XM[3];
  for i:=1 to iv do
  begin
    a:=ABS(XM[i+3]-XM[i+2]);
    b:=ABS(XM[i+1]-XM[i]);
    if a+b<>0 then goto 100;
    Z[i]:=(XM[i+2]+XM[i+1])/2.0;
    goto 200;
100: Z[i]:=(a*XM[i+1]+b*XM[i+2])/(a+b);
200: end;
  {Find relevant table interval}
  i:=0;
300: i:=i+1;
  if xx>X[i] then goto 300;
  i:=i-1;
  {Begin interpolation}
  b:=X[i+1]-X[i];
  a:=xx-X[i];
  yy:=Y[i]+Z[i]*a+(3.0*XM[i+2]-2.0*Z[i]-Z[i+1])*a*a/b;
  yy:=yy+(Z[i]+Z[i+1]-2.0*XM[i+2])*a*a*a/(b*b);
fin: End;

PROCEDURE Trapez1(VAR xi1:double); FORWARD;
PROCEDURE Trapez2(VAR xi2:double); FORWARD;

{*******************************************************
*       General integration subroutine (ITEG)          *
* Interpolation by Akima (or other). Integration by    *
* enhanced trapezoidal rule with Richardson extrapola- *
* tion to give cubic accuracy.  The integration range  *
* is (x1,x2). It is assumed that x1<x2, and that there *
* is at least one table value to the left of x1, and   *
* three to the right of x2. The result is returned in  *
* zz. An error check is returned in z1 (z1=0-->error). *
*******************************************************}
PROCEDURE Integrale;
Label 100,200,fin;
Var xi1,xi2,x3 : double;
Begin
  zz:=0 ; z1:=0;
  {Check to see if end points are correct}
  if x1<X[1] then goto fin;
  if x2>X[iv-3] then goto fin;
  {If x1>x2 then switch and set flag}
  if x1<x2 then goto 100;
  x3:=x1 ; x1:=x2 ; x2:=x3 ; z1:=1;
100: if x2=x1 then goto fin;
  {Start trapezoidal integrations
   First integration to get xi1 }
  Trapez1(xi1);
  {Second round to get xi2}
  Trapez2(xi2);
  {Richardson extrapolation}
  zz:=4.0*xi2/3.0-xi1/3.0;
  {Check to see if the end points have been reversed}
  if z1=0 then goto 200;
  zz:=-zz ; x2:=x1 ; x1:=x3;
  {Reset error flag}
200: z1:=1; 
  {zz is the integral desired}
fin: End;


{*****************************************************
* Routine for the first trapezoidal integration, xi1 *
* n1 keeps track of the number of intervals.         *
*****************************************************}
PROCEDURE Trapez1(VAR xi1:double);
Label 100,150,200,fin;
Var d:double; j1:integer;
Begin
  xi1:=0 ; n1 := 0 ; xx:=x1;
  {Call interpolation subroutine}
  Interpol_Akima;
  {Is there at least one table interval?}
  if x2>X[i+1] then goto 100;
  {If not, integral is simple}
  n1:=n1+1 ; d:=yy ; xx:=x2;
  {Find end point yy value}
  Interpol_Akima;
  xi1:=(yy+d)*(x2-x1)/2.0;
  goto fin;
  {At least one table interval must be summed over}
100: j1:=i;
  xi1:=xi1+(yy+Y[i+1])*(X[i+1]-xx)/2.0;
  {Any more intervals? If not, finish integral with end point}
150: if x2<X[j1+3] then goto 200;
  {Otherwise, keep summing}
  n1:=n1+1;
  xi1:=xi1+(Y[j1+1]+Y[j1+3])*(X[j1+3]-X[j1+1])/2.0;
  j1:=j1+2;
  goto 150;
200: xx:=x2;
  {Find last yy value}
  Interpol_Akima;
  xi1:=xi1+(yy+Y[j1+1])*(x2-X[j1+1])/2.0;
  n1:=n1+1;
fin: End;

{***** integration for xi2 *****}
PROCEDURE Trapez2(VAR xi2:double);
Label 600,650,700,fin;
Var d:double; j1:integer;
Begin
  xi2:=0 ; xx:=x1;
  Interpol_Akima;
  d:=yy;
  if x2>X[i+1] then goto 600;
  xx:=x1+(x2-x1)/2.0;
  Interpol_Akima;
  xi2:=xi2+(d+yy)*(x2-x1)/4.0;
  goto fin;
600: xx:=x1+(X[i+1]-x1)/2.0;
  j1:=i;
  Interpol_Akima;
  xi2:=xi2+(yy+d)*(xx-x1)/2.0;
  xi2:=xi2+(yy+Y[j1+1])*(X[j1+1]-xx)/2.0;
650: if x2<X[j1+2] then goto 700;
  xi2:=xi2+(Y[j1+1]+Y[j1+2])*(X[j1+2]-X[j1+1])/2.0;
  j1:=j1+1;
  goto 650;
700: xx:=x2-(x2-X[j1+1])/2.0;
  Interpol_Akima;
  d:=yy;
  xi2:=xi2+(Y[j1+1]+d)*(x2-xx)/2.0;
  xx:=x2;
  Interpol_Akima;
  xi2:=xi2+(d+yy)*(x2-X[j1+1])/4.0;
fin: End;


{main program}
BEGIN
  iv:=14;  {Number of pooints in table}
  clrscr;  {Clear screen}
  writeln;
  writeln(' Integration of SIN(X)');
  writeln;
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
  {---------------------------------------------------------------}
  write(' Input start point  end point: '); read(x1,x2);
  writeln;

  Integrale;  {Call integral procedure}

  writeln(' Integral from ',x1:6:3,' to ',x2:6:3,' equals ',zz:10:7);
  writeln(' Error check: ',z1:1);
  writeln; ReadKey; DoneWinCrt
END.

{End of file integra.pas}