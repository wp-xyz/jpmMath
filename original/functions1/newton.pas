{***************************************************
*   Program to demonstrate Newton interpolation    *
*     of Function SIN(X) in double precision       *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version By J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  Newton interpolation of SIN(X):                 *
*                                                  *
*    Input X: 0.5                                  *
*                                                  *
*    SIN(X) = 0.47961461 (exact value: 0.47942554) *
*    Error estimate: -0.00000028                   *
*                                                  *
***************************************************}
PROGRAM Newton;
Uses WinCrt;

Label   100;

CONST
        SIZE = 15;

TYPE
        Tab = Array[0..SIZE] of DOUBLE;

VAR
        X,Y,Y1,Y2,Y3 : Tab;
        e,xx,yy : DOUBLE;
        iv, n : integer;


{*******************************************************
* Newton divided differences interpolation subroutine  *
* ---------------------------------------------------- *
* Calculates cubic interpolation for a given table.    *
* iv is the total number of table values. X(i), Y(i)   *
* are the coordinate table values, Y(i) being the      *
* dependant variable. The X(i) may be arbitrarily spa- *
* ced. xx is the interpolation point which is assumed  *
* to be in the interval  with at least one table value *
* to the left, and 3 to the right. If this is violated,*
* n will be set to zero. It is assumed that the table  *
* values are in ascending X(i) order. e is the error   *
* estimate. The returned value is yy.                  *
*******************************************************}
PROCEDURE Interpol_Newton;
Label 100,200,300,fin;
Var i:integer; a,b,c:DOUBLE;
Begin
  {Check to see if interpolation point is correct}
  n:=1;
  if xx>=X[1] then goto 100;
  n:=0 ; goto fin;
100: if xx<=X[iv-3] then goto 200;
  n:=0 ; goto fin;
  {Generate divided differences}
200: for i:=1 to iv-1 do
    Y1[i]:=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
  for i:=1 to iv-2 do
    Y2[i]:=(Y1[i+1]-Y1[i])/(X[i+1]-X[i]);
  for i:=1 to iv-3 do
    Y3[i]:=(Y2[i+1]-Y2[i])/(X[i+1]-X[i]);
  {Find relevant table interval}
  i:=0;
300: i:=i+1;
  if xx>X[i] then goto 300;
  i:=i-1;
  {Begin interpolation}
  a:=xx-X[i];
  b:=a*(xx-X[i+1]);
  c:=b*(xx-X[i+2]);
  yy:=Y[i]+a*Y1[i]+b*Y2[i]+c*Y3[i];
  {Calculate next term in the expansion for an error estimate}
  e:=c*(xx-X[i+3])*yy/24.0;
fin: end;


{main^program}
BEGIN
  iv:=14;
  clrscr;
  writeln;
  writeln(' Newton interpolation of SIN(X):');
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

  Interpol_Newton;

  if n<>0 then
  begin
    writeln;
    writeln('   SIN(X) = ',yy:10:8,' (exact value: ',sin(xx):10:8,')');
    writeln('   Error estimate: ',e:11:8);
    goto 100
  end;
  writeln;
  DoneWinCrt

END.

{End of file newton.pas}