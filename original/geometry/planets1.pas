{****************************************************
*              Position of planets                  *
* ------------------------------------------------- *
* This program calculates the ephemeris of a planet *
* in our solar system to adjust an equatorial tele- *
* scope.                                            *
* ------------------------------------------------- *
* Ref.: "Mathematiques par l'informatique indivi-   *
*        duelle - Programmes en BASIC, MASSON,      *
*        Paris, 1982, p 100 - 105" [BIBLI 06].      *
* ------------------------------------------------- *
* Inputs:                                           *
*   Date: day,month,year                            *
*   Hour UT: hour                                   *
*   Planet number: 1 to 8                           *
*   (Mercury: 1 Venus : 2 Mars   : 4 Jupiter: 5     *
*    Saturn : 6 Uranus: 7 Neptune: 8)               *
* Outputs:                                          *
*   Ascent in hours (0 to 24)                       *
*   Declination in deg. (-90 to 90 North or South)  *
*                                                   *
* SAMPLE RUN:                                       *
* (Find ascent and declination of planet Mars on    *
*  March 10th, 1982 at 6h UT)                       *
*                                                   *
* Date (D,M,Y): 10 3 1982                           *
* Hour UT: 6                                        *
* Mercury: 1 Venus : 2 Mars   : 4 Jupiter: 5        *
* Saturn : 6 Uranus: 7 Neptune: 8                   *
* Planet number: 4                                  *
*                                                   *
* Ascent     : 13 H  8 MN                           *
* Declination:  3 ° 45 MN S                         *
*                                                   *
* This program allows building the following table: *
*                                                   *
* Date: March 10th, 1982 at 6H UT.                  *
*                                                   *
*       Planet     Ascent     Declination           *
*       Mercury    21 H 51    14 ° 45 mn S          *
*       Venus      20 H 26    14 ° 57 mn S          *
*       Mars       13 H 08     3 ° 45 mn S          *
*       Jupiter    14 H 32    13 ° 30 mn S          *
*       Saturn     13 H 22     5 ° 42 mn S          *
*                                                   *
*              Pascal Version By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)           *
****************************************************}
PROGRAM PLANETS;
Uses WinCrt;

Type Temp = Array[1..9] of double;

Var
     A: Array[1..9,1..8] of double;
     B: Array[1..12] of double;
     temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8: Temp;
     aa,bb,c,d,g,h,t,x,y,z: double;
     i,j,m,n: integer;
     As:char;

Procedure planet_coordinates(j:Integer);
Var xl,xm,o,p,q,e,xi,xa,xj,r,u,v: double;
Begin
{calculate planetary constants}
  XL:=A[1,j]+A[2,j]*t;
  O:=A[3,j]+A[4,j]*t;
  P:=A[5,j]+A[6,j]*t;
  E:=A[7,j];
  XI:=A[8,j];
  XA:=A[9,j];
{solve Kepler's equation}
  xm:=XL-O; u:=xm;
  for j:=1 to 10 do u:=xm+E*sin(u);
{formula (3) of reference book}
 r:=XA*(cos(u)-E);
 v:=XA*SQRT(1.0-E*E)*sin(u);
 O:=O-P; xm:=sin(O); O:=cos(O);
 q:=sin(P); P:=cos(P);
 xj:=sin(XI); XI:=cos(XI);
 x:=(P*O-XI*q*xm)*r+(-P*xm-XI*q*O)*v;
 y:=(q*O+XI*P*xm)*r+(-q*xm+XI*P*O)*v;
 z:=xj*xm*r+xj*O*v
End;

Function ATAN2(y,x:Double):Double;
Var z:Double;
Begin
  if ABS(x) < 1e-12 then
  begin
    if y>0 then
      z:=PI/2.0
    else if ABS(y) < 1e-12 then
      z:=-PI/2.0
    else
      z:=0.0
  end
  else
  begin
    z:=ARCTAN(y/x);
    if x<0.0 then z:=z+PI
  end;
  ATAN2:=z
End;

{main program}
BEGIN

  {Init planetary constants (fill table A by columns) }
  temp1[1]:=4.01166; temp1[2]:=0.071425454; temp1[3]:=1.32493; temp1[4]:=0.000000742289;
  temp1[5]:=0.823045; temp1[6]:=0.000000566185; temp1[7]:=0.205615; temp1[8]:=0.122225;
  temp1[9]:=0.387099;
  temp2[1]:=3.60861; temp2[2]:=0.027963119; temp2[3]:=2.271616; temp2[4]:=0.00000065572;
  temp2[5]:=1.32291; temp2[6]:=0.000000436681; temp2[7]:=0.006816; temp2[8]:=0.0592301;
  temp2[9]:=0.723332;
  temp3[1]:=1.72727; temp3[2]:=0.0172028; temp3[3]:=1.76688; temp3[4]:=0.000000818559;
  temp3[5]:=0.0; temp3[6]:=0.0; temp3[7]:=0.016751; temp3[8]:=0.0;
  temp3[9]:=1.0;
  temp4[1]:=2.17756; temp4[2]:=0.0091467658; temp4[3]:=5.83378; temp4[4]:=0.000000879297;
  temp4[5]:=0.851616; temp4[6]:=0.000000371232; temp4[7]:=0.093309; temp4[8]:=0.0322939;
  temp4[9]:=1.5236;
  temp5[1]:=4.68279; temp5[2]:=0.00145099; temp5[3]:=0.2289; temp5[4]:=0.000000857;
  temp5[5]:=1.73578; temp5[6]:=0.000000482933; temp5[7]:=0.048376; temp5[8]:=0.0228418;
  temp5[9]:=5.202799;
  temp6[1]:=4.8567; temp6[2]:=0.00058484; temp6[3]:=1.5974; temp6[4]:=0.000000412;
  temp6[5]:=1.96856; temp6[6]:=0.000000417308; temp6[7]:=0.054311; temp6[8]:=0.0435026;
  temp6[9]:=9.552098;
  temp7[1]:=4.3224; temp7[2]:=0.000205424; temp7[3]:=2.9523; temp7[4]:=0.000000762;
  temp7[5]:=1.2825; temp7[6]:=0.000000238237; temp7[7]:=0.047319; temp7[8]:=0.013482;
  temp7[9]:=19.21694;
  temp8[1]:=1.5223; temp8[2]:=0.000105061; temp8[3]:=0.7637; temp8[4]:=0.000000393;
  temp8[5]:=2.28102; temp8[6]:=0.00000052517; temp8[7]:=0.008262; temp8[8]:=0.0310536;
  temp8[9]:=30.12912;

{init calendar constants}
  B[1]:=0.0;   B[2]:=31.0; B[3]:=59.0;  B[4]:=25.0;  B[5]:=90.0;   B[6]:=25.0;
  B[7]:=120.0; B[8]:=25.0; B[9]:=151.0; B[10]:=25.0; B[11]:=181.0; B[12]:=25.0;

for i:=1 to 9 do
begin
  A[i,1]:=temp1[i];
  A[i,2]:=temp2[i];
  A[i,3]:=temp3[i];
  A[i,4]:=temp4[i];
  A[i,5]:=temp5[i];
  A[i,6]:=temp6[i];
  A[i,7]:=temp7[i];
  A[i,8]:=temp8[i]
end;

{print_mat;}

writeln;
{input hour, date (DD,MM,YY) and planet number}
write(' Date (D,M,Y): '); readln(j, m, aa);
write(' Hour UT: '); readln(h);
writeln(' Mercury: 1 Venus : 2 Mars   : 4 Jupiter: 5');
writeln(' Saturn : 6 Uranus: 7 Neptune: 8');
write(' Planet number: '); readln(n);

{calculate time t}
t:=365.25*(aa-1901)+B[m]+j;
t:=INT(t) + h/24.0;
{calculate earth coordinates}
planet_coordinates(3);                   {planet #3 coordinates}
g:=x; h:=y;                             {save earth coordinates}
{calculate coordinates of planet #n}
planet_coordinates(n);
{calculate geocentric equatorial coordinates}
x:=x-g; y:=y-h;
t:=y*0.917484-z*0.397772; 
z:=y*0.397772+z*0.917484;
y:=t;
{calculate ascent and declination}
aa:=ATAN2(y,x);
d:=ATAN2(z,SQRT(x*x+y*y));
{conversion}
aa:=aa*12.0/PI;
if aa<0.0 then aa:=24.0+aa;
h:=INT(aa); m:=Round(INT(60*(aa-h)));
d:=d*180.0/PI;
As:='N'; if d<0.0 then As:='S';
d:=ABS(d); bb:=INT(d); c:=INT(60*(d-bb));
{print results}
writeln;
writeln(' Ascent     : ',h:2:0,' H ',m:2,' MN');
writeln(' Declination: ',INT(bb):2:0,' ° ',INT(c):2:0,' MN ',As);
writeln;
Readkey; Donewincrt

END.

{end of file planets1.pas}