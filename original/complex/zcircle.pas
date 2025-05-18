{******************************************************
* Program to demonstrate the zero searching algorithm *
* --------------------------------------------------- *
*  Reference; BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*               Pascal Version By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
* What is the initial guess for x and y (x0,y0):      *
*                                                     *
*    X0 = 0                                           *
*    Y0 = 0                                           *
*                                                     *
* What is the radius of the first search circle: 5    *
*                                                     *
* By what fraction this circle will be reduced on     *
* each iteration: .5                                  *
*                                                     *
* How many evaluation points per quadrant: 4          *
*                                                     *
* Maximum number of iterations: 10                    *
*                                                     *
* The approximate solution is:                        *
* Z =  0.000576 +  0.999857 I                         *
*                                                     *
* Number of iterations:  10                           *
*                                                     *
******************************************************}
PROGRAM DEMO_ZCIRCLE;
Uses WinCrt;

Var     {global variables}
        e,w,x0,y0,xx,yy : DOUBLE;
        k,m,n   : INTEGER;
        z       : STRING[3];


{*******************************************
  Complex function(x,y) subroutine         }
PROCEDURE Eval(x,y:DOUBLE;VAR u,v:DOUBLE);
Begin    
  u:=x*x-y*y+1.0;
  v:=2.0*x*y
End;
{******************************************}


{************************************************
*        Complex root search subroutine         *
* --------------------------------------------- *
* This routine searches for the complex roots   *
* of an analytical function by encircling the   *
* zero and estimating where it is. The circle   *
* is subsequently tightened by a factor e, and  *
* a new estimate made.                          *
* The inputs to the subroutine are;             *
*     w      initial radius of search circle    *
*     x0,y0  the initial guess                  *
*     e      factor by which circle is reduced  *
*     n      maximum number of iterations       *
*     m      evaluation points per quadrant     *
* The routine returns Z=X+IY (X,Y), and the     *
* number of iterations, k.                      *
************************************************}
PROCEDURE Zcircle;
Label 50,100,200,250,300,310,320,330,340,350,400,450,460,470,500,510,520,521,fin;
Var i,j,m3,m4 : INTEGER;
    X,Y,U,V : Array[1..80] of DOUBLE;
    a,b1,b2,x1,x2,xm1,xm2 : DOUBLE;
    x3,y3,x4,y1,y2,y4 : DOUBLE;
    uu,vv : DOUBLE;

  {Local procedures}
  Procedure Set_to_zero;
  Begin
    xx:=0 ; yy:=0
  End;

  Procedure Auxilliary;
  Begin
    i:=i+1; j:=j-1;
    if i>m then i:=i-m;
    if j<1 then j:=j+m
  End;

  Procedure Interpolation;
  Begin
    j:=i-1;
    if j<1 then j:=j+m;
    {Regula Falsi interpolation for the zero}
    xx:=(X[i]*U[j]+X[j]*U[i])/(U[i]+U[j]);
    yy:=(Y[i]*V[j]+Y[j]*V[i])/(V[i]+V[j])
  End;

Begin  {Zcircle}
  m:=m*4 ; k:=1;
  a:=2*PI/m;
50: for i:=1 to m do
  begin
    X[i]:=w*COS(a*(i-1))+x0;
    Y[i]:=w*SIN(a*(i-1))+y0
  end;
  {Determine the corresponding U[i] and V[i] }
  for i:=1 to m do
  begin
    xx:=X[i] ; yy:=Y[i];
    Eval(xx,yy,uu,vv);
    U[i]:=uu;
    V[i]:=vv
  end;
  {Find the position at which uu changes sign
   in the counterclockwise direction         }
  i:=1 ; uu:=U[i];
100: Auxilliary;
  if uu*U[i]<0 then goto 200;
  {Guard against infinite loop}
  if i=1 then Set_to_zero;
  goto 100;
  {Transition found}
200: xm1:=i;
  {Search for the other transition, starting
   on the other side of the circle          }
  i:=Round(xm1+m DIV 2);
  if i>m then i:=i-m;
  j:=i;
  uu:=U[i];
  {Flip directions alternately}
250: Auxilliary;
  if uu*U[i]<0 then goto 300;
  if uu*U[j]<0 then goto 310;
  {Guard against infinite loop}
  if i=xm1+(xm2/2) then Set_to_zero;
  if j=xm1+(xm2/2) then Set_to_zero;
  goto 250;
  {Transition found}
300: m3:=i;
  goto 320;
  {Transition found}
310: if j=m then j:=0;
  m3:=j+1;
  {xm1 and m3 have been determined, now for xm2 and m4,
   now for the vv transitions                         }
320: i:=Round(xm1+m DIV 4);
  if i>m then i:=i-m;
  j:=i; vv:=V[i];
330: Auxilliary;
  if vv*V[i]<0 then goto 340;
  if vv*V[j]<0 then goto 350;
  {Guard against infinite loop}
  if i=xm1+m DIV 4 then Set_to_zero;
  if j=xm1+m DIV 4 then Set_to_zero;
  goto 330;
  {xm2 has been found}
340: xm2:=i;
  goto 400;
350: if j=m then j:=0;
  xm2:=j+1;
  {xm2 has been found, now for m4}
400: i:=Round(xm2)+(m DIV 2);
  if i>m then i:=i-m;
  j:=i; vv:=V[i];
450: Auxilliary;
  if uu*V[i]<0 then goto 460;
  if vv*V[j]<0 then goto 470;
  {Guard against infinite loop}
  if i=xm2+m DIV 2 then Set_to_zero;
  if j=xm2+m DIV 2 then Set_to_zero;
  goto 450;
460: m4:=i;
  goto 500;
470: if j=m then j:=0;
  m4:=j+1;
  {All the intersections have been determined
   Interpolate to find the four (x,y) coordinates }
500: i:=Round(xm1);
  Interpolation;
  x1:=xx ; y1:=yy ; i:=Round(xm2);
  Interpolation;
  x2:=xx ; y2:=yy ; i:=m3;
  Interpolation;
  x3:=xx ; y3:=yy ; i:=m4;
  Interpolation;
  x4:=xx ; y4:=yy;
  {Calculate the intersection of the lines
   Guard against a divide by zero         }
  if x1<>x3 then goto 510;
  xx:=x1; yy:=(y1+y3)/2.0;
  goto 520;
510: xm1:=(y3-y1)/(x3-x1);
  if x2<>x4 then goto 520;
  xm2:=1e8;
  goto 521;
520: xm2:=(y2-y4)/(x2-x4);
521: b1:=y1-xm1*x1;
  b2:=y2-xm2*x2;
  xx:=-(b1-b2)/(xm1-xm2);
  yy:=(xm1*b2+xm2*b1)/(xm1+xm2);
  {is another iteration in order ? }
  if k=n then goto fin;
  x0:=xx ; y0:=yy ; k:=k+1 ; w:=w*e;
  goto 50;
fin: End; {Zcircle}


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' What is the initial guess for x and y (x0,y0): ');
  writeln;
  write('    X0 = '); read(x0);
  write('    Y0 = '); read(y0);
  writeln;
  write(' What is the radius of the first search circle: '); read(w);
  writeln;
  write(' By what fraction this circle will be reduced on each iteration: '); read(e);
  writeln;
  write(' How many evaluation points per quadrant: '); read(m);
  writeln;
  write(' Maximum number of iterations: '); read(n);
  writeln;
  
  Zcircle;     {Call ZCircle subroutine}

  writeln;
  if yy>=0 then z:=' + ';
  if yy< 0 then z:=' - ';
  writeln(' The approximate solution is:');
  writeln(' Z = ',xx:9:6,z,yy:9:6,' I');
  writeln;
  writeln(' Number of iterations: ',k);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file zcircle.pas}















