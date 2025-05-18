{****************************************************
*        Program to demonstrate the complex         *
*             root counting subroutine              *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*              Pascal Version By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)           *
* ------------------------------------------------- *
* Example: Find the number of complex roots of      *
*          F(z) = z^2 + 1                           *
*                                                   *
* SAMPLE RUN:                                       *
* Where is the center of the search circle (x0,y0): *
*                                                   *
*    X0 = 0                                         *
*    Y0 = 0                                         *
*                                                   *
* What is the radius of this circle: 4              *
* How many evaluation points per quadrant: 4        *
*                                                   *
* Number of complete cycles found:  2               *
* Residual:  0                                      *
*                                                   *
****************************************************}
PROGRAM Num_Complex_Roots;
Uses WinCrt;

VAR
        i,m,nn : INTEGER;
        N : Array[1..80] of INTEGER;
        a,w,x0,y0 : DOUBLE;


{*******************************************
  Complex function(x,y) subroutine         }
PROCEDURE Eval(x,y:DOUBLE;VAR u,v:DOUBLE);
Begin    
  u:=x*x-y*y+1.0;
  v:=2.0*x*y
End;
{******************************************}


{************************************************
*       Complex root counting subroutine        *
* --------------------------------------------- *
* This routine calculates the number of complex *
* roots within a circle of radius w centered on *
* (x0,y0) by counting (u,v) transitions around  *
* the circumference. The input parameters are:  *
*     w      radius of the circle               *
*     x0,y0  center of the circle               *
*     m      evaluation points per quadrant     *
* The routine returns nn, the number of roots   *
* found, and a, where a<>0 indicates a failure. *
************************************************}
PROCEDURE RootNum;
Label 100;
Var i:INTEGER;
    x,y,u,v:DOUBLE;
Begin
  a:=PI/(2.0*m);
  {Establish the N[i] array}
  for i:=1 to 4*m do
  begin
    x:=w*COS(a*(i-1))+x0;
    y:=w*SIN(a*(i-1))+y0;
    Eval(x,y,u,v);  {Call complex function subroutine}
    if u>=0 then if v>=0 then N[i]:=1;
    if u< 0 then if v>=0 then N[i]:=2;
    if u< 0 then if v< 0 then N[i]:=3;
    if u>=0 then if v< 0 then N[i]:=4
  end;
  {Count complete cycles counterclockwise}
  nn:=N[1]; a:=0.0;
  for i:=2 to 4*m do
  begin
    if nn=N[i] then goto 100;
    if nn<>4 then if nn=N[i]+1 then a:=a-1;
    if nn= 1 then if N[i]=4 then a:=a-1;
    if nn=4  then if N[i]=1 then a:=a+1;
    if nn+1=N[i] then a:=a+1;
    nn:=N[i];
100: end;
  {Complete circle}
  if nn<>4 then if nn=N[1]+1 then a:=a-1;
  if nn= 1 then if N[1]=4 then a:=a-1;
  if nn=4  then if N[1]=1 then a:=a+1;
  if nn+1=N[1] then a:=a+1;
  a:=ABS(a); nn:=Round(a/4);
  a:=a-4.0*INT(a/4)
End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Where is the center of the search circle (x0,y0): ');
  writeln;
  write('    X0 = '); read(x0);
  write('    Y0 = '); read(y0);
  writeln;
  write(' What is the radius of this circle: '); read(w);
  writeln;
  write(' How many evaluation points per quadrant: '); read(m);
  if m>20 then m:=20;
  writeln;
  
  RootNum;   {Call RootNum subroutine}

  writeln;
  writeln(' Number of complete cycles found: ',nn);
  writeln;
  if a<>0.0 then
    writeln(' Residual: ',a)
  else
    writeln(' Residual: 0');
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file rootnum.pas}













