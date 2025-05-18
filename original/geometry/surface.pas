{****************************************************
* This program calculates the internal surface of a *
* polygon defined by n points of coordinates:       *
*         x1, y1, x2, y2,...xn, yn.                 *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* npoints=13                                        *
* x1=3    y1=3                                      *
* x2=12   y2=3                                      *
* x3=12   y3=6                                      *
* x4=11   y4=6                                      *
* x5=11   y5=7                                      *
* x6=15   y6=7                                      *
* x7=15   y7=9                                      *
* x8=8    y8=9                                      *
* x9=8    y9=6                                      *
* x10=4   y10=6                                     *
* x11=4   y11=9                                     *
* x12=3   y12=9                                     *
* x13=3   y13=3                                     *
*                                                   *
* Surface is: 4.70000000000000E+0001                *
*                                                   *
****************************************************}
Program Surface;
Uses WinCrt;

Const   NMAX=50;

Var  s: REAL;
     i,npoints: INTEGER;
     buffer: Array[1..NMAX,1..NMAX] of REAL;

BEGIN
  s:=0;
  npoints:=13;                    {13 points}
  {init buffer x(i), y(i) }
  for i:=1 to npoints do
  begin
    buffer[1,1]:= 3; buffer[1,2]:=3;
    buffer[2,1]:=12; buffer[2,2]:=3;
    buffer[3,1]:=12; buffer[3,2]:=6;
    buffer[4,1]:=11; buffer[4,2]:=6;
    buffer[5,1]:=11; buffer[5,2]:=7;
    buffer[6,1]:=15; buffer[6,2]:=7;
    buffer[7,1]:=15; buffer[7,2]:=9;
    buffer[8,1]:= 8; buffer[8,2]:=9;
    buffer[9,1]:= 8; buffer[9,2]:=6;
    buffer[10,1]:=4; buffer[10,2]:=6;
    buffer[11,1]:=4; buffer[11,2]:=9;
    buffer[12,1]:=3; buffer[12,2]:=9;
    buffer[13,1]:=3; buffer[13,2]:=3;
  end;
  {calculate internal surface}
  for i:=1 to npoints-1 do
    s:=s+(buffer[i,1]+buffer[i+1,1])*(buffer[i,2]-buffer[i+1,2]);
  s:=abs(s)/2.0;
  {print surface}
  writeln;
  writeln(' Surface is: ',s);
  writeln;
  Readkey; Donewincrt
END.

{end of file surface.pas}