{*******************************************************************************
*                         THE   CYCLOIDAL  CURVES                              *
* ---------------------------------------------------------------------------- *
* Explanations:                                                                *
*                                                                              *
* Here are some definitions concerning the cycloidal family:                   *
*                                                                              *
* Cycloids: this curve is described as the trajectory in 2D plane of a given   *
* circle point, when the circle rolls on a straight line. An extended or       *
* shortened cycloid is called trochoid. The parametric equations of a cycloid  *
* are:                                                                         *
*                   X = R.(T - sin(T))                                         *
*                   Y = R.(1 - cos(T))                                         *
*                                                                              *
* Here we do not study these classical curves but their sisters:               *
*                                                                              *
* Epicycloid: this curve is described as the trajectory in 2D plane of a       *
* given circle point, when the circle rolls on an outter fixed circle.         *
*                                                                              *
* Hypocycloid: this curve is described as the trajectory in 2D plane of a      *
* given circle point, when the circle rolls on an inner fixed circle.          *
*                                                                              *
* An extended or shortened Epicycloid (Hypocycloid) is called Epitrochoid      *
* (Hypotrochoid).                                                              *
*                                                                              *
* Here the general equations of these curves:                                  *
*                                                                              *
*          X = (R1+R2).cos(T) - L.R2.cos(((R1+R2)/R2).T)                       *
*          Y = (R1+R2).sin(T) - L.R2.sin(((R1+R2)/R2).T)                       *
*                                                                              *
*          R1 = radius of the fixed circle                                     *
*          R2 = radius of the mobile circle                                    *
*           (R2>0 --> Epicycloid, R2<0 --> Hypocycloid)                        *
*          L is a parameter, such as:                                          *
*            L > 1, we get an extended curve,                                  *
*            L = 1, we get a normal curve,                                     *
*            0 < L < 1, we get a shortened curve.                              *
*                                                                              *
* Some examples:    R1=10, R2=4,  L=1                                          *
*                   R1=10, R2=-4, L=1                                          *
*                   R1=10, R2=6,  L=2                                          *
*                   R1=10, R2=-6, L=1                                          *
*                   R1=10, R2=4.5, L=0.6                                       *
*                   R1=10, R2=-4.5, L=0.6                                      *
*                   R1=5,  R2=-3, L=0.5                                        *
*                   R1=17, R2=2,  L=1                                          *
*                   R1=10, R2=17, L=1                                          *
*                   R1=17, R2=2,  L=1                                          *
* ---------------------------------------------------------------------------- *
* Reference:                                                                   *
*    From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0 de    *
*    R. Dony - MASSON 1990 page 113" [BIBLI 12].                               *
*                                                                              *
*                                          TPW Release By J-P Moreau, Paris.   *
*                                                 (www.jpmoreau.fr)            *
*******************************************************************************}
Program Cyclo;  {Draw cycloidal curves}

Uses WinCrtMy, WinProcs, Strings, Type_def, CrtGr2D;

Var
    aux1,aux2,aux3: REAL;
    R1,R2,L: REAL;
    MaxRadius: REAL;
    Period: REAL;
    Generators: char;


    Function FX(t:real): real;
    begin
      FX:=aux1*cos(t) - aux3*cos(aux2*t)
    end;

    Function FY(t:real): real;
    begin
      FY:=aux1*sin(t) - aux3*sin(aux2*t)
    end;

    Procedure SeekPeriod;
    Var Xstart, Ystart: real;
        CurrentX, CurrentY: real;
        deltax,deltay: real;
    Begin
      Period:=0;
      Xstart:=FX(0);
      Ystart:=FY(0);
      Repeat
        Period:=Period + 2*PI;
        CurrentX:=FX(Period);
        CurrentY:=FY(Period);
        deltax:=abs(CurrentX-Xstart);
        deltay:=abs(CurrentY-Ystart);
      Until ((deltax<1E-6) and (deltay<1E-6)) or (Period>200*PI)
    End;

    Procedure Data;
    Begin
      ClrScr;
      writeln;
      writeln('     Cycloidal Curves');
      writeln('     ----------------');
      writeln;
      writeln('  R1=radius of fixed circle');
      writeln('  R2=radius of mobile circle');
      writeln;
      writeln('  If R2 > 0, we get an Epicycloid');
      writeln('  If R2 < 0, we get an Hypocycloid');
      writeln;
      writeln('  If L > 1, we get an extended curve');
      writeln('  If L = 1, we get a normal curve');
      writeln('  If 0 < L < 1, we get a shortened curve');
      writeln;
      write('       R1 = '); readln(R1);
      write('       R2 = '); readln(R2);
      Repeat
        write('       L = '); readln(L)
      Until L>=0;
      aux1:=R1+R2; aux2:=aux1/R2; aux3:=L*R2;
      writeln;
      write('  Do you want to draw the generator circles (y/n): '); Generators:=ReadKey
    End;

    Procedure SeekWindow;
    Begin
      MaxRadius:=0;
      if L>=1 then
        MaxRadius:=R1+R2 + abs(L*R2)
      else if R2 > 0 then
        MaxRadius:=R1+2*R2
      else
        MaxRadius:=R1
    End;

    Procedure DrawCurve;
    Var a,step: real;
    Begin
      a:=0;
      step:=PI/48.0;
      MoveXY(CrtDc,FX(0),FY(0));
      Repeat
        a:=a+step;
        LineXY(CrtDc,FX(a),FY(a));
      Until a>=Period
    End;

    Procedure GeneratorCircles;
    Begin
      Circle1(CrtDc,0,0,R1,False);
      Circle1(CrtDc,R1+R2,0,R2,False)
    End;

    Procedure DisplayParam;
    Var s1,s2,s3,s4,s5: String;
        ch:Array[0..79] of char;
    Begin
      Str(R1:5:2,s1);
      Str(R2:5:2,s2);
      Str(L:5:2,s3);
      Str(Period:8:5,s4);
      Str(MaxRadius:8:5,s5);
      Cloture(100,MaxX-100,100,111);
      StrPCopy(ch,'R1='+s1+'    R2='+s2+'    L='+s3+'    Period='+s1+'    MaxRadius='+s5);
      TextOut(CrtDc,80,10,ch,strlen(ch))
    End;


{main program}
BEGIN
  WinCrtInit(' Cycloidal Curves');
  Repeat
    Data;
    SeekPeriod;
    SeekWindow;
    ClrScr;
    Fenetre(-MaxRadius,MaxRadius,-MaxRadius,MaxRadius);
    Cloture(100,MaxX-100,100,MaxY-50);
    Bordure(CrtDc);
    Cloture(150,MaxX-150,110,MaxY-60);
    if Upcase(Generators)='Y' then GeneratorCircles;
    DrawCurve;
    DisplayParam;
    SortieGraphique
  Until rep='n';
  DoneWinCrt
END.

{end of file cyclo1.pas}