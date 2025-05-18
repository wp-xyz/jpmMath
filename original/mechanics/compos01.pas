{***********************************************************************
*  This program calculates the stresses in a unidirectional layer of a *
*  composite material, knowing deformations exx, eyy, gxy, and angle   *
*  theta of fibers in x direction. We consider here that we have a     *
*  plane stress state.                                                 *
*                                                                      *
*  The stresses sxx, syy and sxy are calculated:                       *
*                                                                      *
*    1. in plane (x,y),                                                *
*    2. in the main axes (L,W) of the layer.                           *
*                                                                      *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                                                                      *
*   STRESSES IN A UNIDIRECTIONAL COMPOSITE LAYER                       *
*                                                                      *
*      Angle of fibers:  30.000 °                                      *
*      eps. xx        :  10.000 mm                                     *
*      eps. yy        :  -5.000 mm                                     *
*      gam. xy        :  20.000 mm                                     *
*                                                                      *
*   REDUCED STIFFNESS MATRIX (MAIN AXES):                              *
*                                                                      *
*     41.0509     3.2841     0.0000                                    *
*      3.2841    10.2627     0.0000                                    *
*      0.0000     0.0000     4.5000                                    *
*                                                                      *
*   REDUCED STIFFNESS MATRIX (X, Y AXES):                              *
*                                                                      *
*     28.3391     8.2989     9.5611                                    *
*      8.2989    12.9450     3.7706                                    *
*      9.5611     3.7706     9.5148                                    *
*                                                                      *
*   Stresses in axes x, y:                                             *
*                                                                      *
*    433.1189 MPa                                                      *
*     93.6746 MPa                                                      *
*    267.0540 MPa                                                      *
*                                                                      *
*   Stresses in main axes:                                             *
*                                                                      *
*    579.5334 MPa                                                      *
*    -52.7399 MPa                                                      *
*    -13.4567 MPa                                                      *
*                                                                      *
* -------------------------------------------------------------------- *
*                                Pascal Release By J-P Moreau, Paris.  *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
PROGRAM COMPOS01;

Uses WinCrt;

Const
      NMAX=3;
Type
      MAT = Array[1..NMAX,1..NMAX] of Double;
      VEC = Array[1..NMAX] of Double;

Var
      EL,ET,NULT,GLT,TEMP,TH,EXX,EYY,GXY: Double;
      I:Integer;

      Q, Q1   : MAT;
      S, S1, E: VEC;


  Procedure MATNUL(N:Integer; A:MAT);
  Var I,J:Integer;
  Begin
    For I:=1 to N do
      for J:=1 to N do
        A[I,J]:=0.0
  End;

  Procedure MATMULT(N:Integer; A:MAT; B:VEC; Var C:VEC);  
  var SUM:Double; I,K:Integer;
  Begin
    for I:=1 to N do
    begin
      SUM:= 0.0;
      for K:=1 to N do
        SUM:=SUM+A[I,K]*B[K];
      C[I]:=SUM
    end
  End;

  Procedure MATPRINT(title:String; N:Integer; A:MAT);
  Var I,J:Integer; 
  Begin
    writeln;
    writeln(title);
    writeln;
    for I:=1 to N do
    begin
      for j:=1 to N do write(' ',A[I,J]:10:4);
      writeln
    end
  End;

  FUNCTION dcos4(t:Double): Double; 
  Var a: Double;
  Begin
    a:=cos(t);
    dcos4 := a*a*a*a
  End;

  FUNCTION dsin4(t:Double): Double;
  Var a: Double;
  Begin
    a:=sin(t);
    dsin4 := a*a*a*a
  End;


{main program}
BEGIN
{ material constants (GPa) } 

      EL := 40.0;
      ET := 10.0;
      NULT :=  0.32;
      GLT  :=  4.5;

{ initial data } 

      TH  := 30.0;
      EXX := 0.01;
      EYY :=-5E-03;
      GXY := 2E-02;

      writeln;
      writeln('  STRESSES IN A UNIDIRECTIONAL COMPOSITE LAYER');
      writeln;
      writeln('    Angle of fibers: ',TH:7:3,' °');
      writeln('    eps. xx        : ',1000*EXX:7:3,' mm');
      writeln('    eps. yy        : ',1000*EYY:7:3,' mm');   
      writeln('    gam. xy        : ',1000*GXY:7:3,' mm');   

      ReadKey;

{   reduced stiffness matrix in main axes (L,W) }
      MATNUL(3,Q);
      Q[1,1] := EL/(1.0-((ET/EL)*NULT*NULT));
      Q[2,2] := (ET/EL)*Q[1,1];
      Q[1,2] := NULT*Q[2,2];
      Q[2,1] := Q[1,2];
      Q[3,3] := GLT;

      MATPRINT('  REDUCED STIFFNESS MATRIX (MAIN AXES):',3,Q);

{   reduced stiffness matrix in axes (x,y) }
      MATNUL(3,Q1);
      TH := 30.0*PI/180.0;
      TEMP := sin(TH)*sin(TH)*cos(TH)*cos(TH);
      Q1[1,1]:=Q[1,1]*dcos4(TH)+Q[2,2]*dsin4(TH)+2.0*(Q[1,2]+2.0*Q[3,3])*TEMP;
      Q1[1,2]:=(Q[1,1]+Q[2,2]-4.0*Q[3,3])*TEMP+Q[1,2]*(dcos4(TH)+dsin4(TH));
      Q1[2,1]:=Q1[1,2];
      TEMP := sin(TH)*cos(TH)*cos(TH)*cos(TH);
      Q1[1,3]:=(Q[1,1]-Q[1,2]-2.0*Q[3,3])*TEMP;
      TEMP := sin(TH)*sin(TH)*sin(TH)*cos(TH);
      Q1[1,3]:=Q1[1,3]+(Q[1,2]-Q[2,2]+2.0*Q[3,3])*TEMP;
      Q1[3,1]:=Q1[1,3];
      TEMP := sin(th)*sin(th)*cos(th)*cos(th);
      Q1[2,2]:=Q[1,1]*dsin4(th)+2.0*(Q[1,2]+2.0*Q[3,3])*TEMP+Q[2,2]*dcos4(th);
      TEMP := sin(th)*sin(th)*sin(th)*cos(th);
      Q1[2,3]:=(Q[1,1]-Q[1,2]-2.0*Q[3,3])*TEMP;
      TEMP := cos(th)*cos(th)*cos(th)*sin(th);
      Q1[2,3]:=Q1[2,3]+(Q[1,2]-Q[2,2]+2.0*Q[3,3])*TEMP;
      Q1[3,2]:=Q1[2,3];
      TEMP := sin(th)*sin(th)*cos(th)*cos(th);
      Q1[3,3]:=(Q[1,1]+Q[2,2]-2.0*(Q[1,2]+Q[3,3]))*TEMP+Q[3,3]*(dsin4(th)+dcos4(th));

      MATPRINT('  REDUCED STIFFNESS MATRIX (X, Y AXES):',3,Q1);

{   stresses in axes (x, y): }
      E[1]:=EXX;
      E[2]:=EYY;
      E[3]:=GXY;
      MATMULT(3,Q1,E,S);
      writeln;
      writeln('  Stresses in axes x, y:');
      writeln;
      for I:=1 to NMAX do
        writeln('  ',1000*S[I]:10:4,' MPa');

{   stresses in main axes: } 
      TEMP:=sin(th)*cos(th);
      Q[1,1]:=cos(th)*cos(th);
      Q[1,2]:=sin(th)*sin(th);
      Q[1,3]:=2.0*TEMP;
      Q[2,1]:=Q[1,2];
      Q[2,2]:=Q[1,1];
      Q[2,3]:=-Q[1,3];
      Q[3,1]:=-TEMP;
      Q[3,2]:=TEMP;
      Q[3,3]:=cos(th)*cos(th)-sin(th)*sin(th);
      MATMULT(3,Q,S,S1);
      writeln;
      writeln('  Stresses in main axes:');
      writeln;
      for I:=1 to NMAX do
        writeln('  ',1000*S1[I]:10:4,' MPa');
      writeln;
      ReadKey;
      DoneWinCrt
END.  {of main program}

{end of file Compos01.pas}