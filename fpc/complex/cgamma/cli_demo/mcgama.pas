{*************************************************************
*       Calculate Function Gamma with a complex argument     *
* ---------------------------------------------------------- *
* EXPLANATION:                                               *
* Purpose: This program computes the gamma function G(z)     *
*          or Ln[G(z)] for a complex argument using          *
*          subroutine CGAMA                                  *
* Input :  x  --- Real part of z                             *
*          y  --- Imaginary part of z                        *
*          KF --- Function code                              *
*          KF=0 for Ln[G(z)]                                 *
*          KF=1 for G(z)                                     *
* Output:  GR --- Real part of Ln[G(z)] or G(z)              *
*          GI --- Imaginary part of Ln[G(z)] or G(z)         *
* Examples:                                                  *
*    x         y           Re[G(z)]           Im[G(z)]       *
*  --------------------------------------------------------  *
*   2.50      5.00     .2267360319E-01    -.1172284404E-01   *
*   5.00     10.00     .1327696517E-01     .3639011746E-02   *
*   2.50     -5.00     .2267360319E-01     .1172284404E-01   *
*   5.00    -10.00     .1327696517E-01    -.3639011746E-02   *
*                                                            *
*    x         y          Re[LnG(z)]         Im[LnG(z)]      *
*  --------------------------------------------------------  *
*   2.50      5.00    -.3668103262E+01     .5806009801E+01   *
*   5.00     10.00    -.4285507444E+01     .1911707090E+02   *
*   2.50     -5.00    -.3668103262E+01    -.5806009801E+01   *
*   5.00    -10.00    -.4285507444E+01    -.1911707090E+02   *
* ---------------------------------------------------------- *
* SAMPLE RUNS:                                               *
*                                                            *
* Please enter KF, x and y: 1 2.5 5.0                        *
*                                                            *
*    x         y           Re[G(z)]           Im[G(z)]       *
*  --------------------------------------------------------  *
*   2.50      5.00       0.0226736032      -0.0117228440     *
*                                                            *
* Please enter KF, zx and zy: 0 2.5 5.0                      *
*                                                            *
*    x         y          Re[LnG(z)]         Im[LnG(z)]      *
*  --------------------------------------------------------  *
*   2.50      5.00      -3.6681032624      5.8060098006      *
*                                                            *
* ---------------------------------------------------------- *
* REFERENCE:                                                 *
*    "Fortran Routines for Computation of Special Functions, *
*     jin.ece.uiuc.edu/routines/routines.html".              *
*                                                            *
*                       Pascal Release By J-P Moreau, Paris. *
*                                (www.jpmoreau.fr)           *
*************************************************************}
  PROGRAM MCGAMA;

    Var X,X0,Y,GR,GI: Double;
        KF: Integer;

    {auxiliary functions}
    Function SINH(x:double):double;
    Var expx: double;
    Begin
      expx := exp(x);
      SINH := 0.5*(expx-1.0/expx)
    End;

    Function COSH(x:double):double;
    Var expx: double;
    Begin
      expx := exp(x);
      COSH := 0.5*(expx+1.0/expx)
    End;

    Function Power(x:Double; n:Integer): Double;
    {Calculate x power n}
    var
      i, m: integer;
    begin
      Result := 1.0;
      if n=0 then
        exit;

      m :=  n;
      if n < 0 then m := -n;
      for i:=1 to m do
        result :=x*result;
      if n < 0 then Result := 1.0/result;
    end;


  Procedure CGAMA(X,Y:Double; KF:Integer; Var GR,GI:Double);
{ ===========================================================
        Purpose: Compute the gamma function G(z) or Ln[G(z)]
                 for a complex argument
        Input :  x  --- Real part of z
                 y  --- Imaginary part of z
                 KF --- Function code
                        KF=0 for Ln[G(z)]
                        KF=1 for G(z)
        Output:  GR --- Real part of Ln[G(z)] or G(z)
                 GI --- Imaginary part of Ln[G(z)] or G(z)
  =========================================================== }
  Label Return;
  Var
      A: Array[1..10] of Double;
      G0,GR1,GI1,SR,SI,T,TH,TH1,TH2,X1,Y1,Z1,Z2: Double;
      J,K,NA:Integer;
  Begin
    {Initialize table A}
    A[1]:= 8.333333333333333E-02; A[2]:= -2.777777777777778E-03;
    A[3]:= 7.936507936507937E-04; A[4]:= -5.952380952380952E-04;
    A[5]:= 8.417508417508418E-04; A[6]:= -1.917526917526918E-03;
    A[7]:= 6.410256410256410E-03; A[8]:= -2.955065359477124E-02;
    A[9]:= 1.796443723688307E-01; A[10]:=-1.39243221690590;

    IF (Y = 0.0) AND (X = INT(X)) AND (X <= 0.0) THEN
    begin
       GR:=1E+300;  {arbitrary big number}
       GI:=0.0;
       goto RETURN
    end
    ELSE IF X < 0.0 THEN
    begin
       X1:=X;
       Y1:=Y;
       X:=-X;
       Y:=-Y
    end;
    X0:=X;
    IF X <= 7.0 THEN
    begin
      NA:=Round(7.0-X);
      X0:=X+NA
    end;
    Z1:=SQRT(X0*X0+Y*Y);
    TH:=ArcTan(Y/X0);
    GR:=(X0-0.5)*Ln(Z1)-TH*Y-X0+0.5*Ln(2.0*PI);
    GI:=TH*(X0-0.5)+Y*Ln(Z1)-Y;
    For K:=1 to 10 do
    begin
      T:=Power(Z1,(1-2*K));
      GR:=GR+A[K]*T*COS((2.0*K-1.0)*TH);
      GI:=GI-A[K]*T*SIN((2.0*K-1.0)*TH)
    end;
    IF X <= 7.0 THEN
    begin
      GR1:=0.0;
      GI1:=0.0;
      For J:=0 to NA-1 do
      begin
        GR1:=GR1+0.5*Ln(Sqr(X+J)+Y*Y);
        GI1:=GI1+ArcTan(Y/(X+J))
      end;
      GR:=GR-GR1;
      GI:=GI-GI1
    end;
    IF X1 < 0.0 THEN
    begin
      Z1:=SQRT(X*X+Y*Y);
      TH1:=ArcTan(Y/X);
      SR:=-SIN(PI*X)*COSH(PI*Y);
      SI:=-COS(PI*X)*SINH(PI*Y);
      Z2:=SQRT(SR*SR+SI*SI);
      TH2:=ArcTan(SI/SR);
      IF SR < 0.0 Then TH2:=PI+TH2;
      GR:=Ln(PI/(Z1*Z2))-GR;
      GI:=-TH1-TH2-GI;
      X:=X1;
      Y:=Y1
    end;
    IF KF = 1 THEN
    begin
      G0:=EXP(GR);
      GR:=G0*COS(GI);
      GI:=G0*SIN(GI)
    end;
  Return:
  End;


  {main program}
  BEGIN
    WRITELN;
    WRITE('  Please enter KF, x and y: ');
    READLN(KF, X, Y);
    WRITELN;
    IF KF = 1 THEN
      WRITELN('        x         y           Re[G(z)]           Im[G(z)]')
    ELSE
      WRITELN('        x         y          Re[LnG(z)]          Im[LnG(z)]');
    WRITELN('    ---------------------------------------------------------');

    CGAMA(X,Y,KF,GR,GI);

//    WRITELN(' ',X:10:2, Y:10:2, GR:19:10, GI:19:10);
    WRITELN(' ',X:10:2, Y:10:2, GR:25:15, GI:25:15);
    WRITELN;

    Write('Press ENTER to close...');
    ReadLn;
  END.

{end of file mcgama.pas}
