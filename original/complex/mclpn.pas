{**********************************************************************
*      Calculate the Legendre Polynomials for a Complex Argument      *
* ------------------------------------------------------------------- *
* EXPLANATION:                                                        *
*                                                                     *
*      ==========================================================     *
*      Purpose: This program computes the Legendre polynomials        *
*               Pn(z) and Pn'(z) for a complex argument using         *
*               subroutine CLPN                                       *
*      Input :  x --- Real part of z                                  *
*               y --- Imaginary part of z                             *
*               n --- Degree of Pn(z), n = 0,1,...,N                  *
*      Output:  Complex CPN(n) --- Pn(z)                              *
*               Complex CPD(n) --- Pn'(z)                             *
*                                                                     *
*      Example: z = 3.0 + 2.0 i                                       *
*                                                                     *
*      n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn'(z)]   Im[Pn'(z)]       *
*     -----------------------------------------------------------     *
*      0   .100000D+01   .000000D+00   .000000D+00   .000000D+00      *
*      1   .300000D+01   .200000D+01   .100000D+01   .000000D+00      *
*      2   .700000D+01   .180000D+02   .900000D+01   .600000D+01      *
*      3  -.270000D+02   .112000D+03   .360000D+02   .900000D+02      *
*      4  -.539000D+03   .480000D+03  -.180000D+03   .790000D+03      *
*      5  -.461700D+04   .562000D+03  -.481500D+04   .441000D+04      *
*      ==========================================================     *
*                                                                     *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
*                                                                     *
*   Please enter Nmax, x and y (z=x+iy): 5 3 2                        *
*    x =  3.0,  y =  2.0                                              *
*                                                                     *
*    n    Re[Pn(z)]     Im[Pn(z)]      Re[Pn'(z)]   Im[Pn'(z)]        *
*   -----------------------------------------------------------       *
*    0      1.000000      0.000000       0.000000     0.000000        *
*    1      3.000000      2.000000       1.000000     0.000000        *
*    2      7.000000     18.000000       9.000000     6.000000        *
*    3    -27.000000    112.000000      36.000000    90.000000        *
*    4   -539.000000    480.000000    -180.000000   790.000000        *
*    5  -4617.00D000    562.000000   -4815.000000  4410.000000        *
*                                                                     *
* ------------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,  *
*             jin.ece.uiuc.edu/routines/routines.html".               *
*                                                                     *
*                                Pascal Release By J-P Moreau, Paris. *
*                                         (www.jpmoreau.fr)           *
**********************************************************************}
PROGRAM MCLPN;
Uses WinCrt;

      Const NMAX = 100;

      Type  COMPLEX = Array[1..2] of Double;
            pCVEC = ^CVEC;
             CVEC = Array[0..NMAX] of COMPLEX;

      Var   CPN, CPD: pCVEC;
            K,N: Integer;
            X,Y: Double;


      { x^n }
      Function IPower(x:double; n:integer): double;
      {Calculate x^n}
      var i,m : integer;
      result :double;
      begin
        result := 1.0;
        if n=0 then
        begin
          IPower:=result;
          exit
        end;
        m:=  n;
        if n<0 then m:=-n;
        for i:=1 to m do result :=x*result;
        IPower :=result;
        if n<0 then IPower:=1.0/result
      end;

      {Z=Z1/Z2}
      Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
      Var D:double;
      Begin
        D:=Z2[1]*Z2[1]+Z2[2]*Z2[2];
        if D<1e-12 then exit;
        Z[1]:=(Z1[1]*Z2[1]+Z1[2]*Z2[2])/D;
        Z[2]:=(Z1[2]*Z2[1]-Z1[1]*Z2[2])/D
      End;

      {Z=Z1*Z2}
      Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
      Begin
        Z[1]:=Z1[1]*Z2[1] - Z1[2]*Z2[2];
        Z[2]:=Z1[1]*Z2[2] + Z1[2]*Z2[1]
      End;


      Procedure CLPN(N:Integer; X,Y: Double; Var CPN,CPD:pCVEC);
{     =========================================================
        Purpose: Compute Legendre polynomials Pn(z) and
                 their derivatives Pn'(z) for a complex
                 argument
        Input :  x --- Real part of z
                 y --- Imaginary part of z
                 n --- Degree of Pn(z), n = 0,1,2,...
        Output:  Complex CPN(n) --- Pn(z)
                 Complex CPD(n) --- Pn'(z)
      ========================================================= }
      Var   CP0,CP1,CPF,TMP,TMP1,Z: COMPLEX;
      Begin

        Z[1]:=X; Z[2]:=Y;
        CPN^[0][1]:=1.0; CPN^[0][2]:=0.0;
        CPN^[1][1]:=Z[1]; CPN^[1][2]:=Z[2];
        CPD^[0][1]:=0.0; CPD^[0][2]:=0.0;
        CPD^[1][1]:=1.0; CPD^[1][2]:=0.0;
        CP0[1]:=1.0; CP0[2]:=0.0;
        CP1[1]:=Z[1]; CP1[2]:=Z[2];
        For K:=2 to N do
        begin
           {CPF=(2.0*K-1.0)/K*Z*CP1-(K-1.0)/K*CP0 }
           CMUL(Z,CP1,TMP);
           TMP[1]:=(2.0*K-1.0)/K*TMP[1]; TMP[2]:=(2.0*K-1.0)/K*TMP[2];
           CPF[1]:=TMP[1]-(K-1.0)/K*CP0[1];
           CPF[2]:=TMP[2]-(K-1.0)/K*CP0[2];
           CPN^[K][1]:=CPF[1]; CPN^[K][2]:=CPF[2];
           IF (ABS(X) = 1.0) AND (Y = 0.0) THEN
           begin
             {CPD(K)=0.5*X^(K+1)*K*(K+1.0) }
             CPD^[K][1]:=0.5*IPower(X,K+1)*K*(K+1.0);
             CPD^[K][2]:=0.0
           end
           ELSE
           begin
             {CPD(K)=K*(CP1-Z*CPF)/(1.0-Z*Z) }
             CMUL(Z,CPF,TMP);
             TMP[1]:=K*(CP1[1]-TMP[1]); TMP[2]:=K*(CP1[2]-TMP[2]);
             CMUL(Z,Z,TMP1); TMP1[1]:=1.0-TMP1[1]; TMP1[2]:=-TMP1[2];
             CDIV(TMP,TMP1,CPD^[K])
           end;
           CP0[1]:=CP1[1]; CP0[2]:=CP1[2];
           CP1[1]:=CPF[1]; CP1[2]:=CPF[2]
        end
      End;


      {main program}
      BEGIN

        New(CPN); New(CPD);

        Writeln; 
        Write('  Please enter Nmax, x and y (z:=x+iy): '); Readln(N, X, Y);
        Writeln('   x =', X:5:1,', y =', Y:5:1);
        Writeln;

        CLPN(N,X,Y,CPN,CPD);

        Writeln('   n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn''(z)]   Im[Pn''(z)]');
        Writeln('  -----------------------------------------------------------');

        For K:=0 to N do
          Writeln(K:4,CPN^[K][1]:14:6,CPN^[K][2]:14:6,CPD^[K][1]:14:6,CPD^[K][2]:14:6);
 
        Writeln;
        ReadKey;
        Dispose(CPN); Dispose(CPD);
        DoneWinCrt

      END.

{end of file mclpn.pas}