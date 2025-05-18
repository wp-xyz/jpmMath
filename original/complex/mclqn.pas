{**********************************************************************
*       Calculate the Legendre polynomials for a complex argument     *
* ------------------------------------------------------------------- *
* EXPLANATION:                                                        *
*                                                                     *
*       ==========================================================    *
*       Purpose: This program computes the Legendre polynomials       *
*                Qn(z) and Qn'(z) for a complex argument using        *
*                subroutine CLQN                                      *
*       Input :  x --- Real part of z                                 *
*                y --- Imaginary part of z                            *
*                n --- Degree of Qn(z), n = 0,1,...                   *
*       Output:  CQN(n) --- Qn(z)                                     *
*                CQD(n) --- Qn'(z)                                    *
*       Examples:                                                     *
*                                                                     *
*       z = 0.5 + 0.5 i                                               *
*       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]     *
*      -----------------------------------------------------------    *
*       0   .402359D+00   .553574D+00   .800000D+00   .400000D+00     *
*       1  -.107561D+01   .477967D+00   .602359D+00   .115357D+01     *
*       2  -.136636D+01  -.725018D+00  -.242682D+01   .183390D+01     *
*       3   .182619D+00  -.206146D+01  -.622944D+01  -.247151D+01     *
*       4   .298834D+01  -.110022D+01  -.114849D+01  -.125963D+02     *
*       5   .353361D+01   .334847D+01   .206656D+02  -.123735D+02     *
*                                                                     *
*       z = 3.0 + 2.0 i                                               *
*       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]     *
*      -----------------------------------------------------------    *
*       0   .229073D+00  -.160875D+00  -.250000D-01   .750000D-01     *
*       1   .896860D-02  -.244805D-01   .407268D-02   .141247D-01     *
*       2  -.736230D-03  -.281865D-02   .190581D-02   .155860D-02     *
*       3  -.264727D-03  -.227023D-03   .391535D-03   .314880D-04     *
*       4  -.430648D-04  -.443187D-05   .527190D-04  -.305592D-04     *
*       5  -.481362D-05   .265297D-05   .395108D-05  -.839883D-05     *
*       ==========================================================    *
*                                                                     *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
*                                                                     *
*  Please enter Nmax, x and y (z=x+iy): 5 3 2                         *
*    x =  3.0,  y =  2.0                                              *
*                                                                     *
*    n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]        *
*   -----------------------------------------------------------       *
*    0    0.22907268   -0.16087528   -0.02500000    0.07500000        *
*    1    0.00896860   -0.02448047    0.00407268    0.01412472        *
*    2   -0.00073623   -0.00281865    0.00190581    0.00155860        *
*    3   -0.00026473   -0.00022702    0.00039153    0.00003149        *
*    4   -0.00004306   -0.00000443    0.00005272   -0.00003056        *
*    5   -0.00000481    0.00000265    0.00000395   -0.00000840        *
*                                                                     *
* ------------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,  *
*             jin.ece.uiuc.edu/routines/routines.html".               *
*                                                                     *
*                                Pascal Release By J-P Moreau, Paris. *
*                                         (www.jpmoreau.fr)           *
**********************************************************************}
      PROGRAM MCLQN;
      Uses WinCrt;

      Const NMAX = 100;

      Type  COMPLEX = Array[1..2] of Double;
            pCVEC = ^CVEC;
             CVEC = Array[0..NMAX] of COMPLEX;

      Var   CQN, CQD: pCVEC;
            K,N: Integer;
            X,Y: Double;


      {ABSOLUTE VALUE OF A COMPLEX NUMBER Z=X+I*Y }
      FUNCTION CABS(Z:COMPLEX): Double;
      Label 10, 20;
      Var
          U,V,Q,S: Double;
      Begin
        U := ABS(Z[1]);
        V := ABS(Z[2]);
        S := U+V;
{--------------------------------------------------------------------
      S*1.0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
      TRUE FLOATING ZERO
--------------------------------------------------------------------}
        S := S*1.0;
        if S = 0.0 then goto 20;
        if U > V  then goto 10;
        Q := U/V;
        CABS:=V*sqrt(1.0+Q*Q);
        exit;
10:     Q := V/U;
        CABS:=U*sqrt(1.0+Q*Q);
        exit;
20:     CABS:=0.0
      End;

      {Z=Z1/Z2}
      Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
      Var D:double;
      Begin
        D:=Z2[1]*Z2[1]+Z2[2]*Z2[2];
        if D<1e-12 then exit;
        Z[1]:=(Z1[1]*Z2[1]+Z1[2]*Z2[2])/D;
        Z[2]:=(Z1[2]*Z2[1]-Z1[1]*Z2[2])/D
      End;

      Procedure CLOG(ZA:COMPLEX; var ZB:COMPLEX; var IERR:integer);
      {***BEGIN PROLOGUE  CLOG
            DOUBLE PRECISION COMPLEX LOGARITHM B:=CLOG(A)
            IERR:=0,NORMAL RETURN      IERR:=1, Z:=CMPLX(0.0,0.0)
       ***ROUTINES CALLED  ZABS
       ***END PROLOGUE  CLOG    }
      Label 10,20,30,40,50,60,return;
      Var AR,AI,BR,BI,ZM, DTHETA, DPI, DHPI: double;
      Begin

        DPI := 3.141592653589793238462643383E+0000;
        DHPI:= 1.570796326794896619231321696E+0000;

        IERR:=0;
        AR:=ZA[1]; AI:=ZA[2];
        IF AR = 0.0 THEN GOTO 10;
        IF AI = 0.0 THEN GOTO 20;
        DTHETA := ArcTan(AI/AR);
        IF DTHETA <= 0.0 THEN GOTO 40;
        IF AR < 0.0 THEN DTHETA := DTHETA - DPI;
        GOTO 50;
10:     IF AI = 0.0 THEN GOTO 60;
        BI := DHPI;
        BR := Ln(ABS(AI));
        IF AI < 0.0 THEN BI := -BI;
        GOTO RETURN;
20:     IF AR > 0.0 THEN GOTO 30;
        BR := Ln(ABS(AR));
        BI := DPI;
        GOTO RETURN;
30:     BR := Ln(AR);
        BI := 0.0;
        GOTO RETURN;
40:     IF AR < 0.0 THEN DTHETA := DTHETA + DPI;
50:     ZM := CABS(ZA);
        BR := Ln(ZM);
        BI := DTHETA;
        GOTO RETURN;
60:     IERR:=1;
return: ZB[1]:=BR; ZB[2]:=BI;
      End; {CLOG}

      {Z=Z1*Z2}
      Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
      Begin
        Z[1]:=Z1[1]*Z2[1] - Z1[2]*Z2[2];
        Z[2]:=Z1[1]*Z2[2] + Z1[2]*Z2[1]
      End;

      Procedure CLQN(N:Integer; X,Y: Double; CQN,CQD:pCVEC);
{     =====================================================
        Purpose: Compute the Legendre functions Qn(z) and
                 their derivatives Qn'(z) for a complex
                 argument
        Input :  x --- Real part of z
                 y --- Imaginary part of z
                 n --- Degree of Qn(z), n = 0,1,2,...
        Output:  CQN(n) --- Qn(z)
                 CQD(n) --- Qn'(z)
      ===================================================== }
      Label Return;
      Var 
          CONE,CQ0,CQ1,CQF0,CQF1,CQF2,Z,TMP,TMP1,TMP2: COMPLEX;
          Err,K,KM,LS:Integer;
      Begin

        CONE[1]:=1.0; CONE[2]:=0.0;

        Z[1]:=X; Z[2]:=Y;
        IF (Z[1]=1.0) AND (Z[2]=0.0) THEN
        begin
           For K:=0 to N do
           begin
             CQN^[K][1]:=1.0E+200; CQN^[K][2]:=0.0;
             CQD^[K][1]:=1.0E+200; CQD^[K][2]:=0.0
           end;
           goto Return
        end;
        LS:=1;
        IF CABS(Z) > 1.0 Then LS:=-1;
        {CQ0:=0.5*CLOG(LS*(1.0+Z)/(1.0-Z)) }
        TMP[1]:=LS*(CONE[1]+Z[1]);
        TMP[2]:=LS*(CONE[2]+Z[2]);
        TMP1[1]:=CONE[1]-Z[1];
        TMP1[2]:=CONE[2]-Z[2];
        CDIV(TMP,TMP1,TMP2);
        CLOG(TMP2,CQ0,Err);
        CQ0[1]:=CQ0[1]*0.5; CQ0[2]:=CQ0[2]*0.5;
        {CQ1=Z*CQ0-1.0 }
        CMUL(Z,CQ0,CQ1);
        CQ1[1]:=CQ1[1]-CONE[1];
        CQN^[0][1]:=CQ0[1]; CQN^[0][2]:=CQ0[2];
        CQN^[1][1]:=CQ1[1]; CQN^[1][2]:=CQ1[2];
        IF CABS(Z) < 1.0001 THEN
        begin
           CQF0[1]:=CQ0[1]; CQF0[2]:=CQ0[2];
           CQF1[1]:=CQ1[1]; CQF1[2]:=CQ1[2];
           For K:=2 to N do
           begin
              {CQF2=((2.0*K-1.0)*Z*CQF1-(K-1.0)*CQF0)/K }
              CMUL(Z,CQF1,TMP);
              TMP[1]:=TMP[1]*(2.0*K-1.0);
              TMP[2]:=TMP[2]*(2.0*K-1.0);
              TMP1[1]:=(K-1.0)*CQF0[1];
              TMP1[2]:=(K-1.0)*CQF0[2];
              CQF2[1]:=(TMP[1]-TMP1[1])/K;
              CQF2[2]:=(TMP[2]-TMP1[2])/K;
              CQN^[K][1]:=CQF2[1]; CQN^[K][2]:=CQF2[2];
              CQF0[1]:=CQF1[1]; CQF0[2]:=CQF1[2];
              CQF1[1]:=CQF2[1]; CQF1[2]:=CQF2[2]
           end
        end
        ELSE
        begin
           IF CABS(Z) > 1.1 THEN
              KM:=40+N
           ELSE
           begin
              TMP[1]:=Z[1]-1.0; TMP[2]:=Z[2];
              KM:=(40+N)*Round(INT(-1.0-1.8*Ln(CABS(TMP))))
           end;
           CQF2[1]:=0.0; CQF2[2]:=0.0;
           CQF1[1]:=1.0; CQF1[2]:=0.0;
           For K:=KM Downto 0 do
           begin
              {CQF0=((2*K+3.0)*Z*CQF1-(K+2.0)*CQF2)/(K+1.0) }
              CMUL(Z,CQF1,TMP);
              TMP[1]:=TMP[1]*(2.0*K+3.0);
              TMP[2]:=TMP[2]*(2.0*K+3.0);
              TMP1[1]:=(K+2.0)*CQF2[1];
              TMP1[2]:=(K+2.0)*CQF2[2];
              CQF0[1]:=(TMP[1]-TMP1[1])/(1.0+K);
              CQF0[2]:=(TMP[2]-TMP1[2])/(1.0+K);
              IF K <= N then
              begin
                CQN^[K][1]:=CQF0[1]; CQN^[K][2]:=CQF0[2];
              end;
              CQF2[1]:=CQF1[1]; CQF2[2]:=CQF1[2];
              CQF1[1]:=CQF0[1]; CQF1[2]:=CQF0[2]
           end;
           For K:=0 to N do
           begin
             {CQN(K)=CQN(K)*CQ0/CQF0 }
             CMUL(CQN^[K],CQ0,TMP); CDIV(TMP,CQF0,CQN^[K])
           end
        end;
        {CQD(0)=(CQN(1)-Z*CQN(0))/(Z*Z-1.0) }
        CMUL(Z,CQN^[0],TMP); TMP[1]:=CQN^[1][1]-TMP[1]; TMP[2]:=CQN^[1][2]-TMP[2];
        CMUL(Z,Z,TMP1); TMP1[1]:=TMP1[1]-1.0; CDIV(TMP,TMP1,CQD^[0]);
        For K:=1 to N do
        begin
          {CQD(K)=(K*Z*CQN(K)-K*CQN(K-1))/(Z*Z-1.0) }
          CMUL(Z,CQN^[K],TMP); TMP[1]:=K*TMP[1]; TMP[2]:=K*TMP[2];
          TMP2[1]:=K*CQN^[K-1][1]; TMP2[2]:=K*CQN^[K-1][2];
          TMP[1]:=TMP[1]-TMP2[1]; TMP[2]:=TMP[2]-TMP2[2];
          CDIV(TMP,TMP1,CQD^[K])
        end;
Return:End;


      {main program}
      BEGIN

        New(CQN); New(CQD);
	Writeln;
        Write(' Please enter Nmax, x and y (z=x+iy): ');
        Readln(N,X,Y);
        Writeln('   x =',X:5:1,',  y =',Y:5:1);
        Writeln;

        CLQN(N,X,Y,CQN,CQD);

        Writeln('   n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn''(z)]    Im[Qn''(z)]');
        Writeln('  -----------------------------------------------------------');
        For K:=0 to N do
          Writeln(K:4,CQN^[K][1]:14:8,CQN^[K][2]:14:8,CQD^[K][1]:14:8,CQD^[K][2]:14:8);

        Writeln;
        ReadKey;

        Dispose(CQN); Dispose(CQD);
        DoneWinCrt

      END.

{end of file mclqn.pas}