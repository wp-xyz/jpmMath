{*************************************************************************
*   Calculate the associated Legendre Functions of the Second Kind and   *
*   their First Derivatives for a Complex Argument                       *
* ---------------------------------------------------------------------- *
* EXPLANATION:                                                           *
*                                                                        *
*      ============================================================      *
*      Purpose: This program computes the associated Legendre            *
*               functions Qmn(z) and their derivatives Qmn'(z) for       *
*               a complex argument using subroutine CLQMN                *
*      Definition: Qmn(z)=(-1)^m*(1-z*z)^(m/2)*dm/dzm[Qn(z)]             *
*                  Q0(z)=1/2*LOG[(1+z)/(1-z)]     ( for |z|<1 )          *
*                  Qmn(z)=(z*z-1)^(m/2)*dm/dzm[Qn(z)]                    *
*                  Q0(z)=1/2*LOG[(z+1)/(z-1)]     ( for |z|>1 )          *
*      Input :  x --- Real part of z                                     *
*               y --- Imaginary part of z                                *
*               m --- Order of Qmn(z)  ( m = 0,1,2,... )                 *
*               n --- Degree of Qmn(z) ( n = 0,1,2,... )                 *
*      Output:  CQM(m,n) --- Qmn(z)                                      *
*               CQD(m,n) --- Qmn'(z)                                     *
*      Examples:                                                         *
*               n = 5, x = 0.5, y = 0.2                                  *
*                                                                        *
*      m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]       *
*     -------------------------------------------------------------      *
*      0    .987156D+00   .354345D+00    .324023D+01  -.447297D+01       *
*      1   -.240328D+01   .436861D+01    .281158D+02   .171437D+02       *
*      2   -.245853D+02  -.138072D+02   -.106283D+03   .913792D+02       *
*      3    .102723D+03  -.651233D+02   -.362578D+03  -.429802D+03       *
*      4    .155510D+03   .357712D+03    .196975D+04  -.287414D+02       *
*      5   -.167357D+04  -.680954D+03   -.193093D+04  -.925757D+03       *
*                                                                        *
*               n = 5, x = 2.5, y = 1.0                                  *
*                                                                        *
*      m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]       *
*     -------------------------------------------------------------      *
*      0   -.274023D-04  -.227141D-04    .809834D-04   .210884D-04       *
*      1    .165620D-03   .136108D-03   -.489095D-03  -.124400D-03       *
*      2   -.118481D-02  -.948832D-03    .349090D-02   .825057D-03       *
*      3    .982179D-02   .753264D-02   -.288271D-01  -.596384D-02       *
*      4   -.927915D-01  -.669521D-01    .270840D+00   .451376D-01       *
*      5    .985601D+00   .656737D+00   -.285567D+01  -.332533D+00       *
*      ============================================================      *
*                                                                        *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
*                                                                        *
*  Please enter m, n, x and y: 5 5 0.5 0.2                               *
*  n = 5, m = 5, x = 0.5, y = 0.2                                        *
*                                                                        *
*    m   n     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]     *
*   -----------------------------------------------------------------    *
*    0   5       0.987156      0.354345       3.240234     -4.472973     *
*    1   5      -2.403283      4.368612      28.115787     17.143746     *
*    2   5     -24.585261    -13.807231    -106.282948     91.379180     *
*    3   5     102.723445    -65.123308    -362.578436   -429.802309     *
*    4   5     155.510087    357.711653    1969.747136    -28.741415     *
*    5   5   -1673.568876   -680.953539   -1930.929669   -925.757384     *
*                                                                        *
* ---------------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,     *
*             jin.ece.uiuc.edu/routines/routines.html".                  *
*                                                                        *
*                                   Pascal Release By J-P Moreau, Paris. *
*                                            (www.jpmoreau.fr)           *
*************************************************************************} 
PROGRAM MCLQMN;
Uses WinCrt;

  Const NMAX = 40;

  Type
        COMPLEX = Array[1..2] of Double;
        pCTab = ^CTab;
         CTab = Array[0..NMAX,0..NMAX] of COMPLEX;

  Var
        CQM, CQD: pCTab;
        J,M,N: Integer;
        X,Y  : Double;


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
       ***ROUTINES CALLED  CABS
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

      Procedure CSQRT(ZA:COMPLEX; var ZB:COMPLEX);
      {***BEGIN PROLOGUE  CSQRT
      !
      !   DOUBLE PRECISION COMPLEX SQUARE ROOT, CSQRT(ZA,ZB)
      !
      !***ROUTINES CALLED  CABS
      !***END PROLOGUE  CSQRT }
      Label 10,20,30,40,50,60,70,return;
      Var AR,AI,BR,BI,ZM, DTHETA, DPI, DRT: double;
      Begin
        DRT:=7.071067811865475244008443621E-0001;
        DPI:=3.141592653589793238462643383E+0000;
        AR:=ZA[1]; AI:=ZA[2];
        ZM := CABS(ZA);
        ZM := SQRT(ZM);
        IF AR = 0.0 THEN GOTO 10;
        IF AI = 0.0 THEN GOTO 20;
        DTHETA := ARCTAN(AI/AR);
        IF DTHETA <= 0.0 THEN GOTO 40;
        IF AR < 0.0 THEN DTHETA := DTHETA - DPI;
        GOTO 50;
10:     IF AI > 0.0 THEN GOTO 60;
        IF AI < 0.0 THEN GOTO 70;
        BR := 0.0;
        BI := 0.0;
        GOTO RETURN;
20:     IF AR > 0.0 THEN GOTO 30;
        BR := 0.0;
        BI := SQRT(ABS(AR));
        GOTO RETURN;
30:     BR := SQRT(AR);
        BI := 0.0;
        GOTO RETURN;
40:     IF AR < 0.0 THEN DTHETA := DTHETA + DPI;
50:     DTHETA := DTHETA*0.5;
        BR := ZM*COS(DTHETA);
        BI := ZM*SIN(DTHETA);
        GOTO RETURN;
60:     BR := ZM*DRT;
        BI := ZM*DRT;
        GOTO RETURN;
70:     BR := ZM*DRT;
        BI := -ZM*DRT;
return: ZB[1]:=BR; ZB[2]:=BI
      End; {CSQRT}


      Procedure CLQMN(M,N:Integer; X,Y: Double; Var CQM,CQD:pCTab);
{     ========================================================
        Purpose: Compute the associated Legendre functions of
                 the second kind, Qmn(z) and Qmn'(z), for a
                 complex argument
        Input :  x  --- Real part of z
                 y  --- Imaginary part of z
                 m  --- Order of Qmn(z)  ( m = 0,1,2,... )
                 n  --- Degree of Qmn(z) ( n = 0,1,2,... )
        Output:  COMPLEX CQM(m,n) --- Qmn(z)
                 COMPLEX CQD(m,n) --- Qmn'(z)
      =======================================================}
      Label Return;
      Var CQ0,CQ1,CQ10,CQF,CQF0,CQF1,CQF2,TMP,TMP1,Z,ZQ,ZS: COMPLEX;
          Err,I,K,KM,LS: Integer;
          XC: Double;
      Begin
        Z[1]:=X; Z[2]:=Y;
        IF (ABS(X) = 1.0) AND (Y = 0.0) THEN
        begin
           For I:=0 to M do
             For J:=0 to N do
             begin
               CQM^[I,J][1]:=1E+200; CQM^[I,J][2]:=0.0;
               CQD^[I,J][1]:=1E+200; CQD^[I,J][2]:=0.0
             end;
           goto Return
        end;
        XC:=CABS(Z);
        IF (Z[2] = 0.0) OR (XC < 1.0) Then LS:=1;
        IF XC > 1.0 Then LS:=-1;
        {ZQ=CSQRT(LS*(1.0-Z*Z)) }
        CMUL(Z,Z,TMP); TMP[1]:=LS*(1.0-TMP[1]); TMP[2]:=-LS*TMP[2];
        CSQRT(TMP,ZQ);
        {ZS=LS*(1.0-Z*Z) }
        ZS[1]:=TMP[1]; ZS[2]:=TMP[2];
        {CQ0=0.5*CLOG(LS*(1.0+Z)/(1.0-Z)) }
        TMP[1]:=LS*(1.0+Z[1]); TMP[2]:=LS*Z[2];
        TMP1[1]:=1.0-Z[1]; TMP1[2]:=-Z[2];
        CDIV(TMP,TMP1,TMP);
        CLOG(TMP,CQ0,Err); CQ0[1]:=0.5*CQ0[1]; CQ0[2]:=0.5*CQ0[2];
        IF XC < 1.0001 THEN
        begin
           CQM^[0,0][1]:=CQ0[1]; CQM^[0,0][2]:=CQ0[2];
           {CQM(0,1)=Z*CQ0-1.0 }
           CMUL(Z,CQ0,TMP);
           CQM^[0,1][1]:=TMP[1]-1.0; CQM^[0,1][2]:=TMP[2];
           {CQM(1,0)=-1.0/ZQ }
           TMP[1]:=-1.0; TMP[2]:=0.0; CDIV(TMP,ZQ,CQM^[1,0]);
           {CQM(1,1)=-ZQ*(CQ0+Z/(1.0-Z*Z)) }
           CMUL(Z,Z,TMP); TMP[1]:=1.0-TMP[1]; TMP[2]:=-TMP[2];
           CDIV(Z,TMP,TMP);
           TMP[1]:=CQ0[1]+TMP[1]; TMP[2]:=CQ0[2]+TMP[2];
           CMUL(ZQ,TMP,CQM^[1,1]);
           CQM^[1,1][1]:=-CQM^[1,1][1]; CQM^[1,1][2]:=-CQM^[1,1][2];
           For I:=0 to 1 do
             For J:=2 to N do
             begin
               {CQM(I,J)=((2.0*J-1.0)*Z*CQM(I,J-1)-(J+I-1.0)*CQM(I,J-2))/(J-I) }
               CMUL(Z,CQM^[I,J-1],TMP); TMP[1]:=(2.0*J-1.0)*TMP[1]; TMP[2]:=(2.0*J-1.0)*TMP[2];
               TMP1[1]:=(J+I-1.0)*CQM^[I,J-2][1]; TMP1[2]:=(J+I-1.0)*CQM^[I,J-2][2];
               CQM^[I,J][1]:=(TMP[1]-TMP1[1])/(J-I);
               CQM^[I,J][2]:=(TMP[2]-TMP1[2])/(J-I)
             end;
           For J:=0 to N do
             For I:=2 to M do
             begin
               {CQM(I,J)=-2.0*(I-1.0)*Z/ZQ*CQM(I-1,J)-LS*(J+I-1.0)*(J-I+2.0)*CQM(I-2,J) }
               CDIV(Z,ZQ,TMP); TMP[1]:=-2.0*(I-1.0)*TMP[1]; TMP[2]:=-2.0*(I-1.0)*TMP[2];
               CMUL(TMP,CQM^[I-1,J],TMP);                      {Now TMP=-2.0*(I-1.0)*Z/ZQ*CQM(I-1,J) }
               TMP1[1]:=LS*(J+I-1.0)*(J-I+2.0)*CQM^[I-2,J][1];
               TMP1[2]:=LS*(J+I-1.0)*(J-I+2.0)*CQM^[I-2,J][2]; {Now TMP1=LS*(J+I-1.0)*(J-I+2.0)*CQM(I-2,J) }
               CQM^[I,J][1]:=TMP[1]-TMP1[1]; CQM^[I,J][2]:=TMP[2]-TMP1[2]
             end       
        end
        ELSE
        begin
           IF XC > 1.1 THEN
              KM:=40+M+N
           ELSE
              KM:=(40+M+N)*Round(INT(-1.0-1.8*Ln(XC-1.0)));
           CQF2[1]:=0.0; CQF2[2]:=0.0;
           CQF1[1]:=1.0; CQF1[2]:=0.0;
           For K:=KM Downto 0 do
           begin
             {CQF0:=((2*K+3.0)*Z*CQF1-(K+2.0)*CQF2)/(K+1.0) }
             CMUL(Z,CQF1,TMP); TMP[1]:=(2*K+3.0)*TMP[1]; TMP[2]:=(2*K+3.0)*TMP[2];
             TMP1[1]:=(K+2.0)*CQF2[1]; TMP1[2]:=(K+2.0)*CQF2[2];
             CQF0[1]:=(TMP[1]-TMP1[1])/(K+1.0); CQF0[2]:=(TMP[2]-TMP1[2])/(K+1.0);
             IF K <= N then
             begin
               CQM^[0,K][1]:=CQF0[1]; CQM^[0,K][2]:=CQF0[2]
             end;
             CQF2[1]:=CQF1[1]; CQF2[2]:=CQF1[2];
             CQF1[1]:=CQF0[1]; CQF1[2]:=CQF0[2]
           end;
           For K:=0 to N do
           begin
             {CQM(0,K)=CQ0*CQM(0,K)/CQF0 }
             CMUL(CQ0,CQM^[0,K],TMP);
             CDIV(TMP,CQF0,CQM^[0,K])
           end;
           CQF2[1]:=0.0; CQF2[2]:=0.0;
           CQF1[1]:=1.0; CQF1[2]:=0.0;
           For K:=KM Downto 0 do
           begin
             {CQF0=((2*K+3.0)*Z*CQF1-(K+1.0)*CQF2)/(K+2.0) }
             CMUL(Z,CQF1,TMP); TMP[1]:=(2*K+3.0)*TMP[1]; TMP[2]:=(2*K+3.0)*TMP[2];
             TMP1[1]:=(K+1.0)*CQF2[1]; TMP1[2]:=(K+1.0)*CQF2[2];
             CQF0[1]:=(TMP[1]-TMP1[1])/(K+2.0); CQF0[2]:=(TMP[2]-TMP1[2])/(K+2.0);
             IF K <= N then
             begin
               CQM^[1,K][1]:=CQF0[1]; CQM^[1,K][2]:=CQF0[2]
             end;
             CQF2[1]:=CQF1[1]; CQF2[2]:=CQF1[2];
             CQF1[1]:=CQF0[1]; CQF1[2]:=CQF0[2]
           end; 
           {CQ10=-1.0/ZQ }
           TMP[1]:=-1.0; TMP[2]:=0.0; CDIV(TMP,ZQ,CQ10);
           For K:=0 to N do
           begin
             {CQM(1,K)=CQ10*CQM(1,K)/CQF0 }
             CMUL(CQ10,CQM^[1,K],TMP);
             CDIV(TMP,CQF0,CQM^[1,K])
           end;
           For J:=0 to N do
           begin
             CQ0[1]:=CQM^[0,J][1]; CQ0[2]:=CQM^[0,J][2];
             CQ1[1]:=CQM^[1,J][1]; CQ1[2]:=CQM^[1,J][2];
             For I:=0 to M-2 do
             begin
               {CQF=-2.0*(I+1)*Z/ZQ*CQ1+(J-I)*(J+I+1.0)*CQ0 }
               CDIV(Z,ZQ,TMP); TMP[1]:=-2.0*(I+1)*TMP[1]; TMP[2]:=-2.0*(I+1)*TMP[2];
               CMUL(TMP,CQ1,TMP);               {Now TMP=-2.0*(I+1)*Z/ZQ*CQ1}
               TMP1[1]:=(J-I)*(J+I+1.0)*CQ0[1];
               TMP1[2]:=(J-I)*(J+I+1.0)*CQ0[2]; {Now TMP1=(J-I)*(J+I+1.0)*CQ0}
               CQF[1]:=TMP[1]+TMP1[1]; CQF[2]:=TMP[2]+TMP1[2];
               CQM^[I+2,J][1]:=CQF[1]; CQM^[I+2,J][2]:=CQF[2];
               CQ0[1]:=CQ1[1]; CQ0[2]:=CQ1[2];
               CQ1[1]:=CQF[1]; CQ1[2]:=CQF[2]
             end
           end
        end;
        {CQD(0,0)=LS/ZS }
        TMP[1]:=1.0*LS; TMP[2]:=0.0; CDIV(TMP,ZS,CQD^[0,0]);
        For J:=1 to N do
        begin
          {CQD(0,J)=LS*J*(CQM(0,J-1)-Z*CQM(0,J))/ZS }
          CMUL(Z,CQM^[0,J],TMP);
          TMP[1]:=LS*J*(CQM^[0,J-1][1]-TMP[1]);
          TMP[2]:=LS*J*(CQM^[0,J-1][2]-TMP[2]);
          CDIV(TMP,ZS,CQD^[0,J])
        end;
        For J:=0 to N do
          For I:=1 to M do
          begin
            {CQD(I,J):=LS*I*Z/ZS*CQM(I,J)+(I+J)*(J-I+1.0)/ZQ*CQM(I-1,J) }
            CDIV(Z,ZS,TMP); TMP[1]:=LS*I*TMP[1]; TMP[2]:=LS*I*TMP[2];
            CMUL(TMP,CQM^[I,J],TMP);     {Now TMP=LS*I*Z/ZS*CQM(I,J) }
            TMP1[1]:=(I+J)*(J-I+1.0); TMP1[2]:=0.0;
            CDIV(TMP1,ZQ,TMP1);
            CMUL(TMP1,CQM^[I-1,J],TMP1); {Now TMP1=(I+J)*(J-I+1.0)/ZQ*CQM(I-1,J) }
            CQD^[I,J][1]:=TMP[1]+TMP1[1];
            CQD^[I,J][2]:=TMP[2]+TMP1[2]
          end;
Return: End; {CLQMN}


      {main program}
      BEGIN

        New(CQM); New(CQD);

	Writeln;
        Write('  Please enter m, n, x and y: '); Readln(M, N, X, Y);
        Writeln('  m =', M:2,', n =', N:2,', x =', X:4:1,' y = ', Y:4:1);
	Writeln;

        CLQMN(M,N,X,Y,CQM,CQD);

        Writeln('    m   n     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn''(z)]   Im[Qmn''(z)]');
        Writeln('  -----------------------------------------------------------------');
        For J:=0 to M do
          Writeln(' ',J:4, N:4, ' ',CQM^[J,N][1]:14:6, CQM^[J,N][2]:14:6,
                                ' ',CQD^[J,N][1]:14:6, CQD^[J,N][2]:14:6);

        Writeln;
        ReadKey;
        Dispose(CQM); Dispose(CQD);
        DoneWinCrt

      END.

{end of file mclqmn.pas}