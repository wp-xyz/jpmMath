{*************************************************************************
* Calculate the Asociated Legendre Functions and their First Derivatives *
* for a Complex Argument                                                 *
* ---------------------------------------------------------------------- *
* EXPLANATION:                                                           *
*                                                                        *
*      ============================================================      *
*      Purpose: This program computes the associated Legendre            *
*               functions Pmn(z) and their derivatives Pmn'(z) for       *
*               a complex argument using subroutine CLPMN                *
*      Input :  x --- Real part of z                                     *
*               y --- Imaginary part of z                                *
*               m --- Order of Pmn(z),  m = 0,1,2,...,n                  *
*               n --- Degree of Pmn(z), n = 0,1,2,...,N                  *
*      Output:  CPM(m,n) --- Pmn(z)                                      *
*               CPD(m,n) --- Pmn'(z)                                     *
*      Examples:                                                         *
*               n = 5, x = 0.5, y = 0.2                                  *
*                                                                        *
*      m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]       *
*     -------------------------------------------------------------      *
*      0    .252594D+00  -.530293D+00   -.347606D+01  -.194250D+01       *
*      1    .333071D+01   .135206D+01    .117643D+02  -.144329D+02       *
*      2   -.102769D+02   .125759D+02    .765713D+02   .598500D+02       *
*      3   -.572879D+02  -.522744D+02   -.343414D+03   .147389D+03       *
*      4    .335711D+03  -.389151D+02   -.226328D+03  -.737100D+03       *
*      5   -.461125D+03   .329122D+03    .187180D+04   .160494D+02       *
*                                                                        *
*               n = 5, x = 2.5, y = 1.0                                  *
*                                                                        *
*      m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]       *
*     -------------------------------------------------------------      *
*      0   -.429395D+03   .900336D+03   -.350391D+02   .193594D+04       *
*      1   -.216303D+04   .446358D+04   -.208935D+03   .964685D+04       *
*      2   -.883477D+04   .174005D+05   -.123703D+04   .381938D+05       *
*      3   -.273211D+05   .499684D+05   -.568080D+04   .112614D+06       *
*      4   -.565523D+05   .938503D+05   -.167147D+05   .219713D+06       *
*      5   -.584268D+05   .863328D+05   -.233002D+05   .212595D+06       *
*      ============================================================      *
*                                                                        *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
*                                                                        *
*   Please enter m, n, x and y: 5 5 0.5 0.2                              *
*   m = 5, n = 5, x = 0.5, y = 0.2                                       *
*                                                                        *
*    m   n    Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]      *
*  -----------------------------------------------------------------     *
*    0   5      0.252594     -0.530293      -3.476063     -1.942500      *
*    1   5      3.330709      1.352057      11.764333    -14.432917      *
*    2   5    -10.276875     12.575850      76.571250     59.850000      *
*    3   5    -57.287858    -52.274368    -343.414490    147.389427      *
*    4   5    335.711250    -38.915100    -226.327500   -737.100000      *
*    5   5   -461.124511    329.122362    1871.802220     16.049431      *
*                                                                        *
* ---------------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,     *
*             jin.ece.uiuc.edu/routines/routines.html".                  *
*                                                                        *
*                                   Pascal Release By J-P Moreau, Paris. *
*                                            (www.jpmoreau.fr)           *
*************************************************************************} 
PROGRAM MCLPMN;
Uses WinCrt;

  Const NMAX = 40;

  Type
        COMPLEX = Array[1..2] of Double;
        pCTab = ^CTab;
         CTab = Array[0..NMAX,0..NMAX] of COMPLEX;

  Var
        CPM, CPD: pCTab;
        J,M,N: Integer;
        X,Y  : Double;


      Function Power(x:Double; n:Integer): Double;
      {Calculate x power n}
      var i,m : integer;
          result :double;
      begin
        result := 1.0;
        if n=0 then
        begin
          Power:=result;
          exit;
        end;
        m:=  n;
        if n<0 then m:=-n;
        for i:=1 to m do result :=x*result;
        Power :=result;
        if n<0 then Power:=1.0/result;
      end;

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
      !***ROUTINES CALLED  ZABS
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
      End; {ZSQRT}


      Procedure CLPMN(M,N:Integer; X,Y:Double; Var CPM,CPD:pCTab);
{     ============================================================
        Purpose: Compute the associated Legendre functions Pmn(z)   
                 and their derivatives Pmn'(z) for a complex 
                 argument
        Input :  x  --- Real part of z
                 y  --- Imaginary part of z
                 m  --- Order of Pmn(z),  m = 0,1,2,...,n
                 n  --- Degree of Pmn(z), n = 0,1,2,...,N
                 mm --- Physical dimension of CPM and CPD
        Output:  CPM(m,n) --- Pmn(z)
                 CPD(m,n) --- Pmn'(z)
      ===========================================================}
      Label Return;
      Var TMP,TMP1,Z,ZQ,ZS: COMPLEX;
          I,J,LS: Integer;
      Begin
        Z[1]:=X; Z[2]:=Y;
        For I:=0 to N do
          For J:=0 to M do
          begin
           CPM^[J,I][1]:=0.0; CPM^[J,I][2]:=0.0;
           CPD^[J,I][1]:=0.0; CPD^[J,I][2]:=0.0
          end;
        CPM^[0,0][1]:=1.0;
        IF (ABS(X) = 1.0) AND (Y = 0.0) THEN
        begin
           For I:=1 to N do
           begin
             CPM^[0,I][1]:=Power(X,I); CPM^[0,I][2]:=0.0;
             CPD^[0,I][1]:=0.5*I*(I+1)*Power(X,(I+1));
             CPD^[0,I][2]:=0.0
           end;
           For J:=1 to N do
             For I:=1 to M do
               IF I = 1 THEN
               begin
                 CPD^[I,J][1]:=1E+200;
                 CPD^[I,J][2]:=0.0
               end
               ELSE IF I = 2 THEN
               begin
                 CPD^[I,J][1]:=-0.25*(J+2)*(J+1)*J*(J-1)*Power(X,(J+1));
                 CPD^[I,J][2]:=0.0
               end;
           goto Return
        end;
        LS:=1;
        IF CABS(Z) > 1.0 Then LS:=-1;
        {ZQ=CSQRT(LS*(1.0-Z*Z)) }
        CMUL(Z,Z,TMP);
        TMP[1]:=LS*(1.0-TMP[1]); TMP[2]:=-LS*TMP[2];
        CSQRT(TMP,ZQ);
        {ZS=LS*(1.0D0-Z*Z) }
        ZS[1]:=TMP[1]; ZS[2]:=TMP[2];
        For I:=1 to M do
        begin
          {CPM(I,I)=-LS*(2.0*I-1.0)*ZQ*CPM(I-1,I-1) }
          CMUL(ZQ,CPM^[I-1,I-1],TMP);
          CPM^[I,I][1]:=-LS*(2.0*I-1.0)*TMP[1];
          CPM^[I,I][2]:=-LS*(2.0*I-1.0)*TMP[2]
        end;
        For I:=0 to M do
        begin
          {CPM(I,I+1)=(2.0*I+1.0)*Z*CPM(I,I) }
          CMUL(Z,CPM^[I,I],TMP);
          CPM^[I,I+1][1]:=(2.0*I+1.0)*TMP[1];
          CPM^[I,I+1][2]:=(2.0*I+1.0)*TMP[2]
        end;
        For I:=0 to M do
          For J:=I+2 to N do
          begin
            {CPM(I,J)=((2.0*J-1.0)*Z*CPM(I,J-1)-(I+J-1.0)*CPM(I,J-2))/(J-I) }
            CMUL(Z,CPM^[I,J-1],TMP);
            CPM^[I,J][1]:=((2.0*J-1.0)*TMP[1]-(I+J-1.0)*CPM^[I,J-2][1])/(J-I);
            CPM^[I,J][2]:=((2.0*J-1.0)*TMP[2]-(I+J-1.0)*CPM^[I,J-2][2])/(J-I)
          end;
        CPD^[0,0][1]:=0.0; CPD^[0,0][2]:=0.0;
        For J:=1 to N do
        begin
          {CPD(0,J)=LS*J*(CPM(0,J-1)-Z*CPM(0,J))/ZS }
          CMUL(Z,CPM^[0,J],TMP);
          TMP[1]:=LS*J*(CPM^[0,J-1][1]-TMP[1]);
          TMP[2]:=LS*J*(CPM^[0,J-1][2]-TMP[2]);
          CDIV(TMP,ZS,CPD^[0,J])
        end;
        For I:=1 to M do
          For J:=I to N do
          begin
            {CPD(I,J)=LS*I*Z*CPM(I,J)/ZS+(J+I)*(J-I+1.0)*CPM(I-1,J)/ZQ }
            CMUL(Z,CPM^[I,J],TMP);
            TMP[1]:=TMP[1]*LS*I; TMP[2]:=TMP[2]*LS*I;
            CDIV(TMP,ZS,TMP);   {Now TMP=LS*I*Z*CPM(I,J)/ZS}
            TMP1[1]:=(J+I)*(J-I+1.0)*CPM^[I-1,J][1];
            TMP1[2]:=(J+I)*(J-I+1.0)*CPM^[I-1,J][2];
            CDIV(TMP1,ZQ,TMP1); {Now TMP1=(J+I)*(J-I+1.0)*CPM(I-1,J)/ZQ}
            CPD^[I,J][1]:=TMP[1]+TMP1[1]; CPD^[I,J][2]:=TMP[2]+TMP1[2]
          end;
Return: End;


   {main program}
   BEGIN

     New(CPM); New(CPD);
     writeln;
     write('  Please enter m, n, x and y: '); readln(M,N,X,Y);
     writeln('  m = ',m,', n = ',n,', x = ',x:5:1,', y = ',y:5:1);
     writeln;

     CLPMN(M,N,X,Y,CPM,CPD);

     writeln('   m   n    Re[Pmn(z)]     Im[Pmn(z)]    Re[Pmn''(z)]    Im[Pmn''(z)]');
     writeln(' -------------------------------------------------------------------');
     For J:=0 to M do
       writeln('   ',J,'   ',n, CPM^[J,N][1]:14:6,' ',CPM^[J,N][2]:14:6,' ',
                                CPD^[J,N][1]:14:6,' ',CPD^[J,N][2]:14:6);
     writeln;

     ReadKey;
     Dispose(CPM); Dispose(CPD);
     DoneWinCrt

   END.

{end of file mclpmn.pas}