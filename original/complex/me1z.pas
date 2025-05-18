{********************************************************************
*  Calculate the complex exponential integral E1(z) with a complex  *
*  argument                                                         *
* ----------------------------------------------------------------- *
* EXPLANATION:                                                      *
*                                                                   *
*      =========================================================    *
*      Purpose: This program computes the complex exponential       *
*               integral E1(z) using subroutine E1Z                 *
*      Example:                                                     *
*                    z            Re[E1(z)]       Im[E1(z)]         *
*               -----------------------------------------------     *
*                3.0    2.0    -.90959209D-02   -.69001793D-02      *
*                3.0   -2.0    -.90959209D-02    .69001793D-02      *
*               -3.0    2.0    -.28074891D+01    .59603353D+01      *
*               -3.0   -2.0    -.28074891D+01   -.59603353D+01      *
*               25.0   10.0    -.29302080D-12    .40391222D-12      *
*               25.0  -10.0    -.29302080D-12   -.40391222D-12      *
*              -25.0   10.0     .27279957D+10   -.49430610D+09      *
*              -25.0  -10.0     .27279957D+10    .49430610D+09      *
*      =========================================================    *
*                                                                   *
* ----------------------------------------------------------------- *
* SAMPLE RUNS:                                                      *
*                                                                   *
*  Please enter x and y (z =x+iy): 3 2                              *
*                                                                   *
*         z           Re[E1(z)]             Im[E1(z)]               *
*   ----------------------------------------------------------      *
*    3.0    2.0    -9.09592087479E-0003  -6.90017926221E-0003       *
*                                                                   *
*  Please enter x and y (z =x+iy): 25 10                            *
*                                                                   *
*         z           Re[E1(z)]             Im[E1(z)]               *
*   ----------------------------------------------------------      *
*   25.0   10.0    -5.53483313811E-0014   7.83427777745E-0014       *
*                                                                   *
* ----------------------------------------------------------------- *
* REFERENCE:"Fortran Routines for Computation of Special Functions, *
*            jin.ece.uiuc.edu/routines/routines.html".              *
*                                                                   *
*                              Pascal Release By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)           *
********************************************************************}
      PROGRAM ME1Z;
      Uses WinCrt;

      Type COMPLEX = Array[1..2] of Double;
        
      Var  CE1, Z: COMPLEX;
           x,y: Double;


      {ABSOLUTE VALUE OF A COMPLEX NUMBER Z=X+I*Y }
      FUNCTION CABS(Z:COMPLEX): Double;
      Var
          X,Y,W: Double;
      Begin
        X:=ABS(Z[1]);
        Y:=ABS(Z[2]);
        IF X = 0.0 THEN
          W:=Y
        ELSE
          IF Y = 0.0 THEN
            W:=X
          ELSE
            IF X > Y THEN
              W:=X*SQRT(1.0+(Y/X)*(Y/X))
            ELSE
              W:=Y*SQRT(1.0+(X/Y)*(Y/X));
        CABS:=W
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

      {z1=exp(z) }
      Procedure CEXP(Z:Complex; Var Z1:Complex);
      Var tmp:double;
      Begin
        if exp(Z[1]) > 1e16 then tmp:=1e16 else tmp:=exp(Z[1]);
        Z1[1]:=tmp*cos(Z[2]);
        Z1[2]:=tmp*sin(Z[2])
      End;

      Procedure CLOG(ZA:COMPLEX; var ZB:COMPLEX; var IERR:integer);
      {***BEGIN PROLOGUE  ZLOG
            DOUBLE PRECISION COMPLEX LOGARITHM B:=CLOG(A)
            IERR:=0,NORMAL RETURN      IERR:=1, Z:=CMPLX(0.0,0.0)
       ***ROUTINES CALLED  ZABS
       ***END PROLOGUE  ZLOG    }
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

      Procedure E1Z(Z:COMPLEX; Var CE1:COMPLEX);
{       ====================================================
        Purpose: Compute complex exponential integral E1(z)
        Input :  z   --- Argument of E1(z)
        Output:  CE1 --- E1(z)
        ==================================================== }
      Label 15;
      Var A0,EL,X: Double;
          CK,CONE,CR,CT,CT0,TMP,TMP1: COMPLEX;
          IERR,K:Integer;
      Begin
        EL:=0.5772156649015328;
        X:=Z[1];
        A0:=CABS(Z);
        IF A0 = 0.0 THEN
        begin
           CE1[1]:=1.0E+200; CE1[2]:=0.0
        end
        ELSE IF ((A0 <= 10.0) OR (X < 0.0) AND (A0 <20.0)) THEN
        begin
           CE1[1]:=1.0; CE1[2]:=0.0;
           CR[1]:=1.0; CR[2]:=0.0;
           For K:=1 to 150 do
           begin
             {CR=-CR*K*Z/Sqr(K+1.0) }
             CMUL(CR,Z,CR);
             CR[1]:=-K*CR[1]/Sqr(K+1.0);
             CR[2]:=-K*CR[2]/Sqr(K+1.0);
             CE1[1]:=CE1[1]+CR[1];
             CE1[2]:=CE1[2]+CR[2];
             IF CABS(CR) <= CABS(CE1)*1.0E-15 Then GOTO 15
           end;
15:        {CE1=-EL-CLOG(Z)+Z*CE1}
           CLOG(Z,TMP,IERR); CMUL(Z,CE1,TMP1);
           CE1[1]:=-EL-TMP[1]+TMP1[1];
           CE1[2]:=   -TMP[2]+TMP1[2];
        end
        ELSE
        begin
           CT0[1]:=0.0; CT0[2]:=0.0;
           For K:=120 Downto 1 do
           begin
             {CT0=K/(1.0+K/(Z+CT0)) }
             CK[1]:=1.0*K; CK[2]:=0.0;
             TMP[1]:=1.0+CK[1]; TMP[2]:=0.0;
             TMP1[1]:=Z[1]+CT0[1];
             TMP1[2]:=Z[2]+CT0[2];
             CDIV(TMP,TMP1,TMP);
             CDIV(CK,TMP,CT0)
           end;
           {CT=1.0/(Z+CT0) }
           CONE[1]:=1.0; CONE[2]:=0.0;
           TMP[1]:=Z[1]+CT0[1]; TMP[2]:=Z[2]+CT0[2];
           CDIV(CONE,TMP,CT);
           TMP[1]:=-Z[1]; TMP[2]:=-Z[2];
           CEXP(TMP,TMP1);
           CMUL(TMP1,CT,CE1);
           IF (X <= 0.0) AND (Z[2] = 0.0) Then
           begin
             {CE1=CE1-PI*(0.0D0,1.0D0) }
             TMP[1]:=0.0; TMP[2]:=1.0;
             CE1[1]:=CE1[1]-PI*TMP[1];
             CE1[2]:=CE1[2]-PI*TMP[2]
           end
        end
      End;


    {main program}
    BEGIN

      WRITELN;
      WRITE(' Please enter x and y (z =x+iy): '); READLN(X, Y);

      Z[1]:=X; Z[2]:=Y;

      E1Z(Z, CE1);

      WRITELN;
      WRITELN('       z           Re[E1(z)]             Im[E1(z)]');
      WRITELN(' ----------------------------------------------------------');
      WRITELN(' ',X:5:1,'  ',Y:5:1,'   ', CE1[1]:20, '  ',CE1[2]:20);
      WRITELN;

      ReadKey;
      DoneWinCrt

    END.

{ end of file me1z.pas}