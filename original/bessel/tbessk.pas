{**********************************************************
* This program tests the subroutine BESSK to calculate the*
* modified Bessel function of the third kind of order N   *
* for any real positive argument X.                       *
* ------------------------------------------------------- *
* SAMPLE RUN:                                             *
*                                                         *
*  X = 1.20970000000000E+0000                             *
*                                                         *
*  N = 0                                                  *
*  Y = 3.14324491956902E-0001                             *
*                                                         *
*  N = 1                                                  *
*  Y = 4.28050751380124E-0001                             *
*                                                         *
*  N = 2                                                  *
*  Y = 1.02202185722122E+0000                             *
*                                                         *
*  N = 3                                                  *
*  Y = 3.80747327670449E+0000                             *
*                                                         *
*  N = 4                                                  *
*  Y = 1.99067367949966E+0001                             *
*                                                         *
* ------------------------------------------------------- *
* Reference:   From Numath Library By Tuan Dang Trong     *
* in Fortran 77 [BIBLI 18].                               *
*                                                         *
*                      TPW Version By J-P Moreau, Paris.  *
*                             (www.jpmoreau.fr)           *
**********************************************************}     
PROGRAM TBESSK;
Uses WinCrt;

Const BIG: Double = 1.0E30;

Var
      X, Y: Double;
      i, N: Integer;


    FUNCTION BESSK0(X:Double): Double; Forward;
    FUNCTION BESSK1(X:Double): Double; Forward;
    FUNCTION BESSI0(X:Double): Double; Forward;
    FUNCTION BESSI1(X:Double): Double; Forward;


    FUNCTION BESSK(N:Integer; X:Double): Double;
{ ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! -----------------------------------------------------------------------}
    Label Return;
    Var TOX,BK,BKM,BKP: Double;
        J: Integer;
    Begin
      IF N = 0 THEN
      begin
        BESSK := BESSK0(X);
        goto Return
      end;
      IF N = 1 THEN
      begin
        BESSK := BESSK1(X);
        goto Return
      end;
      IF X = 0.0 THEN
      begin
        BESSK := BIG;  {arbitrary big value}
        goto Return
      end;
      TOX := 2.0/X;
      BK  := BESSK1(X);
      BKM := BESSK0(X);
      For J:=1 to N-1 do
      begin
        BKP := BKM + J*TOX*BK;
        BKM := BK;
        BK  := BKP
      end;
      BESSK := BK;
Return: End;

    FUNCTION BESSK0(X:Double): Double;
{ ----------------------------------------------------------------------
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ---------------------------------------------------------------------}
    Label Return;
    Var   Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7: Double;
    Begin
      P1:=-0.57721566; P2:= 0.42278420; P3:=0.23069756; P4:=0.3488590E-1;
      P5:= 0.262698E-2; P6:=0.10750E-3; P7:=0.74E-5;
      Q1:= 1.25331414; Q2:=-0.7832358E-1; Q3:=0.2189568E-1; Q4:=-0.1062446E-1;
      Q5:= 0.587872E-2; Q6:=-0.251540E-2; Q7:=0.53208E-3;
      IF X = 0.0 THEN
      begin
        BESSK0:=BIG;  {arbitrary big value}
        goto Return
      end;
      IF X <= 2.0 THEN
      begin
        Y:=X*X/4.0;
        AX:=-Ln(X/2.0)*BESSI0(X);
        BESSK0:=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      end
      ELSE
      begin
        Y:=2.0/X;
        AX:=Exp(-X)/Sqrt(X);
        BESSK0:=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      end;
    Return: End;

    FUNCTION BESSK1(X:Double): Double;
{ -------------------------------------------------------------------------
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ------------------------------------------------------------------------}
    Label Return;
    Var Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7: Double;
    Begin
      P1:=1.0; P2:=0.15443144; P3:=-0.67278579; P4:=-0.18156897;
      P5:=-0.1919402E-1; P6:=-0.110404E-2; P7:=-0.4686E-4;
      Q1:=1.25331414; Q2:=0.23498619; Q3:=-0.3655620E-1; Q4:=0.1504268E-1;
      Q5:=-0.780353E-2; Q6:=0.325614E-2; Q7:=-0.68245E-3;
      IF X = 0.0 THEN
      begin
        BESSK1:=BIG;
        goto Return
      end;
      IF X <= 2.0 THEN
      begin
        Y:=X*X/4.0;
        AX:=Ln(X/2.0)*BESSI1(X);
        BESSK1:=AX+(1.0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      end
      ELSE
      begin
        Y:=2.0/X;
        AX:=Exp(-X)/Sqrt(X);
        BESSK1:=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      end;
    Return: End;

{   Bessel Function of the 1st kind of order zero }

    FUNCTION BESSI0(X:Double): Double;
    Var Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX: Double;
    Begin
      P1:=1.0; P2:=3.5156229; P3:=3.0899424; P4:=1.2067429;
      P5:=0.2659732; P6:=0.360768E-1; P7:=0.45813E-2;
      Q1:=0.39894228; Q2:=0.1328592E-1; Q3:=0.225319E-2; Q4:=-0.157565E-2;
      Q5:=0.916281E-2; Q6:=-0.2057706E-1; Q7:=0.2635537E-1;
      Q8:=-0.1647633E-1; Q9:=0.392377E-2;
      IF ABS(X) < 3.75 THEN
      begin
        Y:=Sqr(X/3.75);
        BESSI0:=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      end
      ELSE
      begin
        AX:=ABS(X);
        Y:=3.75/AX;
        BX:=Exp(AX)/Sqrt(AX);
        AX:=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        BESSI0:=AX*BX
      end;
    End;

{   Bessel Function of the 1st kind of order one }

    FUNCTION BESSI1(X:Double): Double;
    Var  Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX: Double;
    Begin
      P1:=0.5; P2:=0.87890594; P3:=0.51498869; P4:=0.15084934;
      P5:=0.2658733E-1; P6:=0.301532E-2; P7:=0.32411E-3;
      Q1:=0.39894228; Q2:=-0.3988024E-1; Q3:=-0.362018E-2;
      Q4:=0.163801E-2;Q5:=-0.1031555E-1; Q6:= 0.2282967E-1;
      Q7:=-0.2895312E-1; Q8:=0.1787654E-1; Q9:=-0.420059E-2;
      IF ABS(X) < 3.75 THEN
      begin
        Y:=Sqr(X/3.75);
        BESSI1:=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      end
      ELSE
      begin
        AX:=ABS(X);
        Y:=3.75/AX;
        BX:=Exp(AX)/Sqrt(AX);
        AX:=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        BESSI1:=AX*BX
      end
    End;


 {main program}
 BEGIN

   N:=5;           {Number of integer orders}
   X:=1.2097;      {value of argument}

   writeln;
   writeln(' X =', X);
   writeln;

   For i:=0 to N-1 do
   begin
     Y:=BESSK(i,X);
     writeln(' N = ', i);
     writeln(' Y =', Y);
     writeln
   end;

   ReadKey;
   DoneWinCrt

 END.

{End of file Tbessk.pas}