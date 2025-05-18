{************************************************************************
*   Eigenvalues and Eigenvectors of a general complex square matrix     *
*               using the QR method for complex matrices                *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*                                                                       *
* Example #1:                                                           *
*                                                                       *
* Input matrix (real symmetric for comparison with Jacobi):             *
*        1,  2,  3, -7,   12                                            *
*        2,  4,  7,  3,   -1                                            *
*        3,  7, 10,  8,    4                                            *
*       -7,  3,  8, -0.75,-9                                            *
*       12, -1,  4, -9,   10                                            *
*                                                                       *
*  Input Matrix (real part):                                            *
*    1.000000    2.000000    3.000000   -7.000000   12.000000           *
*    2.000000    4.000000    7.000000    3.000000   -1.000000           *
*    3.000000    7.000000   10.000000    8.000000    4.000000           *
*   -7.000000    3.000000    8.000000   -0.750000   -9.000000           *
*   12.000000   -1.000000    4.000000   -9.000000   10.000000           *
*                                                                       *
* Error code = 0                                                        *
*                                                                       *
* Eigenvalues:                                                          *
*    4.63349527981456E-0001   0.00000000000000E+0000                    *
*   -7.77457972947830E+0000   0.00000000000000E+0000                    *
*   -1.04865451655156E+0001   0.00000000000000E+0000                    *
*    1.82918206016855E+0001   0.00000000000000E+0000                    *
*    2.37559547653270E+0001   0.00000000000000E+0000                    *
*                                                                       *
*  Eigenvectors (real part in lines):                                   *
*    0.208812    1.000000   -0.478027   -0.275000   -0.216915           *
*    1.000000   -0.326100    0.209620   -0.147689   -0.815422           *
*    0.467172   -0.007947   -0.507827    1.000000    0.264432           *
*    0.093053    0.594911    1.000000    0.455666    0.050740           *
*    0.705820   -0.012677    0.131486   -0.527499    1.000000           *
*                                                                       *
* Example #2:                                                           *
*                                                                       *
* Input matrix (real part):                                             *
*        1,  5,  3, -7,   12                                            *
*        2,  4,  7,  3,   -1                                            *
*        3,  7, 10,  8,    4                                            *
*       -7,  3,  8, -0.75,-9                                            *
*       12, -1,  4, -9,   10                                            *
* Input matrix (imag. part):                                            *
*        1,  1,  1,  1,    1                                            *
*        1,  1,  1,  1,    1                                            *
*        1,  1,  1,  1,    1                                            *
*        1,  1,  1,  1,    1                                            *
*        1,  1,  1,  1,    1                                            *
*                                                                       *
*  Input Matrix (real  part):                                           *
*    1.000000    5.000000    3.000000   -7.000000   12.000000           *
*    2.000000    4.000000    7.000000    3.000000   -1.000000           *
*    3.000000    7.000000   10.000000    8.000000    4.000000           *
*   -7.000000    3.000000    8.000000   -0.750000   -9.000000           *
*   12.000000   -1.000000    4.000000   -9.000000   10.000000           *
*                                                                       *
*  Input Matrix (imag. part):                                           *
*    1.000000    1.000000    1.000000    1.000000    1.000000           *
*    1.000000    1.000000    1.000000    1.000000    1.000000           *
*    1.000000    1.000000    1.000000    1.000000    1.000000           *
*    1.000000    1.000000    1.000000    1.000000    1.000000           *
*    1.000000    1.000000    1.000000    1.000000    1.000000           *
*                                                                       *
* Error code = 0                                                        *
*                                                                       *
* Eigenvalues:                                                          *
*    8.88713348723696E-0001   3.47070002209858E-0002                    *
*   -8.28426842306568E+0000   2.49816236188664E-0002                    *
*   -1.03636733606943E+0001   8.79708207660509E-0001                    *
*    1.88947919026401E+0001   3.15558227567470E+0000                    *
*    2.31144365323962E+0001   9.05020892824933E-0001                    *
*                                                                       *
*  Eigenvectors (real  part in lines):                                  *
*    0.411502    1.000000   -0.532508   -0.217945   -0.417001           *
*    1.000000   -0.324189    0.269140   -0.277587   -0.872694           *
*    0.451233   -0.007289   -0.515266    1.000000    0.272732           *
*    0.179884    0.608599    1.000000    0.456406    0.069197           *
*    0.682657   -0.106663   -0.004456   -0.615170    1.000000           *
*                                                                       *
*  Eigenvectors (imag. part in lines):                                  *
*    0.000770   -0.000000   -0.034750    0.009062   -0.005181           *
*   -0.000000    0.001676    0.051622   -0.074119   -0.037648           *
*    0.043492   -0.070147   -0.050515   -0.000000   -0.066367           *
*   -0.198515    0.075502   -0.000000    0.271030   -0.314928           *
*    0.044071    0.193806    0.292159    0.172780   -0.000000           *
*                                                                       *
* --------------------------------------------------------------------- *
* Reference:   From Numath Library By Tuan Dang Trong in Fortran 77     *
*              [BIBLI 18].                                              *
*                                                                       *
*                               TPW Release 1.1 By J-P Moreau, Paris    *
*                                         (www.jpmoreau.fr)             *
* --------------------------------------------------------------------- *
* Release 1.1: added procedure MTPRINT and Example #2 (12/04/2001).     *
************************************************************************}
PROGRAM TEST_CEIGEN;

Uses WinCrt;     {specific Borland TPW }

Const NM=50;         {maximum size}

Type
      pMAT = ^MAT;
      MAT = Array[0..NM,0..NM] of Double;
      pVEC = ^VEC;
      VEC = Array[1..NM] of Double;

Var
      IERR: Integer;   {error code (must be zero) }
      N   : Integer;   {size of complex matrix}
      AR,AI: pMAT;     {pointers to real, imag. parts of complex matrix A}
      ZR,ZI: pMAT;     {complex matrix to store complex eigenvectors}
      WR,WI: pVEC;     {real, imag. parts of complex eigenvalues}
      WORK: pVEC;      {working space, size = N}
      AVEC: Boolean;   {TRUE = eigenvectors asked for} 

      I,J:  Integer;   {loop index}

      {print a matrix with caption}
      Procedure MPRINT(N,MMAX,M:Integer;Name:String;A:pMAT);
      Begin
        Writeln;
        Writeln(Name);
        For I:=1 to N do
        begin
          For J:=1 to N do write(A^[I,J]:12:6);
          writeln
        end
      End;

      {print transpose matrix with caption}
      Procedure MTPRINT(N,MMAX,M:Integer;Name:String;A:pMAT);
      Begin
        Writeln;
        Writeln(Name);
        For I:=1 to N do
        begin
          For J:=1 to N do write(A^[J,I]:12:6);
          writeln
        end
      End;

      {for debug only: print a vector with caption
      Procedure VPRINT(N:Integer;Name:String;V:pVEC);
      Begin
        Writeln(Name);
        For I:=1 to N do
        begin
          write(V^[I]:12:6);
          if (I MOD 5) =0 then Writeln
        end
      End; }

      {Headers of procedures defined below}
      Procedure CBALAN(NM,N:Integer; Var AR,AI: pMAT;
                       Var LOW,UPP:Integer; Var D:pVEC); Forward;
      Procedure CUNITH(NM,N:Integer; LOW,UPP:Integer; Var AR,AI:pMAT;
                       Var DR, DI:pVEC); Forward;
      Procedure COMQR2(NM,N,LOW,UPP:Integer;Var GR, GI:pVEC; Var HR, HI:pMAT;
                       Var WR,WI:pVEC;AVEC:Boolean; Var ZR,ZI:pMAT;
                       Var IERR:Integer); Forward;
      Procedure CBALBK(NM,N,LOW,UPP,M:Integer;D:pVEC; Var Z,T:pMAT);
                       Forward;
      Procedure NORMAL(NM,N:Integer; Var WR,WI:pVEC; AVEC:Boolean;
                       Var ZR,ZI:pMAT; TR,TI:pVEC); Forward;
      Function  ABSCD (AR,AI:Double):Double; Forward;
      Procedure DCSQRT(X,Y:Double; Var A,B:Double); Forward;
      Procedure DIVCD (AR,AI,BR,BI:Double; Var ZR,ZI:Double); Forward;
      Procedure EXCHG1(M,NM,N:Integer; Var A,B:pMAT;D:pVEC;J,K,L:Integer); Forward;
      Function  DCABS (XX,YY:Double): Double; Forward;

      {Utility function}
      Function IMIN(a,b:Integer):Integer;
      Begin
        if a<b then IMIN := a
        else IMIN := b
      End;


      {main procedure}
      Procedure CEIGEN(Var AR,AI:pMAT;NM,N:Integer;AVEC:Boolean;
                       Var WR,WI:pVEC; Var ZR,ZI:pMAT; Var WORK:pVEC;
                       Var IERR: Integer);
{ --------------------------------------------------------------------
!     This Procedure calculates the eigenvalues and eigenvectors
!     of a general complex matrix. The used method begins with a
!     reduction by unitary transformations to a Hessenberg form.
!     The QR algorithm is then applied to get the eigenvalues,
!     finally the corresponding eigenvectors are found by a simple
!     substitution.
! --------------------------------------------------------------------
!     CALLING MODE:  CEIGEN (AR,AI,LMA,N,AVEC,WR,WI,ZR,ZI,WORK,IERR);
!
!     INPUTS:
!     AR    =   TABLE(LMA,N) OF REAL COEF. OF MATRIX A              R*8
!     AI    =   TABLE(LMA,N) OF IMAGINARY COEF. OF MATRIX A         R*8
!     NM    =   1ST DIMENSION OF MATRICES AR,AI,ZR & ZI
!               AS DECLARED IN CALLING PROGRAM                      I*4
!     N     =   2ND DIMENSION OF MATRICES AR,AI,ZR & ZI             I*4
!     AVEC  =   CODE TO CALCULATE OR NOT EIGENVECTORS           Boolean
!           =   FALSE: NO EIGENVECTORS
!           =   TRUE:  WITH EIGENVECTORS
!
!     OUTPUTS:
!     WR,WI =   TABLE(N) OF EIGENVALUES                             R*8
!     ZR,ZI =   TABLE(LMA,N) OF EIGENVECTORS                        R*8
!     IERR  =   ERROR CODE                                          I*4
!           =   0,  NO ERROR
!               < OU = N,  NO CONVERGENCE
!           =   9999  DIMENSION ERROR IN TABLES A OR ZR, ZI
!
!     WORKING SPACE:
!     WORK  =   TABLE(N)                                            R*8
!
!     REFERENCE:
!     J.H.WILKINSON & C.REINSCH
!     LINEAR ALGEBRA   SPRINGER-VERLAG 1971
! --------------------------------------------------------------------}
    Var  
        LOW,N1,N2,UPP: INTEGER;
        TMP1,TMP2: pVEC;
    Begin
      IF N <= NM THEN
      begin
        New(TMP1); New(TMP2);
        N1:=N+1;
        N2:=N+N1;

        CBALAN(NM,N,AR,AI,LOW,UPP,WORK);
        CUNITH(NM,N,LOW,UPP,AR,AI,TMP1,TMP2);
        COMQR2(NM,N,LOW,UPP,TMP1,TMP2,AR,AI,WR,WI,AVEC,ZR,ZI,IERR);

        IF IERR = 0 THEN
        begin
          CBALBK(NM,N,LOW,UPP,N,WORK,ZR,ZI);
          NORMAL(NM,N,WR,WI,AVEC,ZR,ZI,TMP1,TMP2)
        end
      end
      ELSE
        IERR:= 9999;
      Dispose(TMP1); Dispose(TMP2)
    End; {CEIGEN}

    Procedure CUNITH(NM,N:Integer; LOW,UPP:Integer; Var AR,AI:pMAT;Var DR,DI:pVEC);
{     UNITARY TRANSFORMATIONS TO REDUCE A COMPLEX SQUARE MATRIX
      TO A HESSENBERG FORM }
    Label Return;
    Var
      XR,XI,G,F,H,S: Double;
      J, M: Integer;
    Begin  
      If UPP-1 < LOW+1 then goto Return;
      For M:=LOW+1 to UPP-1 do
      begin
        DR^[M]:=0.0;
        DI^[M]:=0.0;
        S:=0.0;
        H:=0.0;
        For I:=M to UPP do
          S:=S+ABS(AR^[I,M-1])+ABS(AI^[I,M-1]);
        IF S <> 0.0 THEN
        begin
          For I:=UPP Downto M do
          begin
            DR^[I]:=AR^[I,M-1]/S;
            DI^[I]:=AI^[I,M-1]/S;
            H:=H+DR^[I]*DR^[I]+DI^[I]*DI^[I]
          end;
          G:=SQRT(H);
          F:=ABSCD(DR^[M],DI^[M]);
          IF F = 0.0 THEN
          begin
            DR^[M]:=G;
            AR^[M,M-1]:=S
          end
          ELSE
          begin
            H:=H+F*G;
            G:=G/F;
            DR^[M]:=(1.0+G)*DR^[M];
            DI^[M]:=(1.0+G)*DI^[M]
          end;
          For J:=M to N do
          begin
            XR:=0.0;
            XI:=0.0;
            For I:=UPP Downto M do
            begin
              XR:=XR+DR^[I]*AR^[I,J]+DI^[I]*AI^[I,J];
              XI:=XI+DR^[I]*AI^[I,J]-DI^[I]*AR^[I,J]
            end;
            XR:=XR/H;
            XI:=XI/H;
            For I:=M to UPP do
            begin
              AR^[I,J]:=AR^[I,J]-XR*DR^[I]+XI*DI^[I];
              AI^[I,J]:=AI^[I,J]-XR*DI^[I]-XI*DR^[I]
            end
          end;
          For I:=1 to UPP do
          begin
            XR:=0.0;
            XI:=0.0;
            For J:=UPP Downto M do
            begin
              XR:=XR+DR^[J]*AR^[I,J]-DI^[J]*AI^[I,J];
              XI:=XI+DR^[J]*AI^[I,J]+DI^[J]*AR^[I,J]
            end;
            XR:=XR/H;
            XI:=XI/H;
            For J:=M to UPP do
            begin
              AR^[I,J]:=AR^[I,J]-XR*DR^[J]-XI*DI^[J];
              AI^[I,J]:=AI^[I,J]+XR*DI^[J]-XI*DR^[J]
            end
          end;
          DR^[M]:=DR^[M]*S;
          DI^[M]:=DI^[M]*S;
          AR^[M,M-1]:=-G*AR^[M,M-1];
          AI^[M,M-1]:=-G*AI^[M,M-1]
        end {if S<>0}
      end; {M loop}
Return: End; {CUNITH}

    Procedure COMQR2(NM,N,LOW,UPP:Integer;Var GR, GI:pVEC; Var HR, HI:pMAT;
                     Var WR,WI:pVEC; AVEC:Boolean; Var ZR,ZI:pMAT;
                     Var IERR:Integer);
{ ---------------------------------------------------------------------
      QR ALGORITH APPLIED TO A COMPLEX HESSENBERG MATRIX
      (translated from Fortran 77).
  --------------------------------------------------------------------- }
    Label 5,7,10,15,20,35,45,50, 55, 60;
    Var
        SR,SI,TR,TI,YR,YI,XR,XI,VR,VI,MACHEP,NORM,EPS,EPS1: Double;
        I,J,K,L,M,IN0,NA,ITS:Integer;  {IN is a reserved word in Pascal}
    Begin
{     SEEK FLOATING POINT RELATIVE PRECISION FOR TYPE REAL*8 }
      EPS:=1.0;
      Repeat
        EPS:=EPS*0.5;
        EPS1:=1.0+EPS
      Until EPS1 <= 1.0;
      MACHEP:=2.0*EPS;

      IERR:=0;
      IF (AVEC) THEN
        For I:=1 to N do
        begin
          For J:=1 to N do
          begin
            ZR^[I,J]:=0.0;
            ZI^[I,J]:=0.0
          end;
          ZR^[I,I]:=1.0
        end;

{     ACCUMULATE UNITARY TRANSFORMATIONS }
      For I:=UPP-1 Downto LOW+1 do
      begin
         IF ((HR^[I,I-1] = 0.0) AND (HI^[I,I-1] = 0.0)) OR
            ((GR^[I] = 0.0) AND (GI^[I] = 0.0)) THEN GOTO 5;
         NORM:=HR^[I,I-1]*GR^[I]+HI^[I,I-1]*GI^[I];
         For K:=I+1 to UPP do
         begin
           GR^[K]:=HR^[K,I-1];
           GI^[K]:=HI^[K,I-1]
         end;
         For J:=I to UPP do
         begin
           SR:=0.0;
           SI:=0.0;
           For K:=I to UPP do
           begin
             SR:=SR+GR^[K]*ZR^[K,J]+GI^[K]*ZI^[K,J];
             SI:=SI+GR^[K]*ZI^[K,J]-GI^[K]*ZR^[K,J]
           end;
           SR:=SR/NORM;
           SI:=SI/NORM;
           For K:=I to UPP do
           begin
             ZR^[K,J]:=ZR^[K,J]+SR*GR^[K]-SI*GI^[K];
             ZI^[K,J]:=ZI^[K,J]+SR*GI^[K]+SI*GR^[K]
           end
         end;
5:    end;

{     CREATE REAL SUBDIAGONAL }
      L:=LOW+1;
      For I:=L to UPP do
      begin
         M:=IMIN(I+1,UPP);
         IF HI^[I,I-1] = 0.0 THEN GOTO 7;
         NORM:=ABSCD(HR^[I,I-1],HI^[I,I-1]);
         YR:=HR^[I,I-1]/NORM;
         YI:=HI^[I,I-1]/NORM;
         HR^[I,I-1]:=NORM;
         HI^[I,I-1]:=0.0;
         For J:=I to N do
         begin
           SI:= YR*HI^[I,J]-YI*HR^[I,J];
           HR^[I,J]:=YR*HR^[I,J]+YI*HI^[I,J];
           HI^[I,J]:=SI
         end;
         For J:=1 to M do
         begin
           SI:= YR*HI^[J,I]+YI*HR^[J,I];
           HR^[J,I]:=YR*HR^[J,I]-YI*HI^[J,I];
           HI^[J,I]:=SI
         end;
         IF (AVEC) THEN
           For J:=LOW to UPP do
           begin
             SI:= YR*ZI^[J,I]+YI*ZR^[J,I];
             ZR^[J,I]:=YR*ZR^[J,I]-YI*ZI^[J,I];
             ZI^[J,I]:=SI
           end;
7:    end;

{     STORE ROOTS ALREADY ISOLATED BY PROCEDURE BALANC }
      For I:=1 to N do
        IF (I < LOW) OR (I > UPP) THEN
        begin
          WR^[I]:=HR^[I,I];
          WI^[I]:=HI^[I,I]
        end;

      IN0:=UPP;
      TR:=0.0;
      TI:=0.0;

{     SEEK EIGENVALUES }
10:   IF IN0 < LOW THEN GOTO 50;
      ITS:=0;
      NA:=IN0-1;

{     WE LOOK AT REAL PART OF SUBDIAGONAL ELEMENT }
15:   For L:=IN0 Downto LOW do
      begin
         {here index 0 is temporarily used by Fortran! }
         SR:=MACHEP*(ABS(HR^[L-1,L-1])+ABS(HI^[L-1,L-1])
                    +ABS(HR^[L,L])+ABS(HI^[L,L]));
         IF (ABS(HR^[L,L-1]) <= SR) OR (L = LOW) THEN GOTO 20
      end;

{     SHIFTING }
20:   IF L = IN0 THEN GOTO 45;
      IF ITS = 30 THEN GOTO 55;
      IF (ITS = 10) OR (ITS = 20) THEN
      begin
{     EXTRA SHIFTING }
        SR:=ABS(HR^[IN0,NA])+ABS(HR^[NA,IN0-2]);
        SI:=0.0
      end
      ELSE
      begin
         SR:=HR^[IN0,IN0];
         SI:=HI^[IN0,IN0];
         XR:=HR^[NA,IN0]*HR^[IN0,NA];
         XI:=HI^[NA,IN0]*HR^[IN0,NA];

         IF(XR <> 0.0) OR (XI <> 0.0) THEN
         begin
           YR:=(HR^[NA,NA]-SR)*0.5;
           YI:=(HI^[NA,NA]-SI)*0.5;
           DCSQRT(YR*YR-YI*YI+XR,2.0*YR*YI+XI,VR,VI);
           IF YR*VR+YI*VI < 0.0 THEN
           begin
             VR:=-VR;
             VI:=-VI
           end;
           DIVCD(XR,XI,YR+VR,YI+VI,VR,VI);
           SR:=SR-VR;
           SI:=SI-VI
         end
      end;

      For I:=LOW to IN0 do
      begin
        HR^[I,I]:=HR^[I,I]-SR;
        HI^[I,I]:=HI^[I,I]-SI
      end;

      TR:=TR+SR;
      TI:=TI+SI;
      ITS:=ITS+1;

{     TRIANGULARY DECOMPOSITION H := Q * R }
      For I:=L+1 to IN0 do
      begin
         SR:=HR^[I,I-1];
         VR:=HR^[I-1,I-1];
         VI:=HI^[I-1,I-1];
         NORM:=SQRT(VR*VR+VI*VI+SR*SR);
         XR:=VR/NORM;
         XI:=VI/NORM;
         WR^[I-1]:=XR;
         WI^[I-1]:=XI;
         HR^[I-1,I-1]:=NORM;
         HI^[I-1,I-1]:=0.0;
         HR^[I,I-1]:=0.0;
         HI^[I,I-1]:=SR/NORM;
         For J:=I to N do
         begin
           YR:=HR^[I-1,J];
           YI:=HI^[I-1,J];
           VR:=HR^[I,J];
           VI:=HI^[I,J];
           HR^[I-1,J]:=XR*YR+XI*YI+HI^[I,I-1]*VR;
           HI^[I-1,J]:=XR*YI-XI*YR+HI^[I,I-1]*VI;
           HR^[I,J]  :=XR*VR-XI*VI-HI^[I,I-1]*YR;
           HI^[I,J]  :=XR*VI+XI*VR-HI^[I,I-1]*YI
         end
      end;

      SI:=HI^[IN0,IN0];

      IF SI = 0.0 THEN GOTO 35;
      NORM:=ABSCD(HR^[IN0,IN0],SI);
      SR:=HR^[IN0,IN0]/NORM;
      SI:=SI/NORM;
      HI^[IN0,IN0]:=0.0;
      HR^[IN0,IN0]:=NORM;
      IF IN0 <> N THEN
        For J:=IN0+1 to N do
        begin
          YR:=HR^[IN0,J];
          YI:=HI^[IN0,J];
          HR^[IN0,J]:=SR*YR+SI*YI;
          HI^[IN0,J]:=SR*YI-SI*YR
        end;

{     INVERSE OPERATION }
35:   For J:=L+1 to IN0 do
      begin
        XR:=WR^[J-1];
        XI:=WI^[J-1];
        For I:=1 to J do
        begin
          YR:=HR^[I,J-1];
          YI:=0.0;
          VR:=HR^[I,J];
          VI:=HI^[I,J];
          IF I <> J THEN
          begin
            YI:=HI^[I,J-1];
            HI^[I,J-1]:=XR*YI+XI*YR+HI^[J,J-1]*VI
          end;
          HR^[I,J-1]:=XR*YR-XI*YI+HI^[J,J-1]*VR;
          HR^[I,J]  :=XR*VR+XI*VI-HI^[J,J-1]*YR;
          HI^[I,J]  :=XR*VI-XI*VR-HI^[J,J-1]*YI
        end;
        IF(AVEC) THEN
          For I:=LOW to UPP do
          begin
            YR:=ZR^[I,J-1];
            YI:=ZI^[I,J-1];
            VR:=ZR^[I,J];
            VI:=ZI^[I,J];
            ZR^[I,J-1]:=XR*YR-XI*YI+HI^[J,J-1]*VR;
            ZI^[I,J-1]:=XR*YI+XI*YR+HI^[J,J-1]*VI;
            ZR^[I,J]  :=XR*VR+XI*VI-HI^[J,J-1]*YR;
            ZI^[I,J]  :=XR*VI-XI*VR-HI^[J,J-1]*YI
          end
      end;

      IF SI = 0.0 THEN GOTO 15;
      For I:=1 to IN0 do
      begin
         YR:=HR^[I,IN0];
         YI:=HI^[I,IN0];
         HR^[I,IN0]:=SR*YR-SI*YI;
         HI^[I,IN0]:=SR*YI+SI*YR
      end;
      IF (AVEC) THEN
        For I:=LOW to UPP do
        begin
          YR:=ZR^[I,IN0];
          YI:=ZI^[I,IN0];
          ZR^[I,IN0]:=SR*YR-SI*YI;
          ZI^[I,IN0]:=SR*YI+SI*YR
        end;
      GOTO 15;

{     ONE ROOT HAS BEEN FOUND }
45:   HR^[IN0,IN0]:=HR^[IN0,IN0]+TR;
      HI^[IN0,IN0]:=HI^[IN0,IN0]+TI;
      WR^[IN0]:=HR^[IN0,IN0];
      WI^[IN0]:=HI^[IN0,IN0];
      IN0:=NA;
      GOTO 10;

{     ALL ROOTS HAVE BEEN FOUND
      START CALCULATION OF EIGENVECTORS BY SUBSTITUTION }
50:   IF (AVEC) THEN
      begin
      NORM:=0.0;
      For I:=1 to N do
        For J:=I to N do
          NORM:=NORM+ABS(HR^[I,J])+ABS(HI^[I,J]);
      IF NORM = 0.0 THEN GOTO 60;
      For IN0:=N Downto 2 do
      begin
         XR:=WR^[IN0];
         XI:=WI^[IN0];
         HR^[IN0,IN0]:=1.0;
         HI^[IN0,IN0]:=0.0;
         For I:=IN0-1 Downto 1 do
         begin
            VR:=HR^[I,IN0];
            VI:=HI^[I,IN0];
            IF I <> IN0-1 THEN
              For J:=I+1 to IN0-1 do
              begin
                VR:=VR+HR^[I,J]*HR^[J,IN0]-HI^[I,J]*HI^[J,IN0];
                VI:=VI+HR^[I,J]*HI^[J,IN0]+HI^[I,J]*HR^[J,IN0]
              end;
            YR:=XR-WR^[I];
            YI:=XI-WI^[I];
            IF(YR = 0.0) AND (YI = 0.0) THEN YR:=NORM*MACHEP;
            DIVCD(VR,VI,YR,YI,TR,TI);
            HR^[I,IN0]:=TR;
            HI^[I,IN0]:=TI
         end
      end;
      For I:=1 to N-1 do
      begin
         IF (I < LOW) OR (I > UPP) THEN
           For J:=I+1 to N do
           begin
             ZR^[I,J]:=HR^[I,J];
             ZI^[I,J]:=HI^[I,J]
           end
      end;

{     END CALCULATION OH HESSENBERG EIGENVECTORS }
      For J:=N Downto LOW+1 do
       begin
         M:=IMIN(J-1,UPP);
         For I:=LOW to UPP do
         begin
           VR:=ZR^[I,J];
           VI:=ZI^[I,J];
           For K:=LOW to M do
           begin
             VR:=VR+ZR^[I,K]*HR^[K,J]-ZI^[I,K]*HI^[K,J];
             VI:=VI+ZR^[I,K]*HI^[K,J]+ZI^[I,K]*HR^[K,J]
           end;
           ZR^[I,J]:=VR;
           ZI^[I,J]:=VI
         end
       end
      end; {if AVEC}

{     END CALCULATION OF EIGENVECTORS OF INPUT MATRIX }
      GOTO 60;

{     NO CONVERGENCE AFTER 30 ITERATIONS! }
55:   IERR := IN0;
60: End;  {COMQR2}


    Procedure CBALAN(NM,N:Integer; Var AR,AI: pMAT;
                       Var LOW,UPP:Integer; Var D:pVEC);
{     BALANCE A COMPLEX MATRIX (VERSION DERIVED FROM BALANC).
      No effect for a symmetrical matrix).   }
    Label 10,11,20,22,30,40,50;
    Var
      B2,C,F,G,R,S: Double;
      NOCONV:Boolean;
      B,J,J1,K,K1,L,L1: Integer;
    Begin
      B:=16;
      B2:=B*B;
      L:=1;
      K:=N;
10:   For J:=K Downto 1 do
      begin
         J1:=J;
         K1:=K;
         L1:=L;
         For I:=1 to K do
         begin
           IF I <> J THEN
             IF (AR^[J,I] <> 0.0) OR (AI^[J,I] <> 0.0) THEN GOTO 11
         end;
         EXCHG1(K,NM,N,AR,AI,D,J1,K1,L1);
         K:=K-1;
         GOTO 10;
11:   end;
20:   For J:=L to K do
      begin
         J1:=J;
         K1:=K;
         L1:=L;
         For I:=L to K do
         begin
            IF I <> J THEN
              IF (AR^[I,J] <> 0.0) OR (AI^[I,J] <> 0.0) THEN GOTO 22
         end;

         EXCHG1(L,NM,N,AR,AI,D,J1,K1,L1);
         L:=L+1;
         GOTO 20;
22:   end;
      LOW:=L;
      UPP:=K;
      For I:=L to K do D^[I]:=1.0;
30:   NOCONV:=FALSE;
      For I:=L to K do
      begin
        C:=0.0;
        R:=0.0;
        For J:=L to K do
          IF J <> I THEN
          begin
            C:=C+ABS(AR^[J,I])+ABS(AI^[J,I]);
            R:=R+ABS(AR^[I,J])+ABS(AI^[I,J])
          end;
        G:=R/B;
        F:=1.0;
        S:=C+R;
40:     IF C < G THEN
        begin 
          F:=F*B;
          C:=C*B2;
          GOTO 40
        end;
        G:=R*B;
50:     IF C >= G THEN
        begin
          F:=F/B;
          C:=C/B2;
          GOTO 50
        end;
        IF (C+R)/F < 0.95*S THEN
        begin
          G:=1.0/F;
          D^[I]:=D^[I]*F;
          NOCONV:=TRUE;
          For J:=L to N do
          begin
            AR^[I,J]:=AR^[I,J]*G;
            AI^[I,J]:=AI^[I,J]*G
          end;
          For J:=1 to K do
          begin
            AI^[J,I]:=AI^[J,I]*F;
            AR^[J,I]:=AR^[J,I]*F
          end
        end
      end;
      IF (NOCONV) THEN GOTO 30
    End; {CBALAN}

    Procedure CBALBK(NM,N,LOW,UPP,M:Integer;D: pVEC; Var Z,T:pMAT);
{     This Procedure finds the eigenvectors of a complex matrix,
      as output of Procedure CBALAN.  }
    Label 10;
    Var S:Double;
        I,J,K,L:Integer;
    Begin
      For I:=LOW to UPP do
      begin
        S:=D^[I];
        For J:=1 to M do
        begin
          Z^[I,J]:=Z^[I,J]*S;
          T^[I,J]:=T^[I,J]*S
        end
      end;
      For L:=1 to N do
      begin
        I:=L;
        IF (I >= LOW) AND (I <= UPP) THEN GOTO 10;
        IF I < LOW THEN I:=LOW-L;
        K:=Round(D^[I]);
        IF K <> I THEN
          For J:=1 to M do
          begin
            S:=Z^[I,J];
            Z^[I,J]:=Z^[K,J];
            Z^[K,J]:=S;
            S:=T^[I,J];
            T^[I,J]:=T^[K,J];
            T^[K,J]:=S
          end;
10:   end
    End;

    Procedure NORMAL(NM,N:Integer; Var WR,WI:pVEC; AVEC:Boolean;
                     Var ZR,ZI:pMAT; TR,TI:pVEC);
{     Sort eigenvalues by absolute values in ascending order and
      normalize.  }
    Label 5;
    Var
      ZM,Z1,Z2,VR,VI: Double;
      IM,J,K:Integer;
    Begin
{     SORT SOLUTIONS IN ASCENDING ORDER }
      For J:=2 to N do
      begin
         VR:=WR^[J];
         VI:=WI^[J];
         IF (AVEC) THEN
           For K:=1 to N do
           begin
             TR^[K]:=ZR^[K,J];
             TI^[K]:=ZI^[K,J]
           end;

         For I:=J-1 Downto 1 do
         begin
            IF ABSCD(WR^[I],WI^[I]) <= ABSCD(VR,VI) THEN GOTO 5;
            WR^[I+1]:=WR^[I];
            WI^[I+1]:=WI^[I];
            IF (AVEC) THEN
              For K:=1 to N do
              begin
                ZR^[K,I+1]:=ZR^[K,I];
                ZI^[K,I+1]:=ZI^[K,I]
              end
         end;
         I:=0;
5:       WR^[I+1]:=VR;
         WI^[I+1]:=VI;
         IF (AVEC) THEN
           For K:=1 to N do
           begin
             ZR^[K,I+1]:=TR^[K];
             ZI^[K,I+1]:=TI^[K]
           end
      end;

{     NORMALIZE WITH RESPECT TO BIGGEST ELEMENT }
      IF (AVEC) THEN
        For J:=N Downto 1 do
        begin
          ZM:= 0.0;
          For I:=1 to N do
          begin
            VR:=ZR^[I,J];
            VI:=ZI^[I,J];
            Z1:=ABSCD (VR,VI);
            IF Z1 >= ZM THEN
            begin
              IM := I;
              ZM := Z1
            end
          end;
          Z1 := ZR^[IM,J];
          Z2 := ZI^[IM,J];
          For I := 1 to N do
          begin
            DIVCD(ZR^[I,J],ZI^[I,J],Z1,Z2,VR,VI);
            ZR^[I,J] := VR;
            ZI^[I,J] := VI
          end
        end;
    End; {NORMAL}

{ Utility routines for complex numbers }
    FUNCTION ABSCD(AR,AI:Double):Double;
{     ABSOLUTE VALUE OF A COMPLEX NUMBER C:=AR+I*AI
      ABSCD,AR,AI REAL*8
      ABSCD:=SQRT(AR**2+AI**2)  }
    Var
        XR,XI,W: Double;
    Begin
      XR:=ABS(AR);
      XI:=ABS(AI);
      IF XR <= XI THEN
      begin
        W:=XR;
        XR:=XI;
        XI:=W
      end;
      IF XI = 0.0 THEN
        ABSCD:=XR
      ELSE
        ABSCD:=XR*SQRT(1.0+(XI/XR)*(XI/XR))
    End;

    Procedure DCSQRT(X,Y: Double; Var A,B:Double);
{     SQUARE ROOT OF A COMPLEX NUMBER  A+I*B := RAC(X+I*Y) }
    Begin
      IF (X = 0.0) AND (Y = 0.0) THEN
      begin
        A:=0.0;
        B:=0.0
      end
      ELSE
      begin
        A:=SQRT((ABS(X)+DCABS(X,Y))*0.5);
        IF X >= 0.0 THEN
          B:=Y/(A+A)
        ELSE
          IF Y < 0.0 THEN
            B:=-A
          ELSE
          begin
            B:=A;
            A:=Y/(B+B)
          end
      end
    End;

{     ABSOLUTE VALUE OF A COMPLEX NUMBER X+I*Y }

    FUNCTION DCABS(XX,YY:Double): Double;
    Var
        X,Y,W: Double;
    Begin
      X:=ABS(XX);
      Y:=ABS(YY);
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
      DCABS:=W
    End;

    Procedure DIVCD (AR,AI,BR,BI:Double; Var ZR,ZI:Double);
{     EFFECT DIVISION Z:=(ZR+I*ZI):=(AR+I*AI)/(BR+I*BI)
 -----DO NOT USE IF BR:=BI:=0.  }
    Var
        YR,YI,W: Double;
    Begin
      YR:=BR;
      YI:=BI;
      IF ABS(YR) > ABS(YI) THEN
      begin
        W:=YI/YR;
        YR:=W*YI+YR;
        ZR:=(AR+W*AI)/YR;
        ZI:=(AI-W*AR)/YR
      end
      ELSE
      begin
        W:=YR/YI;
        YI:=W*YR+YI;
        ZR:=(W*AR+AI)/YI;
        ZI:=(W*AI-AR)/YI
      end
    End;

{ Exchange elements of a complex matrix }
    Procedure EXCHG1(M,NM,N:Integer; Var A,B:pMAT;D:pVEC;J,K,L:Integer);
    Var F: Double;
    Begin
      D^[M]:=J;
      IF J <> M THEN
      begin
        For I:=1 to K do
        begin
          F:=A^[I,J];
          A^[I,J]:=A^[I,M];
          A^[I,M]:=F;
          F:=B^[I,J];
          B^[I,J]:=B^[I,M];
          B^[I,M]:=F
        end;
        For I:=L to N do
        begin
          F:=A^[J,I];
          A^[J,I]:=A^[M,I];
          A^[M,I]:=F;
          F:=B^[J,I];
          B^[J,I]:=B^[M,I];
          B^[M,I]:=F
        end
      end
    End;


{main program}
BEGIN

  N:=5;   {size of input complex matrix}

  {dynamic allocations}
  New(AR);New(AI); New(ZR); New(ZI); New(WR); New(WI); New(WORK);

  For I:=0 to NM do
    For J:=0 to NM do
    begin
      AR^[I,J] := 0.0;
      AI^[I,J] := 0.0
    end;

  For I:=1 to NM do WORK^[I] := 0.0;

  {data section #1}
  AR^[1,1] := 1.0; AR^[1,2] := 2.0; AR^[1,3] :=  3.0; AR^[1,4] := -7.0; AR^[1,5] := 12.0;
  AR^[2,1] := 2.0; AR^[2,2] := 4.0; AR^[2,3] :=  7.0; AR^[2,4] :=  3.0; AR^[2,5] := -1.0;
  AR^[3,1] := 3.0; AR^[3,2] := 7.0; AR^[3,3] := 10.0; AR^[3,4] :=  8.0; AR^[3,5] :=  4.0;
  AR^[4,1] :=-7.0; AR^[4,2] := 3.0; AR^[4,3] :=  8.0; AR^[4,4] :=-0.75; AR^[4,5] := -9.0;
  AR^[5,1] :=12.0; AR^[5,2] :=-1.0; AR^[5,3] :=  4.0; AR^[5,4] := -9.0; AR^[5,5] := 10.0;

  {data section #2
  AR^[1,1] := 1.0; AR^[1,2] := 5.0; AR^[1,3] :=  3.0; AR^[1,4] := -7.0; AR^[1,5] := 12.0;
  AR^[2,1] := 2.0; AR^[2,2] := 4.0; AR^[2,3] :=  7.0; AR^[2,4] :=  3.0; AR^[2,5] := -1.0;
  AR^[3,1] := 3.0; AR^[3,2] := 7.0; AR^[3,3] := 10.0; AR^[3,4] :=  8.0; AR^[3,5] :=  4.0;
  AR^[4,1] :=-7.0; AR^[4,2] := 3.0; AR^[4,3] :=  8.0; AR^[4,4] :=-0.75; AR^[4,5] := -9.0;
  AR^[5,1] :=12.0; AR^[5,2] :=-1.0; AR^[5,3] :=  4.0; AR^[5,4] := -9.0; AR^[5,5] := 10.0;

  For I:=1 to N do
    For J:=1 to N do
      AI^[I,J] := 1.0;  }

  MPRINT(N,NM,N,'  Input Matrix (real part):', AR);
 {#2 only:
  MPRINT(N,NM,N,'  Input Matrix (imag. part):', AI); }

  {call main subroutine}
  CEIGEN(AR,AI,NM,N,TRUE,WR,WI,ZR,ZI,WORK,IERR);

  writeln;
  writeln('  Error code = ', IERR);
  writeln;

  if IERR = 0 then  {no error}
  begin
    writeln('  Eigenvalues:');
    for I:=1 to N do
      writeln('   ',WR^[i],'  ',WI^[i]);

      MTPRINT(N,NM,N,'  Eigenvectors (real part in lines):', ZR);
     {#2 only:
      MTRINT(N,NM,N,'  Eigenvectors (imag. part in lines):', ZI);  }

  end;

  {ending section}
  ReadKey;
  {free memory}
  Dispose(AR);Dispose(AI); Dispose(ZR); Dispose(ZI); Dispose(WR);
  Dispose(WI); Dispose(WORK);
  {close application (Borland TPW only) }
  DoneWinCrt

END. {of main program

end of file tceigen.pas}