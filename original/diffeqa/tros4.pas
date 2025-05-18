{***************************************************************
* NUMERICAL SOLUTION OF A STIFF SYSTEM OF FIRST ORDER ORDINARY *
* DIFFERENTIAL EQUATIONS Y'=F(X,Y) BY ROSENBROCK METHOD.       *
* ------------------------------------------------------------ *
* SAMPLE RUN:                                                  *
* Example #1:                                                  *
* (Solve set of differential equations (N=2):                  *
*  F(1) = Y(1) * Y(2) + COS(X) - HALF * SIN(TWO * X)           *
*  F(2) = Y(1) * Y(1) + Y(2) * Y(2) - (ONE + SIN(X))           *
*  Find values of F(1), F(2) at X=1.5].                        *
*                                                              *
*  SOLUTION AT X= 1.50000000000000E+0000                       *
*  Y(1) =   1.23600333231324E+0000                             *
*  Y(2) =  -1.04930010453074E-0001                             *
*                                                              *
*  LAST STEP SIZE = 2.83434084350562E-0004                     *
*  ERROR CODE = 1                                              *
*                                                              *
* Example #2:                                                  *
* (Solve set of differential equations (N=5):                  *
*  F(1) = Y(2)                                                 *
*  F(2) = Y(3)                                                 *
*  F(3) = Y(4)                                                 *
*  F(4) = Y(5)                                                 *
*  F(5) = (45.0 * Y(3) * Y(4) * Y(5) -                         *
*          40.0 * Y(4) * Y(4) * Y(4)) / (NINE * Y(3) * Y(3))   *
*  Find values of F(1), F(2), ..., F(5) at X=1.5).             *
*                                                              *
*  SOLUTION AT X= 1.50000000000000E+0000                       *
*  Y(1) =   4.36396121915058E+0000                             *
*  Y(2) =   4.00000034474812E+0000                             *
*  Y(3) =   2.82842769366278E+0000                             *
*  Y(4) =   6.62548799912012E-0007                             *
*  Y(5) =  -3.77123765984714E+0000                             *
*                                                              *
*  LAST STEP SIZE =  3.97038928337456E-0002                    *
*  ERROR CODE = 1                                              *
* ------------------------------------------------------------ *
* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77   *
*       [BIBLI 18].                                            * 
*                                                              *
*                       TPW Release 1.1 By J-P Moreau, Paris   *
*                                (www.jpmoreau.fr)             *
* ------------------------------------------------------------ *
* Release 1.1: corected bug in GRK4T (11/26/04).               *
****************************************************************
  LIST OF USED SUBROUTINES (HERE INCLUDED]:
  ========================================
  ROS4, RO4COR, SHAMP, GRK4A, GRK4T, VELDD, VELDS, LSTAB,
  DECA, DECB, SOL, SOLB.
----------------------------------------------------------------
Note: the banded matrix branch has not been tested here, but is
fully implemented. } 
PROGRAM TROS4;

Uses WinCrt;

Const
      NMX  = 25;
      HALF = 0.5;
      ONE  = 1.0;
      TWO  = 2.0;
      NINE = 9.0;
      ZERO = 0.0;

Type
      pVECT = ^VECT;
      VECT = Array[1..NMX] of Double;
      pIVECT = ^IVECT;
      IVECT = Array[1..NMX] of Integer;

      MAT5 = Array[1..5,1..5] of Double;

Var
      X,XEND, H,RTOL,ATOL: Double;
      I,IDID,N: Integer;
      Y, WORK: pVECT;
      IWORK: pIVECT;

      {variables for statistics (optional use) }
      NFCN,NSTEP,NJAC,NACCPT,NREJCT,NDEC,NSOL: LongInt;

      IDFX,IFCN,IJAC,IMAS,IOUT,ITOL,MLJAC,MLMAS: Integer;
      LE1,LIWORK,LJAC,LMAS,LWORK,MUJAC,MUMAS: Integer;


{define example #1}
Procedure FCN(N:Integer;X:Double;Y:pVECT; Var F:pVect);
Begin
  F^[1] := Y^[1] * Y^[2] + COS(X) - HALF * SIN(TWO * X);
  F^[2] := Y^[1] * Y^[1] + Y^[2] * Y^[2] - (ONE + SIN(X))
End;

{define example #2
Procedure FCN(N:Integer;X:Double;Y:pVECT; Var F:pVect);
Begin
  F^[1] := Y^[2];
  F[2] := Y^[3];
  F[3] := Y^[4];
  F[4] := Y^[5];
  F[5] := (45.0 * Y^[3] * Y^[4] * Y^[5] -
           40.0 * Y^[4] * Y^[4] * Y^[4]) / (NINE * Y^[3] * Y^[3])
End; }


Function IMAX(a,b:Integer): Integer;
Begin
  if a > b then IMAX := a
           else IMAX := b
End;

Function IMIN(a,b:Integer): Integer;
Begin
  if a < b then IMIN := a
           else IMIN := b
End;

Function MAX(a,b:Double): Double;
Begin
  if a > b then MAX := a
           else MAX := b
End;

Function MIN(a,b:Double): Double;
Begin
  if a < b then MIN := a
           else MIN := b
End;

{calculate y^x }
Function Power(y,x:Double): Double;
Begin
  IF x<0 THEN EXIT;
  Power:=Exp(x*Ln(y))
End;

Function SIGN(a,b:Double): Double;
Begin
  if b<0 then SIGN:=-ABS(a)
         else SIGN:=ABS(a)
End;

Procedure RO4COR(N:Integer;Var X:Double; Var Y: pVECT;
                 Var XEND:Double; HMAX: Double; Var H: Double;
                 RTOL,ATOL:Double; ITOL,IJAC,MLJAC,MUJAC: Integer;
                 IDFX,MLMAS,MUMAS,IOUT:Integer; Var IDID: Integer;
                 NMAX:LongInt; UROUND:Double; METH: Integer;
                 FAC1,FAC2,FACREJ:Double; AUTNMS,IMPLCT,BANDED:Boolean;
                 LFJAC,LE,LDMAS:Integer; YNEW,DY1,DY,AK1,AK2,AK3,
                 AK4,FX,FJAC1,EE1,FMAS1:pVECT; IP:pIVECT); Forward;


{*********************************************************************}
Procedure ROS4 (N,IFCN:Integer;Var X:Double;Var Y:pVECT;
                Var XEND,H:Double; RTOL,ATOL:Double;
                ITOL,IJAC,MLJAC,MUJAC,IDFX:Integer;
                IMAS,MLMAS,MUMAS,IOUT:Integer;
                WORK:pVECT;LWORK:Integer;IWORK:pIVECT;
                LIWORK:Integer; Var IDID:Integer);
{ ---------------------------------------------------------------------
{     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
{     SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
{     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4
{     WITH STEP SIZE CONTROL.
{     C.F. SECTION IV.7
{
{     AUTHORS: E. HAIRER AND G. WANNER
{              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
{              CH-1211 GENEVE 24, SWITZERLAND
{              E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
{
{     THIS CODE IS PART OF THE BOOK:
{         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
{         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
{         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
{         SPRINGER-VERLAG [1990]
{
{     VERSION OF OCTOBER 12, 1990
{
{     INPUT PARAMETERS
{     ----------------
{     N           DIMENSION OF THE SYSTEM
{
{     FCN         NAME OF PROCEDURE COMPUTING THE
{                 VALUE OF F(X,Y):
{                    Procedure FCN(N,X,Y,F)
{                    X,Y[N],F[N]: Double;
{                    F[1]:=...   ETC.
{                 (Here removed from parameter list}
{
{     IFCN        GIVES INFORMATION ON FCN:
{                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
{                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
{
{     X           INITIAL X-VALUE
{
{     Y[N]        INITIAL VALUES FOR Y
{
{     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
{
{     H           INITIAL STEP SIZE GUESS;
{                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
{                 H:=1.D0/(NORM OF F'), USUALLY 1e-2 OR 1e-3, IS GOOD.
{                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
{                 ADAPTS ITS STEP SIZE. STUDY THE CHOSEN VALUES FOR A FEW
{                 STEPS IN PROCEDURE SOLOUT, WHEN YOU ARE NOT SURE.
{                 (IF H=0.0, THE CODE PUTS H=1e-6).
{
{     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
{                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
{                 (multiple tolerances not implemented here).
{
{     ITOL        SWITCH FOR RTOL AND ATOL:
{                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
{                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
{                     Y[I] BELOW RTOL*ABS[Y^[I]]+ATOL
{                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
{                     THE CODE KEEPS THE LOCAL ERROR OF Y[I] BELOW
{                     RTOL[I]*ABS[Y^[I]]+ATOL[I].
{
{     JAC         NAME OF THE PROCEDURE WHICH COMPUTES
{                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
{                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1.
{                 FOR IJAC=1, THIS PROCEDURE MUST HAVE THE FORM:
{                    Procedure JAC(N,X,Y,DFY,LDFY)
{                    X,Y[N],DFY[LDFY,N]:Double;
{                    DFY[1,1]:= ...
{                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
{                 FURNISHED BY THE CALLING PROGRAM.
{                 IF MLJAC = N, THE JACOBIAN IS SUPPOSED TO
{                    BE FULL AND THE PARTIAL DERIVATIVES ARE
{                    STORED IN DFY AS
{                       DFY(I,J) := PARTIAL F[I] / PARTIAL Y^[J]
{                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
{                    THE PARTIAL DERIVATIVES ARE STORED
{                    DIAGONAL-WISE AS
{                       DFY(I-J+MUJAC+1,J) := PARTIAL F[I] / PARTIAL Y^[J].
{
{     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
{                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
{                       DIFFERENCES, PROCEDURE JAC IS NEVER CALLED.
{                    IJAC=1: JACOBIAN IS SUPPLIED BY PROCEDURE JAC.
{
{     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
{                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
{                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
{                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
{                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
{                       THE MAIN DIAGONAL).
{
{     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
{                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
{                 NEED NOT BE DEFINED IF MLJAC=N.
{
{     DFX         NAME OF THE PROCEDURE WHICH COMPUTES
{                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
{                 THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
{                 OTHERWISE, THIS PROCEDURE MUST HAVE THE FORM:
{                    Procedure DFX(N,X,Y,FX)
{                    X,Y[N],FX[N]:Double;
{                    FX[1]:= ...
{
{     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX:
{                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE
{                       DIFFERENCES, PROCEDURE DFX IS NEVER CALLED.
{                    IDFX=1: DF/DX IS SUPPLIED BY PROCEDURE DFX.
{
{     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
{     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
{
{     MAS         NAME OF PROCEDURE COMPUTING THE MASS-MATRIX M.
{                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
{                 MATRIX AND NEEDS NOT TO BE DEFINED.
{                 IF IMAS=1, THE PROCEDURE MAS IS OF THE FORM:
{                    Procedure MAS(N,AM,LMAS)
{                    AM(LMAS,N):Double;
{                    AM[1,1]:= ....
{                    IF MLMAS = N, THE MASS-MATRIX IS STORED
{                    AS FULL MATRIX LIKE
{                         AM[I,J] := M[I,J]
{                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
{                    DIAGONAL-WISE AS
{                         AM[I-J+MUMAS+1,J] := M[I,J].
{
{     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
{                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
{                       MATRIX, PROCEDURE MAS IS NEVER CALLED.
{                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
{
{     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
{                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
{                    ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
{                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
{                    MATRIX (>:= NUMBER OF NON-ZERO DIAGONALS BELOW
{                    THE MAIN DIAGONAL).
{                 MLMAS IS SUPPOSED TO BE <= MLJAC.
{
{     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>:= NUMBER OF NON-
{                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
{                 NEED NOT BE DEFINED IF MLMAS=N.
{                 MUMAS IS SUPPOSED TO BE <= MUJAC.
{
{     SOLOUT      NAME OF PROCEDURE PROVIDING THE
{                 NUMERICAL SOLUTION DURING INTEGRATION.
{                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
{                 IT MUST HAVE THE FORM:
{                    Procedure SOLOUT (NR,XOLD,X,Y,N,IRTRN)
{                    X,Y[N]:Double;
{                    ....
{                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
{                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
{                    THE FIRST GRID-POINT).
{                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
{                    IS SET <0, ROS4 RETURNS TO THE CALLING PROGRAM.
{
{     IOUT        GIVES INFORMATION ON THE PROCEDURE SOLOUT:
{                    IOUT=0: Procedure IS NEVER CALLED
{                    IOUT=1: Procedure IS USED FOR OUTPUT
{
{     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
{                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
{                 "LWORK" MUST BE AT LEAST
{                             N*[LJAC+LMAS+LE1+8]+5
{                 WHERE
{                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
{                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JACOBIAN)
{                 AND
{                    LMAS=0              IF IMAS=0
{                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
{                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
{                 AND
{                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
{                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JACOBIAN).
{
{                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
{                 MASS-MATRIX IS THE INDENTITY (IMAS:=0), THE MINIMUM
{                 STORAGE REQUIREMENT IS:
{                             LWORK := 2*N*N+8*N+5.
{
{     LWORK       DECLARED LENGHT OF ARRAY "WORK".
{
{     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
{                 "LIWORK" MUST BE AT LEAST N+2.
{
{     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
{
{ ----------------------------------------------------------------------
{
{     SOPHISTICATED SETTING OF PARAMETERS
{     -----------------------------------
{              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
{              WELL. THEY MAY BE DEFINED BY SETTING WORK[1],..,WORK[5]
{              AS WELL AS IWORK[1],IWORK[2] DIFFERENT FROM ZERO.
{              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
{
{    IWORK[1]  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
{              THE DEFAULT VALUE (FOR IWORK[1]=0) IS 100000.
{
{    IWORK[2]  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
{              IF IWORK[2] = 1,  METHOD OF SHAMPINE
{              IF IWORK[2] = 2,  METHOD GRK4T OF KAPS-RENTROP
{              IF IWORK[2] = 3,  METHOD GRK4A OF KAPS-RENTROP
{              IF IWORK[2] = 4,  METHOD OF VAN VELDHUIZEN [GAMMA:=1/2]
{              IF IWORK[2] = 5,  METHOD OF VAN VELDHUIZEN ["D-STABLE"]
{              IF IWORK[2] = 6,  AN L-STABLE METHOD
{              THE DEFAULT VALUE (FOR IWORK[2]=0) IS IWORK[2]=2.
{
{    WORK[1]   UROUND, THE ROUNDING UNIT, DEFAULT 1e-16.
{
{    WORK[2]   MAXIMAL STEP SIZE, DEFAULT XEND-X.
{
{    WORK[3], WORK[4]   PARAMETERS FOR STEP SIZE SELECTION
{              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
{                 WORK[3] <= HNEW/HOLD <= WORK[4]
{              DEFAULT VALUES: WORK[3]=0.2, WORK[4]=6.0
{
{    WORK[5]   AVOID THE HUMP: AFTER TWO CONSECUTIVE STEP REJECTIONS
{              THE STEP SIZE IS MULTIPLIED BY WORK[5]
{              DEFAULT VALUES: WORK[5]=0.1
{
{-----------------------------------------------------------------------
{
{     OUTPUT PARAMETERS
{     -----------------
{     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
{                 (AFTER SUCCESSFUL RETURN X=XEND)
{
{     Y[N]        SOLUTION AT X
{
{     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
{
{     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
{                   IDID=1,  COMPUTATION SUCCESSFUL,
{                   IDID=-1, COMPUTATION UNSUCCESSFUL.
{
{ ---------------------------------------------------------
{ *** *** *** *** *** *** *** *** *** *** *** *** ***
{          DECLARATIONS
{ *** *** *** *** *** *** *** *** *** *** *** *** *** }
Label return;
Var
      AUTNMS,IMPLCT,JBAND,ARRET: Boolean;
{ --------------------------------------------------------------------
{ --- THESE COMMON VARIABLES CAN BE USED FOR STATISTICS:
{ ---    NFCN      NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
{                  EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
{ ---    NJAC      NUMBER OF JACOBIAN EVALUATIONS [EITHER ANALYTICALLY
{                  OR NUMERICALLY]
{ ---    NSTEP     NUMBER OF COMPUTED STEPS
{ ---    NACCPT    NUMBER OF ACCEPTED STEPS
{ ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
{                  HAS BEEN ACCEPTED)
{ ---    NDEC      NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
{ ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
{ -------------------------------------------------------------------}
      LDE,LDJAC,LDMAS,LDMAS2,METH: Integer;
      NMAX: LongInt;
      FAC1,FAC2,FACREJ,HMAX,UROUND: Double;
      IEYNEW,IEDY1,IEDY,IEAK1,IEAK2,IEAK3,IEAK4,
      IEFX,IEJAC,IEMAS, IEE, ISTORE: Integer;
      IEIP: Integer;

      TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9,TMP10,TMP11: pVECT;
      ITMP: pIVECT;

Begin
{Initialize temporary vectors}
      New(TMP1); New(TMP2); New(TMP3); New(TMP4); New(TMP5); New(TMP6);
      New(TMP7); New(TMP8); New(TMP9); New(TMP10); New(TMP11); New(ITMP);
{ *** *** *** *** *** *** ***
{    SETTING THE PARAMETERS
{ *** *** *** *** *** *** *** }
      NFCN:=0;
      NJAC:=0;
      NSTEP:=0;
      NACCPT:=0;
      NREJCT:=0;
      NDEC:=0;
      NSOL:=0;
      ARRET:=FALSE;
{ -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- }
      IF IWORK^[1] = 0 THEN
         NMAX:=100000
      ELSE
      begin
         NMAX:=IWORK^[1];
         IF NMAX <=0 THEN
         begin
            WRITELN(' WRONG INPUT IWORK[1]=', IWORK^[1]);
            ARRET:=TRUE
         end
      end;
{ -------- METH   COEFFICIENTS OF THE METHOD }
      IF IWORK^[2] = 0 THEN
         METH:=2
      ELSE
      begin
         METH:=IWORK^[2];
         IF (METH <= 0) OR (METH >= 7) THEN
         begin
            WRITELN(' CURIOUS INPUT IWORK[2]=', IWORK^[2]);
            ARRET:=TRUE
         end
      end;
{ -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 }
      IF WORK^[1] = ZERO THEN
         UROUND:=1e-16
      ELSE
      begin
         UROUND:=WORK^[1];
         IF (UROUND <= 1e-14) OR (UROUND >= ONE) THEN
         begin
           WRITELN(' COEFFICIENTS HAVE 16 DIGITS, UROUND=', WORK^[1]);
           ARRET:=TRUE
         end
      end;
{ -------- MAXIMAL STEP SIZE }
      IF WORK^[2] = ZERO THEN
         HMAX:=XEND-X
      ELSE
         HMAX:=WORK^[2];
{ -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION }
      IF WORK^[3] = ZERO THEN
         FAC1:=5.0
      ELSE
         FAC1:=ONE/WORK^[3];
      IF WORK^[4] = ZERO THEN
         FAC2:=ONE/6.0
      ELSE
         FAC2:=ONE/WORK^[4];
{ -------  FACREJ    FOR THE HUMP }
      IF WORK^[5] = ZERO THEN
         FACREJ:=0.1
      ELSE
         FACREJ:=WORK^[5];
{ --------- CHECK IF TOLERANCES ARE O.K.  }
      IF ITOL = 0 THEN
          IF (ATOL <= ZERO) OR (RTOL <= 10.0*UROUND) THEN
          begin
            WRITELN(' TOLERANCES ARE TOO SMALL');
            ARRET:=TRUE
          end
      ELSE
      begin
         {Multiple tolerances not implemented here
          For I:=1 to N do
          IF (ATOL[I] <= ZERO) OR (RTOL[I] <= 10.0*UROUND) THEN
          begin
            WRITELN(' TOLERANCES[',I,'] ARE TOO SMALL');
            ARRET:=TRUE
          end }
      end;
{ *** *** *** *** *** *** *** *** *** *** *** *** ***
{         COMPUTATION OF ARRAY ENTRIES
{ *** *** *** *** *** *** *** *** *** *** *** *** ***
{ ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ?          }
      AUTNMS:=(IFCN=0);
      IMPLCT:=(IMAS<>0);
      JBAND:=(MLJAC<>N);
      ARRET:=FALSE;
{ -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ------
{ -- JACOBIAN        }
      IF (JBAND) THEN
         LDJAC:=MLJAC+MUJAC+1
      ELSE
         LDJAC:=N;
{ -- MATRIX E FOR LINEAR ALGEBRA }
      IF (JBAND) THEN
         LDE:=2*MLJAC+MUJAC+1
      ELSE
         LDE:=N;
{ -- MASS MATRIX   }
      IF (IMPLCT) THEN
      begin
          IF MLMAS <> N THEN
              LDMAS:=MLMAS+MUMAS+1
          ELSE
              LDMAS:=N;
{ ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC" }
          IF (MLMAS > MLJAC) OR (MUMAS > MUJAC) THEN
          begin
            WRITELN(' BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"');
            ARRET:=TRUE
          end
      end
      ELSE
          LDMAS:=0;

      LDMAS2:=IMAX(1,LDMAS);

{ ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- }
      IEYNEW:=6;
      IEDY1:=IEYNEW+N;
      IEDY:=IEDY1+N;
      IEAK1:=IEDY+N;
      IEAK2:=IEAK1+N;
      IEAK3:=IEAK2+N;
      IEAK4:=IEAK3+N;
      IEFX :=IEAK4+N;
      IEJAC:=IEFX +N;
      IEMAS:=IEJAC+N*LDJAC;
      IEE  :=IEMAS+N*LDMAS;
{ ------ TOTAL STORAGE REQUIREMENT ----------- }
      ISTORE:=IEE+N*LDE-1;
      IF ISTORE > LWORK THEN
      begin
        WRITELN(' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=', ISTORE);
        ARRET:=TRUE
      end;
{ ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- }
      IEIP:=3;
{ --------- TOTAL REQUIREMENT --------------- }
      ISTORE:=IEIP+N-1;
      IF ISTORE > LIWORK THEN
      begin
        WRITELN(' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=', ISTORE);
        ARRET:=TRUE
      end;
{ ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID:=-1 }
      IF (ARRET) THEN
      begin
         IDID:=-1;
         goto return
      end;
{ --- prepare arguments of RO4COR --- }
      For I:=1 to N do
        if I+IEYNEW-1 <= N then
          TMP1^[I] := WORK^[I+IEYNEW-1]
        else TMP1^[I] := ZERO;
      For I:=1 to N do
        if I+IEDY1-1 <= N then
          TMP2^[I] := WORK^[I+IEDY1-1]
        else TMP2^[I] := ZERO;
      For I:=1 to N do
        if I+IEDY-1 <= N then
          TMP3^[I] := WORK^[I+IEDY-1]
        else TMP3^[I] := ZERO;
      For I:=1 to N do
        if I+IEAK1-1 <= N then
          TMP4^[I] := WORK^[I+IEAK1-1]
        else TMP4^[I] := ZERO;
      For I:=1 to N do
        if I+IEAK2-1 <= N then
          TMP5^[I] := WORK^[I+IEAK2-1]
        else TMP5^[I] := ZERO;
      For I:=1 to N do
        if I+IEAK3-1 <= N then
          TMP6^[I] := WORK^[I+IEAK3-1]
        else TMP6^[I] := ZERO;
      For I:=1 to N do
        if I+IEAK4-1 <= N then
          TMP7^[I] := WORK^[I+IEAK4-1]
        else TMP7^[I] := ZERO;
      For I:=1 to N do
        if I+IEFX-1 <= N then
          TMP8^[I] := WORK^[I+IEFX-1]
        else TMP8^[I] := ZERO;
      For I:=1 to N do
        if I+IEJAC-1 <= N then
          TMP9^[I] := WORK^[I+IEJAC-1]
        else TMP9^[I] := ZERO;
      For I:=1 to N do
        if I+IEE-1 <= N then
          TMP10^[I] := WORK^[I+IEE-1]
        else TMP10^[I] := ZERO;
      For I:=1 to N do
        if I+IEMAS-1 <= N then
          TMP11^[I] := WORK^[I+IEMAS-1]
        else TMP11^[I] := ZERO;
      For I:=1 to N do
        if I+IEIP-1 <= N then
          ITMP^[I] := IWORK^[I+IEIP-1]
        else ITMP^[I] := 0;

{ -------- CALL TO CORE INTEGRATOR ------------ }
      RO4COR(N,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IJAC,
             MLJAC,MUJAC,IDFX,MLMAS,MUMAS,IOUT,IDID,
             NMAX,UROUND,METH,FAC1,FAC2,FACREJ,AUTNMS,
             IMPLCT,JBAND,LDJAC,LDE,LDMAS2,TMP1,TMP2,
             TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9,TMP10,
             TMP11, ITMP);

{ ----------- RETURN SECTION ----------- }
RETURN: Dispose(TMP1); Dispose(TMP2); Dispose(TMP3); Dispose(TMP4); Dispose(TMP5);
        Dispose(TMP6); Dispose(TMP7); Dispose(TMP8); Dispose(TMP9); Dispose(TMP10);
        Dispose(TMP11); Dispose(ITMP)
End; {ROS4}

Procedure SHAMP(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;
Procedure GRK4T(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;
Procedure GRK4A(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;
Procedure VELDD(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;
Procedure VELDS(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;
Procedure LSTAB(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double); Forward;

Procedure DECA(N, NDIM:Integer; Var A:MAT5; Var IP:pIVECT;
               Var IER:Integer); Forward;
Procedure SOL (N,NDIM:Integer;A:MAT5; Var B:pVECT; IP:pIVECT); Forward;
Procedure DECB(N, NDIM:Integer; Var A:MAT5; ML, MU:Integer; Var IP:pIVECT;
               Var IER:Integer); Forward;
Procedure SOLB(N,NDIM:Integer;A:MAT5;ML,MU:Integer;Var B:pVECT; IP:pIVECT); Forward;


{  --------- ... AND HERE IS THE CORE INTEGRATOR  ---------- }

Procedure RO4COR(N:Integer;Var X:Double; Var Y: pVECT;
                 Var XEND:Double; HMAX: Double; Var H: Double;
                 RTOL,ATOL:Double; ITOL,IJAC,MLJAC,MUJAC: Integer;
                 IDFX,MLMAS,MUMAS,IOUT:Integer; Var IDID: Integer;
                 NMAX:LongInt; UROUND:Double; METH: Integer;
                 FAC1,FAC2,FACREJ:Double; AUTNMS,IMPLCT,BANDED:Boolean;
                 LFJAC,LE,LDMAS:Integer; YNEW,DY1,DY,AK1,AK2,AK3,
                 AK4,FX, FJAC1,EE1,FMAS1:pVECT; IP:pIVECT);
{ ----------------------------------------------------------
{     CORE INTEGRATOR FOR ROS4
{     PARAMETERS SAME AS IN ROS4 WITH WORKSPACE ADDED
{ ----------------------------------------------------------
{         DECLARATIONS
{ ----------------------------------------------------------}
Label 1,2, 12, 14, 79, 80, Return;
Var
      I,J,K,L,MBB,MBDIAG,MBJAC,MDIAG,MDIFF,MLE,MUE: Integer;
      REJECT,RJECT2: Boolean;
      FJAC, E, FMAS: MAT5;
      A21,A31,A32,C21,C31,C32,C41,C42,C43: Double;
      B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4: Double;
      DELT,HMAXN,HOPT,POSNEG,XDELT,XXOLD,YSAFE: Double;
      IRTRN,LBEG,LEND,MD,MUJACJ,MUJACP,NSING:Integer;
      FAC,HC21,HC31,HC32,HC41,HC42,HC43,SUM: Double;
      I1,I2,IB,INFO,J1,J2,MADD: Integer;
      ERR, HD1,HD2,HD3,HD4, HNEW,S,SK: Double;

Begin
{ ------- restore parameters FJAC, E, FMAS ---------- }
      For J:=1 to N do
        For I:=1 to N do
        begin
          FJAC[I,J] := FJAC1^[I+J-1];
          E[I,J] := EE1^[I+J-1];
          FMAS[I,J] := FMAS1^[I+J-1]
        end;
{ ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) THEN MAS(N,FMAS,LDMAS);
  ---- PREPARE BANDWIDTHS ----- }
      IF (BANDED) THEN
      begin
          MLE:=MLJAC;
          MUE:=MUJAC;
          MBJAC:=MLJAC+MUJAC+1;
          MBB:=MLMAS+MUMAS+1;
          MDIAG:=MLE+MUE+1;
          MBDIAG:=MUMAS+1;
          MDIFF:=MLE+MUE-MUMAS
      end;
{ *** *** *** *** *** *** ***
{  INITIALISATIONS
{ *** *** *** *** *** *** *** }
      POSNEG:=SIGN(ONE,XEND-X);

      IF METH = 1 then SHAMP(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);
      IF METH = 2 then GRK4T(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);
      IF METH = 3 then GRK4A(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);
      IF METH = 4 then VELDS(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);
      IF METH = 5 then VELDD(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);
      IF METH = 6 then LSTAB(A21,A31,A32,C21,C31,C32,C41,C42,C43,
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4);

{ --- INITIAL PREPARATIONS --- }
      HMAXN:=MIN(ABS(HMAX),ABS(XEND-X));
      H:=MIN(MAX(1e-10,ABS(H)),HMAXN);
      H:=SIGN(H,POSNEG);
      REJECT:=FALSE;
      NSING:=0;
      IRTRN:=1;
      XXOLD:=X;

      IF (IOUT <>0) THEN {SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN)};
      IF IRTRN < 0 then  GOTO 79;
{ --- BASIC INTEGRATION STEP --- }
1:    IF (NSTEP > NMAX) OR (X+0.1*H = X) OR (ABS(H) <= UROUND) THEN GOTO 79;
      IF (X-XEND)*POSNEG+UROUND > ZERO THEN
      begin
        H:=HOPT;
        IDID:=1;
        goto Return
      end;
      HOPT:=H;
      IF (X+H-XEND)*POSNEG > ZERO THEN H:=XEND-X;

      FCN(N,X,Y,DY1);

      NFCN:=NFCN+1;
{ *** *** *** *** *** *** ***
{  COMPUTATION OF THE JACOBIAN
{ *** *** *** *** *** *** *** }
      NJAC:=NJAC+1;
      IF IJAC = 0 THEN
      begin
{ --- COMPUTE JACOBIAN MATRIX NUMERICALLY --- }
          IF (BANDED) THEN
          begin
{ --- JACOBIAN IS BANDED --- }
            MUJACP:=MUJAC+1;
            MD:=IMIN(MBJAC,N);
            For K:=1 to MD do
            begin
              J:=K;
 12:          AK2^[J]:=Y^[J];
              AK3^[J]:=SQRT(UROUND*MAX(1e-5,ABS(Y^[J])));
              Y^[J]:=Y^[J]+AK3^[J];
              J:=J+MD;
              IF J <= N THEN GOTO 12;
              FCN(N,X,Y,AK1);
              J:=K;
              LBEG:=IMAX(1,J-MUJAC);
 14:          LEND:=IMIN(N,J+MLJAC);
              Y^[J]:=AK2^[J];
              MUJACJ:=MUJACP-J;
              For L:=LBEG to LEND do
                FJAC[L+MUJACJ,J]:=(AK1^[L]-DY1^[L])/AK3^[J];
              J:=J+MD;
              LBEG:=LEND+1;
              IF J <= N THEN GOTO 14
            end
          end
          ELSE
          begin
{ --- JACOBIAN IS FULL --- }
            For I:=1 to N do
            begin
              YSAFE:=Y^[I];
              DELT:=SQRT(UROUND*MAX(1e-5,ABS(YSAFE)));
              Y^[I]:=YSAFE+DELT;
              FCN(N,X,Y,AK1);
              For J:=1 to N do
                FJAC[J,I]:=(AK1^[J]-DY1^[J])/DELT;
              Y^[I]:=YSAFE
            end; 
            MLJAC:=N
          end
      end
      ELSE
{ --- COMPUTE JACOBIAN MATRIX ANALYTICALLY ---
          JAC[N,X,Y,FJAC,LFJAC]; }
      IF (NOT AUTNMS) THEN
          IF IDFX = 0 THEN
          begin
{ --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X --- }
            DELT:=SQRT(UROUND*MAX(1e-5,ABS(X)));
            XDELT:=X+DELT;
            FCN(N,XDELT,Y,AK1);
            For J:=1 to N do
              FX^[J]:=(AK1^[J]-DY1^[J])/DELT
          end  
          ELSE
{ --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X ---
            CALL DFX[N,X,Y,FX]};
{ *** *** *** *** *** *** ***
{  COMPUTE THE STAGES
{ *** *** *** *** *** *** *** }
2:    NDEC:=NDEC+1;
      HC21:=C21/H;
      HC31:=C31/H;
      HC32:=C32/H;
      HC41:=C41/H;
      HC42:=C42/H;
      HC43:=C43/H;
      FAC:=ONE/(H*GAMMA);
      IF (IMPLCT) THEN
      begin
          IF (BANDED) THEN
          begin
{ --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX)  }
            For J:=1 to N do
            begin
              I1:=IMAX(1,MUJAC+2-J);
              I2:=IMIN(MBJAC,N+MUJAC+1-J);
              For I:=I1 to I2 do
                E[I+MLE,J]:=-FJAC[I,J]
            end;
            For J:=1 to N do
            begin
              I1:=IMAX(1,MUMAS+2-J);
              I2:=IMIN(MBB,N+MUMAS+1-J);
              For I:=I1 to I2 do
              begin
                IB:=I+MDIFF;
                E[IB,J]:=E[IB,J]+FAC*FMAS[I,J]
              end
            end;
            DECB(N,LE,E,MLE,MUE,IP,INFO);
            IF INFO <> 0 THEN GOTO 80;
            IF (AUTNMS) THEN
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
  ---   3] THE DIFFERENTIAL EQUATION IS AUTONOMOUS         }
                For I:=1 to N do  AK1^[I]:=DY1^[I];
                SOLB(N,LE,E,MLE,MUE,AK1,IP);
                For I:=1 to N do
                  YNEW^[I]:=Y^[I]+A21*AK1^[I];
                FCN(N,X,YNEW,DY);
                For I:=1 to N do
                  YNEW^[I]:=HC21*AK1^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK2^[I]:=SUM+DY^[I]
                end;
                SOLB(N,LE,E,MLE,MUE,AK2,IP);
                For I:=1 to N do
                  YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                FCN(N,X,YNEW,DY);
                For I:=1 to N do
                  YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK3^[I]:=SUM+DY^[I]
                end;
                SOLB(N,LE,E,MLE,MUE,AK3,IP);
                For I:=1 to N do
                  YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK4^[I]:=SUM+DY^[I]
                end;
                SOLB(N,LE,E,MLE,MUE,AK4,IP)
            end
            ELSE
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
  ---   3] THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS    }
                HD1:=H*D1;
                HD2:=H*D2;
                HD3:=H*D3;
                HD4:=H*D4;
                For I:=1 to N do
                  AK1^[I]:=DY1^[I]+HD1*FX^[I];
                SOLB(N,LE,E,MLE,MUE,AK1,IP);
                For I:=1 to N do
                  YNEW^[I]:=Y^[I]+A21*AK1^[I];
                FCN(N,X+C2*H,YNEW,DY);
                For I:=1 to N do
                  YNEW^[I]:=HC21*AK1^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK2^[I]:=SUM+DY^[I]+HD2*FX^[I]
                end;
                SOLB(N,LE,E,MLE,MUE,AK2,IP);
                For I:=1 to N do
                  YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                FCN(N,X+C3*H,YNEW,DY);
                For I:=1 to N do
                  YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK3^[I]:=SUM+DY^[I]+HD3*FX^[I]
                end;
                SOLB(N,LE,E,MLE,MUE,AK3,IP);
                For I:=1 to N do
                  YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                For I:=1 to N do
                begin
                  SUM:=ZERO;
                  J1:=IMAX(1,I-MLMAS);
                  J2:=IMIN(N,I+MUMAS);
                  For J:=J1 to J2 do
                    SUM:=SUM+FMAS[I-J+MBDIAG,J]*YNEW^[J];
                  AK4^[I]:=SUM+DY^[I]+HD4*FX^[I]
                end; 
                SOLB(N,LE,E,MLE,MUE,AK4,IP)
            end
          end
          ELSE
          begin
              IF MLMAS <> N THEN
              begin
{ --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX)  }
                MADD:=MUMAS+1;
                For J:=1 to N do
                  For I:=1 to N do
                    E[I,J]:=-FJAC[I,J];
                For J:=1 to N do
                begin
                  I1:=IMAX(1,J-MUMAS);
                  I2:=IMIN(N,J+MLMAS);
                  For I:=I1 to I2 do
                    E[I,J]:=E[I,J]+FAC*FMAS[I-J+MADD,J]
                end;
                DECA(N,LE,E,IP,INFO);
                IF INFO <> 0 THEN GOTO 80;
                IF (AUTNMS) THEN
                begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
  ---   3] THE DIFFERENTIAL EQUATION IS AUTONOMOUS   }
                    For I:=1 to N do AK1^[I]:=DY1^[I];
                    SOL(N,LE,E,AK1,IP);
                    For I:=1 to N do
                      YNEW^[I]:=Y^[I]+A21*AK1^[I];
                    FCN(N,X,YNEW,DY);
                    For I:=1 to N do
                      YNEW^[I]:=HC21*AK1^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK2^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK2,IP);
                    For I:=1 to N do
                      YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                    FCN(N,X,YNEW,DY);
                    For I:=1 to N do
                      YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK3^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK3,IP);
                    For I:=1 to N do
                      YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK4^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK4,IP)
                end
                ELSE
                begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
  ---   3] THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS         }
                    HD1:=H*D1;
                    HD2:=H*D2;
                    HD3:=H*D3;
                    HD4:=H*D4;
                    For I:=1 to N do
                      AK1^[I]:=DY1^[I]+HD1*FX^[I];
                    SOL(N,LE,E,AK1,IP);
                    For I:=1 to N do
                      YNEW^[I]:=Y^[I]+A21*AK1^[I];
                    FCN(N,X+C2*H,YNEW,DY);
                    For I:=1 to N do
                      YNEW^[I]:=HC21*AK1^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD2*FX^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK2^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK2,IP);
                    For I:=1 to N do
                      YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                      FCN(N,X+C3*H,YNEW,DY);
                    For I:=1 to N do
                      YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD3*FX^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK3^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK3,IP);
                    For I:=1 to N do
                      YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                    For I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD4*FX^[I];
                      J1:=IMAX(1,I-MLMAS);
                      J2:=IMIN(N,I+MUMAS);
                      For J:=J1 to J2 do
                        SUM:=SUM+FMAS[I-J+MADD,J]*YNEW^[J];
                      AK4^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK4,IP)
                end
              end
              ELSE
              begin
{ --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A FULL OR BANDED MATRIX)  }
                  IF MLJAC = N THEN
                  begin
                    For J:=1 to N do
                      for I:=1 to N do
                        E[I,J]:=FMAS[I,J]*FAC-FJAC[I,J]
                  end
                  ELSE
                  begin
                    MADD:=MUJAC+1;
                    For J:=1 to N do
                      for I:=1 to N do
                        E[I,J]:=FMAS[I,J]*FAC;
                    For J:=1 to N do
                    begin
                      I1:=IMAX(1,J-MUJAC);
                      I2:=IMIN(N,J+MLJAC);
                      For I:=I1 to I2 do
                        E[I,J]:=E[I,J]-FJAC[I-J+MADD,J]
                    end
                  end;
                  DECA(N,LE,E,IP,INFO);
                  IF INFO <> 0 THEN GOTO 80;
                  IF (AUTNMS) THEN
                  begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B IS NOT BANDED
  ---   3] THE DIFFERENTIAL EQUATION IS AUTONOMOUS  }
                    for I:=1 to N do AK1^[I]:=DY1^[I];
                    SOL(N,LE,E,AK1,IP);
                    for I:=1 to N do
                      YNEW^[I]:=Y^[I]+A21*AK1^[I];
                    FCN(N,X,YNEW,DY);
                    for I:=1 to N do
                      YNEW^[I]:=HC21*AK1^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      For J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK2^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK2,IP);
                    for I:=1 to N do
                      YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                    FCN(N,X,YNEW,DY);
                    for I:=1 to N do
                      YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      For J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK3^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK3,IP);
                    for I:=1 to N do
                      YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I];
                      For J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK4^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK4,IP)
                  end
                  ELSE
                  begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
  ---   2] THE MATRIX B IS NOT BANDED
  ---   3] THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS  }
                    HD1:=H*D1;
                    HD2:=H*D2;
                    HD3:=H*D3;
                    HD4:=H*D4;
                    for I:=1 to N do
                      AK1^[I]:=DY1^[I]+HD1*FX^[I];
                    SOL(N,LE,E,AK1,IP);
                    for I:=1 to N do
                      YNEW^[I]:=Y^[I]+A21*AK1^[I];
                    FCN(N,X+C2*H,YNEW,DY);
                    for I:=1 to N do
                      YNEW^[I]:=HC21*AK1^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD2*FX^[I];
                      For J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK2^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK2,IP);
                    for I:=1 to N do
                      YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
                    FCN(N,X+C3*H,YNEW,DY);
                    for I:=1 to N do
                      YNEW^[I]:=HC31*AK1^[I]+HC32*AK2^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD3*FX^[I];
                      For J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK3^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK3,IP);
                    for I:=1 to N do
                      YNEW^[I]:=HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
                    for I:=1 to N do
                    begin
                      SUM:=DY^[I]+HD4*FX^[I];
                      for J:=1 to N do
                        SUM:=SUM+FMAS[I,J]*YNEW^[J];
                      AK4^[I]:=SUM
                    end;
                    SOL(N,LE,E,AK4,IP)
                  end
              end
          end
      end
      ELSE
      begin
          IF (BANDED) THEN
          begin
{ --- THE MATRIX E (B:=IDENTITY, JACOBIAN A BANDED MATRIX)  }
            for J:=1 to N do
            begin
              I1:=IMAX(1,MUJAC+2-J);
              I2:=IMIN(MBJAC,N+MUJAC+1-J);
              for I:=I1 to I2 do
                E[I+MLE,J]:=-FJAC[I,J];
              E[MDIAG,J]:=E[MDIAG,J]+FAC
            end;
            DECB(N,LE,E,MLE,MUE,IP,INFO);
            IF INFO <> 0 THEN GOTO 80;
            IF (AUTNMS) THEN
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
  ---   2] THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
  ---   3] THE DIFFERENTIAL EQUATION IS AUTONOMOUS   }
              for I:=1 to N do  AK1^[I]:=DY1^[I];
              SOLB(N,LE,E,MLE,MUE,AK1,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A21*AK1^[I];
              FCN(N,X,YNEW,DY);
              for I:=1 to N do
                AK2^[I]:=DY^[I]+HC21*AK1^[I];
              SOLB(N,LE,E,MLE,MUE,AK2,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
              FCN(N,X,YNEW,DY);
              for I:=1 to N do
                AK3^[I]:=DY^[I]+HC31*AK1^[I]+HC32*AK2^[I];
              SOLB(N,LE,E,MLE,MUE,AK3,IP);
              for I:=1 to N do
                AK4^[I]:=DY^[I]+HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
              SOLB(N,LE,E,MLE,MUE,AK4,IP)
            end
            ELSE
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
  ---   2] THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
  ---   3] THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS     }
              HD1:=H*D1;
              HD2:=H*D2;
              HD3:=H*D3;
              HD4:=H*D4;
              for I:=1 to N do
                AK1^[I]:=DY1^[I]+HD1*FX^[I];
              SOLB(N,LE,E,MLE,MUE,AK1,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A21*AK1^[I];
              FCN(N,X+C2*H,YNEW,DY);
              for I:=1 to N do
                AK2^[I]:=DY^[I]+HD2*FX^[I]+HC21*AK1^[I];
              SOLB(N,LE,E,MLE,MUE,AK2,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
              FCN(N,X+C3*H,YNEW,DY);
              for I:=1 to N do
                AK3^[I]:=DY^[I]+HD3*FX^[I]+HC31*AK1^[I]+HC32*AK2^[I];
              SOLB(N,LE,E,MLE,MUE,AK3,IP);
              for I:=1 to N do
                AK4^[I]:=DY^[I]+HD4*FX^[I]+HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
              SOLB(N,LE,E,MLE,MUE,AK4,IP)
            end
          end
          ELSE
          begin
{ --- THE MATRIX E (B:=IDENTITY, JACOBIAN A FULL MATRIX)  }
            for J:=1 to N do
            begin
              for I:=1 to N do E[I,J]:=-FJAC[I,J];
              E[J,J]:=E[J,J]+FAC
            end;  
            DECA(N,LE,E,IP,INFO);
            IF INFO <> 0 THEN GOTO 80;
            IF (AUTNMS) THEN
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
  ---   2] THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
  ---   3] THE DIFFERENTIAL EQUATION IS AUTONOMOUS       }
              for I:=1 to N do  AK1^[I]:=DY1^[I];
              SOL(N,LE,E,AK1,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A21*AK1^[I];
              FCN(N,X,YNEW,DY);
              for I:=1 to N do
                AK2^[I]:=DY^[I]+HC21*AK1^[I];
              SOL(N,LE,E,AK2,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
              FCN(N,X,YNEW,DY);
              for I:=1 to N do
                AK3^[I]:=DY^[I]+HC31*AK1^[I]+HC32*AK2^[I];
              SOL(N,LE,E,AK3,IP);
              for I:=1 to N do
                AK4^[I]:=DY^[I]+HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
              SOL(N,LE,E,AK4,IP)
            end
            ELSE
            begin
{ --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
  ---   1] THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
  ---   2] THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
  ---   3] THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS   }
              HD1:=H*D1;
              HD2:=H*D2;
              HD3:=H*D3;
              HD4:=H*D4;
              for I:=1 to N do
                AK1^[I]:=DY1^[I]+HD1*FX^[I];
              SOL(N,LE,E,AK1,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A21*AK1^[I];
              FCN(N,X+C2*H,YNEW,DY);
              for I:=1 to N do
                AK2^[I]:=DY^[I]+HD2*FX^[I]+HC21*AK1^[I];
              SOL(N,LE,E,AK2,IP);
              for I:=1 to N do
                YNEW^[I]:=Y^[I]+A31*AK1^[I]+A32*AK2^[I];
              FCN(N,X+C3*H,YNEW,DY);
              for I:=1 to N do
                AK3^[I]:=DY^[I]+HD3*FX^[I]+HC31*AK1^[I]+HC32*AK2^[I];
              SOL(N,LE,E,AK3,IP);
              for I:=1 to N do
                AK4^[I]:=DY^[I]+HD4*FX^[I]+HC41*AK1^[I]+HC42*AK2^[I]+HC43*AK3^[I];
              SOL(N,LE,E,AK4,IP)
            end
          end
      end;
      NSOL:=NSOL+4;
      NFCN:=NFCN+2;
{ *** *** *** *** *** *** ***
{  ERROR ESTIMATION
{ *** *** *** *** *** *** *** }
      NSTEP:=NSTEP+1;
{ ------------ NEW SOLUTION --------------- }
      for I:=1 to N do
        YNEW^[I]:=Y^[I]+B1*AK1^[I]+B2*AK2^[I]+B3*AK3^[I]+B4*AK4^[I];
{ ------------ COMPUTE ERROR ESTIMATION ---------------- }
      ERR:=ZERO;
      for I:=1 to N do
      begin
        S:=E1*AK1^[I]+E2*AK2^[I]+E3*AK3^[I]+E4*AK4^[I];
        IF ITOL = 0 THEN
          SK:=ATOL+RTOL*MAX(ABS(Y^[I]),ABS(YNEW^[I]))
        ELSE
          {SK:=ATOL[I]+RTOL[I]*DMAX1[DABS[Y^[I]],DABS[YNEW[I]]]};
        ERR:=ERR+SQR(S/SK)
      end;
      ERR:=SQRT(ERR/N);
{ --- COMPUTATION OF HNEW
{ --- WE REQUIRE 0.2<=HNEW/H<=6.0  }
      FAC:=MAX(FAC2,MIN(FAC1,Power(ERR,(0.25/0.9))));
      HNEW:=H/FAC;
{ *** *** *** *** *** *** ***
{  IS THE ERROR SMALL ENOUGH ?
{ *** *** *** *** *** *** *** }
      IF ERR <= ONE THEN
      begin
{ --- STEP IS ACCEPTED --- }
         NACCPT:=NACCPT+1;
         for I:=1 to N do  Y^[I]:=YNEW^[I];
         XXOLD:=X;
         X:=X+H;
         IF IOUT <> 0 THEN {SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN) };
         IF IRTRN < 0 THEN GOTO 79;
         IF ABS(HNEW) > HMAXN THEN HNEW:=POSNEG*HMAXN;
         IF (REJECT) THEN HNEW:=POSNEG*MIN(ABS(HNEW),ABS(H));
         REJECT:=FALSE;
         RJECT2:=FALSE;
         H:=HNEW;
         GOTO 1
      end
      ELSE
      begin
{ --- STEP IS REJECTED --- }
         IF (RJECT2) THEN HNEW:=H*FACREJ;
         IF (REJECT) THEN RJECT2:=TRUE;
         REJECT:=TRUE;
         H:=HNEW;
         IF NACCPT >= 1 THEN  NREJCT:=NREJCT+1;
         GOTO 2
      end;
{ --- EXIT --- }
80:   WRITELN(' MATRIX E IS SINGULAR, INFO = ',INFO);
      NSING:=NSING+1;
      IF NSING >= 5 THEN GOTO 79;
      H:=H*HALF;
      GOTO 2;
79:   WRITELN;
      WRITELN(' EXIT OF ROS4 AT X=', X,'   H:=', H);
      IDID:=-1;
Return: End; {RO4COR}

      Procedure SHAMP(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                      E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
      Begin
         A21:=2.0;
         A31:=48.0/25.0;
         A32:=6.0/25.0;
         C21:=-8.0;
         C31:=372.0/25.0;
         C32:=12.0/5.0;
         C41:=-112.0/125.0;
         C42:=-54.0/125.0;
         C43:=-2.0/5.0;
         B1:=19.0/9.0;
         B2:=1.0/2.0;
         B3:=25.0/108.0;
         B4:=125.0/108.0;
         E1:=17.0/54.0;
         E2:=7.0/36.0;
         E3:=0.0;
         E4:=125.0/108.0;
         GAMMA:=0.5;
         C2:= 0.1000000000000000e+01;
         C3:= 0.6000000000000000e+00;
         D1:= 0.5000000000000000e+00;
         D2:=-0.1500000000000000e+01;
         D3:= 0.2420000000000000e+01;
         D4:= 0.1160000000000000e+00
      End;

      Procedure GRK4A(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                      E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
      Begin
        A21:= 0.1108860759493671e+01;
        A31:= 0.2377085261983360e+01;
        A32:= 0.1850114988899692e+00;
        C21:=-0.4920188402397641e+01;
        C31:= 0.1055588686048583e+01;
        C32:= 0.3351817267668938e+01;
        C41:= 0.3846869007049313e+01;
        C42:= 0.3427109241268180e+01;
        C43:=-0.2162408848753263e+01;
        B1:= 0.1845683240405840e+01;
        B2:= 0.1369796894360503e+00;
        B3:= 0.7129097783291559e+00;
        B4:= 0.6329113924050632e+00;
        E1:= 0.4831870177201765e-01;
        E2:=-0.6471108651049505e+00;
        E3:= 0.2186876660500240e+00;
        E4:=-0.6329113924050632e+00;
        GAMMA:= 0.3950000000000000e+00;
        C2:= 0.4380000000000000e+00;
        C3:= 0.8700000000000000e+00;
        D1:= 0.3950000000000000e+00;
        D2:=-0.3726723954840920e+00;
        D3:= 0.6629196544571492e-01;
        D4:= 0.4340946962568634e+00
      End;

      Procedure GRK4T(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                      E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
      Begin
        A21:= 0.2000000000000000e+01;
        A31:= 0.4524708207373116e+01;
        A32:= 0.4163528788597648e+01;
        C21:=-0.5071675338776316e+01;
        C31:= 0.6020152728650786e+01;
        C32:= 0.1597506846727117e+00;
        C41:=-0.1856343618686113e+01;
        C42:=-0.8505380858179826e+01;
        C43:=-0.2084075136023187e+01;
        B1:= 0.3957503746640777e+01;
        B2:= 0.4624892388363313e+01;
        B3:= 0.6174772638750108e+00;
        B4:= 0.1282612945269037e+01;
        E1:= 0.2302155402932996e+01;
        E2:= 0.3073634485392623e+01;
        E3:=-0.8732808018045032e+00;
        E4:=-0.1282612945269037e+01;
        GAMMA:= 0.2310000000000000e+00;
        C2:= 0.4620000000000000e+00;
        C3:= 0.8802083333333334e+00;
        D1:= 0.2310000000000000e+00;
        D2:=-0.3962966775244303e-01;
        D3:= 0.5507789395789127e+00;
        D4:=-0.5535098457052764e-01
     End;

      Procedure VELDS(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                      E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
{ --- METHOD GIVEN BY VAN VELDHUIZEN --- }
      Begin
        A21:= 0.2000000000000000e+01;
        A31:= 0.1750000000000000e+01;
        A32:= 0.2500000000000000e+00;
        C21:=-0.8000000000000000e+01;
        C31:=-0.8000000000000000e+01;
        C32:=-0.1000000000000000e+01;
        C41:= 0.5000000000000000e+00;
        C42:=-0.5000000000000000e+00;
        C43:= 0.2000000000000000e+01;
        B1:= 0.1333333333333333e+01;
        B2:= 0.6666666666666667e+00;
        B3:=-0.1333333333333333e+01;
        B4:= 0.1333333333333333e+01;
        E1:=-0.3333333333333333e+00;
        E2:=-0.3333333333333333e+00;
        E3:=-0.0000000000000000e+00;
        E4:=-0.1333333333333333e+01;
        GAMMA:= 0.5000000000000000e+00;
        C2:= 0.1000000000000000e+01;
        C3:= 0.5000000000000000e+00;
        D1:= 0.5000000000000000e+00;
        D2:=-0.1500000000000000e+01;
        D3:=-0.7500000000000000e+00;
        D4:= 0.2500000000000000e+00
     End;

     Procedure VELDD(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                     E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
{ --- METHOD GIVEN BY VAN VELDHUIZEN --- }
     Begin
        A21:= 0.2000000000000000e+01;
        A31:= 0.4812234362695436e+01;
        A32:= 0.4578146956747842e+01;
        C21:=-0.5333333333333331e+01;
        C31:= 0.6100529678848254e+01;
        C32:= 0.1804736797378427e+01;
        C41:=-0.2540515456634749e+01;
        C42:=-0.9443746328915205e+01;
        C43:=-0.1988471753215993e+01;
        B1:= 0.4289339254654537e+01;
        B2:= 0.5036098482851414e+01;
        B3:= 0.6085736420673917e+00;
        B4:= 0.1355958941201148e+01;
        E1:= 0.2175672787531755e+01;
        E2:= 0.2950911222575741e+01;
        E3:=-0.7859744544887430e+00;
        E4:=-0.1355958941201148e+01;
        GAMMA:= 0.2257081148225682e+00;
        C2:= 0.4514162296451364e+00;
        C3:= 0.8755928946018455e+00;
        D1:= 0.2257081148225682e+00;
        D2:=-0.4599403502680582e-01;
        D3:= 0.5177590504944076e+00;
        D4:=-0.3805623938054428e-01
     End;

     Procedure LSTAB(Var A21,A31,A32,C21,C31,C32,C41,C42,C43,B1,B2,B3,B4,
                     E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4:Double);
{ --- AN L-STABLE METHOD --- }
     Begin
        A21:= 0.2000000000000000e+01;
        A31:= 0.1867943637803922e+01;
        A32:= 0.2344449711399156e+00;
        C21:=-0.7137615036412310e+01;
        C31:= 0.2580708087951457e+01;
        C32:= 0.6515950076447975e+00;
        C41:=-0.2137148994382534e+01;
        C42:=-0.3214669691237626e+00;
        C43:=-0.6949742501781779e+00;
        B1:= 0.2255570073418735e+01;
        B2:= 0.2870493262186792e+00;
        B3:= 0.4353179431840180e+00;
        B4:= 0.1093502252409163e+01;
        E1:=-0.2815431932141155e+00;
        E2:=-0.7276199124938920e-01;
        E3:=-0.1082196201495311e+00;
        E4:=-0.1093502252409163e+01;
        GAMMA:= 0.5728200000000000e+00;
        C2:= 0.1145640000000000e+01;
        C3:= 0.6552168638155900e+00;
        D1:= 0.5728200000000000e+00;
        D2:=-0.1769193891319233e+01;
        D3:= 0.7592633437920482e+00;
        D4:=-0.1049021087100450e+00
      End;

    Procedure DECA(N, NDIM:Integer; Var A:MAT5; Var IP:pIVECT;
                     Var IER:Integer);
{VERSION REAL DOUBLE PRECISION OF DEC (In pascal, DEC is a reserved word).}
    Label 20,50,70, 80, return;
    Var
      NM1,K,KP1,M,I,J: Integer;
      T: Double;
{-----------------------------------------------------------------------
{  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
{  INPUTS:
{     N = ORDER OF MATRIX.
{     NDIM = DECLARED DIMENSION OF ARRAY  A .
{     A = MATRIX TO BE TRIANGULARIZED.
{  OUTPUTS:
{     A[I,J], I <= J = UPPER TRIANGULAR FACTOR, U.
{     A[I,J], I > J  = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
{     IP[K],  K < N  = INDEX OF K-TH PIVOT ROW.
{     IP[N] = (-1)^(NUMBER OF INTERCHANGES) OR O (^ for Power).
{     IER   = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
{             SINGULAR AT STAGE K.
{  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
{  DETERM[A] := IP[N]*A[1,1]*A[2,2]*...*A[N,N].
{  IF IP[N]=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
{
{  REFERENCE..
{     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
{     C.A.C.M. 15 [1972], P. 274.
{-----------------------------------------------------------------------}
    Begin
      IER := 0;
      IP^[N] := 1;
      IF N  = 1 THEN GOTO 70;
      NM1 := N - 1;
      For K := 1 to NM1 do
      begin
        KP1 := K + 1;
        M := K;
        for I := KP1 to N do
          IF ABS(A[I,K]) > ABS(A[M,K]) THEN  M := I;
        IP^[K] := M;
        T := A[M,K];
        IF M = K THEN GOTO 20;
        IP^[N] := -IP^[N];
        A[M,K] := A[K,K];
        A[K,K] := T;
20:     IF T = ZERO THEN GOTO 80;
        T := ONE / T;
        For I := KP1 to N do  A[I,K] := -A[I,K]*T;
        For J := KP1 to N do
        begin
          T := A[M,J];
          A[M,J] := A[K,J];
          A[K,J] := T;
          IF T = ZERO THEN GOTO 50;
          For I := KP1 to N do
            A[I,J] := A[I,J] + A[I,K]*T;
50:     end
      end;
70:   K := N;
      IF A[N,N] = ZERO THEN GOTO 80;
      goto RETURN;
80:   IER := K;
      IP^[N] := 0;
RETURN: End; {DECA}

    Procedure SOL (N,NDIM:Integer;A:MAT5; Var B:pVECT; IP:pIVECT);
{ VERSION REAL DOUBLE PRECISION }
    Label 50;
    Var
      NM1,K,KP1,M,I,KB,KM1: Integer;
      T: Double;
{-----------------------------------------------------------------------
{  SOLUTION OF LINEAR SYSTEM, A*X = B.
{  INPUTS:
{    N  = ORDER OF MATRIX.
{    NDIM = DECLARED DIMENSION OF ARRAY A.
{    A  = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
{    B  = RIGHT HAND SIDE VECTOR.
{    IP = PIVOT VECTOR OBTAINED FROM DECA.
{  DO NOT USE IF DECA HAS SET IER <> 0.
{  OUTPUT:
{    B  = SOLUTION VECTOR, X.
{-----------------------------------------------------------------------}
    Begin
      IF N = 1 THEN GOTO 50;
      NM1 := N - 1;
      For K := 1 to NM1 do
      begin
        KP1 := K + 1;
        M := IP^[K];
        T := B^[M];
        B^[M] := B^[K];
        B^[K] := T;
        For I := KP1 to N do
          B^[I] := B^[I] + A[I,K]*T
      end;
      For KB := 1 to NM1 do
      begin
        KM1 := N - KB;
        K := KM1 + 1;
        B^[K] := B^[K]/A[K,K];
        T := -B^[K];
        For I := 1 to KM1 do
          B^[I] := B^[I] + A[I,K]*T
      end;
50:   B^[1] := B^[1]/A[1,1]
    End;

    Procedure DECB(N, NDIM:Integer; Var A:MAT5; ML, MU:Integer; Var IP:pIVECT;
                   Var IER:Integer);
    Label 7,20,30,45, 60,70,80, Return;
    Var T: Double;
        I,IJK,J,JK,JU,K,KP1,M,MD,MD1,MDL,MM,NM1:Integer;
{----------------------------------------------------------------------------
{  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
{  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
{  INPUTS:
{     N       ORDER OF THE ORIGINAL MATRIX A.
{     NDIM    DECLARED DIMENSION OF ARRAY  A.
{     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
{                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
{                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
{                ML+1 THROUGH 2*ML+MU+1 OF A.
{     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
{     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
{  OUTPUTS:
{     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
{             THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
{     IP      INDEX VECTOR OF PIVOT INDICES.
{     IP[N]   (-1)^(NUMBER OF INTERCHANGES) OR O (^ for Power).
{     IER     = 0 IF MATRIX A IS NONSINGULAR, OR = K IF FOUND TO BE
{             SINGULAR AT STAGE K.
{  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
{  DETERM[A] := IP[N]*A[MD,1]*A[MD,2]*...*A[MD,N]  WITH MD:=ML+MU+1.
{  IF IP[N]=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
{
{  REFERENCE..
{     THIS IS A MODIFICATION OF
{     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
{     C.A.C.M. 15 [1972], P. 274.
{-----------------------------------------------------------------------}
    Begin
      IER := 0;
      IP^[N] := 1;
      MD := ML + MU + 1;
      MD1 := MD + 1;
      JU := 0;
      IF ML = 0 THEN GOTO 70;
      IF N  = 1 THEN GOTO 70;
      IF N < MU+2 THEN  GOTO 7;
      For J := MU+2 to N do
        For I := 1 to ML do
          A[I,J] := ZERO;
7:    NM1 := N - 1;
      For K := 1 to NM1 do
      begin
        KP1 := K + 1;
        M := MD;
        MDL := IMIN(ML,N-K) + MD;
        For I := MD1 to MDL do
          IF ABS(A[I,K]) > ABS(A[M,K]) THEN  M := I;
        IP^[K] := M + K - MD;
        T := A[M,K];
        IF M = MD THEN GOTO 20;
        IP^[N] := -IP^[N];
        A[M,K] := A[MD,K];
        A[MD,K] := T;
20:     IF T = ZERO THEN GOTO 80;
        T := ONE / T;
        For I := MD1 to MDL do  A[I,K] := -A[I,K]*T;
        JU := IMIN(IMAX(JU,MU+IP^[K]),N);
        MM := MD;
        IF JU < KP1 THEN GOTO 60;
        For J := KP1 to JU do
        begin
          M := M - 1;
          MM := MM - 1;
          T := A[M,J];
          IF M = MM THEN GOTO 30;
          A[M,J] := A[MM,J];
          A[MM,J] := T;
30:       IF T = ZERO THEN GOTO 45;
          JK := J - K;
          For I := MD1 to MDL do
          begin
            IJK := I - JK;
            A[IJK,J] := A[IJK,J] + A[I,K]*T
          end;
45:     end;
60:   end;
70:   K := N;
      IF A[MD,N] = ZERO THEN GOTO 80;
      goto Return;
80:   IER := K;
      IP^[N] := 0;
Return: End; {DECB}

    Procedure SOLB(N,NDIM:Integer;A:MAT5;ML,MU:Integer;Var B:pVECT; IP:pIVECT);
    Label 25, 50;
    Var T: Double;
        IMD,K,KB,KMD,LM,MD,M,MD1,MDL,MDM,NM1:Integer;
{-----------------------------------------------------------------------
{  SOLUTION OF LINEAR SYSTEM, A*X = B.
{  INPUTS:
{    N      ORDER OF MATRIX A.
{    NDIM   DECLARED DIMENSION OF ARRAY A.
{    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
{    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
{    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
{    B      RIGHT HAND SIDE VECTOR.
{    IP     PIVOT VECTOR OBTAINED FROM DECB.
{  DO NOT USE IF DECB HAS SET IER  <> 0.
{  OUTPUT:
{    B      SOLUTION VECTOR, X.
{-----------------------------------------------------------------------}
    Begin
      MD := ML + MU + 1;
      MD1 := MD + 1;
      MDM := MD - 1;
      NM1 := N - 1;
      IF ML = 0 THEN GOTO 25;
      IF N = 1 THEN GOTO 50;
      For K := 1 to NM1 do
      begin
        M := IP^[K];
        T := B^[M];
        B^[M] := B^[K];
        B^[K] := T;
        MDL := IMIN(ML,N-K) + MD;
        For I := MD1 to MDL do
        begin
          IMD := I + K - MD;
          B^[IMD] := B^[IMD] + A[I,K]*T
        end
      end;
25:   For KB := 1 to NM1 do
      begin
        K := N + 1 - KB;
        B^[K] := B^[K] / A[MD,K];
        T := -B^[K];
        KMD := MD - K;
        LM := IMAX(1,KMD+1);
        For I := LM to MDM do
        begin
          IMD := I - KMD;
          B^[IMD] := B^[IMD] + A[I,K]*T
        end
      end;
50:   B^[1] := B^[1]/A[MD,1]
    End;

{main program}
BEGIN

  {Initialize parameters (see ROS4)  }
  N:=2;           {DIMENSION OF THE SYSTEM (N=5 for Example #2) }
  IFCN:=1;        {FCN(N,X,Y,F) MAY DEPEND ON X}
  X:=ZERO;        {INITIAL X-VALUE }
  XEND:=1.5;      {FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) }
  H:=0.001;       {INITIAL STEP SIZE GUESS (0.01 FOR #2) }
  RTOL:=1e-8;     {RELATIVE ERROR TOLERANCE (HERE SCALAR) }
  ATOL:=1e-10;    {ABSOLUTE ERROR TOLERANCE (HERE SCALAR) }
  ITOL:=0;        {BOTH RTOL AND ATOL ARE SCALARS}
  IJAC:=0;        {JACOBIAN IS COMPUTED INTERNALLY BY FINITE
                   DIFFERENCES, Procedure JAC IS NEVER CALLED }
  MLJAC:=N;       {JACOBIAN IS A FULL MATRIX. THE LINEAR ALGEBRA
                   IS DONE BY FULL-MATRIX GAUSS-ELIMINATION }
  IDFX:=0;        {DF/DX IS COMPUTED INTERNALLY BY FINITE
                   DIFFERENCES, Procedure DFX IS NEVER CALLED}
  IMAS:=0;        {M IS SUPPOSED TO BE THE IDENTITY
                   MATRIX, Procedure MAS IS NEVER CALLED }
  MLMAS:=N;       {MLMAS=N: THE FULL MATRIX CASE. THE LINEAR ALGEBRA
                   IS DONE BY FULL-MATRIX GAUSS-ELIMINATION }
  IOUT:=0;        {Procedure SOLOUT IS NEVER CALLED}
  LE1:=N;         {IF MLJAC=N (FULL JACOBIAN) }
  LJAC:=N;        {IF MLJAC=N (FULL JACOBIAN) }
  LMAS:=0;        {IF IMAS=0 }

  LIWORK:= N+2;                    {DECLARED LENGTH OF ARRAY IWORK}
  LWORK := N*(LJAC+LMAS+LE1+8)+5;  {DECLARED LENGTH OF ARRAY LWORK}

  {dynamic allocations}
  New(Y); New(WORK); New(IWORK);

  For I:=1 to NMX do
  begin
    WORK^[I]:=ZERO;  {This triggers default values [see ROS4] }
    IWORK^[I]:=0
  end;

  Y^[1]:=0.5;      {INITIAL VALUES FOR Y }
  Y^[2]:=0.5;      {In #2, Y^[1] := Y^[2] := ... := Y^[5] := 1.0}

  Writeln;
  Writeln(' Computing...');

  {call Rosenbrock Procedure with appropriate parameters
   (here, FCN has been removed from input parameters). }
     ROS4(N,IFCN,X,Y,XEND,H,
          RTOL,ATOL,ITOL,
          IJAC,MLJAC,MUJAC,IDFX,
          IMAS,MLMAS,MUMAS,
          IOUT,WORK,LWORK,IWORK,LIWORK,IDID);

  {print results}
  ClrScr;
  writeln;
  writeln(' SOLUTION AT X=', X);
  for I:=1 to N do
    writeln(' Y(',I,') =  ', Y^[I]);
  writeln;
  writeln(' LAST STEP SIZE =', H);
  writeln(' ERROR CODE = ', IDID);
  writeln;
  ReadKey;
  Dispose(Y); Dispose(WORK); Dispose(IWORK);
  DoneWinCrt

END.  {of main program

end of file tros4.pas}