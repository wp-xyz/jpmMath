{******************************************************
*   LU decomposition routines called by test_clu.pas  *
*                                                     *
*              TPW UNIT version by J-P Moreau, Paris  *
*              (for a general complex linear system). *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* Reference:                                          *
*                                                     *
* "Numerical Recipes by W.H. Press, B. P. Flannery,   *
*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
*  University Press, 1986".                           *
*                                                     *
******************************************************}
UNIT CLU;

INTERFACE

Uses WinCrt;

CONST  ISIZE = 1024;

TYPE
       {complex number}
       COMPLEX = Record
         R, I: REAL;  {algebraic form}
       End;

       {complex vector}
       pCVEC = ^CVEC;
       CVEC =  Array[0..ISIZE] of COMPLEX;
       {integer vector}
       pIVEC = ^IVEC;
       IVEC =  Array[0..ISIZE] of Integer;


PROCEDURE LUDCMP(VAR A:pCVEC; N:INTEGER; VAR INDX:pIVEC; VAR D:INTEGER; VAR CODE:INTEGER);

PROCEDURE LUBKSB(VAR A:pCVEC; N:INTEGER; INDX: pIVEC; VAR B:pCVEC);

Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);


IMPLEMENTATION


  Function CABS(Z:Complex): Real;
  Begin
    CABS := sqrt(Z.R*Z.R+Z.I*Z.I)
  End;

  Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R*Z2.R - Z1.I*Z2.I;
    Z.I := Z1.R*Z2.I + Z1.I*Z2.R
  End;

  Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
  Var d:Real; C:Complex;
  Begin
    d := Z2.R*Z2.R+Z2.I*Z2.I;
    if d<1E-10 then
      writeln(' Complex Divide by zero !')
    else
    begin
      C.R:=Z2.R; C.I:=-Z2.I;
      CMUL(Z1,C,Z);
      Z.R:=Z.R/d; Z.I:=Z.I/d
    end
  End;


{****************************************************************
 * Given a complex N x N matrix A, this routine replaces it by  *
 * the LU decomposition of a rowwise permutation of itself. A   *
 * and N are input. INDX is an output vector which records the  *
 * row permutation effected by the partial pivoting; D is out-  *
 * put as -1 or 1, depending on whether the number of row inter-*
 * changes was even or odd, respectively. This routine is used  *
 * in combination with LUBKSB to solve linear equations or to   *
 * invert a matrix. Return code is 1, if matrix is singular.    *
 ****************************************************************}
 PROCEDURE LUDCMP;   {Release with a complex matrix stored in a complex vector}
 Const NMAX=100; TINY=2.2e-16;
 VAR
        I,IMAX,J,K,N1: INTEGER;
        VV: pCVEC;
        AMAX, CONE, DUM, SUM, TMP: COMPLEX;

 BEGIN

 New(VV);

 CONE.R:=1.0; CONE.I:=0.0;

 D:=1; CODE:=0; N1:=N+1;

 for  I:=1 to N do
 begin
   AMAX.R:=0.0; AMAX.I:=0.0;
   for J:=1 to N do
     IF CABS(A^[I*N1+J]) > CABS(AMAX) then AMAX:=A^[I*N1+J];

   IF CABS(AMAX) < TINY THEN
   begin
     CODE := 1;
     exit
   end;
   {VV^[I] := 1.0 / AMAX}
   CDIV(CONE,AMAX,VV^[I])
 end; { i loop }

 for J:=1 to N do
 begin
   for I:=1 to J-1 do
   begin
     SUM := A^[I*N1+J];
     for K:=1 to I-1 do
     begin
       {SUM = SUM - A^[I*N1+K]*A^[K*N1+J] }
       CMUL(A^[I*N1+K],A^[K*N1+J],TMP);
       SUM.R:=SUM.R - TMP.R;
       SUM.I:=SUM.I - TMP.I
     end;
     A^[I*N1+J] := SUM
   END; { i loop }
   AMAX.R:=0.0; AMAX.I:=0.0;
   for I:=J to N do
   begin
     SUM := A^[I*N1+J];
     for  K:=1 to J-1 do
     begin
       {SUM = SUM - A^[I*N1+K]*A^[K*N1+J] }
       CMUL(A^[I*N1+K],A^[K*N1+J],TMP);
       SUM.R:=SUM.R - TMP.R;
       SUM.I:=SUM.I - TMP.I
     end;
     A^[I*N1+J] := SUM;
     {DUM = VV^[I]*ABS(SUM) }
     DUM.R:=VV^[I].R * CABS(SUM);
     DUM.I:=VV^[I].I * CABS(SUM);
     IF CABS(DUM) >= CABS(AMAX) THEN
     BEGIN
       IMAX := I;
       AMAX := DUM
     END
   END; { i loop }  
   
   IF J <> IMAX THEN
   begin
     for K:=1 to N do
     begin
       DUM := A^[IMAX*N1+K];
       A^[IMAX*N1+K] := A^[J*N1+K];
       A^[J*N1+K] := DUM
     end; { k loop }
     D := -D;
     VV^[IMAX] := VV^[J]
   end;

   INDX^[J] := IMAX;
   IF CABS(A^[J*N1+J]) < TINY then
   begin
     A^[J*N1+J].R := TINY;
     A^[J*N1+J].I := 0.0
   end;

   IF J <> N THEN
   begin
     {DUM = 1.0 / A^[J*N1+J] }
     CDIV(CONE,A^[J*N1+J],DUM);
     for I:=J+1 to N do
     begin
       {A^[I*N1+J] = A^[I*N1+J]*DUM }
       CMUL(A^[I*N1+J],DUM,TMP);
       A^[I*N1+J] := TMP
     end
   end 
 END; { j loop }

 Dispose(VV)

 END; { subroutine LUDCMP }


{******************************************************************
 * Solves the set of N complex linear equations A . X := B.  Here *
 * A is input, not as the matrix A but rather as its LU decomposi-*
 * tion, determined by the routine LUDCMP. INDX is input as the   *
 * permutation vector returned by LUDCMP. B is input as the right-*
 * hand side vector B, and returns with the solution vector X. A, *
 * N and INDX are not modified by this routine and can be used for*
 * successive calls with different right-hand sides. This routine *
 * is also efficient for plain matrix inversion.                  *
 ******************************************************************
 Release with a complex matrix stored in a complex vector }
 PROCEDURE LUBKSB;
 VAR SUM, TMP: COMPLEX;
     I,II,J,LL,N1: INTEGER;

 BEGIN
 II := 0; N1:=N+1;

 for I:=1 to N do
 begin
   LL := INDX^[I];
   SUM := B^[LL];
   B^[LL] := B^[I];
   IF II <>0 THEN
     for J:=II to I-1 do
     begin
       {SUM = SUM - A^[I*N1+J]*B^[J] }
        CMUL(A^[I*N1+J],B^[J],TMP);
        SUM.R:=SUM.R - TMP.R;
        SUM.I:=SUM.I - TMP.I
     end
   ELSE IF CABS(SUM) <> 0.0 THEN
     II := I;
   B^[I] := SUM
 end; { i loop }

 for I:=N downto 1 do
 begin
   SUM := B^[I];
   IF I < N THEN
     for J:=I+1 to N do
     begin
       {SUM = SUM - A^[I*N1+J]*B^[J] }
        CMUL(A^[I*N1+J],B^[J],TMP);
        SUM.R:=SUM.R - TMP.R;
        SUM.I:=SUM.I - TMP.I
     end;
   {B^[I] = SUM / A^[I*N1+I] }
   CDIV(SUM, A^[I*N1+I], B^[I])
 end; { i loop }


 END; { subroutine LUBKSB }


END. {of unit lu.pas

End of file clu.pas}
