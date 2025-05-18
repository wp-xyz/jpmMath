{******************************************************
*    LU decomposition routines called by test_lu.pas  *
*                                                     *
*              TPW UNIT version by J-P Moreau, Paris  *
*                       (www.jpmoreau.fr)             *
* --------------------------------------------------- *
* Reference:                                          *
*                                                     *
* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
*  University Press, 1986" [BIBLI 08].                *
******************************************************}
UNIT LU;

INTERFACE

USES Basis;

PROCEDURE LUDCMP(A:pVECT; N:INTEGER; VAR INDX:pIVECT; VAR D:INTEGER; VAR CODE:INTEGER);

PROCEDURE LUBKSB(A:pVECT; N:INTEGER; INDX: pIVECT; VAR B:pVECT);

IMPLEMENTATION

{  ***************************************************************
   * Given an N x N matrix A, this routine replaces it by the LU *
   * decomposition of a rowwise permutation of itself. A and N   *
   * are input. INDX is an output vector which records the row   *
   * permutation effected by the partial pivoting; D is output   *
   * as -1 or 1, depending on whether the number of row inter-   *
   * changes was even or odd, respectively. This routine is used *
   * in combination with LUBKSB to solve linear equations or to  *
   * invert a matrix. Return code is 1, if matrix is singular.   *
   ***************************************************************  }
 PROCEDURE LUDCMP;
 Const NMAX=100; TINY=1.5e-16;
 VAR
        AMAX,DUM, SUM : DOUBLE;
        I,IMAX,J,K,N1 : INTEGER;
        VV : pVECT;

 BEGIN

 New(VV);

 D:=1; CODE:=0; N1:=N+1;

 for  I:=1 to N do
 begin
   AMAX:=0.0;
   for J:=1 to N do
     IF ABS(A^[I*N1+J]) > AMAX then AMAX:=ABS(A^[I*N1+J]);

   IF(AMAX < TINY) THEN
   begin
     CODE := 1;
     exit
   end;
   VV^[I] := 1.0 / AMAX
 end; { i loop }

 for J:=1 to N do
 begin
   for I:=1 to J-1 do
   begin
     SUM := A^[I*N1+J];
     for K:=1 to I-1 do
       SUM := SUM - A^[I*N1+K]*A^[K*N1+J];
     A^[I*N1+J] := SUM
   END; { i loop }
   AMAX := 0.0;
   for I:=J to N do
   begin
     SUM := A^[I*N1+J];
     for  K:=1 to J-1 do
       SUM := SUM - A^[I*N1+K]*A^[K*N1+J];
     A^[I*N1+J] := SUM;
     DUM := VV^[I]*ABS(SUM);
     IF DUM >= AMAX THEN
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
   IF ABS(A^[J*N1+J]) < TINY then
     A^[J*N1+J] := TINY;

   IF J <> N THEN
   begin
     DUM := 1.0 / A^[J*N1+J];
     for I:=J+1 to N do
       A^[I*N1+J] := A^[I*N1+J]*DUM
   end 
 END; { j loop }

 Dispose(VV)

 END; { subroutine LUDCMP }


{  ******************************************************************
   * Solves the set of N linear equations A . X := B.  Here A is    *
   * input, not as the matrix A but rather as its LU decomposition, *
   * determined by the routine LUDCMP. INDX is input as the permuta-*
   * tion vector returned by LUDCMP. B is input as the right-hand   *
   * side vector B, and returns with the solution vector X. A, N and*
   * INDX are not modified by this routine and can be used for suc- *
   * cessive calls with different right-hand sides. This routine is *
   * also efficient for plain matrix inversion.                     *
   ******************************************************************  }
 PROCEDURE LUBKSB;
 VAR SUM : DOUBLE;
     I,II,J,LL,N1 : INTEGER;

 BEGIN
 II := 0; N1:=N+1;

 for I:=1 to N do
 begin
   LL := INDX^[I];
   SUM := B^[LL];
   B^[LL] := B^[I];
   IF II <>0 THEN
     for J:=II to I-1 do
       SUM := SUM - A^[I*N1+J]*B^[J]
   ELSE IF SUM <>0.0 THEN
     II := I;
   B^[I] := SUM
 end; { i loop }

 for I:=N downto 1 do
 begin
   SUM := B^[I];
   IF I < N THEN
     for J:=I+1 to N do
       SUM := SUM - A^[I*N1+J]*B^[J];
   B^[I] := SUM / A^[I*N1+I]
 end; { i loop }


 END; { subroutine LUBKSB }


END. {of unit lu.pas

End of file lu.pas}