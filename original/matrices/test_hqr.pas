{***********************************************************************
*    Eigenvalues of a real nonsymmetric square matrix using the HQR    *
*    method of unit Linpack.                                           *
*                                                                      *
*                             Pascal version by J-P Moreau, Paris      *
*                                      (www.jpmoreau.fr)               *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                                                                      *
* Input matrix (symmetric for comparison with Jacobi):                 *
*                                                                      *
*        1,  2,  3, -7,   12                                           *
*        2,  4,  7,  3,   -1                                           *
*        3,  7, 10,  8,    4                                           *
*       -7,  3,  8, -0.75,-9                                           *
*       12, -1,  4, -9,   10                                           *
*                                                                      *
*                                                                      *
* Output file (test_hqr.lst):                                          *
*                                                                      *
* =============================================================        *
*   TESTING THE QR ALGORITHM TO FIND THE REAL OR COMPLEX               *
*     EIGENVALUES OF A REAL SQUARE NON SYMMETRIC MATRIX                *
*            IN DOUBLE PRECISION (16 digits                            *
* =============================================================        *
*                                                                      *
* GIVEN MATRIX A:                                                      *
*                                                                      *
*     1.00    2.00    3.00   -7.00   12.00                             *
*     2.00    4.00    7.00    3.00   -1.00                             *
*     3.00    7.00   10.00    8.00    4.00                             *
*    -7.00    3.00    8.00   -0.75   -9.00                             *
*    12.00   -1.00    4.00   -9.00   10.00                             *
*                                                                      *
* MATRIX A AFTER BALANC:                                               *
*                                                                      *
* 1.0000E+0000  2.0000E+0000  3.0000E+0000 -7.0000E+0000  1.2000E+0001 * 
* 2.0000E+0000  4.0000E+0000  7.0000E+0000  3.0000E+0000 -1.0000E+0000 *
* 3.0000E+0000  7.0000E+0000  1.0000E+0001  8.0000E+0000  4.0000E+0000 *
*-7.0000E+0000  3.0000E+0000  8.0000E+0000 -7.5000E-0001 -9.0000E+0000 *
* 1.2000E+0001 -1.0000E+0000  4.0000E+0000 -9.0000E+0000  1.0000E+0001 *
*                                                                      *
* MATRIX A IN HESSENBERG FORM:                                         *
*                                                                      *
* 1.0000E+0000  1.7167E+0001 -9.7385E+0000  2.6744E+0000  3.0000E+0000 *
* 1.2000E+0001  1.6083E+0001 -9.3222E+0000 -1.0084E-0001  4.0000E+0000 *
* 0.0000E+0000  3.3194E+0000 -1.1372E+0001  4.7395E+0000  1.0333E+0001 *
* 0.0000E+0000  0.0000E+0000 -1.1556E+0001  9.8935E+0000  1.5715E+0001 *
* 0.0000E+0000  0.0000E+0000  0.0000E+0000  8.5067E+0000  8.6452E+0000 *
*                                                                      *
*                                                                      *
*                                                                      *
* EIGEN VALUES, REAL AND IMAGINARY PARTS:                              *
*                                                                      *
* -1.04865451655156E+0001   0.00000000000000E+0000                     *
*                                                                      *
* -7.77457972947830E+0000   0.00000000000000E+0000                     *
*                                                                      *
*  2.37559547653270E+0001   0.00000000000000E+0000                     *
*                                                                      *
*  1.82918206016855E+0001   0.00000000000000E+0000                     *
*                                                                      *
*  4.63349527981453E-0001   0.00000000000000E+0000                     *
*                                                                      *
* End of file Test_hqr.lst.                                            *
*                                                                      *
***********************************************************************}
PROGRAM Test_Hessenberg_QR;
USES WinCrt,Type_def,LinPack;

CONST zero = 0.0;
      n    =   5;

CONST        {Test matrix of size 5} 
        Matrix: ARRAY[1..n,1..n] OF REAL_AR =
                (( 1,  2,  3, -7,   12),
                 ( 2,  4,  7,  3,   -1),
                 ( 3,  7, 10,  8,    4),
                 (-7,  3,  8, -0.75,-9), 
                 (12, -1,  4, -9,   10));


VAR
    OUT                          : TEXT;
    IprntSw                      : BOOLEAN;
    Annee,Mois,Jour,Jour_semaine : WORD;
    i,j,low,hi                   : INTEGER;
    A                            : Square_Matrix;
    d,wr,wi                      : Real_Vector;



BEGIN {main program}

IprntSw := TRUE;

IF (IprntSw) THEN
   Assign(OUT,'TEST_HQR.LST')
   ELSE
   AssignCrt(OUT);
Rewrite(OUT);

WriteLn(OUT);
WriteLn(OUT,'=============================================================');
WriteLn(OUT,'  TESTING THE QR ALGORITHM TO FIND THE REAL OR COMPLEX');
WriteLn(OUT,'    EIGENVALUES OF A REAL SQUARE NON SYMMETRIC MATRIX');
WriteLn(OUT,'           IN DOUBLE PRECISION (16 digits)');
WriteLn(OUT,'=============================================================');
WriteLn(OUT);

FOR i := 1 TO n DO
  FOR j := 1 TO n DO A[i,j] := Matrix[i,j];

WriteLn(OUT,'GIVEN MATRIX A:');
WriteLn(OUT);

FOR i := 1 TO n DO
  BEGIN
    FOR j := 1 TO n DO Write (OUT,A[i,j]:8:2);
    WriteLn(OUT);
  END;
WriteLn(OUT);

low := 1;
hi  := n;
Balanc(A,n, d, low, hi);

WriteLn;
WriteLn(' Balanc done.');
WriteLn;
WriteLn(OUT,'MATRIX A AFTER BALANC:');
WriteLn(OUT);
FOR i := 1 TO n DO
  BEGIN
    FOR j := 1 TO n DO Write (OUT,A[i,j]:13,' ');
    WriteLn(OUT);
  END;
WriteLn(OUT);

ElmHes(A,n);
WriteLn(' ElmHes done.');
WriteLn;

WriteLn(OUT,'MATRIX A IN HESSENBERG FORM:');
WriteLn(OUT);
FOR i := 1 TO n DO
  BEGIN
    IF (i < 3 ) THEN
      BEGIN
        FOR j := 1 TO n DO Write(OUT,A[i,j]:13,' ');
        WriteLn(OUT);
      END
      ELSE BEGIN
        FOR j := 1 TO i-2 DO Write(OUT, zero:13,' ');
        FOR j := i-1 TO n DO Write(OUT,A[i,j]:13,' ');
        WriteLn(OUT);
      END;
END;

HQR_MR(A,n,wr,wi);

for i:=1 to 3 do writeln(OUT);

WriteLn(' HQR done.');
WriteLn;
WriteLn(OUT,'EIGEN VALUES, REAL AND IMAGINARY PARTS:');
WriteLn(OUT);

FOR i := 1 TO n DO
  BEGIN
    WriteLn(OUT,' ',wr[i],'  ',wi[i]);
    WriteLn(OUT);
  END;
WriteLn(OUT,' End of file Test_hqr.lst.');
WriteLn(' Results in file test_hqr.lst.');
readkey;
if IprntSw then close(OUT);
Donewincrt
END.

{End of file Test_hqr.pas}