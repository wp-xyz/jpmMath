{*******************************************************
*                   UNIT LINPACK                       *
* ---------------------------------------------------- *
* Reduced english version of linpack.pas by J-P Dumont *
* to be used by demonstration program Test_hqr.pas.    *
*                                                      *
*                   TPW version by J-P Moreau, Paris   *
*                           (www.jpmoreau.fr)          *
*******************************************************}
UNIT Linpack;

INTERFACE
USES Type_def;

CONST Maxc=7; Nsol=7; Ndiag=7;

TYPE
  Square_Matrix   = ARRAY[1..Maxc, 1..Maxc]  OF REAL_AR;  {declared double in Type_def}
  Rect_Matrix     = ARRAY[1..Nsol, 1..Maxc]  OF REAL_AR;
  Matrix_n_mPlus1 = ARRAY[1..Maxc, 0..Ndiag] OF REAL_AR;
  Matrix_3x3      = ARRAY[1..3,1..3]         OF REAL_AR;
  Real_Vector     = ARRAY[1..Maxc]           OF REAL_AR;
  Integer_Vector  = ARRAY[1..Maxc]           OF INTEGER;


PROCEDURE Balanc
                 ( VAR A      : Square_Matrix;
                       n      : INTEGER;
                   VAR d      : Real_Vector;
                   VAR low,hi : INTEGER);

PROCEDURE ElmHes
                 ( VAR A      : Square_Matrix;
                       n      : INTEGER);

PROCEDURE HQR_MR
                 ( VAR A      : Square_Matrix;
                       n      : INTEGER;
               VAR wr,wi      : Real_Vector);


IMPLEMENTATION

  Function Min(x,y:INTEGER):INTEGER;
  Begin
    if x<=y then Min:=x
    else Min:=y
  End;

{--------------------------Documentation-------------------------------
   Parameters - A : Square_Matrix; n : INTEGER;
                VAR d : Real_Vector; VAR low,hi : INTEGER
***********************************************************************
* Balanc: Version 1 of 05/15/1993 translated from ALGOL by J-P Dumont *
* ------------------------------------------------------------------- *
* Ref.: "Handbook FOR automatic computation Volume II Linear Algebra  *
*        J.H. Wilkinson et C.Reinsch, Springer Verlag, 1971"          *
* ------------------------------------------------------------------- *
*   This procedure balances the elements of a non symmetric square    *
* matrix. Given a matrix n x n stored in a table of size Maxc x Maxc, *
* that matrix is replaced by a balanced matrix having the same eigen- *
* values. A symmetric matrix is not affected by this procedure.       *
* ------------------------------------------------------------------- *
* Inputs:                                                             *
*       n       Order of non symmetric square matrix A (n x n)        *
*       A       Table n*n storing the  elements of A                  *
* ------------------------------------------------------------------- *
* Outputs:                                                            *
*       A       The original A matrix is replaced by the balanced     *
*               matrix. However, it is possible to recover exactly    *
*               the original matrix.                                  *
*  low,hi       Two integers such as A[i,j] = 0 if                    *
*               (1) i > j  AND                                        *
*               (2) j =1,...,Low-1 OU i = hi+1,...,n                  *
*       D       Vector [1..n] containing enough information to trace  *
*               the done permutations and the used scale factors.     *
**********************************************************************}
PROCEDURE Balanc;
LABEL Iteration,L1,L2,L3,L4;

CONST b = 2;  {Floating point basis of used CPU}

VAR   i,j,k,l        : INTEGER;
      b2,c,f,g,r,s   : REAL_AR;
      noconv         : BOOLEAN;
{----------------------------------------------------------------------}
  PROCEDURE Exc(m : INTEGER);  {internal procedure of Balanc}
  VAR i : INTEGER;
  BEGIN
    d[m] := 1.0*j;
    IF (j<> m ) THEN
      BEGIN
        FOR i := 1 TO k DO
          BEGIN
            f := A[i,j];
            A[i,j] := A[i,m];
            A[i,m] := f;
          END; { FOR i }
        FOR i := l TO n DO
          BEGIN
            f := A[j,i];
            A[j,i] := A[m,i];
            A[m,i] := f;
          END; { FOR i}
      END; { j<>m }
  END; {Exc}
{----------------------------------------------------------------------}
BEGIN {Balanc}
b2 := b*b;
l := 1;
k := n;
{Search lines isolating an eigenvalue and shift downwards}
L1 :
  FOR j := k DOWNTO 1 DO
    BEGIN
      r := 0;
      FOR i := 1 TO j-1 DO r := r + Abs(A[j,i]);
      FOR i := j+1 TO k DO r := r + Abs(A[j,i]);
      IF (r = 0)  THEN
        BEGIN
          Exc(k);
          DEC(k);
          GOTO L1
        END;
    END; { FOR j }
{Search columns isolating an eigenvalue and shift to the left}
L2 :
  FOR j := l TO k DO
    BEGIN
      c := 0;
      FOR i := l TO j-1 DO c := c + Abs(A[i,j]);
      FOR i := j+1 TO k DO c := c + Abs(A[i,j]);
      IF (c = 0) THEN
        BEGIN
          Exc(l);
          INC(l);
          GOTO L2
        END;
    END; { FOR j }
{Now balance submatrix from ligne l to line k}
  low := l;
  hi  := k;
  FOR i := 1 TO k DO d[i] := 1;
Iteration:
  noconv := FALSE;
  FOR i := l TO k DO
    BEGIN
      c := 0;
      r := c;
      FOR j := l TO i-1 DO
        BEGIN
          c := c + Abs(A[j,i]);
          r := r + Abs(A[i,j]);
        END; { FOR j }
      FOR j := i+1 TO k DO
        BEGIN
          c := c + Abs(A[j,i]);
          r := r + Abs(A[i,j]);
        END; { FOR j }
      g := r/b;
      f := 1;
      s := c+r;
L3:
      IF (c < g) THEN
        BEGIN
          f := f*b;
          c := c*b2;
          GOTO L3;
        END;
      g := r*b;
L4:
      IF (c>=g) THEN
        BEGIN
          f := f/b;
          c := c/b2;
          GOTO L4;
        END;

{Balancing the elements of submatrix}
      IF ( (c+r)/f < 0.95*s ) THEN
        BEGIN
          g := 1/f;
          d[i] := d[i] * f;
          noconv := TRUE;
          FOR j := l TO n DO A[i,j] := A[i,j] * g;
          FOR j := 1 TO k DO A[j,i] := A[j,i] * f;
        END; {IF}

  END; {i loop}
  IF (noconv) THEN GOTO Iteration;
END; {of procedure Balanc}


{-------------------------Documentation-------------------------------}
{               VAR A: Square_Matrix; n: INTEGER);                    }
{ Reduction of a non symmetric real matrix to Hessenberg form by the  }
{ elimination method. Matrix A (n x n), stored in a table of size     }
{ Maxc x Maxc, is replaced by a superior Hessenberg matrix having the }
{ same eigenvalues. It is recommanded to call previously the Balanc   }
{ procedure. In output, the Hessenberg matrix has elements A[i,j] with}
{ i<=j+1. The elements for i>j+1, that in theory equal zero, are      }
{ actually filled with random (not used) values.                      } 
{---------------------------------------------------------------------}
{ Version 1.0 of 05/15/1993 - J.P.DUMONT - Extracted from:            }
{---------------------------------------------------------------------}
{ Reference:                                                          }
{        William H.PRESS, Brian P.FLANNERY, Saul A.TEUKOLSKY AND      }
{        William T.VETTERLING                                         }
{              N U M E R I C A L  R E C I P E S                       }
{              The Art OF Scientific Computing                        }
{              CAMBRIDGE UNIVERSITY PRESS 1987                        }
{/////////////////////////////////////////////////////////////////////}
PROCEDURE ElmHes;
VAR m,j,i : INTEGER;
    y,x   : REAL_AR;

BEGIN
  IF ( n > 2 ) THEN
  BEGIN
    FOR m := 2 TO n-1 DO
    BEGIN
      x := 0.0;
      i := m;
      FOR j := m TO n DO
      BEGIN
        IF (Abs(A[j,m-1]) > Abs(x)) THEN
        BEGIN
          x := A[j,m-1];
          i := j;
        END; {IF Abs}
      END; {FOR j:= m TO n}
      IF i<> m THEN
      BEGIN
        FOR j := m-1 TO n DO
        BEGIN
          y := A[i,j];
          A[i,j] := A[m,j];
          A[m,j] := y;
        END; {FOR j}
        FOR j := 1 TO n DO
        BEGIN
          y := A[j,i];
          A[j,i] := A[j,m];
          A[j,m] := y;
        END; {FOR j}
      END; {IF i <> m}
      IF (x <> 0.0) THEN
      BEGIN
        FOR i := m+1 TO n DO
        BEGIN
          y := A[i,m-1];
          IF y <> 0.0 THEN
          BEGIN
            y := y/x;
            A[i,m-1] := y;
            FOR j := m TO n DO A[i,j] := A[i,j] - y*A[m,j];
            FOR j := 1 TO n DO A[j,m] := A[j,m] + y*A[j,i];
          END; {IF y}
        END; {FOR i}
      END; {IF x}
    END; {FOR m}
  END; {if n>2}
END; {of PROCEDURE ElmHes}


{-------------------------Documentation--------------------------------
                   VAR A : Square_Matrix; n: INTEGER;
                   VAR wr,wi : Real_Vector);}
{/////////////////////////////////////////////////////////////////////}
{               QR algorithm for real Hessenberg matrices             }
{                                                                     }
{   This procedure finds all the eigenvalues of an upper triangle     }
{ Hessenberg matrix of size n x n, stored in a table of size Maxc x   }
{ Maxc. The input matrix A has been previously produced by the Helmes }
{ procedure. In output this A matrix is destroyed. Tables wr and wi   }
{ respectively return real parts and imaginary parts of the eigen-    }
{ values.                                                             }         
{/////////////////////////////////////////////////////////////////////}
{ Reference:                                                          }
{        William H.PRESS, Brian P.FLANNERY, Saul A.TEUKOLSKY AND      }
{        William T.VETTERLING                                         }
{              N U M E R I C A L  R E C I P E S                       }
{              The Art OF Scientific Computing                        }
{              CAMBRIDGE UNIVERSITY PRESS 1987                        }
{/////////////////////////////////////////////////////////////////////}
PROCEDURE HQR_MR;
LABEL 2,3,4;

CONST itsmx = 30;

VAR i,its,j,k,l,m,mmin,nn       : INTEGER;
    anorm,p,q,r,s,t,u,v,w,x,y,z : REAL_AR;
{---------------------------------------------------------------------}
FUNCTION Sign(a,b : REAL_AR) : REAL_AR;
BEGIN
IF (b <0.0) THEN Sign := - Abs(a)
            ELSE Sign :=   Abs(a);
END;
{---------------------------------------------------------------------}
BEGIN
 anorm := Abs(A[1,1]);
 FOR i := 2 TO n DO
  BEGIN
  FOR j := i-1 TO n DO
    BEGIN
    anorm := anorm + Abs(A[i,j]);
    END; {FOR j}
  END; {FOR i}
 nn := n;
 t := 0.0;
 WHILE (nn >= 1) DO
   BEGIN
   its := 0;
2: FOR l := nn DOWNTO 2 DO
     BEGIN
       s := Abs(A[l-1,l-1])+Abs(A[l,l]);
       IF (s = 0.0) THEN s := anorm;
       IF ((Abs(A[l,l-1])+s) = s) THEN GOTO 3
     END; {FOR l}
   l := 1;
3: writeln(' nn=',nn,'  s=',s);
   x := A[nn,nn];
   IF (l = nn) THEN
     BEGIN
       wr[nn] := x+t;
       wi[nn] := 0.0;
       DEC(nn);
     END
     ELSE BEGIN
       y := A[nn-1,nn-1];
       w := A[nn,nn-1]*a[nn-1,nn];
       IF (l = nn - 1 ) THEN
         BEGIN
           p := 0.5*(y-x);
           q := sqr(p)+w;
           z := sqrt(abs(q));
           x := x+t;
           IF (q >= 0.0) THEN
             BEGIN
               z := p + Sign(z,p);
               wr[nn] := x+z;
               wr[nn-1] := wr[nn];
               IF (z <> 0.0) THEN wr[nn] := x - w/z;
               wi[nn]   := 0.0;
               wi[nn-1] := 0.0
             END
             ELSE BEGIN
               wr[nn] := x+p;
               wr[nn-1] := wr[nn];
               wi[nn] := z;
               wi[nn-1] := -z;
             END; {ELSE}
           DEC(nn,2);
         END
         ELSE BEGIN
           IF (its = itsmx) THEN
             BEGIN
               WriteLn('Pause in HQR procedure');
               WriteLn('Too many iterations!');
               Readln;
             END; {IF its}
           IF (its = 10) OR (its = 20) THEN
             BEGIN
               t := t + x;
               FOR i := 1 TO nn DO A[i,i] := A[i,i] - x;
               s := Abs(A[nn,nn-1])+Abs(A[nn-1,nn-2]);
               x := 0.75*s;
               y := x;
               w := - 0.4375*sqr(s);
             END;
           INC(its);
           FOR m := nn-2 DOWNTO 1 DO
             BEGIN
               z := A[m,m];
               r := x - z;
               s := y - z;
               p := (r*s-w)/A[m+1,m]+ A[m,m+1];
               q := A[m+1,m+1] - z - r - s;
               r := A[m+2,m+1];
               s := Abs(p)+Abs(q)+Abs(r);
               p := p/s;
               q := q/s;
               r := r/s;
               IF ( m = 1 ) THEN GOTO 4;
               u := Abs(A[m,m-1])*(Abs(q)+Abs(r));
               v := Abs(p)*(Abs(A[m-1,m-1])+Abs(z)+Abs(A[m+1,m+1]));
               IF((u+v) = v) THEN GOTO 4;
             END; {FOR}
4:           FOR i := m+2 TO nn DO
               BEGIN
                 A[i,i-2] := 0.0;
                 IF (i <> (m+2)) THEN A[i,i-3] := 0.0;
               END; { FOR i }
             FOR k := m TO nn-1 DO
               BEGIN
                 IF (k <> m) THEN
                   BEGIN
                     p := A[k,k-1];
                     q := A[k+1,k-1];
                     r := 0.0;
                     IF (k <> (nn-1)) THEN r := A[k+2,k-1];
                     x := Abs(p)+Abs(q)+Abs(r);
                     IF (x <> 0.0) THEN
                       BEGIN
                         p := p/x;
                         q := q/x;
                         r := r/x;
                       END; { IF x}
                   END; { IF k }
                 s := Sign(sqrt(sqr(p)+sqr(q)+sqr(r)),p);
                 IF (s <> 0.0) THEN
                   BEGIN
                     IF (k = m ) THEN
                       BEGIN
                         IF ( l <> m) THEN A[k,k-1] := -A[k,k-1];
                       END
                       ELSE A[k,k-1] := -s*x;
                     p := p+s;
                     x := p/s;
                     y := q/s;
                     z := r/s;
                     q := q/p;
                     r := r/p;
                     FOR j := k TO nn DO
                       BEGIN
                       p := A[k,j]+q*A[k+1,j];
                       IF (k <> (nn-1)) THEN
                         BEGIN
                           p := p+r*A[k+2,j];
                           A[k+2,j] := A[k+2,j] -p*z;
                         END; {IF k}
                       A[k+1,j] := A[k+1,j] - p*y;
                       A[k,j] := A[k,j] - p*x;
                       END; { FOR j }
                       mmin := min(nn,k+3);
                       FOR i := 1 TO mmin DO
                         BEGIN
                           p := x*A[i,k]+ y*A[i,k+1];
                           IF (k <> (nn-1)) THEN
                             BEGIN
                               p := p+z*A[i,k+2];
                               A[i,k+2] := A[i,k+2] -p*r;
                             END;
                           A[i,k+1] := A[i,k+1] -p*q;
                           A[i,k] := A[i,k] -p;
                         END; {FOR i}
                       END;
                     END;
                   GOTO 2;
                   END;
                 END;
               END;
END; {of PROCEDURE HQR}

END.

{End of file Linpack.pas}