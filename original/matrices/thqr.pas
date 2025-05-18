{=====================================================================
*            Calculation of eigenvalues and eigenvectors             *
*               of a real non symmetric square matrix                *
*                         by QR algorithm                            * 
* ------------------------------------------------------------------ *
* SAMPLE RUN:                                                        *
*                                                                    *
* Input file Hqr5.dat contains:                                      *
*                                                                    *
* 5                                                                  *
*  1  2  3 -7    12                                                  *
*  2  4  7  3    -1                                                  *
*  3  7 10  8     4                                                  *
* -7  3  8 -0.75 -9                                                  *
* 12 -1  4 -9    10                                                  *
*                                                                    *
* Output file Hqr5.lst contains:                                     *
*                                                                    *
* --------------------------------------------------------------     *
* Eigenvalues and Eigenvectors by QR algorithm                       *
* --------------------------------------------------------------     *
* Dimension of the input matrix = 5                                  *
*                                                                    *
* Input matrix:                                                      *
*  1.000000  2.000000  3.000000 -7.000000 12.000000                  *
*  2.000000  4.000000  7.000000  3.000000 -1.000000                  *
*  3.000000  7.000000 10.000000  8.000000  4.000000                  *
* -7.000000  3.000000  8.000000 -0.750000 -9.000000                  *
* 12.000000 -1.000000  4.000000 -9.000000 10.000000                  *
*                                                                    *
* Normalized Eigenvectors:                                           *
*  0.467172  1.000000  0.705820  0.093053  0.208812                  *
* -0.007947 -0.326100 -0.012677  0.594911  1.000000                  *
* -0.507827  0.209620  0.131486  1.000000 -0.478027                  *
*  1.000000 -0.147689 -0.527499  0.455666 -0.275000                  *
*  0.264432 -0.815422  1.000000  0.050740 -0.216915                  *
*                                                                    *
* Eigenvalues:             Iterations:                               *
*-10.486545 +   0.000000 i      0                                    *
* -7.774580 +   0.000000 i      4                                    *
* 23.755955 +   0.000000 i     -4                                    *
* 18.291821 +   0.000000 i      5                                    *
*  0.463350 +   0.000000 i     -5                                    *
*                                                                    *
* Check sum =  2.78359265162582E-0013                                *
* (must be approximately 0).                                         *
* --------------------------------------------------------------     *
*                                                                    *
* Reference: "Numerical Algorithms with C By G. Engeln-Mueller and   *
*             F. Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
*                                                                    *
*                              Pascal Version By J-P Moreau, Paris.  *
*                                       (www.jpmoreau.fr)            *
*====================================================================}
PROGRAM Test_HQR;
Uses WinCrt, Type_def, Feigen0; {for procedure eigen()               }

Label    10;

VAR
         mat: Square_Matrix;   {input  matrix                        } 
         a:   Square_Matrix;   {copy of the input matrix             }
         ev:  Square_Matrix;   {Eigenvectors (if vec = TRUE)         }
         wr:  Real_Vector;     {Eigenvalues (Real part)              }
         wi:  Real_Vector;     {Eigenvalues (Imaginary parts)        }

         n:   Integer;         {Size of matrix mat                   }
         cnt: Integer_Vector;  {Iteration counter                    }
         rc:  Integer;         {Return Code                          }
         vec: Boolean;         {flag for eigenvectors (=False->none) }
     ev_norm: Boolean;        {flag for normalization of Eigenvectors}
                              {False = no normalization              } 
         i,j,k: Integer;
         v, w, norm: REAL_AR;

         fp_in, fp_out: TEXT;

BEGIN

  vec:=True; ev_norm:=True;

  Assign(fp_in,'Hqr5.dat'); Reset(fp_in);
  Assign(fp_out,'Hqr5.lst'); Rewrite(fp_out);

  Writeln(fp_out,'--------------------------------------------------------------');
  Writeln(fp_out,' Eigenvalues and Eigenvectors by QR algorithm');
  Writeln(fp_out,'--------------------------------------------------------------');

  Readln(fp_in, n);

  if n < 1 then
  begin
    Writeln(fp_out,' Dimension must be > 0');
    goto 10
  end;

  for i:=0 to n-1 do
  begin
    for j:=0 to n-2 do read(fp_in,mat[i,j]);
    readln(fp_in,mat[i,n-1])
  end;

  close(fp_in);

  writeln(fp_out,' Dimension of the input matrix = ', n);
  writeln(fp_out);
  writeln(fp_out,' Input matrix:');

  for i:=0 to n-1 do
  begin
    for j:=0 to n-1 do write(fp_out,mat[i,j]:10:6);
    writeln(fp_out)
  end;

  for i:=0 to n-1 do
    for j:=0 to n-1 do
      a[i,j]:=mat[i,j];

  writeln(fp_out);

  eigen(vec, ev_norm, n, mat, ev, wr, wi, cnt, rc);

  if rc<>0 then
  begin
    writeln(fp_out,' Error in procedure Eigen, rc=',rc);
    goto 10
  end;

  {If vec = True, print eigenvectors}
  if (vec) then
  begin
    if (ev_norm) then
      writeln(fp_out,' Normalized Eigenvectors (in columns):')
    else
      writeln(fp_out,' Not normalized Eigenvectors (in columns):');
    for i:=0 to n-1 do
    begin
      for j:=0 to n-1 do write(fp_out,ev[i,j]:10:6);
      writeln(fp_out)
    end
  end;

  writeln(fp_out);
  writeln(fp_out,' Eigenvalues:             Iterations:');
  for i := 0 to n-1 do
    writeln(fp_out, wr[i]:10:6,' + ',wi[i]:10:6,' i   ',cnt[i]:4);
  
  {Check result: sum of L1 norms of Matrix*Eigenvector - Eigenvalue*Eigenvector
   (this must be nearly 0).  }
  if (vec) then
  begin
    norm := ZERO;
    for k := 0 to n-1 do
    begin
      if wi[k] = ZERO then
      begin
        for i := 0 to n-1 do
        begin
          w := ZERO;
          for j := 0 to n-1 do
            w := w + a[i,j] * ev[j,k];
          w := w - wr[k] * ev[i,k];
          norm := norm + ABS (w)
        end
      end
      else
      begin
        for i := 0 to n-1 do
        begin
          w := ZERO;
          for j := 0 to n-1 do
            w := w + a[i,j] * ev[j,k];
          w := w - wr[k] * ev[i,k] - wi[k] * ev[i,k+1];
          v := ZERO;
          for j := 0 to n-1 do
            v := v + a[i,j] * ev[j,k+1];
          v := v - wr[k] * ev[i,k+1] + wi[k] * ev[i,k];
          norm := norm + 2.0 * SQRT (v*v + w*w)
        end;
        Inc(k)
      end
    end;
    writeln(fp_out);
    writeln(fp_out,' Check sum = ', norm);
    writeln(fp_out,' (must be approximately 0).');
  end;

10:Writeln(fp_out,'--------------------------------------------------------------');
  close(fp_out);
  writeln;
  writeln(' Results in hqr5.lst.');
  writeln;
  Readkey; DoneWinCrt

END.

{ end of file thqr.pas}