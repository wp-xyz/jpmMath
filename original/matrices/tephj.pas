{**********************************************************************
*     Calculate eigenvalues and eigenvectors of a Square Hermitian    *
*     Matrix By Jacobi's Method                                       *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
* File Matc4.dat contains:                                            *
* 4                                                                   *
* 12   0   -1   -1   2   0   3   3                                    *
* -1   1   12    0   1  -1  -2   0                                    *
*  2   0    1    1   8   0  -1  -1                                    *
*  3  -3   -2    0  -1   1   8   0                                    *
*                                                                     *
*  Precision: 1E-10                                                   *
*  Max. number of oterations: 25                                      *
*                                                                     *
*  Eigenvalue 1:  4.00000000000000E+0000                              *
*  Eigenvector:                                                       *
* -5.00000000000000E-0001 -  5.00000000000000E-0001 I                 *
*  1.24580090105524E-0018 +  0.00000000000000E+0000 I                 *
*  5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
*  1.00000000000000E+0000 +  0.00000000000000E+0000 I                 *
*  Eigenvalue 2:  8.00000000000000E+0000                              *
*  Eigenvector:                                                       *
* -2.70626288044616E-0017 +  0.00000000000000E+0000 I                 *
* -5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
*  1.00000000000000E+0000 +  0.00000000000000E+0000 I                 *
* -5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
*  Eigenvalue 3:  1.20000000000000E+0001                              *
*  Eigenvector:                                                       *
*  5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
*  1.00000000000000E+0000 +  0.00000000000000E+0000 I                 *
*  5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
*  2.96615611457399E-0016 +  0.00000000000000E+0000 I                 *
*  Eigenvalue 4:  1.60000000000000E+0001                              *
*  Eigenvector:                                                       *
*  1.00000000000000E+0000 +  0.00000000000000E+0000 I                 *
* -5.00000000000000E-0001 +  5.00000000000000E-0001 I                 *
* -1.11681259237675E-0016 +  0.00000000000000E+0000 I                 *
*  5.00000000000000E-0001 -  5.00000000000000E-0001 I                 *
*                                                                     *
* ------------------------------------------------------------------- *
* Reference: "ALGEBRE Algorithmes et programmes en Pascal             *
*             By Jean-Louis Jardrin - Dunod BO-PRE 1988" [BIBLI 10]   *
*                                                                     *
*                              TPW Release 1.1 By J-P Moreau, Paris.  *
*                                       (www.jpmoreau.fr)             *
*                                                                     *
* Release 1.1 (Nov. 30th, 2007): Added procedure NORMAL.              * 
**********************************************************************}                                       
Program TEPHJ;

Uses WinCrt1;

Const SIZE = 50;

Type
      Complex = Record
        R, I: Double
      End;
      pCMat= ^CMat;
      CMat = Array[1..SIZE,1..SIZE] of Complex;
      pVec = ^Vec;
       Vec = Array[1..SIZE] of Double;

Var  I,IT,J,M,N: Integer;
     dta: Double;
     A,VX: pCMat;
     R: pVec;

     fp:TEXT;  {input text file}

     {Absolute value of a complex number}
     Function CABS(C:Complex):Double;
     Begin
       CABS:=Sqrt(Sqr(C.R)+Sqr(C.I))
     End;

     {Add two complex numbers}
     Procedure CADD(c1,c2:Complex; Var c3:Complex);
     Begin
       C3.R:=C1.R+C2.R; C3.I:=C1.I+C2.I
     End;

     {Substract two complex numbers}
     Procedure CDIF(c1,c2:Complex; Var c3:Complex);
     Begin
       C3.R:=C1.R-C2.R; C3.I:=C1.I-C2.I
     End;

     {Multiply two complex numbers}
     Procedure CMUL(c1,c2:Complex; Var c3:Complex);
     Begin
       C3.R:=C1.R*C2.R - C1.I*C2.I;
       C3.I:=C1.R*C2.I + C1.I*C2.R
     End;

     {Return conjugate of a complex number}
     Procedure CONJ(C:Complex; Var C1:Complex);
     Begin
       C1.R:=C.R; C1.I:=-C.I
     End;

     {Multiply a complex number by a real number}
     Procedure CPRO(alpha:Double; C:Complex; Var C1:Complex);
     Begin
       C1.R:=alpha*C.R; C1.I:=alpha*C.I
     End;

{***********************************************************************
*   Compute all eigenvalues/eigenvectors of a square hermitian matrix  *
*   using Jacobi's method.                                             *
* -------------------------------------------------------------------- *
* Inputs:                                                              *
*         dta: required precision (double)                             *
*         M  : maximum number of iterations (integer)                  *
*         N  : size of matrix                                          *
*         A  : pointer to hermitian (complex) matrix                   *
* Outputs:                                                             *
*         it : flag for convergence (integer) =-1, no convergence      *
*         R  : vector(1..N) of eigenvalues (double)                    * 
*         VX : matrix storing complex eigenvectors in columns          *
***********************************************************************}
Procedure EPHJ(dta:Double; M,N:Integer; A:pCMat; Var it:Integer; R:pVec;
               VX:pCMat);
Var I,J,K,K0,L,L0: Integer;
    delta,s,s0,t0,t1,w0: Double;
    c0,c1,c2,c3,u0,u1,z0,z1: Complex;
Begin
  z0.R:=0.0; z0.I:=0.0; z1.R:=1.0; z1.I:=0.0;
  For I:=1 to N do
    For J:=1 to N do
      If I=J then VX^[I,J]:=z1 else VX^[I,J]:=z0;
  it:=-1; L:=1;
  Repeat
    s:=0.0;
    For I:=1 to N-1 do
      For J:=I+1 to N do
      begin
        t0:=CABS(A^[I,J]);
        if t0 > s then
        begin
          s:=t0; K0:=I; L0:=J
        end
      end;
    if s=0.0 then it:=1
    else
    begin
      delta:=Sqr(A^[L0,L0].R-A^[K0,K0].R)+4.0*Sqr(CABS(A^[K0,L0]));
      t0:=A^[L0,L0].R-A^[K0,K0].R + SQRT(delta);
      t1:=A^[L0,L0].R-A^[K0,K0].R - SQRT(delta);
      If Abs(t0) >= Abs(t1) then w0:=t0 else w0:=t1;
      s0:=Abs(w0)/SQRT(Sqr(w0)+4.0*Sqr(CABS(A^[K0,L0])));
      t0:=2.0*s0/w0; CPRO(t0,A^[K0,L0],C0);
      CONJ(C0,C1);
      For I:=1 to K0-1 do
      begin
        u0:=A^[I,K0];
        CMUL(c0,u0,c2); CPRO(s0,A^[I,L0],c3); CADD(c2,c3,A^[I,K0]);
        CMUL(c1,A^[I,L0],c2); CPRO(s0,u0,c3); CDIF(c2,c3,A^[I,L0])
      end;
      For K:=K0+1 to L0-1 do
      begin
        u0:=A^[K0,K];
        CMUL(c1,u0,c2); CONJ(A^[K,L0],u1); CPRO(s0,u1,c3);
        CADD(c2,c3,A^[K0,K]);
        CMUL(c1,A^[K,L0],c2); CONJ(u0,u1); CPRO(s0,u1,c3);
        CDIF(c2,c3,A^[K,L0])
      end;
      For J:=L0+1 to N do
      begin
        u0:=A^[K0,J];
        CMUL(c1,u0,c2); CPRO(s0,A^[L0,J],c3); CADD(c2,c3,A^[K0,J]);
        CMUL(c0,A^[L0,J],c2); CPRO(s0,u0,c3); CDIF(c2,c3,A^[L0,J])
      end;
      t0:=A^[K0,K0].R;
      t1:=4.0*Sqr(s0*CABS(A^[K0,L0]))/w0;
      A^[K0,K0].R:=Sqr(CABS(c0))*t0 + t1+Sqr(s0)*A^[L0,L0].R;
      A^[L0,L0].R:=Sqr(s0)*t0 - t1+Sqr(CABS(c0))*A^[L0,L0].R;
      A^[K0,L0]:=z0;
      For I:=1 to N do
      begin
        u0:=VX^[I,K0];
        CMUL(c0,u0,c2); CPRO(s0,VX^[I,L0],c3); CADD(c2,c3,VX^[I,K0]);
        CMUL(c1,VX^[I,L0],c2); CPRO(s0,u0,c3); CDIF(c2,c3,VX^[I,L0])
      end;
      t0:=0.0;
      For I:=1 to N-1 do
        For J:=I+1 to N do
          t0:=t0+Sqr(CABS(A^[I,J]));
      s:=2.0*t0;
      if s < dta then it:=1 else Inc(L)
    end {else}
  Until (L>M) or (it=1);
  If it=1 then
    For I:=1 to N do R^[I]:=A^[I,I].R
End; {of EPHJ}

    Procedure NORMAL(N:Integer; R:pVec; Z:pCMat);
{--------------------------------------------------------------
    Sort eigenvalues by absolute values in ascending order and
    normalize.
--------------------------------------------------------------}
    Label 5;
    Var
      Tr,Ti: pVec;
      ZM,Z1,Z2,VR: Double;
      V: Complex;
      IM,J,K:Integer;

      Procedure CDIV(AR,AI,BR,BI:Double; Var ZR,ZI:Double);
{     EFFECT DIVISION Z=(ZR+I*ZI):=(AR+I*AI)/(BR+I*BI)
 -----DO NOT USE IF BR:=BI=0.  }
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

    Begin

      New(Tr); New(Ti);

{     SORT SOLUTIONS IN ASCENDING ORDER }
      For J:=2 to N do
      begin
         VR:=R^[J];
         For K:=1 to N do
         begin
           Tr^[K]:=Z^[K,J].R;
           Ti^[K]:=Z^[K,J].I
         end;

         For I:=J-1 Downto 1 do
         begin
           IF ABS(R^[I]) <= ABS(VR) THEN GOTO 5;
           R^[I+1]:=R^[I];
           For K:=1 to N do
           begin
             Z^[K,I+1].R:=Z^[K,I].R;
             Z^[K,I+1].I:=Z^[K,I].I
           end
         end;
         I:=0;
5:       R^[I+1]:=VR;
         For K:=1 to N do
         begin
           Z^[K,I+1].R:=Tr^[K];
           Z^[K,I+1].I:=Ti^[K]
         end
      end;

{     NORMALIZE WITH RESPECT TO BIGGEST ELEMENT }
      For J:=N Downto 1 do
      begin
        ZM:= 0.0;
        For I:=1 to N do
        begin
          V.R:=Z^[I,J].R;
          V.I:=Z^[I,J].I;
          Z1:=CABS(V);
          IF Z1 >= ZM THEN
          begin
            IM := I;
            ZM := Z1
          end
        end;
        Z1 := Z^[IM,J].R;
        Z2 := Z^[IM,J].I;
        For I := 1 to N do
        begin
          CDIV(Z^[I,J].R,Z^[I,J].I,Z1,Z2,V.R,V.I);
          Z^[I,J].R := V.R;
          Z^[I,J].I := V.I
        end
      end;
      Dispose(Tr); Dispose(Ti)
    End; {NORMAL}


{Read from input file complex hermitian matrix}
Procedure LME;
Var I,J: Integer;
Begin
  Assign(fp,'matc4.dat'); Reset(fp);
  Readln(fp,N);
  For I:=1 to N do
  begin
    For J:=1 to N-1 do Read(fp, A^[I,J].R, A^[I,J].I);
    Read(fp,A^[I,N].R); Readln(fp,A^[I,N].I)
  end;
  Close(fp)
End;         


{main program}
BEGIN

  New(A); New(VX); New(R);

  LME;    {read complex matrix from text file}

  ClrScr;
  Writeln('  ** Compute Eigenvalues and Eigenvectors **');
  Writeln('         of a Square Hermitian Matrix');
  Writeln('             By Jacobi''s Method');
  Writeln;
  Write('  Precision: '); Readln(dta);
  Write('  Max. number of iterations: '); Readln(M);
  Writeln;

  EPHJ(dta, M, N, A, it, R, VX);

  NORMAL(N, R, VX);

  If it=-1 then
    Writeln('  No convergence.')
  else                             {display results}
    For J:=1 to N do
    begin
      Writeln('  Eigenvalue ',J,': ', R^[J]);
      Writeln('  Eigenvector:');
      For I:=1 to N do
      begin
        Write(' ',VX^[I,J].R);
        If VX^[I,J].I >= 0.0 then write(' + ')
                             else write(' - ');
        Writeln(Abs(VX^[I,J].I),' I')
      end;
      ReadKey
    end;

  Dispose(A); Dispose(VX); Dispose(R);
  DoneWinCrt

END.

{end of file Tephj.pas}