{*************************************************************
*   Solving a complex homogeneous linear system by Gauss     *
*   Method with full pivoting.                               *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
*                        | -1  i -i |                        *
* Solve AX = 0 with  A = | -i -1  1 |                        *
*                        |  i  1 -1 |                        *
*                                                            *
* Basic Family of system solutions:                          *            
*  Solution 1                                                *
*   X1 =  0.00000000000000E+0000 +  1.00000000000000E+0000 I *
*   X2 =  1.00000000000000E+0000 +  0.00000000000000E+0000 I *
*   X3 =  0.00000000000000E+0000 +  0.00000000000000E+0000 I *
*                                                            *
*  Solution 2                                                *
*   X1 =  0.00000000000000E+0000 -  1.00000000000000E+0000 I *
*   X2 =  0.00000000000000E+0000 +  0.00000000000000E+0000 I *
*   X3 =  1.00000000000000E+0000 +  0.00000000000000E+0000 I *
*                                                            *
* ---------------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal        *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988".          *
*************************************************************}
Program TRSHCGT;

Uses WinCrt;

Const NMAX = 20;

  Type
       Complex = Record
         R, I: Real
       End;

       Matc = Array[1..NMAX,1..NMAX] of Complex;
       Vecc = Array[1..NMAX] of Complex;
       Veci = Array[1..NMAX] of Integer;

  Var
       i,l,N,M0,R0: Integer;
       eps: Real;

       A, VX: Matc;


  Function CABS(Z:Complex): Real;
  Begin
    CABS := sqrt(Z.R*Z.R+Z.I*Z.I)
  End;

  Procedure CADD(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R + Z2.R;
    Z.I := Z1.I + Z2.I
  End;

  Procedure CDIF(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R - Z2.R;
    Z.I := Z1.I - Z2.I
  End;

  Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R*Z2.R - Z1.I*Z2.I;
    Z.I := Z1.R*Z2.I + Z1.I*Z2.R
  End;

  Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
  Var d:Real; C:Complex;
  Begin
    d := Z2.R*Z2.R + Z2.I*Z2.I;
    if d<2E-12 then
      {writeln(' Complex Divide by zero !') }
    else
    begin
      C.R:=Z2.R; C.I:=-Z2.I;
      CMUL(Z1,C,Z);
      Z.R:=Z.R/d; Z.I:=Z.I/d
    end
  End;

  Procedure CSwap(Var Z1,Z2: Complex);
  Var tmp:Complex;
  Begin
    tmp.R:=Z2.R; tmp.I:=Z2.I;
    Z2.R:=Z1.R; Z2.I:=Z1.I;
    Z1.R:=tmp.R; Z1.I:=tmp.I
  End;

  {*********************************************************
  * Procedure RSHCGT solves the complex homogeneous linear *
  * system AX = 0 by Gauss method with full pivoting.      *
  * ------------------------------------------------------ *
  * INPUTS:                                                *
  *         eps:  required precision                       *
  *          N :  size of complex matrix A (NxN)           *
  * OUTPUTS:                                               *
  *          R0:  rank of A in output                      *
  *          M0:  dimension of system solutions space      *
  *               (if R0 <> N)                             *
  *          VX:  contains in output the M0 solution       *
  *               vectors (stored in columns 1..M0) when   *
  *               R0<>N, or the unique solution vector     *
  *               (stored in column 1) when R0=N.          * 
  *********************************************************}   
  Procedure RSHCGT(eps:real; N:integer; A:Matc; Var R0,M0:integer;
                     Var VX:Matc);
  Var i,j,k:integer;
      k0,l,l0,n0:integer;
      C0,C1:Complex;
      P0,Q0,S:Complex;
      Z0,Z1:Complex;
      I0:Array[1..NMAX] of integer;
  Begin
    Z0.R:=0.0; Z0.I:=0.0;
    Z1.R:=1.0; Z1.I:=0.0;
    For K:=1 to N do I0[k]:=K;
    R0:=0; K:=1;
    Repeat
      P0:=A[k,k]; l0:=k; k0:=k;
      For i:=k to N do
        For j:=k to N do
          if CABS(A[i,j]) > CABS(P0) then
          begin
            P0:=A[i,j]; l0:=i; k0:=j
          end;
      if CABS(P0)<eps then
      begin
        R0:=k-1;
        M0:=N-R0
      end
      else
        if k=N then
        begin
          R0:=N; M0:=0
        end
        else
        begin
          if l0<>k then
            For j:=k to N do CSwap(A[k,j],A[l0,j]);
          if k0<>k then
            For i:=1 to N do CSwap(A[i,k],A[i,k0]);
          N0:=I0[k]; I0[k]:=I0[k0]; I0[k0]:=N0
        end;
      For i:=k+1 to N do
      begin
        CDIV(A[i,k],P0,Q0);
        For j:=k+1 to N do
        begin
          C0:=A[i,j];
          CMUL(Q0,A[k,j],C1);
          CDIF(C0,C1,A[i,j])
        end
      end;
      Inc(k)
    Until (R0<>0) or (k>N);

    if R0=N then
      For i:=1 to N do VX[i,1]:=Z0
    else
    begin
      For l:=1 to M0 do
      begin
        For i:=N Downto R0+1 do
          if i=R0+l then VX[i,l]:=Z1
                    else VX[i,l]:=Z0;
        For i:=R0 Downto 1 do
        begin
          S:=Z0;
          For j:=i+1 to N do
          begin
            C0:=S;
            CMUL(A[i,j],VX[j,l],C1);
            CADD(C0,C1,S)
          end;
          CDIV(S,A[i,i],C0);
          VX[i,l].R:=-C0.R;
          VX[i,l].I:=-C0.I
        end
      end;
      For i:=1 to N-1 do
      begin
        j:=i+1;
        While j<=N do
          if I0[j]=i then
          begin
            For l:=1 to M0 do CSwap(VX[j,l],VX[i,l]);
            I0[j]:=I0[i]; I0[i]:=i;
            j:=N+1
          end
          else
            Inc(j)
      end
    end

  End; {RSHCGT}                          
  

{main program}
BEGIN

  N:=3;   {size of complex matrix A}

  A[1,1].R:=-1.0; A[1,1].I:= 0.0;
  A[1,2].R:= 0.0; A[1,2].I:= 1.0;
  A[1,3].R:= 0.0; A[1,3].I:=-1.0;

  A[2,1].R:= 0.0; A[2,1].I:=-1.0;
  A[2,2].R:=-1.0; A[2,2].I:= 0.0;
  A[2,3].R:= 1.0; A[2,3].I:= 0.0;

  A[3,1].R:= 0.0; A[3,1].I:= 1.0;
  A[3,2].R:= 1.0; A[3,2].I:= 0.0;
  A[3,3].R:=-1.0; A[3,3].I:= 0.0;

  eps:=1E-10;

  RSHCGT(eps,N,A,R0,M0,VX);

  Writeln;
  if R0=N then
  begin
    writeln(' System solution:');
    For i:=1 to N do writeln('  X',I,' = ',VX[i,1].R)
  end
  else
  begin
    writeln(' Basic Family of system solutions:');
    For l:= 1 to M0 do
    begin
      writeln('  Solution ', l);
      For i:=1 to N do
      begin
        write('   X',I,' = ', VX[i,l].R);
        if VX[i,l].I >= 0.0 then write(' + ')
                            else write(' - ');
        writeln(ABS(VX[i,l].I),' I')
      end
    end
  end;       

  ReadKey;
  DoneWinCrt

END.

{end of file rshcgt.pas}