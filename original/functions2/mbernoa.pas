{********************************************************************
* This program computes Bernoulli number Bn using subroutine BERNOA *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*  Compute Bernoulli number Bn for n = 0,1,...,10.                  *
*                                                                   *
*  Please enter Nmax: 10                                            *
*                                                                   *
*    n            Bn                                                *
*  --------------------------                                       *
*    0        1.000000000000                                        *
*    1       -0.500000000000                                        *
*    2        0.166666666667                                        *
*    4       -0.033333333333                                        *
*    6        0.023809523810                                        *
*    8       -0.033333333333                                        *
*   10        0.075757575758                                        *
*  --------------------------                                       *
*                                                                   *
* ----------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,*
*             jin.ece.uiuc.edu/routines/routines.html".             *
*                                                                   *
*                              Pascal Release By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)           *
********************************************************************}
  PROGRAM MBERNOA;
  USES WINCRT;

  Const NMAX = 100;

  Type
       pVEC = ^VEC;
        VEC = Array[0..NMAX] of Double;

  Var
      B: pVEC;
      K,N: Integer;


  Procedure BERNOA(N:Integer; BN:pVEC);
{ ======================================
   Purpose: Compute Bernoulli number Bn
   Input :  n --- Serial number
   Output:  BN(n) --- Bn
  ====================================== }
  Var J,K,M: Integer;
      R,S: Double;
  Begin
    BN^[0]:=1.0;
    BN^[1]:=-0.5;
    For M:=2 to N do
    begin
      S:=-(1.0/(M+1.0)-0.5);
      For K:=2 to M-1 do
      begin
        R:=1.0;
        For J:=2 to K do R:=R*(J+M-K)/J;
        S:=S-R*BN^[K]
      end;
      BN^[M]:=S
    end; 
    M:=3;
    Repeat
      BN^[M]:=0.0;
      Inc(M,2);
    Until M > N;
  End;


  {main program}
  BEGIN
    New(B);
    writeln;
    write('  Please enter Nmax: ');
    readln(N);

    BERNOA(N,B);

    writeln;
    writeln('   n            Bn');
    writeln(' --------------------------');
    writeln('   0  ',B^[0]:20:12);
    writeln('   1  ',B^[1]:20:12);
    K:=2;
    Repeat
      writeln(K:4,'  ',B^[K]:20:12);
      Inc(K,2);
    Until K > N;
    writeln(' --------------------------');
    writeln;
    ReadKey;
    Dispose(B);
    DoneWinCrt
  END.

{end of file mbernoa.pas}