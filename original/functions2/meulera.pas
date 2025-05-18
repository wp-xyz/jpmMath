{********************************************************************
*   This program computes Euler number En using subroutine EULERA   *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*  Compute Euler number En for n = 0,1,...,10.                      *
*                                                                   *
*  Please enter Nmax: 10                                            *
*                                                                   *
*    n            En                                                *
*  --------------------------                                       *
*    0                  1.00                                        *
*    2                 -1.00                                        *
*    4                  5.00                                        *
*    6                -61.00                                        *
*    8               1385.00                                        *
*   10             -50521.00                                        *
*  --------------------------                                       *
*                                                                   *
* ----------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,*
*             jin.ece.uiuc.edu/routines/routines.html".             *
*                                                                   *
*                              Pascal Release By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)           *
********************************************************************}
  PROGRAM MEULERA;
  Uses WinCrt;

  Const NMAX = 200;

  Type
       pVEC = ^VEC;
        VEC = Array[0..NMAX] of Double;

  Var
       E: pVEC;
       K,N: Integer;


  Procedure EULERA(N:Integer; EN:pVEC);
{   ======================================
      Purpose: Compute Euler number En
      Input :  n --- Serial number
      Output:  EN(n) --- En
    ====================================== }
  Var J,K,M:Integer; R,S:Double;
  Begin
    EN^[0]:=1.0;
    For M:=1 to N Div 2 do
    begin
      S:=1.0;
      For K:=1 to M-1 do
      begin
        R:=1.0;
        For J:=1 to 2*K do
          R:=R*(2.0*M-2.0*K+J)/J;
        S:=S+R*EN^[2*K]
      end;
      EN^[2*M]:=-S
    end
  End;


  {main program}
  BEGIN
    New(E);
    writeln;
    write('  Please enter Nmax: '); Readln(N);

    EULERA(N,E);

    writeln;
    writeln('   n            En');
    writeln(' --------------------------');
    K:=0;
    Repeat
      writeln(K:4,'  ',E^[K]:20:2);
      Inc(K,2)
    Until K>N;
    writeln(' --------------------------');
    ReadKey;
    Dispose(E);
    DoneWinCrt

  END.

{end of file meulera.pas}