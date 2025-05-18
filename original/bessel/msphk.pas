{******************************************************************
*        Purpose: This program computes the modified spherical    *
*                Bessel functions kn(x) and kn'(x) using          *
*                 subroutine SPHK                                 *
*        Input :  x --- Argument of kn(x)  ( x > 0 )              *
*                 n --- Order of kn(x) ( n = 0 to 250 )           *
*        Output:  SK(n) --- kn(x)                                 *
*                 DK(n) --- kn'(x)                                *
*        Example: x= 10.0                                         *
*                   n          kn(x)               kn'(x)         *
*                 --------------------------------------------    *
*                   0     .7131404291D-05    -.7844544720D-05     * 
*                   1     .7844544720D-05    -.8700313235D-05     *
*                   2     .9484767707D-05    -.1068997503D-04     * 
*                   3     .1258692857D-04    -.1451953914D-04     *
*                   4     .1829561771D-04    -.2173473743D-04     *
*                   5     .2905298451D-04    -.3572740841D-04     * 
* --------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special         *
*             Functions,                                          *
*             jin.ece.uiuc.edu/routines/routines.html".           *
*                                                                 *
*                            Pascal Release By J-P Moreau, Paris. *
*                                     (www.jpmoreau.fr)           *
******************************************************************}
Program MSPHK;

Uses WinCrt;

Type VEC = Array[0..250] of double;

Var
     SK, DK: VEC;
     k, n, nm, ns: Integer;
     X: double;

Procedure SPHK(N:integer; X:double; Var NM:integer; Var SK, DK: VEC);
{       =====================================================
!       Purpose: Compute modified spherical Bessel functions
!                of the second kind, kn(x) and kn'(x)
!       Input :  x --- Argument of kn(x)  ( x ò 0 )
!                n --- Order of kn(x) ( n = 0,1,2,... )
!       Output:  SK(n) --- kn(x)
!                DK(n) --- kn'(x)
!                NM --- Highest order computed
!       ===================================================== }
Label 10, 20;
var K:integer; F, F0, F1:double;
Begin
        NM:=N;
        if X < 1E-60 then
        begin
           for K:=0 to N do
           begin
              SK[K]:=1E+100;
              DK[K]:=-1E+100;
           end;
           goto 20;
	end;
        SK[0]:=0.5*PI/X*Exp(-X);
        SK[1]:=SK[0]*(1.0+1.0/X);
        F0:=SK[0];
        F1:=SK[1];
	for K:=2 to N do
        begin
           F:=(2.0*K-1.0)*F1/X+F0;
           SK[K]:=F;
           if abs(F) > 1E+100 then goto 10;
           F0:=F1;
           F1:=F
        end; 
10:     NM:=K-1;
        DK[0]:=-SK[1];
        for K:=1 to NM do
          DK[K]:=-SK[K-1]-(K+1.0)/X*SK[K];
20: End;


{main program}
BEGIN
        InitWinCrt;
        writeln;
        write(' Please enter n and x: ');
        readln(n, X);

        n:=n+1;

	if n <= 10 then
           ns:=1
        else
        begin
           write(' Please enter order step Ns: ');
           readln(ns)
        end;
        
	SPHK(n,X,nm,SK,DK);
        
	writeln;
        writeln('  n          jn(x)                 jn''(x)');
        writeln('----------------------------------------------------');
        k:=0;
        Repeat 
          writeln(k:3,' ',SK[k],' ',DK[k]);
          k:=k+ns
        until k = nm+1;
        ReadKey;
        DoneWinCrt
END.

{end of file msphk.pas}