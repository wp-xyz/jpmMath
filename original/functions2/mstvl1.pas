{*************************************************************
*      Purpose: This program computes the modified Struve    *
*               function L1(x) using subroutine STVL1        *
*      Input :  x   --- Argument of L1(x) ( x >= 0 )         *
*      Output:  SL1 --- L1(x)                                *
*      Example:                                              *
*                    x        L1(x)                          * 
*                -----------------------                     *
*                   0.0   .00000000D+00                      *
*                   5.0   .23728216D+02                      *
*                  10.0   .26703583D+04                      *
*                  15.0   .32812429D+06                      *
*                  20.0   .42454973D+08                      *
*                  30.0   .76853204D+12                      *
*                  40.0   .14707396D+17                      *
*                  50.0   .29030786D+21                      *
* ---------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special    *
*             Functions, jin.ece.uiuc.edu/routines/routines  *
*             .html".                                        *
*                                                            *
*                       Pascal Release By J-P Moreau, Paris. *
*                                 (www.jpmoreau.fr)          *
*************************************************************} 
Program MSTVL1;

Uses WinCrt;

Var X, SL1: Double;


Procedure STVL1(X: double; Var SL1: double);
{ ================================================
!  Purpose: Compute modified Struve function L1(x)
!  Input :  x   --- Argument of L1(x) ( x ò 0 )
!  Output:  SL1 --- L1(x)
!  ================================================ }
Label   15, 25, 35;
Var     A1,BI1,R,S: Double;
	K,KM: Integer;
Begin
        R:=1.0;
        if (X <= 20.0) then
        begin
           S:=0.0;
           for K:=1 to 60 do
           begin
              R:=R*X*X/(4.0*K*K-1.0);
              S:=S+R;
              if (abs(R/S) < 1.0E-12)  then goto 15
           end;
15:        SL1:=2.0/PI*S
        end
        else
        begin
           S:=1.0;
           KM:=Round(0.50*X);
           if (X > 50.0) then KM:=25;
           for K:=1 to KM do
           begin
              R:=R*(2.0*K+3.0)*(2.0*K+1.0)/(X*X);
              S:=S+R;
              if (abs(R/S) < 1.0E-12) then goto 25
           end;
25:        SL1:=2.0/PI*(-1.0+1.0/(X*X)+3.0*S/(X*X*X*X));
           A1:=exp(X)/sqrt(2.0*PI*X);
           R:=1.0;
           BI1:=1.0;
           for K:=1 to 16 do
           begin
              R:=-0.125*R*(4.0-(2.0*K-1.0)*(2.0*K-1.0))/(K*X);
              BI1:=BI1+R;
              if (abs(R/BI1) < 1.0E-12) then goto 35
           end;
35:        SL1:=SL1+A1*BI1
        end
End;

{main program}
BEGIN

        writeln;  
        write(' Please enter x: ');
        readln(X);
        writeln;
        writeln('   x          L1(x)\n');
        writeln('------------------------------');
        
	STVL1(X, SL1);
        
	writeln(X:5:1,' ',SL1);

        ReadKey;
        DoneWinCrt

END.

{end of file mstvl1.pas}