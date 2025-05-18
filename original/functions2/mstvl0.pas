{**************************************************************
!*      Purpose: This program computes modified Struve        *
!*               function L0(x) using subroutine STVL0        *
!*      Input :  x   --- Argument of L0(x) ( x ò 0 )          *
!*      Output:  SL0 --- L0(x)                                *
!*      Example:                                              *
!*                  x        L0(x)                            *
!*              ------------------------                      *
!*                 0.0   .00000000D+00                        *
!*                 5.0   .27105917D+02                        *
!*                10.0   .28156522D+04                        *
!*                15.0   .33964933D+06                        *
!*                20.0   .43558283D+08                        *
!*                30.0   .78167230D+12                        *
!*                40.0   .14894775D+17                        *
!*                50.0   .29325538D+21                        *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
!*                                 (www.jpmoreau.fr)          *
!*************************************************************}  
Program MSTVL0;

Uses WinCrt;

Var  X, SL0: Double;


Procedure STVL0(X: double; Var SL0: double);
{       ================================================
!       Purpose: Compute modified Struve function L0(x)
!       Input :  x   --- Argument of L0(x) ( x ò 0 )
!       Output:  SL0 --- L0(x)
!       ================================================ }
Label   15, 25, 35;
Var
        A0,A1,BI0,R,S: Double;
        K, KM: Integer;
Begin	    
        S:=1.0;
        R:=1.0;
        if X <= 20.0 then
        begin
           A0:=2.0*X/PI;
           for K:=1 to 60 do
           begin
              R:=R*Sqr(X/(2.0*K+1.0));
              S:=S+R;
              if abs(R/S) < 1.0E-12 then goto 15
           end;
15:        SL0:=A0*S
        end
        else
        begin
           KM:=Round(0.5*(X+1.0));
           if (X >= 50.0) then KM:=25;
           for K:=1 to KM do
           begin
              R:=R*Sqr((2.0*K-1.0)/X);
              S:=S+R;
              if (abs(R/S) < 1.0E-12) then goto 25
           end;
25:        A1:=exp(X)/sqrt(2.0*PI*X);
           R:=1.0;
           BI0:=1.0;
           for K:=1 to 16 do
           begin
              R:=0.125*R*Sqr(2.0*K-1.0)/(K*X);
              BI0:=BI0+R;
              if (abs(R/BI0) < 1.0e-12) then goto 35
           end;
35:        BI0:=A1*BI0;
           SL0:=-2.0/(PI*X)*S+BI0
        end
End;


BEGIN

        writeln;
	write(' Please enter x: ');
        readln(X);
        writeln;
        writeln('   x          L0(x)           ');
        writeln('------------------------------');
        
	STVL0(X, SL0);
       
	writeln(X:5:1,' ',SL0);
        readkey;
        DoneWinCrt

END.

{end of file mstvl0.pas}