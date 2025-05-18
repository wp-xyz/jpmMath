{**************************************************************
!*      Purpose:  This program computes the modified Struve   *
!*                function Lv(x) for an arbitrary order v     * 
!*                using subroutine STVLV                      *
!*      Input :   v   --- Order of Lv(x)  ( |v| ó 20 )        *
!*                x   --- Argument of Lv(x) ( x ò 0 )         *
!*      Output:   SLV --- Lv(x)                               *
!*      Example:  x = 10.0                                    * 
!*                  v          Lv(x)                          *
!*                ------------------------                    *
!*                 0.5     .27785323D+04                      * 
!*                 1.5     .24996698D+04                      *
!*                 2.5     .20254774D+04                      *
!*                 3.5     .14816746D+04                      *
!*                 4.5     .98173460D+03                      *
!*                 5.5     .59154277D+03                      *
!*                ------------------------                    *
!*                                                            *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
!*                                (www.jpmoreau.fr)           *
!*************************************************************}  
Program MSTVLV;

Uses WinCrt;

Label 10;

Var   SLV, V, VMIN, VMAX, DV, X: double;


Function Power(y,x: double): double;
Begin
  IF x<0 THEN EXIT;
  Power:=Exp(x*Ln(y))
End;


Function GAMMA(X: double): double;

{       ==================================================
!       Purpose: Compute gamma function gamma(x)
!       Input :  x  --- Argument of gamma(x)
!                       ( x is not equal to 0,-1,-2)
!       Output:  gamma(x)
!       ================================================== }
Var
        GA, GR, R, Z: double;
        G: Array[1..26] of double;
        K, M, M1: integer;
Begin
        if X=int(X) then
	   if (X>0.0) then
           begin
              GA:=1.0;
              M1:=Round(X-1);
              for K:=2 to M1 do  GA:=GA*K
           end
           else
              GA:=1E+100
	 else
         begin
	    if abs(X)>1.0 then
            begin
              Z:=abs(X);
              M:=Round(Z);
              R:=1.0;
              for K:=1 to M do  R:=R*(Z-K);
              Z:=Z-M
            end
            else
             Z:=X;

           G[1]:=1.0; G[2]:=0.5772156649015329;
	       G[3]:=-0.6558780715202538; G[4]:=-0.420026350340952e-1;
	       G[5]:=0.1665386113822915; G[6]:=-0.421977345555443e-1;
	       G[7]:=-0.96219715278770e-2; G[8]:=0.72189432466630e-2;
	       G[9]:=-0.11651675918591e-2; G[10]:=-0.2152416741149e-3;
           G[11]:=0.1280502823882e-3; G[12]:=-0.201348547807e-4;
	       G[13]:=-0.12504934821e-5; G[14]:=0.11330272320e-5;
	       G[15]:=-0.2056338417e-6; G[16]:=0.61160950e-8;
	       G[17]:=0.50020075e-8; G[18]:=-0.11812746e-8;
	       G[19]:=0.1043427e-9; G[20]:=0.77823e-11;
	       G[21]:=-0.36968e-11; G[22]:=0.51e-12;
	       G[23]:=-0.206e-13; G[24]:=-0.54e-14;
	       G[25]:=0.14e-14; G[26]:=0.1e-15;

           GR:=G[26];
           for K:=25 Downto 1 do  GR:=GR*Z+G[K];
           GA:=1.0/(GR*Z);
           if abs(X)>1.0 then
           begin
              GA:=GA*R;
              if X < 0.0 then GA:=-PI/(X*GA*sin(PI*X))
           end
        end;
        GAMMA := GA
End;


Procedure STVLV(V, X: double; Var SLV: double);
{       ======================================================
!       Purpose:  Compute modified Struve function Lv(x) with
!                 an arbitrary order v
!       Input :   v   --- Order of Lv(x)  ( |v| ó 20 )
!                 x   --- Argument of Lv(x) ( x ò 0 )
!       Output:   SLV --- Lv(x)
!       Routine called: GAMMA to compute the gamma function
!       ====================================================== }
Label   15, 30, fin;
Var
        BIV, BIV0, GA, GB, R, R1, R2, S, S0, SA, U, U0, V0, VA, VB, VT: double;
	BF, BF0, BF1: double;
	K, L, N: Integer;
Begin
        if (X = 0.0) then
        begin
           if (V > -1.0) or (Round(V)-V = 0.5) then
              SLV:=0.0
           else if (V < -1.0) then
              SLV:=power(-1, int(0.5-V)-1)*1.0e+300
           else if (V = -1.0) then
              SLV:=2.0/PI;
	   goto fin
        end;
        if (X <= 40.0) then
        begin
           V0:=V+1.5;
           GA := GAMMA(V0);
           S:=2.0/(sqrt(PI)*GA);
           R1:=1.0;
           for K:=1 to 100 do
           begin
              VA:=K+1.5;
              GA := GAMMA(VA);
              VB:=V+K+1.5;
              GB := GAMMA(VB);
              R1:=R1*(0.5*X)*(0.5*X);
              R2:=R1/(GA*GB);
              S:=S+R2;
              if abs(R2/S) < 1.0e-12 then goto 15
	   end;
15:        SLV:=power(0.5*X,V+1.0)*S
        end
        else
        begin
           SA:=-1.0/PI*power(0.5*X, V-1.0);
           V0:=V+0.5;
           GA := GAMMA(V0);
           S:=-sqrt(PI)/GA;
           R1:=-1.0;
           for K:=1 to 12 do
           begin
              VA:=K+0.5;
              GA := GAMMA(VA);
              VB:=-K+V+0.5;
              GB := GAMMA(VB);
              R1:=-R1/power(0.5*X,2);
              S := S + R1*GA/GB
	   end;
           S0:=SA*S;
           U:=abs(V);
           N:=Round(U);
           U0:=U-N;
           for L:=0 to 1 do
           begin
              VT:=U0+L;
              R:=1.0;
              BIV:=1.0;
              for K:=1 to 16 do
              begin
                 R:=-0.125*R*(4.0*VT*VT-power(2.0*K-1.0,2))/(K*X);
                 BIV := BIV + R;
                 if abs(R/BIV) < 1.0e-12 then goto 30;
	      end;
30:           if (L = 0) then BIV0:=BIV
	   end;
           BF0:=BIV0;
           BF1:=BIV;
           for K:=2 to N do
           begin
              BF:=-2.0*(K-1.0+U0)/X*BF1+BF0;
              BF0:=BF1;
              BF1:=BF
           end;
           if (N = 0) then BIV:=BIV0;
           if (N > 1) then BIV:=BF;
           SLV:=exp(X)/sqrt(2.0*PI*X)*BIV+S0
        end;
fin: End;


{main program}
BEGIN

        writeln;
        write(' Please enter vmin, vmax, dv and x: ');
        readln(VMIN, VMAX, DV, X);

        writeln;
        writeln('   v           Lv(x)     ');
        writeln(' ------------------------');
        
        V:=VMIN;                                  

10:     STVLV(V,X,SLV);
        
        writeln(V:5:1, SLV:18:8);
	V := V + DV;
	if V<=VMAX then goto 10;

        writeln(' ------------------------');
        ReadKey;

        DoneWinCrt

END.

{end of file mstvlv.pas}