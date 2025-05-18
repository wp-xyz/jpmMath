{**************************************************************
!*      Purpose:  This program computes Struve function Hv(x) *
!*                for an arbitrary order using subroutine     *
!*                STVHV                                       * 
!*      Input :   v  --- Order of Hv(x)  ( -8.0 to 12.5 )     *
!*                x  --- Argument of Hv(x) ( x > 0 )          *
!*      Output:   HV --- Hv(x)                                *
!*      Example:  x = 10.0                                    *
!*                  v           Hv(x)                         *
!*                -----------------------                     *
!*                  .5     .46402212D+00                      * 
!*                 1.5     .14452322D+01                      *                      
!*                 2.5     .31234632D+01                      *
!*                 3.5     .53730255D+01                      *
!*                 4.5     .72083122D+01                      *
!*                 5.5     .76851132D+01                      *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
!*                                (www.jpmoreau.fr)           *
!*************************************************************}
Program MSTVHV;

Uses WinCrt;

Label 10;

Var   HV, V, VMIN, VMAX, DV, X: Double;

{y power x}
Function XPower(y,x:double): double;
BEGIN
  IF x<0 THEN EXIT;
  XPower:=Exp(x*Ln(y))
END;

{x power n}
Function Power(x:double; n:integer): double;
var i,m : integer;
    result :double;
begin
  result := 1.0;
  if n=0 then
  begin
    Power:=result;
    exit;
  end;
  m:=  n;
  if n<0 then m:=-n;
  for i:=1 to m do result :=x*result;
  Power :=result;
  if n<0 then Power:=1.0/result;
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


Procedure STVHV(V: double; X: double; Var HV: double);
{       =====================================================
!       Purpose: Compute Struve function Hv(x) with an
!                arbitrary order v
!       Input :  v  --- Order of Hv(x)  ( -8.0 ó v ó 12.5 )
!                x  --- Argument of Hv(x) ( x ò 0 )
!       Output:  HV --- Hv(x)
!       Routine called: GAMMA to compute the gamma function
!       ===================================================== }
Label   15, fin;
Var     GA,GB,PU0,PU1,QU0,QU1,R1,R2,S,S0,SA,T0,T1,U,U0,V0,VA,VB,VT: double;
	BF,BF0,BF1,BY0,BY1,BYV,SR: double;
        K, L, N: integer;
Begin
        if X = 0.0 then
        begin
           if (V > -1.0) or (Round(V)-V = 0.5) then
              HV:=0.0
           else if (V < -1.0) then
	      HV:=power(-1, Round(0.5-V)-1)*1.0E+300
           else if (V = -1.0) then
              HV:=2.0/PI;
           goto fin
        end;
        if (X <= 20.0) then
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
              R1:=-R1*power(0.5*X,2);
              R2:=R1/(GA*GB);
              S:=S+R2;
              if (abs(R2) < abs(S)*1.0e-12) then goto 15;
	       end;
15:        HV:=Xpower(0.5*X, V+1.0)*S;
        end
        else
        begin
           SA:=Xpower(0.5*X, V-1.0)/PI;
           V0:=V+0.5;
           GA := GAMMA(V0);
           S:=sqrt(PI)/GA;
           R1:=1.0;
           for K:=1 to 12 do
           begin
	          VA:=K+0.5;
              GA := GAMMA(VA);
              VB:=-K+V+0.5;
              GB := GAMMA(VB);
              R1:=R1/power(0.5*X,2);
              S:=S+R1*GA/GB;
           end;
           S0:=SA*S;
           U:=abs(V);
           N:=Round(U);
           U0:=U-N;
           for L:=0 to 1 do
           begin
              VT:=4.0*power(U0+L,2);
              R1:=1.0;
              PU1:=1.0;
              for K:=1 to 12 do
              begin
                 R1:=-0.0078125*R1*(VT-power(4.0*K-3.0,2))*(VT-power(4.0*K-1.0,2))/((2.0*K-1.0)*K*X*X);
                 PU1:=PU1+R1
	          end;
              QU1:=1.0;
              R2:=1.0;
              for K:=1 to 12 do
              begin
                 R2:=-0.0078125*R2*(VT-power(4.0*K-1.0,2))*(VT-power(4.0*K+1.0,2))/((2.0*K+1.0)*K*X*X);
                 QU1:=QU1+R2
	          end;
              QU1:=0.125*(VT-1.0)/X*QU1;
              if (L = 0) then
              begin
                 PU0:=PU1;
                 QU0:=QU1
              end
	       end;
           T0:=X-(0.5*U0+0.25)*PI;
           T1:=X-(0.5*U0+0.75)*PI;
           SR:=sqrt(2.0/(PI*X));
           BY0:=SR*(PU0*sin(T0)+QU0*cos(T0));
           BY1:=SR*(PU1*sin(T1)+QU1*cos(T1));
           BF0:=BY0;
           BF1:=BY1;
           for K:=2 to N do
           begin
              BF:=2.0*(K-1.0+U0)/X*BF1-BF0;
              BF0:=BF1;
              BF1:=BF
           end;
           if (N = 0) then BYV:=BY0;
           if (N = 1) then BYV:=BY1;
           if (N > 1) then BYV:=BF;
           HV:=BYV+S0
        end;
fin: End;


{main program}       
BEGIN

        writeln;
        write(' Please enter vmin, vmax, dv and x: ');
        readln(VMIN,VMAX,DV,X);

        writeln;
        writeln('    v         Hv(x)   ');
        writeln(' ------------------------------');
        
        V:=VMIN;                                  

10:     STVHV(V,X,HV);
        
        writeln(V:5:1,'  ',HV);
 	V := V + DV;
	if V<=VMAX then goto 10;

        writeln(' ------------------------------');

        ReadKey;
        DoneWinCrt

END.

{end of file mstvhv.pas}