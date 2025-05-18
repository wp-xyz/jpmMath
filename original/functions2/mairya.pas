{**********************************************************************
!*       Purpose: This program computes Airy functions and their      *
!*                derivatives using subroutine AIRYA                  *
!*       Input:   xstart, xend --- Arguments of Airy function         * 
!*                xstep --- x increment                               *
!*       Output:  AI --- Ai(x)                                        * 
!*                BI --- Bi(x)                                        *
!*                AD --- Ai'(x)                                       *
!*                BD --- Bi'(x)                                       *
!*       Example:                                                     *
!*                xstart =0 xend = 30 xstep = 10                      *
!*                                                                    *
!*   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)      *
!*  ----------------------------------------------------------------  *
!*   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00   *
!*  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10   *
!*  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26   *
!*  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48   *
!*                                                                    *
!*   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)     *
!*  ----------------------------------------------------------------  *
!*   0       .35502805      .61492663     -.25881940      .44828836   *
!*  10       .04024124     -.31467983      .99626504      .11941411   *
!*  20      -.17640613     -.20013931      .89286286     -.79142903   *
!*  30      -.08796819     -.22444694     1.22862060     -.48369473   *
!* ------------------------------------------------------------------ *
!* REFERENCE: "Fortran Routines for Computation of Special Functions, *
!*             jin.ece.uiuc.edu/routines/routines.html".              *
!*                                                                    *
!*                               Pascal Release By J-P Moreau, Paris. *
!*                                        (www.jpmoreau.fr)           *
!*********************************************************************}        
Program Mairya;

Uses WinCrt1;

Label 5, 10;

Var
    AI,BI,AD,BD,X,XSTART,XEND,XSTEP: Double;

{calculate y power x}
Function Power(y,x:Double): Double;
Begin
  IF y<=0 THEN
    Power:=0.0
  else
    Power:=Exp(x*Ln(y))
End;

Procedure AJYIK(X: Double; Var VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2: Double);
{       =======================================================
!       Purpose: Compute Bessel functions Jv(x) and Yv(x),
!                and modified Bessel functions Iv(x) and
!                Kv(x), and their derivatives with v=1/3,2/3
!       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
!                      Kv(x) ( x ò 0 )
!       Output:  VJ1 --- J1/3(x)
!                VJ2 --- J2/3(x)
!                VY1 --- Y1/3(x)
!                VY2 --- Y2/3(x)
!                VI1 --- I1/3(x)
!                VI2 --- I2/3(x)
!                VK1 --- K1/3(x)
!                VK2 --- K2/3(x)
!       ======================================================= }
Label   20, 50, 60, 90, fin;
Var
	A0,B0,CK,GN1,GN2,GP1,GP2,QX,PX,R,RP,RP2,RQ,SK,UJ1,UJ2,UU0,VV,VJL,VL,VV0,X2,XK: Double;
	C0,GN,PV1,PV2,SUM,VIL,VSL: Double;
	K, K0, L: Integer;
Begin
        if X = 0.0 then
        begin
           VJ1:=0.0;
           VJ2:=0.0;
           VY1:=-1.0E+100;
           VY2:=1.0E+100;
           VI1:=0.0;
           VI2:=0.0;
           VK1:=-1.0E+100;
           VK2:=-1.0E+100;
           goto fin
        end;
        RP2:=0.63661977236758;
        GP1:=0.892979511569249;
        GP2:=0.902745292950934;
        GN1:=1.3541179394264;
        GN2:=2.678938534707747;
        VV0:=0.444444444444444;
        UU0:=1.1547005383793;
        X2:=X*X;
        K0:=12;
        if (X >= 35.0) then K0:=10;
        if (X >= 50.0) then K0:=8;
        if (X <= 12.0) then
	  for L:=1 to 2 do
          begin
              VL:=L/3.0;
              VJL:=1.0;
              R:=1.0;
              for K:=1 to 40 do
              begin
                 R:=-0.25*R*X2/(K*(K+VL));
                 VJL:=VJL+R;
                 if (abs(R) < 1.0E-15) then goto 20;
	      end;
20:           A0:=power(0.5*X,VL);
              if (L = 1) then VJ1:=A0/GP1*VJL;
              if (L = 2) then VJ2:=A0/GP2*VJL
          end
        else
            for L:=1 to 2 do
            begin
              VV:=VV0*L*L;
              PX:=1.0;
              RP:=1.0;
              for K:=1 to K0 do
              begin
                 RP:=-0.78125E-2*RP*(VV-power(4.0*K-3.0,2.0))*(VV-power(4.0*K-1.0,2.0))/(K*(2.0*K-1.0)*X2);
                 PX:=PX+RP;
              end; 
              QX:=1.0;
              RQ:=1.0;
              for K:=1 to K0 do
              begin
                 RQ:=-0.78125E-2*RQ*(VV-power(4.0*K-1.0,2.0))*(VV-power(4.0*K+1.0,2.0))/(K*(2.0*K+1.0)*X2);
                 QX:=QX+RQ;
              end;

              QX:=0.125*(VV-1.0)*QX/X;
              XK:=X-(0.5*L/3.0+0.25)*PI;
              A0:=sqrt(RP2/X);
              CK:=cos(XK);
              SK:=sin(XK);

              if (L = 1) then
              begin
                 VJ1:=A0*(PX*CK-QX*SK);
                 VY1:=A0*(PX*SK+QX*CK)
              end
              else if (L = 2) then
              begin
                 VJ2:=A0*(PX*CK-QX*SK);
                 VY2:=A0*(PX*SK+QX*CK)
              end
            end;
		   
        if (X <= 12.0) then
        begin
           for L:=1 to 2 do
           begin
              VL:=L/3.0;
              VJL:=1.0;
              R:=1.0;
	      for K:=1 to 40 do
              begin
                 R:=-0.25*R*X2/(K*(K-VL));
                 VJL:=VJL+R;
                 if (abs(R) < 1.0E-15) then goto 50;
              end;
50:           B0:=power(2.0/X,VL);
              if (L = 1) then UJ1:=B0*VJL/GN1;
              if (L = 2) then UJ2:=B0*VJL/GN2
           end;
           PV1:=PI/3.0;
           PV2:=PI/1.5;
           VY1:=UU0*(VJ1*cos(PV1)-UJ1);
           VY2:=UU0*(VJ2*cos(PV2)-UJ2)
        end;
        if (X <= 18.0) then 
	   for L:=1 to 2 do
           begin
              VL:=L/3.0;
              VIL:=1.0;
              R:=1.0;
	      for K:=1 to 40 do
              begin
                 R:=0.25*R*X2/(K*(K+VL));
                 VIL := VIL + R;
                 if (abs(R) < 1.0E-15) then goto 60
              end;
60:           A0:=power(0.5*X, VL);
              if (L = 1)  then VI1:=A0/GP1*VIL;
              if (L = 2)  then VI2:=A0/GP2*VIL
           end
	 else
         begin
           C0:=exp(X)/sqrt(2.0*PI*X);
           for L:=1 to 2 do
           begin
              VV:=VV0*L*L;
              VSL:=1.0;
              R:=1.0;
              for K:=1 to K0 do
              begin
                 R:=-0.125*R*(VV-power(2.0*K-1.0,2.0))/(K*X);
                 VSL := VSL + R
              end;
              if (L = 1) then VI1:=C0*VSL;
              if (L = 2) then VI2:=C0*VSL
           end
        end;
        if (X <= 9.0) then
           for L:=1 to 2 do
           begin
              VL:=L/3.0;
              if (L = 1) then GN:=GN1;
              if (L = 2) then GN:=GN2;
              A0:=power(2.0/X, VL)/GN;
              SUM:=1.0;
              R:=1.0;
              for K:=1 to 60 do
              begin
                 R:=0.25*R*X2/(K*(K-VL));
                 SUM := SUM + R;
                 if (abs(R) < 1.0E-15) then goto 90
              end;
90:           if (L = 1) then VK1:=0.5*UU0*PI*(SUM*A0-(VI1));
              if (L = 2) then VK2:=0.5*UU0*PI*(SUM*A0-(VI2))
           end
	 else
         begin
           C0:=exp(-X)*sqrt(0.5*PI/X);
           for L:=1 to 2 do
           begin
              VV:=VV0*L*L;
              SUM:=1.0;
              R:=1.0;
              for K:= 1 to K0 do
              begin
                 R:=0.125*R*(VV-power(2.0*K-1.0,2.0))/(K*X);
                 SUM := SUM + R
              end; 
              if (L = 1) then VK1:=C0*SUM;
              if (L = 2) then VK2:=C0*SUM
           end
        end;
fin: End;

Procedure AIRYA(X: Double; Var AI, BI, AD, BD: Double);

{       ======================================================
!       Purpose: Compute Airy functions and their derivatives
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       Routine called:
!                AJYIK for computing Jv(x), Yv(x), Iv(x) and
!                Kv(x) with v=1/3 and 2/3
!       ====================================================== }
Var
        C1,C2,PIR,SR3,VI1,VI2,VJ1,VJ2,VK1,VK2,VY1,VY2,XA,XQ,Z: Double;
Begin
        XA:=abs(X);
        PIR:=0.318309886183891;
        C1:=0.355028053887817;
        C2:=0.258819403792807;
        SR3:=1.732050807568877;

        Z:=Power(XA,1.5)/1.5;
        XQ:=sqrt(XA);

        AJYIK(Z,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2);

        if (X = 0.0) then
        begin
           AI:=C1;
           BI:=SR3*C1;
           AD:=-C2;
           BD:=SR3*C2
        end 
        else if (X > 0.0) then
        begin
           AI:=PIR*XQ/SR3*VK1;
           BI:=XQ*(PIR*VK1+2.0/SR3*VI1);
           AD:=-XA/SR3*PIR*VK2;
           BD:=XA*(PIR*VK2+2.0/SR3*VI2)
        end
        else
        begin
           AI:=0.5*XQ*(VJ1-VY1/SR3);
           BI:=-0.5*XQ*(VJ1/SR3+VY1);
           AD:=0.5*XA*(VJ2+VY2/SR3);
           BD:=0.5*XA*(VJ2/SR3-VY2)
        end
End;


{main program}
BEGIN
        InitWinCrt;
        writeln;
        write(' Please enter x start, x end, x step: ');
        readln(XSTART, XEND, XSTEP);

	X := XSTART;

5:      AIRYA(X,AI,BI,AD,BD);

	if (X = XSTART) then
        begin
           writeln;
	   writeln('    x        Ai(x)                    Bi(x)                  Ai''(x)      ');
           writeln('    Bi''(x)');      
	   writeln('  -------------------------------------------------------------------------')
	end;       
        writeln(X:3:0,' ',AI,' ',BI,' ',AD);
        writeln(' ',BD);
	X := X + XSTEP;
	if (X <= XEND) then goto 5;

	X := XSTART;

10:     AIRYA(-X,AI,BI,AD,BD);

        if (X = XSTART) then
        begin
           writeln;
	   writeln('    x       Ai(-x)                   Bi(-x)                 Ai''(-x)      ');
           writeln('    Bi''(x)');
           writeln('  -------------------------------------------------------------------------')
        end;
        writeln(X:3:0,' ',AI,' ',BI,' ',AD);
        writeln(' ',BD);
	X := X + XSTEP;
	if (X <= XEND) then goto 10;
        writeln;
        Readkey;
        DoneWinCrt
END.

{end of file mairya.pas}