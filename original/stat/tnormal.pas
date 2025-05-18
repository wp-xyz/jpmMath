{********************************************************************
* This program calculates the standardized normal law probabilities *
* (mean=0, standard deviation=1) for -37.0 <= X <= 37.0             *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*   X=5.         Number of standard deviations                      *
*   P1           Probability to the left of X                       *
*   P2           Probability to the right of X                      *
*   Q            Probability Density for X                          *
* (Use successively function alnorm and subroutines NORMP and       *
*  NPROB).                                                          *
*                                                                   *
*   P1= 9.99999713348428E-0001                                      *
*   P2= 2.86651571867398E-0007                                      *
*                                                                   *
*   P1= 9.99999713348428E-0001                                      *
*   P2= 2.86651571892333E-0007                                      *
*   Q=  1.48671951473430E-0006                                      *
*                                                                   *
*   P1= 9?99999713348428E-0001                                      *
*   P2= 2.86651571867398E-0007                                      *
*   Q=  1.48671951467306E-0006                                      *
*                                                                   *
* ----------------------------------------------------------------- *
* Ref.: Journal of Applied Statistics (1973) vol22 no.3.            *
*                                                                   *
*                             Pascal Release By J-P Moreau, Paris.  *
*                                      (www.jpmoreau.fr)            *
********************************************************************}
Program TNormal;
Uses WinCrt;

Const
      zero: Double = 0.0;
      one: Double = 1.0;
      half: Double = 0.5;
      ltone: Double = 7.0;
      utzero: Double = 18.66;
      con: Double = 1.28;

Var
      P1, P2, Q, X: Double;
      UP: Boolean;


{ This file includes the Applied Statistics algorithm AS 66 for calculating
  the tail area under the normal curve, and two alternative routines which
  give higher accuracy.   The latter have been contributed by Alan Miller of
  CSIRO Division of Mathematics & Statistics, Clayton, Victoria.   Notice
  that each function or routine has different call arguments.             }

    Function alnorm(x:double; upper:boolean): Double;

{       Algorithm AS66 Applied Statistics (1973) vol22 no.3

        Evaluates the tail area of the standardised normal curve
        from x to infinity if upper is .true. or
        from minus infinity to x if upper is .false.              }

    Label 10,20,30,40;
    Var
        z,y,Temp: double;
        p,q,r,a1,a2,a3,b1,b2,c1,c2,c3,c4,c5,c6: double;
        d1,d2,d3,d4,d5: double;
        up: boolean;

    Begin
{*** machine dependent constants }
      p:=0.398942280444; q:=0.39990348504; r:=0.398942280385;   
      a1:=5.75885480458; a2:=2.62433121679; a3:=5.92885724438;  
      b1:=-29.8213557807; b2:=48.6959930692;
      c1:=-3.8052E-8; c2:=3.98064794E-4; c3:=-0.151679116635;
      c4:=4.8385912808; c5:=0.742380924027; c6:=3.99019417011;  
      d1:=1.00000615302; d2:=1.98615381364; d3:=5.29330324926;  
      d4:=-15.1508972451; d5:=30.789933034;

      up:=upper;
      z:=x;
      if z >= zero then goto 10;
      up:= not up;
      z:=-z;
10:   if ((z <= ltone) or (up=TRUE)) and (z <= utzero) then goto 20;
      Temp:=zero;
      goto 40;
20:   y:=half*z*z;
      if z > con then goto 30;
      Temp:=half-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
      goto 40;
30:   Temp:=r*Exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
40:   if up then alnorm:=Temp
            else alnorm:=one-Temp
    End;


    Procedure NORMP(Z:Double; Var P, Q, PDF:Double);

{	Normal distribution probabilities accurate to 1.e-15.
 	Z = no. of standard deviations from the mean.
	P, Q = probabilities to the left & right of Z.   P + Q = 1.
        PDF = the probability density.
        Based upon algorithm 5666 for the error function, from:
        Hart, J.F. et al, 'Computer Approximations', Wiley 1968

        Programmer: Alan Miller

        Latest revision - 30 March 1986       }

    Label Return;
    Var P0,P1,P2,P3,P4,P5,P6: Double;
        Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7: Double;
        CUTOFF,EXPNTL,ROOT2PI: Double;
        ZABS:Double;

    Begin

	P0:=220.2068679123761; P1:=221.2135961699311; P2:=112.0792914978709;
      	P3:= 33.91286607838300; P4:= 6.373962203531650;
      	P5:=  0.7003830644436881; P6:= 0.3526249659989109E-01;
        Q0:=440.4137358247522; Q1:=793.8265125199484; Q2:=637.3336333788311;
      	Q3:=296.5642487796737; Q4:= 86.78073220294608; Q5:=16.06417757920695;
        Q6:=  1.755667163182642; Q7:=0.8838834764831844E-01;
        CUTOFF:= 7.071; ROOT2PI:=2.506628274631001;

	ZABS := ABS(Z);

{	|Z| > 37.0 }

	IF ZABS > 37.0 THEN
        begin
	  PDF := zero;
	  IF Z >  zero THEN
          begin
	    P := one;
	    Q := zero
          end
	  ELSE
          begin
	    P := zero;
	    Q := one
	  end;
	  goto RETURN
	end;

{	|Z| <= 37.0 }

	EXPNTL := Exp(-half*ZABS*ZABS);
	PDF := EXPNTL/ROOT2PI;

{	|Z| < CUTOFF = 10/sqrt(2)  }

	IF ZABS <  CUTOFF THEN
	  P := EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS +
     	       P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS  +
     	       Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + Q0)

{	|Z| >= CUTOFF }

	ELSE
	  P := PDF/(ZABS + one/(ZABS + 2.0/(ZABS + 3.0/(ZABS + 4.0/
     		(ZABS + 0.65)))));

	IF Z < zero THEN
	  Q := one - P
	ELSE
        begin
	  Q := P;
	  P := one - Q
	end;
RETURN: End;


      Procedure NPROB(Z:Double; Var P,Q,PDF:Double);

{       P, Q = PROBABILITIES TO THE LEFT AND RIGHT OF Z
        FOR THE STANDARD NORMAL DISTRIBUTION.
        PDF  = THE PROBABILITY DENSITY FUNCTION

        REFERENCE: ADAMS,A.G. AREAS UNDER THE NORMAL CURVE,
        ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969.

        LATEST REVISION - 23 JANUARY 1981                 }

      Label 10,20,30, Return;

      Var A0,A1,A2,A3,A4,A5,A6,A7: Double;
          B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11: Double;
          Y, ZABS:Double;
      Begin

        A0:=half; A1:=0.398942280444; A2:=0.399903438504;
        A3:=5.75885480458; A4:=29.8213557808;
        A5:=2.62433121679; A6:=48.6959930692; A7:=5.92885724438;

        B0:=0.398942280385; B1:=3.8052E-8; B2:=1.00000615302;
        B3:=3.98064794E-4; B4:=1.98615381364; B5:=0.151679116635;
        B6:=5.29330324926; B7:=4.8385912808; B8:=15.1508972451;
        B9:=0.742380924027; B10:=30.789933034; B11:=3.99019417011;

        ZABS := ABS(Z);
        IF ZABS > 12.7 Then GOTO 20;
        Y := A0*Z*Z;
        PDF := Exp(-Y)*B0;
        IF ZABS > 1.28 Then GOTO 10;

{       Z BETWEEN -1.28 AND +1.28 } 

        Q := A0-ZABS*(A1-A2*Y/(Y+A3-A4/(Y+A5+A6/(Y+A7))));
        IF Z < zero Then GOTO 30;
        P := one-Q;
        GOTO RETURN;

{       ZABS BETWEEN 1.28 AND 12.7 }

10:     Q := PDF/(ZABS-B1+B2/(ZABS+B3+B4/(ZABS-B5+B6/(ZABS+B7-B8/
        (ZABS+B9+B10/(ZABS+B11))))));
        IF Z < zero Then GOTO 30;
        P := one-Q;
        GOTO RETURN;

{       Z FAR OUT IN TAIL }

20:     Q := zero;
        PDF := zero;
        IF Z < zero Then GOTO 30;
        P := one;
        goto RETURN;

{       NEGATIVE Z, INTERCHANGE P AND Q }

30:     P := Q;
        Q := one-P;

RETURN: End;


{main program}
BEGIN

  X:=5.0;

  UP:=False;
  P1:= alnorm(X,UP);
  UP:=True;
  P2:= alnorm(X,UP);

  writeln;
  writeln('  P1=', P1);
  writeln('  P2=', P2);

  NORMP(X, P1, P2, Q);
  writeln;
  writeln('  P1=', P1);
  writeln('  P2=', P2);
  writeln('  Q= ', Q);

  NPROB(X, P1, P2, Q);
  writeln;
  writeln('  P1=', P1);
  writeln('  P2=', P2);
  writeln('  Q= ', Q);
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file Tnormal.pas}