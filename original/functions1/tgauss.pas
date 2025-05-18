{********************************************************
* Integration by Gauss method of a real function F=F(X) *
* or F=F(X,Y) or F=F(X,Y,Z). The integral is calculated *
* by using from 2 to 10 Gauss points.                   *
* ----------------------------------------------------- *
* Ref.: "Mécanique des vibrations linéaires By          *
*        M. Lalanne, P. Berthier and J. Der Hagopian,   *
*        Masson, Paris, 1980" [BIBLI 16].               *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* (Integrate F=sin(x) from x=0 to x=1).                 *
*                                                       *
* INTEGRATION OF A REAL FUNCTION BY GAUSS               *
*        F(X), F(X,Y) or F(X,Y,Z)                       *
*                                                       *
* Number of variables (1 to 3): 1                       *
*                                                       *
* How many Gauss points (2 to 10): 4                    *
*                                                       *
* Minimum value of X: 0                                 *
* Maximum value of X: 1                                 *
*                                                       *
* Value of integral =  4.59697693864200E-0001           *
*                                                       *
*                    TPW Version By J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)           *
********************************************************}
Program Test_Gauss;
Uses WinCrt;

Label   10,20,30,40;

Var
        A, H: Array[1..10] of REAL;
        x,x1,x2,y,y1,y2,z,z1,z2: REAL;
        i,j,k,n,n1,n2: INTEGER;
        a1,a2,a3,b1,b2,b3,xi9: REAL;

{define here function F}
Function F(x,y,z:REAL): REAL;
begin
  F := SIN(x)
end;

BEGIN
  writeln;
  writeln(' INTEGRATION OF A REAL FUNCTION BY GAUSS');
  writeln('        F(X), F(X,Y) or F(X,Y,Z)');
  writeln;
  write(' Number of variables (1 to 3): '); readln(n2);
  writeln;
  write(' How many Gauss points (2 to 10): '); readln(n);
  writeln;
  n1 := n - 1;
  Case n1 of
    1:begin
        A[1] := -0.57735026919;
        H[1] := 1.0
      end;
    2:begin
        A[1] := -0.774596669241;
        A[2] := 0.0;
        H[1] := 0.555555555556;
        H[2] := 0.888888888889
      end;
    3:begin
        A[1] := -0.8611363115939999;
        A[2] := -0.339981043585;
        H[1] := 0.347854845137;
        H[2] := 0.652145154863
      end;
    4:begin
        A[1] := -0.906179845939;
        A[2] := -0.538469310106;
        A[3] := 1e-12;
        H[1] := 0.236926885056;
        H[2] := 0.478628670499;
        H[3] := 0.568888888889
      end;
    5:begin
        A[1] := -0.9324695142029999;
        A[2] := -0.661209386466;
        A[3] := -0.238619186083;
        H[1] := 0.171324492379;
        H[2] := 0.360761573048;
        H[3] := 0.467913934573
      end;   
    6:begin
        A[1] := -0.949107912343;
        A[2] := -0.741531185599;
        A[3] := -0.405845151377;
        A[4] := 0.0;
        H[1] := 0.129484966169;
        H[2] := 0.279705391489;
        H[3] := 0.381830050505;
        H[4] := 0.417959183673
      end;
    7:begin
        A[1] := -0.949107912343;
        A[2] := -0.741531185599;
        A[3] := -0.405845151377;
        A[4] := 0.0;
        H[1] := 0.129484966169;
        H[2] := 0.279705391489;
        H[3] := 0.381830050505;
        H[4] := 0.417959183673
      end;
    8:begin
        A[1] := -0.960289856497;
        A[2] := -0.796666477414;
        A[3] := -0.525532409916;
        A[4] := -0.183434642496;
        H[1] := 0.10122853629;
        H[2] := 0.222381034453;
        H[3] := 0.313706645878;
        H[4] := 0.362683783378
      end;
    9:begin
        A[1] := -0.968160239508;
        A[2] := -0.836031107327;
        A[3] := -0.6133714327000001;
        A[4] := -0.324253423404;
        A[5] := 0.0;
        H[1] := 0.0812743883616;
        H[2] := 0.180648160695;
        H[3] := 0.260610696403;
        H[4] := 0.31234707704;
        H[5] := 0.330239355001
      end;
    10:begin
        A[1] := -0.973906528517;
        A[2] := -0.865063366689;
        A[3] := -0.679409568299;
        A[4] := -0.433395394129;
        A[5] := -0.148874338982;
        H[1] := 0.0666713443087;
        H[2] := 0.149451349151;
        H[3] := 0.219086362516;
        H[4] := 0.26926671931;
        H[5] := 0.295524224715
      end                
  End;
  FOR i := 1 TO n div 2 do
  begin
    j := n + 1 - i;
    A[j] := -A[i];
    H[j] := H[i]
  end;
  writeln;
  write(' Minimum value of X: '); readln(x1);
  write(' Maximum value of X: '); readln(x2);
  writeln;
  IF n2 - 1 > 0 THEN
  begin
    write(' Minimum value of Y: '); readln(y1);
    write(' Maximum value of Y: '); readln(y2);
    writeln
  end;
  IF n2 - 2 > 0 THEN
  begin
    write(' Minimum value of Z: '); readln(z1);
    write(' Maximum value of Z: '); readln(z2);
  end;
  xi9 := 0.0;
  FOR i := 1 TO n do
  begin
    IF n2 - 1 = 0 THEN GOTO 10;
    FOR j := 1 TO n do
    begin
      IF n2 - 2 = 0 THEN GOTO 10;
	FOR k := 1 TO n do
        begin
10:       a1 := -(x1 - x2) / 2;
	  b1 := (x1 + x2) / 2;
	  x := a1 * A[i] + b1;
	  IF n2 - 1 = 0 THEN GOTO 20;
	  a2 := -(y1 - y2) / 2;
	  b2 := (y1 + y2) / 2;
	  y := a2 * A[j] + b2;
	  IF n2 - 2 = 0 THEN GOTO 20;
	  a3 := -(z1 - z2) / 2;
	  b3 := (z1 + z2) / 2;
	  z := a3 * A[j] + b3;
20:       Case n2 of
            1:begin
                xi9 := xi9 + F(x,y,z) * H[i] * a1;
	        GOTO 40
              end;
            2:begin
                xi9 := xi9 + F(x,y,z) * H[i] * H[j] * a1 * a2;
	        GOTO 30
              end;
            3:xi9 := xi9 + F(x,y,z) * H[i] * H[j] * H[k] * a1 * a2 * a3
          end
	end;
30:   end;
40: end; 
  writeln;
  writeln(' Value of integral = ', xi9);
  writeln;
  Readkey; DoneWinCrt
END.

{end of file tgauss.pas}