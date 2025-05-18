{*****************************************************************
*    Complex Bessel Function of the 1st Kind of integer order    *
* -------------------------------------------------------------- *
* SAMPLE RUN:                                                    *
*                                                                *
* Complex Bessel Function of the 1st Kind of integer order       *
*                                                                *
* Input complex argument (real imaginary): 1 2                   *
* Input integer order: 1                                         *
*                                                                *
* Function value:  1.2918475193703E+0000  1.01048836515558E+0000 *
*                                                                *
*                                                                *
*                          TPW Release 1.2 By J-P Moreau, Paris. *
*                                    (www.jpmoreau.fr)           *
* -------------------------------------------------------------- *
* Release 1.1: Corrected bug in Function ZLn (ArcTan replaced by *
*              Procedure ATAN2) 11/10/2005.                      *
* Release 1.2: ZPower replaced by IZpower (integer exponant).    *
*              Limitations increased in CDIV, MAXK=20            *
*****************************************************************}
Program Test_CBESSJ;

Uses WinCrt;

Const
      MAXK = 20;    {09/21/2009}

      HALF = 0.5;
      ONE  = 1.0;
      FPF  = 5.5;

Type Complex = Array[1..2] of Double;

Var  nu: Integer;
     z,z1: Complex;


  {Z=Z1/Z2}
  Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
  Var D:double;
  Begin
    D:=Z2[1]*Z2[1]+Z2[2]*Z2[2];
    if D<1e-20 then exit;
    Z[1]:=(Z1[1]*Z2[1]+Z1[2]*Z2[2])/D;
    Z[2]:=(Z1[2]*Z2[1]-Z1[1]*Z2[2])/D
  End;

  {Z=Z1*Z2}
  Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z[1]:=Z1[1]*Z2[1] - Z1[2]*Z2[2];
    Z[2]:=Z1[1]*Z2[2] + Z1[2]*Z2[1]
  End;

  {compute Z^N }
  Procedure IZPower(z: COMPLEX; n: integer; Var z1:COMPLEX);
  Var temp,temp1: COMPLEX; i: integer;
  Begin
    if n=0 then
    begin
      z1[1]:=1.0;
      z1[2]:=0.0
    end
    else if n=1 then
    begin
      z1[1]:=z[1];
      z1[2]:=z[2]
    end
    else
    begin
      temp1[1]:=z[1]; temp1[2]:=z[2];
      for i:=2 to n do
      begin
	CMUL(temp1,z,temp);
	temp1[1]:=temp[1];
	temp1[2]:=temp[2]
      end;
      z1[1]:=temp[1];
      z1[2]:=temp[2]
    end
  End;

  Function Fact(k:Integer):Double;
  Var i:Integer;
      f:Double;
  Begin
    F:=1.0;
    for i:=2 to k do f:=f*i;
    Fact:=f
  End;

{******************************************
*           FUNCTION  GAMMA(X)            *
* --------------------------------------- *
* Returns the value of Gamma(x) in double *
* precision as EXP(LN(GAMMA(X))) for X>0. *
******************************************}
Function Gamma(xx:double):double;
Var
  cof:Array[1..6] of double;
  stp,x,tmp,ser:double;
  j:integer;
Begin
  cof[1]:=76.18009173;
  cof[2]:=-86.50532033;
  cof[3]:=24.01409822;
  cof[4]:=-1.231739516;
  cof[5]:=0.120858003e-2;
  cof[6]:=-0.536382e-5;
  stp:=2.50662827465;
  
  x:=xx-ONE;
  tmp:=x+FPF;
  tmp:=(x+HALF)*LN(tmp)-tmp;
  ser:=one;
  for j:=1 to 6 do
  begin
    x:=x+ONE;
    ser:=ser+cof[j]/x
  end;
  Gamma := EXP(tmp+LN(stp*ser))
End;


Procedure CBESSJ(z:Complex;nu:Integer;Var z1:Complex);
{---------------------------------------------------
                        inf.     (-z^2/4)^k
    Jnu(z) = (z/2)^nu x Sum  ------------------
                        k=0  k! x Gamma(nu+k+1)
   (nu must be >= 0).
----------------------------------------------------}
Var k:Integer;
    sum,tmp,tmp1:Complex;
Begin
  sum[1]:=0.0; sum[2]:=0.0;
  For k:=0 to MAXK do
  begin
    {calculate (-z^2/4)^k }
    CMUL(z,z,tmp);
    tmp[1]:=-tmp[1]; tmp[2]:=-tmp[2];
    tmp1[1]:=4.0; tmp1[2]:=0.0;
    CDIV(tmp,tmp1,tmp);
    IZPower(tmp,k,tmp);
    {divide by k! }
    tmp1[1]:=Fact(k);
    CDIV(tmp,tmp1,tmp);
    {divide by Gamma(nu+k+1) }
    tmp1[1]:=Gamma(1.0*(nu+k+1));
    CDIV(tmp,tmp1,tmp);
    {actualize sum}
    sum[1]:=sum[1]+tmp[1];
    sum[2]:=sum[2]+tmp[2];
  end;
  {calculate (z/2)^nu }
  tmp[1]:=2.0; tmp[2]:=0.0;
  CDIV(z,tmp,tmp);
  IZPower(tmp,nu,tmp);
  {multiply (z/2)^nu by sum }
  CMUL(tmp,sum,z1)
End;


BEGIN
  Writeln;
  Writeln(' Complex Bessel Function of the 1st Kind of integer order');
  writeln;
  write(' Input complex argument (real imaginary): ');
  readln(z[1],z[2]);
  write(' Input integer order: '); readln(nu);
  writeln;

  CBESSJ(z,nu,z1);

  writeln(' Function value: ', z1[1],' ',z1[2]);
  writeln;
  ReadKey;
  DoneWinCrt

END.

{end of file cbessj.pas}