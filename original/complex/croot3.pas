{**********************************************************
*          Solve complex equation of degree 3:            *
*             a z^3 + b z^2 + c z + d = 0                 *
* ------------------------------------------------------- *
* SAMPLE RUN:                                             *
*             a = (1,0) b = (2,-10) c = (-4,1) d = (5,2)  *
*                                                         *
* Equation az^3 + bz^2 + cz + d = 0                       *
* Solutions z:                                            *
* z = ( -0.4703,  0.6685)                                 *
* F(z) = ( -0.0000,  0.0000)                              *
* z = (  0.6018, -0.2911)                                 *
* F(z) = (  0.0000, -0.0000)                              *
* z = ( -2.1315,  9.6226)                                 *
* F(z) = ( -0.0000,  0.0000)                              *
*                                                         *
*                                     J-P Moreau, Paris.  *
*                                     (www.jpmoreau.fr)   *
***********************************************************
  NOTE:
  z = Z - (b/3a)
  p = (c/a) - (b²/3a²)
  q = (2b^3/27a^3) + (d/a) - (bc)/3a²)

  Now the equation is reduced to Z^3 + pZ + q = 0

  with the roots:

  Z[1..3] = cubroot(-q/2 + sqrt(q²/4 + p^3/27)) + cubroot(-q/2 - sqrt(q²/4 + p^3/27))
-------------------------}
Program CRoot3;

Uses WinCrt;

Type
     COMPLEX = Record
       x,y: Double;
     End;

Var
    a,b,c,d, det, p, q, sdet, u, u1, u2, v, v1, v2, w, z1: COMPLEX;
    z: Array[1..9] of COMPLEX;
    i: Integer;


{return ABS(z) }
Function CABS(z:COMPLEX): Double;
Begin
  CABS := sqrt(z.x*z.x + z.y*z.y)
End;

{ z = z1 * z2 }
Procedure CMUL(z1,z2:COMPLEX; Var z:COMPLEX);
Begin
  z.x := z1.x*z2.x - z1.y*z2.y;
  z.y := z1.x*z2.y + z1.y*z2.x
End;

{ Z = Z1 / Z2 }
Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
Var d:Real; C:Complex;
Begin
  d := Z2.x*Z2.x + Z2.y*Z2.y;
  if d<1E-15 then
    writeln(' Complex Divide by zero !')
  else
  begin
    C.x:=Z2.x; C.y:=-Z2.y;
    CMUL(Z1,C,Z);
    Z.x:=Z.x/d; Z.y:=Z.y/d
  end
End;

{ z1 = SQRT(z) }
Procedure CSQRT(z:COMPLEX; Var z1:COMPLEX);
Var r: double;
Begin
  r:=sqrt(z.x*z.x+z.y*z.y);
  z1.x := sqrt((r+z.x)/2);
  z1.y := sqrt((r-z.x)/2);
  if z.y<0 then z1.y:=-z1.y
End;

{return y^x}
Function Power(y,x:double): double;
Begin
  if y>0 then
    Power:=Exp(x*Ln(y))
  else if y=0 then
    Power:=0.0
  else
    Power:=-Exp(x*Ln(-y))
End;

Procedure ATAN(Numerator,denominator:double;
               Var Phase:double);
{Return a phase between -PI and +PI}
Begin
  if Abs(denominator) < 1E-15 then
  begin
    if (Numerator < 0) then Phase:=-PI/2
                       else Phase:=PI/2;
    exit
  end
  else
  begin
    Phase:=Arctan(Numerator/denominator);
    if denominator < 0 then Phase := Phase + Pi;
  end;
End; {ATAN}

{ z1, z2, z3 = z^(1/3) }
Procedure CUBRT(z:COMPLEX; Var z1,z2,z3:COMPLEX);
Var r,t,tt: double;
Begin
  r:=sqrt(z.x*z.x+z.y*z.y);
  r:=Power(r,1/3);
  ATan(z.y,z.x,t);
  t:=t/3.0;
  z1.x:=r*cos(t);
  z1.y:=r*sin(t);
  tt:=t-(2*PI/3);
  z2.x:=r*cos(tt);
  z2.y:=r*sin(tt);
  tt:=t+(2*PI/3);
  z3.x:=r*cos(tt);
  z3.y:=r*sin(tt)
End;

{ z1 = a*z^3 + b*z2 + c*z + d }
Procedure F(z:COMPLEX; Var z1:COMPLEX);
Begin
  CMUL(a,z,u); CMUL(z,u,v); CMUL(z,v,w);
  z1.x:=w.x; z1.y:=w.y;
  CMUL(b,z,u); CMUL(u,z,v);
  z1.x:=z1.x+v.x; z1.y:=z1.y+v.y;
  CMUL(c,z,u);
  z1.x:=z1.x+u.x+d.x; z1.y:=z1.y+u.y+d.y
End;

{main program}
BEGIN

{ seek three complex roots of complex equation:
  a z^3 + b z^2 + c z + d = 0                } 

  a.x :=  1.0; a.y :=   0.0;
  b.x :=  2.0; b.y := -10.0;
  c.x := -4.0; c.y :=   1.0;
  d.x :=  5.0; d.y :=   2.0;

  {calculate p = c/a - b²/3a² }
  CDIV(c,a,u); CMUL(b,b,v); CMUL(a,a,w);
  w.x:=3.0*w.x; w.y:=3.0*w.y; CDIV(v,w,z1);
  p.x := u.x - z1.x; p.y := u.y - z1.y;

  {calculate q = 2b^3/27a^3 + d/a - bc/3a² } 
  CMUL(b,v,w); w.x:=2.0*w.x; w.y:=2.0*w.y;
  CMUL(a,a,u); CMUL(u,a,v);
  v.x:=27.0*v.x; v.y:=27.0*v.y;
  CDIV(w,v,q);
  CDIV(d,a,u);
  q.x:=q.x+u.x; q.y:=q.y+u.y;
  CMUL(b,c,u); CMUL(a,a,v);
  v.x:=3.0*v.x; v.y:=3.0*v.y;
  CDIV(u,v,w);
  q.x:=q.x-w.x; q.y:=q.y-w.y;

  {calculate det = q²/4 + p^3/27 }
  CMUL(q,q,u); u.x:=u.x/4.0; u.y:=u.y/4.0;
  CMUL(p,p,v); CMUL(p,v,w); w.x:=w.x/27.0; w.y:=w.y/27.0;
  det.x:=u.x + w.x; det.y:=u.y + w.y;

  CSQRT(det,sdet);  {now sdet contains sqrt(q²/4 + p^3/27) }

  v.x := -q.x/2.0 + sdet.x; v.y := -q.y/2.0 + sdet.y;

  CUBRT(v,u,u1,u2);   {3 cubic roots}

  w.x := -q.x/2.0 - sdet.x; w.y := -q.y/2.0 - sdet.y;

  CUBRT(w,v,v1,v2);   {3 cubic roots}

  z[1].x := u.x + v.x; z[1].y := u.y + v.y;
  z[2].x := u.x + v1.x; z[2].y := u.y + v1.y;
  z[3].x := u.x + v2.x; z[3].y := u.y + v2.y;
  z[4].x := u1.x + v.x; z[4].y := u1.y + v.y;
  z[5].x := u1.x + v1.x; z[5].y := u1.y + v1.y;
  z[6].x := u1.x + v2.x; z[6].y := u1.y + v2.y;
  z[7].x := u2.x + v.x; z[7].y := u2.y + v.y;
  z[8].x := u2.x + v1.x; z[8].y := u2.y + v1.y;
  z[9].x := u2.x + v2.x; z[9].y := u2.y + v2.y;

{ Note: only 3 z[i] are correct solutions}

  u.x:=3.0*a.x; u.y:=3.0*a.y;
  CDIV(b,u,v);

  For i:=1 to 9 do
  begin
    z[i].x := z[i].x - v.x;
    z[i].y := z[i].y - v.y
  end;

  writeln;
  writeln(' Equation az^3 + bz^2 + cz + d = 0');
  writeln(' Solutions z:');

  For i:=1 to 9 do
  begin
    F(z[i],z1);
    if CABS(z1) < 0.001 then  {detect if solution is correct}
    begin
      writeln(' z = (',z[i].x:8:4,',',z[i].y:8:4,')');
      writeln(' F(z) = (',z1.x:8:4,',',z1.y:8:4,')')
    end
  end;

  writeln;
  ReadKey;
  DoneWinCrt

END.

{end of file croot3.pas}