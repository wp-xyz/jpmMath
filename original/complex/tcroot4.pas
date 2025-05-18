{**************************************************************
*            Solve complex equation of degree 4:              *
*            a z^4 + b z^3 + c z^2 + d z + e = 0              *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*    a = (1,0) b = (2,-10) c = (-4,1) d = (5,2) e = (3,-7.5)  *
*                                                             *
* Equation az^4 + bz^3 + cz^2 + d z + e = 0                   *
* Solutions z1..z4:                                           *
*                                                             *
* z1 = (-2.13984854250220E+0000, 9.62514245773279E+0000)      *
* F(z1) = (-2.94946289614018E-0011,-4.11546352552250E-0011)   *
* z2 = ( 6.90192380071045E-0001,-7.06379865429986E-0001)      *
* F(z2) = ( 4.27213819875760E-0013,-1.74571468392060E-0012)   *
* z3 = (-8.75877413322122E-0001, 2.23149298527118E-0001)      *
* F(z3) = ( 1.45794487593776E-0012,-1.03821395924797E-0011)   *
* z4 = ( 3.25533575753275E-0001, 8.58088109170080E-0001)      *
* F(z4) = (-4.44166925461786E-0012,-7.68629604408488E-0012)   *
*                                                             *
* ----------------------------------------------------------- *
* Reference: Mathématiques et statitiques - Programmes en     *
*            BASIC. Editions du P.S.I., 1981.                 *
*                                                             *
*            Adapted to complex Domain By J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)                 *
***************************************************************
 Explanations:
 ------------
 Let us solve complex equation z^4 + a*z^3 + b*z^2 + c*z + d = 0  (1)

 Let us set z = z1 - a/4 and equation (1) becomes:

     (z1²+2ckz1+l)(z1²-2ckz1+m) = 0       (2)

 where:   l + m -4ck² = q  |            q = 3a²/8
          2 ck (m-1)  = r  | (3) with:  r = c-(ab/2)+(a^3/8)
                 lm   = s  |            s = d-(ac/4)+(a²b/16)-(3a^4/256)

 ck² = zz is a solution of complex cubic equation:

          zz^3 + aa zz^2 + bb z + cc = 0   (4)

 with:    aa = q/2, bb = (q²-4s)/16, cc = -r²/64

 So we solve equation (4) using Croot3 (3 roots zz1, zz2, zz3), then
 we calculate complex l and m by solving system (3), finally, we solve
 both 2nd degree equations of product (2). We obtain three sets of four
 complex roots, only one of which is valid, i.e. annulates equation (1). 

-----------------------------------------------------------------------}
Program Test_CRoot4;

Uses WinCrt1;

Type
     COMPLEX = Record
       x,y: Double;
     End;

Var
    a,b,c,d,e,z1,z2,z3,z4,zz: COMPLEX;


{ return Abs(z) }
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

{ z = z1 / z2 }
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

{return y^x } 
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
{Returna a phase between -PI and +PI}
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

{ z1 = z^1/3) }
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

  {return a*z^4+b*z^3+c*z^2+d^z+e }
  Procedure F4(a,b,c,d,e,z:Complex; Var z1:Complex);
  Var u,v,w: Complex;
  begin
    CMUL(z,z,u); CMUL(u,u,v); CMUL(a,v,z1);
    CMUL(u,z,v); CMUL(b,v,w);
    z1.x:=z1.x+w.x; z1.y:=z1.y+w.y;
    CMUL(c,u,v);
    z1.x:=z1.x+v.x; z1.y:=z1.y+v.y;
    CMUL(d,z,u);
    z1.x:=z1.x+u.x+e.x; z1.y:=z1.y+u.y+e.y
  end;

{return complex roots of equation a*z^2+b*z+c = 0 }
Procedure Equa2c(a,b,c: COMPLEX; VAR s1, s2: COMPLEX);
Var u,v,w:COMPLEX; r: Double;
Begin
  u.x := b.x * b.x - b.y * b.y - 4 * a.x * c.x + 4 * a.y * c.y;
  u.y := 2 * b.x * b.y - 4 * a.x * c.y - 4 * a.y * c.x;
  r := sqrt(u.x*u.x + u.y*u.y);
  v.x := sqrt((r+u.x)/2);
  v.y := sqrt((r-u.x)/2);
  if u.y<0 then v.y:=-v.y;
  w.x := (-b.x - v.x)/2;
  w.y := (-b.y - v.y)/2;
  u.x := (-b.x + v.x)/2;
  u.y := (-b.y + v.y)/2;
  r := a.x * a.x + a.y * a.y;
  s1.x := (a.x * w.x + a.y * w.y)/r;
  s1.y := (a.x * w.y - a.y * w.x)/r;
  s2.x := (a.x * u.x + a.y * u.y)/r;
  s2.y := (a.x * u.y - a.y * u.x)/r;
End;

{solve complex equation a*z^3+b*z^2+c*z+d=0 and return
 3 solutions z1,z2,z3)  }
Procedure Croot3(a,b,c,d:Complex; Var z1,z2,z3:Complex);
{ z = Z - (b/3a)
  p = (c/a) - (b²/3a²)
  q = (2b^3/27a^3) + (d/a) - (bc)/3a²)

  Now the equation is reduced to Z^3 + pZ + q = 0

  with the roots:

  Z[1..3] = cubroot(-q/2 + sqrt(q²/4 + p^3/27)) + cubroot(-q/2 - sqrt(q²/4 + p^3/27))  }

Var p,q,det,sdet,u,u1,u2,v,v1,v2,w,zz: Complex;
    z: Array[1..9] of COMPLEX;  i,index:integer;

  { z1 = z^3 + p*z + q }
  Procedure F1(z:COMPLEX; Var z1:COMPLEX);
  Begin
    CMUL(z,z,u); CMUL(z,u,v); 
    CMUL(p,z,w);
    z1.x:=v.x + w.x + q.x;
    z1.y:=v.y + w.y + q.y
  End;

Begin
  {calculate p}
  CDIV(c,a,u); CMUL(b,b,v); CMUL(a,a,w);
  w.x:=3.0*w.x; w.y:=3.0*w.y; CDIV(v,w,z1);
  p.x := u.x - z1.x; p.y := u.y - z1.y;

  {calculate q}
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

  CMUL(q,q,u); u.x:=u.x/4.0; u.y:=u.y/4.0;
  CMUL(p,p,v); CMUL(p,v,w); w.x:=w.x/27.0; w.y:=w.y/27.0;
  det.x:=u.x + w.x; det.y:=u.y + w.y;

  CSQRT(det,sdet);  {now sdet contains sqrt(q²/4 + p^3/27)) }

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

  index:=1;
  For i:=1 to 9 do
  begin
    F1(z[i],zz);
    if CABS(zz) < 0.001 then  {detect if solution is correct}
    begin
      if index=1 then
      begin
        z1.x:=z[i].x; z1.y:=z[i].y; Inc(index)
      end
      else if index=2 then
      begin
        z2.x:=z[i].x; z2.y:=z[i].y; Inc(index)
      end
      else
      begin
        z3.x:=z[i].x; z3.y:=z[i].y
      end
    end
  end;

  u.x:=3.0*a.x; u.y:=3.0*a.y;
  CDIV(b,u,v);

  z1.x := z1.x - v.x; z1.y := z1.y - v.y;
  z2.x := z2.x - v.x; z2.y := z2.y - v.y;
  z3.x := z3.x - v.x; z3.y := z3.y - v.y

End; {croot3}

{solve complex equation a*z^4+b*z^3+c*z^2+d*z+e=0 and return
 4 solutions z1,z2,z3,z4)  }
Procedure Croot4(a,b,c,d,e:Complex; Var z1,z2,z3,z4:Complex);
Label 100,200,300,400;
Var a1,aa,b1,bb,c1,cc, q,r,s,t,u,v,w: Complex;
    ck,l,m: Complex;
    z:Array[1..3] of Complex;
    zz:Array[1..12] of Complex;
    i: integer;

  {return z^4+a*z^3+b*z^2+c^z+d }
  Procedure F(a,b,c,d,z:Complex; Var z1:Complex);
  begin
    CMUL(z,z,u); CMUL(u,u,v); z1:=v;
    CMUL(u,z,v); CMUL(a,v,w);
    z1.x:=z1.x+w.x; z1.y:=z1.y+w.y;
    CMUL(b,u,v);
    z1.x:=z1.x+v.x; z1.y:=z1.y+v.y;
    CMUL(c,z,u);
    z1.x:=z1.x+u.x+d.x; z1.y:=z1.y+u.y+d.y
  end;

Begin
  t:=a;
  CDIV(b,t,a); CDIV(c,t,b); CDIV(d,t,c); CDIV(e,t,d);
{ Now equation is: z^4+a*z^3+b*z^2+c*z+d = 0

  q = b - (3.0 * a * a / 8.0)  }
  CMUL(a,a,u); u.x:=u.x*3.0/8.0; u.y:=u.y*3.0/8.0;
  q.x:=b.x - u.x; q.y:=b.y - u.y;

{ r = c - (a * b / 2) + (a * a * a / 8.0)  }
  CMUL(a,b,u);  r.x:=-u.x/2.0; r.y:=-u.y/2.0;
  CMUL(a,a,u); CMUL(a,u,v); v.x:=v.x/8.0; v.y:=v.y/8.0;
  r.x:=r.x+c.x + v.x; r.y:=r.y+c.y + v.y;

{ s = d - (a * c / 4) + (a * a * b / 16.0) - (3 * a * a * a * a / 256.0) }
  CMUL(a,c,u); u.x:=u.x/4.0; u.y:=u.y/4.0;
  CMUL(a,a,v); CMUL(b,v,w); v.x:=w.x/16.0; v.y:=w.y/16.0;
  CMUL(a,a,t); CMUL(t,t,w); w.x:=w.x*3.0/256.0; w.y:=w.y*3.0/256.0;
  s.x:=d.x - u.x + v.x -w.x; s.y:=d.y - u.y + v.y -w.y;

{ Défine coefficients of cubic equation z^3+aa*z^2+bb*z+cc=0 (1)
  aa = q / 2.0  }
  aa.x:=q.x/2.0; aa.y:=q.y/2.0;
{ bb = (q * q - 4 * s) / 16.0  }
  CMUL(q,q,u);
  bb.x:=(u.x-4.0*s.x)/16.0; bb.y:=(u.y-4.0*s.y)/16.0;
{ cc = -(r * r / 64.0)  }
  CMUL(r,r,u);
  cc.x:=-u.x/64.0; cc.y:=-u.y/64.0;

{ Calculate complex roots equation (1) }
  a1.x:=1.0; a1.y:=0.0;
  IF (CABS(r) > 1E-15) OR (CABS(bb) > 1E-15) THEN GOTO 100;

{ Particular case when equation (1) is of 2nd order }
  b1:=aa; c1:=bb;
  Equa2c(a1,b1,c1,z[1],z[2]);
  z[3].x:=0.0; z[3].y:=0.0;
  GOTO 200;

100: Croot3(a1,aa,bb,cc,z[1],z[2],z[3]);

200: For i:=1 to 3 do
  begin
    CSQRT(z[i],ck);
    {Calculate L and M if k=0}
    IF CABS(ck) = 0 THEN
    begin
      {r = SQRT(q * q - 4 * s)}
      CMUL(q,q,u); v.x:=4.0*s.x;  v.y:=4.0*s.y;
      r.x := u.x - v.x; r.y := u.y - v.y;
      GOTO 300
    end;
    {q = q + (4 * z[i]) }
    u.x:=4.0*z[i].x; u.y:=4.0*z[i].y;
    q.x := q.x + u.x; q.y := q.y + u.y;
    {r = r / (2 * ck)  }
    u.x:=2.0*ck.x; u.y:=2.0*ck.y;
    CDIV(r,u,r);
300:{l = (q - r) / 2.0, m = (q + r) / 2.0 }
    l.x:=(q.x-r.x)/2.0; l.y:=(q.y-r.y)/2.0;
    m.x:=(q.x+r.x)/2.0; m.y:=(q.y+r.y)/2.0;

    {Solving two equations of degree 2}
    b1.x := 2.0 * ck.x; b1.y := 2.0 * ck.y;
    {1st equation}
    Equa2c(a1,b1,l,u,v);

    b1.x:=-b1.x; b1.y:=-b1.y;
    {2nd equation}
    Equa2c(a1,b1,m,w,t);

    {Transferring solutions in zz(i) }
    Case i of
      1: begin zz[1]:=u; zz[2]:=v; zz[3]:=w; zz[4]:=t end;
      2: begin zz[5]:=u; zz[6]:=v; zz[7]:=w; zz[8]:=t end;
      3: begin zz[9]:=u; zz[10]:=v; zz[11]:=w; zz[12]:=t end
    end; {Note: only 4 solutions are correct}
  end; {i loop}

  u.x:=a.x/4.0; u.y:=a.y/4.0;
  For i:=1 to 12 do  {shift zz by -a/4}
  begin
    zz[i].x:=zz[i].x-u.x; zz[i].y:=zz[i].y-u.y
  end;
  {select 4 correct solutions}
  i:=1;
400: F(a,b,c,d,zz[i],t);
  if CABS(t) < 0.0001 then
    Case i of
      1: begin z1:=zz[1]; z2:=zz[2]; z3:=zz[3]; z4:=zz[4] end;
      5: begin z1:=zz[5]; z2:=zz[6]; z3:=zz[7]; z4:=zz[8] end;
      9: begin z1:=zz[9]; z2:=zz[10]; z3:=zz[11]; z4:=zz[12] end
    end;
  Inc(i,4);
  if i<=9 then goto 400;

  {now complex roots are in z1,z2,z3,z4}   

End;

{main program}
BEGIN

  writeln;
  writeln(' Equation az^4 + bz^3 + cz^2 + d z + e = 0');
  writeln(' Solutions z1..z4:');

  a.x :=  1.0; a.y :=   0.0;
  b.x :=  2.0; b.y := -10.0;
  c.x := -4.0; c.y :=   1.0;
  d.x :=  5.0; d.y :=   2.0;
  e.x :=  3.0; e.y :=  -7.5;

  Croot4(a,b,c,d,e, z1,z2,z3,z4);

  writeln;
  writeln(' z1 = (',z1.x,',',z1.y,')');
  F4(a,b,c,d,e,z1,zz);
  writeln(' F(z1) = (',zz.x,',',zz.y,')');
  writeln(' z2 = (',z2.x,',',z2.y,')');
  F4(a,b,c,d,e,z2,zz);
  writeln(' F(z2) = (',zz.x,',',zz.y,')');
  writeln(' z3 = (',z3.x,',',z3.y,')');
  F4(a,b,c,d,e,z3,zz);
  writeln(' F(z3) = (',zz.x,',',zz.y,')');
  writeln(' z4 = (',z4.x,',',z4.y,')');
  F4(a,b,c,d,e,z4,zz);
  writeln(' F(z4) = (',zz.x,',',zz.y,')');
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file tcroot4.pas}