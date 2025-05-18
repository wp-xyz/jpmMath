{*********************************************
* Elementary operations with complex numbers *
* ------------------------------------------ *
* See demo program tcomplex.pas.             *
*                                            *
*       Pascal version by J-P Moreau, Paris. *
*               (www.jpmoreau.fr)            *
*********************************************}
UNIT COMPLEX2;

INTERFACE

CONST   INF  = 1.2E16;           {big number}
        TINY = 1E-16;            {small number}

{Note: PI is supposed compiler defined}

TYPE
        Complex = RECORD
                    a,b:DOUBLE; {cartesian form}
                    r,t:DOUBLE  {polar form}
                  END;

  Procedure AssignXY(VAR Z:Complex; x,y:DOUBLE);
  Procedure AssignRT(VAR Z:Complex; r,t:DOUBLE);
  Procedure DisplayComplex(Z:Complex;i,j:Integer);
  Procedure DisplayComplexR(Z:Complex;i,j:Integer);
  Procedure AddComplex(VAR Z2:Complex; Z,Z1:Complex);
  Procedure SubComplex(VAR Z2:Complex; Z,Z1:Complex);
  Procedure ChsComplex(VAR Z1:Complex; Z:Complex);
  Procedure MulComplex(VAR Z2:Complex; Z,Z1:Complex);
  Procedure DivComplex(VAR Z2:Complex; Z,Z1:Complex);
  Procedure ExpComplex(VAR Z1:Complex; Z:Complex);
  Procedure LnComplex(VAR Z1:Complex; Z:Complex);
  Procedure PowerComplex(VAR Z2:Complex; Z,Z1:Complex);
  Procedure CosComplex(VAR Z1:Complex; Z:Complex);
  Procedure SinComplex(VAR Z1:Complex; Z:Complex);
  Procedure TanComplex(VAR Z1:Complex; Z:Complex);
  Procedure ChComplex(VAR Z1:Complex; Z:Complex);
  Procedure ShComplex(VAR Z1:Complex; Z:Complex);
  Procedure ThComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArcCosComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArcSinComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArcTanComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArgChComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArgShComplex(VAR Z1:Complex; Z:Complex);
  Procedure ArgThComplex(VAR Z1:Complex; Z:Complex);


IMPLEMENTATION

  {define complex number by x and y}
  Procedure AssignXY(VAR Z:Complex; x,y:DOUBLE);
  Begin
    Z.a:=x; Z.b:=y;
    Z.r:=sqrt(x*x+y*y);
    if x=0 then
      if y>0 then Z.t:=PI/2
      else if y=0 then Z.t:=-PI/2 else Z.t:=0
    else
    begin
      Z.t:=ArcTan(y/x);
      if x<0 then
        if y>=0 then Z.t:=Z.t+PI else Z.t:=Z.t-PI;
      if Z.t>PI then Z.t:=Z.t-2*PI;
      if Z.t<-PI then Z.t:=Z.t+2*PI;
    end
  End;

  {define complex number by r and t in radians}
  Procedure AssignRT(VAR Z:Complex; r,t:DOUBLE);
  Begin
    Z.r:=r; Z.t:=t;
    if Z.t>PI then Z.t:=Z.t-2*PI;
    if Z.t<-PI then Z.t:=Z.t+2*PI;
    Z.a:=r*cos(t);
    Z.b:=r*sin(t)
  End;

  {display complex number with x and y}
  Procedure DisplayComplex(Z:Complex;i,j:Integer);
  Begin
    write('(',Z.a:i:j,',',Z.b:i:j,')');
  End;

  {display complex number with radius and phase in radians}
  Procedure DisplayComplexR(Z:Complex;i,j:Integer);
  Begin
    write('(',Z.r:i:j,',',Z.t:i:j,')');
  End;

  {add two complex numbers: Z2=Z+Z1}
  Procedure AddComplex(VAR Z2:Complex; Z,Z1:Complex);
  Begin
    Z2.a:=Z.a+Z1.a; Z2.b:=Z.b+Z1.b;
    AssignXY(Z2,Z2.a,Z2.b)
  End;

  {substract two complex numbers: Z2=Z-Z1}
  Procedure SubComplex(VAR Z2:Complex; Z,Z1:Complex);
  Begin
    Z2.a:=Z.a-Z1.a; Z2.b:=Z.b-Z1.b;
    AssignXY(Z2,Z2.a,Z2.b)
  End;

  {change sign of a complex number}
  Procedure ChsComplex(VAR Z1:Complex; Z:complex);
  Begin
    Z1.a:=-Z.a; Z1.b:=-Z.b;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  {multiply two complex numbers: Z2=Z*Z1}
  Procedure MulComplex(VAR Z2:Complex; Z,Z1:Complex);
  Begin
    Z2.r:=Z.r*Z1.r; Z2.t:=Z.t+Z1.t;
    AssignRT(Z2,Z2.r,Z2.t)
  End;

  {divide two complex numbers: Z2=Z/Z1}
  Procedure DivComplex(VAR Z2:Complex; Z,Z1:Complex);
  Begin
    if Z1.r < TINY then Z2.r := INF
                   else Z2.r := Z.r/Z1.r;
    Z2.t:=Z.t-Z1.t;
    AssignRT(Z2,Z2.r,Z2.t)
  End;

  {exponential complex function: Z1=Exp(Z) }
  Procedure ExpComplex(VAR Z1:Complex; Z:Complex);
  Var temp:double;
  Begin
    if exp(Z.a) > INF then temp:=INF else temp:=exp(Z.a); 
    Z1.a:=temp*cos(Z.b); Z1.b:=temp*sin(Z.b);
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=LN(Z) }
  Procedure LnComplex(VAR Z1:Complex; Z:Complex);
  Begin
    if Z.r<=0 then Z1.a:=-INF
              else Z1.a:=LN(Z.r);
    Z1.b:=Z.t;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  {Complex power ZZ=Z1^Z2 }
  Procedure PowerComplex(VAR Z2:Complex; Z,Z1:Complex);
  Var temp:Complex;
  Begin
    LnComplex(temp,Z);
    MulComplex(temp,Z1,temp);
    ExpComplex(Z2,temp)
  End;

  { Z1=COS(Z) }
  Procedure CosComplex(VAR Z1:Complex; Z:Complex);
  Begin
    Z1.a:=(EXP(-z.b)*COS(z.a) + EXP(z.b)*COS(-z.a))/2;
    Z1.b:=(EXP(-z.b)*SIN(z.a) + EXP(z.b)*SIN(-z.a))/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=SIN(Z) }
  Procedure SinComplex(VAR Z1:Complex; Z:Complex);
  Begin
    Z1.a:= (EXP(-z.b)*SIN(z.a) - EXP(z.b)*SIN(-z.a))/2;
    Z1.b:=-(EXP(-z.b)*COS(z.a) - EXP(z.b)*COS(-z.a))/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=TAN(Z) }
  Procedure TanComplex(VAR Z1:Complex; Z:Complex);
  Var  z2,z3: COMPLEX;
  Begin
    SinComplex(z2,Z);
    CosComplex(z3,Z);
    DivComplex(Z1,z2,z3)
  End;

  { Z1=CH(Z) }
  Procedure ChComplex(VAR Z1:Complex; Z:Complex);
  Begin
    Z1.a:=(EXP(z.a)*COS(z.b) + EXP(-z.a)*COS(-z.b))/2;
    Z1.b:=(EXP(z.a)*SIN(z.b) + EXP(-z.a)*SIN(-z.b))/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=SH(Z) }
  Procedure ShComplex(VAR Z1:Complex; Z:Complex);
  Begin
    Z1.a:= (EXP(z.a)*COS(z.b) - EXP(-z.a)*COS(-z.b))/2;
    Z1.b:=-(EXP(z.a)*SIN(z.b) - EXP(-z.a)*SIN(-z.b))/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=TH(Z) }
  Procedure ThComplex(VAR Z1:Complex; Z:Complex);
  Var  z2,z3: COMPLEX;
  Begin
    ShComplex(z2,Z);
    ChComplex(z3,Z);
    DivComplex(Z1,z2,z3)
  End;

  { Z1=ARCCOS(Z) }
  Procedure ArcCosComplex(VAR Z1:Complex; Z:Complex);
  Var temp,u:Complex;
  Begin
    temp.a:=1-z.a*z.a+z.b*z.b;
    temp.b:=-2*z.a*z.b;
    AssignXY(temp,temp.a,temp.b);
    AssignXY(u,0.5,0);
    PowerComplex(temp,temp,u);
    temp.a:=z.a-temp.b;
    temp.b:=z.b+temp.a;
    AssignXY(temp,temp.a,temp.b);
    LnComplex(temp,temp);
    Z1.a:=temp.b; Z1.b:=-temp.a;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=ARCSIN(Z) }
  Procedure ArcSinComplex(VAR Z1:Complex; Z:Complex);
  Var temp,u:Complex;
  Begin
    temp.a:=1-z.a*z.a+z.b*z.b;
    temp.b:=-2*z.a*z.b;
    AssignXY(temp,temp.a,temp.b);
    AssignXY(u,0.5,0);
    PowerComplex(temp,temp,u);
    temp.a:=temp.a-z.b;
    temp.b:=temp.b+z.a;
    AssignXY(temp,temp.a,temp.b);
    LnComplex(temp,temp);
    Z1.a:=temp.b; Z1.b:=-temp.a;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=ARCTAN(Z) }
  Procedure ArcTanComplex(VAR Z1:Complex; Z:Complex);
  Var  z2,z3: COMPLEX;
  Begin
    z2.a:=-Z.a; z2.b:=1-Z.b;
    z3.a:=Z.a; z3.b:=1+Z.b;
    AssignXY(z2,z2.a,z2.b);
    AssignXY(z3,z3.a,z3.b);
    DivComplex(z3,z2,z3);
    LnComplex(z3,z3);
    Z1.a:=z3.b/2; Z1.b:=-z3.a/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;

  { Z1=ARGCH(Z) }
  Procedure ArgChComplex(VAR Z1:Complex; Z:Complex);
  Var temp,u:Complex;
  Begin
    temp.a:=-1+Z.a*Z.a-Z.b*Z.b;
    temp.b:=2*Z.a*Z.b;
    AssignXY(temp,temp.a,temp.b);
    AssignXY(u,0.5,0);
    PowerComplex(temp,temp,u);
    temp.a:=temp.a+Z.a;
    temp.b:=temp.b+Z.b;
    AssignXY(temp,temp.a,temp.b);
    LnComplex(Z1,temp)
  End;

  { Z1=ARGSH(Z) }
  Procedure ArgShComplex(VAR Z1:Complex; Z:Complex);
  Var temp,u:Complex;
  Begin
    temp.a:=1+Z.a*Z.a-Z.b*Z.b;
    temp.b:=2*Z.a*Z.b;
    AssignXY(temp,temp.a,temp.b);
    AssignXY(u,0.5,0);
    PowerComplex(temp,temp,u);
    temp.a:=temp.a+Z.a;
    temp.b:=temp.b+Z.b;
    AssignXY(temp,temp.a,temp.b);
    LnComplex(Z1,temp)
  End;

  { Z1=ARGTH(Z) }
  Procedure ArgThComplex(VAR Z1:Complex; Z:Complex);
  Var  z2,z3: COMPLEX;
  Begin
    z2.a:=1+Z.a; z2.b:=Z.b;
    z3.a:=1-Z.a; z3.b:=-Z.b;
    AssignXY(z2,z2.a,z2.b);
    AssignXY(z3,z3.a,z3.b);
    DivComplex(z3,z2,z3);
    LnComplex(z3,z3);
    Z1.a:=z3.a/2; Z1.b:=z3.b/2;
    AssignXY(Z1,Z1.a,Z1.b)
  End;


END.

{end of file complex2.pas}