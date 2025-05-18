{Procedures and Functions used by program ucomplex.pas}
UNIT COMPLEX1;

Interface

TYPE

  {complex number}
  COMPLEX = Record
              x,y: REAL;  {algebraic form}
              r,t: REAL;  {polar form}
            End;

  Procedure RectPol(VAR n:COMPLEX);
  Procedure PolRect(VAR n: COMPLEX);
  Function  ZSum(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZMinus(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZMult(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZDiv(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZExp(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZLn(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZSqr(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZSqrt(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZInv(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZPower(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZSh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZCh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZSin(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZCos(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZTan(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZTh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArcsin(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArccos(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArctan(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArgsh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArgch(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Function  ZArgth(z1: COMPLEX; VAR z: COMPLEX):Boolean;


Implementation

  {rectangular to polar conversion}
  Procedure RectPol(VAR n: COMPLEX);
  Begin
    with n do
    begin
      r:=SQRT(SQR(x)+SQR(y));
      if x=0 then
        if y>0 then t:=PI/2
        else if y<0 then t:=-PI/2 else t:=0
      else
      begin
        t:=ARCTAN(y/x);
        if x<0 then
          if y>=0 then t:=t+PI else t:=t-PI
      end
    end
  End;

  {polar to rectangular conversion}
  Procedure PolRect(VAR n: COMPLEX);
  Begin
    if n.t>PI then n.t:=n.t-2*PI;
    if n.t<-PI then n.t:=n.t+2*PI;
    n.x:=n.r*cos(n.t); n.y:=n.r*sin(n.t)
  End;

  {add two complex numbers}
  Function ZSum(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=z1.x+z2.x;
    z.y:=z1.y+z2.y;
    RectPol(z);
    ZSum:=TRUE
  End;

  {subtract two complex numbers}
  Function ZMinus(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=z1.x-z2.x;
    z.y:=z1.y-z2.y;
    RectPol(z);
    ZMinus:=TRUE
  End;

  {multiply two complex numbers}
  Function ZMult(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.r:=z1.r*z2.r;
    z.t:=z1.t+z2.t;
    PolRect(z);
    ZMult:=TRUE
  End;

  {divide two complex numbers}
  Function ZDiv(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    ZDiv:=FALSE;
    if z.r=0 then exit
    else
    begin
      z.r:=z1.r/z2.r;
      z.t:=z1.t-z2.t;
      PolRect(z);
      ZDiv:=TRUE
    end
  End;

  {compute exp(z) }
  Function ZExp(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=exp(z1.x)*cos(z1.y);
    z.y:=exp(z1.x)*sin(z1.y);
    RectPol(z);
    ZExp:=TRUE
  End;

  {compute Ln(z) }
  Function ZLn(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    if z1.r<=0 then ZLn:=FALSE
    else
    begin
      z.x:=LN(z1.r);
      z.y:=z1.t;
      RectPol(z);
      ZLn:=TRUE
    end
  End;

  {compute Z^2 }
  Function ZSqr(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=SQR(z1.x)-SQR(z1.y);
    z.y:=2*z1.x*z1.y;
    RectPol(z);
    ZSqr:=TRUE
  End;

  {compute SQRT(Z) }
  Function ZSqrt(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=SQRT(z1.r)*COS(z1.t/2);
    z.y:=SQRT(z1.r)*SIN(z1.t/2);
    RectPol(z);
    ZSqrt:=TRUE
  End;

  {compute 1/Z }
  Function ZInv(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    if z1.r=0 then ZInv:=FALSE
    else
    begin
      z.r:=1 / z1.r;
      z.t:= -z1.t;
      PolRect(z);
      ZInv:=TRUE
    end
  End;

  {compute Z1^Z2}
  Function ZPower(z1,z2: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z3: COMPLEX;
  Begin
    if ZLn(z1,z3) then
      if ZMult(z3,z2,z3) then
        if ZExp(z3,z) then ZPower:=TRUE
                      else ZPower:=FALSE
  End;

  {compute Sh(Z) }
  Function ZSh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=(EXP(z1.x)*COS(z1.y) - EXP(-z1.x)*COS(-z1.y))/2;
    z.y:=(EXP(z1.x)*SIN(z1.y) - EXP(-z1.x)*SIN(-z1.y))/2;
    RectPol(z);
    ZSh:=TRUE
  End;

  {compute Ch(Z) }
  Function ZCh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=(EXP(z1.x)*COS(z1.y) + EXP(-z1.x)*COS(-z1.y))/2;
    z.y:=(EXP(z1.x)*SIN(z1.y) + EXP(-z1.x)*SIN(-z1.y))/2;
    RectPol(z);
    ZCh:=TRUE
  End;

  {compute Sin(Z) }
  Function ZSin(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:= (EXP(-z1.y)*SIN(z1.x) - EXP(z1.y)*SIN(-z1.x))/2;
    z.y:=-(EXP(-z1.y)*COS(z1.x) - EXP(z1.y)*COS(-z1.x))/2;
    RectPol(z);
    ZSin:=TRUE
  End;

  {compute Cos(Z) }
  Function ZCos(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Begin
    z.x:=(EXP(-z1.y)*COS(z1.x) + EXP(z1.y)*COS(-z1.x))/2;
    z.y:=(EXP(-z1.y)*SIN(z1.x) + EXP(z1.y)*SIN(-z1.x))/2;
    RectPol(z);
    ZCos:=TRUE
  End;

  {compute Tan(Z) }
  Function ZTan(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var  z3,z4: COMPLEX;
  Begin
    if ZSin(z1,z3) then
      if ZCos(z1,z4) then
        if ZDiv(z3,z4,z) then ZTan:=TRUE
                         else ZTan:=FALSE
  End;

  {compute Th(Z) }
  Function ZTh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var  z3,z4: COMPLEX;
  Begin
    if ZSh(z1,z3) then
      if ZCh(z1,z4) then
        if ZDiv(z3,z4,z) then ZTh:=TRUE
                         else ZTh:=FALSE
  End;

  Function ZArcsin(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z3: COMPLEX;
  Begin
    ZArcsin:=FALSE;
    z3.x:=1-SQR(z1.x)+SQR(z1.y);
    z3.y:=-2*z1.x*z1.y;
    RectPol(z3);
    if Not ZSqrt(z3,z3) then exit;
    z3.x:=z3.x-z1.y;
    z3.y:=z3.y+z1.x;
    RectPol(z3);
    if Not ZLn(z3,z3) then exit;
    z.x:=z3.y;
    z.y:=-z3.x;
    RectPol(z);
    ZArcsin:=TRUE
  End;

  Function ZArccos(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z3: COMPLEX;
  Begin
    ZArccos:=FALSE;
    z3.x:=1-SQR(z1.x)+SQR(z1.y);
    z3.y:=-2*z1.x*z1.y;
    RectPol(z3);
    if Not ZSqrt(z3,z3) then exit;
    z3.x:=z1.x-z3.y;
    z3.y:=z1.y+z3.x;
    RectPol(z3);
    if Not ZLn(z3,z3) then exit;
    z.x:=z3.y;
    z.y:=-z3.x;
    RectPol(z);
    ZArccos:=TRUE
  End;

  Function ZArctan(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z2,z3: COMPLEX;
  Begin
    ZArctan:=FALSE;
    z2.x:=-z1.x;
    z2.y:=1-z1.y;
    z3.x:=z1.x;
    z3.y:=1+z1.y;
    RectPol(z2); RectPol(z3);
    if Not ZDiv(z2,z3,z3) then exit;
    if Not ZLn(z3,z3) then exit;
    z.x:=z3.y/2;
    z.y:=-z3.x/2;
    RectPol(z);
    ZArctan:=TRUE
  End;

  Function ZArgsh(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z3: COMPLEX;
  Begin
    ZArgsh:=FALSE;
    z3.x:=1+SQR(z1.x)-SQR(z1.y);
    z3.y:=2*z1.x*z1.y;
    RectPol(z3);
    if Not ZSqrt(z3,z3) then exit;
    z3.x:=z3.x+z1.x;
    z3.y:=z3.y+z1.y;
    RectPol(z3);
    if Not ZLn(z3,z) then exit;
    ZArgsh:=TRUE
  End;

  Function ZArgch(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z3: COMPLEX;
  Begin
    ZArgch:=FALSE;
    z3.x:=-1+SQR(z1.x)-SQR(z1.y);
    z3.y:=2*z1.x*z1.y;
    RectPol(z3);
    if Not ZSqrt(z3,z3) then exit;
    z3.x:=z3.x+z1.x;
    z3.y:=z3.y+z1.y;
    RectPol(z3);
    if Not ZLn(z3,z) then exit;
    ZArgch:=TRUE
  End;

  Function ZArgth(z1: COMPLEX; VAR z: COMPLEX):Boolean;
  Var z2,z3: COMPLEX;
  Begin
    ZArgth:=FALSE;
    z2.x:=1+z1.x;
    z2.y:=z1.y;
    z3.x:=1-z1.x;
    z3.y:=-z1.y;
    RectPol(z2); RectPol(z3);
    if Not ZDiv(z2,z3,z3) then exit;
    if Not ZLn(z3,z3) then exit;
    z.x:=z3.x/2;
    z.y:=z3.y/2;
    RectPol(z);
    ZArgth:=TRUE
  End;


END.

{end of file complex1.pas}