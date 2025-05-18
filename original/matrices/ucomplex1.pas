{used by program carpol1.pas}
UNIT UCOMPLEX1;

Interface

TYPE
      Complex = Record
        r,i:Double
      End;

Procedure CADD(c1,c2:Complex; Var c3:Complex);
Procedure CDIF(c1,c2:Complex; Var c3:Complex);
Procedure CMUL(c1,c2:Complex; Var c3:Complex);
Procedure CPRO(alpha:Double;C:Complex;VAR c1:Complex);


Implementation

Procedure CADD(c1,c2:Complex; Var c3:Complex);
Begin
  c3.r:=c1.r+c2.r; c3.i:=c1.i+c2.i
End;

Procedure CDIF(c1,c2:Complex; Var c3:Complex);
Begin
  c3.r:=c1.r-c2.r; c3.i:=c1.i-c2.i
End;

Procedure CMUL(c1,c2:Complex; Var c3:Complex);
Begin
  c3.r:=c1.r*c2.r-c1.i*c2.i;
  c3.i:=c1.r*c2.i+c1.i*c2.r
End;

Procedure CPRO(alpha:Double;c:Complex;VAR c1:Complex);
Begin
  c1.r:=alpha*c.r; c1.i:=alpha*c.i
End;


END.

{end of file ucomplex1.pas}

