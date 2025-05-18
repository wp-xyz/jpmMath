{procedures used by program apollo.pas}
Unit Maths_2D;

Interface
Uses Type_Def;


  Procedure ToModeDeg;
  Function IStr(x: longint; format: byte): string;
  Procedure real_arSwap( VAR a,b: real_ar);

Implementation

  CONST  modeDegres: boolean = false;

  Procedure ToModeDeg;
  Begin
    modeDegres := true
  End;

  Function IStr;
  Var  aux: string;
  Begin
    Str(x:format,aux);
    IStr:=aux
  End;

  Procedure real_arSwap;
  Var  c: real_ar;
  Begin
    c:=a; a:=b; b:=c
  End;

END.

{end of file maths_2d.pas}