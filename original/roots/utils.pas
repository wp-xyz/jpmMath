{====================================
 Utility functions for project unnes

 Release 1.1: Added subroutine mmulvt
 ====================================}
Unit Utils;

Interface

Uses WinCrt, Strings;

Const SIZE = 50;

Const zero = 0.0; one = 1.0; two = 2.0; ten = 10.0;

Type
      pMat = ^Mat;
      Mat = Array[0..SIZE,0..SIZE] of Double;
      pVec = ^Vec;
      Vec = Array[0..SIZE] of Double;
Var
      nfetot:Integer;
      fp_out: TEXT;

      Function Dot_Product(n:integer; a,b:pVec):double;
      Function IMAX(a,b:integer): integer;
      Function IMIN(a,b:integer): integer;
      Procedure Line0;
      Procedure Line1;
      Function log(x:double):double;
      Function MAX(a,b:double): double;
      Function MIN(a,b:double): double;
      Procedure mmulv(N:integer; A:pMat; B, C:pVec);
      Procedure mmulvt(N:integer; A:pMat; B, C:pVec);
      Procedure Msg(text:Pchar);
      Function Power(x:double; n:integer): double;
      Function Power1(y,x:double):double;
      Procedure PrintVec(n:integer; Name:pChar; V:pVec);
      Function Sign(a,b:double): double;
     

Implementation

    Function Dot_Product(n:integer; a,b:pVec):double;
    Var i:integer; sum:double;
    begin
      sum:=zero;
      For i:=1 to n do sum:=sum+a^[i]*b^[i];
      Dot_Product:=sum
    end;

    Function IMAX(a,b:integer): integer;
    begin
      if a>=b then IMAX:=a else IMAX:=b
    end;

    Function IMIN(a,b:integer): integer;
    begin
      if a<=b then IMIN:=a else IMIN:=b
    end;

    Procedure Line0;
    Var k:Word;
    Begin
      Write(fp_out,'   *');
      For k:=1 to 71 do Write(fp_out,' ');
      Writeln(fp_out,'*')
    End;

    Procedure Line1;
    Var k:Word;
    Begin
      Write(fp_out,'   *');
      For k:=1 to 71 do Write(fp_out,'-');
      Writeln(fp_out,'*')
    End;

    FUNCTION log(x:double):double;
    begin
      if x>1E-12 then
        log:= ln(x)/ln(10)
      else
        log:= -1E12
    end;

    Function MAX(a,b:double): double;
    begin
      if a>=b then MAX:=a else MAX:=b
    end;

    Function MIN(a,b:double): double;
    begin
      if a<=b then MIN:=a else MIN:=b
    end;

    Procedure Msg(text:Pchar);
    var i,l:word;
    begin
      Write(fp_out,'   *');
      Write(fp_out,text);
      l:=strlen(text);
      For i:=1 to 73-l-2 do Write(fp_out,' ');
      Writeln(fp_out,'*');
    end;

    Function Power(x:double; n:integer): double;
    {calculate x power n}
    var i:integer; result:double;
    begin
      result :=one;
      if n=0 then
      begin
        power:=result;
        exit
      end
      else
        for i:=1 to n do result := x * result;
      Power :=result
    end;

    Function Power1(y,x:double):double;
    {calculate y power x}
    begin
      IF x<0 THEN Exit;
      Power1:=Exp(x*Ln(y))
    end;

    Function Sign(a,b:double): double;
    begin
      if b>=0 then Sign:=a
              else Sign:=-a
    end;

    Procedure mmulv(N:integer; A:pMat; B, C:pVec);
    {**********************************************                                     
    * MULTIPLICATION OF A REAL MATRIX BY A VECTOR *
    * ------------------------------------------- *                                     
    * INPUTS:    A  MATRIX N*N                    *                                     
    *            B  VECTOR N*1                    *                                     
    *            N  INTEGER                       *
    * ------------------------------------------- *
    * OUTPUT:    C  VECTOR N*1  PRODUCT A*B       *                                     
    **********************************************}
    VAR 
      SUM: DOUBLE;  
      I,K: integer;
    BEGIN                                              
      for I:=1 to N do
      begin                                                                  
        SUM:=zero;                                                                
        for K:=1 to N do SUM:=SUM+A^[I][K]*B^[K];
        C^[I]:=SUM
      end                                                              
    END;

     Procedure mmulvt(N:integer; A:pMat; B, C:pVec);
    {**********************************************                                     
    * MULTIPLICATION OF THE TRANSPOSE OF A REAL   *
    * MATRIX BY A VECTOR                          *
    * ------------------------------------------- *                                     
    * INPUTS:    A  MATRIX N*N                    *                                     
    *            B  VECTOR N*1                    *                                     
    *            N  INTEGER                       *
    * ------------------------------------------- *
    * OUTPUT:    C  VECTOR N*1  PRODUCT A*B       *                                     
    **********************************************}
    VAR 
      SUM: DOUBLE;  
      I,K: integer;
    BEGIN                                              
      for I:=1 to N do
      begin                                                                  
        SUM:=zero;                                                                
        for K:=1 to N do SUM:=SUM+A^[K][I]*B^[K];                                               
        C^[I]:=SUM
      end                                                              
    END;

    {for debug only}
    Procedure PrintVec(n:integer; Name:pChar; V:pVec);
    Var i:Word;
    Begin
      Write(Name);
      For i:=1 to n do Write(' ',V^[i]);
      Writeln
    End;

END.