=============================================
 Utility functions for project unnes
 -----------------------------------
 In this version of utils.pas, output is sent
 to screen instead of output text file. This 
 unit can be used by other programs, such as
 tequdif.pas.
 ============================================}
Unit Utils1;

Interface

Uses WinCrt1, Strings;

Const SIZE = 50;

Const zero = 0.0; one = 1.0; two = 2.0; ten = 10.0;

Type
      pMat1 = ^Mat1;
      Mat1 = Array[0..SIZE,0..SIZE] of Double;
      pVec1 = ^Vec1;
      Vec1 = Array[0..SIZE] of Double;
      pIVec1 = ^IVec1;
      IVec1 = Array[0..SIZE] of Integer;
Var
      nfetot:Integer;
      fp_out: TEXT;

      Function Dot_Product(n:integer; a,b:pVec1):double;
      Function IMAX(a,b:integer): integer;
      Function IMIN(a,b:integer): integer;
      Procedure Line0;
      Procedure Line1;
      Function log(x:double):double;
      Procedure matprt (nrowpr, ncolpr: integer; a: pMat1);
      Function MAX(a,b:double): double;
      Function MIN(a,b:double): double;
      Procedure mmulv(N:integer; A:pMat1; B, C:pVec1);
      Procedure mmulvt(N:integer; A:pMat1; B, C:pVec1);
      Procedure Msg(text:String);
      Procedure IMsg(text:String; N:integer);
      Procedure RMsg(text:String; X:double);
      Function Power(x:double; n:integer): double;
      Function Power1(y,x:double):double;
      Procedure PrintVec(n:integer; Name:String; V:pVec1);
      Function Sign(a,b:double): double;
     

Implementation

    Function Dot_Product(n:integer; a,b:pVec1):double;
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
      Write('   *');
      For k:=1 to 71 do Write(' ');
      Writeln('*')
    End;

    Procedure Line1;
    Var k:Word;
    Begin
      Write('   *');
      For k:=1 to 71 do Write('-');
      Writeln('*')
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

    Procedure Msg(text:String);
    var i,l:word;
    begin
      Write('   *');
      Write(text);
      l:=length(text);
      For i:=1 to 73-l-2 do Write(' ');
      Writeln('*');
    end;

    Procedure RMsg(text:String; X:double);
    var i,l:word; s,s1:string;
    begin
      Write('   *');
      if (Abs(X)<100000) and (Abs(X)>1E-5) then Str(X:12:6,s1)
      else Str(X:13:-6,s1);
      s:=text+' '+s1;
      Write(s);
      l:=length(s);
      For i:=1 to 73-l-2 do Write(' ');
      Writeln('*');
    end;

    Procedure IMsg(text:String; N:integer);
    var i,l:word; s,s1:string;
    begin
      Write('   *');
      Str(N,s1);
      s:=text+' '+s1;
      Write(s);
      l:=length(s);
      For i:=1 to 73-l-2 do Write(' ');
      Writeln('*');
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

    Procedure mmulv(N:integer; A:pMat1; B, C:pVec1);
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

     Procedure mmulvt(N:integer; A:pMat1; B, C:pVec1);
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

Procedure matprt (nrowpr, ncolpr: integer; a: pMat1);
{-----------------------------------------------------------------------------
!    FEB. 6, 1991
!
!    THIS SUBROUTINE PRINTS RECTANGULAR BLOCKS STARTING WITH ELEMENT A^[1,1]
!    OF SIZE NROWPR BY NCOLPR FOR MATRIX A ^[WHICH HAS DECLARED SIZE NROWA BY
!    NCOLA].  THE MATRIX IS PRINTED AS A BLOCK FOR SIZES UP TO 5X5 OR BY
!    COLUMNS IF IT IS LARGER.
!
!    NROWPR IS THE NUMBER OF ROWS TO BE PRINTED
!    NCOLPR IS THE NUMBER OF COLUMNS TO BE PRINTED
!
!    IF MATRIX PRINTING IS TO BE SUPPRESSED THEN LOGICAL
!    VARIABLE MATSUP MUST BE SET TO TRUE BEFORE THE CALL TO NNES.
!----------------------------------------------------------------------------}
Label return;
Var
  i, j, k, limit: integer;
  s:string[80]; s1:string[20]; s2:string[4];
Begin

{ FOR NCOLPR <= 5 WRITE MATRIX AS A WHOLE }

  s:=''; s1:=''; s2:='';
  Line0;
  IF ncolpr <= 5 THEN
  begin
    s:='           ';
    For k:= 1 to ncolpr do
    begin
      Str(k:3,s2); s:=s+s2; s:=s+'        ';
    end;
    Msg(s); Line0;
    For i:=1 to nrowpr do
    begin
      s:='';
      Str(i:3,s2); s:=s+s2; s:=s+'     ';
      For k:= 1 to ncolpr do
      begin
        Str(a^[i,k]:10:3,s1); s:=s+s1; s:=s+' '
      end;
      Msg(s);
    end;
    Line0
  end
  ELSE
  begin

{   LIMIT IS THE NUMBER OF GROUPS OF 5 COLUMNS }
    limit:=ncolpr Div 5;

{   WRITE COMPLETE BLOCKS FIRST (LEFTOVERS LATER) }
    For j:=1 to limit do
    begin
      s:='           ';
      For k:= 1+(j-1)*5 to 5+(j-1)*5 do
      begin
        Str(k:3,s2); s:=s+s2; s:=s+'        ';
      end;
      Msg(s); Line0;
      For i:=1 to nrowpr do
      begin
        s:='';
        Str(i:3,s2); s:=s+s2; s:=s+'     ';
        For k:= 1+(j-1)*5 to 5+(j-1)*5 do
        begin
          if abs(a^[i,k]) < 100000 then
          begin
            Str(a^[i,k]:10:3,s1); s:=s+s1; s:=s+' '
          end
          else
          begin
            Str(a^[i,k]:10:-3,s1); s:=s+s1; s:=s+' '
          end
        end;
        Msg(s)
      end;
      Line0
    end;

{   WRITE REMAINING ELEMENTS }
    s:='           ';
    For k:= 5*limit+1 to ncolpr do
    begin
      Str(k:3,s2); s:=s+s2; s:=s+'        ';
    end;
    Msg(s); Line0;
    For i:=1 to nrowpr do
    begin
      s:='';
      Str(i:3,s1); s:=s+s1; s:=s+'     ';
      For k:= 5*limit+1 to ncolpr do
      begin
        if abs(a^[i,k]) < 100000.0 then
        begin
          Str(a^[i,k]:10:3,s1); s:=s+s1; s:=s+' '
        end
        else
        begin
          Str(a^[i,k]:10:-3,s1); s:=s+s1; s:=s+' '
        end
      end;
      Msg(s)
    end;
    Line0
  end;

return: End; {matprt}


Procedure PrintVec(n:integer; Name:String; V:pVec1);
Var i:Word;
Begin
  Msg(Name);
  For i:=1 to n do RMsg('', V^[i]);
  Line0
End;

END.