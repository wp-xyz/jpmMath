{**************************************************************
*   Calculate the formal derivative of a real function F(x)   *
* ----------------------------------------------------------- *
* Ref.: *Mathematiques en Turbo Pascal by Alain Reverchon and *
*        Marc Ducamp, Armand Colin Editor, Paris, 1987".      *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*  Formal derivation of a real function F(x):                 *
*  Function f(x) to derivate: sin(2*x)                        *
*  The formal derivative is f(x) = 2*COS(2*X)                 *
*                                                             *
*  Formal derivation of a real function F(x):                 *
*  Function f(x) to derivate: sqr(1-x)                        *
*  The formal derivative is f(x) = -0.5/SQR(1-X)              *
*                                                             *
*                   TPW English Version By J-P Moreau, Paris. *
**************************************************************}
Program Test_FunctionDeriv;

Uses Wincrt, UFunct, Simplif2;

VAR
     f, g: FONC;
     s   : CHAINE;
     ch  : MAXSTRING;

Function DerivFunction(f:FONC; VAR g:FONC): Boolean;
Var g1,g2: FONC;
Begin
  DerivFunction := FALSE;
  case f.typ of
    CTE: begin
           g.typ:=CTE;
           g.value:=0
         end;
    VARIABLE: begin
                g.typ:=CTE;
                g.value:=1
              end;
    OPUNAIRE:
    begin
      if Not DerivFunction(f.left^,g1) then exit;
      f1:=f.left^;
      f2:=g1;
      Case f.oper1 of
        FMOINS:
          if Not CompileFunction('-U''', g) then exit;
        FPLUS:
          if Not CompileFunction('U''', g) then exit;
        FEXP:
          if Not CompileFunction('U''*EXP(U)', g) then exit;
        FLOG:
          if Not CompileFunction('U''/U', g) then exit;
        FSIN:
          if Not CompileFunction('U''*COS(U)', g) then exit;
        FCOS:
          if Not CompileFunction('-U''*SIN(U)', g) then exit;
        FTAN:
          if Not CompileFunction('U''/(COS(U)^2)', g) then exit;
        FSH:
          if Not CompileFunction('U''*CH(U)', g) then exit;
        FCH:
          if Not CompileFunction('U''*SH(U)', g) then exit;
        FTH:
          if Not CompileFunction('U''/(CH(U)^2)', g) then exit;
        FSQR:
          if Not CompileFunction('U''/2/SQR(U)', g) then exit;
        FARCTAN:
          if Not CompileFunction('U''/(1+U^2)', g) then exit;
        FARGTH:
          if Not CompileFunction('U''/(1-U^2)', g) then exit;
        FARCSIN:
          if Not CompileFunction('U''/SQR(1-U^2)', g) then exit;
        FARCCOS:
          if Not CompileFunction('-U''/SQR(1-U^2)', g) then exit;
        FARGSH:
          if Not CompileFunction('U''/SQR(U^2+1)', g) then exit;
        FARGCH:
          if Not CompileFunction('U''/SQR(U^2-1)', g) then exit;
        FABS,FINT,FFRAC: exit;  {no derivative!}
      End;
      DestroyFunction(g1)
    End;            
    OPBINAIRE:
    begin
      if Not DerivFunction(f.left^, g1) then exit;
      if Not DerivFunction(f.right^, g2) then exit;
      f1:=f.left^;
      f2:=g1;
      f3:=f.right^;
      f4:=g2;
      Case f.oper2 of
        FSUM:
          if Not CompileFunction('U''+V''', g) then exit;
        FSUB:
          if Not CompileFunction('U''-V''', g) then exit;
        FMUL:
          if Not CompileFunction('U''*V + U*V''', g) then exit;
        FDIV:
          if Not CompileFunction('U''*V-U*V'')/(V^2)', g) then exit;
        FPUI:
          if f3.typ = CTE then
            if f3.value=0 then
            begin
              g.typ:=CTE;
              g.value:=0
            end
            else
              if f3.value=1 then
                if Not DerivFunction(f3, g) then exit
              else
                if Not CompileFunction ('U''*V*(U^(V-1))', g) then exit
          else
            if Not CompileFunction
                  ('(U^V)*(V''*LOG(U)+V*U''/U)', g) then exit
      end;
      DestroyFunction(g1);
      DestroyFunction(g2)
    end
  end;
  SimplifyFunction (g);
  DerivFunction:=TRUE
End;

{main program}
BEGIN
  writeln;
  writeln(' Formal derivation of a real function F(x):');
  writeln;
  write(' Function f(x) to derivate: '); readln(s);
  if Not EnterFunction(s, f) then exit;
  if Not DerivFunction(f, g) then
    writeln(' Impossible to derivate !!!')
  else
  begin
    DeCompileFunction(g, ch);
    writeln;
    writeln(' The formal derivative is f(x) = ', ch)
  end;
  DestroyFunction(f);
  DestroyFunction(g);
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file deriv.pas}