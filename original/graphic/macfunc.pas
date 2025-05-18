{*************************************************************
*    Approximating a function F(x) By Maclaurin's Series     *
* ---------------------------------------------------------- *
* This is the general Maclaurin's formula:                   *
*                                                            *
*                f'(0)  f"(0)  f"'(0)        fn-1(0)         *
*  f(x) = f(0) + ---- + ---- + ----- + ... + ------- + ...   *
*                 1!     2!     3!           (n-1)!          *
*                                                            *
* In the case of sin(x) between x=0 and x=2*PI, the formula  *
* becomes:                                                   *
*                                                            *
*  sin(x) = x - x^3/3! + x^5/5! - x^7/7! + x^9/9! ...        *
*                                                            *
* This small program shows that with 9 terms, we have a fair *
* approximation of sin(x) between 0 and 2*PI.                *
*                                                            *
* With option 2, you can approximate exp(x):                 *
*                                                            *
*  exp(x) = 1 + x/1! + x^2/2! + ... + x^n/n! + ...           *
* ---------------------------------------------------------- *
* REFERENCE:                                                 *
*     "Graphisme dans le plan et dans l'espace avec Turbo    *
*      Pascal 4.0 By R. Dony - MASSON, Paris 1990".          *
*                                                            *
*                         TPW Version By J-P Moreau, Paris.  *
*************************************************************}
Program MacFunc;

Uses WinCrtMy, Type_def, CrtGr2D;

Const
      B = 6.29;
      step = 0.05;

Var
      Example, NbTerms: Byte;
      Y: Array[1..256] of real_ar;

  Function Func(x:real_ar):real_ar;
  begin
    if Example=1 then
      Func := sin(x)
    else
      Func := exp(x)
  end;

  Function Fact(N:integer):real_ar;
  Var tmp:real_ar; i:integer;
  Begin
    if N < 1 then
    begin
      Fact:=1;
      exit
    end;
    tmp:=1.0;
    For i:=1 to N do tmp:=tmp*i;
    Fact:=tmp
  End;

  Function Power(X:real_ar; N:integer):real_ar;
  Var tmp:real_ar; i:integer;
  Begin
    if N < 1 then
    begin
      Power:=1;
      exit
    end;
    tmp:=1.0;
    For i:=1 to N do tmp:=tmp*X;
    Power:=tmp
  End;

  Procedure Data;
  begin
    ClrScr;
    Writeln(' Approximating SIN(X) or EXP(X) By a Maclaurin''s Series');
    Writeln(' ======================================================');
    GotoXY(1,4);
    Write(' # example (1=sin(x) 2=exp(x): '); Readln(Example);
    GotoXY(1,6);
    Write(' Input amount of terms ......: '); Readln(Nbterms)
  end;

  Procedure DrawFunc;
  var i:real_ar;
  begin
    if Example=1 then I:=0.0
                 else I:=-4.0;
    While I <= B do
    begin
      MoveXY(CrtDc,I, Func(I));
      LineXY(CrtDc,1.01*I, 1.01*Func(I));
      I := I + step
    end
  end;

  Procedure DrawFuncMaclaurin(Nbterms:Byte);
  Var term,cpt: Integer;
      i,D:integer;
      X,T: real_ar;
  Begin
    FillChar(Y,SizeOf(Y),0);
    D:=-1;
    For term:=1 to Nbterms do
    begin
      cpt:=1;
      if Example=1 then
      begin
        Inc(D,2);
        X:=0.0
      end
      else
      begin
        Inc(D,1);
        X:=-4.0
      end;
      T:=Power(X,D)/Fact(D);
      if Example=1 then
        if Odd(term) then
          Y[cpt]:=Y[cpt]+T
        else
          Y[cpt]:=Y[cpt]-T
      else
        Y[cpt]:=Y[cpt]+T;
      MoveXY(CrtDc,X,Y[cpt]);
      While X <= B do
      begin
        X := X + step;
        T:=Power(X,D)/Fact(D);
        Inc(cpt);
        if Example=1 then
          if Odd(term) then
            Y[cpt]:=Y[cpt]+T
          else
            Y[cpt]:=Y[cpt]-T
        else
          Y[cpt]:=Y[cpt]+T;
        LineXY(CrtDc,X,Y[cpt])
      end;
      GotoXY(10,2);
      Writeln(term,' term(s).');
      if term<NbTerms then ReadKey
    end
  End;


  {main program}
  BEGIN
    WinCrtInit(' MACFUNC');
    Repeat
      Data;
      ClrScr;
      if Example=1 then
        Fenetre(0, 6, -2, 2)
      else
        Fenetre(-4,2,-2,8);
      Cloture(50,MaxX-10,95,MaxY-10);
      Axes(CrtDc);
      Gradue(CrtDc,1,1);
      Grille(CrtDc,0.5,0.5);
      Bordure(CrtDc);
      DrawFunc;
      DrawFuncMaclaurin(Nbterms);
      SortieGraphique
    Until rep='n';
    DoneWinCrt
  END.

{end of file macfunc.pas}