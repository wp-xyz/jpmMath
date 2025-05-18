{**************************************************************
*      Draw a real function F(x) by giving its equation       *
* ----------------------------------------------------------- *
* Ref.: "Mathématiques en Turbo Pascal By Alain Reverchon and *
*        Marc Ducamp, Editions Eyrolles, Paris, 1990"         *
*        [BIBLI 03].                                          *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*                                                             *
* Enter equation of Function f(x): sin(x)/x                   *
*                                                             *
*  Begin x: -4*pi                                             *
*  End x  : 4*pi                                              *
*                                                             *
* Number of points: 512                                       *
* Caption: SIN X / X                                          *
*                                                             *
* (A graph is displayed to screen).                           * 
*                                                             *
*                           TPW version by J-P Moreau, Paris. *
* ----------------------------------------------------------- *
* Note: this program demonstrates the capability of the unit  *
*       UFunct.pas to compile a user defined real function    *
*       F(x), given by its equation, and to evaluate a current*
*       point.  Then the curve is drawn by using the unit     *
*       graph_2d.                                             *
*                                                             *
* Accepted signs and operators:                               *
* ----------------------------                                *
* 0,1,2,3,4,5,6,7,8,9,(,),-,+,*,/,^,SIN,COS,TAN,ARCSIN,ARCCOS,*
* ARCTAN,SH,CH,TH,ARGSH,ARGCH,ARGTH,ABS,INT,FRAC,FACT,EXP,LN, *
* SQR (^ for power, SQR for square root).                     *
**************************************************************}
Program Test_Unit_Function;
Uses WincrtMy, Type_def, UFunct, Graph_2d;

VAR
     f: FONC;
     s: CHAINE;
     Y: RV;      {pointeur to real vector}  
     i,ndata: Integer;
     dx,x,x1,x2: real_ar;
     title: Array[0..40] of char;

{main program}
BEGIN

  WinCrtInit(' Draw a real function F(x) defined by its equation');
  writeln;
  write(' Enter equation of Function f(x): '); readln(s);
  if Not EnterFunction(s, f) then exit;
  writeln;
  write('  Begin x: '); readln(s); x1:=RealValue(s);
  write('  End x  : '); readln(s); x2:=RealValue(s);
  writeln;
  write(' Number of points: '); readln(ndata);
  write(' Caption: '); readln(title);

  {Put ndata ordinates of F(x) in table Y}
  New(Y);
  dx:=(x2-x1)/(ndata-1);
  x:=x1-dx;
  for i:=1 to ndata do
  begin
    x:=x+dx;
    Y^[i]:=Evaluate(f,x)
  end;

  {draw curve F(x) from x1 to x2}
  Clrscr;
  CourbeXY(CrtDc,ndata,10,Y,x1,x2);
  Legendes(CrtDc,title,'X','Y');
  SortieGraphique;

  {end program}
  DestroyFunction(f);
  Dispose(Y);
  DoneWinCrt
END.

{end of file gfunct.pas}