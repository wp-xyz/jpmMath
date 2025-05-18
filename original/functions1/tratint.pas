{********************************************************
*       Polynomial Interpolation or Extrapolation       *
*              of a Discrete Function F(x)              *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Example: Function sin(x) - 2*cos(x) is given by 12   *
*          points from x=0 to x=1.1.                    *
*          Extrapolate for x=1.255).                    *
*                                                       *
*  For X             =  1.25500000000000E+0000          *
*  Estimated Y value =  3.29402327208658E-0001          *
*  Estimated Error   = -1.32858449311977E-0010          *
*  Exact Y value     =  3.29402327220005E-0001          *
*                                                       *
* ----------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific  *
*             Computing by W.H. Press, B.P. Flannery,   *
*             S.A. Teukolsky and W.T. Vetterling,       *
*             Cambridge University Press, 1986"         *
*             [BIBLI 08].                               *
*                                                       *
*                    TPW Release By J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)           *
********************************************************}   
PROGRAM TEST_POLINT;

Uses WinCrt;

Const
      NMAX=25;
      TINY=1e-25;

Type
      pVEC = ^VEC;
       VEC = Array[1..NMAX] of Double;

Var

      X, Y: pVEC;
      XX, YY, XERR: Double;
      N: Integer;


FUNCTION FCT(X:Double):Double;
Begin
  FCT:=SIN(X) - 2.0*COS(X)
END;

Procedure RATINT(XA,YA:pVEC; N:Integer; X: Double; Var Y,DY:Double);
{****************************************************
*    Interpolation or Extrapolation of a Discrete   *
*       Function By a Quotient of Polynomials       *
* ------------------------------------------------- *
* INPUTS:                                           *
*   XA:    Table of abscissas (N)                   *
*   YA:    Table of ordinates (N)                   *
*    N:    Number of points                         *
*    X:    Interpolation abscissa value             *
* OUTPUTS:                                          *
*    Y:    Returned estimation of function for X    *
*   DY:    Estimated error for Y                    *
****************************************************}
Var
    C, D: pVEC;
    DD,H,HH,T,W: Double;
    I,M,NS:Integer;
Begin

  New(C); New(D);
  NS:=1;
  HH:=ABS(X-XA^[1]);
  For I:=1 to N do
  begin
    H:=ABS(X-XA^[I]);
    IF H = 0.0 THEN
    begin
      Y:=YA^[I];
      DY:=0.0;
      Halt
    end
    else if H < HH then
    begin
      NS:=I;         {index of closest table entry}
      HH:=H
    end;
    C^[I]:=YA^[I];
    D^[I]:=YA^[I]+TINY {TINY is to prevent a zero-over-zero}
  end;                 {condition.}
  Y:=YA^[NS];
  NS:=NS-1;
  For M:=1 to N-1 do
  begin
    For I:=1 to N-M do
    begin
      W:=C^[I+1]-D^[I];
      H:=XA^[I+M]-X;          {H<>0 (tested before) }
      T:=(XA^[I]-X)*D^[I]/H;
      DD:=T-C^[I+1];
      IF DD = 0.0 then
      begin
        writeln(' Pole at requested value of X.');
        ReadKey;
        DoneWinCrt
      end;
      DD:=W/DD;
      D^[I]:=C^[I+1]*DD;
      C^[I]:=T*DD
    end;
    IF 2*NS < N-M THEN
      DY:=C^[NS+1]
    ELSE
    begin
      DY:=D^[NS];
      NS:=NS-1
    end;
    Y:=Y+DY
  end;
  Dispose(C);
  Dispose(D)

End;


{main program}
BEGIN

  New(X); New(Y);

  N := 12;     {Number of points}

  {define tables X and Y }
  X^[1] := 0.0; Y^[1]:=FCT(X^[1]);
  X^[2] := 0.1; Y^[2]:=FCT(X^[2]);
  X^[3] := 0.2; Y^[3]:=FCT(X^[3]);
  X^[4] := 0.3; Y^[4]:=FCT(X^[4]);
  X^[5] := 0.4; Y^[5]:=FCT(X^[5]);
  X^[6] := 0.5; Y^[6]:=FCT(X^[6]);
  X^[7] := 0.6; Y^[7]:=FCT(X^[7]);
  X^[8] := 0.7; Y^[8]:=FCT(X^[8]);
  X^[9] := 0.8; Y^[9]:=FCT(X^[9]);
  X^[10] := 0.9; Y^[10]:=FCT(X^[10]);
  X^[11] := 1.0; Y^[11]:=FCT(X^[11]);
  X^[12] := 1.1; Y^[12]:=FCT(X^[12]);

  {define interpolation abscissa }
  XX := 1.255;

  {call interpolation procedure }
  RATINT(X,Y,N,XX,YY,XERR);
  
  {print results }
  writeln;
  writeln(' For X             = ', XX);
  writeln(' Estimated Y value = ', YY);
  writeln(' Estimated Error   = ', XERR);
  writeln(' Exact Y value     = ', FCT(XX));
  writeln;
  ReadKey;

  Dispose(X); Dispose(Y);
  DoneWinCrt  

END.

{end of file tpolint.pas}