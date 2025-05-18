{********************************************************
*       Polynomial Interpolation or Extrapolation       *
*              of a discrete Function F(x)              *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Example: Function sin(x) - 2*cos(x) is given by 12   *
*          points from x=0 to x=1.1.                    *
*          Extrapolate for x=1.255).                    *
*                                                       *
*  For X             =  1.25500000000000E+0000          *
*  Estimated Y value =  3.29402327224072E-0001          *
*  Estimated Error   = -8.30468527028353E-0011          *
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

Procedure POLINT(XA,YA:pVEC; N:Integer; X: Double; Var Y,DY:Double);
{****************************************************
*     Polynomial Interpolation or Extrapolation     *
*            of a discrete Function                 *
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
    DEN,DIF,DIFT,HO,HP,W: Double;
    I,M,NS:Integer;
Begin

  New(C); New(D);
  NS:=1;
  DIF:=ABS(X-XA^[1]);
  For I:=1 to N do
  begin
    DIFT:=ABS(X-XA^[I]);
    IF DIFT < DIF THEN
    begin
      NS:=I;                {index of closest table entry}
      DIF:=DIFT
    end;
    C^[I]:=YA^[I];          {Initialize the C's and D's}
    D^[I]:=YA^[I]
  end;
  Y:=YA^[NS];               {Initial approximation of Y}
  NS:=NS-1;
  For M:=1 to N-1 do
  begin
    For I:=1 to N-M do
    begin
      HO:=XA^[I]-X;
      HP:=XA^[I+M]-X;
      W:=C^[I+1]-D^[I];
      DEN:=HO-HP;
      IF DEN = 0.0 then
      begin
        writeln;
        writeln(' Error: two identical abscissas.');
        ReadKey;
        DoneWinCrt  {stop windows program}
      end;
      DEN:=W/DEN;
      D^[I]:=HP*DEN;        {Update the C's and D's}
      C^[I]:=HO*DEN
    end;
    IF 2*NS < N-M THEN     {After each column in the tableau XA is completed,}
      DY:=C^[NS+1]         {we decide which correction, C or D, we want to   }
    ELSE                   {add to our accumulating value of Y, i.e. which   }
    begin                  {path to take through the tableau, forking up or  }
      DY:=D^[NS];          {down. We do this in such a way as to take the    }
      NS:=NS-1             {most "straight line" route through the tableau to}
    end;                   {its apex, updating NS accordingly to keep track  }
    Y:=Y+DY                {of where we are. This route keeps the partial    }
  end;                     {approximations centered (insofar as possible) on }
  Dispose(C);              {the target X.The last DY added is thus the error }
  Dispose(D)               {indication.                                      }

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
  POLINT(X,Y,N,XX,YY,XERR);
  
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