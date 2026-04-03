{*******************************************************
*   Program to demonstrate the Akima spline fitting    *
*       of Function SIN(X) in double precision         *
* ---------------------------------------------------- *
*   Reference: BASIC Scientific Subroutines, Vol. II   *
*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                      *
*             Pascal Version by J.-P. Moreau, Paris.   *
*                       (www.jpmoreau.fr)              *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* Akima spline fitting of SIN(X):                      *
*                                                      *
*   X   SIN(X) HANDBOOK   AKIMA INTERPOLATION   ERROR  *
* ---------------------------------------------------- *
* 0.00    0.0000000          0.0000000       0.0000000 *
* 0.05    0.0499792          0.0500402      -0.0000610 *
* 0.10    0.0998334          0.0998435      -0.0000101 *
* 0.15    0.1494381          0.1494310       0.0000072 *
* 0.20    0.1986693          0.1986459       0.0000235 *
* 0.25    0.2474040          0.2474157      -0.0000118 *
* 0.30    0.2955202          0.2955218      -0.0000016 *
* 0.35    0.3428978          0.3428916       0.0000062 *
* 0.40    0.3894183          0.3894265      -0.0000081 *
* 0.45    0.4349655          0.4349655       0.0000000 *
* 0.50    0.4794255          0.4794204       0.0000051 *
* 0.55    0.5226872          0.5226894      -0.0000021 *
* 0.60    0.5646425          0.5646493      -0.0000068 *
* 0.65    0.6051864          0.6051821       0.0000043 *
* 0.70    0.6442177          0.6442141       0.0000035 *
* 0.75    0.6816388          0.6816405      -0.0000017 *
* ---------------------------------------------------- *
*                                                      *
*******************************************************}
program AkimaLaz;

{*******************************************************
*          Akima spline fitting subroutine             *
* ---------------------------------------------------- *
* The input table is X(i), Y(i), where Y(i) is the     *
* dependant variable. The interpolation point is xx,   *
* which is assumed to be in the interval of the table  *
* with at least one table value to the left, and three *
* to the right. The interpolated returned value is yy. *
* The function result is false when an error occored.  *
* It is also assumed that the X(i) are in ascending    *
* order.                                               *
*******************************************************}

function Interpol_Akima(X, Y: Array of Double; xx: Double; var yy: Double): Boolean;
var
  i, n: integer;
  a, aa, aaa, b, bb: Double;
  X0: Double;
  XM: array of Double = nil;
  Z: array of Double = nil;
begin
  Result := true;

  n := Length(X);
  if Length(Y) < n then
    n := Length(Y);

  { Special case xx=0 }
  if xx = 0.0 then
  begin
    yy := 0.0;
    exit;
  end;

  { Check to see if interpolation point is correct }
  if (xx < X[0]) or (xx >= X[n-4]) then
  begin
    Result := false;
    exit;
  end;

  X0 := 2.0*X[0] - X[1];

  { Calculate Akima coefficients, a and b }
  SetLength(XM, n+3);
  SetLength(Z, n);
  for i := 0 to n-1 do
    {Shift i to i+2}
    XM[i+2] := (Y[i+1]-Y[i]) / (X[i+1]-X[i]);
  XM[n+1] := 2.0*XM[n] - XM[n-1];
  XM[n+2] := 2.0*XM[n+1] - XM[n];
  XM[1] := 2.0*XM[2] - XM[3];
  XM[0] := 2.0*XM[1] - XM[2];

  for i := 0 to n-1 do
  begin
    a := abs(XM[i+3] - XM[i+2]);
    b := abs(XM[i+1] - XM[i]);
    if a + b = 0 then
    begin
      Z[i] := (XM[i+2] + XM[i+1]) / 2.0;
      break;
    end;
    Z[i] := (a*XM[i+1] + b*XM[i+2]) / (a + b);
  end;

  {Find relevant table interval}
  i := 0;
  while (xx > X[i]) do
    inc(i);
  dec(i);

  {Begin interpolation}
  if i = -1 then
  begin
    b := X[i+1] - X0;
    a := xx - X0;
  end else
  begin
    b := X[i+1] - X[i];
    a := xx - X[i];
  end;

  a := xx - X[i];
  aa := a*a;
  aaa := aa*a;
  bb := b*b;
  yy := Y[i] + Z[i]*a + (3.0*XM[i+2]-2.0*Z[i]-Z[i+1])*aa/b;
  yy := yy + (Z[i]+Z[i+1]-2.0*XM[i+2])*aaa/(bb);
end;

{main program}
var
  X: Array[0..13] of Double;
  Y: Array[0..13] of Double;
  xx, yy: Double;
  i: Integer;
begin
  WriteLn;
  WriteLn(' Akima spline fitting of SIN(X):');

  {Input sine table
  -----------------------------------------------------------------
   Sine table values from  Handbook of mathematical functions
   by M. Abramowitz and I.A. Stegun, NBS, june 1964
  ----------------------------------------------------------------}
  X[0]:=0.000; Y[0]:=0.00000000;
  X[1]:=0.125; Y[1]:=0.12467473;
  X[2]:=0.217; Y[2]:=0.21530095;
  X[3]:=0.299; Y[3]:=0.29456472;
  X[4]:=0.376; Y[4]:=0.36720285;
  X[5]:=0.450; Y[5]:=0.43496553;
  X[6]:=0.520; Y[6]:=0.49688014;
  X[7]:=0.589; Y[7]:=0.55552980;
  X[8]:=0.656; Y[8]:=0.60995199;
  X[9]:=0.721; Y[9]:=0.66013615;
  X[10]:=0.7853981634; Y[10]:=0.7071067812;
  X[11]:=0.849; Y[11]:=0.75062005;
  X[12]:=0.911; Y[12]:=0.79011709;
  X[13]:=0.972; Y[13]:=0.82601466;

  { print header }
  WriteLn;
  WriteLn(' X   SIN(X) HANDBOOK   AKIMA INTERPOLATION   ERROR  ');
  WriteLn('----------------------------------------------------');

  { Main loop }
  xx := 0.0;
  for i:=1 to 16 do
  begin
    if Interpol_Akima(X, Y, xx, yy) then
      WriteLn(xx:4:2, '    ', sin(xx):9:7, '          ', yy:9:7,'      ', (sin(xx)-yy):10:7)
    else
    begin
      WriteLn(xx:4:2,'    ', 'Error calculating the Akima interpolation');
      break;
    end;
    xx := xx + 0.05
  end;

  { Print footer }
  WriteLn('----------------------------------------------------');
  WriteLn;
  Write('Press ENTER to close...');
  ReadLn;
end.
{ End of file akima.pas}
