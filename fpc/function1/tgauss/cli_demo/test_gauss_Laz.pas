{********************************************************
* Integration by Gauss method of a real function F=F(X) *
* or F=F(X,Y) or F=F(X,Y,Z). The integral is calculated *
* by using from 2 to 10 Gauss points.                   *
* ----------------------------------------------------- *
* Ref.: "Mécanique des vibrations linéaires By          *
*        M. Lalanne, P. Berthier and J. Der Hagopian,   *
*        Masson, Paris, 1980" [BIBLI 16].               *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* (Integrate F=sin(x) from x=0 to x=1).                 *
*                                                       *
* INTEGRATION OF A REAL FUNCTION BY GAUSS               *
*        F(X), F(X,Y) or F(X,Y,Z)                       *
*                                                       *
* Number of variables (1 to 3): 1                       *
*                                                       *
* How many Gauss points (2 to 10): 4                    *
*                                                       *
* Minimum value of X: 0                                 *
* Maximum value of X: 1                                 *
*                                                       *
* Value of integral =  4.59697693864200E-0001           *
*                                                       *
*                    TPW Version By J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)           *
********************************************************}
Program test_gauss_Laz;

{$mode objFPC}{$H+}

uses
  jpmTypes, jpmIntegration;

{ Define here functions F1, F2, F3 of 1,2,3 variables}
Function F1(x: Float): Float;
begin
  Result := sin(X);
end;

Function F2(x, y: Float): Float;
begin
  Result := x * y; //sin(x) * cos(x);
end;

Function F3(x,y,z: Float): Float;
begin
  Result := x * y * z; //sin(x) * cos(x) * x;
end;

var
  nVar, n, n1, n2: Integer;
  x1, x2, y1, y2, z1, z2: Float;
  Integral: Float;
begin
  WriteLn;
  WriteLn(' INTEGRATION OF A REAL FUNCTION BY GAUSS');
  WriteLn('        F(X), F(X,Y) or F(X,Y,Z)');
  WriteLn;
  Repeat
    Write(' Number of variables (1 to 3): '); ReadLn(nVar);
    if (nVar > 3) or (nVar < 1) then
      WriteLn('  --> Only values 1, 2 or 3 allowed.');
  until (nVar in [1, 2, 3]);

  WriteLn;
  Write(' Minimum value of X: '); ReadLn(x1);
  Write(' Maximum value of X: '); ReadLn(x2);
  WriteLn;
  iF nVar > 1 then
  begin
    Write(' Minimum value of Y: '); ReadLn(y1);
    Write(' Maximum value of Y: '); ReadLn(y2);
    Writeln;
  end;
  if nVar > 2 then
  begin
    Write(' Minimum value of Z: '); ReadLn(z1);
    Write(' Maximum value of Z: '); ReadLn(z2);
  end;

  n1 := 2;
  n2 := 10;
  WriteLn(' The number of Gauss points (polynomial order) is varied between ', n1, ' and ', n2);
  for n := n1 to n2 do
  begin
    case nVar of
      1: Integral := GaussIntegral(@F1, n, x1, x2);
      2: Integral := GaussIntegral(@F2, n, x1, x2, y1, y2);
      3: Integral := GaussIntegral(@F3, n, x1, x2, y1, y2, z1, z2);
    end;
    WriteLn( ' ', n:2, ' Gauss points --> value of integral = ', Integral);
  end;

  WriteLn;
  Write('Press ENTER to close...');
  ReadLn;
end.

{end of file tgauss.pas}
