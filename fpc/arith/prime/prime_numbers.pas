{******************************************************************
*                         PRIME NUMBERS                           *
* --------------------------------------------------------------- *
* This small program tests if a given integer number is a prime   *
* number or not.                                                  *
* --------------------------------------------------------------- *
* Ref.: "Mathématiques par l'informatique individuelle (1)        *
*         par H. Lehning et D. Jakubowicz, Masson, Paris, 1982"   *
*         [BIBLI 06].                                             *
*                                                                 *
*                                  Pascal version By J-P Moreau.  *
*                                        (www.jpmoreau.fr)        *
* --------------------------------------------------------------- *
* Modernized version for FreePascal/Lazarus:                      *
* - Make cross-patform by removing WinCRT                         *
* - Remove goto instructions and labels                           *
* - Use mod operator to check divisibility --> clearer code.      *
*******************************************************************
See also program factors.pas}
program prime_numbers;

const
  EPS = 1e-6;

var
  d: Real;
  i, m, n: Integer;

procedure Pause;
begin
  WriteLn;
  Write('Press ENTER to close...');
  ReadLn;
end;

procedure NoPrime;
begin
  WriteLn(' Not a prime number.');
  Pause;
end;

begin
  WriteLn;
  Write(' Enter number to test: ');
  ReadLn(n);
  WriteLn;

  {Test if multiple of 2}
  {
  d := n - 2*INT(n/2);
  if abs(d)<eps then
  }
  if n mod 2 = 0 then
  begin
    NoPrime;
    exit;
  end;

  {Test if multiple of 3}
  {
  d := N - 3*INT(n/3);
  if abs(d)<eps then
  }
  if n mod 3 = 0 then
  begin
    NoPrime;
    exit;
  end;

  {Test if multiple of 6i-1 and 6i+1 from i=6 to sqrt(N)+1}
  i := 6;
  m := trunc(sqrt(n)) + 1;
  while i <= m do
  begin
  {
    d := n - (i-1)*INT(n/(i-1));
    if abs(d)<eps then
  }
    if n mod (i-1) = 0 then
    begin
      NoPrime;
      exit;
    end;
    {
    d := n - (i+1)*INT(n/(i+1));
    if abs(d)<eps then
    }
    if n mod (i+1) = 0 then
    begin
      NoPrime;
      exit;
    end;
    i := i + 6;
  end;
  WriteLn(' Prime number.');
  Pause;
end.

{end of file prime.pas}
