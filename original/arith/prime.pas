{********************************************
*              PRIME NUMBERS                *
* ----------------------------------------- *
* This small program tests if a given inte- *
* ger number is a prime number or not.      *
* ----------------------------------------- *
* Ref.: "Mathématiques par l'informatique   *
*        individuelle (1) par H. Lehning    *
*        et D. Jakubowicz, Masson, Paris,   *
*        1982" [BIBLI 06].                  *
*                                           *
*            Pascal version By J-P Moreau.  *
*                  (www.jpmoreau.fr)        *
*********************************************
See also program factors.pas}
PROGRAM Prime_numbers;
Uses WinCrt;

Label 100, fin;

Var   d,eps,i,m,n : REAL;

BEGIN
  eps:=1e-6;
  clrscr;
  writeln;
  write(' ? '); readln(n);
  writeln;
  {Test if multiple of 2}
  d:=n-2*INT(n/2);
  if abs(d)<eps then goto 100;
  {Test if multiple of 3}
  d:=N-3*INT(n/3);
  if abs(d)<eps then goto 100;
  {Test if multiple of 6i-1 and 6i+1 from i=6 to sqrt(N)+1}
  i:=6; m:=INT(SQRT(n))+1;
  While i<=m do
  begin
    d:=n-(i-1)*INT(n/(i-1));
    if abs(d)<eps then goto 100;
    d:=n-(i+1)*INT(n/(i+1));
    if abs(d)<eps then goto 100;
    i:=i+6
  end;
  writeln(' Prime number.'); goto fin;
100:writeln(' Not a prime number.');
fin:Readkey; DoneWinCrt
END.

{end of file prime.pas}