{********************************************
*    Factorization of an integer number     *
* ----------------------------------------- *
* Sample run:                               *
*                                           *
* ? 394616                                  *
*                                           *
* 2 2 2 107 461                             *
*                                           *
* ----------------------------------------- *
* Ref.: "Mathématiques par l'informatique   *
*        individuelle (1) par H. Lehning    *
*        et D. Jakubowicz, Masson, Paris,   *
*        1982" [BIBLI 06].                  *
*                                           *
*            Pascal version by J-P Moreau.  *
*                  (www.jpmoreau.fr)        *
********************************************}
PROGRAM Factors;
Uses WinCrt;

Label 50,100,150,200;

VAR   d,eps,m,n,i : REAL;


BEGIN
  eps:=1e-6;
  {Enter integer number to be factorized}
  clrscr;
  writeln;
  write(' ? '); readln(n);
  writeln;
  write(' ');
  {Test if multiple of 2}
50:D:=N-2*round(N/2);
  if abs(D)<eps then begin write('2 '); N:=N/2; goto 50 end;
  {Test if multiple of 3}
100:D:=N-3*round(N/3);
  if abs(D)<eps then begin write('3 '); N:=N/3; goto 100 end;
  {Test of divisors 6i-1 and 6i+1 up to sqrt(N)
   Prime numbers are of the form 6i-1 or 6i+1 }
  m:=round(sqrt(n))+1;
  i:=6;
  while i<m+1 do
  begin
150:D:=N-(i-1)*Round(N/(i-1));
    if abs(D)<eps then begin write(i-1:4:0,' '); N:=N/(i-1); goto 150 end;
200:D:=N-(i+1)*Round(N/(i+1));
    if abs(D)<eps then begin write(i+1:4:0,' '); N:=N/(i+1); goto 200 end;
    i:=i+6
  end;
  if N>1 then writeln(N:4:0);
  ReadKey; DoneWinCrt
END.

{end of file factors.pas}