{********************************************
*     Elementary operations on fractions    *
* ----------------------------------------- *
* This program allows user performing       *
* operations on fractions such as:          *
*         3/4 x 8/9 - 4/5                   *
*                                           *
* Sample run:                               *
*                                           *
* 3 4                                       *
* *                                         *
* 8 9                                       *
* +                                         *
* -4 5                                      *
* =                                         *
* Result is -2/15                           *
* ----------------------------------------- *
* Ref.: "Mathématiques par l'informatique   *
*        individuelle (1) par H. Lehning    *
*        et D. Jakubowicz, Masson, Paris,   *
*        1982" [BIBLI 06].                  *
*                                           *
*           Pascal version by J-P Moreau.   *
*                 (www.jpmoreau.fr)         *
********************************************}
PROGRAM Fractions;
Uses Wincrt;

Label 10,20,fin;

Var   b,g,h,p,q,r,s,t,u : INTEGER;
      cs : CHAR;

      Procedure P300;   {substract}
      Begin
        r:=-r;
        p:=p*s+q*r;
        q:=q*s          {common denominator for all operations}
      End;

      Procedure P310;   {add}
      Begin
        p:=p*s+q*r;
        q:=q*s          {common denominator for all operations}
      End;

      Procedure P330;   {Divide}
      Begin
        b:=r;
        r:=s; s:=b;
        p:=p*r;
        q:=q*s          {common denominator for all operations}
      End;

      Procedure P360;   {Multiply}
      Begin
        p:=p*r;
        q:=q*s          {common denominator for all operations}
      End;

      Procedure P390;   {simplify result}
      Label 10,20;
      Begin
        t:=p;
        u:=q;
        10:b:=t-u*Round(t/u);
        if b=0 then goto 20;
        t:=u; u:=b;
        goto 10;
        20:p:=p DIV u; q:=q DIV u;
        if q<0 then begin q:=-q; p:=-p end
      End;

{main program}
BEGIN
  clrscr;
  writeln;
  {*** enter first fraction p/q (p=0 to exit) ***}
10:readln(p,q); if (p=0) or (q=0) then goto fin;
  {*** enter operator:
       = to view current result p/q
       s to save current result
       +,-,*,/ arithmetic operator
  ***}
20:readln(cs);
  if cs='=' then begin writeln(' Result is ',p,'/',q); goto 20 end; 
  if cs='s' then begin g:=p; h:=q; goto 10 end;
  {*** enter next fraction r/s (if s=0, take current stored result) ***}
  readln(r,s); if r=0 then goto fin;
  if s=0 then begin r:=g; s:=h end;
  {computing result}
  if cs='+' then P310;
  if cs='-' then P300;
  if cs='*' then P360;
  if cs='/' then P330;
  P390;
  goto 20;
fin: readkey;
  DoneWinCrt
END.


{end of file fraction.pas}