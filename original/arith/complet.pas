{***********************************************************
*   This small program finds coincidences of letters in    *
*   two given words of up to 26 letters.                   *
* -------------------------------------------------------- *
* SAMPLE RUN:                                              *
*                                                          *
* Input first word : Bonjour                               *
* Input second word: Monsieur                              *
*                                                          *
* The letter o is common to the two words.                 *
* The letter n is common to the two words.                 *
* The letter u is common to the two words.                 *
* The letter r is common to the two words.                 *
*                                                          *
* -------------------------------------------------------- *
* Reference: "Exercices en Turbo Pascal By Claude Delannoy *
*             Eyrolles, 1997".                             *
***********************************************************}
Program Common_letters;

Uses WinCrt;

Const Lmax = 26;           {max. length of a word}

Var   mot1: String[Lmax];  {1st word}
      mot2: String[Lmax];  {2nd word}
      found: Boolean;      {2 letters identical if TRUE}
      i,j: integer;        {loop counters}


BEGIN

  {Read two words}
  writeln;
  write(' Input first word : '); readln(mot1);
  write(' Input second word: '); readln(mot2);
  writeln;

  {Compare the two words}
  For i:=1 to length(mot1) do
  begin
    found:=False;
    j:=1;
    While (Not found) and (j<=length(mot2)) do
    begin
      if mot1[i]=mot2[j] then
      begin
        writeln(' The letter ',mot1[i],' is common to the two words.');
        mot2[j]:=' ';     {to deactivate the found coincidence}
        found:=True
      end;
      Inc(j)
    end
  end;
  ReadKey;
  DoneWinCrt

END.

{end of file comlet.pas}