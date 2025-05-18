{***********************************************************
*     This small program displays a Pascal's triangle      *
*     up to 15 lines.                                      *
* -------------------------------------------------------- *
* SAMPLE RUN:                                              *
*                                                          *
* How many lines: 10                                       *
*                                                          *
* p      0    1    2    3    4    5    6    7    8    9    *
* n                                                        *
* -----------------------------------------------------    *
* 0 --   1                                                 *
* 1 --   1    1                                            *
* 2 --   1    2    1                                       *
* 3 --   1    3    3    1                                  *
* 4 --   1    4    6    4    1                             *
* 5 --   1    5   10   10    5    1                        *
* 6 --   1    6   15   20   15    6    1                   *
* 7 --   1    7   21   35   35   21    7    1              *
* 8 --   1    8   28   56   70   56   28    8    1         *
* 9 --   1    9   36   84  126  126   84   36    9    1    *
* -------------------------------------------------------- *
* Reference: "Exercices en Turbo Pascal By Claude Delannoy *
*             Eyrolles, 1997".                             *
***********************************************************}

Program Paslcal_Triangle;
Uses WinCrt1;

Const  Nmax = 15;                     {max. number of lines}

Var
     t  : array[0..Nmax-1] of Integer;  {one line of triangle}
     nl : integer;                      {desired number of lines}
     i,j: integer;                      {current line/column}


BEGIN

  {read number of lines & display caption}
  writeln;
  write(' How many lines: '); readln(nl);
  if nl>Nmax then nl:=Nmax;
  writeln;
  writeln;
  write(' p   ');
  for i:=0 to nl-1 do write(i:5);
  writeln;
  writeln(' n');
  for i:=0 to nl do write('-----');
  writeln;

  {create and display each line}
  for i:=0 to nl-1 do
  begin
    t[i]:=1;
    for j:=i-1 downto 1 do t[j]:=t[j]+t[j-1];
    write(i:2,' --');
    for j:=0 to i do write(t[j]:5);
    writeln
  end;

  for i:=0 to nl do write('-----');
  writeln;

  ReadKey;
  DoneWinCrt    

END.

{end of file pastri.pas}