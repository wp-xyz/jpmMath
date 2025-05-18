(*****************************************************************************
*          This demo program encodes a text file by the letter               *
*          transposing method                                                *
* -------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                *
* File test.pas contains:                                                    *
*                                                                            *
* Program Relativistic_Mass;                                                 *
*                                                                            *
* Uses WinCrt1;                                                              *
*                                                                            *
* Var                                                                        *
*                                                                            *
* rest_mass, rela_mass,         { kg }                                       *
* voltage,                      { volt }                                     *
* speed,                        { m/s }                                      *
* e,                            { electron charge in coulomb }               *
* c: double;                    { speed of light in m/s }                    *
*                                                                            *
* Begin                                                                      *
*                                                                            *
* e:=1.602e-19; c:=2.9979e8; rest_mass:=9.109e-31;                           *
*                                                                            *
* writeln;                                                                   *
* write(' Give electron gun voltage in volts: ');                            *
*                                                                            *
* readln(voltage);                                                           *
*                                                                            *
* rela_mass := (voltage*e + rest_mass*c*c) / (c*c);                          *
*                                                                            *
* speed := c*sqrt(1.0 - Sqr(rest_mass/rela_mass));                           *
*                                                                            *
* writeln;                                                                   *
* writeln(' Relativistic mass (kg) and speed (m/s):');                       *
* writeln(' ', rela_mass, ' ', speed);                                       *
* writeln;                                                                   *
*                                                                            *
* ReadKey;                                                                   *
* DoneWinCrt                                                                 *
*                                                                            *
* End.                                                                       *
*                                                                            *
*     TRANSPOSITION CIPHER                                                   *
* Give name of input text file: test.pas                                     *
* Give name of output text file: test.txt                                    *
* Give flag (e) to encode (d) decode: e                                      *
* Give transposition distance: 12                                            *
*                                                                            *
* File test.txt now contains:                                                *
*                                                                            *
* tMvcstii_iaaPeo ramgRrlsi                                                  *
* C                                                                          *
* t1;r                                                                       *
* n                                                                          *
* Wss                                                                        *
* s                                                                          *
* U                                                                          *
* e; Vtem ss,ar_lsar                                                         *
*                                                                            *
*                                                                            *
* rea                                                                        *
*   kg{}                                                                     *
* _ a s, s m                oetagl,v  te}s                                   *
*                                                                            *
* p el v {      od            ,           e        ,                         *
* }/s m                                                                      *
* {   l    { e e            cib ooulcmn  tgoa chnrre}e        ; l            *
* u dc:  o                                                                   *
* b slefd oe pi           {g                                                 *
*                                                                            *
*                                                                            *
* negiB                                                                      *
*                                                                            *
*                                                                            *
*                                                                            *
* h  /n mist}                                                                *
* 19;2c:= .99- 2e6=1.:0 e7se:09.1=9s-a9_8s re;tem3ern                        *
*  ;wlit1r                                                                   *
*                                                                            *
*                                                                            *
* w;itc run gontveee'eGiv  (lot                                              *
* :                                                                          *
* ');                                                                        *
* s                                                                          *
* llvane ig to a e                                                           *
* ;                                                                          *
* )                                                                          *
* g t oe(dlnavrlr(eoetagl*v  e:asmas_ l=+*c)c/ ( *c)c seat_mssr*;:( r*sqct=1 * 
* e                                                                          *
* p  s                                                                       *
* e                                                                          *
* d._aaes/rslm_t0e-(Sqr r smw                                                *
* i;elnt                                                                     *
* r  a                                                                       *
* s                                                                          *
* );                                                                         *
* )                                                                          *
* s  lcttvisiia ew i(elnt'rRm /p(ed emssdaas)(kg  sn)ern,' '( let:r)         *
*  ;w'il                                                                     *
* p;ed)e                                                                     *
* s ,a m ss,a'_'  ;ReadKey                                                   *
* w                                                                          *
* i;elnt                                                                     *
* r                                                                          *
*                                                                            *
*                                                                            *
* .End                                                                       *
*                                                                            *
*  t CDineWon r                                                              *
*                                                                            *
*     TRANSPOSITION CIPHER                                                   *
* Give name of input text file: test.txt                                     *
* Give name of output text file: test1.pas                                   *
* Give flag (e) to encode (d) decode: d                                      *
* Give transposition distance: 12                                            *
*                                                                            *
* File test1.pas now contains: same text as test.pas.                        *
*                                                                            *
* -------------------------------------------------------------------------- *
* Reference: "Advanced Turbo C By Herbert Schildt, Borland-Osborne/          *
*             McGraw-Hill, 1987."                                            *
*                                                                            *
*                                        TPW Release By J-P Moreau, Paris.   *
*                                               (www.jpmoreau.fr)            *
*****************************************************************************)
Program Encode;

Uses WinCrt, Strings;

Label fin;

Const SIZE = 4096;

Var
    dist: Word;          {user-specified transposition distance}
    input:String[40];    {name of input text file (to encode/decode}
    output:String[40];   {name of output text file (to encode/decode}
    flag: Char;          { (e)ncode or (d)ecode }


Procedure code(input, output: string; dist: Word);
Var
    done: Boolean;
    temp: char;
    t: integer;
    s: Array[1..SIZE] of Char;

    fp1, fp2: TEXT;

Begin

  Assign(fp1, input); Reset(fp1);

  Assign(fp2, output); Rewrite(fp2);

  done:=False;

  Repeat
    for t:=1 to dist*2 do
    begin
      read(fp1,temp);
      s[t]:=temp;
      if EOF(fp1) then done:=True
    end;
    for t:=1 to dist do
    begin
      temp:=s[t];
      s[t]:=s[t+dist];
      s[t+dist]:=temp;
      Inc(t);
      temp:=s[t];
      s[t]:=s[dist*2-t];
      s[dist*2-t]:=temp
    end;
    for t:=1 to dist*2 do write(fp2,s[t])
  Until done;
  close(fp1);
  close(fp2)
End; {code}


Procedure decode(input, output: string; dist: Word);
Var
    done: Boolean;
    temp: char;
    t: integer;
    s: Array[1..SIZE] of Char;

    fp1, fp2: TEXT;

Begin

  Assign(fp1, input); Reset(fp1);

  Assign(fp2, output); Rewrite(fp2);

  done:=False;

  Repeat
    for t:=1 to dist*2 do
    begin
      read(fp1,temp);
      s[t]:=temp;
      if EOF(fp1) then done:=True
    end;
    for t:=1 to dist do
    begin
      Inc(t);
      temp:=s[t];
      s[t]:=s[dist*2-t];
      s[dist*2-t]:=temp;
      Dec(t);
      temp:=s[t];
      s[t]:=s[t+dist];
      s[t+dist]:=temp;
      Inc(t)
    end;
    for t:=1 to dist*2 do write(fp2, s[t])
  Until done;
  close(fp1);
  close(fp2);
End; {decode}


{main program}
BEGIN

  Writeln;
  Writeln('    TRANSPOSITION CIPHER');
  Writeln;
  Write(' Give name of input text file: ');
  Readln(input);

  Writeln;
  Write(' Give name of output text file: ');
  Readln(output);

  Writeln;
  Write(' Give flag (e) to encode (d) decode: ');
  Readln(flag);

  if (flag<>'e') and (flag<>'d') then
  begin
    Writeln(' Usage (example): toto.in toto.out e');
    goto fin
  end;

  Writeln;
  Write(' Give transposition distance: ');
  Readln(dist);
  Writeln;

  if flag='e' then
    code(input,output,dist)
  else
    decode(input,output,dist);

  Writeln(' Program terminated.'); 

fin: ReadKey;

  DoneWinCrt

END.

{end of file encode.pas}