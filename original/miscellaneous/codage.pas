(*************************************************************************
*          Encoding an ASCII text using a random encoding table          *
* ---------------------------------------------------------------------- *
* The text to encode is in file codage.dat.                              *
* The random encoding table is in file code.lst.                         *
* The output encoded text is in file codage.lst.                         *
*                                                                        *
* The text can be decoded by  decodage.exe that must have access to      *
* code.lst and codage.lst files.                                         *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
* File codage.dat contains:                                              *
*                                                                        *
* Program ASCII;                                                         *
* Uses WincrtMy;                                                         *
*                                                                        *
* VAR                                                                    *
*     i: word;                                                           *
*                                                                        *
* BEGIN                                                                  *
*   for i:=160 to 255 do                                                 *
*     write(i:4,' ',chr(i));                                             *
*   Readkey;                                                             *
*   DoneWinCrt                                                           *
* END.                                                                   *
*                                                                        *
* The output file code.lst contains:                                     *
*                                                                        *
* t!="1#h$y%X&''<(Q)?*W+7,g-Y.s/c061;2K3/425F6i758E9@:.;}<L=!>0?w@(ATBHC *
* uDpE[FZGRHVInJkK#LfMaNGOJP+QvR]S\TAUrVIWeXOY>Z4["\$]j^B_l`8amb{c,dSe|f *
* %gUhxi*jMkDl_mzndoop9q&r`s3t-uNv~wqxCy:z^{)|P}b~                       *
*                                                                        *
* (! is coded t, " is coded =, # is coded 1, etc.)                       *
*                                                                        *
* The output file codage.lst contains;                                   *
*                                                                        *
* J&d%&8_ (]HVV.                                                         *
* A`S` Ixz{&3fC.                                                         *
*                                                                        *
* r(v                                                                    *
*     x@ ~d&,.                                                           *
*                                                                        *
* TpZVa                                                                  *
*   |d& x@L6Fc 3d ;22 ,d                                                 *
*     ~&x3S<x@/7' '7{U&<xQQ.                                             *
*   vS8,MSC.                                                             *
*   udzSIxzH&3                                                           *
* pauY                                                                   *
*                                                                        *
*                                    TPW Release By J-P Moreau, Paris.   *
*                                           (www.jpmoreau.fr)            *
*************************************************************************)
Program Codage;

Uses WinCrt;

Const MINCAR=33;     {ASCII characters from 33 to 126 only} 
      MAXCAR=126;
      MAXSTR=80;     {lines of length 80} 

Type  Tab = Array[MINCAR..MAXCAR] of char;

Var
      A, Z: Tab;
      alea,fait,i,j: integer;  
      ligne: String[MAXSTR];
      codee: String[MAXSTR];

      fp1,fp2,fp3: TEXT;

BEGIN

  Randomize;     {initialize random generator}

  Assign(fp1,'codage.dat'); Reset(fp1);
  Assign(fp2,'codage.lst'); Rewrite(fp2);
  Assign(fp3,'code.lst'); Rewrite(fp3);
    
  for i:=MINCAR to MAXCAR do  A[i]:=Chr(i);  
  
  {generate a random encoding table}        
  alea:=Random(MAXCAR-MINCAR+1);

  Z[MINCAR]:=Chr(MINCAR+alea);  
  for i:=MINCAR+1 to MAXCAR do 
    repeat
      alea:=Random(MAXCAR-MINCAR+1);
      Z[i]:=Chr(MINCAR+alea);
      fait:=1;
      for j:=MINCAR to i-1 do
        if Z[j]=Z[i] then fait:=0
    Until fait=1;

  {encode input text line by line and save}
  while Not(eof(fp1)) do
  begin
    Readln(fp1,ligne);
    codee:=ligne;
    for i:=1 to MAXSTR do
      for j:=MINCAR to MAXCAR do
      	if ligne[i]=A[j] then codee[i]:=Z[j];
    writeln(fp2, codee)
  end;
  close(fp1);                     
  writeln(fp2);
  close(fp2);

  {Print random encoding table in file code.lst} 
  for i:=MINCAR to MAXCAR do
  begin
    write(fp3, Z[i],A[i]);
    if ((i-MINCAR+1) mod 35)=0 then writeln(fp3)
  end;  
  writeln(fp3);  
  close(fp3);
  writeln;
  writeln(' Program terminated. Results in file codage.lst');
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file codage.pas}
