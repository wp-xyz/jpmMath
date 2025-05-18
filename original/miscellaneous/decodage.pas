(*************************************************************************
*                  Decode an encoded ASCII Test                          *
* ---------------------------------------------------------------------- * 
* The ASCII text must have been encoded by codage.exe that generates two *
* output files:                                                          *
* 1. code.lst: a random encoding table                                   *
* 2. codage.lst: the encoded text.                                       *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
* File codage.lst contains:                                              *
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
* File code.lst contains:                                                *
*                                                                        *
* t!="1#h$y%X&''<(Q)?*W+7,g-Y.s/c061;2K3/425F6i758E9@:.;}<L=!>0?w@(ATBHC *
* uDpE[FZGRHVInJkK#LfMaNGOJP+QvR]S\TAUrVIWeXOY>Z4["\$]j^B_l`8amb{c,dSe|f *
* %gUhxi*jMkDl_mzndoop9q&r`s3t-uNv~wqxCy:z^{)|P}b~                       *
*                                                                        *
* (! is coded t, " is coded =, # is coded 1, etc.)                       *
*                                                                        *
* The output file decodage.lst contains;                                 *
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
*                                    TPW Release By J-P Moreau, Paris.   *
*                                           (www.jpmoreau.fr)            *
*************************************************************************)
Program Decodage;

Uses WinCrt;

Const MINCAR=33;     {ASCII characters from 33 to 126 only} 
      MAXCAR=126;
      MAXSTR=80;     {lines of length 80} 

Type  Tab = Array[MINCAR..MAXCAR] of char;

Var
      A, Z: Tab;
      alea,fait,i,j: integer;  
      ligne: String[MAXSTR];
      decodee: String[MAXSTR];

      fp1,fp2,fp3: TEXT;

BEGIN

  Assign(fp1,'codage.lst'); Reset(fp1);
  Assign(fp2,'decodage.lst'); Rewrite(fp2);
  Assign(fp3,'code.lst'); Reset(fp3);

  for i:=MINCAR to MAXCAR do  A[i]:=Chr(i);
         
  {read random encoding table}   
  writeln(' Random encoding table:');
  Readln(fp3, ligne);  
  Writeln(ligne);
  for i:=1 to 35 do Z[MINCAR+i-1] := ligne[2*(i-1)+1]; 
  Readln(fp3, ligne);
  Writeln(ligne);
  for i:=1 to 35 do Z[MINCAR+i+34] := ligne[2*(i-1)+1];    
  Readln(fp3, ligne);
  Writeln(ligne);
  for i:=1 to 24 do Z[MINCAR+i+69] := ligne[2*(i-1)+1];
  close(fp3);
  
  {decode line by line and save}  
  while Not(eof(fp1)) do
  begin
    readln(fp1, ligne);
    decodee:=ligne;
    for i:=1 to MAXSTR do
      for j:=MINCAR to MAXCAR do
	if ligne[i]=Z[j] then decodee[i]:=A[j];    
    writeln(fp2, decodee)
  end;
  close(fp1);                     
  writeln(fp2);
  close(fp2);
  writeln;
  writeln(' Program terminated. Results in file decodage.lst');
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file decodage.pas}
