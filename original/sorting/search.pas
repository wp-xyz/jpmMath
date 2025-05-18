{*****************************************************
*    Program to demonstrate the use of searching     *
*    subroutines of unit FSearch.                    *
* -------------------------------------------------- *
* SAMPLE RUN:                                        *
* (An ordered list of 16 names is given in text file *
*  'search.dat').                                    *
*                                                    *
* # of items in list = 16                            *
* What name would you like to look for? David        *
* Linear serach for (f)irst, (a)ll, or (b)inary? f   *
* David found at position 5.                         *
* Try again (y/n)? y                                 *
* What name would you like to look for? Grace        *
* Linear serach for (f)irst, (a)ll, or (b)inary? a   *
* Grace at position 9.                               *
* Grace at position 10.                              *
* Grace at position 11.                              *
* 3 occurence(s) found.                              *
* Try again (y/n)? y                                 *
* What name would you like to look for? Grace        *
* Linear serach for (f)irst, (a)ll, or (b)inary? b   *
* Grace found at position 10.                        *
* Try again (y/n)? n                                 *
*                                                    *
* -------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90, By David R.*
*        Brooks, Springer Verlag, 1997".             *
*                                                    *
*                 TPW Release By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)           *
*****************************************************}
PROGRAM Search;

Uses WinCrt, FSearch;

Var
    YesNo, choice: Char;
    size, where, how_many: Integer;


  Procedure GetList(var size:Integer);
  Var fp:TEXT;
  Begin
    Assign(fp,'search.dat'); Reset(fp);
    size:=1;
    While Not EoF(fp) do
    begin
      Readln(fp, a[size]);
      Inc(size)
    end;
    size:=size-1;
    Close(fp);
    writeln;
    writeln(' # of items in list = ',  size)

  End;     


{main program}
BEGIN

  ClrScr;
  GetList(size);   {read names from input file}
  YesNo:='y';

  While YesNo = 'y' do
  begin
    write(' What name would you like to look for? '); readln(target);
    write(' Linear serach for (f)irst, (a)ll, or (b)inary? '); readln(choice);
	
    CASE choice of
      'a': begin
	     FindAll(size,how_many);
	     IF how_many>0 THEN
	       writeln(' ',how_many,' occurence(s) found.')
             ELSE
	       writeln(' ',target,' not found.')
           end;
      'b': begin
	     Binary(1,size,where);
	     IF where>0 THEN
	       writeln(' ',target,' found at position ', where,'.')
             ELSE
	       writeln(' ',target,' not found.')
           end;
      'f': begin 
	     FindFirst(size,where);
	     IF where>0 THEN
	       writeln(' ',target,' found at position ', where,'.')
             ELSE
	       writeln(' ',target,' not found.')
           end
    End;

    write(' Try again (y/n)? '); readln(YesNo)

  end;

  DoneWincrt

END.

{end of file search.pas}