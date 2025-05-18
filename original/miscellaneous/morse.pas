{******************************************************************************
*                      Decode  a  text  in  morse                             *
* --------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                 *
* File morse.dat contains:                                                    *
*.--- . ...-.- ... ..- .. ... ...-.- -.-. --- -. - . -. - ...-.- -.. . ...-.- *
*...- --- ..- ... ...-.- ...- --- .. .-. -.-.- ...-.-                         *
*.--- . .- -. -....- .--. .. . .-. .-. . ...-.- -- --- .-. . .- ..- ...-.-    *
*...---                                                                       *
*                                                                             *
* Screen output:                                                              *
*.--- . ...-.- ... ..- .. ... ...-.- -.-. --- -. - . -. - ...-.- -.. . ...-.- *
* je suis content de                                                          *
*...- --- ..- ... ...-.- ...- --- .. .-. -.-.- ...-.-                         *
* vous voir.                                                                  *
*.--- . .- -. -....- .--. .. . .-. .-. . ...-.- -- --- .-. . .- ..- ...-.-    *
* jean-pierre moreau                                                          *
*...---                                                                       *
*                                                                             *
* Program terminated.                                                         *
*                                                                             *
*                                  TPW Release By Jean-Pierre Moreau, Paris.  *
*                                              (www.jpmoreau.fr)              *
******************************************************************************}
Program Decode_Morse;

Uses WinCrt;

Label 10;

type
     Str6 = string[6];

var
     morse: array[0..45] of Str6;
     alpha: array[0..45] of char;

     ligne, decodee: String;

     mot: Str6;

    fp1:TEXT;

    i, count, j, jj: integer;


BEGIN

{initialize table of morse signs}
  morse[0]:='.';
  morse[1]:='-';
  morse[2]:='..';
  morse[3]:='.-';
  morse[4]:='-.';
  morse[5]:='--';
  morse[6]:='...';
  morse[7]:='..-';
  morse[8]:='.-.';
  morse[9]:='.--';
  morse[10]:='-..';
  morse[11]:='--.';
  morse[12]:='-.-';
  morse[13]:='---';
  morse[14]:='....';
  morse[15]:='...-';
  morse[16]:='..-.';
  morse[17]:='.-..';
  morse[18]:='.-.-';
  morse[19]:='.--.';
  morse[20]:='..--';
  morse[21]:='.---';
  morse[22]:='-...';
  morse[23]:='-.-.';
  morse[24]:='-..-';
  morse[25]:='--..';
  morse[26]:='---.';
  morse[27]:='--.-';
  morse[28]:='-.--';
  morse[29]:='.----';
  morse[30]:='..---';
  morse[31]:='...--';
  morse[32]:='....-';
  morse[33]:='.....';
  morse[34]:='-....';
  morse[35]:='--...';
  morse[36]:='---..';
  morse[37]:='----.';
  morse[38]:='-----';
  morse[39]:='-.-.-';
  morse[40]:='--..--';
  morse[41]:='---...';
  morse[42]:='-....-';
  morse[43]:='--..--';
  morse[44]:='...-.-';
  morse[45]:='#13';

{initialize table of alphanumeric characters}
  alpha[0]:='e';
  alpha[1]:='t';
  alpha[2]:='i';
  alpha[3]:='a';
  alpha[4]:='n';
  alpha[5]:='m';
  alpha[6]:='s';
  alpha[7]:='u';
  alpha[8]:='r';
  alpha[9]:='w';
  alpha[10]:='d';
  alpha[11]:='g';
  alpha[12]:='k';
  alpha[13]:='o';
  alpha[14]:='h';
  alpha[15]:='v';
  alpha[16]:='f';
  alpha[17]:='l';
  alpha[18]:='*';
  alpha[19]:='p';
  alpha[20]:='*';
  alpha[21]:='j';
  alpha[22]:='b';
  alpha[23]:='c';
  alpha[24]:='x';
  alpha[25]:='z';
  alpha[26]:='*';
  alpha[27]:='q';
  alpha[28]:='y';
  alpha[29]:='1';
  alpha[30]:='2';
  alpha[31]:='3';
  alpha[32]:='4';
  alpha[33]:='5';
  alpha[34]:='6';
  alpha[35]:='7';
  alpha[36]:='8';
  alpha[37]:='9';
  alpha[38]:='0';
  alpha[39]:='.';
  alpha[40]:=',';
  alpha[41]:=':';
  alpha[42]:='-';
  alpha[43]:='?';
  alpha[44]:=' ';
  alpha[45]:=Chr(13);

  Assign(fp1,'morse.dat'); Reset(fp1);

  while Not eof(fp1) do
  begin
    readln(fp1,ligne);
    writeln(ligne);

    decodee:=''; mot:='';

    for i:=1 to length(ligne) do
      if (ligne[i]='.') or (ligne[i]='-') then
	mot:=mot + ligne[i]
      else
      begin
	if mot='...---' then goto 10; {exit program}
        for j:=0 to 45 do
	  if mot = morse[j] then
	    decodee:=decodee + alpha[j];
        mot:=''
      end;
    writeln(' ',decodee)
  end;

10: Close(fp1);
  writeln(' Program terminated.');
  writeln;

  ReadKey;
  DoneWinCrt

End.

{end of file morse1.cpp

 j   e         s   u  i   s          c    o  n  t e n  t         d  e
.--- . ...-.- ... ..- .. ... ...-.- -.-. --- -. - . -. - ...-.- -.. . ...-.-
  v   o   u   s          v    o   i  r  .
...- --- ..- ... ...-.- ...- --- .. .-. -.-.-
 j   e a  n     -    p   i  e  r   r  e        m   o   r  e a   u
.--- . .- -. -....- .--. .. . .-. .-. . ...-.- -- --- .-. . .- ..-
...---  }



