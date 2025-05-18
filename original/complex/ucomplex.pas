{*********************************************
*            COMPLEX  CALCULATOR             *
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal By    *
*        M. Ducamp and A. Reverchon (vol 2), *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* SAMPLE RUN:                                *
*                                            *
* Calculate: sqrt((5-i)*Ln(3+2i))            *
*                                            *
* Enter operation: xa    Real part: 5        *
*                        Imaginary part: -1  *
* Enter operation: sto0                      *
* Enter operation: xa    Real part: 3        *
*                        Imaginary part: 2   *
* Enter operation: sto1                      *
* Enter operation: ln                        *
* Enter operation: sto2                      *
* Enter operation: *                         *
* Enter operation: sto3                      *
* Enter operation: sqr                       *
* Enter operation: sto4                      *
*                                            *
* (result is: 2.6640487 + i  0.3110939)      *
*                                            *
*          TPW Version By J-P Moreau, Paris. *
*                  (www.jpmoreau.fr)         *
*********************************************}
{
 The screen should look like that:

        COMPLEX NUMBERS CALCULATOR IN INVERSE POLISH NOTATION
 -------------------------------------------------------------------
 |  |        ALGEBRAIC FORM        |          POLAR FORM           |
 | T|  0.0000000                   |                               |
 | Z|  0.0000000                   |                               |
 | Y|  0.0000000                   |                               |
 | X|  2.6640487 + i  0.3110939    |  2.6821511  EXP+i  0.1162483  |
 -------------------------------------------------------------------
  Enter operation:
 -------------------------------------------------------------------
 |M0|  5.0000000 - i  1.0000000    |  5.0990195  EXP-i  0.1973956  |
 |M1|  3.0000000 + i  2.0000000    |  3.6055513  EXP+i  0.5880026  |
 |M2|  1.2824747 + i  0.5880026    |  1.4108467  EXP+i  0.4298922  |
 |M3|  7.0003760 + i  1.6575383    |  7.1939348  EXP+i  0.2324967  |
 |M4|  2.6640487 + i  0.3110939    |  2.6821511  EXP+i  0.1162483  |
 |M5|  0.0000000                   |                               |
 |M6|  0.0000000                   |                               |
 |M7|  0.0000000                   |                               |
 |M8|  0.0000000                   |                               |
 |M9|  0.0000000                   |                               |
 -------------------------------------------------------------------
    +    -     *    /    1/X    X^2    SQR     ^     EXP    LN
    SIN  COS   TAN  SH   CH     TH    ARCSIN ARCCOS ARCTAN ARGCH
   ARGSH ARGTH XA   XP   CLX    CLS    CLM    STOx   RCLx   Q
 -------------------------------------------------------------------

 Command explanations:
    XA:  enter a complex number x+iy, by giving x and y.
    XP:  enter a complex number r*exp(i*t), by giving r and t.
    CLX: stack register X is set to 0 and stack is shifted downwards.
         This allows to correct an input mistake.
    CLS: All four stack registers are set to 0.
    STOx: x from 0 to 9. stores X content in memory number x.
    RCLx: x from 0 to 9, recalls content of memory x in X register.
          The stack is previously shifted upwards.
    CLM: All ten memories are set to 0.
    Q:   Quit program.

}
PROGRAM Complex_calculator;
Uses WinCrt,WinProcs,Strings,Complex1;

Const  NBSTACK = 4;
       NBMEM   = 9;

{The COMPLEX Record is defined in complex1.pas}

Var    stack : Array[1..NBSTACK] of COMPLEX;
       memory: Array[0..NBMEM] of COMPLEX;
       oper  : String;
       i,j,k : integer;


  Procedure DownStack(n:integer);
  var i: integer;
  Begin
    if n<1 then n:=1;
    for i:=n to pred(NBSTACK) do stack[i]:=stack[i+1];
    fillchar(stack[NBSTACK], sizeof(stack[NBSTACK]), 0)
  End;

  Procedure UpStack;
  var i: integer;
  Begin
    for i:=NBSTACK downto 2 do stack[i]:=stack[i-1];
    fillchar(stack[1], sizeof(stack[1]), 0)
  End;

  {write a complex number (algebraic form + polar form) }
  Procedure DisplayNumber(n:COMPLEX; l:integer);
  var i: integer;
  Begin
    gotoxy(10,l); write(n.x:10:7);
    gotoxy(20,l);
    if n.y<>0 then
    begin
      if n.y>0 then write(' +') else write(' -');
      write(' i ',abs(n.y):10:7)
    end
    else write(' ':17);
    gotoxy(41,l);
    if n.r<>0 then write(n.r:10:7)
              else write(' ':10);
    gotoxy(51,l);
    if n.t<>0 then
    begin
      write(' EXP');
      if n.t>0 then write('+') else write('-');
      write('i ',abs(n.t):10:7)
    end
    else write(' ':17)
  End;

  Procedure InitDisplay;
  var i: integer;
  Begin
    Clrscr;
    gotoxy(12,1);
    write('COMPLEX NUMBERS CALCULATOR IN INVERSE POLISH NOTATION');
    gotoxy(5,2); for i:=1 to 67 do write('-');
    gotoxy(5,3);
    write('|  |        ALGEBRAIC FORM        ');
    write('|          POLAR FORM           |');
    gotoxy(5,4); write('| T|',' ':30,'|',' ':31,'|');
    gotoxy(5,5); write('| Z|',' ':30,'|',' ':31,'|');
    gotoxy(5,6); write('| Y|',' ':30,'|',' ':31,'|');
    gotoxy(5,7); write('| X|',' ':30,'|',' ':31,'|');
    for i:=1 to 4 do DisplayNumber(stack[1],8-i);
    gotoxy(5,8); for i:=1 to 67 do write('-');
    gotoxy(5,10); for i:=1 to 67 do write('-');
    for i:=0 to 9 do
    begin
      gotoxy(5,11+i);
      write('|M',i:1,'|',' ':30,'|',' ':31,'|')
    end;
    for i:=0 to 9 do DisplayNumber(memory[i],11+i);
    gotoxy(5,21); for i:=1 to 67 do write('-');
    gotoxy(5,22); write('  +    -     *    /    1/X    X^2    SQR     ^     EXP    LN  ');
    gotoxy(5,23); write('  SIN  COS   TAN  SH   CH     TH    ARCSIN ARCCOS ARCTAN ARGCH');
    gotoxy(5,24); write(' ARGSH ARGTH XA   XP   CLX    CLS    CLM    STOx   RCLx    Q  ');
    gotoxy(5,25); for i:=1 to 67 do write('-')
  End;

  {**********************************************************
  * Open a CRT window with title with possibility of text   *
  * and/or graph with the condition to uses WinCrtMy unit   *
  * instead of WinCrt (else CrtWindow is not visible).      *
  * Here the graphic possibility is not used.               *
  **********************************************************}
  PROCEDURE WinCrtInit(Name:PChar);
  BEGIN
    WindowOrg.X:= 70;    {upper left corner position}
    WindowOrg.Y:=100;    {and sizes of CRT window   }
    WindowSize.X:=620;
    WindowSize.Y:=425;
    StrCopy(WindowTitle,Name);  {window title}
    InitWinCrt;                {call standard Borland procedure}
  END;

  {execute command calling elementary routines of unit complex1.pas}
  Procedure Compute(oper: STRING);
  Var inter: COMPLEX;
  Begin
    if oper='XA' then
    begin
      UpStack;
      gotoxy(40,9); clreol;
      write('Real part: '); readln(stack[1].x);
      gotoxy(40,9); clreol;
      write('Imaginary part: '); readln(stack[1].y);
      RectPol(stack[1]);
      exit
    end;
    if oper='XP' then
    begin
      UpStack;
      gotoxy(40,9); clreol;
      write('Modulus: '); readln(stack[1].r);
      gotoxy(40,9); clreol;
      write('Argument: '); readln(stack[1].t);
      PolRect(stack[1]);
      exit
    end;
    if oper='CLX' then
    begin
      DownStack(1);
      exit
    end;
    if oper='CLS' then
    begin
      fillchar(stack,sizeof(stack),0);
      exit
    end;
    if oper='CLM' then
    begin
      fillchar(memory,sizeof(memory),0);
      InitDisplay;
      exit
    end;
    if copy(oper,1,3)='STO' then
    begin
      val(copy(oper,4,1),i,j);
      if j<>0 then MessageBeep(0)
      else begin
        memory[i]:=stack[1];
        DisplayNumber(memory[i],11+i)
      end;
      exit
    end;
    if copy(oper,1,3)='RCL' then
    begin
      val(copy(oper,4,1),i,j);
      if j<>0 then MessageBeep(0)
      else begin
        UpStack;
        stack[1]:=memory[i];
        DisplayNumber(memory[i],11+i)
      end;
      exit
    end;
    if oper='+' then
      if ZSum(stack[2],stack[1],stack[1]) then
        begin DownStack(2); exit end;
    if oper='-' then
      if ZMinus(stack[2],stack[1],stack[1]) then
        begin DownStack(2); exit end;
    if oper='*' then
      if ZMult(stack[2],stack[1],stack[1]) then
        begin DownStack(2); exit end;
    if oper='/' then
      if ZDiv(stack[2],stack[1],stack[1]) then
        begin DownStack(2); exit end;
    if oper='^' then
      if ZPower(stack[2],stack[1],stack[1]) then
        begin DownStack(2); exit end;
    if oper='EXP' then if ZExp(stack[1],stack[1]) then exit;
    if oper='LN'  then if ZLn(stack[1],stack[1]) then exit;
    if oper='1/X' then if ZInv(stack[1],stack[1]) then exit;
    if oper='X^2' then if ZSqr(stack[1],stack[1]) then exit;
    if oper='SQR' then if ZSqrt(stack[1],stack[1]) then exit;
    if oper='SIN' then if ZSin(stack[1],stack[1]) then exit;
    if oper='COS' then if ZCos(stack[1],stack[1]) then exit;
    if oper='TAN' then if ZTan(stack[1],stack[1]) then exit;
    if oper='SH'  then if ZSh(stack[1],stack[1]) then exit;
    if oper='CH'  then if ZCh(stack[1],stack[1]) then exit;
    if oper='TH'  then if ZTh(stack[1],stack[1]) then exit;
    if oper='ARCSIN' then if ZArcsin(stack[1],stack[1]) then exit;
    if oper='ARCCOS' then if ZArccos(stack[1],stack[1]) then exit;
    if oper='ARCTAN' then if ZArctan(stack[1],stack[1]) then exit;
    if oper='ARGSH'  then if ZArgsh(stack[1],stack[1]) then exit;
    if oper='ARGCH'  then if ZArgch(stack[1],stack[1]) then exit;
    if oper='ARGTH'  then if ZArgth(stack[1],stack[1]) then exit;
    MessageBeep(0)
  End;


{main program}
BEGIN
  WinCrtInit(' COMPLEX CALCULATOR');       {open application main window}
  fillchar(stack,sizeof(stack), 0);              {fill stack with zeroes}
  fillchar(memory,sizeof(memory), 0);           {fill memory with zeroes}
  InitDisplay;                                          {initial display}
  repeat                                                      {main loop}
    gotoxy(6,9); write('Enter operation: ');
    clreol; readln(oper);                                 {enter command}
    for i:=1 to length(oper) do              {put oper in upcase letters}
      oper[i]:=Upcase(oper[i]);
    if oper<>'Q' then
    begin
      Compute(oper);                                    {execute command}
      for i:=1 to 4 do DisplayNumber(stack[i],8-i)       {display result}
    end
  until oper='Q';
  DoneWinCrt                                          {close application}
END.

{end of file ucomplex.pas}