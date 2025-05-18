{*************************************************************************
*     A PChar-Based Calculator with big integers using Unit PCalcFun     *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
* >2 ^ 100                                                               *
*                                                                        *
* POWER: 1,267,650,600,228,229,401,496,703,205,376                       *
*                                                                        *
* >200 !                                                                 *
*                                                                        *
* FACTORIAL: 30,414,093,201,713,378,043,612,608,166,064,768,844,377,     *
*            641,568,960,512,000,000,000,000                             *
*                                                                        *
* ---------------------------------------------------------------------- *
* * Ref.: "Turbo Pascal for Windows - Techniques and Utilities By        *
*          Neil J. Rubenking, PC Magazine & ZD Press, 1992."             *
*                                                                        *
*************************************************************************}
Program PCalculator;

Uses WinCrt1, WinTypes, Strings, PCalcFun;

Const Max = 5000;   {Maximum length of integer result!}

Var
    op, opRand, rem: Array[0..Max] of Char;
    operation: Char;
    finished: Boolean;
    Status: word;


    Function HeapFunc(size:word): integer; FAR;
    Begin
      HeapFunc:=1
    End;

    Function GetParams: Boolean;
    { Purpose: returns True if parameters are correctly passed
      on the input line -- and assigns them to the correct
      variables if so.}
    Var
        S, P, Pnm1, Pnm2: PChar;
        Len,posn: Byte;
    Begin
      GetMem(S, 2*Max+4);
      GetParams:=False;
      Write('>');
      Len:=ReadBuf(S, 2*Max+3);
      if Len=2 then
      begin
        Finished:=True;
        FreeMem(S, 2*Max+4);
        exit
      end;
      S[Len-2]:=#0;
      P:=StrScan(S,' ');
      if P<>Nil then
      begin
        posn:=P-S;
        operation:=S[Succ(posn)];
        pnm1:=StrScan(S,'#');
        if pnm1=Nil then
        begin
          StrLCopy(op, S, posn);
          StrCopy(opRand, P+3)
        end
        else if pnm1 > P then
        begin
          StrCopy(opRand, StripComma(op));
          StrLCopy(op, S, posn)
        end
        else
        begin
          StripComma(op);
          pnm2:=StrRScan(S,'#');
          if pnm2 > P then StrCopy(opRand, op)
          else StrCopy(opRand, S+posn+3)
        end;
        FreeMem(S, 2*Max+3)
      end
      else
      begin
        Writeln(' Enter n1 op n2, where op is + - * / or ^');
        Writeln(' or    n1 ! for factorial...');
        FreeMem(S, 2*Max+3);
        exit
      end;

      if Not AllNums(op) then
      begin
        writeln(' op is not a valid number!');
        exit
      end;

      Case operation of
        '!' : ;
        '+','-','*','/','^' :
        begin
          if opRand[0]=#0 then
          begin
            writeln('  The operation ',operation,' requires a second operand.');
            exit
          end;
          if Not AllNums(opRand) then
          begin  
            writeln(' opRand is not a valid number!');
            exit
          end
        end
        else
        begin
          writeln(' Valid operations are: + - * / ^ and !');
          exit
        end
      end;
      GetParams:=True
    End; {GetParams}                          

 {main program}
 BEGIN
   HeapError:=@HeapFunc;
   Repeat
     if GetParams then
     begin
       Case operation of
         '+':
         begin
           write('  SUM: ');
           add(op,opRand,Max,Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max))
         end;
         '-':
         begin
           write('  DIFF: ');
           sub(op,opRand,Max,Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max))
         end;
         '*':
         begin
           write('  PRODUCT: ');
           prod(op,opRand,Max,Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max))
         end;
         '/':
         begin
           write('  QUOTIENT: ');
           divide(op,opRand,rem,Max,Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max));
           write('Remainder: ');
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(rem,Max));
         end;
         '!':
         begin
            write('  FACTORIAL: ');
           Fact(op, Max, Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max));
         end;
         '^':
         begin
           write('  POWER: ');
           Power(op,opRand,Max,Status);
           if Status <> cer_Ok then writeln('*ERROR*')
           else writeln(AddComma(op,Max))
         end;
       end
     end
   Until finished;
   ReadKey;
   DoneWinCrt


 END.

 {end of file Pcalc.pas}






