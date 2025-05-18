 {***************************************
 * Errors handling unit used by Ef1.pas *
 ***************************************} 
 UNIT ERRORHAN;

 INTERFACE
 Uses WinCrt1, Errors;

 Procedure ErrorHandler(ErrorNumber,ErrorAddress:Integer);

 IMPLEMENTATION

  Type Hex_Type = String[2];

  FUNCTION  Byte_to_Hex(What:Byte):Hex_Type;
  Var I:Byte;
      Rslt:Hex_Type;
      Rslt_Len:Byte absolute Rslt;
  Begin
      Rslt_Len := 2;
      For I:=2 downto 1 do
          Begin
          If (What and 15) < 10 then
             Rslt[I]:=Char((What and 15)+$30)
          else
             Rslt[I]:=Char((What and 15)+55);
          What:=What shr 4;
          End;
      Byte_to_Hex:=Rslt;
  End;    (* of Byte_to_Hex *)

 Procedure ErrorHandler(ErrorNumber,ErrorAddress:Integer);
 { This procedure tells the user the occured error type }
 Type ErrorA=Array[1..1] of Error_Type;
 Const Error_Type:Array[1..2] of String[10] = ('I/O ','Execution ');
 Var C:Char;
     ErrorP:^ErrorA;
     I,J,K: Integer;
     Tipe,Why:Byte;

 BEGIN
   Tipe:=Hi(ErrorNumber);
   Why :=Lo(ErrorNumber);
   Case Tipe of

     0:    Begin
           WriteLn(' Normal exit ');
           Readln
           End;

     1,
     2:    begin
                ClrScr;
                If Tipe =1 then
                   ErrorP:=Addr(IO_Errors)
                else
                   ErrorP:=Addr(Run_Errors);
                I:=1;
                While (ErrorP^[I].Error_ID <> Why) and (Why < $FF) do
                      I:=I+1;
                ClrEol;
                If ErrorP^[I].Error_ID = Why then
                   Writeln('   '+ErrorP^[I].Error_Msg)
                else if (Why < $FF) then
                   WriteLn(Error_Type[Tipe]+'Error # ',Byte_to_Hex(Why));
                ClrEol;
                WriteLn(' occured at address: ',
                          Byte_to_Hex(Hi(ErrorAddress)),
                          Byte_to_Hex(Lo(ErrorAddress)));
                WriteLn;
                WriteLn(' Please, report this error to your selling agent.');
                WriteLn(' Try to gather all informations concerning the ');
                WriteLn(' context of this error.');
                WriteLn;
                Write(' Any key to continue...');
                Read(C);
                End;
            End; (* of case *)
   Halt(Tipe+1);
 End; (* ErrorHandler *)

 END.
