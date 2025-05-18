{**************************************************************************
*                  Symbolic Parser with Polynomials                       *
*                 Example : '(A+B)^2' ==> A^2+2AB+B^2                     *
*                 TPW release 2.1 By J-P Moreau, Paris                    *
*                   (translated from Basic release)                       *
*                          (www.jpmoreau.fr)                              *
* ----------------------------------------------------------------------- *
* SAMPLE RUN:                                                             *
*                                                                         *
*                   SYMBOLIC PARSER for POLYNOMIALS                       *
*                                                                         *
*                   Example: (A+B)^2 ==> A^2+2AB+B^2                      *
*                                                                         *
* Inputs (0=SCREEN  1=algebra.dat file): 0                                *
* Input string to expand:                                                 *
* (A+B+C)^3                                                               *
*                                                                         *
* Detail analysis (0=NO  1=YES): 1                                        *
* Outputs (0=SCREEN  1=algebra.lst file): 0                               *
*                                                                         *
* SIMPLIFICATION:                                                         *
* (A+B+C)^3                                                               *
* POWER:                                                                  *
* (A^3+3A^2B+3A^2C+3AB^2+6ABC+3AC^2+B^3+3B^2C+3BC^2+C^3)                  *
* ADDITION:                                                               *
* A^3+3A^2B+3A^2C+3AB^2+6ABC+3AC^2+B^3+3B^2C+3BC^2+C^3                    *
*                                                                         *
* FINAL SIMPLIFICATION AND RESULT:                                        *
*                                                                         *
* A^3+3A^2B+3A^2C+3AB^2+6ABC+3AC^2+B^3+3B^2C+3BC^2+C^3                    *
*                                                                         *
* Evaluate another string (0=NO  1=YES): 1                                *
*                                                                         *
* (Assuming the input file algabra.dat contains the line:                 *
* (A+B+C+D)^4+(A+B+C-D)^4 )                                               *
*                                                                         *
* Inputs (0=SCREEN  1=algebra.dat file): 1                                *
* Input string to evaluate:                                               *
* (A+B+C+D)^4+(A+B+C-D)^4                                                 *
*                                                                         *
* Detail analysis (0=NO  1=YES): 1                                        *
* Outputs (0=SCREEN  1=algebra.lst file): 1                               *
*                                                                         *
* The file algebra.lst contains:                                          *
*                                                                         *
* STRING TO EVALUATE:                                                     *
*  (A+B+C+D)^4+(A+B+C-D)^4                                                *
*                                                                         *
* SIMPLIFICATION:                                                         *
*  (A+B+C+D)^4+(A+B+C-D)^4                                                *
* POWER:                                                                  *
*  (A^4+4A^3B+4A^3C+4A^3D+6A^2B^2+12A^2BC+12A^2BD+6A^2C^2+12A^2CD+6A^2D^2 *
*  +4AB^3+12AB^2C+12AB^2D+12ABC^2+24ABCD+12ABD^2+4AC^3+12AC^2D+12ACD^2+4A *
*  D^3+B^4+4B^3C+4B^3D+6B^2C^2+12B^2CD+6B^2D^2+4BC^3+12BC^2D+12BCD^2+4BD^ *
*  3+C^4+4C^3D+6C^2D^2+4CD^3+D^4)+(A+B+C-D)^4                             *
* ADDITION:                                                               *
*  A^4+4A^3B+4A^3C+4A^3D+6A^2B^2+12A^2BC+12A^2BD+6A^2C^2+12A^2CD+6A^2D^2+ *
*  4AB^3+12AB^2C+12AB^2D+12ABC^2+24ABCD+12ABD^2+4AC^3+12AC^2D+12ACD^2+4AD *
*  ^3+B^4+4B^3C+4B^3D+6B^2C^2+12B^2CD+6B^2D^2+4BC^3+12BC^2D+12BCD^2+4BD^3 *
*  +C^4+4C^3D+6C^2D^2+4CD^3+D^4+(A+B+C-D)^4                               *
* POWER:                                                                  *
*  A^4+4A^3B+4A^3C+4A^3D+6A^2B^2+12A^2BC+12A^2BD+6A^2C^2+12A^2CD+6A^2D^2+ *
*  4AB^3+12AB^2C+12AB^2D+12ABC^2+24ABCD+12ABD^2+4AC^3+12AC^2D+12ACD^2+4AD *
*  ^3+B^4+4B^3C+4B^3D+6B^2C^2+12B^2CD+6B^2D^2+4BC^3+12BC^2D+12BCD^2+4BD^3 *
*  +C^4+4C^3D+6C^2D^2+4CD^3+D^4+(A^4+4A^3B+4A^3C-4A^3D+6A^2B^2+12A^2BC-12 *
*  A^2BD+6A^2C^2-12A^2CD+6A^2D^2+4AB^3+12AB^2C-12AB^2D+12ABC^2-24ABCD+12A *
*  BD^2+4AC^3-12AC^2D+12ACD^2-4AD^3+B^4+4B^3C-4B^3D+6B^2C^2-12B^2CD+6B^2D *
*  ^2+4BC^3-12BC^2D+12BCD^2-4BD^3+C^4-4C^3D+6C^2D^2-4CD^3+D^4)            *
* ADDITION:                                                               *
*  A^4+4A^3B+4A^3C+4A^3D+6A^2B^2+12A^2BC+12A^2BD+6A^2C^2+12A^2CD+6A^2D^2+ *
*  4AB^3+12AB^2C+12AB^2D+12ABC^2+24ABCD+12ABD^2+4AC^3+12AC^2D+12ACD^2+4AD *
*  ^3+B^4+4B^3C+4B^3D+6B^2C^2+12B^2CD+6B^2D^2+4BC^3+12BC^2D+12BCD^2+4BD^3 *
*  +C^4+4C^3D+6C^2D^2+4CD^3+D^4+A^4+4A^3B+4A^3C-4A^3D+6A^2B^2+12A^2BC-12A *
*  ^2BD+6A^2C^2-12A^2CD+6A^2D^2+4AB^3+12AB^2C-12AB^2D+12ABC^2-24ABCD+12AB *
*  D^2+4AC^3-12AC^2D+12ACD^2-4AD^3+B^4+4B^3C-4B^3D+6B^2C^2-12B^2CD+6B^2D^ *
*  2+4BC^3-12BC^2D+12BCD^2-4BD^3+C^4-4C^3D+6C^2D^2-4CD^3+D^4              *
*                                                                         *
* FINAL SIMPLIFICATION AND RESULT:                                        *
*                                                                         *
*  2A^4+8A^3B+8A^3C+12A^2B^2+24A^2BC+12A^2C^2+12A^2D^2+8AB^3+24AB^2C+24AB *
*  C^2+24ABD^2+8AC^3+24ACD^2+2B^4+8B^3C+12B^2C^2+12B^2D^2+8BC^3+24BCD^2+2 *
*  C^4+12C^2D^2+2D^4                                                      *
*                                                                         *
* ----------------------------------------------------------------------- *
* Reference: "Calcul symbolique et informatique. Du calcul numérique au   *
*             calcul littéral (programmes en BASIC de A. DESHAYES),       *
*             Masson Paris, 1985" [BIBLI 02].                             *
**************************************************************************}
  PROGRAM ALGEBRA;
  Uses WinCrt,Strings;                   {Labels inherited from Basic}
  LABEL debut,130,220,250,290,300,340,350,370,460,470,490,510,540,560;
  LABEL 1000,1465,1470,1580,8000,res,fin;
  CONST MaxSize = 270;          {Maximum number of terms in expansion}
  TYPE  Str20 = string[20];               {Maximum length of one term}
        pTab = ^Tab;
        Tab = array[1..MaxSize] of Str20;
  VAR
      CHs : string;                                 {String to expand}
      RRs : array[0..770] of char;                      {Final answer}
      ALFA: array[1..26] of char;            {Table of variables A..Z}
      
      CPs,Rs,Ps,SIs: pTab;
      Es  : array[1..26] of Str20;
      IE  : array[1..26] of word;
      CP ,P,SI: array[1..MaxSize] of single;
      IPI,IPO : array[1..MaxSize] of integer;
      
      Ds,ECs,IEs,Ts,Ys,PRs,NBs,ENs,VLs,Ms,M1s: Str20;
      err,I,IC,IC1,IE1,IEN,IEX,INB,ISI,IVL,LE1,NB,NC,NL,RE  : integer;
      IAP,IDI,IDM,IFM,IH,IHC,IHO,IHP,INO,IPF,IPP,JPO,LE : integer;
      I4,IDC,IMD,IMG,IO1,IR,ISL,KS,NS,RF,RT : integer;

      Cs,SUs : char;
      EC,EX,M: single;
      E1s    : string[2];
      F1     : TEXT;
      Ok     : Boolean;


    Function Verif(VAR N:integer): boolean;
    Begin
      if N > MaxSize then
        Verif:=FALSE
      else
        Verif:=TRUE
    End; {Verif}

    Procedure Message; 
    {Message if a table index is greater than MaxSize.
     To be used after a call to Verif.}
    Begin
      if RT=0 then
        writeln('  Too many terms!')
      else
      begin
        writeln(F1,'  Too many terms!');
        Close(F1)
      end
    End;

    Function SIMPLI(IS1:integer): integer;
    {Simplify polynomial. Example: 3AB...-AB ==> 2AB }
    Var  IS2,NB : integer;
    Begin
      SIMPLI:=0;
      if Not Verif(IS1) then
      begin
        SIMPLI:=IS1;
        exit
      end;
      IF IS1 = 1 then
      begin
        NS:=1;
        exit
      end;
      IS2:=IS1-1;
      for NB:=1 to IS2 do
        if SI[NB] <> 0.0 then
        for KS:=NB+1 to IS1 do
          IF SIs^[KS] = SIs^[NB] THEN
          begin
            SI[NB]:=SI[NB]+SI[KS];
            SI[KS]:=0.0
          end;

      NS:=0;
      for I:=1 to IS1 do
        IF SI[I] <> 0.0 then
        begin
          NS:=NS+1;
          if Not Verif(NS) then
          begin
            SIMPLI:=NS;
            exit
          end;
          SIs^[NS]:=SIs^[I];
          SI[NS]:=SI[I]
        end;

      IF NS = 0 then
      begin
        NS:=1;
        SIs^[1]:=''
      end;
    End;  {Simpli}

    Procedure MONOME(Var Ts: Str20);
    {Simplify monom, ex.: BAB ==> AB^2 }
    Label 80,200;
    Var  I,I2,IVL : integer;
         IEs,SSs : Str20;
         SCs  : char;
    Begin
      for I:=1 to 26 do IE[I]:=0;
      I2:=1;
      SCs:=Ts[I2];
   80:IF I2 > length(Ts) then GOTO 200;
      for I:=1 to 26 do
        IF SCs=ALFA[I] then NL:=I;
      IE[NL]:=IE[NL]+1;
      I2:=I2+1;
      SCs:=Ts[I2];

      IF (SCs='^') and (I2<=length(Ts)) then
      begin
        I2:=I2+1;
        {modif. for version 2.0}
        if (Ts[I2+1] in ['0'..'9']) and (I2<length(Ts)) then
          SSs:=Ts[I2]+Ts[I2+1]
        else
          SSs:=Ts[I2];
        VAL(SSs,IVL,err);
        IE[NL]:=IE[NL]+IVL-1;
        I2:=I2+Length(SSs);   
        SCs:=Ts[I2]
      end;
      GOTO 80;

  200:Ts:='';
      for I2:=1 to 26 do
      begin
        if IE[I2] <> 0 then
        begin
          Ts:=Ts+ALFA[I2];
          IF IE[I2] > 1 then
          begin
            Str(IE[I2],IEs);
            Ts:=Ts+'^'+IEs
          end;
          IE[I2]:=0
        end
      end
    End;  { Monome }

    Function COPIE: integer;
    Label 9620,9840;
    Var  I,IHR,IZZ : integer;
         CO : single;
    Begin
      COPIE:=0;
      IF IDC <= 0 then
      begin
        for I:=NS downto 1 do
        begin
          IR:=IR+1;
          CPs^[IR]:=SIs^[I];
          CP[IR]:=SI[I]
        end;
        JPO:=IR+1;
        CPs^[JPO]:='(';
        CP[JPO]:=IPF;
        CP[IPF]:=JPO;
        IDM:=JPO-1;
        IAP:=IO1+1;
 9620:  IF CPs^[IAP]='$' then
        begin
          IAP:=IAP+1;
          GOTO 9620
        end;
        IF IO1=IDI then IDI:=IDI+IDC;
        IF IDC=0 then exit;
        for I:=JPO+1 to IO1 do
        begin
          CPs^[I]:='$';
          CP[I]:=0.0
        end;
        exit
      end;

      IF IO1<>IDI then
      begin
        IZZ:=IO1+1;
        IHR:=IHC-1;
        for I:=1 to IHR do
        begin
          CO:=CP[IPI[I]];
          IF CO > JPO then CP[IPI[I]]:=CO+IDC
        end;
        I:=IDI+IDC;
        if Not Verif(I) then
        begin
          COPIE:=I;
          exit
        end;
        for I:=IDI downto IZZ do
        begin
          CPs^[I+IDC]:=CPs^[I];
          IF CPs^[I]=')' then
          begin
            CP[I+IDC]:=CP[I]+IDC;
            GOTO 9840
          end;
          IF (CPs^[I]='(') and (CP[I] > IPF) then
          begin
            CP[I+IDC]:=CP[I]+IDC;
            GOTO 9840
          end;
          CP[I+IDC]:=CP[I];
 9840:  end; { for I }
      end;   { if }

      IDI:=IDI+IDC;
      if Not Verif(NS) then
      begin
        COPIE:=NS;
        exit
      end;
      for I:=NS downto 1 do
      begin
        IR:=IR+1;
        CPs^[IR]:=SIs^[I];
        CP[IR]:=SI[I]
      end;

      JPO:=IR+1;
      CPs^[JPO]:='(';
      CP[JPO]:=IPF;
      CP[IPF]:=JPO;
      IDM:=JPO-1;
      IAP:=JPO+1;
      While CPs^[IAP]='$' do IAP:=IAP+1;
    End;  { Copie }

    Function IMPR: boolean;
    {print partial or final results}
    Label 500;
    Var  I,IMO,INO,ISP,LR,LR1,LR2 : integer;
         Zs  : array[0..10] of char;
         Zcp : array[0..40] of char;
    Begin
      IMPR:=TRUE;
      ISP:=0;
      {Empty string RRs} 
      for i:=0 to SizeOf(RRs)-1 do RRs[i]:=#0;

      for I:=IDI downto 3 do
      begin
        IF CPs^[I] = '$' then GOTO 500;
        IMO:=0; INO:=0;
        IF (CPs^[I]>='A') and (CPs^[I]<='Z') then IMO:=1;
        if CPs^[I]='' then INO:=1;
        IF ISP=1 then
        begin
          IF ((IMO=1) or (INO=1)) and (CP[I]>0.0) and (StrLen(RRs)>0) then
            StrLCat(RRs,'+',Sizeof(RRs)-1);
          if CPs^[I]='(' then StrLCat(RRs,'+',Sizeof(RRs)-1);
        end;

        if ((IMO=1) or (INO=1)) and (Round(Abs(CP[I]))<>1)
          and (CP[I]<>0.0) then
        begin
          Str(Round(CP[I]),Zs);
          StrLCat(RRs,Zs,Sizeof(RRs)-1);
        end;

        if (Round(CP[I])=-1) and (CPs^[I]='') then
          StrLCat(RRs,'-1',Sizeof(RRs)-1);

        if (Round(CP[I])=-1) and (CPs^[I]<>'') then
          StrLCat(RRs,'-',Sizeof(RRs)-1);

        if (Round(CP[I])=1) and (CPs^[I]='') then
          StrLCat(RRs,'1',Sizeof(RRs)-1);

        if CP[I]<>0.0 then
        begin
          StrPCopy(Zcp,CPs^[I]);
          StrLCat(RRS,Zcp,Sizeof(RRS)-1)
        end;

        {Only for partial results}
        if (CP[I]=0.0) and (CPs^[I][1] in ['*','+','-']) then
        begin
          StrPCopy(Zcp,CPs^[I]);
          StrLCat(RRS,Zcp,Sizeof(RRS)-1)
        end;

        ISP:=0;
        IF (IMO=1) or (INO=1) or (CPs^[I]=')') then ISP:=1;
  500:end;

      LR:=StrLen(RRs);
      if RT=0 then     {print to screen}
      begin
      if LR=0 then write('  0')
      else
        write('  ');for i:=0 to 68 do WRITE(RRs[i]); writeln(RRs[69]);
      if LR > 69 then
      begin
        write('  ');for i:=70 to 138 do WRITE(RRs[i]); writeln(RRs[139]);
      end;
      if LR > 139 then
      begin
        write('  ');for i:=140 to 208 do WRITE(RRs[i]); writeln(RRs[209]);
      end;
      if LR > 209 then
      begin
        write('  ');for i:=210 to 278 do WRITE(RRs[i]); writeln(RRs[279]);
      end;
      if LR > 279 then
      begin
        write('  ');for i:=280 to 348 do WRITE(RRs[i]); writeln(RRs[349]);
      end;
      if LR > 349 then
      begin
        write('  ');for i:=350 to 418 do WRITE(RRs[i]); writeln(RRs[419]);
      end;
      if LR > 419 then
      begin
        write('  ');for i:=420 to 488 do WRITE(RRs[i]); writeln(RRs[489]);
      end;
      if LR > 489 then
      begin
        write('  ');for i:=490 to 558 do WRITE(RRs[i]); writeln(RRs[559]);
      end;
      if LR > 559 then
      begin
        write('  ');for i:=560 to 628 do WRITE(RRs[i]); writeln(RRs[629]);
      end;
      if LR > 629 then
      begin
        write('  ');for i:=630 to 698 do WRITE(RRs[i]); writeln(RRs[699]);
      end;
      if LR > 699 then
      begin
        write('  ');for i:=700 to 768 do WRITE(RRs[i]); writeln(RRs[769]);
      end;
      IF LR > 769 then
      begin
        WRITELN('  CHAINE RESULTAT TROP LONGUE (>',LR,') !');
        IMPR:=FALSE
      end
      end
      else   {print to text file ALGEBRA.LST}
      begin  {fill lines with blanks up to 70 characters}
      if LR<70 then
      begin
        LR2:=69-LR;
        for i:=LR to LR+LR2 do RRs[i]:=' ';
      end
      else if LR<770 then
      begin
        LR1:=LR DIV 70;
        LR2:=LR-LR1*70;
        for i:=LR to LR+69-LR2 do RRs[i]:=' '
      end;
      if LR=0 then write(F1,'  0')
      else
        write(F1,'   ');for i:=0 to 68 do WRITE(F1,RRs[i]); writeln(F1,RRs[69]);
      if LR > 69 then
      begin
        write(F1,'   ');for i:=70 to 138 do WRITE(F1,RRs[i]); writeln(F1,RRs[139]);
      end;
      if LR > 139 then
      begin
        write(F1,'   ');for i:=140 to 208 do WRITE(F1,RRs[i]); writeln(F1,RRs[209]);
      end;
      if LR > 209 then
      begin
        write(F1,'   ');for i:=210 to 278 do WRITE(F1,RRs[i]); writeln(F1,RRs[279]);
      end;
      if LR > 279 then
      begin
        write(F1,'   ');for i:=280 to 348 do WRITE(F1,RRs[i]); writeln(F1,RRs[349]);
      end;
      if LR > 349 then
      begin
        write(F1,'   ');for i:=350 to 418 do WRITE(F1,RRs[i]); writeln(F1,RRs[419]);
      end;
      if LR > 419 then
      begin
        write(F1,'   ');for i:=420 to 488 do WRITE(F1,RRs[i]); writeln(F1,RRs[489]);
      end;
      if LR > 489 then
      begin
        write(F1,'   ');for i:=490 to 558 do WRITE(F1,RRs[i]); writeln(F1,RRs[559]);
      end;
      if LR > 559 then
      begin
        write(F1,'   ');for i:=560 to 628 do WRITE(F1,RRs[i]); writeln(F1,RRs[629]);
      end;
      if LR > 629 then
      begin
        write(F1,'   ');for i:=630 to 698 do WRITE(F1,RRs[i]); writeln(F1,RRs[699]);
      end;
      if LR > 699 then
      begin
        write(F1,'   ');for i:=700 to 768 do WRITE(F1,RRs[i]); writeln(F1,RRs[769]);
      end;
      IF LR > 769 then
      begin
        WRITELN(F1,'   RESULT STRING TOO LONG (>769)!');
        IMPR:=FALSE
      end
      end
    End;

    Function PARE: integer;
    {MULTIPLY PARENTHESES}
    Label 5010,5054,5075;
    Var  ID1,ID2,IF1,IF2,IO2,IS1,IUM,L1,L2,N1,N2 : integer;
         I,J,MA,MB,NM : integer;
    Begin
      PARE:=0; Ok:=TRUE;
      for i:=1 to MaxSize do SIs^[i]:='';
      IO2:=JPO;
      IF2:=IPF;
      IF1:=IAP+1;

 5010:IF (CPs^[IF1]='$') and (IF1 < IDI) then
      begin
        IF1:=IF1+1;
        GOTO 5010
      end;
      IO1:=Round(CP[IF1]);
      ID1:=IO1-1;
      L1:=IF1+1;
      ID2:=IO2-1;
      L2:=IF2+1;
      IF ID1 = L1 then
      begin
        N1:=1;
        GOTO 5054
      end;
      IUM:=0;

      for I:=ID1 downto L1 do
      begin
        IUM:=IUM+1;
        if Not Verif(IUM) then
        begin
          PARE:=IUM;
          exit
        end;
        SIs^[IUM]:=CPs^[I];
        SI[IUM]:=CP[I]
      end;

      IS1:=IUM;
      I:=SIMPLI(IS1);
      if I>0 then
      begin
        PARE:=I;
        exit
      end;
      N1:=NS;
      IUM:=L1-1;
      for I:=N1 downto 1 do
      begin
        IUM:=IUM+1;
        if Not Verif(IUM) then
        begin
          PARE:=IUM;
          exit
        end;
        CPs^[IUM]:=SIs^[I];
        CP[IUM]:=SI[I]
      end;

 5054:IF ID2 = L2 then
      begin
        N2:=1;
        GOTO 5075
      end;
      IUM:=0;
      for I:=ID2 downto L2 do
      begin
        IUM:=IUM+1;
        if Not Verif(IUM) then
        begin
          PARE:=IUM;
          exit
        end;
        SIs^[IUM]:=CPs^[I];
        SI[IUM]:=CP[I]
      end;

      IS1:=IUM;
      I:=SIMPLI(IS1);
      if I>0 then
      begin
        PARE:=I;
        exit
      end;
      N2:=NS;
      IUM:=L2-1;
      for I:=N2 downto 1 do
      begin
        IUM:=IUM+1;
        if Not Verif(IUM) then
        begin
          PARE:=IUM;
          exit
        end;
        CPs^[IUM]:=SIs^[I];
        CP[IUM]:=SI[I]
      end;

 5075:MA:=N1+L1-1;
      MB:=N2+L2-1;
      NM:=0;
      for I:=MA downto L1 do
        for J:=MB downto L2 do
        begin
          Ts:=CPs^[I]+CPs^[J];

          MONOME(Ts);

          NM:=NM+1;
          if Not Verif(NM) then
          begin
            PARE:=NM;
            exit
          end;
          SIs^[NM]:=Ts;
          SI[NM]:=CP[I]*CP[J]
        end;

      IS1:=NM;
      I:=SIMPLI(IS1);
      if I>0 then
      begin
        PARE:=I;
        exit
      end;
      IR:=IF2;
      IDC:=NS-IO1+L2;

      I:=COPIE;
      if I>0 then
      begin
        PARE:=I;
        exit
      end;
      IF RE = 1 then
      begin
        if RT=0 then
          Writeln('     MULTIPLY PARENTHESES:')
        else
          Writeln(F1,'     MULTIPLY PARENTHESES:');
        if Not IMPR then
        begin
          Ok:=FALSE;
          PARE:=-1
        end
      end;
    End;  { PARE }

    Function PUIS: integer;    {Integer power} 
    Label 9080,9480,9485;
    Var I,IS1,ITP,IU,IUP,IX1,IZP,JP,KP : integer;
    Begin
      PUIS:=0;  Ok:=TRUE;
      {Current exponent is put into EX}
      if EC=0 then
      begin
        EX:=CP[IPP-1];
        IO1:=JPO;
      end;
      if EX=0.0 then
      begin
        SIs^[1]:='';
        SI[1]:=1.0;
        NS:=1;
        GOTO 9480
      end;
      IF EX=1.0 then
      begin
        IF EC=0.0 then GOTO 9485;
        IUP:=0;
        for I:=IDM downto IFM do
        begin
          IUP:=IUP+1;
          if Not Verif(IUP) then
          begin
            PUIS:=IUP;
            exit
          end;
          SIs^[IUP]:=CPs^[I];
          SI[IUP]:=CP[I]
        end;
        NS:=IDM-IFM+1;

        GOTO 9480
      end;

      ITP:=IDM-IFM+1;
      IF ITP=1 then
      begin
        NS:=1;
        Ps^[1]:=CPs^[IFM];
        P[1]:=CP[IFM];
        GOTO 9080
      end;
      I:=0;
      for JP:=IDM downto IFM do
      begin
        I:=I+1;
        if Not Verif(I) then
        begin
          PUIS:=I;
          exit
        end;
        SIs^[I]:=CPs^[JP];
        SI[I]:=CP[JP]
      end;

      IS1:=ITP;
      I:=SIMPLI(IS1);
      if I>0 then
      begin
        PUIS:=I;
        exit
      end;

      IX1:=IFM+NS-1;
      IUP:=IX1+1;
      if Not Verif(IUP) then
      begin
        PUIS:=IUP;
        exit
      end;
      for I:=1 to NS do
      begin
        IUP:=IUP-1;
        CPs^[IUP]:=SIs^[I];
        CP[IUP]:=SI[I];
        Ps^[I]:=SIs^[I];
        P[I]:=SI[I]
      end;

 9080:IEX:=Round(EX)-1;
      for KP:=1 to IEX do
      begin
        IZP:=0;
        for I:=IX1 downto IFM do
          for IU:=1 to NS do
          begin
            IZP:=IZP+1;
            if Not Verif(IZP) then
            begin
              PUIS:=IZP;
              exit
            end;
            Ts:=CPs^[I]+Ps^[IU];

            MONOME(Ts);

            SIs^[IZP]:=Ts;
            SI[IZP]:=CP[I]*P[IU]
          end;

        IS1:=IZP;
        I:=SIMPLI(IS1);
        if I>0 then
        begin
          PUIS:=I;
          exit
        end;

        if Not Verif(NS) then
        begin
          PUIS:=NS;
          exit
        end;
        for I:=1 to NS do
        begin
          Ps^[I]:=SIs^[I];
          P[I]:=SI[I]
        end;    
      end;  {for kp}

 9480:IDC:=NS-IO1+IPF+1;

      IR:=IPF;
      COPIE;

      IF EC=1.0 then IFM:=IPF+1;
 9485:IF EC=0.0 then
      begin
        CPs^[IPP]:='$';
        CP[IPP]:=0.0;
        CPs^[IPP-1]:='$';
        CP[IPP-1]:=0.0;
        IPP:=IPP-2
      end;

      IF RE = 1 then
      begin
        if RT=0 then
          Writeln('     POWER:')
        else
          Writeln(F1,'     POWER:');
        if Not IMPR then
        begin
          Ok:=FALSE;
          PUIS:=-1
        end
      end;
      EC:=0.0
    End;  { PUIS }

    Procedure ADDI;
    {Suppress parentheses preceded by a plus sign}
    begin
      CP[JPO]:=0.0;
      CP[IPF]:=0.0;
      CPs^[JPO]:='$';
      CPs^[IPF]:='$';
      IF CPs^[IAP]='+' then
      begin
        CPs^[IAP]:='$';
        CP[IAP]:=0.0
      end;
      IF RE = 1 then
      begin
        if RT=0 then
          Writeln('     ADDITION:')
        else
          Writeln(F1,'     ADDITION:');
        IMPR
      end;
    End;  { Addi }

    Procedure SOUS;
    {Suppress parentheses preceded by a minus sign}
    begin
      for I:=IDM downto IFM do
        CP[I]:=-CP[I];
      CP[JPO]:=0.0;
      CP[IPF]:=0.0;
      CPs^[JPO]:='$';
      CPs^[IPF]:='$';
      CPs^[IAP]:='$';
      CP[IAP]:=0.0;
      IF RE = 1 then
      begin
        if RT=0 then
          Writeln('     SUBSTRACTION:')
        else
          Writeln(F1,'     SUBSTRACTION:');
        IMPR
      end;
    End;  { Sous }

    Procedure OPERA(CHs:string;Var IC:integer;Var Cs:char;Var NC:integer;
                    Var Ds:Str20);
    {Seek a derivation operator, £}
    Label 670,680,700;
    Begin
      IF Cs='£' then
      begin
  670:  Ds:='';
  680:  IF Cs <> '(' then
        begin
          Ds:=Ds+Cs;
          IC:=IC+1;
          Cs:=CHs[IC];
          GOTO 680
        end;
        NC:=NC+1;
        CPs^[NC]:=Ds
      end;
  700:IF Cs = '(' then
      begin
        NC:=NC+1;
        CPs^[NC]:='(';
        IC:=IC+1;
        Cs:=CHs[IC];
        GOTO 700
      end;
      IF Cs='£' then GOTO 670
    End;  {Opera}

    Procedure DERI;
    Begin
      {Derivation not implemented here}
    End;

    PROCEDURE WinCrtInit(Nom:PChar);
    {Open a CRT window (specific to MS Windows environment) }
    BEGIN
      WindowOrg.X:=10;
      WindowOrg.Y:=20;
      WindowSize.X:=620;
      WindowSize.Y:=440;
      StrCopy(WindowTitle,Nom);
      InitWinCrt;
    END;

  {Main Program}
  BEGIN
      New(CPs); New(Rs); New(Ps); New(SIs);
      RF:=0;
      WinCrtInit(' ALGEBRA');
      Writeln;
      Writeln('                     SYMBOLIC PARSER FOR POLYNOMIALS');
      Writeln;
      Writeln('                     Example: (A+B)^2 ==> A^2+2AB+B^2');
      Writeln;
      for I:=1 to 26 do ALFA[I]:=Chr(I+64);
debut:if RF=1 then
      begin 
        Clrscr;
        writeln
      end;
      Write('  Inputs (0=Screen  1=File algebra.dat): ');
      Readln(RT);
      if RT=0 then
      begin
        Writeln('  Input string to expand: ');
        Write('  '); readln(CHs);
        Writeln
      end
      else
      begin
        Assign(F1,'algebra.dat');
        reset(F1);
        Read(F1,CHs);
        Close(F1);
        Writeln('  String to expand: ');
        IC:=length(CHs);
        if IC<71 then
        begin
          Write('  ');for i:=1 to IC-1 do write(CHs[i]);writeln(CHs[IC]);
        end
        else if IC<141 then
        begin
          Write('  ');for i:=1 to 69 do write(CHs[i]);writeln(CHs[70]);
          WRITE('  ');for i:=71 to IC-1 do write(CHs[i]);writeln(CHs[IC]);
        end
        else if IC<211 then
        begin
          Write('  ');for i:=1 to 69 do write(CHs[i]);writeln(CHs[70]);
          WRITE('  ');for i:=71 to 139 do write(CHs[i]);writeln(CHs[140]);
          WRITE('  ');for i:=141 to IC-1 do write(CHs[i]);writeln(CHs[IC]);
        end;
        Writeln
      end;
      Write('  Show computing details (0=No  1=Yes): '); readln(RE);
      Write('  Outputs (0=Screen  1=File algebra.lst): ');
      Readln(RT);
      if RT=1 then
      begin
        Assign(F1,'algebra.lst');
        rewrite(F1);
        Writeln(F1);
        Writeln(F1,'     STRING TO EXPAND:');
        IC:=length(CHs);
        if IC<71 then
        begin
          Write(F1,'   ');for i:=1 to IC-1 do write(F1,CHs[i]);writeln(F1,CHs[IC]);
        end
        else if IC<141 then
        begin
          Write(F1,'   ');for i:=1 to 69 do write(F1,CHs[i]);writeln(F1,CHs[70]);
          WRITE(F1,'   ');for i:=71 to IC-1 do write(F1,CHs[i]);writeln(F1,CHs[IC]);
        end
        else if IC<211 then
        begin
          Write(F1,'   ');for i:=1 to 69 do write(F1,CHs[i]);writeln(F1,CHs[70]);
          WRITE(F1,'   ');for i:=71 to 139 do write(F1,CHs[i]);writeln(F1,CHs[140]);
          WRITE(F1,'   ');for i:=141 to IC-1 do write(F1,CHs[i]);writeln(F1,CHs[IC]);
        end;
        Writeln(F1)
      end;

      {Preliminary loop to add a * sign between two
       parentheses )(,  or a number and (.  }
      for ic:=1 to length(CHs) do
      begin
        Cs:=CHs[ic];
        if (Cs=')') and (CHs[ic+1]='(') then
          Chs:=Copy(Chs,1,ic)+'*'+Copy(CHs,ic+1,length(CHs));
        if (Cs IN ['0'..'9']) and (CHs[ic+1]='(') then
          Chs:=Copy(Chs,1,ic)+'*'+Copy(CHs,ic+1,length(CHs));
      end;

      IC:=1;
      Cs:=CHs[IC];
      NC:=0;

      OPERA(CHs,IC,Cs,NC,Ds);

      { boucle analyse de la chaîne }
  130:IF IC > length(CHs) then goto 1000;
      if Cs='-' then ISI:=-1 else ISI:=1;
      if (Cs='+') or (Cs='-') then
      begin
        IC:=IC+1;
        Cs:=CHs[IC]
      end;
      if Cs='(' then goto 220;
      if Cs<>'£' then goto 250;
      {Opening parenthesis}
  220:PRs:=CHs[IC-1];
      NC:=NC+1;
      CPs^[NC]:= PRs;

      OPERA(CHs,IC,Cs,NC,Ds);

      Goto 130;
  250:if Cs<>')' then goto 460;
      {Closing parenthesis}
  290:if Cs<>')' then goto 300;
      NC:=NC+1;
      CPs^[NC]:=')';
      IC:=IC+1;
      Cs:=CHs[IC];
      GOTO 290;
  300:SUs:=CHs[IC+1];
      IF (SUs='(') or (SUs='£') then
      begin
        IC:=IC+1;
        Cs:=SUs;                       
        GOTO 220
      end;
      IF Cs<>'^' then GOTO 370;
    {Integer Power}
      NC:=NC+1;
      CPs^[NC]:='^';
      EX:=1.0;
  340:IF Cs<>'^' then GOTO 350;
        NC:=NC-1;
        if CHs[IC+2] in ['0'..'9'] then
          E1s:=CHs[IC+1]+CHs[IC+2]  {2 digits exponent}
        else
          E1s:=CHs[IC+1];           {1 digit exponent}
        VAL(E1s,IE1,err);
        EX:=EX*IE1;
        IC:=IC+length(E1s);
        Cs:=CHs[IC];
      GOTO 340;

  350:NC:=NC+1;
      CP[NC]:=EX;
      SUs:=CHs[IC+1];
      IF (SUs='(') or (SUs='£') then   {£ is for derivation}
      begin
        IC:=IC+1;
        Cs:=SUs;
        GOTO 220
      end;

  370:IF Cs<>'*' then GOTO 130;
{Begin MULTIPLY MONOMS}
      NC:=NC+1;
      CPs^[NC]:='*';
      NC:=NC+1;
      CPs^[NC]:='';
      CP[NC]:=1.0;
      GOTO 470;
{End MULTIPLY MONOMS}
  460:NC:=NC+1;
      CP[NC]:=ISI;
      CPs^[NC]:='';
{CONVERT/COMPACT MONOM}
  470:IF (CHs[IC]<'0') or (CHs[IC]>'9')
        or (IC > length(CHs)) then GOTO 510;
      NBs:='';
      {Modif. for version 2.0}
      if length(E1s)=2 then IC1:=IC-1 else IC1:=IC;
      {Read an integer of unknown length}
      While (CHs[IC1] in ['0'..'9']) and (IC1 <= length(CHs)) do
      begin
        NBs:=NBs+CHs[IC1];
        Inc(IC1);
      end;
      VAL(NBs,NB,err);
      IC:=IC+length(NBs);
      Cs:=CHs[IC];
  490:IF Cs='^' then
      begin
        IC:=IC+1;
        ENs:=CHs[IC];
        VAL(ENs,IEN,err);

      {NB power IEN}
        INB:=NB;
        NB:=1;
        for i:=1 to IEN do NB:=NB*INB;

        IC:=IC+1;
        Cs:=CHs[IC];
        GOTO 490
      end;
      CP[NC]:=CP[NC]*NB;
  510:IF (Cs < 'A') or (Cs > 'Z') or (IC > length(CHs)) then GOTO 560;
      for i:=1 to 26 do
        if ALFA[i]=Cs then NL:=i;
      IE[NL]:=IE[NL]+1;
      IC:=IC+1;
      Cs:=CHs[IC];
      IEX:=1;

  540:IF Cs='^' then
      begin
        IC:=IC+1;
        if CHs[IC+1] in ['0'..'9'] then
          ECs:=CHs[IC]+CHs[IC+1]
        else
          ECs:=CHs[IC];
        Val(ECs,IVL,err);
        IEX:=IEX*IVL;
        IC:=IC+length(ECs);
        Cs:=CHs[IC];
        Goto 540
      end;
      IE[NL]:=IE[NL]+IEX-1;
  560:if Cs='*' then
      begin
        Inc(IC);
        Cs:=CHs[IC];
        if Cs='-' then
        begin
          CP[NC]:=-CP[NC];
          Inc(IC);
          Cs:=CHs[IC]
        end
      end;
      if (Cs>='A') and (Cs<='Z') then goto 470;
      if (Cs>='0') and (Cs<='9') then goto 470;
      for i:=1 to 26 do
      begin
        if IE[i]<>0 then
        begin
          CPs^[NC]:=CPs^[NC]+ALFA[I];
          if IE[i]>1 then
          begin
            Str(IE[i],IEs);
            CPs^[NC]:=CPs^[NC]+'^'+IEs
          end;
          IE[i]:=0
        end { if IE }
      end; { for i }
      {end closing parenthesis}
      goto 130;
      {end string analysis loop}

      {Simplify monoms, if necessary}
 1000:IDI:=NC+2;
      if Not Verif(IDI) then
      begin
        Message;   {too many terms}
        goto fin
      end;
      IH:=0;
      INO:=0;
      LE:=IDI DIV 2;
      for NC:= 1 to LE do
      begin
        ECs:=CPs^[NC];
        CPs^[NC]:=CPs^[IDI-NC+1];
        CPs^[IDI-NC+1]:=ECs;
        EC:=CP[NC];
        CP[NC]:=CP[IDI-NC+1];
        CP[IDI-NC+1]:=EC
      end; { for NC }
      for i:=3 to IDI do
      begin
        if CPs^[i]=')' then
        begin
          Inc(IH);
          IPI[IH]:=i;
          Inc(INO);
          IPO[INO]:=i
        end;
        if CPs^[i]='(' then
        begin
          IHO:=IPO[INO];
          CP[IHO]:=i;
          CP[i]:=IHO;
          Dec(INO)
        end
      end; { for i }
      if RE=1 then
      begin
        if RT=0 then
          writeln('     SIMPLIFICATION MONOMS:')
        else
          writeln(F1,'     SIMPLIFICATION MONOMS:');
        IMPR
      end;

      {Operations on monoms}
      IHP:=IH;  IHC:=IH;   EC:=0.0;
1465: if IHP=0 then goto res;
1470:   IPF:=IPI[IHC];
        JPO:=Round(CP[IPF]);
        IDM:=JPO-1;
        IFM:=IPF+1;
        IPP:=IPF-1;
        IAP:=JPO+1;
        if (CPs^[IPP]='^') and (CPs^[IPP-1]='(') then
        begin
          Dec(IHC);
          Goto 1470
        end;
        if CPs^[IPP]='^' then
          if PUIS<>0 then
          begin
            if Ok=TRUE then
            begin
              Message;   {Power failed}
              goto fin
            end
            else
              RE:=0
          end;
        if CPs^[IAP]='^' then
        begin
          EX:=0.0;
          for i4:=IFM to IDM do EX:=EX+CP[i4];
          IFM:=IAP+2;
          IDM:=Round(CP[IAP+1])-1;
          IO1:=Round(CP[IAP+1]);
          Dec(IHP);
          EC:=1.0;
          if PUIS<>0 then
          begin
            if Ok=TRUE then
            begin
              Message;    {Power failed}
              goto fin
            end
            else
              RE:=0
          end
        end;  { if CPs^[IAP }
        Cs:=CPs^[IAP][1];
        if Cs='£' then DERI;
        Ms:='';
        M:=1.0; IMG:=0; IMD:=0;
        if (CPs^[IAP]='*') and (CPs^[IAP+1]<>')') and (CPs^[IAP+1]<>'$') then
          IMG:=1;
        Cs:=CPs^[IPP-1][1];
        if (CPs^[IPP]='*') and (CPs^[IPP-1]<>'(') and (Cs<>'£') then
          IMD:=1;
        if (IMG<>0) or (IMD<>0) then goto 8000;
 1580:  if CPs^[IAP]='*' then
        begin
          if PARE<>0 then
          begin
            if Ok=TRUE then
            begin
              Message;     {PARE failed}
              goto fin
            end
            else
              RE:=0
          end;
          Dec(IHP)
        end;
        if CPs^[IPP]='*' then
        begin
          Dec(IHC);
          Goto 1470
        end;
        if CPs^[IAP]='-' then SOUS else ADDI;
        Dec(IHP); Dec(IHC);
        Goto 1465;

      {Multiply monoms}
 8000:if IMG=1 then
      begin
        Ms:=Ms+CPs^[IAP+1];
        M:=M*CP[IAP+1];
        CPs^[IAP]:='$';
        CP[IAP]:=0.0;
        CPs^[IAP+1]:='$';
        CP[IAP+1]:=0.0;
        IAP:=IAP+2
      end; {if IMG=1}

      if IMD=1 then
      begin
        Ms:=Ms+CPs^[IPP-1];
        M:=M*CP[IPP-1];
        CPs^[IPP]:='$';
        CP[IPP]:=0.0;
        CPs^[IPP-1]:='$';
        CP[IPP-1]:=0.0;
        IPP:=IPP-2
      end; {if IMD=1}

      for i:=IDM Downto IFM do
      begin
        if CP[i]<>0.0 then
        begin
          CP[i]:=CP[i]*M;
          Ts:=CPs^[i]+Ms;
          MONOME(Ts);
          CPs^[i]:=Ts
        end
      end; { for i:=IDM }
      if RE=1 then
      begin
        if RT=0 then
          writeln('     MULTIPLICATION MONOMS:')
        else
          writeln(F1,'     MULTIPLICATION MONOMS:');
        if Not IMPR then goto fin
      end;
      Goto 1580;

      {section simplifications & results}
 res: if RT=0 then
      begin
        Writeln;
        Writeln('  FINAL SIMPLIFICATION AND RESULTS:');
        Writeln
      end
      else
      begin
        Writeln(F1);
        Writeln(F1,'   FINAL SIMPLIFICATION AND RESULTS:');
        Writeln(F1)
      end;
      ISL:=IDI-1;
      for i:=3 to ISL do
        if CP[i]<>0.0 then
          for KS:=i+1 to IDI do
            if CPs^[KS]=CPs^[i] then
            begin
              CP[i]:=CP[i]+CP[KS];
              CP[KS]:=0.0
            end;
      IMPR;
      if RT=1 then Close(F1);
fin:  Writeln;
      Write('  Continue (0=No  1=Yes): ');
      Readln(RF);
      if RF=1 then
        goto debut
      else
        begin
        Dispose(CPs); Dispose(Rs); Dispose(Ps); Dispose(SIs);
        DoneWinCrt
      end

 END.

 {end of file algebra.pas}