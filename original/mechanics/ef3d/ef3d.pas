{*************************************************************
*           Program EF3D - Release 2.0, december 1991        *
*         (C) Copyright Jean-Pierre DUMONT 1987,1991         *
*                       Jean-Pierre MOREAU 2006              *
* ---------------------------------------------------------- *
* Release History:                                           *
* 1) DOS                                                     *
* Calculate eigenfrequencies of circular beams               *
* 6/9/87   First 1D-version (2ddl by node: v and theta)      *
* 13/9/87  Imposed DOF elimination                           *
*          Data from text file                               *
*          Added local masses for 1D-version                 *
* 18/9/87  Beam element with shearing (1D-version)           *
* 27/9/87  3D-Beam element                                   *
* 3/10/87  Bar & torsion bar elements                        *
* 6/10/87  Spring elements between two DOF                   *
* 7/10/87  Data reorganization                               *
* 14/10/87 Arbitrary node numbering                          *
* 18/10/87 Storage of matrices K,M,Phi in heap               *
* 25/10/87 Title, Error handling in Reduc1, Tql2             *
* 13/11/87 ErrorHandler, table of modes                      *
* 10/06/88 Saving results to disk, etc.                      *
* ---------------------------------------------------------- *
* 2) WINDOWS                                                 *
* May 1998   French Release 2.0 By J-P Moreau                *
* May 2006   English Release 3.0 for Internet By J-P Moreau  *
* ---------------------------------------------------------- *
* Required Units:                                            *
*            Errors, ErrorHan, Declara1, AlgebLin,           *
*            Utilit, Typedef, Graph_2D.                      *
*************************************************************}
PROGRAM EF3D;

USES WinCrt1,WinTypes,WinProcs,Strings,       {Borland Pascal}
     Errors,ErrorHan,Declara1,AlgebLin,           {J-P Dumont}
     Utilit,graph_2D;                             {J-P Moreau}

{Note: WinCrt1 is a modified version of Borland WinCrt unit}

 Label debut;

 Const
     Tol =5.E-15;

     Menu_Item:  Array[1..7] of String[20]=
                 ('  General Parameters',
                  '  Table of modes    ',
                  '  Shape of modes    ',
                  '  Draw a mode       ',
                  '  Orthogonality Test',
                  '  Exit              ',
                  '  Other Problem     ');

     Messg : String = 'To continue, strike any key...';

 VAR   {global variables}

     ES,EI,GJ,PhiY,VCol,Longueur,elmass,Roll_Inertia : pVect;  {vectors}
     d,dl,e,omega2,frequence, U,V, masse_generalisee : pVect;

     K,M,Phi          : pMatc;   {square matrices}
     Xcoord           : Vecteur_Position;
     ElType           : Vecteur_Byte_court;
     Unknown          : Vecteur_Byte_long;
     idl              : ddl_status;
     elnode           : connexion;
     Numero           : Classe;
     General          : Vecteur_Ent;
     sortie           : Array[1..6] of arg_reel;

     i,j,l,n,Nel,Nodes,Mode,Lines,iglob,index          : Integer;
     ddl,Free,Prescribed,Noeuds_a_sortir               : Integer;
     Duree,biggest,produit                             : arg_reel;
     Error,Resultat,NewSW,OrthoSW,SaveSw,StatSw,ModeSw : Boolean;

     Donnees                 : TEXT;
     Dump                    : File of Gen_Data_Record;
     Dump1                   : File of Vecteur_Colonne;
     Buff                    : Array[0..$FF] of Byte;
     Titre                   : Enregistrement;
     Probleme,Fichier        : Nom_Fichier;
     St                      : String[2];
     jj                      : byte;         {main menu option}
     Exist                   : boolean;
     General_Data            : Gen_Data_Record;
     Msg                     : Array[0..60] of char;
     answ                    : Char;

{************************************************}
PROCEDURE Pause(S:Pchar);
begin
  gotoxy(40,28);
  TextOut(CrtDC,15,402,S,strlen(S));
  Readkey
end;
{************************************************}
 PROCEDURE Read_An_Integer(Var ix:Integer);
 begin
   repeat
     {$I-}
     Read(ix);
     {$I+}
   until(IOresult=0);
 end;
(************************************************)
 FUNCTION Power(x:arg_reel;n:integer): arg_reel;
 { Calculate x power n }
 var i,m : integer;
     result :arg_reel;
 begin
   result :=1.0;
   if n=0 then
   begin
     Power:=result;
     exit
   end;
   m:=  n;
   if n<0 then m:=-n;
   for i:=1 to m do result :=x*result;
   Power :=result;
   if n<0 then Power:=1.0/result
 end;
(*************************************************************)
 FUNCTION log10 (x : arg_reel) : arg_reel;
 begin
   if x>0 then log10:= ln(x)/ln(10)
          else if x=0 then writeln ('Null Argument for Log10')
                      else writeln ('Log10 Argument =', x);
 end;
{*************************************************************
 Put in upper-case letters                                   }
 FUNCTION Majuscule(Mot:Mot_Cle):Mot_Cle;
 VAR i : Byte;
     Long : Byte Absolute Mot;
 BEGIN
   FOR i:=1 to Long do Mot[i]:=Upcase(Mot[i]);
   Majuscule:=Mot;
 END;
(************************************************************)
 FUNCTION Answer: Boolean;
 Var Xhold,Yhold : Integer;
     Car         : Char;
 Begin
   Xhold:=WhereX;
   Yhold:=WhereY;
   Answer:=True;
   Repeat
     GotoXY(Xhold,Yhold);
     ClrEol;
     Car:=UpCase(Readkey);
     if (Car='N') then Answer:=False
   until Car in ['Y','N','O']
 End;
(************************************************************)
 PROCEDURE ClearScreen;
 Begin
   ClrScr
 End;
{************************************************************
  draw a double frame of thickness 3 pixels  }
  PROCEDURE Cadre(P:HDC;a,b,c,d:INTEGER);
  begin
    Rectangle(P,a,b,c,d);
    Rectangle(P,a+3,b+3,c-3,d-3)
  end;
(************************************************************)
 PROCEDURE Erreur(i:Byte);
 BEGIN
   Writeln; writeln('  ERROR:'); writeln;
   Write('     ');
   Case i of
   1:  WriteLn('Number of elements is > ',Maxddl:5);
   2:  WriteLn('Number of independant unknowns is > ',Maxc);
   3:  WriteLn('Number of nodes is > ',Maxndp);
   4:  WriteLn('Spring added to unknown node!');
   5:  WriteLn('Mass added to unknown node!');
   6:  WriteLn('First node of element is unknown');
   7:  WriteLn('Second node of element is unknown');
   8:  WriteLn('Mass matrix is singular!');
   9:  WriteLn('No convergence of eigenvalue');
   10: WriteLn(': Element kind unknown!');
   11: WriteLn('Impossible to open file');
   12: WriteLn('Option not implemented');
   13: WriteLn('Exit after 3 tries');
   14: Writeln('Wrong parameter!');
   15..255 : WriteLn('Fatal error!');
   end;
   Readln;
   Halt
 End;
(************************************************************)
 PROCEDURE Open_Input_File;
 Var i,taille,Essai   : Byte;
     Filename         : Nom_Fichier;

 BEGIN
   ClrScr; FileName:='';
   Essai:=0;
   Repeat
     Clrscr; GotoXY(1,2);
     Write('     Name of data file ? ');
     Readln(FileName);
     if copy(filename,length(filename)-3,1) <> '.' then
       filename:=filename+'.DAT';   {DAT suffix by default}
     Assign(Donnees,FileName);
     (*$I-*) Reset(Donnees) (*$I+*);
     Exist:=(IOresult=0);
     if Not Exist then
     begin
       Writeln;
       Write('     File: ',FileName,' does not exist, retry!');
       essai:=essai+1;
       MessageBeep(0); readln
     end
   Until (exist) or (Essai>2);
   if (Essai>2) then Erreur(13);
   Taille:=Length(FileName);
   if copy(filename,taille-3,1)='.' then probleme:=copy(filename,1,taille-4)
                                    else probleme:=filename;
 END;  {Open_Input_File}
(************************************************************)
PROCEDURE Open_Save_Read(Chaine:Suffixe);
Var Filename : Nom_Fichier;

BEGIN
  FileName:=Probleme+'.'+chaine;
  {$I-}
  if chaine='GEN' then
                  begin
                    Assign(Dump,FileName);
                    Reset(Dump)
                  end
                  else
                  begin
                    Assign(Dump1,FileName);
                    Reset(Dump1)
                  end;
  {$I+}
  Exist:=(IOresult=0);
  if (Not Exist) then
  begin
    WriteLn('     File: ',FileName,' does not exist.');
    Readkey;
    Erreur(11);
  end; (* if *)
END;  {Open_Save_Read}
(************************************************************)
PROCEDURE Open_Save_Write(Chaine:Suffixe);
Var OK       : Boolean;
    Essai    : Byte;
    Chr      : Char;
    X,Y      : Integer;
    Filename : Nom_Fichier;

BEGIN
  FileName:=Probleme+'.'+chaine;
  {$I-}
  if chaine='GEN' then
                  begin
                    Assign(Dump,FileName);
                    Reset(Dump)
                  end
                  else
                  begin
                    Assign(Dump1,FileName);
                    Reset(Dump1)
                  end;
  {$I+}
  Exist:=(IOresult=0);
  X:=WhereX;Y:=WhereY;
  if (Exist) then
  begin
    if chaine='MOD' then WriteLn;
    WriteLn('     File: ',FileName,' already exists');
    Write('     Overwrite (Y/N) ? ');
    OK:=Answer;
    GotoXY(1,WhereY);ClrEol;
    GotoXY(1,WhereY-1);ClrEol;
  end; (* if *)
  if (Not Exist) or (OK) then
  begin
    if chaine='GEN' then Rewrite(Dump)     {All ok}
                    else Rewrite(Dump1);
    exit;
  end
  else
  Write('     Please, give another name:');
  essai:=0;
  X:=WhereX;
  Y:=WhereY;
  Repeat
    ReadLn(FileName);
    Assign(Dump,FileName);
    (*$I-*) Reset(Dump) (*$I+*);
    Exist:=(IOresult=0);
    Ok:=true;
    if (Exist) then
    begin
      WriteLn('     File: ',FileName,' already exists');
      Write('     Overwrite (Y/N) ? ');
      OK:=Answer;
      gotoXY(1,WhereY);ClrEol;
      gotoXY(1,WhereY-1);ClrEol;
      essai:=essai+1;
      GotoXY(X,Y);
      ClrEol;
    end
  Until OK or (Essai >2);
  if (Essai>2) then Erreur(13);
  GotoXY(1,Y);ClrEol;
  if chaine='GEN' then Rewrite(Dump)
                  else Rewrite(Dump1);
END;     {Open_Save_Write
************************************************************}
PROCEDURE Points_Nodaux(Var n:Integer;
                        Var Numero: Classe;
                        Var x:vecteur_Position;
                        Var idl:ddl_status);

Var i,j,npt : Integer;
{ Numero = arbitrary node number given by user}
Begin
ReadLn(Donnees,n);
If (n>maxndp) then Erreur(3);
Free:=0;
Prescribed:=0;
npt:=0;
FOR i:=1 to n do
    Begin
    ReadLn(Donnees,Numero[i],x[1][i],x[2][i],x[3][i],
           idl[1][i],idl[2][i],idl[3][i],idl[4][i],idl[5][i],idl[6][i]);
    npt:=npt+1;
    For j:=1 to 6 do
        begin
        if (idl[j][i]=0) then Free:=Free+1             (* idl =0 libre *)
                     else Prescribed:=Prescribed+1;    (* idl =1 Prescribed *)
        end; (* j *)

    end; (* i *)
Read(Donnees,j);
if(j<>999) then Erreur(2);
END; (* Points_Nodaux *)
(*****************************************************)
FUNCTION Num_Interne(i:Integer):Byte;
Var k     : byte;
BEGIN
For k:=1 to Nodes do
    begin
    if (Numero[k]=i) then
       begin
       Num_Interne:=k;
       exit;
       end;
    end;
Erreur(5);
End;

(***********************************************************)
PROCEDURE Elements(Var Nel:integer;
                   Var ES,EI,GJ,Longueur,elmass,Roll_Inertia,Phiy :pVect;
                   Var x     : Vecteur_Position;
                   Var np    : Connexion;
                   Var ElType: Vecteur_Byte_court);

Var i,i1,i2,j,isave  : Integer;
    Typ              : Byte;
    Blanc            : Char;
    Type_Element     : Mot_Cle;

    Dext,Dint,Rho,Long,Section,Inertie,E,G,Ky : arg_reel;

PROCEDURE Numero_Interne(Var i1,i2 :Integer);
Var k     : byte;
Label skip,fin;
BEGIN
  For k:=1 to Nodes do
  begin
  if (Numero[k]=i1) then
  begin
    i1:=k;
    goto skip
    end
  end;
  Erreur(6);
skip: For k:=1 to Nodes do
  begin
    if (Numero[k]=i2) then
    begin
      i2:=k;
      goto fin;
    end
  end;
  Erreur(7);
Fin:if (i2<i1) then
  begin
    isave:=i1;
    i1:=i2;
    i2:=isave;
  end
End;

BEGIN
  i:=0;     {Number of elements in model}
  Write('     The model includes: ');

  Repeat
    Typ:=0;
    Read(Donnees,Nel);

    if (Nel<>999) THEN
    Begin
      ReadLn(Donnees,blanc,Type_Element);
      Type_Element:=Majuscule(Type_Element);
      GotoXY(26,WhereY);
      WriteLn(Nel:10,'  elements (',Type_Element,')');
      If (Type_Element='BEAM') THEN typ:=1;
      If (Type_Element='BEAMS') THEN typ:=2;
      If (Type_Element='BARRE') THEN typ:=3;
      If (Type_Element='TORS') THEN typ:=4;
      Case Typ of
               {Beam element}
        1:     For j:=1 to Nel do begin
               i:=i+1;
               if(i>Maxc) then Erreur(1);
               ElType[i]:=1;
               ReadLn(Donnees,i1,i2,E,G,Dext,Dint,Rho);
               Numero_interne(i1,i2);
               Inertie:=Pi*(Power(Dext,4)-Power(Dint,4))/64;
               EI^[i]:=E*Inertie;
               GJ^[i]:=2*G*Inertie;
               Long:=sqr(x[1][i2]-x[1][i1])+sqr(x[2][i2]-x[2][i1])
                      +sqr(x[3][i2]-x[3][i1]);
               Long:=sqrt(Long);
               Longueur^[i]:=Long;
               Section:=0.25*Pi*(Dext*Dext-Dint*Dint);
               ES^[i]:=E*Section;
               elmass^[i]:=Rho*Section*Long;
               Roll_inertia^[i]:=2*Inertie*Rho*Long;
               np[1][i]:=i1;
               np[2][i]:=i2;
               end;

               {Beams element}
         2:    For j:=1 to Nel do begin
               i:=i+1;
               if(i>Maxc) then Erreur(1);
               ElType[i]:=2;
               ReadLn(Donnees,i1,i2,E,G,Dext,Dint,Rho,Ky);
               Numero_interne(i1,i2);

               { Ky = 0.9 for a full circular section.
                 Ky = 0.5 for a thin circular section. }

               Inertie:=Pi*(Power(Dext,4)-Power(Dint,4))/64;
               EI^[i]:=E*Inertie;
               GJ^[i]:=2*G*Inertie;   (* Inertie de torsion *)
               Long:=sqr(x[1][i2]-x[1][i1])+sqr(x[2][i2]-x[2][i1])
                      +sqr(x[3][i2]-x[3][i1]);
               Long:=sqrt(Long);
               Longueur^[i]:=Long;
               Section:=0.25*Pi*(Dext*Dext-Dint*Dint);
               ES^[i]:=E*Section;
               elmass^[i]:=Rho*Section*Long;
               Roll_inertia^[i]:=2*Inertie*Rho*Long;
               if (Ky>0) then
               PhiY^[i]:=12*EI^[i]/(Ky*G*Section*Long*Long)
               else PhiY^[i]:=0;
               np[1][i]:=i1;
               np[2][i]:=i2;
               end;

               {Bar element}
         3:    For j:=1 to Nel do begin
               i:=i+1;
               if(i>Maxc) then Erreur(1);
               ElType[i]:=3;
               ReadLn(Donnees,i1,i2,E,Dext,Dint,Rho);
               Numero_interne(i1,i2);
               Long:=sqr(x[1][i2]-x[1][i1])+sqr(x[2][i2]-x[2][i1])
                      +sqr(x[3][i2]-x[3][i1]);
               Long:=sqrt(Long);
               Longueur^[i]:=Long;
               Section:=0.25*Pi*(Dext*Dext-Dint*Dint);
               ES^[i]:=E*Section;
               elmass^[i]:=Rho*Section*Long;
               np[1][i]:=i1;
               np[2][i]:=i2;
               end;

               {Torsion bar element}
         4:    For j:=1 to Nel do begin
               i:=i+1;
               if(i>Maxc) then Erreur(1);
               ElType[i]:=4;
               ReadLn(Donnees,i1,i2,G,Dext,Dint,Rho);
               Numero_interne(i1,i2);
               Inertie:=Pi*(Power(Dext,4)-Power(Dint,4))/64;
               GJ^[i]:=2*G*Inertie;
               Long:=sqr(x[1][i2]-x[1][i1])+sqr(x[2][i2]-x[2][i1])
                      +sqr(x[3][i2]-x[3][i1]);
               Long:=sqrt(Long);
               Longueur^[i]:=Long;
               Section:=0.25*Pi*(Dext*Dext-Dint*Dint);
               Roll_inertia^[i]:=2*Inertie*Rho*Long;
               np[1][i]:=i1;
               np[2][i]:=i2;
               end;

        else   begin    {unknown element}
                 Write(Type_element,' ');
                 Erreur(10)
               end
      end  (* Case *)
    end  (* if *)
  until (Nel=999);
  Nel:=i;
  GotoXY(34,WhereY);
  WriteLn('-----------------');
  GotoXY(19,WhereY);
  WriteLn('Total: ',Nel:10,' elements');
End;

(**********************************************************)
PROCEDURE Partitionner(Nodes      :Integer;
                       Var Unknown:Vecteur_Byte_long;
                       Var Global :Vecteur_Ent);

Var i,ig,ddl,suppress,last:Integer;

BEGIN

{ ig = global index for DOF = (node-1)*6 + 1
  unknown[ig]= Unknown number after partition corresponding to ig.
  Global[inc]= Number ig corresponding to unknown number inc.       }

Last:=0;
Suppress:=Free; (* Free défini dans Points_Nodaux, visible ici *)
For i:=1 to Nodes do
    Begin
    ig:=(i-1)*6;
    For ddl:=1 to 6 do
        begin
        ig:=ig+1;
        if (idl[ddl][i]=0) then
                     begin
                     Last:=Last+1;
                     Unknown[ig]:=Last;
                     Global[Last]:=ig;
                     end
                 else
                     begin
                     suppress:=suppress+1;
                     Unknown[ig]:=suppress;
                     Global[suppress]:=ig;
                     end;
     end; (* ddl *)
end; (* i *)
End; (* Partitionner *)
(**********************************************************)
PROCEDURE Change_Base(alpha,beta,gamma:arg_reel;
                      Var Lambda:Matrix_3x3);
BEGIN
Lambda[1][1]:=alpha;
Lambda[2][1]:=beta;
Lambda[3][1]:=gamma;
if ((Abs(alpha)<Macheps) and (Abs(beta)<Macheps)) then
   begin
   Lambda[1][2]:=0;
   Lambda[2][2]:=1;
   Lambda[3][2]:=0;
   if(Gamma>0) then Lambda[1][3]:=-1
               else Lambda[1][3]:=1;
   Lambda[2][3]:=0;
   Lambda[3][3]:=0;
   end
else
   begin
   Lambda[1][2]:=-Beta;
   Lambda[2][2]:=alpha;
   Lambda[3][2]:=0;

   Lambda[1][3]:=-alpha*gamma;
   Lambda[2][3]:=-Beta*Gamma;
   Lambda[3][3]:=alpha*alpha+Beta*beta;
   end;
End;

{**********************************************************}
PROCEDURE Assemble_K(Var K:pMatc;
                     Var ES,EI,GJ,Longueur,PhiY:pVect;
                     Nel:Integer;
                     Elnode:Connexion;
                     ElType:Vecteur_Byte_court);

Var elk :Block12x12;
    Lambda,LambdaT,Bloc,LambdaK   :Matrix_3x3;
    long,Phi,Flex,tension,torsion,
    alpha,beta,gamma,
    a0,a1,a2,a3,a4                :arg_reel;
    i,j,ielm,n1,n2,i1,i2,ib,jb,
    ielk,jelk,
    ddl,iblk,jblk,ideb,jdeb       :Integer;
    g1,inc1                       :Array[1..12] of Integer;


begin
a0:=12;
For ielm:=1 to Nel do
    begin

    For i:=1 to 12 do
        For j:=1 to 12 do elk[j][i]:=0.;

    Case Eltype[ielm] of
    1 : Begin
        long:=Longueur^[ielm];
        tension:=ES^[ielm]/long;
        Flex:=EI^[ielm]/Power(Long,3);
        Torsion:=GJ^[ielm]/Long;
        a1:=a0*Flex;
        a2:=6*long*Flex;
        a3:=2*long*long*Flex;

        elk[1][1]:=tension;
        elk[7][1]:=-tension;

        elk[2][2]:=a1;
        elk[6][2]:=a2;
        elk[8][2]:=-a1;
        elk[12][2]:=a2;

        elk[3][3]:=a1;
        elk[5][3]:=-a2;
        elk[9][3]:=-a1;
        elk[11][3]:=-a2;

        elk[4][4]:=torsion;
        elk[10][4]:=-torsion;

        elk[5][5]:=2*a3;
        elk[9][5]:=a2;
        elk[11][5]:=a3;

        elk[6][6]:=2*a3;
        elk[8][6]:=-a2;
        elk[12][6]:=a3;

        elk[7][7]:=tension;

        elk[8][8]:=a1;
        elk[12][8]:=-a2;

        elk[9][9]:=a1;
        elk[11][9]:=a2;

        elk[10][10]:=Torsion;

        elk[11][11]:=2*a3;

        elk[12][12]:=2*a3;
        end;

    2 : Begin
        long:=Longueur^[ielm];
        Phi:=PhiY^[ielm];
        tension:=ES^[ielm]/long;
        Flex:=EI^[ielm]/(Power(Long,3)*(1+Phi));
        Torsion:=GJ^[ielm]/Long;
        a1:=a0*Flex;
        a2:=6*long*Flex;
        a4:=long*long*Flex;
        a3:=(2-Phi)*a4;
        a4:=(4+Phi)*a4;

        elk[1][1]:=tension;
        elk[7][1]:=-tension;

        elk[2][2]:=a1;
        elk[6][2]:=a2;
        elk[8][2]:=-a1;
        elk[12][2]:=a2;

        elk[3][3]:=a1;
        elk[5][3]:=-a2;
        elk[9][3]:=-a1;
        elk[11][3]:=-a2;

        elk[4][4]:=torsion;
        elk[10][4]:=-torsion;

        elk[5][5]:=a4;
        elk[9][5]:=a2;
        elk[11][5]:=a3;

        elk[6][6]:=a4;
        elk[8][6]:=-a2;
        elk[12][6]:=a3;

        elk[7][7]:=tension;

        elk[8][8]:=a1;
        elk[12][8]:=-a2;

        elk[9][9]:=a1;
        elk[11][9]:=a2;

        elk[10][10]:=Torsion;

        elk[11][11]:=a4;

        elk[12][12]:=a4;

        end;
        (* Bar Element tension/compression *)
  3 :   Begin
        long:=Longueur^[ielm];
        tension:=ES^[ielm]/long;
        elk[1][1]:=tension;
        elk[7][1]:=-tension;
        elk[7][7]:=tension;
        end;
        (* Torsion Bar Element *)
  4 :   Begin
        long:=Longueur^[ielm];
        Torsion:=GJ^[ielm]/Long;
        elk[4][4]:=torsion;
        elk[10][4]:=-torsion;
        elk[10][10]:=Torsion;
        end;
end; (* Case *)

    (* Number of  ending nodes *)
    n1:=elnode[1][ielm];
    n2:=elnode[2][ielm];

    (* Number of unknowns in full partitionned matrix *)

    For i:=1 to 6 do
        begin
        i1:=6*(n1-1)+i;g1[i]:=i1;
        inc1[i]:=unknown[i1];
        end;
    For i:=7 to 12 do
        begin
        ddl:=i-6;
        i2:=6*(n2-1)+ddl;g1[i]:=i2;
        inc1[i]:=unknown[i2];
        end;

    (* Calculate the change base matrix *)
    alpha:=(Xcoord[1][n2]-Xcoord[1][n1])/Long;
    Beta:=(Xcoord[2][n2]-Xcoord[2][n1])/Long;
    gamma:=(Xcoord[3][n2]-Xcoord[3][n1])/Long;

    Change_base(alpha,beta,gamma,Lambda);

    Transpose_Matrix_3x3(Lambda,LambdaT);

    (* Partition matrix elk into 3x3-Blocks and loop on Blocks *)

    FOR Iblk :=1 to 4 do
        begin
        ideb:=(Iblk-1)*3+1;
        FOR Jblk :=Iblk to 4 do
            begin
            jdeb:=(Jblk-1)*3+1;
            (* Extract a 3x3-block, of indexes ib,jb *)
            ib:=0; (* line index *)
            For ielk:=ideb to ideb+2 do
                begin
                ib:=ib+1;
                jb:=0; (* column index *)
                For jelk:=jdeb to jdeb+2 do
                    begin
                    jb:=jb+1;
                    Bloc[jb][ib]:=elk[jelk][ielk];
                    end; (* jelk *)
                end; (* ielk *)

            (* Left-Multiply by Lambda transposed *)
            Multiplie_Matrix_3x3(LambdaT,Bloc,LambdaK);
            Multiplie_Matrix_3x3(LambdaK,Lambda,Bloc);

    (* Loop on elements in block *)
    For ib:=1 to 3 do
        For jb:=1 to 3 do
           begin
           ielk:=(Iblk-1)*3+ib;
           jelk:=(Jblk-1)*3+jb;
           if ((inc1[ielk]<=Free)and(inc1[jelk]<=Free)) then
              begin
              i1:=inc1[ielk];
              i2:=inc1[jelk];
              K^[i2][i1]:=K^[i2][i1]+Bloc[jb][ib];
              end; (* if *)
           end; (* jb,ib *)

       end; (* Jblk *)
       end; (* Iblk *)

    end; (* ielm *)

End; (* assemble_K *)

(**********************************************************)
PROCEDURE Assemble_M(Var M:pMatc;
                     Var Elmass,Roll_Inertia,Longueur : pVect;
                     Nel:Integer;
                     Elnode:Connexion;
                     Eltype:Vecteur_Byte_court);

Var elm :Block12x12;
    Lambda,LambdaT,Bloc,LambdaM   :Matrix_3x3;
    long,mass,Inertie,tempo,
    alpha,beta,gamma,
    a0,a1,a2,a3,a4,a5,a6,a7       :arg_reel;
    i,j,ielm,n1,n2,i1,i2,ib,jb,
    ielk,jelk,
    ddl,iblk,jblk,ideb,jdeb       :Integer;
    g1,inc1                       :Array[1..12] of Integer;


begin
a0:=1/420.;
a1:=140.*a0;
a2:=156.*a0;
a3:=22.*a0;
a4:=54.*a0;
a5:=13.*a0;
a6:=4.*a0;
a7:=3*a0;

For ielm:=1 to Nel do
    begin

    For i:=1 to 12 do
        For j:=i to 12 do elm[j][i]:=0;
  Case ElType[ielm] of
  1..2 : Begin                  (* Full beam elements *)
    long:=Longueur^[ielm];
    mass:=elmass^[ielm];
    Inertie:=Roll_Inertia^[ielm];
    Tempo:=a1*mass;
    elm[1][1]:=tempo;
    elm[7][7]:=Tempo;
    elm[7][1]:=0.5*tempo;
    Tempo:=a2*mass;
    elm[2][2]:=Tempo;
    elm[3][3]:=Tempo;
    elm[8][8]:=Tempo;
    elm[9][9]:=Tempo;
    Tempo:=a1*Inertie;
    elm[4][4]:=Tempo;
    elm[10][4]:=0.5*Tempo;
    elm[10][10]:=Tempo;
    Tempo:=a3*mass*Long;
    elm[6][2]:=Tempo;
    elm[5][3]:=Tempo;
    elm[12][8]:=-Tempo;
    elm[11][9]:=-Tempo;
    Tempo:=a4*mass;              (* 54/420 x M *)
    elm[8][2]:=Tempo;
    elm[9][3]:=Tempo;
    Tempo:=a6*mass*Long*Long;    (* 4/420 x M x L2 *)
    elm[5][5]:=Tempo;
    elm[6][6]:=Tempo;
    elm[11][11]:=Tempo;
    elm[12][12]:=Tempo;
    Tempo:=a5*mass*Long;         (* 13/420 x M x L *)
    elm[9][5]:=Tempo;
    elm[8][6]:=Tempo;
    elm[12][2]:=-Tempo;
    elm[11][3]:=-Tempo;
    Tempo:=a7*mass*Long*Long;
    elm[11][5]:=-Tempo;
    elm[12][6]:=-Tempo;
    end;
 3: begin                        (* Bar Element *)
    long:=Longueur^[ielm];
    mass:=elmass^[ielm];
    Tempo:=a1*mass;
    elm[1][1]:=tempo;
    elm[7][7]:=Tempo;
    elm[7][1]:=0.5*tempo;
    end;
 4: begin                      (* Torsion Bar Element *)
    long:=Longueur^[ielm];
    Inertie:=Roll_Inertia^[ielm];
    Tempo:=a1*Inertie;
    elm[4][4]:=Tempo;
    elm[10][4]:=0.5*Tempo;
    elm[10][10]:=Tempo;
    end;
 end; (* Case *)
  (* Nø des noeuds extr‚mit‚s *)
     n1:=elnode[1][ielm];
     n2:=elnode[2][ielm];

  (* Number of unknowns in full partitionned matrix *)

    For i:=1 to 6 do
        begin
        i1:=6*(n1-1)+i;g1[i]:=i1;
        inc1[i]:=unknown[i1];
        end;
    For i:=7 to 12 do
        begin
        ddl:=i-6;
        i2:=6*(n2-1)+ddl;g1[i]:=i2;
        inc1[i]:=unknown[i2];
        end;

    (* Calculate the change base matrix *)
    alpha:=(Xcoord[1][n2]-Xcoord[1][n1])/Long;
    Beta:=(Xcoord[2][n2]-Xcoord[2][n1])/Long;
    gamma:=(Xcoord[3][n2]-Xcoord[3][n1])/Long;

    Change_base(alpha,beta,gamma,Lambda);

    Transpose_Matrix_3x3(Lambda,LambdaT);

    (* Partition matrix elm into 3x3-Blocks and loop on Blocks *)

    FOR Iblk :=1 to 4 do
        begin
        ideb:=(Iblk-1)*3+1;
        FOR Jblk :=Iblk to 4 do
            begin
            jdeb:=(Jblk-1)*3+1;
            (* Extract a 3x3-block, of indexes ib,jb *)
            ib:=0; (* indice de ligne *)
            For ielk:=ideb to ideb+2 do
                begin
                ib:=ib+1;
                jb:=0; (* indice de colonne *)
                For jelk:=jdeb to jdeb+2 do
                    begin
                    jb:=jb+1;
                    Bloc[jb][ib]:=elm[jelk][ielk];
                    end; (* jelk *)
                end; (* ielk *)

            (* Left-Multiply by Lambda transposed *)
            Multiplie_Matrix_3x3(LambdaT,Bloc,LambdaM);
            Multiplie_Matrix_3x3(LambdaM,Lambda,Bloc);

    (* Loop on elements in block *)
    For ib:=1 to 3 do
        For jb:=1 to 3 do
           begin
           ielk:=(Iblk-1)*3+ib;
           jelk:=(Jblk-1)*3+jb;
           if ((inc1[ielk]<=Free)and(inc1[jelk]<=Free)) then
              begin
              i1:=inc1[ielk];
              i2:=inc1[jelk];
              M^[i2][i1]:=M^[i2][i1]+Bloc[jb][ib];
              end; (* if *)
           end; (* jb,ib *)

       end; (* Jblk *)
       end; (* Iblk *)

    end; (* ielm *)

End; (* assemble_M *)
(******************************************************************)
PROCEDURE Add_Local_Masses(M:pMatc);

Var i,j,ddl          : Byte;
    Numero_Noeud,global,Noeud,
    Nombre_Masses    : Integer;
    Mass             : arg_reel;

Begin
ReadLn(Donnees,Nombre_Masses);
if(Nombre_Masses<>999) then
   begin
   For i:=1 to Nombre_Masses do
    begin
    ReadLn(Donnees,Noeud,ddl,Mass);
    Numero_Noeud :=Num_Interne(Noeud);
    global:=(Numero_Noeud-1)*6+ddl;
    j:=unknown[global];
    if j>Free then
       writeLn
     ('* Warning * Mass added to fixed node',Noeud:5,' DOF ',ddl:1,'!');
    M^[j][j]:=M^[j][j]+Mass;
    end;
   end; (* if *)
Close(Donnees);
END;
(******************************************************************)
PROCEDURE Add_Springs(K:pMatc);


Var i,j1,j2,ddl1,ddl2   : Byte;
    Noeud1,Noeud2,
    g1,g2,g3,
    Nombre_Ressorts     : Integer;
    Raideur             : arg_reel;

Begin
ReadLn(Donnees,Nombre_Ressorts);
if(Nombre_Ressorts<>999) then
   begin
   For i:=1 to Nombre_Ressorts do
    begin
    ReadLn(Donnees,Noeud1,ddl1,Noeud2,ddl2,Raideur);
    Noeud1:=Num_Interne(Noeud1);
    Noeud2:=Num_Interne(Noeud2);
    g1:=(Noeud1-1)*6+ddl1;
    g2:=(Noeud2-1)*6+ddl1;
    if(g2<g1) then
              begin
              g3:=g1;
              g1:=g2;
              g2:=g3;
              end;
    j1:=unknown[g1];
    j2:=unknown[g2];
    if (j1>Free)and (j2>Free) then
        begin
        writeLn('* Warning * Spring added between two fixed DOF!');
        writeLn(Noeud1:5,ddl1:2,Noeud2:5,ddl2:2,Raideur:10:3);
        end
        else
        begin
        if(j1<free) then K^[j1][j1]:=K^[j1][j1]+Raideur;
        if(j2<free) then K^[j2][j2]:=K^[j2][j2]+Raideur;
        if(j1<Free)and(j2<Free) then K^[j1][j2]:=K^[j1][j2]-Raideur;
        end;
    end;(* i *)
   end; (* if *)
END;
(******************************************************************)
PROCEDURE Table_of_Modes;
Var Modes_a_sortir,Lines,i: Byte;

  Procedure En_Tete;
  begin
    TextOut(CrtDC,24,14,'Mode',4);
    TextOut(CrtDC,100,14,'Frequency',9);
    TextOut(CrtDC,208,14,'Generalized Mass',16);
    GotoXY(2,3); Writeln('---------------------------------------------');
  end;

Begin
  En_tete;
  Modes_a_sortir:=Free;
  For i:=1 to Free do
  begin
    Write(i:6); Write('     '); disp_real(Frequence^[i]);
    Write('       '); disp_real(Masse_Generalisee^[i]); Writeln;
    Lines:=WhereY;
    If (Lines=28) then
    begin
      Modes_a_sortir:=Modes_a_sortir-18;
      if (Modes_a_sortir>0) then
      begin
        Pause(Msg);
        ClrScr;
        En_Tete;
      end;
    end;
  end; (* i *)
  WriteLn(' ---------------------------------------------');
  Pause(Msg);
  ClearScreen
End;
{**********************************************************}
PROCEDURE Frequences_et_Masses_Generalisees;
{Calculate frequencies and generalized masses of modes}
VAR i,j :Integer;

BEGIN
  New(omega2);
  New(Frequence);
  New(Masse_Generalisee);
  New(U);
  New(V);

  For j:=1 to Free do          {Loop on Modes}
  begin
    omega2^[j]:=d^[j];
    if(omega2^[j]<0) then
    begin
      writeLn('     * ATTENTION * Mode ',j:2,' Negative eigenvalue!');
      omega2^[j]:=0;
    end; (*if*)
    Frequence^[j]:=0.5*sqrt(omega2^[j])/Pi;
    for i:=1 to Free do U^[i]:=Phi^[i][j];
    Matvec(Free,M,U,V,-1);
    Produit_Scalaire(Free,U,V,masse_generalisee^[j]);
  end; (* j *)
END; {Frequences_et_Masses_Generalisees}
(**********************************************************************)
PROCEDURE Shape_of_Modes;
Label again;
Var  s : array[0..20] of char;

    PROCEDURE Ecrire_EnTete;
    Begin
      gotoxy(1,3);
      WriteLn(' Mode #',j:2,'   Eigenvalue=   ',Omega2^[j]:16);
      Writeln('            Frequency=        ',frequence^[j]:16:8,' Hz');
      Writeln('            Generalized Mass= ',masse_generalisee^[j]:16:8);
      WriteLn;
      WriteLn(' Eigenvector:');
      WriteLn(' ------------');
      Lines:=WhereY;
      GotoXY(2,Lines);Write('# Node');
      GotoXY(13,Lines);Write('u');
      GotoXY(24,Lines);Write('v');
      GotoXY(35,Lines);Write('w');
      GotoXY(44,Lines);Write('theta x');
      GotoXY(55,Lines);Write('theta y');
      GotoXY(66,Lines);Writeln('theta z');
      Writeln
    end;

BEGIN     (* Forme_des_Modes *)
  Repeat
    again: ClearScreen;
    WriteLn;
    Writeln(' Input number of desired mode (99 to exit) ');
    GotoXY(6,4); Write(' Mode number? ');
    Read_An_Integer(j);
    if(j=99) then
    begin
      ClearScreen;
      exit
    end;
    if (j > Free) or (j < 1) then goto again;

    ClrScr;
    Str(j:2,St);
    StrPCopy(s,'Shape of Mode '+St);
    TextOut(CrtDC,7,2,s,strlen(s));
    Ecrire_EnTete;
    Noeuds_a_sortir:=Nodes;
    For i:= 1 to Nodes do
    begin
      iglob:=(i-1)*6;
      For ddl:=1 to 6 do
      begin
        iglob:=iglob+1;
        if (idl[ddl][i]=0) then
        begin
          index:=Unknown[iglob];
          sortie[ddl]:=Phi^[index][j];
        end
        else sortie[ddl]:=0;
      end; (* ddl *)
      write(Numero[i]:5);
      For ddl:=1 to 6 do Write(' ',sortie[ddl]:10:4);
      writeLn;
      Lines:=WhereY;
      if (Lines=28) Then
      begin
        Noeuds_a_sortir:=Noeuds_a_sortir-13;
        if (Noeuds_a_sortir>0) then
        begin
          Pause(Msg);
          ClrScr;
          Ecrire_EnTete
        end
      end;
    End; (* i *)
    Pause(Msg);
    ClearScreen;
  until (j=99);
END; (* Shape_of_Modes *)
(**********************************************************************)
PROCEDURE Draw_Mode;
Label again;
Var  s,s1 : array[0..20] of char;
     miniy,maxiy, x,xi,xf : ar_reel;
     nddl : INTEGER;
BEGIN
  Repeat
    again: ClearScreen;
    WriteLn;
    Writeln(' Input number of desired mode (99 to exit) ');
    GotoXY(6,4); Write(' Mode number? ');
    Read_An_Integer(j);
    if j<>99 then
      Repeat
        GotoXY(6,5); Write(' DOF number? ');
        Read_An_Integer(nddl)
      Until (nddl>0) AND (nddl<7);
    if(j=99) then
    begin
      ClearScreen;
      exit
    end;
    if (j > Free) or (j < 1) then goto again;

    ClrScr;
    Str(j:2,St);
    StrPCopy(s,'Shape of Mode '+St);
    Str(nddl:2,St);
    StrPCopy(s1,'ddl '+St);
    Noeuds_a_sortir:=Nodes;

    {seek minimum/maximum values of mode}
    miniy:=0; maxiy:=0;
    For i:= 1 to Nodes do
    begin
      iglob:=(i-1)*6;
      For ddl:=1 to 6 do
      begin
        iglob:=iglob+1;
        if (idl[ddl][i]=0) then
        begin
          index:=Unknown[iglob];
          sortie[ddl]:=Phi^[index][j];
          if sortie[ddl] < miniy then miniy:=sortie[ddl];
          if sortie[ddl] > maxiy then maxiy:=sortie[ddl]
        end
        else sortie[ddl]:=0;
      end; { ddl }
    end; {for i}

    {initialize mode graph}
    xi:=0; xf:=Nodes;
    if nddl<4 then
      InitFenetre(CrtDC,10,xi,xf,-1.0,1.0)
    else
      InitFenetre(CrtDC,10,xi,xf,miniy,maxiy);

    {Draw mode}
    For i:= 1 to Nodes do
    begin
      iglob:=(i-1)*6;
      For ddl:=1 to 6 do
      begin
        iglob:=iglob+1;
        if (idl[ddl][i]=0) then
        begin
          index:=Unknown[iglob];
          sortie[ddl]:=Phi^[index][j];
        end
        else sortie[ddl]:=0;
      end; { ddl }

      x:=int(i);
      if i=1 then MoveXY(CrtDC,x,sortie[nddl])
             else LineXY(CrtDC,x,sortie[nddl]);
    end; {for i}
    Legendes(CrtDC,s,'N° node',s1);
    SortieGraphique
  until (j=99);
END; (* Draw_Mode
*******************************************************************)
PROCEDURE Write_Logo;
var j : integer; ch: char;
Begin
  ClrScr;
  Cadre(CrtDC,5,5,315,450);
  TextOut(CrtDC,35,132,'EF3D FINITE ELEMENTS PROGRAM',28);
  GotoXY(4,11);
  Write('--------------------------------');
  GotoXY(4,12);
  Write('   Release TPW 3.0, May 2006   ');
  GoToXY(3,27);
  Write('(C)Copyright J.-P. Dumont 1987,1991');
  GoToXY(3,28);
  Write('             J.-P. Moreau 2006');
  GotoXY(45,30); write('Any key to continue...');
  ch:=readkey;
End;
(******************************************************************)
PROCEDURE Solve_New_Problem;
VAR  answer: char;
BEGIN

  Points_Nodaux(Nodes,Numero,Xcoord,Idl);
  WriteLn('     Number of unknowns= ',Free);
  WriteLn('     Nombre of imposed displacements= ',Prescribed);
  WriteLn('     Total number of DOF= ',6*Nodes);
  WriteLn;

  Partitionner(Nodes,Unknown,General);

  Elements(Nel,ES,EI,GJ,Longueur,elmass,Roll_Inertia,PhiY,Xcoord,elnode,eltype);
  WriteLn;

  (* Initializations *)
  n:=Free;
  New(M);
  New(Phi);
  New(K);

  For j:= 1 to n do
    for i:=1 to n do
    begin
      K^[i,j]:=0;
      M^[i,j]:=0;
    end;

  Assemble_K(K,ES,EI,GJ,Longueur,PhiY,Nel,elnode,ElType);

  Add_Springs(K);

  Assemble_M(M,elmass,Roll_Inertia,Longueur,Nel,elnode,ElType);

  Add_Local_Masses(M);

  Writeln('     Assembling Ok.');


  (* Make matrices symmetric *)
  For i:= 1 to free do
    for j:=1 to i-1 do
    begin
      K^[j,i]:=K^[i,j];
      M^[j,i]:=M^[i,j];
    end;

  Writeln; Writeln;
  Write('     Do you want to save the results (y/n)? ');
  answer:=Readkey;
  GotoXY(1,WhereY);ClrEol;
  If (answer='y') then
  begin
    {Save general data}
    Open_Save_Write('GEN');
    WITH General_Data do
    begin
      Rec_Problem:=Probleme;
      Rec_Titre:=Titre;
      Rec_Nel:=Nel;
      Rec_Nodes:=Nodes;
      Rec_Free:=Free;
      Rec_Prescribed:=Prescribed;
      Rec_Vecteur_Byte_court:=ElType;
      Rec_Vecteur_Byte_Long:=Unknown;
      Rec_ddl_Status:=idl;
      Rec_Connexion:=ElNode;
      Rec_Class:=Numero;
      Rec_Vecteur_Entier:=General;
      Rec_Vecteur_Colonne1:=ES^;
      Rec_Vecteur_Colonne2:=EI^;
      Rec_Vecteur_Colonne3:=GJ^;
      Rec_Vecteur_Colonne4:=PhiY^;
      Rec_Vecteur_Colonne5:=Longueur^;
      Rec_Vecteur_Colonne6:=ElMass^;
      Rec_Vecteur_Colonne7:=Roll_Inertia^;
      Rec_Vecteur_Position:=Xcoord
    end;
    write(Dump,General_Data);
    Close(Dump);

    {Save stiffness matrix before call Reduc1 which
     destroyes the main diagonal of K }

    Open_Save_Write('STF');
    FOR i:=1 to Free do
    begin
      For j:=1 to Free do VCol^[j]:=K^[j][i];
      write(Dump1,VCol^);
    end; {i}
    Close(Dump1);

    {Save mass matrix}
    Open_Save_Write('MAS');
    FOR i:=1 to Free do
    begin
      For j:=1 to Free do VCol^[j]:=M^[j][i];
      write(Dump1,VCol^);
    end; {i}
    Close(Dump1);
    Writeln;
    Writeln('     End of saving K and M matrices.')
  End; {If}

  Reduc1(Free,K,M,dl,error);
  if (error) then Erreur(8);

  Tred2(Free,tol,K,Phi,d,e);
  Dispose(K);

  Tql2(Free,tol,d,e,Phi,error);
  if error then Erreur(9);

  Rebaka(Free,1,Free,M,dl,Phi);
{
  Normalize to unity eigenvectors on biggest translation
  or, by default, on biggest rotation.
}
  For Mode:=1 to Free do
  begin
    Biggest:=0;
    For i:=1 to Free do
    begin
        ddl:= General[i] mod 6;
        if (ddl=0) then ddl:=6;
        if (ddl<4) and (Abs(Phi^[i][mode])> Abs(Biggest))
                   then Biggest:=Phi^[i][Mode];
    end; (* i *)
    if (Abs(Biggest)=0) then   (* All DOF fixed in translation *)
    begin
        For i:=1 to Free do
        begin
            ddl:= General[i] mod 6;
            if(ddl=0) then ddl:=6;
            if (ddl>3) and (Abs(Phi^[i][mode])> Abs(Biggest)) then
                         Biggest:=Phi^[i][Mode];
            end; (* i *)
        end; (* if *)
    For i:=1 to Free do Phi^[i][Mode]:=Phi^[i][Mode]/Biggest
  end;

  Frequences_et_Masses_Generalisees;

  IF (answer='y') then
  begin
    {Sauvegarde de la matrice modale}
    Open_Save_Write('MOD');
    write(Dump1,Omega2^,Frequence^,Masse_Generalisee^);
    FOR i:=1 to Free do
    begin
      For j:=1 to Free do VCol^[j]:=Phi^[j][i];
      write(Dump1,VCol^)
    end;
    Close(Dump1);
  end;  {If}

END; (* Solve_New_Problem *)
(*********************************************************************)
PROCEDURE Test_Orthogonality;
Var ch: char;
BEGIN
  ClrScr;
  WriteLn;
  Biggest:=0;
  n:=1;
  index:=1;
  For j:=1 to Free do
  begin
    for i:=1 to Free do U^[i]:=Phi^[i][j];
    Matvec(Free,M,U,V,-1);
    For l:=j+1 to Free do
    begin
      for i:=1 to Free do U^[i]:=Phi^[i][l];
      Produit_Scalaire(Free,U,V,produit);
      if Abs(Produit) > Biggest then
      begin
        n:=j;
        index:=l;
        biggest:=Abs(Produit);
      end
    end (* l *)
  end; (* j *)
  TextOut(CrtDC,5,5,'TEST OF ORTHOGONALITY:',22);
  Writeln;
  WriteLn('     Greatest crossed product: ',Biggest:10);
  Pause(Msg);
  ClearScreen
END; (* Test_Orthogonality *)

{************************************************************************
Main program }
BEGIN
  WinCrtInit(' EF3D');
  SetTextColor(CrtDC,RGB(0,0,255));
  ExitProc:=@ErrorHandler;
  Write_Logo;
  SaveSw:=False;
  Statsw:=False;
  ModeSw:=True;
  OrthoSW:=false;
  probleme:='';
debut:
  {initialize column vectors}
  New(dl);New(d); New(e); New(ES); New(EI); New(GJ); New(PhiY); New(VCol);
  New(Longueur); New(elmass); New(Roll_Inertia);
  ClrScr;
  jj:=0;
  if probleme='' then
  begin
    GotoXY(1,2);
    Write('     Is it a new problem (y/n)? ');
    answ:=readkey; if jj=6 then answ:=readkey;
    if answ='y' then Newsw:=TRUE else Newsw:=FALSE;
    Open_Input_File
  end
  else
  begin
    fichier:=probleme+'.DAT';
    Assign(donnees,fichier);
    {$I-} Reset(donnees); {$I+}
    if IOResult <> 0 then Erreur(14)
  end;

  StrPCopy(Msg,'  '+Probleme+' - '+Messg);

CASE Newsw OF

True: {Case of a new problem}
      begin
        Clrscr;
        ReadLn(Donnees,Titre);
        GotoXY(1,2);
        WriteLn('     '+Titre);
        WriteLn;
        Solve_New_Problem
      end;

False: {Case continue former calculation}
      begin

      {Read general data}
      Open_Save_Read('GEN');
      Read(Dump,General_Data);

      WITH General_Data do
      begin
        Probleme:=Rec_Problem;
        Titre:=Rec_Titre;
        Nel:=Rec_Nel;
        Nodes:=Rec_Nodes;
        Free:=Rec_Free;
        Prescribed:=Rec_Prescribed;
        ElType:=Rec_Vecteur_Byte_court;
        Unknown:=Rec_Vecteur_Byte_long;
        idl:=Rec_ddl_Status;
        ElNode:=Rec_Connexion;
        Numero:=Rec_Class;
        General:=Rec_Vecteur_Entier;
        ES^:=Rec_Vecteur_Colonne1;
        EI^:=Rec_Vecteur_Colonne2;
        GJ^:=Rec_Vecteur_Colonne3;
        PhiY^:=Rec_Vecteur_Colonne4;
        Longueur^:=Rec_Vecteur_Colonne5;
        ElMass^:=Rec_Vecteur_Colonne6;
        Roll_Inertia^:=Rec_Vecteur_Colonne7;
        Xcoord:=Rec_Vecteur_Position;
      end; {With}
      Close(Dump);

      {read stiffness matrix}
      Open_Save_Read('STF');
      New(K);
      FOR i:=1 to Free do
      begin
        Read(Dump1,VCol^);
        For j:=1 to Free do K^[j][i]:=VCol^[j]
      end; {i}
      Close(Dump1);

      {read mass matrix}
      Open_Save_Read('MAS');
      New(M);
      FOR i:=1 to Free do
      begin
        Read(Dump1,VCol^);
        For j:=1 to Free do M^[j][i]:=VCol^[j];
      end; {i}
      Close(Dump1);

      {read frequencies and generalized masses}
      Open_Save_Read('MOD');
      New(Omega2);
      New(Frequence);
      New(Masse_Generalisee);
      New(Phi);
      New(U);  { These two tables are required}
      New(V);  { for orthogonality test       }

      Read(Dump1,Omega2^,Frequence^,Masse_Generalisee^);

      {read modal matrix}
      FOR i:=1 to Free do
      begin
        Read(Dump1,VCol^);
        For j:=1 to Free do Phi^[j][i]:=VCol^[j]
      end; {i}
      Close(Dump1);

      WriteLn;
      WriteLn('     End of reading saved results.');
      Pause(Msg);
    End {False}
END; {CASE}

Repeat
   Clrscr;
   Writeln;
   Writeln;
   TextOut(CrtDC,100,15,'MENU OF RESULTS',15);
   Cadre(CrtDC,40,35,280,160);  {frame of menu}
   for i:=1 to 7 do
   begin
     gotoxy(8,i+3);
     writeln(Menu_Item[i],':',i:2)
   end;
   Writeln;
   Write('            Option (1 to 7): '); read(jj);
   Case jj of
     1: Begin    {Données générales }
          Clrscr;
          Writeln;
          TextOut(CrtDC,20,20,'GENERAL DATA:',13);
          Writeln;
          Writeln;
          WriteLn('     Number of unknowns= ',Free);
          WriteLn('     Number of imposed displacements= ',Prescribed);
          WriteLn('     Total number of DOF= ',6*Nodes);
          WriteLn;
          Writeln('     Number of elements= ',Nel);
          Pause(Msg)
        end;
     2: Begin    {Table of modes }
          Clrscr;
          Table_of_Modes
        end;
     3: Begin    {Shape of modes }
          Shape_of_Modes
        end;
     4: Begin    {Draw a mode}
          Draw_mode
        end;
     5: Begin    {Test of orthogonality of modes }
          Test_Orthogonality
        end;
     6: Begin    {Exit section}
          {free all column vectors}
          New(dl);Dispose(d); Dispose(e); Dispose(ES); Dispose(EI); Dispose(GJ);
          Dispose(PhiY); Dispose(VCol); Dispose(Longueur); Dispose(elmass);
          Dispose(Roll_Inertia); Dispose(omega2); Dispose(frequence); Dispose(U);
          Dispose(V); Dispose(masse_generalisee);
          {libération des matrices M, Phi}
          Dispose(M); Dispose(Phi);
          DoneWinCrt
        End;
     7: begin    {Other Problem}
          Probleme:='';
          {free all column vectors}
          New(dl);Dispose(d); Dispose(e); Dispose(ES); Dispose(EI); Dispose(GJ);
          Dispose(PhiY); Dispose(VCol); Dispose(Longueur); Dispose(elmass);
          Dispose(Roll_Inertia); Dispose(omega2); Dispose(frequence); Dispose(U);
          Dispose(V); Dispose(masse_generalisee);
          {libération des matrices M, Phi}
          Dispose(M); Dispose(Phi);
          Goto debut;
        end
   end;
Until jj=6 {normal exit}
END.

{end of file ef1.pas}
