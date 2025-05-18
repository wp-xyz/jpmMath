{*****************************************************
*  General Constants and Types Declarations used by  *
*  Program EF1.pas.                                  *
*****************************************************}
UNIT DECLARA1;

INTERFACE
USES WinTypes;

CONST
  Macheps= 1E-15;
  MaxSize=  4096;
  Size   =  1024;
  npoin  =   512;
  Taille_Bloc=128;
  maxr   =     1;
  maxc   =    89;  { Max. number of unknowns: Turbo=104 ,Turbo-87=89}
  maxndp =    89;  { Max. number of nodal points}
  maxddl =   534;  { 6*Maxc, Max. number of DOF}
  neq    =     1;
  nsol   =     1;
  ndiag  =     2;
  IBM    = false;

{  Maxc is chosen not to exceed the 64 ko limit for Tables:
   with Turbo   : 104*104*6 = 64896 < 65536
        Turbo-87:  90* 90*8 = 64800 < 65536       }

TYPE
  arg_reel        = Double;
  Matrice         = Array[1..maxr, 1..maxc] of arg_reel;
  pMatc           = ^Matrice_carree;
  Matrice_carree  = Array[1..maxc, 1..maxc] of arg_reel;
  Matrice_n_mPlus1= Array[1..neq,0..ndiag] of arg_reel;
  Matrice_nxr     = Array[1..neq,1..nsol] of arg_reel;
  Matrix_3x3      = Array[1..3,1..3] of arg_reel;
  Block12x12      = Array[1..12,1..12] of arg_reel;
  ddl_status      = Array[1..6,1..maxndp] of Byte;
  connexion       = Array[1..2,1..maxc] of Byte;
  Classe          = Array[1..Maxndp] of Integer;

  Vecteur_ligne   = array[1..maxr] of arg_reel;
  pVect           = ^Vecteur_colonne;
  Vecteur_colonne = array[1..maxc] of arg_reel;
  Vecteur_Position= array[1..3,1..maxndp] of arg_reel;
  Vecteur_Reel    = Array[1..Size] of arg_reel;
  Vecteur_Ent     = Array[1..Maxc] of Integer;
  Vecteur_Byte_long = Array[1..Maxddl] of Byte;
  Vecteur_Byte_court= Array[1..Maxc]   of Byte;

  Mot_Cle         = String[6];
  Descriptif      = String[20];
  Enregistrement  = String[80];
  Nom_Fichier     = String[40];
  Suffixe         = String[3];

  Gen_Data_Record = Record
                      Rec_Problem : Nom_Fichier;
                      Rec_Titre   : Enregistrement;
                      Rec_Nel,Rec_Nodes,Rec_Free,Rec_Prescribed : Integer;
                      Rec_Vecteur_Byte_long : Vecteur_Byte_long;
                      Rec_Vecteur_Byte_court: Vecteur_Byte_court;
                      Rec_ddl_Status : ddl_status;
                      Rec_Connexion  : Connexion;
                      Rec_Class      : Classe;
                      Rec_Vecteur_Entier : Vecteur_Ent;
                      Rec_Vecteur_Colonne1 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne2 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne3 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne4 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne5 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne6 : Vecteur_Colonne;
                      Rec_Vecteur_Colonne7 : Vecteur_Colonne;
                      Rec_Vecteur_Position : Vecteur_Position;
                    End;

 VAR  CrtDC : HDC;  {Device Context of main window}

 IMPLEMENTATION
   {No procedures}
 END.

