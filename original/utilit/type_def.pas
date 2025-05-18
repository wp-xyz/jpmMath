{*******************************************************
*             UNIT Type_Def.pas                        *
* ---------------------------------------------------- *
* Definition of types and global variables used by any *
* program with 2D graphs.                              *
*                                   J-P Moreau         *
*                                (www.jpmoreau.fr)     *
*******************************************************}
UNIT Type_Def;

INTERFACE
Uses WinTypes;        {for HDC}

Const
  Macheps = 1.2E-16;  {1E-12 if no coprocessor}
  Size   =  2048;     {this value can be modified by user}

Type
  real_ar         = DOUBLE;   { REAL if no coprocessor }
  complex         = Array[1..2] of REAL;

  Real_Vector     = Array[1..Size] of REAL;
  Integer_Vector  = Array[1..Size] of Integer;
  Complex_Vector  = Array[1..Size] of Complex;

  RV=^Real_Vector;
  CV=^Complex_Vector;

  Descriptif      = String[20];
  FileName        = String[40];
  Title           = String[50];

VAR
  MaxX,MaxY : integer;             { used by par WGRAPH_2D.PAS }
  Log_X,Log_Y,Ech_auto : boolean;
  X_mini,X_maxi,Y_mini,Y_maxi : real_ar;
  CrtDC : HDC;

IMPLEMENTATION

Begin    {By default: no log. scales and automatic mini, maxi}
  Log_X:=false; Log_Y:=false; Ech_auto:=true;
End.

{End of unit type_def.pas}

