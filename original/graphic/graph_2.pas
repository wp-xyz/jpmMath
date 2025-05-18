{***************************************************
*    2D procedures used by program apollo.pas      *
*                                                  *
*                         by Jean-Pierre Moreau    *
*                           (www.jpmoreau.fr)      *
***************************************************}
Unit Graph_2;

Interface
Uses Type_Def,    {for real_ar}
     Maths_2D;    {for real_arSwap }

  Type  coordonnee  = real_ar;
        Point2D     = Record
                        x,y: coordonnee
                      End;
        Transform2D = Array[1..2,1..3] of coordonnee;
        Vecteur2D   = Point2D;

  Const Orig2D  : point2D     = (x:0; y:0);
        Eps     : coordonnee  = 1E-10;
        Erreur2D: boolean     = false;
        Id2D    : Transform2D = ((1,0,0),(0,1,0));

  Function  Norme2D(V: vecteur2D): real_ar;
  Function  Dist2D(ptA,ptB: point2D): real_ar;
  Procedure GetVecteur2D(ptA,ptB: point2D; VAR V: vecteur2D);

Implementation

  {return the norm (length) of a vector of type
   vecteur2D defined by its traces along Ox and Oy}
  Function Norme2D;
  Begin
    Erreur2D:=false;
    With V do
    begin
      x:=abs(x); y:=abs(y);
      if x < y then real_arSwap(x,y);
      if x > eps then norme2D := x*sqrt(1+sqr(y/x))
                 else norme2D := 0
    end
  End;

  {define a vector by its traces (projections) along
   Ox and Oy from the end coordinates. Returns the
   type vecteur2D}
  Procedure GetVecteur2D;
  Begin
    Erreur2D:=false;
    V.x:=ptB.x-ptA.x;
    V.y:=ptB.y-ptA.y
  End;

  {return the distance between two points of type point2D}
  Function Dist2D;
  Var V: vecteur2D;
  Begin
    Erreur2D:=false;
    GetVecteur2D(ptA,ptB,V);
    Dist2D:=Norme2D(V)
  End;

END.

{end of file graph_2.pas}