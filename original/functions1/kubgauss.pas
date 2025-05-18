{ ----------------------- UNIT kubgauss.pas ------------------------- *
 *    "Numerical Algorithms with C, By Gisela Engeln-Muellges          *
 *     and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].             *
 *                                                                     *
 *                           Pascal unit version By J-P Moreau, Paris. *
 *                                       (www.jpmoreau.fr)             *
 ----------------------------------------------------------------------}
UNIT Kubgauss;

Interface

  {visible from calling program}
  Function Kub3GauN (
                     Px, Py: Double;
                     Qx, Qy: Double;
                     Rx, Ry: Double;
                     n, m: Integer;
                 Var Wert: Double
                    ): Integer;

Implementation

Const NMAX = 25;
{***********************************************************************
*       Here the constants for cubature over rectangles using          *
*       Newton-Cotes formulas are provided.                            *
*                                                                      *
*       The structure kub4_ta contains a pair of doubles:              *
*         a node t and the corresponding weight a.                     *
*                                                                      *
*       The arrays K_i (i = 0..7) contain (i+1) such pairs             *
*         according to the order of the cubature formula.              *
*                                                                      *
*       Author:        Uli Eggermann, 02.17.1991 (C version)           *
***********************************************************************}
Type

     kub4_ta = Record      {node t und weight a}
       t, a: Double
     End;

     VEC = Array[0..NMAX] of Double;

     Tripel = Record       {weights, x-, y-coordinates}
       w, x, y: Double
     End;

Var
     K_0: kub4_ta;
     K_1: array[0..1] of kub4_ta;
     K_2: array[0..2] of kub4_ta;
     K_3: array[0..3] of kub4_ta;
     K_4: array[0..4] of kub4_ta;
     K_5: array[0..5] of kub4_ta;
     K_6: array[0..6] of kub4_ta;
     K_7: array[0..7] of kub4_ta;

{**********************************************************************
* Constants for n-point Gaussian cubature:                            *
**********************************************************************}
     Gau1Konst: Tripel;
     Gau2Konst: Array[0..1] of Tripel;
     Gau3Konst: Array[0..2] of Tripel;
     Gau7Konst: Array[0..6] of Tripel;


  {user defined function(x,y) }
  Function func(x,y:Double):Double;
  begin
    func := Exp(-(x*x+y*y))
  end;

Function Kub4GauE (
                    a, b: Double; Nx: Integer;
                    c, d: Double; Ny: Integer;
                    Verf: Integer;
                Var Wert: Double;
                    Schaetzen: Integer;
                Var FehlerSch: Double
                  ): Integer;
{***********************************************************************
* Cubature over rectangles using Newton-Cotes formulas.                *
*                                                                      *
* Integrate the function  f(x,y) using the summed cubature formula     *
* of Gauss for the rectangle (a,b; c,d).                               *
* The edges of the sub-rectangles have the lengths:  (b-a) / Nx   and  *
* (d-c) / Ny .                                                         *
*                                                                      *
* Input parameters:                                                    *
*  REAL a, b, Nx     left, right x-end points, number of intervals     *
*  REAL c, d, Ny     left, right y-end points, number of intervals     *
*  integer Verf      order of method  (0 <= Verf <= 7)                 *
*  REAL f ()         function                                          *
*  REAL Wert         value of integral                                 *
*  integer Schaetzen if nonzero : find estimate                        *
*  REAL FehlerSch    error estimate, obtained if desired by a second   *
*                    cubature pass using half of the step size         *
*                                                                      *
* Return value :                                                       *
*   0:               o.k.                                              *
*   1:               Nx improper                                       *
*   2:               Ny improper                                       *
*   3:               order incorrect                                   *
*   4:               Integration interval has length zero              *
*                                                                      *
* Author             Uli Eggermann, 03.31.1996                         *
***********************************************************************}
Label return;
Var
  i, j, k, kend, u, v, Nab, Ncd: Integer;
  Hab, Hcd, Wert1: Double;
  w,x,y,z: Double;

Begin

  Wert1 := 0.0;
  if Nx < 1 then
  begin
    Kub4GauE := 1;
    goto return
  end;
  if Ny < 1
  then
  begin
    Kub4GauE := 2;
    goto return
  end;
  if (Verf < 0) or (Verf > 7) then
  begin
    Kub4GauE := 3;
    goto return
  end;
  if (a = b) or (c = d) then
  begin
    Kub4GauE := 4;
    goto return
  end;

  if Schaetzen<>0 then kend:=2 else kend:=1;
  for k := 1 to kend do
  begin
    Nab := k * Nx;                            {number of intervals}
    Ncd := k * Ny;
    Hab := 0.5 * (b - a) / Nab;               {half of x-step size}
    Hcd := 0.5 * (d - c) / Ncd;               {half of y-step size}

    Wert := 0.0;                              {initialize}

    for i := 0 to Nab-1 do
      for j := 0 to Ncd-1 do
        for u := 0 to Verf do
          for v := 0 to Verf do
          begin
            Case Verf of
              0:begin
                w := Hab * Hcd * K_0.a * K_0.a;
                x := a + Hab * (K_0.t + 2 * i + 1);
                y := c + Hcd * (K_0.t + 2 * j + 1)
                end;
              1:begin
                w := Hab * Hcd * K_1[u].a * K_1[v].a;
                x := a + Hab * (K_1[u].t + 2 * i + 1);
                y := c + Hcd * (K_1[v].t + 2 * j + 1)
                end;
              2:begin
                w := Hab * Hcd * K_2[u].a * K_2[v].a;
                x := a + Hab * (K_2[u].t + 2 * i + 1);
                y := c + Hcd * (K_2[v].t + 2 * j + 1)
                end;
              3:begin
                w := Hab * Hcd * K_3[u].a * K_3[v].a;
                x := a + Hab * (K_3[u].t + 2 * i + 1);
                y := c + Hcd * (K_3[v].t + 2 * j + 1)
                end;
              4:begin
                w := Hab * Hcd * K_4[u].a * K_4[v].a;
                x := a + Hab * (K_4[u].t + 2 * i + 1);
                y := c + Hcd * (K_4[v].t + 2 * j + 1)
                end;
              5:begin
                w := Hab * Hcd * K_5[u].a * K_5[v].a;
                x := a + Hab * (K_5[u].t + 2 * i + 1);
                y := c + Hcd * (K_5[v].t + 2 * j + 1)
                end;
              6:begin
                w := Hab * Hcd * K_6[u].a * K_6[v].a;
                x := a + Hab * (K_6[u].t + 2 * i + 1);
                y := c + Hcd * (K_6[v].t + 2 * j + 1)
                end;
              7:begin
                w := Hab * Hcd * K_7[u].a * K_7[v].a;
                x := a + Hab * (K_7[u].t + 2 * i + 1);
                y := c + Hcd * (K_7[v].t + 2 * j + 1)
                end
            End;
            z := func(x, y);
            Wert := Wert + w * z
          end;
    if (Schaetzen<>0) and (k = 1) then Wert1 := Wert;    {store value}
  end;

  if Schaetzen<>0 then                                   {estimate}
    FehlerSch := (Wert - Wert1) / 3.0;

  Kub4GauE := 0;

Return: end;

{ -------------------------------------------------------------------- }
Function Kub4GauV (
                   x: VEC; Nx: Integer;
                   y: VEC; Ny: Integer;
                   Verf: Integer;
               Var Wert: Double;
                   Schaetzen: Integer;
               Var FehlerSch: Double
                  ): Integer;
{***********************************************************************
* Cubature over rectangles using  Newton-Cotes formulas.               *
*                                                                      *
* Integrate the function f(x,y) using the summed Gaussian cubature     *
* formula on the rectangle (a,b) x (c,d).                              *
*                                                                      *
* The edge lengths for the sub-rectangles are not identical, however,  *
* as in Kub4GauE, but are given by thr two vectors  X []  and  Y [] .  *
*                                                                      *
* Parameters:                                                          *
*  REAL X []         vector of x-interval end points:                  *
*                      a = X[0] < X[1] < .. < X[Nx] = b                *
*  REAL Y []         vector of y-interval end points:                  *
*                      c = Y[0] < Y[1] < .. < Y[Ny] = d                *
*  integer Verf      order of method  (0 <= Verf <= 7)                 *
*  integer Schaetzen 0: no error estimation                            *
*                    1: error estimation using a second pass of        *
*                       cubature for half the step size                *
*  REAL Wert         value for integral                                *
*  REAL FehlerSch    error estimate of  Wert                           *
*                                                                      *
* Return value :                                                       *
* --------------                                                       *
*   0:               o.k.                                              *
*   1:               order number improper                             *
*   2:               x-interval of length zero                         *
*   3:               y-interval of length zero                         *
*                                                                      *
* Author:         Uli Eggermann, 03.31.1991 (C version)                *
***********************************************************************}
Label return;
Var
    k, kend, Si, Sj, i, j, u, v: Integer;
    Hx, Hy, Wert1, Tx, Ty: Double;
Begin

  Wert1 := 0.0;
  if (Verf < 0) or (Verf > 7) then                { check valid order}
  begin
    Kub4GauV := 1;
    goto return
  end;
  for i := 0 to Nx-1 do                           { X-interval test  }
    if x[i+1] <= x[i] then
    begin
      Kub4GauV := 2;
      goto return
    end;
  for j := 0 to Ny-1 do                           { Y-interval test  }
    if y[j+1] <= y[j] then
    begin
      Kub4GauV := 3;
      goto return
    end;

  {#define KUB KubArr[Verf] }

  if Schaetzen<>0 then kend:=2 else kend:=1;
  for k := 1 to kend do                           { with estimate ?  }
  begin
    Wert:=0.0;                                    { initialize       }
    for i := 0 to Nx-1 do                         { X-intervals      }
    begin
      Hx := (x[i+1] - x[i]) / (2*k);
      Si:=1;                                      { halve X-direction}
      While Si<2*k do
      begin
        Tx := x[i] + Si * Hx;                     { X-interval center}
        for j := 0 to Ny-1 do                     { Y-intervals      }
        begin
          Hy := (y[j+1] - y[j]) / (2*k);          { halve Y-direction}
          Sj:=1;
          While Sj<2*k do
          begin
            Ty := y[j] + Sj * Hy;                 { Y-interval center}
            for u := 0 to Verf do
              for v := 0 to Verf do
              begin
                Case Verf of
                  0: begin
                       Wert := Wert + Hx * K_0.a * Hy * K_0.a *
                       func(Tx + Hx * K_0.t, Ty + Hy * K_0.t)
                     end;
                  1: begin
                       Wert := Wert + Hx * K_1[u].a * Hy * K_1[v].a *
                       func(Tx + Hx * K_1[u].t, Ty + Hy * K_1[v].t)
                     end;
                  2: begin
                       Wert := Wert + Hx * K_2[u].a * Hy * K_2[v].a *
                       func(Tx + Hx * K_2[u].t, Ty + Hy * K_2[v].t)
                     end;
                  3: begin
                       Wert := Wert + Hx * K_3[u].a * Hy * K_3[v].a *
                       func(Tx + Hx * K_3[u].t, Ty + Hy * K_3[v].t)
                     end;
                  4: begin
                       Wert := Wert + Hx * K_4[u].a * Hy * K_4[v].a *
                       func(Tx + Hx * K_4[u].t, Ty + Hy * K_4[v].t)
                     end;
                  5: begin
                       Wert := Wert + Hx * K_5[u].a * Hy * K_5[v].a *
                       func(Tx + Hx * K_5[u].t, Ty + Hy * K_5[v].t)
                     end;
                  6: begin
                       Wert := Wert + Hx * K_6[u].a * Hy * K_6[v].a *
                       func(Tx + Hx * K_6[u].t, Ty + Hy * K_6[v].t)
                     end;
                  7: begin
                       Wert := Wert + Hx * K_7[u].a * Hy * K_7[v].a *
                       func(Tx + Hx * K_7[u].t, Ty + Hy * K_7[v].t)
                     end
                end
              end;
            Inc(Sj,2)
          end
        end;
        Inc(Si,2)
      end
    end;

    if (Schaetzen<>0) and (k = 1) then Wert1 := Wert;   {store value}
  end;

  if Schaetzen<>0 then                                     {estimate}
    FehlerSch := (Wert - Wert1) / 3.0;

  Kub4GauV := 0;

return: end;


{ ------------------------------------------------------------------ }
Function Kub3GauN (
                   Px, Py: Double;
                   Qx, Qy: Double;
                   Rx, Ry: Double;
                   n, m: Integer;
               Var Wert: Double
                  ): Integer;
{***********************************************************************
* Cubature over triangular regions using the n-point Gauss formula     *
*                                                                      *
* Integrate the function f (x,y) over the triangle PQR using the summed*
* n-point Gauss cubature formula on m x m subtriangles.                *
* The subtriangles have edges of length 1/m of the original edge       *
* lengths.                                                             *
*                                                                      *
* Input parameters:                                                    *
*  REAL   Px,Py    coordinates of  P                                   *
*  REAL   Qx,Ry    coordinates of  Q                                   *
*  REAL   Rx,Ry    coordinates of  R                                   *
*  integer n       order of method (= number of points in each sub-    *
*                  triangle)                                           *
*  integer m       number of subtriangles along one edge               *
*                                                                      *
* Output parameter:                                                    *
*  REAL   Wert     value of integral                                   *
*                                                                      *
* Return value :                                                       *
*   0:             o.k.                                                *
*   1:             m  improper                                         *
*   2:             the corners P, Q and R are collinear                *
*   3:             nth order method not implemented                    *
*                                                                      *
* Author:          Uli Eggermann, 8.1.1990 (C version)                 *
***********************************************************************}
Label return;

{ GauNEpsilon serves as a check of collinearity; if the area of P, Q  }
{ and R has area less than Epsilon/2, we judge the three points to be }
{ collinear.                                                          }

Const GauNEpsilon = 0.0001;

Var
    GauNmax: Integer;

    d, i, j, k: Integer;
    Fak, Area, Dx, Dy, X, Y: Double;
    hPQx, hPQy, hPRx, hPRy, Gw, Gx, Gy: Double;

Begin

  GauNmax := 13;

  if m < 1 then                               { m  o.k. ?            }
  begin
    Kub3GauN := 1;
    goto return
  end;

  Area := Px * (Qy - Ry)                      { Test collinearity    }
       + Qx * (Ry - Py)
       + Rx * (Py - Qy);

  if ABS(Area) < GauNEpsilon then
  begin
    Kub3GauN := 2;
    goto return
  end;

  if Not n in [1,2,3,7] then                  { desired method       }
  begin                                       {  implemented ?       }
    Kub3GauN := 3;
    goto return
  end;

  Wert := 0.0;                                { initialize           }
  Area := Area / (m * m);                     { double triangle area }
  hPQx := (Qx - Px) / m;
  hPRx := (Rx - Px) / m;                      { edge vectors for the }
  hPQy := (Qy - Py) / m;                      {   m * m              }
  hPRy := (Ry - Py) / m;                      {   subtriangles       }

  for d := 0 to 1 do                          { types of triangles d }
  begin
    if d<>0 then Fak:=-1.0 else Fak:=1.0;     { d = 1: reflected     }
    for j := d to m-1 do                      { j: along   PR        }
      for i := d to m-j+d-1 do                { i: along   PQ        }
      begin
        Dx := Px + i * hPQx + j * hPRx;       { (Dx,Dy) ist top      }
        Dy := Py + i * hPQy + j * hPRy;       { corner of subtriangle}
        for k := 0 to n-1 do                  { Sum of weighted      }
        begin
          Case n of
            1:begin
                Gw:=Gau1Konst.w;
                Gx:=Gau1Konst.x;
                Gy:=Gau1Konst.y
              end;
            2:begin
                Gw:=Gau2Konst[k].w;
                Gx:=Gau2Konst[k].x;
                Gy:=Gau2Konst[k].y
              end;
            3:begin
                Gw:=Gau3Konst[k].w;
                Gx:=Gau3Konst[k].x;
                Gy:=Gau3Konst[k].y
              end;
            7:begin
                Gw:=Gau7Konst[k].w;
                Gx:=Gau7Konst[k].x;
                Gy:=Gau7Konst[k].y
              end
          end;
          X := Dx + Fak * (Gx * hPQx + Gy * hPRx);
          Y := Dy + Fak * (Gx * hPQy + Gy * hPRy);
          Wert := Wert + Gw * func(X,Y)       { function values      }
        end
      end                                     { (see above for X, Y) }
  end;

  Wert := Wert * Area;                {multiply by double the subarea}
  Kub3GauN := 0;

return: end;


BEGIN

{initialize constants}
  K_0.t:=0.0; K_0.a:=2.0;
  K_1[0].t:= -0.577350269189626; K_1[0].a:=1.0;
  K_1[1].t:=  0.577350269189626; K_1[1].a:=1.0;
  K_2[0].t:= -0.774596669241483; K_2[0].a:=0.5555555555555556;
  K_2[1].t:=  0.0;               K_2[1].a:=0.8888888888888888;
  K_2[2].t:=  0.774596669241483; K_2[2].a:=0.5555555555555556;
  K_3[0].t:= -0.861136311594053; K_3[0].a:=0.347854845137454;
  K_3[1].t:= -0.339981043584856; K_3[1].a:=0.652145154862546;
  K_3[2].t:=  0.339981043584856; K_3[2].a:=0.652145154862546;
  K_3[3].t:=  0.861136311594053; K_3[3].a:=0.347854845137454;
  K_4[0].t:= -0.906179845938664; K_4[0].a:=0.236926885056189;
  K_4[1].t:= -0.538469310105683; K_4[1].a:=0.478628670499366;
  K_4[2].t:=  0.0;               K_4[2].a:=0.5688888888888889;
  K_4[3].t:=  0.538469310105683; K_4[3].a:=0.478628670499366;
  K_4[4].t:=  0.906179845938664; K_4[4].a:=0.236926885056189;
  K_5[0].t:= -0.9324695142031521; K_5[0].a:=0.17132449237917;
  K_5[1].t:= -0.661209386466265; K_5[1].a:=0.360761573048139;
  K_5[2].t:= -0.238619186083197; K_5[2].a:=0.467913934572691;
  K_5[3].t:=  0.238619186083197; K_5[3].a:=0.467913934572691;
  K_5[4].t:=  0.661209386466265; K_5[4].a:=0.360761573048139;
  K_5[5].t:=  0.9324695142031521; K_5[5].a:=0.17132449237917;
  K_6[0].t:= -0.949107912342759; K_6[0].a:=0.12948496616887;
  K_6[1].t:= -0.741531185599394; K_6[1].a:=0.279705391489277;
  K_6[2].t:= -0.405845151377397; K_6[2].a:=0.381830050505119;
  K_6[3].t:=  0.0;               K_6[3].a:=0.417959183673469;
  K_6[4].t:=  0.405845151377397; K_6[4].a:=0.381830050505119;
  K_6[5].t:=  0.741531185599394; K_6[5].a:=0.279705391489277;
  K_6[6].t:=  0.949107912342759; K_6[6].a:=0.12948496616887;
  K_7[0].t:= -0.960289856497536; K_7[0].a:=0.101228536290376;
  K_7[1].t:= -0.7966664774136269; K_7[1].a:=0.222381034453374;
  K_7[2].t:= -0.525532409916329; K_7[2].a:=0.313706645877887;
  K_7[3].t:= -0.18343464249565;  K_7[3].a:=0.362683783378362;
  K_7[4].t:=  0.18343464249565; K_7[4].a:=0.362683783378362;
  K_7[5].t:=  0.525532409916329; K_7[5].a:=0.313706645877887;
  K_7[6].t:=  0.7966664774136269; K_7[6].a:=0.222381034453374;
  K_7[7].t:=  0.960289856497536; K_7[7].a:=0.101228536290376;

  Gau1Konst.w:=0.5; Gau1Konst.x:=1.0/3.0; Gau1Konst.y:=1.0/3.0;
  Gau2Konst[0].w:=0.25; Gau2Konst[0].x:=1.0/6.0; Gau2Konst[0].y:=0.5;
  Gau2Konst[1].w:=0.25; Gau2Konst[1].x:=0.5; Gau2Konst[1].y:=1.0/6.0;
  Gau3Konst[0].w:=1.0/6.0; Gau3Konst[0].x:=1.0/6.0; Gau3Konst[0].y:=1.0/6.0;
  Gau3Konst[1].w:=1.0/6.0; Gau3Konst[1].x:=2.0/3.0; Gau3Konst[1].y:=1.0/6.0;
  Gau3Konst[2].w:=1.0/6.0; Gau3Konst[2].x:=1.0/6.0; Gau3Konst[2].y:=2.0/3.0;
  Gau7Konst[0].w:=0.1125;  Gau7Konst[0].x:=1.0/3.0; Gau7Konst[0].y:=1.0/3.0;
  Gau7Konst[1].w:=0.0661970763942531; Gau7Konst[1].x:=0.4701420641051151; Gau7Konst[1].y:=0.4701420641051151;
  Gau7Konst[2].w:=0.0661970763942531; Gau7Konst[2].x:=0.05971587178976981; Gau7Konst[2].y:=0.4701420641051151;
  Gau7Konst[3].w:=0.0661970763942531; Gau7Konst[3].x:=0.4701420641051151; Gau7Konst[3].y:=0.05971587178976981;
  Gau7Konst[4].w:=0.06296959027241357; Gau7Konst[4].x:=0.1012865073234563; Gau7Konst[4].y:=0.1012865073234563;
  Gau7Konst[5].w:=0.06296959027241357; Gau7Konst[5].x:=0.7974269853530873; Gau7Konst[5].y:=0.1012865073234563;
  Gau7Konst[6].w:=0.06296959027241357; Gau7Konst[6].x:=0.1012865073234563; Gau7Konst[6].y:=0.7974269853530873;

END.

{ ------------------------ END kubgauss.pas -------------------------- }