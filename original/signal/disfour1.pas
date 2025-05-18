{*********************************************************
* This program calculates the Fourier coefficients of a  *
* periodic discrete function F(x) by using the procedure *
* DiscreetFourierHn of unit Fourier with a graph of the  *
* reconstructed signal from Fourier coefficients.        *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* F(x) has the following shape:                          *
*                                                        *
*  1 ---- 2  pi/4                                        *
*    !  !                                                *
*    !  !                                                *
* ------!-------------------------------------------> x  *
*  -pi  !    0         ! pi                              *
*       !              !                                 *
*     3 ---------------- 4 -pi/4                         *
*                                                        *
* Calculate the Fourier coefficients of a periodic       *
* discrete function:                                     *
*                                                        *
* Number of points:  4                                   *
*                                                        *
*  x1 = -3.1415927                                       *
*  y1 = 0.7853982                                        *
*                                                        *
*  x2 = -1.5707963                                       *
*  y2 = 0.7853982                                        *
*                                                        *
*  x3 = -1.5707963                                       *
*  y3 = -0.7853982                                       *
*                                                        *
*  x4 = 3.1415927                                        *
*  y4 = -0.7853982                                       *
*                                                        *
* Lowest  harmonic: 0                                    *
* Highest harmonic: 125                                  *
*                                                        *
* A graph is displayed (reconstructed curve F(x) with    *
* 125 harmonics).                                        *
*                                                        *
*                     TPW version by J-P Moreau, Paris.  *
*                            (www.jpmoreau.fr)           *
**********************************************************
Note: Exact Fourier coefficients for this periodic function are:
      an=0 if n is even, -1/2n if n=4p+1, 1/2n if n=4p+3
      bn=0 if n=4p, -1/2n if n=4p+1, 1/n if n=4p+2, -1/2n if n=4p+3
------------------------------------------------------------------}
Program Test_DiscreetFourier;
Uses WinCrtmy,Type_def,Graph_2D,Fourier;

Var
    h1,h2,i,j,npoints,mpoints: Integer;
    Xi,Yi,Ai,Bi: pTab;
    Y: RV;
    a,b,dx,x,x0,xf: Double;

BEGIN
  New(Xi); New(Yi); New(Ai); New(Bi); New(Y);
  WinCrtInit(' Fourier coefficients of a periodic discreet function');
  writeln;
  write(' Number of points: '); readln(npoints);
  writeln;
  for i:=1 to npoints do
  begin
    write('  x',i,'= '); readln(Xi^[i]);
    write('  y',i,'= '); readln(Yi^[i]);
    writeln
  end;
  write(' Lowest  harmonic: '); readln(h1);
  write(' Highest harmonic: '); readln(h2);
  for i:=h1 to h2 do
  begin
    DiscreetFourierHn(npoints,Xi,Yi,i,a,b);
    Ai^[i+1] := a;
    Bi^[i+1] := b
  end;
  {reconstruct signal from Fourier series coefficients and show result}
  mpoints:=256; x0:=Xi^[1]; xf:=Xi^[npoints]; dx:=(xf-x0)/(mpoints-1);
  x:=x0-dx;
  for i:=1 to mpoints do
  begin
    x:=x+dx; Y^[i]:=0;
    for j:=h1 to h2 do
    begin
      if j=0 then
        Y^[i]:=Y^[i]+Ai^[1]
      else
      begin
        om:=2*PI*j/T;
        Y^[i]:=Y^[i]+Ai^[j+1]*cos(om*x)+Bi^[j+1]*sin(om*x)
      end
    end
  end;
  {draw curve F(x) }
  ClrScr;
  CourbeXY(CrtDc,mpoints,10,Y,x0,xf);
  Legendes(CrtDc,'Reconstructed Signal F(x)','X','Y');
  SortieGraphique;
  Dispose(Xi); Dispose(Yi); Dispose(Ai); Dispose(Bi); Dispose(Y);
  DoneWinCrt
END.

{end of file disfour1.pas}