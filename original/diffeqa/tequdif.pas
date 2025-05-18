{***********************************************************
*     TEST PROGRAM OF PROCEDURE ODEINT OF UNIT EQUDIF      *
*       (Runge-Kutta Method with Time Step Control)        *
* -------------------------------------------------------- *
* SAMPLE RUN:                                              *
* Integrate ODE First Order System:                        *
*   Y1' = Y2                                               *
*   Y2' = (3-C)*Y1+2*Y4-q1                                 *
*   Y3' = Y4                                               *
*   Y4' = -2*Y2-C*Y3                                       *
*   Y5' = Y6                                               *
*   Y6' = -(1+C)*Y5                                        *
* from t1=0 to t2=10 with additional data:                 *
*   Y start = (0,0,0,0,1,0)                                *
*   Starting integration step = 0.01                       *
*   Desired precision = 1E-8                               *
*   Minimum integration step = 0.001                       *
*   C=1, q1=2                                              *
*                                                          *
* At time =    10.0000000000000                            *
* Y(1) = -0.121038E+04                                     *
* Y(2) = -0.908048E+03                                     *
* Y(3) =  0.116298E+04                                     *
* Y(4) =  0.870537E+03                                     *
* Y(5) = -0.496866E-02                                     *
* Y(6) = -0.141420E+01                                     *
*                                                          *
* Final time step = 2.175370624046558E-002                 *
*                                                          *
*                       TPW Release By J-P Moreau, Paris.  *
*                               (www.jpmoreau.fr)          *
***********************************************************}
Program Tequdif;

Uses WinCrt1, Type_def, Equdif, Utils1;

Var
    Y: vecteur_etat;
    i,nbad,nok,nvar: integer;
    h,t1,t2: REAL_AR;
    s,s1:String;

BEGIN
  nvar:=6;
  For i:=1 to nvar do Y[i]:=0.0;
  Y[5]:=1.0;
  t1:=0;
  t2:=10;
  h:=0.01;

  odeint1(Y,nvar,t1,t2,1E-8,h,0.001,nok,nbad);

  writeln;
  Line1; Line0;
  RMsg('  At time = ', t2);
  For i:=1 to nvar do
  begin
    Str(i,s1); s:='  Y('+s1+') = ';
    RMsg(s, Y[i])
  end;
  Line0;
  RMsg('  Final time step =', h);
  Line0; Line1;
  Readln;
  DoneWinCrt
END.

{end of file tequdif.pas}