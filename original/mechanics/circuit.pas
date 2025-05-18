{*********************************************************
*                    PROGRAM CIRCUIT                     *
* Use Runge-Kutta method to solve a LRC circuit or equi- *
* valent damped mass-spring problem.                     *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* (Solve mass-spring system with:                        *
*  F=75 N, M=50 kg, K=100 N (C=0.01), D=0.05)            *
*                                                        *
* Specify [o]scillating or [n]o oscillating term...      *
* o                                                      *
* Give V or F, L or M, R or D, C or 1/K (4L/C-R^2 > 0):  *
* 75 50 0.05 0.01                                        *
* one period at t=   4.442879    s                       *
* Give t_final and number of points:                     *
* 10 10                                                  *
*                                                        *
*     time        q or d    i or speed  analytic q or d  *
*   1.000000      0.7504      1.4993      0.6331         *
*   2.000000      2.2485     -0.0017      1.4628         *
*   3.000000      0.7483     -2.9955      1.0890         *
*   4.000000     -2.2440      0.0055      0.1436         *
*   5.000000      0.7555      5.9850      0.2223         *
*   6.000000      6.7320     -0.0150      1.1913         *
*   7.000000      0.7351    -11.9580      1.4147         *
*   8.000000    -11.2021      0.0378      0.5163         *
*   9.000000      0.7878     23.8922      0.0131         *
*  10.000000     24.6302     -0.0915      0.7537         *
*                                                        *
* (For a better accuracy, increase number of points).    *
*                                                        *
* ------------------------------------------------------ *
* Reference: "Problem Solving with Fortran 90 By David   *
*             R. Brooks, Springer-Verlag New York, 1997" *
*                                                        *
*                  Pascal version by J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)            *
*********************************************************} 
Program Circuit;
Uses WinCrt;
 

{------------------------------------}
      FUNCTION AofT(x,v,D,K,M,F:REAL):REAL;
{ Calculate acceleration for "mass and spring" problem.}
      Begin
        AofT:=-(-F+D*v+K*x)/M
      End;
{---------------------------------------------}
      PROCEDURE MassAndSpring(VAR x,v:REAL; D,K,M,F,dt:REAL);
{ Calculate motion for mass and spring problem with constant force term.}
      Var
          k1_x,k2_x,k3_x,k4_x,k1_v,k2_v,k3_v,k4_v: REAL;
      Begin
{ Runge-Kutta coefficients...}
        k1_x:=v;
        k2_x:=v+AofT(x,v,D,K,M,F)*dt/2.0;
        k3_x:=v+AofT(x,v,D,K,M,F)*dt/2.0;
        k4_x:=v+AofT(x,v,D,K,M,F)*dt;
        k1_v:=AofT(x,v,D,K,M,F);
        k2_v:=AofT(x+v*dt/2.,v+k1_v*dt/2.,D,K,M,F);
        k3_v:=AofT(x+v*dt/2.,v+k2_v*dt/2.,D,K,M,F);
        k4_v:=AofT(x+v*dt,v+k3_v*dt,D,K,M,F);
{ Propagate solution...}
        x:=x+(k1_x+2.0*k2_x+2.0*k3_x+k4_x)*dt/6.0;
        v:=v+(k1_v+2.0*k2_v+2.0*k3_v+k4_v)*dt/6.0
      End; {MassAndSpring}


{main program
*********************************************************
* Use Runge-Kutta method to solve LRC circuit problems. *
* Variable equivalences with mass-and-spring problem:   *
* V   => force F                                        *
* q   => displacement x                                 *
* i   => velocity v                                     *
* L   => mass m                                         *
* R   => damping constant D                             *
* 1/C => spring constant K                              *
********************************************************}
VAR
      i,q,L,C,R,V,t,dt,t_final: REAL;
      k1_i,k2_i,k3_i,k4_i: REAL;
      j,n: INTEGER;
      choice: CHAR;
BEGIN

{ Choose circuit type... }

      Writeln;
      Write(' Specify [o]scillating or [n]o oscillating term... ');
      Readln(choice);
      if (choice='o') or (choice='O') then
      begin
{ Ld^2q/dt^2+Rdq/dt+q/C=V }
        Repeat
          Writeln(' Give V or F, L or M, R or 1/K, C or D (4L/C-R^2 > 0):');
          Write(' '); Readln(V,L,R,C);
          Writeln(' one period at t=',4*pi*L/SQRT(4.0*L/C-R*R),' s');
          Writeln(' Give t_final and number of points:');
          Write(' '); Readln(t_final,n)
        Until 4.0*L/C-R*R > 0.0;
        q:=0.0; i:=C*V*R/2./L; t:=0.0;  {Initial values}
        dt:=t_final/n;
        Writeln;
        Writeln('       time        q or d    i or speed  analytic q or d');
        For j:=1 to n do
        begin
          MassAndSpring(q,i,R,1./C,L,V,dt);
          t:=t+dt;
          Writeln(' ',t:12:6,q:12:4,i:12:4,C*V*(1.0-EXP(-R*t/2.0/L)*COS(SQRT(4.0*L/C-R*R)*t/2.0/L)):12:4)
        end
      end
      else if (choice='n') or (choice='N') then
      begin
{ Ldi/dt+Ri=V (no oscillating term)...}
        Writeln(' Give V or F, L or M, R or 1/K:');
        Write(' '); Readln(V,L,R);
        Writeln;
        Writeln(' time constant at t=',L/R,' s');
        Writeln(' Give t_final and number of points:');
        Write(' '); Readln(t_final,n);
        i:=0.0; t:=0.0;   {Initial values}
        dt:=t_final/n;
        Writeln;
        Writeln('       time      i or speed  analytic i');
        For j:=1 to n do
        begin
{ Runge-Kutta coefficients...}
          k1_i:=(V-R*i             )/L;
          k2_i:=(V-R*(i+k1_i*dt/2.0))/L;
          k3_i:=(V-R*(i+k2_i*dt/2.0))/L;
          k4_i:=(V-R*(i+k3_i*dt)   )/L;
{ Propagate solution...}
          i:=i+(k1_i+2.0*k2_i+2.0*k3_i+k4_i)*dt/6.0;
          t:=t+dt;
          Writeln(' ',t:12:6,i:12:4,V/R*(1.-EXP(-R*t/L)):12:4)
        end
      end
      else
        Writeln(' No such choice.  Try again...');
      Readkey; DoneWinCrt

END.

{end of file circuit.pas}