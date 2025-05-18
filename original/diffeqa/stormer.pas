{************************************************************
*   Differential equation y"=f(x,y,y') by Stormer's method  *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
* Integrate y" = 8yy / (1 + 2x) from x=0 to x=1,            *
* with initial conditions: x(0)=0, y(0)=1 and y'(0)=-2      *
* and compare with exact solution: y = 1 / (1 + 2x)         *
*                                                           *
* Output file (stormer.lst):                                *
*                                                           *
* --------------------------------------------------------- *
*   Differential equation y"=f(x,y,y') by Stormer's method  *
* --------------------------------------------------------- *
*       X          Y        Y exact        Error            *
*    0.0000    1.000000    1.000000      0.00000000         *
*    0.0100    0.980392    0.980392      0.00000000         *
*    0.0200    0.961538    0.961538      0.00000000         *
*    0.0300    0.943396    0.943396      0.00000003         *
*    0.0400    0.925926    0.925926      0.00000008         *
*    0.0500    0.909091    0.909091      0.00000016         *
*    ...       ...         ...           ...                *
*    0.9500    0.344866    0.344828      0.00003818         *
*    0.9600    0.342505    0.342466      0.00003889         *
*    0.9700    0.340176    0.340136      0.00003962         *
*    0.9800    0.337878    0.337838      0.00004035         *
*    0.9900    0.335612    0.335570      0.00004109         *
*    1.0000    0.333375    0.333333      0.00004183         *
*                                                           *
*  End of file.                                             *
*                             Pascal Version By J-P Moreau  *
*                                  (www.jpmoreau.fr)        *
************************************************************}
PROGRAM Stormer;
Uses WinCrtMy;

TYPE
        Table = Array[1..4] of DOUBLE;

VAR
        c, x, y, z   : Table;
        h, a1, a2, a3: DOUBLE;
        xend,yex,er  : DOUBLE;
        k            : INTEGER;
        fp           : TEXT;

   
   {Example: y' := 8 yy / (1 + 2x) }
   Function F(x,y,z:DOUBLE): DOUBLE;
   begin
     F := z
   end;

   Function G(x,y,z:DOUBLE): DOUBLE;
   begin
     G := 8.0*y*y/(1.0+2.0*x)
   end;

   {exact solution} 
   Function Fx(x:DOUBLE): DOUBLE;
   begin
     Fx := 1.0/(1.0+2.0*x)
   end;

   Procedure RK4D2(x,y,z,h:DOUBLE;VAR x1,y1,z1:DOUBLE);
   Var c1,d1,c2,d2,c3,d3,c4,d4:DOUBLE;
   Begin
     c1:=F(x,y,z);
     d1:=G(x,y,z);
     c2:=F(x+h/2.0,y+h/2*c1,z+h/2.0*d1);
     d2:=G(x+h/2.0,y+h/2*c1,z+h/2.0*d1);
     c3:=F(x+h/2.0,y+h/2*c2,z+h/2.0*d2);
     d3:=G(x+h/2.0,y+h/2*c2,z+h/2.0*d2);
     c4:=F(x+h,y+h*c3,z+h*d3);
     d4:=G(x+h,y+h*c3,z+h*d3);
     x1:=x+h;
     y1:=y+h*(c1+2*c2+2*c3+c4)/6.0;
     z1:=z+h*(d1+2*d2+2*d3+d4)/6.0
   End;


BEGIN
    
  h:=0.01;
  a1:=1.08333333333333;
  a2:=-2.0*(a1-1.0);
  a3:=a1-1.0;
  
  Assign(fp,'stormer.lst'); Rewrite(fp);
  writeln(fp,'-----------------------------------------------------------------');
  writeln(fp,'  Differential equation y"=f(x,y,y'') by Stormer''s method       ');
  writeln(fp,'-----------------------------------------------------------------');
  writeln;
  writeln(' Input x begin, x end, y begin, y'' begin:');
  writeln;
  write('    x begin = '); read(x[1]);
  write('    x end   = '); read(xend);
  write('    y begin = '); read(y[1]);
  write('    y'' begin = '); read(z[1]);
  
  yex:=Fx(x[1]); er:=0.0;
  writeln(fp,'       X          Y        Y exact        Error');
  writeln(fp,x[1]:10:4,y[1]:12:6,yex:12:6,'  ',er:14:8);
  {call Runge-Kutta for first 2 steps}
  for k:=1 to 2 do
  begin
    RK4D2(x[k],y[k],z[k],h,x[k+1],y[k+1],z[k+1]);
    yex:=Fx(x[k+1]); er:=abs(yex-y[k+1]);
    writeln(fp,x[k+1]:10:4,y[k+1]:12:6,yex:12:6,'  ',er:14:8);
  end;
  {main Stormer loop}
  while x[3] < xend  do
  begin
    for k:=2 to 4 do
      c[k]:=G(x[5-k],y[5-k],z[5-k]);
    y[4]:=2.0*y[3]-y[2]+h*h*(a1*c[2]+a2*c[3]+a3*c[4]);
    x[4]:=x[3]+h; yex:=Fx(x[4]); er:=abs(yex-y[4]);
    writeln(fp,x[4]:10:4,y[4]:12:6,yex:12:6,'  ',er:14:8);
    for k:=1 to 3 do
    begin
      x[k]:=x[k+1]; y[k]:=y[k+1]
    end
  end;
  writeln(fp);
  writeln(fp,'  End of file.');
  close(fp);
  writeln;
  writeln('  Results in file stormer.lst');
  Readkey; Donewincrt

END.

{end of file stormer.pas}