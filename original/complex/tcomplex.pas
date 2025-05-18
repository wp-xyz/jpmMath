{************************************************
* Test elementary operations on complex numbers *
* --------------------------------------------- *
* Uses: unit complex2.pas.                      *
*                                               *
* SAMPLE RUN:                                   *
*                                               *
* Z1=(1.00,-2.00)                               *
* Z2=(2.50,5.75)                                *
* Z1+Z2=(3.50,3.75)                             *
* Z1+Z2 with modulus and phase:                 *
* (5.1296,0.8199)                               *
* Z1 with modulus and phase:                    *
* (2.2361,-1.1071)                              *
* Z2 with modulus and phase:                    *
* (6.2700,1.1607)                               *
* Z1*Z2=(14.00,0.75)                            *
* Z1*Z2 with modulus and phase:                 *
* (14.0201,0.0535)                              *
* Z1/Z2 with modulus and phase:                 *
* (0.3566,-2.2678)                              *
* Exp(Z1)=(-1.13,-2.47)                         *
* Exp(Z1) with modulus and phase:               *
* (2.7183,-2.0000)                              *
* ArgSh(Z2)=(2.5246,1.1560)                     *
*                                               *
*          Pascal version by J-P Moreau, Paris. *
*                    (www.jpmoreau.fr)          *
************************************************}
Program Test01;
Uses WinCrt,Complex2;

Var  Z1,Z2,Z3: Complex;
     x,y,r,t:Double;

BEGIN
  x:=1; y:=-2;
  AssignXY(Z1,x,y);
  x:=2.5; y:=5.75;
  AssignXY(Z2,x,y);
  AddComplex(Z3,Z1,Z2);
  writeln;
  write(' Z1=');
  DisplayComplex(Z1,4,2);
  writeln;
  write(' Z2=');
  DisplayComplex(Z2,4,2);
  writeln;
  write(' Z1+Z2=');
  DisplayComplex(Z3,4,2);
  writeln;
  writeln(' Z1+Z2 with modulus and phase:'); 
  write(' '); DisplayComplexR(Z3,6,4);
  writeln;
  writeln(' Z1 with modulus and phase:'); 
  write(' '); DisplayComplexR(Z1,6,4);
  writeln;
  writeln(' Z2 with modulus and phase:');
  write(' '); DisplayComplexR(Z2,6,4);
  writeln;
  MulComplex(Z3,Z1,Z2);
  write(' Z1*Z2=');
  DisplayComplex(Z3,4,2);
  writeln;
  writeln(' Z1*Z2 with modulus and phase:'); 
  write(' '); DisplayComplexR(Z3,6,4);
  writeln;
  DivComplex(Z3,Z1,Z2);
  writeln(' Z1/Z2 with modulus and Phase:');
  write(' '); DisplayComplexR(Z3,6,4);
  writeln;
  ExpComplex(Z3,Z1);
  write(' Exp(Z1)=');
  DisplayComplex(Z3,4,2);
  writeln;
  writeln(' Exp(Z1) with modulus and Phase:');
  write(' '); DisplayComplexR(Z3,6,4);
  writeln;
  ArgShComplex(Z3,Z2);
  write(' ArgSh(Z2)=');
  DisplayComplex(Z3,6,4);
  writeln;
  Readkey; Donewincrt

END.

{end of file tcomplex.pas}