{**************************************************************************************************
*                 Inversion of a real square matrix by Householder's method                       *
* ----------------------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                                     *
*                                                                                                 *
* Input file householder.dat contains:                                                            *
*                                                                                                 *
* 4                                                                                               *
* 1.0       0.5        0.333333  0.25                                                             *
* 0.5       0.333333   0.25      0.2                                                              *
* 0.333333  0.25       0.2       0.166667                                                         *
* 0.25      0.2        0.166667  0.142857                                                         *
*                                                                                                 *
* Output file householder.txt contains:                                                           *
*                                                                                                 *
* Input square matrix:                                                                            *
*  1.00000000000000E+000   5.00000000000000E-001   3.33333000000000E-001   2.50000000000000E-001  *
*  5.00000000000000E-001   3.33333000000000E-001   2.50000000000000E-001   2.00000000000000E-001  *
*  3.33333000000000E-001   2.50000000000000E-001   2.00000000000000E-001   1.66667000000000E-001  *
*  2.50000000000000E-001   2.00000000000000E-001   1.66667000000000E-001   1.42857000000000E-001  *
*                                                                                                 *
* Transformation #1                                                                               *
* -1.19315166214903E+000  -6.70492926740837E-001  -4.74932372787676E-001  -3.69835516397145E-001  *
*  0.00000000000000E+000   6.64812024380633E-002   6.57297129201566E-002   5.86884055729692E-002  *
*  0.00000000000000E+000   7.20989795265739E-002   7.71532647936291E-002   7.24593645897091E-002  *
*  0.00000000000000E+000   6.65741012190317E-002   7.45318564600783E-002   7.22012027864846E-002  *
* Transformation #2                                                                               *
* -1.19315166214903E+000  -6.70492926740837E-001  -4.74932372787676E-001  -3.69835516397145E-001  *
*  0.00000000000000E+000  -1.18533219307947E-001  -1.25655520474182E-001  -1.17542173236451E-001  *
*  0.00000000000000E+000   0.00000000000000E+000   2.57161924423547E-003   3.78339450247119E-003  *
*  0.00000000000000E+000   0.00000000000000E+000   5.66533362590667E-003   8.78780895033711E-003  *
* Transformation #3                                                                               *
* -1.19315166214903E+000  -6.70492926740837E-001  -4.74932372787676E-001  -3.69835516397145E-001  *
*  0.00000000000000E+000  -1.18533219307947E-001  -1.25655520474182E-001  -1.17542173236451E-001  *
*  0.00000000000000E+000   0.00000000000000E+000  -6.22167426262024E-003  -9.56580449944889E-003  *
*  0.00000000000000E+000   0.00000000000000E+000   0.00000000000000E+000   1.87201461739764E-004  *
*                                                                                                 *
* Inverted matrix:                                                                                *
*  1.60326813541718E+001  -1.20367368862529E+002   2.40887472621963E+002  -1.40578963337616E+002  *
* -1.20367368862531E+002   1.20413864092798E+003  -2.71000936710905E+003   1.68653440306041E+003  *
*  2.40887472621969E+002  -2.71000936710907E+003   6.50422122896783E+003  -4.21581745593188E+003  *
* -1.40578963337620E+002   1.68653440306043E+003  -4.21581745593189E+003   2.81033136738221E+003  *
*                                                                                                 *
* Determinant: -1.64722238571522E-007                                                             *
*                                                                                                 *
* AP=A*A^-1 matrix:                                                                               *
*  9.99999999999965E-001   5.11590769747272E-013  -9.09494701772928E-013   6.82121026329696E-013  *
* -5.82346670885414E-015   1.00000000000008E+000  -1.15019105351166E-013   7.66608998503671E-014  *
* -3.08260361681079E-015   3.41116024316079E-014   9.99999999999921E-001   5.43454170554014E-014  *
* -2.43208231331948E-015   4.05508959744339E-014  -6.53921361504217E-014   1.00000000000006E+000  *
*                                                                                                 *
* ----------------------------------------------------------------------------------------------- *
* Ref.: "Méthodes de calcul numérique - Tome 2 - By Claude Nowakowski, PSI Editions, 1984".       *
*                                                                                                 *
*                                                 Turbo Pascal Version By J-P Moreau, Paris.      *
**************************************************************************************************}
Program Householder;
Uses WinCrt;
Const
    SIZE=25; SIZE1=50;
Var
    n: Integer;  {size of matrix A}
    A: Array[1..SIZE,1..SIZE1] of double;
    A0:Array[1..SIZE,1..SIZE] of double;
    V: Array[1..SIZE] of double;
    D, S, XA, XB, XG, XL: Double;
    B: Array[1..SIZE] of double;
    AI:Array[1..SIZE,1..SIZE] of double;
    AP:Array[1..SIZE,1..SIZE] of double;
    X: Array[1..SIZE] of double;

    fp1,fp2: TEXT;

    i,j,k: Integer;

    Function Sign(x:double): double;
    Begin
      if x<0.0 then
	    Sign:=-1.0
      else if x>0.0 then
	    Sign:=1.0
	  else
	    Sign:=0.0	
    End;

    Procedure P1000;
	Var
	    i, j, K: Integer; S: Double;
    Begin
        WriteLn(fp2);
        For K := 1 To n - 1 do
	begin
            S := 0;
            For i := K To n do
                S := S + A[i, K] * A[i, K];

	    XL := -Sign(A[K, K]) * Sqrt(S);
            XA := XL * (XL - A[K, K]);
            V[K] := A[K, K] - XL;
            For i := K + 1 To n do
                V[i] := A[i, K];
            For j := K + 1 To 2 * n do
	    begin
                XB := 0.0;
                For i := K To n do
                    XB := XB + V[i] * A[i, j];
                XG := XB / XA;
                For i := K To n do
                    A[i, j] := A[i, j] - XG * V[i];
            end;
            A[K, K] := XL;
            For i := K + 1 To n do
                A[i, K] := 0.0;
            WriteLn(fp2, ' Transformation #', K);
            For i := 1 To n do
                For j := 1 To n do
                    If j = n Then
                        WriteLn(fp2, '  ', A[i, j])
                    Else
                        Write(fp2, '  ', A[i, j])
        end
    End;

    Procedure P2000;
    Label 20, 40;
    Var
        i, j: Integer; S: Double;
    Begin
        For i := n DownTo 1 do
		begin
            j := n; S := 0.0;
20:         If i = j Then GoTo 40;
            S := S - A[i, j] * X[j];
            j := j - 1;
            GoTo 20;
40:         X[i] := (B[i] + S) / A[i, i]
        end
    End;

    Procedure P3000;
	Var
	    i, j: Integer;
		S: Double;
    Begin
        For i := 1 To n do
            For j := 1 To n do
			begin
                S := 0.0;
                For k := 1 To n do
                    S := S + A0[i, k] * AI[k, j];
                AP[i, j] := S
            end
    End;

    Procedure P4000;
	Begin
        D := 1.0;
        For i := 1 To n do
            D := D * A[i, i]
    End;


    {main program}
    BEGIN
        {read matrix to be inverted
         Open input/output text files}
        Assign(fp1,'househol.dat'); Reset(fp1);
        Assign(fp2,'househol.txt'); Rewrite(fp2);
        WriteLn(fp2);
        WriteLn(fp2, ' Input square matrix:');
        Read(fp1, n);
        For i := 1 To n do
            For j := 1 To n do
	    begin
                if j<n then
                  Read(fp1, A[i, j])
                else
                  Readln(fp1, A[i, j]);
                If j = n Then
                    WriteLn(fp2, '  ', A[i, j])
                Else
		    Write(fp2, '  ', A[i, j])
            end;
        Close(fp1);
        For i := 1 To n do
		begin
            For j := 1 To n do
			begin
                A[i, j + n] := 0.0;
                A0[i, j] := A[i, j]
            end;
            A[i, i + n] := 1.0
        end;
        {Transform A into triangular matrix}
        P1000;
        {N linear systems to solve}
        For K := 1 To n do
		begin
            For i := 1 To n do
                B[i] := A[i, K + n];
            {Solve triangular system}
            P2000;
            For i := 1 To n do
                AI[i, K] := X[i]
        end;
        WriteLn(fp2);
        WriteLn(fp2, ' Inverted matrix:');
        For i := 1 To n do
            For j := 1 To n do
                If j = n Then
                    WriteLn(fp2, '  ', AI[i, j])
                Else
                    Write(fp2, '  ', AI[i, j]);
        {Calculate determinant}
        P4000;
		WriteLn(fp2);
        WriteLn(fp2, ' Determinant: ', D);
        {Check AP:=A*A^-1:=Unity Matrix}
        WriteLn(fp2);
        WriteLn(fp2, ' AP=A*A^-1 matrix:');
        P3000;
        For i := 1 To n do
            For j := 1 To n do
                If j = n Then
                    WriteLn(fp2, '  ', AP[i, j])
                Else
				    Write(fp2, '  ', AP[i, j]);

        Close(fp2)
    END.

{end of file househol.pas}