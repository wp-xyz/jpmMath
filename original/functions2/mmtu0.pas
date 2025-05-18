{**********************************************************************
*       Calculate Mathieu Functions and their First Derivatives       *
* ------------------------------------------------------------------- *
* EXPLANATION:                                                        *
*                                                                     *
*      Purpose: This program computes Mathieu functions cem(x,q),     *
*               sem(x,q) and their derivatives using subroutine       *
*               MTU0 ( q = 0 )                                        *
*      Input :  KF  --- Function code                                 *
*                       KF=1 for computing cem(x,q) and cem'(x,q)     *
*                       KF=2 for computing sem(x,q) and sem'(x,q)     *
*               m   --- Order of Mathieu functions                    *
*               q   --- Parameter of Mathieu functions                *
*               x   --- Argument of Mathieu functions (in degrees)    *
*      Output:  CSF --- cem(x,q) or sem(x,q)                          *
*               CSD --- cem'x,q) or sem'x,q)                          *
*      Example: x = 40                                                *
*          m     q    cem(x,q)   cem'(x,q)    sem(x,q)  sem'(x,q)     *
*         --------------------------------------------------------    *
*          0    5.0   .3025683    .9470247                            *
*          1    5.0   .7669652   1.2873097    .2988052   .9606824     *
*          2    5.0   .9102723   -.3463855    .7549264  1.4743128     *
*          5    5.0  -.9810931   -.6328576    .1694850 -4.8676455     *
*          0   25.0   .0515371    .3823737                            *
*          1   25.0   .2074402   1.2646301    .0515365   .3823777     *
*          2   25.0  -.5297051  -2.4292679    .2074275  1.2646996     *
*          5   25.0   .7507159  -3.9047012   1.1881232   .3258081     *
*         --------------------------------------------------------    *
*                                                                     *
* ------------------------------------------------------------------- *
* SAMPLE RUNS:                                                        *
*                                                                     *
*  Please enter KF, m, q and x (in degrees): 1 5 25 40                *
*  KF = 1, m = 5, q = 25.0, x = 40.0                                  *
*                                                                     *
*   x(degs)    cem(x,q)       cem'(x,q)                               *
*   -------------------------------------                             *
*    40.0    0.750715922     -3.904701223                             *
*                                                                     *
*  Please enter KF, m, q and x (in degrees): 2 5 25 40                *
*  KF = 2, m = 5, q = 25.0, x = 40.0                                  *
*                                                                     *
*   x(degs)    sem(x,q)       sem'(x,q)                               *
*   -------------------------------------                             *
*    40.0    1.188123243      0.325808111                             *
*                                                                     *
* ------------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special Functions,  *
*             jin.ece.uiuc.edu/routines/routines.html".               *
*                                                                     *
*                               Pascal Release By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)            *
**********************************************************************}
      PROGRAM MMTU0;
      Uses WinCrt;

      Const NMAX = 251;

      Type  pVEC = ^VEC;
             VEC = Array[1..NMAX] of Double;

      Var   KF,M: Integer;
            CSF,CSD,Q,X : Double;

      {Headers of functions used below}
      Procedure CVA2(KD,M:Integer; Q:Double; Var A: Double); Forward;
      Procedure FCOEF(KD,M:Integer; Q,A:Double; Var FC:pVEC); Forward;
      Procedure CV0(KD,M:Integer; Q:Double; Var A0:Double); Forward;
      Procedure REFINE(KD,M:Integer; Q:Double; Var A:Double; IFLAG:Integer); Forward;
      Procedure CVQM(M:Integer; Q:Double; Var A0:Double); Forward;
      Procedure CVQL(KD,M:Integer; Q:Double; Var A0:Double); Forward;
      Procedure CVF(KD,M:Integer; Q,A:Double; MJ:Integer; Var F:Double); Forward;


      Procedure MTU0(KF,M:Integer; Q,X:Double; Var CSF,CSD:Double);
{       ===============================================================
        Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
                 and their derivatives ( q = 0 )
        Input :  KF  --- Function code
                         KF=1 for computing cem(x,q) and cem'(x,q)
                         KF=2 for computing sem(x,q) and sem'(x,q)
                 m   --- Order of Mathieu functions
                 q   --- Parameter of Mathieu functions
                 x   --- Argument of Mathieu functions (in degrees)
        Output:  CSF --- cem(x,q) or sem(x,q)
                 CSD --- cem'x,q) or sem'x,q)
        Routines called:
             (1) CVA2 for computing the characteristic values
             (2) FCOEF for computing the expansion coefficients
        ===============================================================}
      Label 15, 25;
      Var
            FG:pVEC;
            A,EPS,QM,RD,XR:Double;
            IC,K,KD,KM:Integer;
      Begin

        New(FG);
        EPS:=1E-14;

        IF (KF = 1) AND (M = 2*(M Div 2)) Then KD:=1;
        IF (KF = 1) AND (M <> 2*(M Div 2)) Then KD:=2;
        IF (KF = 2) AND (M <> 2*(M Div 2)) Then KD:=3;
        IF (KF = 2) AND (M = 2*(M Div 2)) Then KD:=4;

        CVA2(KD,M,Q,A);

        IF Q <= 1.0 THEN
           QM:=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM:=17.0+3.1*SQRT(Q)-0.126*Q+0.0037*SQRT(Q)*Q;

        KM:=Round(Int(QM+0.5*M));

        FCOEF(KD,M,Q,A,FG);

        IC:=Round(M Div 2)+1;
        RD:=1.74532925199433E-2;
        XR:=X*RD;
        CSF:=0.0;
        For K:=1 to KM do
        begin
           IF KD = 1 THEN
              CSF:=CSF+FG^[K]*COS((2*K-2)*XR)
           ELSE IF KD = 2 THEN
              CSF:=CSF+FG^[K]*COS((2*K-1)*XR)
           ELSE IF KD = 3 THEN
              CSF:=CSF+FG^[K]*SIN((2*K-1)*XR)
           ELSE IF KD = 4 THEN
              CSF:=CSF+FG^[K]*SIN(2*K*XR);
           IF (K >= IC) AND (ABS(FG^[K]) < ABS(CSF)*EPS) Then GOTO 15
        end;
15:     CSD:=0.0;
        For K:=1 to KM do
        begin
           IF KD = 1 THEN
              CSD:=CSD-(2*K-2)*FG^[K]*SIN((2*K-2)*XR)
           ELSE IF KD = 2 THEN
              CSD:=CSD-(2*K-1)*FG^[K]*SIN((2*K-1)*XR)
           ELSE IF KD = 3 THEN
              CSD:=CSD+(2*K-1)*FG^[K]*COS((2*K-1)*XR)
           ELSE IF KD = 4 THEN
              CSD:=CSD+2.0*K*FG^[K]*COS(2*K*XR);
           IF (K >= IC) AND (ABS(FG^[K]) < ABS(CSD)*EPS) Then GOTO 25
        end;
25:   Dispose(FG)
      End;

      Procedure FCOEF(KD,M:Integer; Q,A:Double; Var FC:pVEC);
{       =====================================================
        Purpose: Compute expansion coefficients for Mathieu
                 functions and modified Mathieu functions
        Input :  m  --- Order of Mathieu functions
                 q  --- Parameter of Mathieu functions
                 KD --- Case code
                        KD=1 for cem(x,q)  ( m = 0,2,4,...)
                        KD=2 for cem(x,q)  ( m = 1,3,5,...)
                        KD=3 for sem(x,q)  ( m = 1,3,5,...)
                        KD=4 for sem(x,q)  ( m = 2,4,6,...)
                 A  --- Characteristic value of Mathieu
                        functions for given m and q
        Output:  FC(k) --- Expansion coefficients of Mathieu
                        functions ( k= 1,2,...,KM )
                        FC(1),FC(2),FC(3),... correspond to
                        A0,A2,A4,... for KD=1 case, A1,A3,
                        A5,... for KD=2 case, B1,B3,B5,...
                        for KD=3 case and B2,B4,B6,... for
                        KD=4 case
        =====================================================}
      Label 45,70,85, Return;
      Var F,F1,F2,F3,QM,S,S0,SP,SS,TMP,U,V: Double;
          I,J,K,KB,KM:Integer;
      Begin

        IF Q <= 1.0 THEN
           QM:=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM:=17.0+3.1*SQRT(Q)-0.126*Q+0.0037*SQRT(Q)*Q;
        KM:=Round(Int(QM+0.5*M));
        IF Q = 0.0 THEN
        begin
           For K:=1 to KM do FC^[K]:=0.0;
           IF KD = 1 THEN
           begin
              FC^[(M+2) Div 2] := 1.0;
              IF M = 0 Then FC^[1]:=1.0/SQRT(2.0)
           end
           ELSE IF KD = 4 THEN
              FC^[M Div 2]:=1.0
           ELSE
              FC^[(M+1) Div 2]:=1.0;
           goto Return
        end;
        KB:=0;
        S:=0.0;
        F:=1.0E-100;
        U:=0.0;
        FC^[KM]:=0.0;
        IF KD = 1 THEN
        begin
           For K:=KM Downto 3 do
           begin
              V:=U;
              U:=F;
              F:=(A-4.0*K*K)*U/Q-V;
              IF ABS(F) < ABS(FC^[K+1]) THEN
              begin
                 KB:=K;
                 FC^[1]:=1.0E-100;
                 SP:=0.0;
                 F3:=FC^[K+1];
                 FC^[2]:=A/Q*FC^[1];
                 FC^[3]:=(A-4.0)*FC^[2]/Q-2.0*FC^[1];
                 U:=FC^[2];
                 F1:=FC^[3];
                 For I:=3 to KB do
                 begin
                   V:=U;
                   U:=F1;
                   F1:=(A-4.0*Sqr(I-1.0))*U/Q-V;
                   FC^[I+1]:=F1;
                   IF I = KB Then F2:=F1;
                   IF I <> KB Then SP:=SP+F1*F1
                 end;
                 SP:=SP+2.0*Sqr(FC^[1])+Sqr(FC^[2])+Sqr(FC^[3]);
                 SS:=S+SP*Sqr(F3/F2);
                 S0:=SQRT(1.0/SS);
                 For J:=1 to KM do
                   IF J <= KB+1 THEN
                     FC^[J]:=S0*FC^[J]*F3/F2
                   ELSE
                     FC^[J]:=S0*FC^[J];
                 GOTO 85
              end
              ELSE
              begin
                FC^[K]:=F;
                S:=S+F*F
              end
           end; {of K loop}

           FC^[2]:=Q*FC^[3]/(A-4.0-2.0*Q*Q/A);
           FC^[1]:=Q/A*FC^[2];
           S:=S+2.0*Sqr(FC^[1])+Sqr(FC^[2]);
           S0:=SQRT(1.0/S);
           For K:=1 to KM do FC^[K]:=S0*FC^[K]
        end
        ELSE IF (KD = 2) OR (KD = 3) THEN
        begin
           For K:=KM Downto 3 do
           begin
              V:=U;
              U:=F;
              F:=(A-Sqr(2.0*K-1.0))*U/Q-V;
              IF ABS(F) >= ABS(FC^[K]) THEN
              begin
                 FC^[K-1]:=F;
                 S:=S+F*F
              end
              ELSE
              begin
                 KB:=K;
                 F3:=FC^[K];
                 GOTO 45
              end
           end;
           {TMP=(-1)^KD}
           if KD=2*(KD Div 2) Then TMP:=1.0
                              else TMP:=-1.0;
           FC^[1]:=Q/(A-1.0-TMP*Q)*FC^[2];
           S:=S+FC^[1]*FC^[1];
           S0:=SQRT(1.0/S);
           For K:=1 to KM do FC^[K]:=S0*FC^[K];
           GOTO 85;
45:        FC^[1]:=1.0E-100;
           {TMP=(-1)^KD}
           if KD=2*(KD Div 2) Then TMP:=1.0
                              else TMP:=-1.0;
           FC^[2]:=(A-1.0-TMP*Q)/Q*FC^[1];
           SP:=0.0;
           U:=FC^[1];
           F1:=FC^[2];
           For I:=2 to KB-1 do
           begin
              V:=U;
              U:=F1;
              F1:=(A-Sqr(2.0*I-1.0))*U/Q-V;
              IF I <> KB-1 THEN
              begin
                FC^[I+1]:=F1;
                SP:=SP+F1*F1
              end
              ELSE
                F2:=F1
           end;
           SP:=SP+Sqr(FC^[1]) + Sqr(FC^[2]);
           SS:=S+SP*Sqr(F3/F2);
           S0:=1.0/SQRT(SS);
           For J:=1 to KM do
           begin
             IF J < KB Then FC^[J]:=S0*FC^[J]*F3/F2;
             IF J >= KB Then FC^[J]:=S0*FC^[J]
           end
        end
        ELSE IF KD = 4 THEN
        begin
           For K:=KM Downto 3 do
           begin
              V:=U;
              U:=F;
              F:=(A-4.0*K*K)*U/Q-V;
              IF ABS(F) >= ABS(FC^[K]) THEN
              begin
                FC^[K-1]:=F;
                S:=S+F*F
              end
              ELSE
              begin
                KB:=K;
                F3:=FC^[K];
                GOTO 70
              end
           end;
           FC^[1]:=Q/(A-4.0)*FC^[2];
           S:=S+FC^[1]*FC^[1];
           S0:=SQRT(1.0/S);
           For K:=1 to KM do FC^[K]:=S0*FC^[K];
           GOTO 85;
70:        FC^[1]:=1.0E-100;
           FC^[2]:=(A-4.0)/Q*FC^[1];
           SP:=0.0;
           U:=FC^[1];
           F1:=FC^[2];
           For I:=2 to KB-1 do
           begin 
             V:=U;
             U:=F1;
             F1:=(A-4.0*I*I)*U/Q-V;
             IF I <> KB-1 THEN
             begin
                FC^[I+1]:=F1;
                SP:=SP+F1*F1
             end
             ELSE
                F2:=F1
           end;
           SP:=SP+Sqr(FC^[1]) + Sqr(FC^[2]);
           SS:=S+SP*Sqr(F3/F2);
           S0:=1.0/SQRT(SS);
           For J:=1 to KM do
           begin
             IF J < KB Then FC^[J]:=S0*FC^[J]*F3/F2;
             IF J >= KB Then FC^[J]:=S0*FC^[J]
           end
        end;
85:     IF FC^[1] < 0.0 THEN
           For J:=1 to KM do FC^[J]:=-FC^[J];
Return: End;

      Procedure CVA2(KD,M:Integer; Q:Double; Var A: Double);
{       ======================================================
        Purpose: Calculate a specific characteristic value of
                 Mathieu functions
        Input :  m  --- Order of Mathieu functions
                 q  --- Parameter of Mathieu functions
                 KD --- Case code
                        KD=1 for cem(x,q)  ( m = 0,2,4,...)
                        KD=2 for cem(x,q)  ( m = 1,3,5,...)
                        KD=3 for sem(x,q)  ( m = 1,3,5,...)
                        KD=4 for sem(x,q)  ( m = 2,4,6,...)
        Output:  A  --- Characteristic value
        Routines called:
              (1) REFINE for finding accurate characteristic
                  values using an iteration method
              (2) CV0 for finding initial characteristic
                  values using polynomial approximation
              (3) CVQM for computing initial characteristic
                  values for q ó 3*m
              (3) CVQL for computing initial characteristic
                  values for q ò m*m
        ======================================================}
      Label 5,15;
      Var I,IFLAG,NDIV,NN: Integer;
          A1,A2,DELTA,Q1,Q2,QQ: Double;
      Begin  
        IF (M <= 12) OR (Q <= 3.0*M) OR (Q > M*M) THEN
        begin
          CV0(KD,M,Q,A);
          IF Q <> 0.0 Then REFINE(KD,M,Q,A,1)
        end
        ELSE
        begin
           NDIV:=10;
           DELTA:=(M-3.0)*M/NDIV;
           IF Q-3.0*M <= M*M-Q THEN
           begin
5:            NN:=Round((Q-3.0*M)/DELTA)+1;
              DELTA:=(Q-3.0*M)/NN;
              Q1:=2.0*M;
              CVQM(M,Q1,A1);
              Q2:=3.0*M;
              CVQM(M,Q2,A2);
              QQ:=3.0*M;
              For I:=1 to NN do
              begin
                 QQ:=QQ+DELTA;
                 A:=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1);
                 IFLAG:=1;
                 IF I = NN Then IFLAG:=-1;
                 REFINE(KD,M,QQ,A,IFLAG);
                 Q1:=Q2;
                 Q2:=QQ;
                 A1:=A2;
                 A2:=A
              end;
              IF IFLAG = -10 THEN
              begin
                 NDIV:=NDIV*2;
                 DELTA:=(M-3.0)*M/NDIV;
                 GOTO 5
              end
           end
           ELSE
           begin
15:           NN:=Round((M*M-Q)/DELTA)+1;
              DELTA:=(M*M-Q)/NN;
              Q1:=M*(M-1.0);
              CVQL(KD,M,Q1,A1);
              Q2:=M*M;
              CVQL(KD,M,Q2,A2);
              QQ:=M*M;
              For I:=1 to NN do
              begin
                QQ:=QQ-DELTA;
                A:=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1);
                IFLAG:=1;
                IF I = NN Then IFLAG:=-1;
                REFINE(KD,M,QQ,A,IFLAG);
                Q1:=Q2;
                Q2:=QQ;
                A1:=A2;
                A2:=A
              end;
              IF IFLAG = -10 THEN
              begin
                NDIV:=NDIV*2;
                DELTA:=(M-3.0)*M/NDIV;
                GOTO 15
              end
           end
        end
      End;

      Procedure REFINE(KD,M:Integer; Q:Double; Var A:Double; IFLAG:Integer);
{       =====================================================
        Purpose: calculate the accurate characteristic value
                 by the secant method
        Input :  m --- Order of Mathieu functions
                 q --- Parameter of Mathieu functions
                 A --- Initial characteristic value
        Output:  A --- Refineed characteristic value
        Routine called:  CVF for computing the value of F for
                         characteristic equation
       ========================================================}
      Label 5, 15, Return;
      Var CA,DELTA,F,F0,F1,EPS,X,X0,X1: Double;
          IT,MJ:Integer;
      Begin
        EPS:=1E-14;
        MJ:=10+M;
        CA:=A;
        DELTA:=0.0;
        X0:=A;
        CVF(KD,M,Q,X0,MJ,F0);
        X1:=1.002*A;
        CVF(KD,M,Q,X1,MJ,F1);
5:      For IT:=1 to 100 do
        begin
          MJ:=MJ+1;
          X:=X1-(X1-X0)/(1.0-F0/F1);
          CVF(KD,M,Q,X,MJ,F);
          IF (ABS(1.0-X1/X) < EPS) OR (F = 0.0) Then GOTO 15;
          X0:=X1;
          F0:=F1;
          X1:=X;
          F1:=F
        end;
15:     A:=X;
        IF DELTA > 0.05 THEN
        begin
           A:=CA;
           IF IFLAG < 0 THEN IFLAG:=-10;
           goto RETURN
        end;
        IF ABS((A-CA)/CA) > 0.05 THEN
        begin
          X0:=CA;
          DELTA:=DELTA+0.005;
          CVF(KD,M,Q,X0,MJ,F0);
          X1:=(1.0+DELTA)*CA;
          CVF(KD,M,Q,X1,MJ,F1);
          GOTO 5
        end;
Return: End;

      Procedure CVF(KD,M:Integer; Q,A:Double; MJ:Integer; Var F:Double);
{       ======================================================
        Purpose: Compute the value of F for characteristic
                 equation of Mathieu functions
        Input :  m --- Order of Mathieu functions
                 q --- Parameter of Mathieu functions
                 A --- Characteristic value
        Output:  F --- Value of F for characteristic equation
        ======================================================}
      Var B,T0,T1,T2: Double;
          IC,J,J0,JF,L,L0: Integer;
      Begin
        B:=A;
        IC:=M Div 2;
        L:=0;
        L0:=0;
        J0:=2;
        JF:=IC;
        IF KD = 1 Then L0:=2;
        IF KD = 1 Then J0:=3;
        IF (KD = 2) OR (KD = 3) Then L:=1;
        IF KD = 4 Then JF:=IC-1;
        T1:=0.0;
        For J:=MJ Downto IC+1 do
          T1:=-Q*Q/(Sqr(2.0*J+L)-B+T1);
        IF M <= 2 THEN
        begin
          T2:=0.0;
          IF (KD = 1) AND (M = 0) Then T1:=T1+T1;
          IF (KD = 1) AND (M = 2) Then T1:=-2.0*Q*Q/(4.0-B+T1)-4.0;
          IF (KD = 2) AND (M = 1) Then T1:=T1+Q;
          IF (KD = 3) AND (M = 1) Then T1:=T1-Q
        end
        ELSE
        begin
          IF KD = 1 Then T0:=4.0-B+2.0*Q*Q/B;
          IF KD = 2 Then T0:=1.0-B+Q;
          IF KD = 3 Then T0:=1.0-B-Q;
          IF KD = 4 Then T0:=4.0-B;
          T2:=-Q*Q/T0;
          For J:=J0 to JF do T2:=-Q*Q/(Sqr(2.0*J-L-L0)-B+T2)
        end;
        F:=Sqr(2.0*IC+L) + T1 + T2 - B
      End;

      Procedure CV0(KD,M:Integer; Q:Double; Var A0:Double);
{       =====================================================
        Purpose: Compute the initial characteristic value of
                 Mathieu functions for m = 12  or q = 300 or
                 q = m*m
        Input :  m  --- Order of Mathieu functions
                 q  --- Parameter of Mathieu functions
        Output:  A0 --- Characteristic value
        Routines called:
              (1) CVQM for computing initial characteristic
                  value for q = 3*m
              (2) CVQL for computing initial characteristic
                  value for q = m*m
        ====================================================}
      Var Q2: Double;
      Begin
        Q2:=Q*Q;
        IF M = 0 THEN
        begin
          IF Q <= 1.0 THEN
            A0:=(((0.0036392*Q2-0.0125868)*Q2+0.0546875)*Q2-0.5)*Q2
          ELSE IF Q <= 10.0 THEN
            A0:=((3.999267E-3*Q-9.638957E-2)*Q-0.88297)*Q +0.5542818
          ELSE
            CVQL(KD,M,Q,A0)
        end
        ELSE IF M = 1 THEN
        begin
          IF (Q <= 1.0) AND (KD = 2) THEN
             A0:=(((-6.51E-4*Q-0.015625)*Q-0.125)*Q+1.0)*Q+1.0
          ELSE IF (Q <= 1.0) AND (KD = 3) THEN
             A0:=(((-6.51E-4*Q+0.015625)*Q-0.125)*Q-1.0)*Q+1.0
          ELSE IF (Q <= 10.0) AND (KD = 2) THEN
             A0:=(((-4.94603E-4*Q+1.92917E-2)*Q-0.3089229)*Q+1.33372)*Q+0.811752
          ELSE IF (Q <= 10.0) AND (KD = 3) THEN
             A0:=((1.971096E-3*Q-5.482465E-2)*Q-1.152218)*Q+1.10427
          ELSE
             CVQL(KD,M,Q,A0)
        end
        ELSE IF M = 2 THEN
        begin
          IF (Q <= 1.0) AND (KD = 1) THEN
             A0:=(((-0.0036391*Q2+0.0125888)*Q2-0.0551939)*Q2+0.416667)*Q2+4.0
          ELSE IF (Q <= 1.0) AND (KD = 4) THEN
             A0:=(0.0003617*Q2-0.0833333)*Q2+4.0
          ELSE IF (Q <= 15.0) AND (KD = 1) THEN
             A0:=(((3.200972E-4*Q-8.667445E-3)*Q-1.829032E-4)*Q+0.9919999)*Q+3.3290504
          ELSE IF (Q <= 10.0) AND (KD = 4) THEN
             A0:=((2.38446E-3*Q-0.08725329)*Q-4.732542E-3)*Q+4.00909
          ELSE
             CVQL(KD,M,Q,A0)
        end
        ELSE IF M = 3 THEN
        begin
           IF (Q <= 1.0) AND (KD = 2) THEN
              A0:=((6.348E-4*Q+0.015625)*Q+0.0625)*Q2+9.0
           ELSE IF (Q <= 1.0) AND (KD = 3) THEN
              A0:=((6.348E-4*Q-0.015625)*Q+0.0625)*Q2+9.0
           ELSE IF (Q <= 20.0) AND (KD = 2) THEN
              A0:=(((3.035731E-4*Q-1.453021E-2)*Q+0.19069602)*Q-0.1039356)*Q+8.9449274
           ELSE IF (Q <= 15.0) AND (KD = 3) THEN
              A0:=((9.369364E-5*Q-0.03569325)*Q+0.2689874)*Q+8.771735
           ELSE
              CVQL(KD,M,Q,A0)
        end
        ELSE IF M = 4 THEN
        begin
           IF (Q <= 1.0) AND (KD = 1) THEN
              A0:=((-2.1E-6*Q2+5.012E-4)*Q2+0.0333333)*Q2+16.0
           ELSE IF (Q <= 1.0) AND (KD = 4) THEN
              A0:=((3.7E-6*Q2-3.669E-4)*Q2+0.0333333)*Q2+16.0
           ELSE IF (Q <= 25.0) AND (KD = 1) THEN
              A0:=(((1.076676E-4*Q-7.9684875E-3)*Q+0.17344854)*Q-0.5924058)*Q+16.620847
           ELSE IF (Q <= 20.0) AND (KD = 4) THEN
              A0:=((-7.08719E-4*Q+3.8216144E-3)*Q+0.1907493)*Q+15.744
           ELSE
              CVQL(KD,M,Q,A0)
        end
        ELSE IF M = 5 THEN
        begin
           IF (Q <= 1.0) AND (KD = 2) THEN
              A0:=((6.8E-6*Q+1.42E-5)*Q2+0.0208333)*Q2+25.0
           ELSE IF (Q <= 1.0) AND (KD = 3) THEN
              A0:=((-6.8E-6*Q+1.42E-5)*Q2+0.0208333)*Q2+25.0
           ELSE IF (Q <= 35.0) AND (KD = 2) THEN
              A0:=(((2.238231E-5*Q-2.983416E-3)*Q+0.10706975)*Q-0.600205)*Q+25.93515
           ELSE IF (Q <= 25.0) AND (KD = 3) THEN
              A0:=((-7.425364E-4*Q+2.18225E-2)*Q+4.16399E-2)*Q+24.897
           ELSE
              CVQL(KD,M,Q,A0)
        end
        ELSE IF (M = 6) THEN
        begin
           IF (Q <= 1.0) THEN
              A0:=(0.4E-6*Q2+0.0142857)*Q2+36.0
           ELSE IF (Q <= 40.0) AND (KD = 1) THEN
              A0:=(((-1.66846E-5*Q+4.80263E-4)*Q+2.53998E-2)*Q-0.181233)*Q+36.423
           ELSE IF (Q <= 35.0) AND (KD = 4) THEN
              A0:=((-4.57146E-4*Q+2.16609E-2)*Q-2.349616E-2)*Q+35.99251
           ELSE
              CVQL(KD,M,Q,A0)
        end
        ELSE IF (M = 7) THEN
        begin
           IF (Q <= 10.0) THEN
              CVQM(M,Q,A0)
           ELSE IF (Q <= 50.0) AND (KD = 2) THEN
              A0:=(((-1.411114E-5*Q+9.730514E-4)*Q-3.097887E-3)*Q+3.533597E-2)*Q+49.0547
           ELSE IF (Q <= 40.0) AND (KD = 3) THEN
              A0:=((-3.043872E-4*Q+2.05511E-2)*Q-9.16292E-2)*Q+49.19035
           ELSE
              CVQL(KD,M,Q,A0)
        end
        ELSE IF M >= 8 THEN
        begin
           IF Q <= 3.*M THEN
              CVQM(M,Q,A0)
           ELSE IF Q > M*M THEN
              CVQL(KD,M,Q,A0)
           ELSE
           begin
              IF (M = 8) AND (KD = 1) THEN
                 A0:=(((8.634308E-6*Q-2.100289E-3)*Q+0.169072)*Q-4.64336)*Q+109.4211
              ELSE IF (M = 8) AND (KD = 4) THEN
                 A0:=((-6.7842E-5*Q+2.2057E-3)*Q+0.48296)*Q+56.59
              ELSE IF (M = 9) AND (KD = 2) THEN
                 A0:=(((2.906435E-6*Q-1.019893E-3)*Q+0.1101965)*Q-3.821851)*Q+127.6098
              ELSE IF (M = 9) AND (KD = 3) THEN
                 A0:=((-9.577289E-5*Q+0.01043839)*Q+0.06588934)*Q+78.0198
              ELSE IF (M = 10) AND (KD = 1) THEN
                 A0:=(((5.44927E-7*Q-3.926119E-4)*Q+0.0612099)*Q-2.600805)*Q+138.1923
              ELSE IF (M = 10) AND (KD = 4) THEN
                 A0:=((-7.660143E-5*Q+0.01132506)*Q-0.09746023)*Q+99.29494
              ELSE IF (M = 11) AND (KD = 2) THEN
                 A0:=(((-5.67615E-7*Q+7.152722E-6)*Q+0.01920291)*Q-1.081583)*Q+140.88
              ELSE IF (M = 11) AND (KD = 3) THEN
                 A0:=((-6.310551E-5*Q+0.0119247)*Q-0.2681195)*Q+123.667
              ELSE IF (M = 12) AND (KD = 1) THEN
                 A0:=(((-2.38351E-7*Q-2.90139E-5)*Q+0.02023088)*Q-1.289)*Q+171.2723
              ELSE IF (M = 12) AND (KD = 4) THEN
                 A0:=(((3.08902E-7*Q-1.577869E-4)*Q+0.0247911)*Q-1.05454)*Q+161.471
           end
        end
      End; {CV0}

      Procedure CVQL(KD,M:Integer; Q:Double; Var A0:Double);
{       ========================================================
        Purpose: Compute the characteristic value of Mathieu
                 functions  for q ò 3m
        Input :  m  --- Order of Mathieu functions
                 q  --- Parameter of Mathieu functions
        Output:  A0 --- Initial characteristic value
        ========================================================}
      Var  C1,CV1,CV2,D1,D2,D3,D4,P1,P2,W,W2,W3,W4,W6: Double;
      Begin
        IF (KD = 1) OR (KD = 2) Then W:=2.0*M+1.0;
        IF (KD = 3) OR (KD = 4) Then W:=2.0*M-1.0;
        W2:=W*W;
        W3:=W*W2;
        W4:=W2*W2;
        W6:=W2*W4;
        D1:=5.0+34.0/W2+9.0/W4;
        D2:=(33.0+410.0/W2+405.0/W4)/W;
        D3:=(63.0+1260.0/W2+2943.0/W4+486.0/W6)/W2;
        D4:=(527.0+15617.0/W2+69001.0/W4+41607.0/W6)/W3;
        C1:=128.0;
        P2:=Q/W4;
        P1:=SQRT(P2);
        CV1:=-2.0*Q+2.0*W*SQRT(Q)-(W2+1.0)/8.0;
        CV2:=(W+3.0/W)+D1/(32.0*P1)+D2/(8.0*C1*P2);
        CV2:=CV2+D3/(64.0*C1*P1*P2)+D4/(16.0*C1*C1*P2*P2);
        A0:=CV1-CV2/(C1*P1)
      End;

      Procedure CVQM(M:Integer; Q:Double; Var A0:Double);
{       =====================================================
        Purpose: Compute the characteristic value of Mathieu
                 functions for q ó m*m
        Input :  m  --- Order of Mathieu functions
                 q  --- Parameter of Mathieu functions
        Output:  A0 --- Initial characteristic value
        =====================================================}
      Var  HM1,HM3,HM5: Double;
      Begin
        HM1:=0.5*Q/(M*M-1.0);
        HM3:=0.25*HM1*HM1*HM1/(M*M-4.0);
        HM5:=HM1*HM3*Q/((M*M-1.0)*(M*M-9.0));
        A0:=M*M+Q*(HM1+(5.0*M*M+7.0)*HM3+(9.0*M*M*M*M+58.0*M*M+29.0)*HM5)
      End;


  {main program}
  BEGIN

    WRITELN;;
    WRITE(' Please enter KF, m, q and x (in degrees): ');
    READLN(KF,M,Q,X);
    WRITELN(' KF = ',KF,', M = ',M,', Q = ',Q:4:1,', X = ',X:4:1);
    WRITELN;
    IF KF = 1 Then WRITELN(' x(degs.)    cem(x,q)       cem''(x,q)');
    IF KF = 2 Then WRITELN(' x(degs.)    sem(x,q)       sem''(x,q)');
    WRITELN(' --------------------------------------');

    MTU0(KF,M,Q,X,CSF,CSD);

    WRITELN('  ',X:5:1,CSF:16:9,CSD:16:9);
    WRITELN;

    ReadKey;
    DoneWinCrt

  END.

{ end of file mmtu0.pas}