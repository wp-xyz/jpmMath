{*****************************************************************
*      MULTI-DIMENSIONAL CURVE FITTING BY SIMPLEX METHOD         *
* -------------------------------------------------------------- *
* Reference:                                                     *
*                                                                *
*  J. P. CHANDLER, COMPUTER SCIENCE DEPARTMENT,                  *
*     OKLAHOMA STATE UNIVERSITY                                  ù
* -------------------------------------------------------------- *
* Explanations:                                                  *
*                                                                *
*  THIS PROGRAM FITS THE MODEL:                                  *
*                                                                *
*               A                                                *
*     Z = ----------------                                       *
*           1 + BX + CY                                          *
*                                                                *
*  WRITTEN AS:                                                   *
*                                                                *
*     FIT(J) = X(1)/(1.0 + X(2)*T(J,1) + X(3)*T(J,2))            *
*                                                                *
*  TO DATA POINTS  (T(J,1),T(J,2),Y(J)).                         *
*                                                                *
*  T(J,1) AND T(J,2) ARE THE TWO INDEPENDENT VARIABLES X AND Y   *
*  IN DATA SPACE.  Y(J) IS THE MEASURED DEPENDENT VARIABLE Z AND *
*  FIT(J) IS THE VALUE TO BE FITTED TO Y(J) BY MINIMIZING THE    *
*  WEIGHTED SUM OF SQUARES,                                      *
*                                                                *
*           J=NPTS                                               *
*                                                                *
*     PHI =  SUM ((FIT(J)-Y(J))/YSIG(J))**2                      *
*                                                                *
*            J=1                                                 *
*                                                                *
*  THE VALUE OF PHI IS STORED IN THE VARIABLE FOBJ.              *
*  THE MODEL Z=F(X,Y) IS IMPLEMENTED IN FUNCTION FUNK.           *
*  THE RESULTS ARE BEST COEFFICIENTS A, B, C, STORED IN X(1),    *
*  X(2), X(3) TO FIT THE DATA PROVIDED IN FILE SIMPTEST.DAT.     *
* -------------------------------------------------------------- *
* SAMPLE RUN:                                                    *
* File Simptest.dat contains:                                    *
* 12                                                             *
* 0.0       1.0       16.6      0.2                              *
* 1.0       0.0       12.4      0.2                              *
* 1.0       1.0       8.2       0.2                              *
* 0.0       2.0       10.1      0.2                              *
* 2.0       0.0       7.3       0.2                              *
* 2.0       2.0       4.7       0.2                              *
* 1.0       3.0       5.1       0.2                              *
* 3.0       1.0       4.3       0.2                              *
* 3.0       3.0       3.0       0.2                              *
* 5.0       1.0       2.5       0.2                              *
* 1.0       5.0       3.4       0.2                              *
* 5.0       5.0       2.1       0.2                              *
*                                                                *
* Output file Simpout.lst contains at the end:                   *
* --/--                                                          ***********
*                                                                          *
* TERMINATED WHEN THE DIMENSIONS OF THE SIMPLEX BECAME AS SMALL            *
* AS THE DELMIN(J).                                                        *
*                                                                          *
*                                                                          *
*                                                                          *
* 544 FUNCTION COMPUTATIONS:                                               *
*                                                                          *
* FINAL VALUE OF FOBJ =       6.5040337830664168                           *
*                                                                          *
* FINAL VALUES OF X(J):       48.2028883       2.8763282       1.9030361   *
*                                                                          *
*                                                                          *
*                       +   3.6498                                         *
*  X(1) =     48.202888                STDEVS =  1.00000                   *
                        -   3.0447                                         *
*                                                                          *
*                                                                          *
*                       +   0.2638                                         *
*  X(2) =      2.876328                STDEVS =  1.00000                   *
*                       -   0.2239                                         *
*                                                                          *
*                                                                          *
*                       +   0.1964                                         *
*  X(3) =      1.903036                STDEVS =  1.00000                   *
*                       -   0.1681                                         *
*                                                                          *
* ------------------------------------------------------------------------ *
*                                      TPW Release By J-P Moreau, Paris.   *
*                                              (www.jpmoreau.fr)           *
****************************************************************************
Note: To change example, adjust function FUNK to your problem and check
      NV, NPTS and size of T[300,2].                                       }
PROGRAM SIMPLEX;  
      
Uses WinCrt;

Type
      pTab20 = ^Tab20;
      Tab20 = Array[1..20] of Double;
      pTab300 = ^Tab300;
      Tab300 = Array[1..300] of Double;
      pMat = ^Mat;
      Mat = Array[1..20,1..21] of Double;
      pMat300 = ^Mat300;
      Mat300 = Array[1..300,1..2] of Double;
      Tab2 = Array[1..2] of Double;

Var   {global variables}
      ERGSAV:pTab20; FIDINT: Tab2;                                            
      
      X, XMAX, XMIN, DELTX, DELMIN: pTab20;
      MASK: Array[1..20] of Integer;
      ERFRAC, ERGUES,FOBJ,STDEVS,TOLFID: Double;
      JXFID,MAXIND,MATRX,NV,NTRACE,NTRACF,NTRSET: Integer;
                         
      K,NFMAX,NFLAT,JVARY,NXTRA,KFLAG,NOREP,KERFL: Integer;                           

      FLAMBD,FNU,RELDIF,RELMIN:Double;
      METHD,KALCP,KORDIF,MAXIT,LEQU,MAXSUB,MAXUPD,NPTS:Integer;

      {common variables for SMPLX, SIBEG...}
      FZ: pTab20;
      HUGE,ALPHA,BETA,GAMMA: Double;
      MAXPT,NVPLUS,NSSW,NF: Integer;                                                         

      Y, YSIG: pTab300;                              
      T: pMat300;
      Z: pMat;

      fp_in,fp_out: TEXT;                                                  

      Procedure SMPLX; Forward;

Procedure  FUNK;
{-------------------------------------------------------------
!  J. P. CHANDLER, COMPUTER SCIENCE DEPT.,                                      
!     OKLAHOMA STATE UNIVERSITY                                                 
!                                                                               
!  THIS VERSION OF SUBROUTINE FUNK COMPUTES THE FITTED VALUES                   
!  AND THE WEIGHTED SUM OF SQUARES FOR THE PROBLEM OF FITTING                   
!  THE MODEL                                                                    
!                                                                               
!     FIT(J) = X(1)/(1.0 + X(2)*T(J,1) + X(3)*T(J,2))                           
!                                                                               
!  TO DATA POINTS  (T(J,1),T(J,2),Y(J)).                                      
!                                                                               
!-------------------------------------------------------------}
Var
   FITM: pTab300;  J:Integer;     
Begin                                                   

{  FOBJ WILL CONTAIN THE WEIGHTED SUM OF SQUARES, PHI.                          
   SIMPLEX REQUIRES THAT FUNK COMPUTE FOBJ.                        
   Here FOBJ is a global variable. }

   New(FITM);   
   FOBJ:=0.0;                                                                
   For J:=1 to NPTS do                                                            
   begin        
     FITM^[J]:=X^[1]/(1.0+X^[2]*T^[J,1]+X^[3]*T^[J,2]);                           
     FOBJ:=FOBJ+Sqr((FITM^[J]-Y^[J])/YSIG^[J])                                  
   end;
   Dispose(FITM)                                                               
End;                                                                       


Procedure FIDO(JXFID:Integer;STDEVS,TOLFID,ERGUES:Double; MAXIND,
                 NTRSET,NTRACF:Integer; Var FIDINT:Tab2);                  
{-----------------------------------------------------------------
!  FIDO 4.5          DECEMBER 1991                                              
!  A.N.S.I. STANDARD FORTRAN 77                                                 
!                                                                               
!  J. P. CHANDLER, COMPUTER SCIENCE DEPARTMENT,                                 
!     OKLAHOMA UNIVERSITY, STILLWATER, OKLAHOMA 74078                           
!                                                                               
!  SUBROUTINE FIDO IS USED IN CONJUNCTION WITH SIMPLEX
!  TO COMPUTE CONFIDENCE HALF-INTERVALS ("ERRORS") FOR                  
!  THE PARAMETERS IN A FITTING PROBLEM, USING THE METHOD OF                     
!  SUPPORT PLANES.                                                              
!  THE QUANTITY BEING MINIMIZED (FOBJ) MUST EITHER BE                           
!  CHI-SQUARE OR IT MUST BE TWICE THE NEGATIVE OF THE NATURAL                   
!  LOGARITHM OF THE LIKELIHOOD FUNCTION.                                        
!                                                                               
! ---------------------------------------------------------------                 
!                                                                               
!  INPUT QUANTITIES.....   JXFID,STDEVS,TOLFID,ERGUES,                     
!                          MAXIND,NTRSET,NTRACF,X(*),                 
!                                                                               
!  OUTPUT QUANTITY......   FIDINT( )                                            
!                                                                               
!                                                                               
!     JXFID     --  THE INDEX OF THE PARAMETER, X(JXFID), FOR                   
!                   WHICH CONFIDENCE HALF-INTERVALS ARE TO BE                   
!                   COMPUTED                                              
!                                                                               
!     STDEVS    --  THE NUMBER OF STANDARD DEVIATIONS TO WHICH                  
!                   THE HALF-INTERVALS ARE TO CORRESPOND                     
!                                                                               
!     TOLFID    --  CONVERGENCE TOLERANCE (USE TOLFID=0.05)                     
!                                                                               
!     ERGUES    --  ROUGH ESTIMATE OF THE LENGTH OF THE                         
!                   HALF-INTERVALS                                           
!                                                                               
!     MAXIND    --  =1 TO COMPUTE THE HALF-INTERVAL WITH THE                    
!                      SAME SIGN AS ERGUES,                                   
!                   =2 TO COMPUTE BOTH HALF-INTERVALS                           
!                                                                               
!     NTRSET    --  VALUE OF NTRACE TO BE USED IN THE                           
!                   FITTING SUBROUTINE (STEPIT OR MARQ,                      
!                   ETC.) WHEN CALLED BY FIDO (Not used here).                                
!                                                                               
!     NTRACF    --  A SWITCH THAT CONTROLS PRINTING...                          
!                                                                               
!                      =-1 TO PRINT NOTHING IN FIDO EXCEPT ANY                  
!                          WARNING OR ERROR MESSAGES,                            
!                                                                               
!                      = 0 FOR INITIAL AND SUMMARY PRINTOUT,                    
!                                                                               
!                      =+1 TO PRINT EVERY FIDO ITERATION                        
!                                                                               
!     X( )      --  THE POSITION OF THE MINIMUM OF FOBJ                         
!                                                                               
!                                                                               
!     FIDINT(J) --  RETURN THE COMPUTED LENGTHS OF THE                          
!                   CONFIDENCE HALF-INTERVALS.
!                                                                               
!  FOR A LEAST SQUARES PROBLEM, FOBJ IS EQUAL TO CHI-SQUARE,                    
!  THE WEIGHTED SUM OF SQUARES...                                               
!                                                                               
!             NPTS                                                              
!     FOBJ =  SUM  ((FIT(JPT)-YDATA(JPT))/YSIGMA(JPT))^2                       
!            JPT=1                                                              
!                                                                               
!  THE STANDARD ERRORS YSIGMA( ) MUST BE CORRECTLY SCALED.                      
!  IF THE YSIGMA( ) ARE NOT KNOWN ACCURATELY, NORMALIZE THEM                    
!  SO THAT THE VALUE OF FOBJ AT THE MINIMUM IS EQUAL TO THE                     
!  NUMBER OF DEGREES OF FREEDOM,                                                
!     N.D.F. = (NO. DATA POINTS) - (NO. ADJUSTABLE PARAMETERS)                  
!                                                                               
!  FIDO IS ESSENTIALLY A ROOT FINDING ROUTINE.  IT USES INVERSE                 
!  QUADRATIC INTERPOLATION OR EXTRAPOLATION, INVERSE LINEAR                     
!  INTERPOLATION, AND BISECTION (INTERVAL HALVING), AS                          
!  SUCCESSIVELY MORE DIFFICULTIES ARE ENCOUNTERED IN FINDING                    
!  THE ROOT.                                                                    
!                                                                               
!  PRIMARY REFERENCES....                                                       
!                                                                               
!     W. T. EADIE ET AL., "STATISTICAL METHODS IN EXPERIMENTAL                  
!        EXPERIMENTAL PHYSICS" (AMERICAN ELSEVIER, 1971),                       
!           CHAPTER 9                                                           
!                                                                               
!     P. R. BEVINGTON, "DATA REDUCTION AND ERROR ANALYSIS IN                    
!        THE PHYSICAL SCIENCES" (MCGRAW-HILL, 1969),                            
!           PAGES 242-245                                                       
!                                                                               
!  OTHER REFERENCES....                                                         
!                                                                               
!     H. SCHEFFE, "THE ANALYSIS OF VARIANCE" (WILEY, 1959)                      
!                                                                               
!     H. STONE, DISCUSSION ON PAPER BY E. M. L. BEALE,                          
!        J. ROY. STATIS. SOC. B., V. 22 (1960), P. 41, PP. 84-5                 
!                                                                               
!     G. E. P. BOX AND G. A. COUTIE, PROC. INST. ELEC. ENGRS.                   
!        V. 103, PART B, SUPPL. NO. 1 (1956).                                   
!                                                                               
!     G. W. BOOTH, G. E. P. BOX, M. E. MULLER, AND                              
!        T. I. PETERSON, "FORECASTING BY GENERALIZED REGRESSION                 
!        METHODS; NON-LINEAR ESTIMATION" (PRINCETON/IBM),                       
!        FEBRUARY 1959, IBM MANUAL                                              
!                                                                               
!     D. W. MARQUARDT, "LEAST-SQUARES ESTIMATION OF NONLINEAR                   
!        PARAMETERS", SHARE DISTRIBUTION 3094                                   
!                                                                               
!     D. W. MARQUARDT, R. G. BENNETT, AND E. J. BURRELL, "LEAST                 
!        SQUARES ANALYSIS OF ELECTRON PARAMAGNETI! RESONANCE                    
!        SPECTRA", J. OF MOLECULAR SPECTROSCOPY 7 (1961) 269                    
!                                                                               
!-------------------------------------------------------------------}                 
Label 20,40,70,100,120,140,150,170,190, 210,240,250,270,290,310,330,350,370;
Label 380,400,430;
Var
   XSAVE,XKA,XKB,XBEST: pTab20;  {pointers to Tables of size 20}
      
   FSAVE:Double;                        
       
   RSMALL,CONVRG,BRACK,FACMAX,FRMIN,RZERO,RHALF,UNITR,FPSGSQ,SGN,
   FB,FC,DX,DC,DIFF,FKA,FKB,FAC,XX,FRAC,FBEST: Double;                                                    
      
   ITER,J,JJ,K,KLOSB,MATSAV,MAXTRY,MLIN,MSKSAV,NACTIV,NTRSAV: Integer;

   Function ZSQRT(ARG:double):Double;
   Begin
     ZSQRT := Sqrt(ARG)
   End;

   Function T(ARG:double):Double;
   Begin
     T := ABS(SQRT(ARG-FSAVE)-STDEVS)
   End;

 Begin

   New(XSAVE); New(XKA); New(XKB); New(XBEST);
{ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                 
!                                                                               
!  SET SOME PARAMETERS.                                                         
!                                                                               
!  RSMALL IS USED IN SETTING A PHONY VALUE IF A VALUE OF FOBJ                   
!  IS FOUND THAT IS LESS THAN THE "GLOBAL MINIMUM" VALUE.                       
}
   RSMALL:=0.0001;
                                                                               
{  CONVRG = SATISFACTORY RATE OF CONVERGENCE }
                                                                               
   CONVRG:=0.5;
                                                                               
{  BRACK = FACTOR FOR BRACKETTING SHOTS }
                                                                               
   BRACK:=2.0;
                                                                               
{  MAXTRY = MAXIMUM NUMBER OF ITERATIONS }
                                                                               
   MAXTRY:=20;
                                                                               
{  FACMAX = MAXIMUM FACTOR FOR QUADRATIC INTERPOLATION }
                                                                               
   FACMAX:=4.0;
                                                                               
{  FRMIN = MINIMUM FRACTION FOR LINEAR INTERPOLATION }
                                                                               
   FRMIN:=0.1;
     
   RZERO:=0.0;
      
   RHALF:=0.5;
      
   UNITR:=1.0;
                                                                               
{  SAVE SOME INPUT QUANTITIES, AND INITIALIZE EVERYTHING }                       
                                                                               
   NTRSAV:=NTRACE;
      
   NTRACE:=NTRSET;
      
   MATSAV:=MATRX;
      
   MATRX:=0;
      
   MSKSAV:=MASK[JXFID];
      
   MASK[JXFID]:=1;
      
   NACTIV:=0;
      
   For J:=1 to NV do
   begin      
     XBEST^[J]:=X^[J];
     XSAVE^[J]:=X^[J];
     IF MASK[J] = 0 Then NACTIV:=NACTIV+1
   end;                                                               

   FUNK;

   FSAVE:=FOBJ;

   FBEST:=FOBJ;
      
   FPSGSQ:=FOBJ+Sqr(STDEVS);
      
   X^[JXFID]:=XSAVE^[JXFID]+STDEVS*ERGUES;
                                                                               
{  LOOP OVER THE FIDINT(JJ) }                                                    
                                                                               
   JJ:=1;
   
20:IF ERGUES <> RZERO Then GOTO 40;                                              

   Writeln(fp_out);
   WRITELN(fp_out, ' ERROR IN INPUT TO FIDO...   ERGUES IS EQUAL TO ZERO.');                                          
   Halt;                                                                      
   
40:SGN:=UNITR;
      
   IF ERGUES < RZERO Then SGN:=-UNITR;
      
   IF NTRACF < 0 Then GOTO 70;                                                  
                                              
   Writeln(fp_out);
   Writeln(fp_out,' SUBROUTINE FIDO:');
   Writeln(fp_out,' JXFID =',JXFID:3,'  STDEVS =',STDEVS:12:5,'  TOLFID = ',TOLFID:12:5);
   Writeln(fp_out,' GLOBAL MINIMUM OF FOBJ =', FSAVE:15:7,' IS AT X(JXFID) = ',XSAVE^[JXFID]:15:7);
   Writeln(fp_out,' FIRST GUESS FOR LENGTH OF HALF-INTERVAL = ',ERGUES:12:5);
   Writeln(fp_out);
   Writeln(fp_out,' NV =', NV:11, '  NACTIV =',NACTIV:11);
   
70:KLOSB:=0;
      
   MLIN:=0;
      
   For K:=1 to NV do XKA^[K]:=XSAVE^[K];
     
   FB:=FSAVE;
      
   FKA:=FSAVE;
                                                                               
{  BEGIN THE ITERATION LOOP FOR LOCATING THE PLANE                              
   PERPENDICULAR TO THE X(JXIFD) AXIS IN WHICH THE MINIMUM                      
   VALUE OF FOBJ IS EQUAL TO    FPSGSQ = FSAVE + STDEVS^2,
   WHERE FSAVE IS THE VALUE OF FOBJ AT THE GLOBAL MINIMUM.  }                      
                                                                               
   For ITER:=1 to MAXTRY do
   begin

     FC:=FB;

     IF NACTIV = 0 THEN
       FUNK                                                   
     ELSE                                                                   
       SMPLX;

     IF FOBJ >= FBEST Then GOTO 100;                                            
         
     FBEST:=FOBJ;
         
     For K:=1 to NV do XBEST^[K]:=X^[K];
   
100: DX:=X^[JXFID]-XSAVE^[JXFID];
         
     DC:=FOBJ-FPSGSQ;

     if NTRACF >= 1 then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' ITERATION ',ITER:2,'  X(',JXFID,') = ',X^[JXFID]);
       Writeln(fp_out,' DISTANCE FROM GLOBAL MINIMUM = ',DX:12:5);
       Writeln(fp_out,' MINIMUM FOBJ IN THIS PLANE = ',FOBJ:15:7);
       Writeln(fp_out,' VALUE SOUGHT = ',FPSGSQ:15:7,'  DIFFERENCE = ',DC:15:5);
     end;
         
     IF FOBJ = FSAVE Then GOTO 140;                                            
         
     IF FOBJ < FSAVE Then GOTO 120;                                            

{  TEST FOR CONVERGENCE                                                        
   IF(ABS(FOBJ-FPSGSQ).LE.TOL) ............................ }
                                                                               
     DIFF:=FOBJ-FPSGSQ;

     IF DIFF < RZERO Then DIFF:=-DIFF;
         
     IF DIFF <= TOLFID Then GOTO 330;                                           
         
     IF FOBJ >= FPSGSQ Then GOTO 170;                                           
         
     GOTO 150;                                                              
  
120: DIFF:=FSAVE-FOBJ;
         
     Writeln(fp_out);
     Write(fp_out,' ***** FOBJ LESS THAN VALUE AT INPUT POINT BY ',DIFF:13:5);
     Write(fp_out,' X(J) AT THIS POINT: ');
     For K:=1 to NV do Write(fp_out,' ',X^[K]);
     Writeln(fp_out);
  
140: FOBJ:=FSAVE+RSMALL*(FPSGSQ-FSAVE);

{  CHECK FOR TIGHTER BRACKETTING OF FPSGSQ FROM BELOW }                          
  
150: IF (X^[JXFID]-XKA^[JXFID])*SGN <= RZERO Then GOTO 190;                       
         
     For K:=1 to NV do  XKA^[K]:=X^[K];
        
     FKA:=FOBJ;
         
     GOTO 190;                                                              

{  CHECK FOR BRACKETTING OF FPSGSQ FROM ABOVE }                                  

170: IF KLOSB = 1 THEN                                                    
            
       IF (X^[JXFID]-XKB^[JXFID])*SGN >= RZERO Then GOTO 190;                    
 
     KLOSB:=1;
         
     For K:=1 to NV do  XKB^[K]:=X^[K];

     FKB:=FOBJ;
  
190: IF T(FOBJ) < T(FC) Then FB:=FOBJ;

{  CHECK THE RATE OF CONVERGENCE.  IF IT IS SATISFACTORY, AND
   IF LINEAR INTERPOLATION HAS BEEN USED, USE IT AGAIN.       }
         
     IF (ITER >= 2) AND (T(FOBJ) > CONVRG*T(FC)) Then GOTO 210;                  
         
     IF MLIN = 1 Then GOTO 250;                                                

{  USE INVERSE QUADRATIC INTERPOLATION }                                         
         
     FAC:=STDEVS/ZSQRT(FOBJ-FSAVE);

{  FAC:=AMIN1(FACMAX,AMAX1(FAC,1./FACMAX)) ................... }
         
     IF FAC < UNITR/FACMAX Then FAC:=UNITR/FACMAX;
         
     IF FAC > FACMAX Then FAC:=FACMAX;
         
     XX:=XSAVE^[JXFID]+FAC*(X^[JXFID]-XSAVE^[JXFID]);

{  CHECK THAT THE PROPOSED POINT IS INSIDE THE BRACKETTED INTERVAL }
         
     IF (XX-XKA^[JXFID])*SGN <= RZERO Then GOTO 210;                             
         
     IF KLOSB = 1 THEN                                                    
            
       IF (XX-XKB^[JXFID])*SGN >= RZERO Then GOTO 240;                          

     For K:=1 to NV do
       IF(K = JXFID) OR (MASK[K] = 0) Then                                    
          X^[K]:=XSAVE^[K] + FAC*(X^[K]-XSAVE^[K]);
         
     GOTO 310;                                                              
  
210: IF KLOSB = 1 Then GOTO 240;                                               

{  CONVERGENCE IS POOR, AND FPSGSQ HAS NOT YET BEEN BRACKETTED.
   TRY TO BRACKET IT. }
         
     FRAC:=BRACK;
         
     IF FOBJ >= FPSGSQ Then FRAC:=UNITR/BRACK;
         
     For K:=1 to NV do
       IF (K = JXFID) OR (MASK[K] = 0) Then X^[K]:=XSAVE^[K]+FRAC*(X^[K]-XSAVE^[K]);
         
     IF NTRACF >= 1 Then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' BRACKETTING SHOT....')
     end;

     GOTO 310;                                                              

240: IF MLIN = 1 Then GOTO 270;
                                                                              
{  TRY LINEAR INTERPOLATION BETWEEN THE TWO BRACKETTING POINTS }
 
250: MLIN:=1;
         
     FRAC:=(FPSGSQ-FKA)/(FKB-FKA);

{  FRAC:=AMAX1(FRMIN,AMIN1(1.0-FRMIN,FRAC)) ............. }
         
     IF FRAC > UNITR-FRMIN Then FRAC:=UNITR-FRMIN;
         
     IF FRAC < FRMIN Then FRAC:=FRMIN;
 
     IF NTRACF >= 1 Then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' LINEAR INTERPOLATION....')
     end;

     GOTO 290;                                                              

{  CONVERGENCE IS POOR, AND LINEAR INTERPOLATION HAS BEEN USED.
   BISECT THE BRACKETTED INTERVAL.  }
  
270: FRAC:=RHALF;
         
     IF NTRACF >= 1 Then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' INTERVAL BISECTED....')
     end;

290: For K:=1 to NV do
       IF (K = JXFID) OR (MASK[K] = 0) Then
         X^[K]:=XKA^[K] + FRAC*(XKB^[K]-XKA^[K]);

310: end;

{  END OF ITERATION LOOP.  THE ITERATION FAILED TO CONVERGE }                    
      
   Writeln(fp_out);
   Writeln(fp_out,' CONVERGENCE FAILURE IN SUBROUTINE FIDO.');
     
   FIDINT[JJ]:=RZERO;
      
   GOTO 350;                                                                 

{  CONVERGENCE ACHIEVED ...  A SATISFACTORY PLANE HAS BEEN LOCATED }

330: FIDINT[JJ]:=X^[JXFID]-XSAVE^[JXFID];
      
   IF NTRACF >= 0 Then
   begin
     Writeln(fp_out);
     Writeln(fp_out,' THE LENGTH OF THE CONFIDENCE HALF-INTERVAL FOR');
     Writeln(fp_out,' X(',JXFID,') IS ',FIDINT[JJ]:13:5,' (STDEVS = ',STDEVS:12:5,')')
   end;
  
350: JJ:=JJ+1;
      
   IF JJ > MAXIND Then GOTO 370;                                                
                                                                             
{  REFLECT X THROUGH XSAVE, AND SEARCH FOR THE OTHER HALF-INTERVAL }
     
   For K:=1 to NV do
     IF (K = JXFID) OR (MASK[K] = 0) Then  X^[K]:=XSAVE^[K]+(XSAVE^[K]-X^[K]);

   ERGUES:=X^[JXFID]-XSAVE^[JXFID];
      
   GOTO 20;                                                                  

{  END OF THE JJ LOOP.  PRINT SUMMARY RESULTS }                                  
                                                                                
370:IF (MAXIND < 2) OR (FIDINT[1] = RZERO) Then GOTO 400;                         
      
   IF (FIDINT[1] < RZERO) AND (FIDINT[2] <= RZERO) Then GOTO 400;                 
      
   IF (FIDINT[1] > RZERO) AND (FIDINT[2] >= RZERO) Then GOTO 400;                 
      
   IF (FIDINT[1] > RZERO) AND (FIDINT[2] < RZERO) Then GOTO 380;                 
      
   XX:=FIDINT[1];
      
   FIDINT[1]:=FIDINT[2];
      
   FIDINT[2]:=XX;

380: XX:=-FIDINT[2];

   { PRINT A SUMMARY OF RESULTS FROM FIDO... }
   Writeln(fp_out,'                       + ',FIDINT[1]:8:4);
   Writeln(fp_out,'  X(',JXFID,') = ',XSAVE^[JXFID]:13:6,'                STDEVS = ',STDEVS:8:5);
   Writeln(fp_out,'                       - ',XX:8:4);

{  RESTORE THE SAVED ENTRY VALUES }                                              
  
400: For J:=1 to NV do X^[J]:=XBEST^[J];
         
   FUNK;                                                              
      
   IF FOBJ >= FSAVE Then GOTO 430;

   Writeln(fp_out);
   Writeln(fp_out,' SUBROUTINE FIDO FOUND A BETTER MINIMUM.');
   Writeln(fp_out,' OLD FOBJ = ',FSAVE:18:10,' NEW FOBJ = ',FOBJ:18:10);
   Write(fp_out,' X(J) = ');
   For J:=1 to NV do Write(fp_out,' ',X^[J]:18:10);
   Writeln(fp_out);

430: MASK[JXFID]:=MSKSAV;
      
   NTRACE:=NTRSAV;
      
   MATRX:=MATSAV;

   Dispose(XSAVE); Dispose(XKA); Dispose(XKB); Dispose(XBEST)

End; {FIDO}                                                                       

Procedure SIBEG; Forward;
Procedure SIFUN; Forward;
Procedure DATSW(NSSW:Integer; Var JUMP:Integer); Forward;

Procedure SMPLX;
{---------------------------------------------------------------
!  SIMPLEX 2.12          DECEMBER 1991                                          
!                                                                               
!  A.N.S.I. STANDARD FORTRAN 77                                                 
!                                                                               
!  COPYRIGHT (C) 1965, 1975, 1991 J. P. CHANDLER                                
!     (PRESENT ADDRESS ...                                                      
!        COMPUTER SCIENCE DEPARTMENT,                                           
!        OKLAHOMA STATE UNIVERSITY,                                             
!        STILLWATER, OKLAHOMA 74078                                             
!        (405)-744-5676                )                                        
!                                                                               
!  SIMPLEX FINDS LOCAL MINIMA OF A SMOOTH FUNCTION OF SEVERAL                   
!  PARAMETERS.  IT WILL ALSO HANDLE SOME NON-SMOOTH FUNCTIONS.                  
!  THE METHOD IS OFTEN VERY SLOW IF THE NUMBER OF PARAMETERS IS                 
!  AT ALL LARGE (GREATER THAN ABOUT SIX TO EIGHT).                              
!                                                                               
!     "A SIMPLEX METHOD FOR FUNCTION MINIMIZATION",                             
!     J. A. NELDER AND R. MEAD,                                                 
!     THE COMPUTER JOURNAL 7 (1965) 308-313                                     
!                                                                               
!  FOR APPLICATIONS, SEE                                                        
!                                                                               
!     "'DIRECT SEARCH' SOLUTION OF                                              
!        NUMERICAL AND STATISTICAL PROBLEMS",                                   
!     ROBERT HOOKE AND T. A. JEEVES, JOURNAL OF THE                             
!     ASSOCIATION FOR COMPUTING MACHINERY 8 (1961) 212-229                      
!                                                                               
!     "THE NELDER-MEAD SIMPLEX PROCEDURE FOR FUNCTION                           
!        MINIMIZATION",                                                         
!     D. M. OLSSON AND L. S. NELSON,                                            
!     TECHNOMETRICS 17 (1975) 45-51                                             
!                                                                               
!                                                                               
!  SIMPLEX 2.9 IS AVAILABLE FROM THE                                            
!     QUANTUM CHEMISTRY PROGRAM EXCHANGE                                        
!     DEPT. OF CHEMISTRY, INDIANA UNIVERSITY                                    
!     BLOOMINGTON, INDIANA 47401                                                
!                                                                               
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                 
!                                                                               
!  INPUT QUANTITIES.....     NV,NTRACE,MASK(*),X(*),                          
!                            XMAX(*),XMIN(*),DELTX(*),                          
!                            DELMIN(*),NFMAX,NFLAT,NXTRA,KW                     
!                                                                               
!  OUTPUT QUANTITIES....  X(*), FOBJ, KFLAG, NOREP                         
!                                                                               
!     NV         --  THE NUMBER OF PARAMETERS, X(J)                             
!                                                                               
!     NTRACE     --  :=0 FOR NORMAL OUTPUT,
!                    :=+1 FOR TRACE OUTPUT,
!                    :=+2 FOR DEBUG OUTPUT,
!                    :=-1 FOR NO OUTPUT EXCEPT ERROR MESSAGES
!                                                                               
!     MATRX      --  USED BY STEPIT PROPER BUT NOT BY SIMPLEX                   
!                                                                               
!     FOBJ       --  THE VALUE OF THE FUNCTION TO BE MINIMIZED                  
!                                                                               
!     MASK(J)    --  NONZERO IF X(J) IS TO BE HELD FIXED                        
!                                                                               
!     X(J)       --  THE J-TH PARAMETER                                         
!                                                                               
!     XMAX(J)    --  THE UPPER LIMIT ON X(J)                                    
!                                                                               
!     XMIN(J)    --  THE LOWER LIMIT ON X(J)                                    
!                                                                               
!     DELTX(J)   --  THE INITIAL STEP SIZE FOR X(J)                             
!                                                                               
!     DELMIN(J)  --  THE LOWER LIMIT (CONVERGENCE TOLERANCE) ON                 
!                       THE STEP SIZE FOR X(J)                                  
!                                                                               
!     ERR(J,K)   --  SCRATCH STORAGE                                            
!                                                                               
!     NFMAX      --  THE MAXIMUM NUMBER OF FUNCTION                             
!                       COMPUTATIONS TO BE ALLOWED                              
!                                                                               
!     NFLAT      --  NONZERO IF THE SEARCH IS TO TERMINATE WHEN                 
!                       ALL TRIAL STEPS GIVE IDENTICAL FUNCTION                 
!                       VALUES.                                                 
!                       THE RECOMMENDED VALUE OF NFLAT IS                       
!                       USUALLY  NFLAT:=1 .
!                                                                               
!     JVARY      --  USED BY STEPIT AND STP BUT NOT BY SIMPLEX                  
!                       (SIMPLEX SETS JVARY TO ZERO                             
!                       PERMANENTLY)                                            
!                                                                               
!     NXTRA      --  NUMBER OF EXTRA POINTS TO BE ADDED TO                      
!                       THE SIMPLEX                                             
!                       (NXTRA.GT.0 CAUSES A MORE THOROUGH                      
!                       SEARCH)                                                 
!                                                                               
!     KFLAG      --  RETURNED .GT. ZERO FOR A NORMAL EXIT,                      
!                       RETURNED .LT. ZERO FOR AN ABNORMAL EXIT                 
!                                                                               
!     NOREP      --  RETURNED .GT. ZERO IF THE FUNCTION WAS NOT                 
!                       REPRODUCIBLE                                            
!                                                                               
!     KERFL      --  NOT USED BY THIS ROUTINE                                   
!                                                                               
! ---------------------------------------------------------------}                 
                                                                            
{  SIMPLEX SHOULD USUALLY BE RUN USING A FLOATING POINT                         
   PRECISION OF AT LEAST TEN SIGNIFICANT DIGITS.  ON MOST                       
   COMPUTERS, THIS REQUIRES THE USE OF DOUBLE PRECISION. }
Label 10,40,50,70,100,130,140,150,160,170,190,220,250,280,310,340,360,380;
Label 390,400,420,440,450,460,470,490,510,530,550,560,570,590,610,640,660;
Label 690,720;                          
Var      
   RZERO,UNITR,RTWO,SCALK,XS,FSAVE,DX,DZ,ZMAX,ZMIN,FSTAR,ZBARJ,ZJK,
   ZJJH,XJ,ZKJ,ZKJL,XK: Double;                                                         

   J,JDIFF,JH,JHSAV,JHTIE,JL,JLOW,JS,JUMP,K,LATER,NOW,KONTR: Integer;
      
   ZBAR, ZSTAR: pTab20;
                                                                               
{  SIMPLEX CALLS NO FUNCTIONS, EITHER EXTERNAL OR INTRINSIC.                    
   THE SUBROUTINES CALLED ARE FUNK, SIBEG, SIFUN, AND DATSW.                    
   SIMPLEX TERMINATES IF SENSE SWITCH NUMBER -NSSW- IS ON.                      
   THE STATEMENT    CALL DATSW(NSSW,JUMP)    RETURNS JUMP:=1 IF
   SENSE SWITCH NUMBER -NSSW- IS ON, AND JUMP:=2 IF IT IS OFF.
   IF NO SENSE SWITCH IS TO BE USED, THE USER SHOULD SUPPLY                     
   A DUMMY SUBROUTINE FOR DATSW.                                                
   THIS SUBROUTINE CONTAINS NO REAL OR DOUBLE PRECISION                         
   CONSTANTS.         }                                                          

Begin      

   New(ZBAR); New(ZSTAR);

   SIBEG;
      
   IF KFLAG < 0 Then GOTO 720;                                                  
      
   FSAVE:=FOBJ;
      
   KERFL:=0;

{  SET FIXED QUANTITIES....                                                     
                                                                               
   METHD:=1 TO USE THE METHOD OF NELDER AND MEAD,
   METHD:=2 TO USE A MODIFIED METHOD.
   METHD:=1 CAN CAUSE THE BEST KNOWN POINT TO BE DISCARDED,
   AND HENCE IS NOT RECOMMENDED.                                                
   METHD:=2 IS RECOMMENDED.         }
     
   METHD:=2;
      
   RZERO:=0.0;
      
   UNITR:=1.0;
      
   RTWO:=2.0;
      
   JL:=NVPLUS;
      
   LATER:=0;
                                                                  
{  BEGIN THE NEXT ITERATION.                                                    
   FIND JH AND JL, THE INDICES OF THE POINTS WITH THE HIGHEST                   
   AND LOWEST FUNCTION VALUES, RESPECTIVELY.  }                                    
  
10:JH:=JL;
      
   For J:=1 to NVPLUS do
   begin      
     IF FZ^[J] > FZ^[JH] Then JH:=J;
     IF FZ^[J] < FZ^[JL] Then JL:=J;
   end;

{  CHECK FOR POSSIBLE TERMINATION }                                              
       
   IF KFLAG <> 0 Then GOTO 690;                                                  
      
   IF(JH = JL) AND (NFLAT > 0) Then GOTO 660;                                   

{  CHECK FOR TIES -- MORE THAN ONE POINT IN THE SIMPLEX HAVING                  
   THE HIGHEST FUNCTION VALUE.  }                                                  
      
   JHTIE:=0;
      
   For J:=1 to NVPLUS do
     IF(J <>JH) AND (FZ^[J] >= FZ^[JH]) Then JHTIE:=7;

   IF (JH <> JL) AND (JHTIE = 0) Then GOTO 70;                                    

{  THERE IS A TIE FOR HIGHEST FUNCTION VALUE, CHOOSE H FAR FROM L. }                                                         
      
   DX:=RZERO;
      
   JHSAV:=JH;
      
   For J:=1 to NVPLUS do
   begin      
     IF(J = JL) OR (FZ^[J] < FZ^[JH]) Then GOTO 50;                              

     For K:=1 to NV do
     begin       
       IF MASK[K] <> 0 Then GOTO 40;                                           
       SCALK:=X^[K];
       IF SCALK = RZERO Then SCALK:=DELTX^[K];
       DZ:=(Z^[K,J]-Z^[K,JL])/SCALK;
       IF DZ < RZERO Then DZ:=-DZ;
       IF DZ <= DX Then GOTO 40;                                               
       JH:=J;
       DX:=DZ;
40:  end;                                                            
50:end;                                                               
      
   IF NTRACE >= 2 Then
   begin
     Writeln(fp_out);
     Writeln(fp_out,' TIE BREAKING....     JHSAV =',JHSAV:3,'  JH =',JH:3)
   end;

IF DX <= RZERO Then GOTO 640;                                                 
   
70: IF NTRACE < 2 Then GOTO 100;                                                 

   Writeln(fp_out);
   Writeln(fp_out,' FZ(',JL,') = ',FZ^[JL]:24:16);
   Write(fp_out,' Z(J,JL) = ');
   For J:=1 to NV do Write(fp_out,' ',Z^[J,JL]:24:16);
   Writeln(fp_out);

   Writeln(fp_out);
   Writeln(fp_out,' FZ(',JH,') = ',FZ^[JH]:24:16);
   Write(fp_out,' Z(J,JH) = ');
   For J:=1 to NV do Write(fp_out,' ',Z^[J,JH]:24:16);
   Writeln(fp_out);

{  ADD EXTRA POINTS TO THE SIMPLEX, IF DESIRED.                                 
   ANY EXTRA POINTS WILL BE SUPERIMPOSED ON THE HIGHEST POINT,                  
   P(JH), SO THAT THEY WILL SPLIT AS SOON AS POSSIBLE.  }                          
  
100:IF LATER <> 0 Then  GOTO 140;                                                  
      
   LATER:=7;
      
   IF NVPLUS+NXTRA > MAXPT Then  NXTRA:=MAXPT-NVPLUS;
      
   IF NXTRA <= 0 Then GOTO 130;                                                  
      
   JLOW:=NVPLUS+1;
      
   NVPLUS:=NVPLUS+NXTRA;
      
   For J:=JLOW to NVPLUS do
   begin      
     FZ^[J]:=FZ^[JH];
     For K:=1 to NV do  Z^[K,J]:=Z^[K,JH]
   end; 
  
130:FSAVE:=FZ^[JL];

{  CALCULATE PBAR, THE CENTROID OF ALL POINTS EXCEPT P(JH).                     
   THE MEAN VALUE IS COMPUTED USING A STABLE UPDATE FORMULA BY                  
   D. H. D. WEST, COMMUNICATIONS OF THE A.C.M. 22 (1979) 532-535.  }                                                        
  
140:For J:=1 to NV do
    begin      

     ZBARJ:=X^[J];
         
     IF MASK[J] <> 0 Then GOTO 190;                                             
         
     ZMAX:=Z^[J,JL];
         
     ZMIN:=Z^[J,JL];
         
     XS:=RZERO;

     For K:=1 to NVPLUS do
     begin
       IF K = JH Then GOTO 150;                                               
       ZJK:=Z^[J,K];
       IF ZJK > ZMAX Then ZMAX:=ZJK;
       IF ZJK < ZMIN Then ZMIN:=ZJK;
       XS:=XS+UNITR;
       ZBARJ:=ZBARJ+(ZJK-ZBARJ)/XS;
150: end;                                                            

{  CHECK THE ROUNDING.                                                          
   ZBARJ MUST LIE IN THE CLOSED INTERVAL (ZMIN,ZMAX). }                           
 
     IF ZBARJ <= ZMAX Then GOTO 160;                                            
         
     ZBARJ:=ZMAX;
         
     GOTO 170;                                                              
  
160: IF ZBARJ >= ZMIN Then GOTO 190;                                            
         
     ZBARJ:=ZMIN;
  
170: NOW:=1;
         
     IF NTRACE >= 1 Then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' ROUNDING OF X(',J,') CORRECTED AT CHECKPOINT ',NOW,' WITH NF =',NF)
     end;

190: ZBAR^[J]:=ZBARJ
  
   end; {J loop}                                                                

   If NTRACE >= 2 Then
   begin
     Writeln(fp_out);
     Write(fp_out,' ZBAR(J) = ');
     For J:=1 to NV do Write(fp_out,' ',ZBAR^[J]:24:16);
     Writeln(fp_out)
   end;

{  ATTEMPT A REFLECTION.  FORM P* .                                             
   THE FORMS USED BELOW FOR REFLECTION, EXPANSION, AND                          
   CONTRACTION ALL HAVE OPTIMAL ROUNDOFF PROPERTIES.   }                            

   JDIFF:=0;

   For J:=1 to NV do
   begin
     XJ:=X^[J];
     IF MASK[J] <> 0 Then GOTO 220;                                             
     XJ:=ZBAR^[J]+ALPHA*(ZBAR^[J]-Z^[J,JH]);
     IF XJ <> Z^[J,JH] Then JDIFF:=7;
     X^[J]:=XJ;
220: ZSTAR^[J]:=XJ
   end;                                                               

   IF JDIFF <> 0 Then GOTO 250;                                                  
      
   IF NTRACE >= 1 Then
     Writeln(fp_out,' ZSTAR = Z(JH),  TREAT AS A CONTRACTION FAILURE');

   GOTO 440;                                                                 

250: SIFUN;                                                         
      
   FSTAR:=FOBJ;
      
   KONTR:=0;
      
   IF FSTAR < FZ^[JL] Then GOTO 280;                                             
  
   IF NTRACE >= 1 Then
   begin
     Writeln(fp_out);
     Writeln(fp_out,' FOBJ =',FSTAR:24:16,'   REFLECTION FAILED,   NF =',NF)
   end;
      
   For J:=1 to NVPLUS do
   begin
{  HERE WE DEVIATE FROM THE FLOWCHART IN THE ARTICLE BY NELDER                  
   AND MEAD.  THEY USE   FSTAR.LE.FZ(J) .                                       
   THAT IS NOT CONSISTENT WITH THE TEXT OF THE ARTICLE, AND CAN                 
   CREATE AN INFINITE LOOP OF REFLECTIONS.  }                                      
         
     IF(J <> JH) AND (FSTAR < FZ^[J]) Then GOTO 340                             
  
   end;                                                               

{  HERE WE DEVIATE FROM THE FLOWCHART IN THE ARTICLE BY NELDER                  
   AND MEAD.  THEY USE   FSTAR.LE.FZ(JH)  (SEE COMMENT ABOVE). }                  

   KONTR:=1;
      
   IF FSTAR < FZ^[JH] Then GOTO 340;                                             
      
   GOTO 380;                                                                 
  
280: IF NTRACE < 1 Then GOTO 310;

   Writeln(fp_out);
   Writeln(fp_out,' FOBJ =',FSTAR:24:16,'   REFLECTION SUCCEEDED,   NF =',NF);
   Write(fp_out,' X = ');
   For J:=1 to NV do Write(fp_out,' ',X^[J]:15:7);
   Writeln(fp_out);

{  THE REFLECTED VALUE FSTAR IS LESS THAN THE LOWEST VALUE,                     
   FZ(JL), IN THE SIMPLEX.                                                      

   ATTEMPT AN EXPANSION.  FORM P** . }
  
310: For J:=1 to NV do
       IF MASK[J] = 0 Then X^[J]:=X^[J]+(GAMMA-UNITR)*(X^[J]-ZBAR^[J]);

   SIFUN;                                                         

{  CHOOSE ONE OF TWO STRATEGIES FOR ACCEPTING THE EXPANSION POINT P** .  }

   IF((METHD = 1) AND (FOBJ < FZ^[JL])) OR ((METHD = 2) AND (FOBJ < FSTAR)) Then GOTO 490;                            
      
   IF NTRACE >= 1 Then
   begin
     Writeln(fp_out);
     Writeln(fp_out,' FOBJ =',FOBJ:24:16,',  EXPANSION FAILED.')
   end;

{  THE REFLECTION FAILED BUT THE REFLECTED POINT WAS BETTER                     
   THAN SOME OTHER POINT, OR ELSE THE REFLECTION SUCCEEDED                      
   BUT THE EXPANSION FAILED.                                                    
                                                                               
   ACCEPT THE REFLECTED POINT.  REPLACE P(JH) BY P* . }                           
  
340:IF NTRACE < 1 Then GOTO 360;                                                 

   Writeln(fp_out);
   Writeln(fp_out,' ACCEPT REFLECTED POINT:');
   For J:=1 to NV do Write(fp_out,' ',ZSTAR^[J]:15:7);
   Writeln(fp_out);

360: For J:=1 to NV do  Z^[J,JH]:=ZSTAR^[J];
      
    FZ^[JH]:=FSTAR;

{  IF THE REFLECTION FAILED AND THE REFLECTED POINT WAS BETTER                  
   THAN ONLY P(JH), CONTRACT THE SIMPLEX.  }                                       
      
    IF KONTR = 0 Then GOTO 530;                                                  

{  ATTEMPT A CONTRACTION.  FORM P** }                                           
  
380: JDIFF:=0;
      
    For J:=1 to NV do
    begin     
      IF MASK[J] <> 0 Then GOTO 400;                                             
         
      ZBARJ:=ZBAR^[J];
         
      ZJJH:=Z^[J,JH];
         
      XJ:=ZBARJ+BETA*(ZJJH-ZBARJ);

{  CHECK THE ROUNDING.                                                          
   XJ MIGHT HAVE BEEN ROUNDED UP TO ZJJH, WHICH COULD CREATE                    
   AN INFINITE LOOP.                                                            
   XJ MUST LIE IN THE CLOSED INTERVAL (ZJJH,ZBARJ), AND                         
   XJ MUST NOT BE EQUAL TO ZJJH UNLESS ZJJH.EQ.ZBARJ .                          
   EQUIVALENTLY, AN OPEN INTERVAL (ZJJH,ZBARJ) MUST EXIST AND                   
   XJ MUST LIE IN IT, OR XJ MUST BE EQUAL TO ZBARJ. }                             
         
      IF ((XJ > ZJJH) AND (XJ < ZBARJ)) OR ((XJ > ZBARJ) AND (XJ < ZJJH)) Then GOTO 390;                                
      IF XJ = ZBARJ Then GOTO 390;                                              
         
      NOW:=2;
         
      IF NTRACE >= 1 Then
      begin
        Writeln(fp_out);
        Writeln(fp_out,' ROUNDING OF X(',J,') CORRECTED AT CHECKPOINT ',NOW,' WITH NF =',NF)
      end;
         
      XJ:=ZBARJ;

390:  IF XJ <> ZJJH Then JDIFF:=7;
         
      X^[J]:=XJ;
  
400: end; {J loop}                                                              
      
   IF JDIFF = 0 Then GOTO 420;                                                  
      
   SIFUN;                                                         
      
   IF FOBJ > FZ^[JH] Then GOTO 420;
      
   IF NTRACE < 1 Then GOTO 420;                                                 

   Writeln(fp_out);
   Writeln(fp_out,' FOBJ =',FOBJ:24:16,'   CONTRACTION SUCCEEDED,   NF =',NF);
      
   GOTO 510;                                                                 

420: IF NTRACE >= 1 Then
     begin
       Writeln(fp_out);
       Writeln(fp_out,' FOBJ =',FOBJ:24:16,',  CONTRACTION FAILED.')
     end;

440: JS:=JL;

{  REPLACE ALL P(J) BY (P(J)+P(JL))/2.0 }                                       
      
   For J:=1 to NVPLUS do
   begin
         
     IF J = JL Then GOTO 470;                                                  

     For K:=1 to NV do
     begin            

       IF MASK[K] <> 0 Then GOTO 460;                                          
            
       ZKJ:=Z^[K,J];
            
       ZKJL:=Z^[K,JL];
            
       XK:=ZKJL+(ZKJ-ZKJL)/RTWO;
 
{  CHECK THE ROUNDING.                                                          
   XK MIGHT HAVE BEEN ROUNDED UP TO ZKJ, WHICH COULD CREATE                     
   AN INFINITE LOOP.                                                            
   XK MUST LIE IN THE CLOSED INTERVAL (ZKJ,ZKJL), AND                           
   XK MUST NOT BE EQUAL TO ZKJ UNLESS ZKJ.EQ.ZKJL .                             
   EQUIVALENTLY, AN OPEN INTERVAL (ZKJ,ZKJL) MUST EXIST AND                     
   XK MUST LIE IN IT, OR XK MUST BE EQUAL TO ZKJL }                             
            
       IF((XK > ZKJ) AND (XK < ZKJL)) OR ((XK > ZKJL) AND (XK < ZKJ)) Then GOTO 450;
       IF XK = ZKJL Then GOTO 450;                                            
            
       NOW:=3;
            
       IF NTRACE >= 1 Then
       begin
         Writeln(fp_out);
         Writeln(fp_out,' ROUNDING OF X(',J,') CORRECTED AT CHECKPOINT ',NOW,' WITH NF =',NF)
       end;
            
       XK:=ZKJL;

450:   Z^[K,J]:=XK;
            
       X^[K]:=XK;
  
460: end; {K loop}                                                            
         
     SIFUN;                                                    
         
     IF FOBJ < FZ^[JS] Then JS:=J;
         
     FZ^[J]:=FOBJ;
  
470: end; {J loop}                                                               
      
   IF JS = JL Then GOTO 530;                                                    
      
   JL:=JS;

   if NTRACE >=1 Then
   begin 
     Writeln(fp_out);
     Writeln(fp_out,' FZ(JS) = ', FZ^[JS]:24:16);
     Write(fp_out,' Z = ');
     For K:=1 to NV do Write(fp_out,' ',Z^[K,JS]:15:7);
     Writeln(fp_out)
   end;

   GOTO 530;                                                                 

490: IF NTRACE < 1 Then GOTO 510;

   Writeln(fp_out);
   Writeln(fp_out,' FOBJ =',FOBJ:24:16,',  EXPANSION SUCCEEDED.');
      
   Write(fp_out,' X = ');
   For J:=1 to NV do Write(fp_out,' ',X^[J]:15:7);
   Writeln(fp_out);

{  REPLACE P(JH) BY P* OR P** }                                                 
  
510: For J:=1 to NV do Z^[J,JH]:=X^[J];
      
   FZ^[JH]:=FOBJ;

{  TEST FOR CONVERGENCE }                                                        
  
530:For J:=1 to NV do
    begin     
      IF MASK[J] <> 0 Then GOTO 550;                                             
      ZMAX:=Z^[J,1];
      ZMIN:=Z^[J,1];
      For K:=2 to NVPLUS do
      begin
        ZJK:=Z^[J,K];
        IF ZJK > ZMAX Then ZMAX:=ZJK;
        IF ZJK < ZMIN Then ZMIN:=ZJK
      end;                                                            
      DZ:=ZMAX-ZMIN;
      IF DZ > DELMIN^[J] Then GOTO 560;                                          
550:end;                                                               
      
    GOTO 640;                                                                 

{  RETURN IF THE SENSE SWITCH IS ON }                                            

560: JUMP:=2;
      
    DATSW(NSSW,JUMP);
      
    IF JUMP <= 1 Then GOTO 570;                                                   

    IF NF > NFMAX Then GOTO 590;                                                 
      
    GOTO 10;                                                                  
 
570: KFLAG:=-3;
      
    IF NTRACE >= -1 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' ABNORMAL TERMINATION...   TERMINATED BY OPERATOR VIA SENSE SWITCH ',NSSW)
    end;
      
    GOTO 610;                                                                 

590: KFLAG:=-2;
      
    IF NTRACE >= -1 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' ABNORMAL TERMINATION...   MORE THAN ',NFMAX,' CALLS TO THE FOBJ SUBROUTINE.')
    end;

610: IF NTRACE < -1 Then GOTO 690;                                                

    Writeln(fp_out);
    For J:=1 to NVPLUS do
    begin
      Writeln(fp_out,' SIMPLEX POINT ',J,'  FOBJ=',FZ^[J]:24:16);
      Write(fp_out,' X(J) = ');
      For K:=1 to NV do Write(fp_out,' ',Z^[K,J]:15:7);
      Writeln(fp_out)
    end;                                                               
      
    GOTO 690;                                                                 
  
640:KFLAG:=1;
      
    IF NTRACE >= 0 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' TERMINATED WHEN THE DIMENSIONS OF THE SIMPLEX BECAME AS SMALL');
      Writeln(fp_out,' AS THE DELMIN(J).')
    end;
      
    GOTO 690;                                                                 

660:KFLAG:=2;
      
    IF NTRACE >= 0 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' TERMINATED WHEN THE FUNCTION VALUES AT ALL POINTS OF THE SIMPLEX');
      Writeln(fp_out,' WERE EXACTLY EQUAL.')
    end;

{  GO BACK AND COMPUTE JL AND JH ONE LAST TIME }                                 
   GOTO 10;                                                                  
{  FINISH UP AND EXIT }                                                          
 
690:For J:=1 to NV do X^[J]:=Z^[J,JL];
      
    SIFUN;                                                         
      
    IF (FOBJ <= FSAVE) AND (FOBJ = FZ^[JL]) Then GOTO 720;                          
      
    NOREP:=NOREP+2;
      
    IF NTRACE >= -1 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' WARNING....  FOBJ IS NOT A REPRODUCIBLE FUNCTION OF X(J).');
      Writeln(fp_out,' NF=',NF,' FSAVE=',FSAVE:24:16,' FZ(JL)=',FZ^[JL]:24:16,' FOBJ=',FOBJ:24:16) 
    end;
  
720:IF NTRACE >= 0 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out);
      Writeln(fp_out);
      Writeln(fp_out,' ',NF,' FUNCTION COMPUTATIONS:');
      Writeln(fp_out);
      Writeln(fp_out,' FINAL VALUE OF FOBJ = ',FOBJ:24:16);
      Writeln(fp_out);
      Write(fp_out,' FINAL VALUES OF X(J): ');
      For J:=1 to NV do Write(fp_out,' ',X^[J]:15:7);
      Writeln(fp_out);
    end;
    Writeln(fp_out);

    Dispose(ZBAR); Dispose(ZSTAR)

End; {SMPLX}


Procedure SIBEG;                                                   
{--------------------------------------------------------------
!  SIBEG 1.3          OCTOBER 1991                                              
!                                                                               
!  A.N.S.I. STANDARD FORTRAN 77                                                 
!                                                                               
!  COPYRIGHT (C) 1965, 1975, 1990 J. P. CHANDLER,                               
!     COMPUTER SCIENCE DEPARTMENT,                                              
!     OKLAHOMA STATE UNIVERSITY                                                 
!                                                                               
!  SIBEG SETS DEFAULT VALUES AND PRINTS INITIAL OUTPUT                          
!  FOR SIMPLEX.                                                                 
!                                                                               
! ------------------------------------------------------------                 
!                                                                               
!  INPUT QUANTITIES.....  FUNK,X(*),XMAX(*),XMIN(*),DELTX(*),                   
!                              DELMIN(*),NV,NTRACE,MASK(*),                     
!                              NFMAX,NFLAT,NXTRA,KW                             
!                                                                               
!  OUTPUT QUANTITIES....  FZ,HUGE,ALPHA,BETA,GAMMA,MAXPT,                       
!                              NVPLUS,NSSW,NF,KFLAG,NOREP,                      
!                              DELTX(*), AND DELMIN(*)                          
!-------------------------------------------------------------}
Label 20,30,50,60,70,80,90,100,180,200,230,260;

Var
    J,JUMP,NVMAX,NACTIV: Integer;
    DELDF,FSAVE,RZERO,XPLUS,XS: Double;

Begin

{ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                 
!                                                                               
!  SET FIXED QUANTITIES ....                                                    
!                                                                               
   NSSW = SENSE SWITCH NUMBER FOR SEARCH TERMINATION
   (IRRELEVANT IF A DUMMY DATSW ROUTINE IS USED)  }                             
                                                                               
   NSSW:=6;
                                                                               
{  HUGE := A VERY LARGE REAL NUMBER
     (DEFAULT VALUE FOR XMAX AND -XMIN, AND                                    
      OUT-OF-BOUNDS VALUE FOR FOBJ)  }                                             
                                                                               
   HUGE:=1E35;
                                                                               
{  NVMAX IS THE MAXIMUM VALUE OF NV.                                            
   NVMAX IS ALSO THE DIMENSION OF X, XMAX, XMIN, DELTX, DELMIN,                 
   MASK, ZBAR, ZSTAR, AND THE FIRST DIMENSION OF ERR AND Z. }                     
                                                                               
   NVMAX:=20;
                                                                               
{  MAXPT IS THE MAXIMUM NUMBER OF POINTS IN THE SIMPLEX.                        
   MAXPT IS ALSO THE DIMENSION OF FZ, AND THE SECOND DIMENSION                  
   OF ERR AND Z.                                                                
   MAXPT MUST BE .GE. (NVMAX+1).  }                                                
                                                                               
   MAXPT:=21;
                                                                               
{  DELDF := DEFAULT VALUE FOR DELTX(J) }
                                                                               
   DELDF:=0.01;
                                                                               
{  ALPHA := REFLECTION PARAMETER }
                                                                               
   ALPHA:=1.0;
                                                                               
{  BETA := CONTRACTION PARAMETER }
                                                                               
   BETA:=0.5;
                                                                               
{  GAMMA := EXPANSION PARAMETER }
                                                                               
   GAMMA:=2.0;
    
   RZERO:=0.0;
 
{  NO REAL OR DOUBLE PRECISION CONSTANTS ARE USED BEYOND THIS                   
   POINT, IN THIS SUBROUTINE }                                                   

   KFLAG:=0;
   NOREP:=0;
      
   KERFL:=0;
                                                                               
{  CHECK SOME INPUT QUANTITIES, AND SET THEM TO DEFAULT VALUES                  
   IF DESIRED }                                                                  
                                                                               
{  MAKE SURE THAT THE SENSE SWITCH IS OFF }                                      
                                                                               
   JUMP:=2;
      
   DATSW(NSSW,JUMP);                                                    
      
   IF JUMP = 2 Then GOTO 30;                                                    
                                                                               
{  THIS IS THE ONLY USAGE OF THE CONSOLE TYPEWRITER }                            
                                                                               
   WRITELN(' TURN OFF SENSE SWITCH ', NSSW);
   
20: DATSW(NSSW,JUMP);                                                    
      
   IF JUMP = 1 Then GOTO 20;
   
30: JVARY:=0;
      
   IF (NV >= 1) AND (NV <=NVMAX) AND (NVMAX <= MAXPT) Then GOTO 50;
      
   KFLAG:=-1;

   Writeln(fp_out);
   Writeln(fp_out,' ERROR IN INPUT TO SUBROUTINE SIMPLEX, NV=',NV,' NVMAX=',NVMAX,' MAXPT=',MAXPT);
      
   IF NTRACE >= -1 Then Writeln(fp_out,' NACTIV=',NACTIV);                        
      
   Halt;                                                                      

50: For J:=1 to NV do
    begin     
      IF MASK[J] <> 0 Then GOTO 100;                                             
      IF DELMIN^[J] < RZERO Then DELMIN^[J]:=-DELMIN^[J];
{  CHECK THAT DELTX(J) IS NOT NEGLIGIBLE }                                       
      XPLUS:=X^[J]+DELTX^[J];
      IF XPLUS = X^[J] Then GOTO 60;                                             
      XPLUS:=X^[J]-DELTX^[J];
      IF XPLUS <> X^[J] Then GOTO 80;                                             
60:   IF X^[J] = RZERO Then GOTO 70;                                             
      DELTX^[J]:=DELDF*X^[J];
      GOTO 80;                                                               
70:   DELTX^[J]:=DELDF;
80:   IF XMAX^[J] > XMIN^[J] Then GOTO 90;                                        
      XMAX^[J]:=HUGE;
      XMIN^[J]:=-HUGE;
90:   IF X^[J] > XMAX^[J] Then X^[J]:=XMAX^[J];
      IF X^[J] < XMIN^[J] Then X^[J]:=XMIN^[J];
100:end;                                                               

    IF NTRACE < 0 Then GOTO 180;                                                 

    Writeln(fp_out);
    Writeln(fp_out,' SIMPLEX MINIMIZATION,  INITIAL VALUES....');                                             
    Write(fp_out,' MASK = '); For J:=1 to NV do Write(fp_out,' ',MASK[J]);
    Writeln(fp_out);
    Write(fp_out,' X    = '); For J:=1 to NV do Write(fp_out,' ',X^[J]);
    Writeln(fp_out);
    Write(fp_out,' XMAX = '); For J:=1 to NV do Write(fp_out,' ',XMAX^[J]);
    Writeln(fp_out);
    Write(fp_out,' XMIN = '); For J:=1 to NV do Write(fp_out,' ',XMIN^[J]);
    Writeln(fp_out);
    Write(fp_out,' DELTX  = '); For J:=1 to NV do Write(fp_out,' ',DELTX^[J]);
    Writeln(fp_out);
    Write(fp_out,' DELMIN = '); For J:=1 to NV do Write(fp_out,' ',DELMIN^[J]);
    Writeln(fp_out);
                                                                             
{  CALCULATE THE INITIAL P(I) AND Y(I) (Z AND FZ,                               
   RESPECTIVELY).                                                               
                                                                               
   NF = NUMBER OF CALLS TO FUNK SO FAR }
 
180: NF:=0;
                                                                               
{  NACT := NUMBER OF ACTIVE X(J) }
                                                                               
    NACTIV:=0;
     
    For J:=1 to NV do
    begin
     IF MASK[J] <> 0 Then GOTO 200;                                             
     NACTIV:=NACTIV+1;
     For K:=1 to NV do Z^[K,NACTIV]:=X^[K];
     Z^[J,NACTIV]:=Z^[J,NACTIV]+DELTX^[J];
     XS:=X^[J];
     X^[J]:=Z^[J,NACTIV];
     SIFUN;                                                      
     X^[J]:=XS;
     FZ^[NACTIV]:=FOBJ;
200: end;                                                             
      
    IF NTRACE >= 0 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' ',NV,' VARIABLES,  ',NACTIV,' ACTIVE.');
      Writeln(fp_out,' NFMAX=',NFMAX,'  NFLAT=',NFLAT,'  NXTRA=',NXTRA);
      Writeln(fp_out,' ALPHA=',ALPHA:5:2,'  BETA=',BETA:5:2,'  GAMMA=',GAMMA:5:2)
    end;
       
    IF NACTIV > 0 Then GOTO 230;                                                 
      
    KFLAG:=-1;
      
    IF NTRACE >= -1 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' WARNING ...  MASK(J).EQ.0 FOR ALL J IN PROCEDURE SMPLX.');
      Writeln(fp_out,' FOBJ WILL BE EVALUATED BUT NOT MINIMIZED.')
    end;

{  SET THE BASE POINT }                                                          
  
230: NVPLUS:=NACTIV+1;
      
    For K:=1 to NV do  Z^[K,NVPLUS]:=X^[K];
      
    SIFUN;                                                         
      
    FZ^[NVPLUS]:=FOBJ;
                                                                               
{  CHECK THE REPRODUCIBILITY OF FOBJ, AND PRINT THE REMAINING LINES }
                                                                               
    FSAVE:=FOBJ;
      
    SIFUN;                                                         
      
    IF FOBJ = FSAVE Then GOTO 260;                                               
      
    NOREP:=1;
      
    IF NTRACE >= -1 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' WARNING....  FOBJ IS NOT A REPRODUCIBLE FUNCTION OF X(J).');
      Writeln(fp_out,' NF=',NF,' FSAVE=',FSAVE:24:16,' FOBJ=',FOBJ:24:16)
    end;   
  
260:IF NTRACE >= 0 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' FOBJ=',FOBJ:24:16);
      Writeln(fp_out);
      Writeln(fp_out,' BEGIN MINIMIZATION....');
    end;
      
End; {SIBEG}                                                                       


Procedure SIFUN;                                                   
{-----------------------------------------------------
!  SIFUN 1.2          DECEMBER 1991                                             
!                                                                               
!  A.N.S.I. STANDARD FORTRAN 77                                                 
!                                                                               
!  COPYRIGHT (C) 1975, 1990 BY J. P. CHANDLER,                                  
!     COMPUTER SCIENCE DEPARTMENT,                                              
!     OKLAHOMA STATE UNIVERSITY                                                 
!                                                                               
!  SIFUN CALLS FUNK TO COMPUTE FOBJ, IF X IS IN BOUNDS.                         
!  SIFUN IS CALLED BY SIMPLEX AND SIBEG.                                        
!-----------------------------------------------------}
Label 20, 30;
Var     
    JF: Integer;
Begin
    For JF:=1 to NV do
    begin     
      IF ((MASK[JF] = 0) AND (X^[JF] > XMAX^[JF])) OR (X^[JF] < XMIN^[JF]) Then                        
        GOTO 20                                                            
    end;
      
    FUNK;                                                                 
      
    NF:=NF+1;
      
    GOTO 30;                                                                  
   
20: FOBJ:=HUGE;
   
30: IF NTRACE >= 2 Then
    begin
      Writeln(fp_out);
      Writeln(fp_out,' TRIAL FOBJ=',FOBJ:24:16);
      Write(fp_out,' X(J) = ');
      For JF:=1 to NV do Write(fp_out,' ',X^[JF]:24:16);
      Writeln(fp_out)
    end

End;                                                                       

Procedure DATSW(NSSW:Integer; Var JUMP:Integer);                                              
{-----------------------------------------------------------
!  DUMMY VERSION OF SUBROUTINE DATSW  --  ALL SWITCHES OFF.                     
!                                                                               
!  J. P. CHANDLER, COMPUTER SCIENCE DEPARTMENT,                                 
!     OKLAHOMA STATE UNIVERSITY                                                 
!----------------------------------------------------------}
Begin     
  JUMP:=2
End;                                                                       

Procedure STSET;                                                          
{------------------------------------------------------------
!  STSET 3.2          DECEMBER 1991                                             
!                                                                               
!  STSET SETS SOME INPUT QUANTITIES TO DEFAULT VALUES, FOR                      
!  SUBROUTINE SIMPLEX.                      
!                                                                               
!  J. P. CHANDLER, DEPARTMENT OF COMPUTER SCIENCE,                              
!     OKLAHOMA STATE UNIVERSITY                                                 
!                                                                               
!  USAGE.....                                                                   
!                                                                               
!  CALL STSET.                                                                  
!  THEN SET SOME INPUT QUANTITIES (NV, AT LEAST) AND RESET ANY                  
!  OF THOSE SET IN STSET (BETTER VALUES OF X(J), ETC.) BEFORE                   
!  CALLING SIMPLEX.
!------------------------------------------------------------}
Var
    RZERO: Double;
    JX,NVMAX: Integer; 
Begin

    HUGE:=1E30;
      
    RZERO:=0.0;

{  NVMAX IS THE MAXIMUM PERMISSIBLE VALUE OF NV, GIVEN THE                      
   PRESENT DIMENSIONS OF ARRAYS.                                                
   NVMAX IS THE DIMENSION OF THE ARRAYS X(*), XMAX(*), XMIN(*),                 
   DELTX(*), DELMIN(*), AND MASK(*).  NVMAX IS ALSO THE FIRST                   
   DIMENSION OF ERR(*,*).  THE SECOND DIMENSION OF ERR(*,*)                     
   IS NVMAX+1 }                                                                  
       
    NVMAX:=20;

{  THE USER MUST SET NV AFTER CALLING STSET }                                    
      
    NV:=-1;
      
    NTRACE:=0;     {=1 for additional output}
      
    NFMAX:=32000;
     
    MAXIT:=50;
      
    MAXSUB:=30;
      
    METHD:=1;
      
    KALCP:=0;
      
    LEQU:=0;
      
    NFLAT:=1;
      
    MATRX:=105;
      
    NXTRA:=0;
      
    FLAMBD:=1.0;
      
    FNU:=10.0;
      
    KORDIF:=1;
      
    RELDIF:=1E-8;
      
    RELMIN:=1E-6;
                                                                               
    For JX:=1 to NVMAX do
    begin     
      X^[JX]:=RZERO;
      XMAX^[JX]:=HUGE;
      XMIN^[JX]:=-HUGE;
      DELTX^[JX]:=RZERO;
      DELMIN^[JX]:=RZERO;
      MASK[JX]:=0
    end;

End;


{main program}
BEGIN

   {Allocate memory}
   New(X); New(XMAX); New(XMIN); New(DELTX);  New(DELMIN);
   New(Y); New(YSIG); New(T); New(ERGSAV); New(Z); New(FZ);

   {Open input/output files}
   Assign(fp_in,'SIMPTEST.DAT'); Reset(fp_in);
   Assign(fp_out,'SIMPOUT.LST'); Rewrite(fp_out);
           
   Writeln(fp_out,' **  SIMPLEX PROGRAM  **');
                                                                           
{  READ IN THE VALUE OF NPTS, THEN READ IN THE DATA POINTS }                    
     
   Readln(fp_in, NPTS);                                                           
   Writeln(fp_out,' NPTS = ',NPTS);   

   For K:=1 to NPTS do
   begin      
     Readln(fp_in,T^[K,1],T^[K,2],Y^[K],YSIG^[K]);
     Writeln(fp_out,' ',K:2,'  ',T^[K,1]:10:3,T^[K,2]:10:3,Y^[K]:10:3,YSIG^[K]:10:3)
   end;                                                               

   Close(fp_in);

{  INITIALIZE FOR THE FIT }                                                      
     
   STSET;
       
   NV:=3;   {Number of constants X(J) }
      
   NTRACE:=0;   {Normal output}
{  NTRACE:=1;   For additional output}
      
   NFMAX:=550;  {Maximum number of FUNK evaluations}
      
   X^[1]:=10.0;    {Initial guess}
   X^[2]:=1.0;
   X^[3]:=1.0;

{  FIT THE MODEL TO THE DATA USING SIMPLEX METHOD}
     
   SMPLX;                                                         
      
{  USE PROCEDURE FIDO TO COMPUTE ERRORS USING SIMPLEX }                               
      
   ERFRAC:=0.1;
      
   For JXFID:=1 to NV do
   begin

     ERGSAV^[JXFID]:=ERFRAC*ABS(X^[JXFID]);
         
     IF ERGSAV^[JXFID] = 0.0 Then ERGSAV^[JXFID]:=ERFRAC
   
   end;                                                               
     
   STDEVS:=1.0;
      
   TOLFID:=0.05;
      
   MAXIND:=2;
      
   NTRACF:=-1;   {0 or 1 to additional print in FIDO}
         
   NTRSET:=-1;

   For JXFID:=1 to NV do
   begin      

     ERGUES:=ERGSAV^[JXFID];
         
     FIDO(JXFID,STDEVS,TOLFID,ERGUES,MAXIND,NTRSET,NTRACF,FIDINT)                            
  
   end;                                                               
     
   Close(fp_out);
   Writeln;
   Writeln('  Results in file Simpout.lst.');
   Writeln('  Program terminated.');
   Writeln;

   {Free memory}
   Dispose(X); Dispose(XMAX); Dispose(XMIN); Dispose(DELTX);  Dispose(DELMIN);
   Dispose(Y); Dispose(YSIG); Dispose(T); Dispose(ERGSAV); Dispose(Z); Dispose(FZ);

   ReadKey;
   DoneWinCrt                                                            

END.

{SAMPLE DATA (File Simptest.dat must contain:)                                                                          
12                                                                           
0.0       1.0       16.6      0.2                                               
1.0       0.0       12.4      0.2                                               
1.0       1.0       8.2       0.2                                               
0.0       2.0       10.1      0.2                                               
2.0       0.0       7.3       0.2                                               
2.0       2.0       4.7       0.2                                               
1.0       3.0       5.1       0.2                                               
3.0       1.0       4.3       0.2                                               
3.0       3.0       3.0       0.2                                               
5.0       1.0       2.5       0.2                                               
1.0       5.0       3.4       0.2                                               
5.0       5.0       2.1       0.2

end of file tsmplx.pas }