{***************************************************************
* THE SALESMAN'S PROBLEM: FIND THE BEST ITINERARY TO VISIT NV  *
* TOWNS WITH THE MINIMAL TOTAL DISTANCE USING THE "SIMULATED   *
* ANNEALING METHOD". THE MAXIMUM NUMBER OF TOWNS IS HERE SET   *
* TO FIFTY.                                                    *
* ------------------------------------------------------------ *
* LISTE OF MAIN VARIABLES:                                     *
*   NV          : NUMBER OF TOWNS TO VISIT                     *
*   CH$(NV,4)   : TABLE OF TOWN NAMES                          *
*   V$(NV,4)    : SHORTEST ITINERARY                           *
*   D(NV,NV)    : MATRIX OF DISTANCES IN KM (OPTION 2)         *
*   ICO(NV,2)   : GEOGRAPHICAL COORDINATES OF TOWNS (OPTION 1) *
*   IT1 ET LL   : INDICES OF STEPS OF RESOLUTION               *
* ------------------------------------------------------------ *
* SAMPLE RUN:                                                  *
* (Find shortest itinerary to visit 38 French towns:           *
*  AMIENS, ANGOULEME, AUXERRE, BAYONNE, BORDEAUX, BOURGES,     *
*  BREST, CAEN, CALAIS, CHERBOURG, CLERMONT-FERRAND, DIJON,    *
*  GRENOBLE, LE HAVRE, LE MANS, LILLE, LIMOGES, LYON,          *
*  MARSEILLE, METZ, MONTPELLIER, MULHOUSE, NANCY, NANTES, NICE,*
*  ORLEANS, PARIS, PAU, PERIGUEUX, POITIERS, REIMS, RENNES,    *
*  ROUEN, SAINT-ETIENNE, STRASBOURG, TOULOUSE, TOURS, TROYES). *
*                                                              *
* Input file "Villes.dat" contains the number of towns, the    *
* names and geographical coordinates in km of the 38 towns (The*
* origin is an arbitrary point 30 km west of Brest and at the  *
* latitude of Calais):                                         * 
*                                                              *
* AMIENS           525 115                                     *
* ANGOULEME        365 585                                     *
* AUXERRE          620 340                                     *
* BAYONNE          228 830                                     *
* BORDEAUX         305 675                                     *
* BOURGES          535 425                                     *
* BREST            030 260                                     *
* CAEN             333 190                                     *
* CALAIS           497 0                                       *
* CHERBOURG        250 140                                     *
* CLERMONT-FERRAND 590 565                                     *
* DIJON            730 395                                     *
* GRENOBLE         795 630                                     *
* LE_HAVRE         370 157                                     *
* LE_MANS          370 325                                     *
* LILLE            578 034                                     *
* LIMOGES          450 565                                     *
* LYON             725 570                                     *
* MARSEILLE        775 835                                     *
* METZ             810 195                                     *
* MONTPELLIER      655 805                                     *
* MULHOUSE         910 340                                     *
* NANCY            810 245                                     *
* NANTES           235 405                                     *
* NICE             920 780                                     *
* ORLEANS          495 330                                     *
* PARIS            530 228                                     *
* PAU              316 844                                     *
* PERIGUEUX        410 630                                     *
* POITIERS         380 480                                     *
* REIMS            655 185                                     *
* RENNES           233 305                                     *
* ROUEN            440 165                                     *
* SAINT-ETIENNE    690 605                                     *
* STRASBOURG       925 245                                     *
* TOULOUSE         460 810                                     *
* TOURS            410 390                                     *
* TROYES           657 286                                     *
*                                                              *
* Distances between towns (km):                                *
*     1 : as the crow flies                                    *
*     2 : by road                                              *
* Your choice (1 or 2): 2                                      *
*                                                              *
* Starting number of iterations: 10                            *
* Maximum number of iterations : 2000                          *
* Increment value of iterations: 10                            *
*                                                              *
* See results in file Villes.lst.                              *
*                                                              *
* NOTE: Since this program depends on random numbers, you are  *
*       never sure to really have the shortest itinerary in    *
*       one pass. With above data, the output file contains:   *
*                                                              *
*    Shortest Itinerary:                                       *
*                                                              *
*  DIJON            -->  TROYES           -->                  *
*  AUXERRE          -->  PARIS            -->                  *
*  ORLEANS          -->  BOURGES          -->                  *
*  CLERMONT-FERRAND -->  SAINT-ETIENNE    -->                  *
*  LYON             -->  GRENOBLE         -->                  *
*  NICE             -->  MARSEILLE        -->                  *
*  MONTPELLIER      -->  TOULOUSE         -->                  *
*  PAU              -->  BAYONNE          -->                  *
*  BORDEAUX         -->  PERIGUEUX        -->                  *
*  ANGOULEME        -->  LIMOGES          -->                  *
*  POITIERS         -->  TOURS            -->                  *
*  LE_MANS          -->  NANTES           -->                  *
*  BREST            -->  RENNES           -->                  *
*  CHERBOURG        -->  CAEN             -->                  *
*  LE_HAVRE         -->  ROUEN            -->                  *
*  AMIENS           -->  CALAIS           -->                  *
*  LILLE            -->  REIMS            -->                  *
*  METZ             -->  NANCY            -->                  *
*  STRASBOURG       -->  MULHOUSE         -->                  *
*                                                              *
*    NITER= 490   DMIN=  5509.00 KM.                           *
*                                                              *
*                                                              *
* This is probably one of the best itineraries possible, but   *
* other solutions may exist.                                   *
* Note: with option 1 (by air), total distance is of course    *
*       shorter.                                               *
* ------------------------------------------------------------ *
*                                                              *
*                            TPW Release By J-P Moreau, Paris. *
*                                    (www.jpmoreau.fr)         *
****************************************************************
NOTE: With F90:  NBITER=1085  DMIN=5599  (see tanneal.f90)
      With C++:  NBITER=1450  DMIN=5630  (see tanneal.cpp)
---------------------------------------------------------------}
Program TANNEAL;

Uses WinCrt;

Const NBRVILLES=50;     {Maximum number of towns}

Type  NomVille = String[16];
      pDist = ^Dist;
      Dist = Array[1..NBRVILLES,1..NBRVILLES] of Double;
      pIDist = ^IDist;
      IDist = Array[1..NBRVILLES,1..NBRVILLES] of Integer;
      pCoord = ^Coord;
      Coord = Array[1..NBRVILLES] of Double;
      pTab=^Tab;
      Tab = Array[1..NBRVILLES] of NomVille;
      pICord = ^ICord;
      ICord = Array[1..NBRVILLES,1..2] of Integer;
      pItab = ^ITab;
      ITab = Array[1..NBRVILLES] of Integer;

Var
      CH, Vs: pTab;
      ICO: pICord;
      D: pDist;
      L: pITab;
      C: pCoord;
      ID:pIDist;

      i,j, NV: Word;
      choice,diter,itermin,nbreiter,nmax: Integer;
      DDMIN,DMAX,DMIN,DT,TP,X,Y: Double;

      fp, fp1: TEXT;  {I/O Text Files}

      {indices of towns}
      Amiens,Angouleme,Auxerre,Bayonne,Bordeaux,Bourges,
      Brest,Caen,Calais,Cherbourg,Clermont,Dijon,Grenoble,
      LeHavre,LeMans,Lille,Limoges,Lyon,Marseille,Metz,
      Montpellier,Mulhouse,Nancy,Nantes,Nice,
      Orleans,Paris,Pau,Perigueux,Poitiers,Reims,Rennes,
      Rouen,StEtienne,Strasbourg,Toulouse,Tours,Troyes:Integer;


{   GENERATE A RANDOM TIME }
    Procedure TIME;
    Label a15,a20,a25;
    Var i5,it,j5,k5,n5:integer;
	r: Double;
    Begin
      for i5:=1 to NV do
      begin
	L^[i5]:=0;
	C^[i5]:=0.0;
      end;
      n5:=NV;
a15:  it:= 1 + Random(n5);
      k5:=0; j5:=1;
      Repeat
        if C^[j5]=0.0 then Inc(k5);
        if k5<>it then goto a20;
        k5:=j5;
        goto a25;
a20:    Inc(j5)
      Until j5 > NV;
a25:  L^[n5]:=k5;
      Dec(n5);
      C^[k5]:=1.0;
      if n5>0 then goto a15
    End; {TIME}

{ Calculate total distance, DT }
    Function DELTAT: Double;
    Var	T: Double; i6: integer;
    Begin
      T:=0.0;
      for i6:=1 to NV-1 do
        T := T + D^[L^[i6],L^[i6+1]];
      DELTAT:=T
    End;

{   Debug only }
    Procedure aff_L(j:integer);
    Var i,n: integer;
    Begin
      writeln(fp1);
      writeln(fp1,'(',j,')'); n:=0;
      for i:=1 to NV do
      begin
        Inc(n);
        write(fp1,L^[i]:3);
        if n=20 then begin writeln(fp1); n:=0 end
      end
    End;

{   Simulated Annealing Method }
    Procedure ANNEAL(nbiter:Integer; Var itmin:Integer);
    Label a,b,a25;
    Var
      DD,EX,R,R1: Double;
      I,I1,I2,I9,II,IR,IT1,J,J1,J2,K,K2,KF,LL,M,M1:Integer;
      V1s:pTab;
    Begin
      New(V1s);
      DT:=DELTAT;
      DMIN:=35000.0;
      TP:=110.0;
      for IT1:=1 to 40 do
      begin
	for LL:=1 to nbiter do
        begin
	  I:= 1 + Random(NV);
          if I<1 then I:=1;
	  Repeat
	    J:= 1 + Random(NV);
            if J<1 then J:=1;
	  Until J<>I;
	  if J < I then
          begin
	    II:=I;
	    I:=J;
	    J:=II
	  end;
	  I1:=I-1;
	  if I1<1 then I1:=I1+NV;
	  I2:=I+1;
	  J1:=J-1;
	  R:=random;
	  if R > 0.5 then goto a;
	  J2:=J+1;
	  if J2>NV then J2:=J2-NV;
a25:      K:= 1 + Random(NV);
	  K2:=K+1;
	  if K2 > NV then K2:=1;
	  if K=I then goto a25;
	  if K=J then goto a25;
	  if K<=I then
          begin
	    II:=K;
	    K:=I;
	    I:=II;
	    II:=K2;
	    K2:=I2;
	    I2:=II
	  end;
	  if K<=J then
          begin
	    II:=K;
	    K:=J;
	    J:=II;
	    II:=K2;
	    K2:=J2;
	    J2:=II
	  end;
	  DD:=D^[L^[I],L^[J2]]+D^[L^[K],L^[I2]]+D^[L^[J],L^[K2]];
	  DD:=DD-(D^[L^[I],L^[I2]]+D^[L^[J],L^[J2]]+D^[L^[K],L^[K2]]);
	  R1:=Random;
          {Guard against underflow or over flow in Exp}
          if DD/TP>60 then
            EX:=0.0
          else if DD/TP<-60 then
            EX:=1.0
          else
            EX:=Exp(-DD/TP);
	  if R1>EX then goto b;
	  for M:=1 to I do  C^[M]:=1.0*L^[M];
	  M:=I+1;
	  for M1:=J2 to K do
          begin
	    C^[M]:=1.0*L^[M1];
	    Inc(M)
	  end;
	  for M1:=I2 to J do
          begin
	    C^[M]:=1.0*L^[M1];
	    Inc(M)
	  end;
	  if K2 <> 1 then
            for M1:=K2 to NV do
            begin
	      C^[M]:=1.0*L^[M1];
	      Inc(M)
	    end;
	  for M1:=1 to NV do L^[M1] := Round(C^[M1]);
{         aff_L(1); }
	  goto b;
a:        if I=J1 then goto b;
	  DD:=D^[L^[I1],L^[J1]]+D^[L^[I],L^[J]]-(D^[L^[I1],L^[I]]+D^[L^[J1],L^[J]]);
	  R1:=Random;
          {Guard against underflow or over flow in Exp}
          if DD/TP>60 then
            EX:=0.0
          else if DD/TP<-60 then
            EX:=1.0
          else
            EX:=Exp(-DD/TP);
	  if R1<=EX then
          begin
{  	    aff_L(2); }
	    KF:= I+Round(0.5+((J1-I) Div 2));
	    for K:=I to KF do
            begin
	      II:=L^[K];
	      L^[K]:=L^[J1+I-K];
	      L^[J1+I-K]:=II
	    end;
{	    aff_L(3); }
	  end;
b:      end; {for LL}
	TP:=TP*0.9;   {lower temperature}
	DT:=DELTAT;
	DMAX:=DT;
	if DT<DMIN then
        begin
	  DMIN:=DT;
	  for I9:=1 to NV do V1s^[L^[I9]] := CH^[L^[I9]]
	end
      end; {for IT1}

      if DMIN < DDMIN then
      begin
        itmin:=nbiter;
	DDMIN:=DMIN;
	for I9:=1 to NV do Vs^[I9] := V1s^[I9]
      end;
      Dispose(V1s)
    End; {Anneal}

{ Initialize road distances (according to Michelin France Map) }
    Procedure InitDist1;    {part 1}
    Begin
{     initialize town indices }
      Amiens    :=  1;
      Angouleme :=  2;
      Auxerre   :=  3;
      Bayonne   :=  4;
      Bordeaux  :=  5;
      Bourges   :=  6;
      Brest     :=  7;
      Caen      :=  8;
      Calais    :=  9;
      Cherbourg := 10;
      Clermont  := 11;
      Dijon     := 12;
      Grenoble  := 13;
      LeHavre   := 14;
      LeMans    := 15;
      Lille     := 16;
      Limoges   := 17;
      Lyon      := 18;
      Marseille := 19;
      Metz      := 20;
      Montpellier:=21;
      Mulhouse  := 22;
      Nancy     := 23;
      Nantes    := 24;
      Nice      := 25;
      Orleans   := 26;
      Paris     := 27;
      Pau       := 28;
      Perigueux := 29;
      Poitiers  := 30;
      Reims     := 31;
      Rennes    := 32;
      Rouen     := 33;
      StEtienne := 34;
      Strasbourg:= 35;
      Toulouse  := 36;
      Tours     := 37;
      Troyes    := 38;
{     distances between Amiens and other towns }
      D^[Amiens,Angouleme]  :=588.0;
      D^[Amiens,Auxerre]    :=302.0;
      D^[Amiens,Bayonne]    :=888.0;
      D^[Amiens,Bordeaux]   :=711.0;
      D^[Amiens,Bourges]    :=357.0;
      D^[Amiens,Brest]      :=617.0;
      D^[Amiens,Caen]       :=246.0;
      D^[Amiens,Calais]     :=154.0;
      D^[Amiens,Cherbourg]  :=366.0;
      D^[Amiens,Clermont]   :=520.0;
      D^[Amiens,Dijon]      :=444.0;
      D^[Amiens,Grenoble]   :=696.0;
      D^[Amiens,LeHavre]    :=178.0;
      D^[Amiens,LeMans]     :=317.0;
      D^[Amiens,Lille]      :=110.0;
      D^[Amiens,Limoges]    :=586.0;
      D^[Amiens,Lyon]       :=592.0;
      D^[Amiens,Marseille]  :=907.0;
      D^[Amiens,Metz]       :=342.0;
      D^[Amiens,Montpellier]:=981.0;
      D^[Amiens,Mulhouse]   :=512.0;
      D^[Amiens,Nancy]      :=353.0;
      D^[Amiens,Nantes]     :=475.0;
      D^[Amiens,Nice]       :=1063.0;
      D^[Amiens,Orleans]    :=262.0;
      D^[Amiens,Paris]      :=132.0;
      D^[Amiens,Pau]        :=910.0;
      D^[Amiens,Perigueux]  :=674.0;
      D^[Amiens,Poitiers]   :=478.0;
      D^[Amiens,Reims]      :=154.0;
      D^[Amiens,Rennes]     :=422.0;
      D^[Amiens,Rouen]      :=122.0;
      D^[Amiens,StEtienne]  :=648.0;
      D^[Amiens,Strasbourg] :=505.0;
      D^[Amiens,Toulouse]   :=837.0;
      D^[Amiens,Tours]      :=366.0;
      D^[Amiens,Troyes]     :=284.0;
{     distances between Angouleme and other towns }
      D^[Angouleme,Auxerre]    :=439.0;
      D^[Angouleme,Bayonne]    :=295.0;
      D^[Angouleme,Bordeaux]   :=118.0;
      D^[Angouleme,Bourges]    :=295.0;
      D^[Angouleme,Brest]      :=536.0;
      D^[Angouleme,Caen]       :=460.0;
      D^[Angouleme,Calais]     :=751.0;
      D^[Angouleme,Cherbourg]  :=559.0;
      D^[Angouleme,Clermont]   :=283.0;
      D^[Angouleme,Dijon]      :=530.0;
      D^[Angouleme,Grenoble]   :=569.0;
      D^[Angouleme,LeHavre]    :=532.0;
      D^[Angouleme,LeMans]     :=334.0;
      D^[Angouleme,Lille]      :=678.0;
      D^[Angouleme,Limoges]    :=103.0;
      D^[Angouleme,Lyon]       :=460.0;
      D^[Angouleme,Marseille]  :=740.0;
      D^[Angouleme,Metz]       :=794.0;
      D^[Angouleme,Montpellier]:=640.0;
      D^[Angouleme,Mulhouse]   :=759.0;
      D^[Angouleme,Nancy]      :=692.0;
      D^[Angouleme,Nantes]     :=240.0;
      D^[Angouleme,Nice]       :=927.0;
      D^[Angouleme,Orleans]    :=339.0;
      D^[Angouleme,Paris]      :=461.0;
      D^[Angouleme,Pau]        :=317.0;
      D^[Angouleme,Perigueux]  :=86.0;
      D^[Angouleme,Poitiers]   :=110.0;
      D^[Angouleme,Reims]      :=630.0;
      D^[Angouleme,Rennes]     :=347.0;
      D^[Angouleme,Rouen]      :=504.0;
      D^[Angouleme,StEtienne]  :=397.0;
      D^[Angouleme,Strasbourg] :=865.0;
      D^[Angouleme,Toulouse]   :=334.0;
      D^[Angouleme,Tours]      :=227.0;
      D^[Angouleme,Troyes]     :=510.0;
{     distances between Auxerre and other towns }
      D^[Auxerre,Bayonne]    :=707.0;
      D^[Auxerre,Bordeaux]   :=530.0;
      D^[Auxerre,Bourges]    :=144.0;
      D^[Auxerre,Brest]      :=695.0;
      D^[Auxerre,Caen]       :=407.0;
      D^[Auxerre,Calais]     :=454.0;
      D^[Auxerre,Cherbourg]  :=537.0;
      D^[Auxerre,Clermont]   :=240.0;
      D^[Auxerre,Dijon]      :=150.0;
      D^[Auxerre,Grenoble]   :=445.0;
      D^[Auxerre,LeHavre]    :=366.0;
      D^[Auxerre,LeMans]     :=288.0;
      D^[Auxerre,Lille]      :=501.0;
      D^[Auxerre,Limoges]    :=331.0;
      D^[Auxerre,Lyon]       :=298.0;
      D^[Auxerre,Marseille]  :=614.0;
      D^[Auxerre,Metz]       :=314.0;
      D^[Auxerre,Montpellier]:=596.0;
      D^[Auxerre,Mulhouse]   :=379.0;
      D^[Auxerre,Nancy]      :=258.0;
      D^[Auxerre,Nantes]     :=450.0;
      D^[Auxerre,Nice]       :=768.0;
      D^[Auxerre,Orleans]    :=150.0;
      D^[Auxerre,Paris]      :=162.0;
      D^[Auxerre,Pau]        :=794.0;
      D^[Auxerre,Perigueux]  :=433.0;
      D^[Auxerre,Poitiers]   :=331.0;
      D^[Auxerre,Reims]      :=202.0;
      D^[Auxerre,Rennes]     :=447.0;
      D^[Auxerre,Rouen]      :=301.0;
      D^[Auxerre,StEtienne]  :=361.0;
      D^[Auxerre,Strasbourg] :=410.0;
      D^[Auxerre,Toulouse]   :=621.0;
      D^[Auxerre,Tours]      :=262.0;
      D^[Auxerre,Troyes]     :=81.0;
{     distances between Bayonne and other towns }
      D^[Bayonne,Bordeaux]   :=177.0;
      D^[Bayonne,Bourges]    :=590.0;
      D^[Bayonne,Brest]      :=798.0;
      D^[Bayonne,Caen]       :=755.0;
      D^[Bayonne,Calais]     :=1046.0;
      D^[Bayonne,Cherbourg]  :=854.0;
      D^[Bayonne,Clermont]   :=545.0;
      D^[Bayonne,Dijon]      :=811.0;
      D^[Bayonne,Grenoble]   :=812.0;
      D^[Bayonne,LeHavre]    :=827.0;
      D^[Bayonne,LeMans]     :=602.0;
      D^[Bayonne,Lille]      :=973.0;
      D^[Bayonne,Limoges]    :=398.0;
      D^[Bayonne,Lyon]       :=722.0;
      D^[Bayonne,Marseille]  :=680.0;
      D^[Bayonne,Metz]       :=1079.0;
      D^[Bayonne,Montpellier]:=516.0;
      D^[Bayonne,Mulhouse]   :=1097.0;
      D^[Bayonne,Nancy]      :=965.0;
      D^[Bayonne,Nantes]     :=502.0;
      D^[Bayonne,Nice]       :=835.0;
      D^[Bayonne,Orleans]    :=634.0;
      D^[Bayonne,Paris]      :=756.0;
      D^[Bayonne,Pau]        :=103.0;
      D^[Bayonne,Perigueux]  :=298.0;
      D^[Bayonne,Poitiers]   :=405.0;
      D^[Bayonne,Reims]      :=891.0;
      D^[Bayonne,Rennes]     :=609.0;
      D^[Bayonne,Rouen]      :=799.0;
      D^[Bayonne,StEtienne]  :=691.0;
      D^[Bayonne,Strasbourg] :=1237.0;
      D^[Bayonne,Toulouse]   :=276.0;
      D^[Bayonne,Tours]      :=522.0;
      D^[Bayonne,Troyes]     :=805.0;
{     distances between Bordeaux and other towns }
      D^[Bordeaux,Bourges]    :=413.0;
      D^[Bordeaux,Brest]      :=621.0;
      D^[Bordeaux,Caen]       :=578.0;
      D^[Bordeaux,Calais]     :=869.0;
      D^[Bordeaux,Cherbourg]  :=677.0;
      D^[Bordeaux,Clermont]   :=368.0;
      D^[Bordeaux,Dijon]      :=634.0;
      D^[Bordeaux,Grenoble]   :=654.0;
      D^[Bordeaux,LeHavre]    :=650.0;
      D^[Bordeaux,LeMans]     :=425.0;
      D^[Bordeaux,Lille]      :=796.0;
      D^[Bordeaux,Limoges]    :=221.0;
      D^[Bordeaux,Lyon]       :=545.0;
      D^[Bordeaux,Marseille]  :=648.0;
      D^[Bordeaux,Metz]       :=902.0;
      D^[Bordeaux,Montpellier]:=484.0;
      D^[Bordeaux,Mulhouse]   :=816.0;
      D^[Bordeaux,Nancy]      :=788.0;
      D^[Bordeaux,Nantes]     :=325.0;
      D^[Bordeaux,Nice]       :=803.0;
      D^[Bordeaux,Orleans]    :=457.0;
      D^[Bordeaux,Paris]      :=579.0;
      D^[Bordeaux,Pau]        :=199.0;
      D^[Bordeaux,Perigueux]  :=121.0;
      D^[Bordeaux,Poitiers]   :=228.0;
      D^[Bordeaux,Reims]      :=714.0;
      D^[Bordeaux,Rennes]     :=432.0;
      D^[Bordeaux,Rouen]      :=622.0;
      D^[Bordeaux,StEtienne]  :=514.0;
      D^[Bordeaux,Strasbourg] :=1060.0;
      D^[Bordeaux,Toulouse]   :=244.0;
      D^[Bordeaux,Tours]      :=345.0;
      D^[Bordeaux,Troyes]     :=628.0;
{     distances between Bourges and other towns }  
      D^[Bourges,Brest]      :=634.0;
      D^[Bourges,Caen]       :=363.0;
      D^[Bourges,Calais]     :=526.0;
      D^[Bourges,Cherbourg]  :=493.0;
      D^[Bourges,Clermont]   :=180.0;
      D^[Bourges,Dijon]      :=259.0;
      D^[Bourges,Grenoble]   :=466.0;
      D^[Bourges,LeHavre]    :=508.0;
      D^[Bourges,LeMans]     :=230.0;
      D^[Bourges,Lille]      :=455.0;
      D^[Bourges,Limoges]    :=187.0;
      D^[Bourges,Lyon]       :=357.0;
      D^[Bourges,Marseille]  :=637.0;
      D^[Bourges,Metz]       :=458.0;
      D^[Bourges,Montpellier]:=540.0;
      D^[Bourges,Mulhouse]   :=488.0;
      D^[Bourges,Nancy]      :=402.0;
      D^[Bourges,Nantes]     :=347.0;
      D^[Bourges,Nice]       :=792.0;
      D^[Bourges,Orleans]    :=106.0;
      D^[Bourges,Paris]      :=236.0;
      D^[Bourges,Pau]        :=734.0;
      D^[Bourges,Perigueux]  :=289.0;
      D^[Bourges,Poitiers]   :=187.0;
      D^[Bourges,Reims]      :=346.0;
      D^[Bourges,Rennes]     :=506.0;
      D^[Bourges,Rouen]      :=323.0;
      D^[Bourges,StEtienne]  :=326.0;
      D^[Bourges,Strasbourg] :=554.0;
      D^[Bourges,Toulouse]   :=561.0;
      D^[Bourges,Tours]      :=150.0;
      D^[Bourges,Troyes]     :=225.0;
{     distances between Brest and other towns }
      D^[Brest,Caen]       :=371.0;
      D^[Brest,Calais]     :=713.0;
      D^[Brest,Cherbourg]  :=415.0;
      D^[Brest,Clermont]   :=750.0;
      D^[Brest,Dijon]      :=854.0;
      D^[Brest,Grenoble]   :=1106.0;
      D^[Brest,LeHavre]    :=479.0;
      D^[Brest,LeMans]     :=402.0;
      D^[Brest,Lille]      :=814.0;
      D^[Brest,Limoges]    :=606.0;
      D^[Brest,Lyon]       :=1002.0;
      D^[Brest,Marseille]  :=1317.0;
      D^[Brest,Metz]       :=920.0;
      D^[Brest,Montpellier]:=1105.0;
      D^[Brest,Mulhouse]   :=1081.0;
      D^[Brest,Nancy]      :=915.0;
      D^[Brest,Nantes]     :=296.0;
      D^[Brest,Nice]       :=1472.0;
      D^[Brest,Orleans]    :=545.0;
      D^[Brest,Paris]      :=596.0;
      D^[Brest,Pau]        :=820.0;
      D^[Brest,Perigueux]  :=622.0;
      D^[Brest,Poitiers]   :=485.0;
      D^[Brest,Reims]      :=731.0;
      D^[Brest,Rennes]     :=245.0;
      D^[Brest,Rouen]      :=495.0;
      D^[Brest,StEtienne]  :=896.0;
      D^[Brest,Strasbourg] :=1077.0;
      D^[Brest,Toulouse]   :=865.0;
      D^[Brest,Tours]      :=484.0;
      D^[Brest,Troyes]     :=733.0;
{     distances between Caen and other towns }
      D^[Caen,Calais]     :=342.0;
      D^[Caen,Cherbourg]  :=120.0;
      D^[Caen,Clermont]   :=617.0;
      D^[Caen,Dijon]      :=542.0;
      D^[Caen,Grenoble]   :=794.0;
      D^[Caen,LeHavre]    :=108.0;
      D^[Caen,LeMans]     :=156.0;
      D^[Caen,Lille]      :=342.0;
      D^[Caen,Limoges]    :=469.0;
      D^[Caen,Lyon]       :=690.0;
      D^[Caen,Marseille]  :=1005.0;
      D^[Caen,Metz]       :=567.0;
      D^[Caen,Montpellier]:=988.0;
      D^[Caen,Mulhouse]   :=769.0;
      D^[Caen,Nancy]      :=540.0;
      D^[Caen,Nantes]     :=283.0;
      D^[Caen,Nice]       :=1160.0;
      D^[Caen,Orleans]    :=257.0;
      D^[Caen,Paris]      :=240.0;
      D^[Caen,Pau]        :=777.0;
      D^[Caen,Perigueux]  :=606.0;
      D^[Caen,Poitiers]   :=410.0;
      D^[Caen,Reims]      :=379.0;
      D^[Caen,Rennes]     :=176.0;
      D^[Caen,Rouen]      :=124.0;
      D^[Caen,StEtienne]  :=745.0;
      D^[Caen,Strasbourg] :=725.0;
      D^[Caen,Toulouse]   :=768.0;
      D^[Caen,Tours]      :=233.0;
      D^[Caen,Troyes]     :=418.0;
{     distances between Calais and other towns }
      D^[Calais,Cherbourg]  :=462.0;
      D^[Calais,Clermont]   :=678.0;
      D^[Calais,Dijon]      :=603.0;
      D^[Calais,Grenoble]   :=855.0;
      D^[Calais,LeHavre]    :=284.0;
      D^[Calais,LeMans]     :=418.0;
      D^[Calais,Lille]      :=112.0;
      D^[Calais,Limoges]    :=746.0;
      D^[Calais,Lyon]       :=751.0;
      D^[Calais,Marseille]  :=1066.0;
      D^[Calais,Metz]       :=467.0;
      D^[Calais,Montpellier]:=1049.0;
      D^[Calais,Mulhouse]   :=636.0;
      D^[Calais,Nancy]      :=523.0;
      D^[Calais,Nantes]     :=668.0;
      D^[Calais,Nice]       :=1221.0;
      D^[Calais,Orleans]    :=420.0;
      D^[Calais,Paris]      :=292.0;
      D^[Calais,Pau]        :=1059.0;
      D^[Calais,Perigueux]  :=828.0;
      D^[Calais,Poitiers]   :=630.0;
      D^[Calais,Reims]      :=278.0;
      D^[Calais,Rennes]     :=518.0;
      D^[Calais,Rouen]      :=218.0;
      D^[Calais,StEtienne]  :=806.0;
      D^[Calais,Strasbourg] :=624.0;
      D^[Calais,Toulouse]   :=995.0;
      D^[Calais,Tours]      :=524.0;
      D^[Calais,Troyes]     :=399.0;
{     distances between Cherbourg and other towns }
      D^[Cherbourg,Clermont]   :=663.0;
      D^[Cherbourg,Dijon]      :=684.0;
      D^[Cherbourg,Grenoble]   :=936.0;
      D^[Cherbourg,LeHavre]    :=228.0;
      D^[Cherbourg,LeMans]     :=276.0;
      D^[Cherbourg,Lille]      :=462.0;
      D^[Cherbourg,Limoges]    :=589.0;
      D^[Cherbourg,Lyon]       :=810.0;
      D^[Cherbourg,Marseille]  :=1125.0;
      D^[Cherbourg,Metz]       :=687.0;
      D^[Cherbourg,Montpellier]:=1108.0;
      D^[Cherbourg,Mulhouse]   :=889.0;
      D^[Cherbourg,Nancy]      :=660.0;
      D^[Cherbourg,Nantes]     :=319.0;
      D^[Cherbourg,Nice]       :=1280.0;
      D^[Cherbourg,Orleans]    :=507.0;
      D^[Cherbourg,Paris]      :=360.0;
      D^[Cherbourg,Pau]        :=876.0;
      D^[Cherbourg,Perigueux]  :=645.0;
      D^[Cherbourg,Poitiers]   :=450.0;
      D^[Cherbourg,Reims]      :=499.0;
      D^[Cherbourg,Rennes]     :=212.0;
      D^[Cherbourg,Rouen]      :=244.0;
      D^[Cherbourg,StEtienne]  :=865.0;
      D^[Cherbourg,Strasbourg] :=845.0;
      D^[Cherbourg,Toulouse]   :=893.0;
      D^[Cherbourg,Tours]      :=356.0;
      D^[Cherbourg,Troyes]     :=399.0;
{     distances between Clermont-Ferrand and other towns }
      D^[Clermont,Dijon]      :=282.0;
      D^[Clermont,Grenoble]   :=286.0;
      D^[Clermont,LeHavre]    :=581.0;
      D^[Clermont,LeMans]     :=507.0;
      D^[Clermont,Lille]      :=605.0;
      D^[Clermont,Limoges]    :=178.0;
      D^[Clermont,Lyon]       :=177.0;
      D^[Clermont,Marseille]  :=457.0;
      D^[Clermont,Metz]       :=546.0;
      D^[Clermont,Montpellier]:=360.0;
      D^[Clermont,Mulhouse]   :=463.0;
      D^[Clermont,Nancy]      :=490.0;
      D^[Clermont,Nantes]     :=460.0;
      D^[Clermont,Nice]       :=612.0;
      D^[Clermont,Orleans]    :=289.0;
      D^[Clermont,Paris]      :=508.0;
      D^[Clermont,Pau]        :=554.0;
      D^[Clermont,Perigueux]  :=250.0;
      D^[Clermont,Poitiers]   :=285.0;
      D^[Clermont,Reims]      :=487.0;
      D^[Clermont,Rennes]     :=528.0;
      D^[Clermont,Rouen]      :=516.0;
      D^[Clermont,StEtienne]  :=146.0;
      D^[Clermont,Strasbourg] :=569.0;
      D^[Clermont,Toulouse]   :=501.0;
      D^[Clermont,Tours]      :=307.0;
      D^[Clermont,Troyes]     :=321.0;
{     distances between Dijon and other towns }
      D^[Dijon,Grenoble]   :=295.0;
      D^[Dijon,LeHavre]    :=506.0;
      D^[Dijon,LeMans]     :=435.0;
      D^[Dijon,Lille]      :=530.0;
      D^[Dijon,Limoges]    :=425.0;
      D^[Dijon,Lyon]       :=191.0;
      D^[Dijon,Marseille]  :=507.0;
      D^[Dijon,Metz]       :=264.0;
      D^[Dijon,Montpellier]:=490.0;
      D^[Dijon,Mulhouse]   :=229.0;
      D^[Dijon,Nancy]      :=210.0;
      D^[Dijon,Nantes]     :=597.0;
      D^[Dijon,Nice]       :=662.0;
      D^[Dijon,Orleans]    :=297.0;
      D^[Dijon,Paris]      :=312.0;
      D^[Dijon,Pau]        :=836.0;
      D^[Dijon,Perigueux]  :=520.0;
      D^[Dijon,Poitiers]   :=446.0;
      D^[Dijon,Reims]      :=286.0;
      D^[Dijon,Rennes]     :=605.0;
      D^[Dijon,Rouen]      :=441.0;
      D^[Dijon,StEtienne]  :=247.0;
      D^[Dijon,Strasbourg] :=335.0;
      D^[Dijon,Toulouse]   :=726.0;
      D^[Dijon,Tours]      :=409.0;
      D^[Dijon,Troyes]     :=152.0;
{     distances between Grenoble and other towns }                                                                            
      D^[Grenoble,LeHavre]    :=758.0;
      D^[Grenoble,LeMans]     :=613.0;
      D^[Grenoble,Lille]      :=782.0;
      D^[Grenoble,Limoges]    :=464.0;
      D^[Grenoble,Lyon]       :=104.0;
      D^[Grenoble,Marseille]  :=281.0;
      D^[Grenoble,Metz]       :=560.0;
      D^[Grenoble,Montpellier]:=300.0;
      D^[Grenoble,Mulhouse]   :=450.0;
      D^[Grenoble,Nancy]      :=505.0;
      D^[Grenoble,Nantes]     :=730.0;
      D^[Grenoble,Nice]       :=334.0;
      D^[Grenoble,Orleans]    :=549.0;
      D^[Grenoble,Paris]      :=564.0;
      D^[Grenoble,Pau]        :=709.0;
      D^[Grenoble,Perigueux]  :=544.0;
      D^[Grenoble,Poitiers]   :=585.0;
      D^[Grenoble,Reims]      :=581.0;
      D^[Grenoble,Rennes]     :=857.0;
      D^[Grenoble,Rouen]      :=693.0;
      D^[Grenoble,StEtienne]  :=143.0;
      D^[Grenoble,Strasbourg] :=548.0;
      D^[Grenoble,Toulouse]   :=536.0;
      D^[Grenoble,Tours]      :=533.0;
      D^[Grenoble,Troyes]     :=447.0;
{     distances between Le Havre and other towns }
      D^[LeHavre,LeMans]     :=264.0;
      D^[LeHavre,Lille]      :=284.0;
      D^[LeHavre,Limoges]    :=577.0;
      D^[LeHavre,Lyon]       :=654.0;
      D^[LeHavre,Marseille]  :=969.0;
      D^[LeHavre,Metz]       :=531.0;
      D^[LeHavre,Montpellier]:=952.0;
      D^[LeHavre,Mulhouse]   :=733.0;
      D^[LeHavre,Nancy]      :=504.0;
      D^[LeHavre,Nantes]     :=500.0;
      D^[LeHavre,Nice]       :=1124.0;
      D^[LeHavre,Orleans]    :=282.0;
      D^[LeHavre,Paris]      :=204.0;
      D^[LeHavre,Pau]        :=876.0;
      D^[LeHavre,Perigueux]  :=645.0;
      D^[LeHavre,Poitiers]   :=450.0;
      D^[LeHavre,Reims]      :=343.0;
      D^[LeHavre,Rennes]     :=284.0;
      D^[LeHavre,Rouen]      :=88.0;
      D^[LeHavre,StEtienne]  :=709.0;
      D^[LeHavre,Strasbourg] :=689.0;
      D^[LeHavre,Toulouse]   :=840.0;
      D^[LeHavre,Tours]      :=305.0;
      D^[LeHavre,Troyes]     :=502.0;
{     distances between Le Mans and other towns }  
      D^[LeMans,Lille]      :=418.0;
      D^[LeMans,Limoges]    :=313.0;
      D^[LeMans,Lyon]       :=509.0;
      D^[LeMans,Marseille]  :=816.0;
      D^[LeMans,Metz]       :=534.0;
      D^[LeMans,Montpellier]:=716.0;
      D^[LeMans,Mulhouse]   :=662.0;
      D^[LeMans,Nancy]      :=503.0;
      D^[LeMans,Nantes]     :=181.0;
      D^[LeMans,Nice]       :=971.0;
      D^[LeMans,Orleans]    :=150.0;
      D^[LeMans,Paris]      :=203.0;
      D^[LeMans,Pau]        :=651.0;
      D^[LeMans,Perigueux]  :=420.0;
      D^[LeMans,Poitiers]   :=224.0;
      D^[LeMans,Reims]      :=345.0;
      D^[LeMans,Rennes]     :=157.0;
      D^[LeMans,Rouen]      :=200.0;
      D^[LeMans,StEtienne]  :=504.0;
      D^[LeMans,Strasbourg] :=749.0;
      D^[LeMans,Toulouse]   :=615.0;
      D^[LeMans,Tours]      :=80.0;
      D^[LeMans,Troyes]     :=326.0;
{     distances between Lille and other towns }  
      D^[Lille,Limoges]    :=673.0;
      D^[Lille,Lyon]       :=678.0;
      D^[Lille,Marseille]  :=993.0;
      D^[Lille,Metz]       :=349.0;
      D^[Lille,Montpellier]:=976.0;
      D^[Lille,Mulhouse]   :=563.0;
      D^[Lille,Nancy]      :=405.0;
      D^[Lille,Nantes]     :=595.0;
      D^[Lille,Nice]       :=1148.0;
      D^[Lille,Orleans]    :=347.0;
      D^[Lille,Paris]      :=219.0;
      D^[Lille,Pau]        :=995.0;
      D^[Lille,Perigueux]  :=784.0;
      D^[Lille,Poitiers]   :=588.0;
      D^[Lille,Reims]      :=205.0;
      D^[Lille,Rennes]     :=565.0;
      D^[Lille,Rouen]      :=218.0;
      D^[Lille,StEtienne]  :=733.0;
      D^[Lille,Strasbourg] :=505.0;
      D^[Lille,Toulouse]   :=922.0;
      D^[Lille,Tours]      :=451.0;
      D^[Lille,Troyes]     :=326.0;
{     distances between Limoges and other towns }
      D^[Limoges,Lyon]       :=355.0;
      D^[Limoges,Marseille]  :=635.0;
      D^[Limoges,Metz]       :=645.0;
      D^[Limoges,Montpellier]:=550.0;
      D^[Limoges,Mulhouse]   :=654.0;
      D^[Limoges,Nancy]      :=589.0;
      D^[Limoges,Nantes]     :=310.0;
      D^[Limoges,Nice]       :=790.0;
      D^[Limoges,Orleans]    :=324.0;
      D^[Limoges,Paris]      :=454.0;
      D^[Limoges,Pau]        :=363.0;
      D^[Limoges,Perigueux]  :=102.0;
      D^[Limoges,Poitiers]   :=121.0;
      D^[Limoges,Reims]      :=533.0;
      D^[Limoges,Rennes]     :=500.0;
      D^[Limoges,Rouen]      :=541.0;
      D^[Limoges,StEtienne]  :=324.0;
      D^[Limoges,Strasbourg] :=760.0;
      D^[Limoges,Toulouse]   :=315.0;
      D^[Limoges,Tours]      :=225.0;
      D^[Limoges,Troyes]     :=412.0;
{     distances between Lyon and other towns }
      D^[Lyon,Marseille]  :=315.0;
      D^[Lyon,Metz]       :=456.0;
      D^[Lyon,Montpellier]:=298.0;
      D^[Lyon,Mulhouse]   :=375.0;
      D^[Lyon,Nancy]      :=401.0;
      D^[Lyon,Nantes]     :=626.0;
      D^[Lyon,Nice]       :=470.0;
      D^[Lyon,Orleans]    :=445.0;
      D^[Lyon,Paris]      :=460.0;
      D^[Lyon,Pau]        :=731.0;
      D^[Lyon,Perigueux]  :=415.0;
      D^[Lyon,Poitiers]   :=462.0;
      D^[Lyon,Reims]      :=477.0;
      D^[Lyon,Rennes]     :=753.0;
      D^[Lyon,Rouen]      :=589.0;
      D^[Lyon,StEtienne]  :=63.0;
      D^[Lyon,Strasbourg] :=480.0;
      D^[Lyon,Toulouse]   :=534.0;
      D^[Lyon,Tours]      :=429.0;
      D^[Lyon,Troyes]     :=343.0;
{     distances between Marseille and other towns }
      D^[Marseille,Metz]       :=771.0;
      D^[Marseille,Montpellier]:=168.0;
      D^[Marseille,Mulhouse]   :=690.0;
      D^[Marseille,Nancy]      :=716.0;
      D^[Marseille,Nantes]     :=973.0;
      D^[Marseille,Nice]       :=187.0;
      D^[Marseille,Orleans]    :=760.0;
      D^[Marseille,Paris]      :=776.0;
      D^[Marseille,Pau]        :=577.0;
      D^[Marseille,Perigueux]  :=654.0;
      D^[Marseille,Poitiers]   :=742.0;
      D^[Marseille,Reims]      :=793.0;
      D^[Marseille,Rennes]     :=1069.0;
      D^[Marseille,Rouen]      :=904.0;
      D^[Marseille,StEtienne]  :=308.0;
      D^[Marseille,Strasbourg] :=796.0;
      D^[Marseille,Toulouse]   :=404.0;
      D^[Marseille,Tours]      :=736.0;
      D^[Marseille,Troyes]     :=658.0;
{     distances between Metz and other towns }  
      D^[Metz,Montpellier]:=754.0;
      D^[Metz,Mulhouse]   :=227.0;
      D^[Metz,Nancy]      :=56.0;
      D^[Metz,Nantes]     :=701.0;
      D^[Metz,Nice]       :=926.0;
      D^[Metz,Orleans]    :=453.0;
      D^[Metz,Paris]      :=331.0;
      D^[Metz,Pau]        :=1192.0;
      D^[Metz,Perigueux]  :=747.0;
      D^[Metz,Poitiers]   :=645.0;
      D^[Metz,Reims]      :=188.0;
      D^[Metz,Rennes]     :=671.0;
      D^[Metz,Rouen]      :=419.0;
      D^[Metz,StEtienne]  :=511.0;
      D^[Metz,Strasbourg] :=163.0;
      D^[Metz,Toulouse]   :=990.0;
      D^[Metz,Tours]      :=557.0;
      D^[Metz,Troyes]     :=208.0;
    End;

    Procedure InitDist2;  {part 2}
    Var I,J: Integer;
    Begin
{     distances between Montpellier and other towns }
      D^[Montpellier,Mulhouse]   :=673.0;
      D^[Montpellier,Nancy]      :=700.0;
      D^[Montpellier,Nantes]     :=809.0;
      D^[Montpellier,Nice]       :=323.0;
      D^[Montpellier,Orleans]    :=743.0;
      D^[Montpellier,Paris]      :=759.0;
      D^[Montpellier,Pau]        :=413.0;
      D^[Montpellier,Perigueux]  :=490.0;
      D^[Montpellier,Poitiers]   :=659.0;
      D^[Montpellier,Reims]      :=776.0;
      D^[Montpellier,Rennes]     :=916.0;
      D^[Montpellier,Rouen]      :=887.0;
      D^[Montpellier,StEtienne]  :=291.0;
      D^[Montpellier,Strasbourg] :=779.0;
      D^[Montpellier,Toulouse]   :=240.0;
      D^[Montpellier,Tours]      :=636.0;
      D^[Montpellier,Troyes]     :=642.0;
{     distances between Mulhouse and other towns }
      D^[Mulhouse,Nancy]      :=171.0;
      D^[Mulhouse,Nantes]     :=824.0;
      D^[Mulhouse,Nice]       :=698.0;
      D^[Mulhouse,Orleans]    :=524.0;
      D^[Mulhouse,Paris]      :=539.0;
      D^[Mulhouse,Pau]        :=1106.0;
      D^[Mulhouse,Perigueux]  :=701.0;
      D^[Mulhouse,Poitiers]   :=675.0;
      D^[Mulhouse,Reims]      :=358.0;
      D^[Mulhouse,Rennes]     :=832.0;
      D^[Mulhouse,Rouen]      :=668.0;
      D^[Mulhouse,StEtienne]  :=430.0;
      D^[Mulhouse,Strasbourg] :=113.0;
      D^[Mulhouse,Toulouse]   :=909.0;
      D^[Mulhouse,Tours]      :=636.0;
      D^[Mulhouse,Troyes]     :=274.0;
{     distances between Nancy and other towns }
      D^[Nancy,Nantes]     :=708.0;
      D^[Nancy,Nice]       :=872.0;
      D^[Nancy,Orleans]    :=365.0;
      D^[Nancy,Paris]      :=300.0;
      D^[Nancy,Pau]        :=1046.0;
      D^[Nancy,Perigueux]  :=728.0;
      D^[Nancy,Poitiers]   :=589.0;
      D^[Nancy,Reims]      :=244.0;
      D^[Nancy,Rennes]     :=662.0;
      D^[Nancy,Rouen]      :=439.0;
      D^[Nancy,StEtienne]  :=457.0;
      D^[Nancy,Strasbourg] :=188.0;
      D^[Nancy,Toulouse]   :=871.0;
      D^[Nancy,Tours]      :=477.0;
      D^[Nancy,Troyes]     :=177.0;
{     distances between Nantes and other towns }  
      D^[Nantes,Nice]       :=1128.0;
      D^[Nantes,Orleans]    :=300.0;
      D^[Nantes,Paris]      :=378.0;
      D^[Nantes,Pau]        :=524.0;
      D^[Nantes,Perigueux]  :=326.0;
      D^[Nantes,Poitiers]   :=189.0;
      D^[Nantes,Reims]      :=513.0;
      D^[Nantes,Rennes]     :=107.0;
      D^[Nantes,Rouen]      :=353.0;
      D^[Nantes,StEtienne]  :=606.0;
      D^[Nantes,Strasbourg] :=859.0;
      D^[Nantes,Toulouse]   :=569.0;
      D^[Nantes,Tours]      :=197.0;
      D^[Nantes,Troyes]     :=488.0;
{     distances between Nice and other towns }
      D^[Nice,Orleans]    :=915.0;
      D^[Nice,Paris]      :=931.0;
      D^[Nice,Pau]        :=732.0;
      D^[Nice,Perigueux]  :=809.0;
      D^[Nice,Poitiers]   :=775.0;
      D^[Nice,Reims]      :=948.0;
      D^[Nice,Rennes]     :=1224.0;
      D^[Nice,Rouen]      :=1059.0;
      D^[Nice,StEtienne]  :=463.0;
      D^[Nice,Strasbourg] :=808.0;
      D^[Nice,Toulouse]   :=559.0;
      D^[Nice,Tours]      :=891.0;
      D^[Nice,Troyes]     :=813.0;
{     distances between Orleans and other towns }  
      D^[Orleans,Paris]      :=130.0;
      D^[Orleans,Pau]        :=656.0;
      D^[Orleans,Perigueux]  :=426.0;
      D^[Orleans,Poitiers]   :=256.0;
      D^[Orleans,Reims]      :=265.0;
      D^[Orleans,Rennes]     :=297.0;
      D^[Orleans,Rouen]      :=217.0;
      D^[Orleans,StEtienne]  :=503.0;
      D^[Orleans,Strasbourg] :=611.0;
      D^[Orleans,Toulouse]   :=574.0;
      D^[Orleans,Tours]      :=112.0;
      D^[Orleans,Troyes]     :=188.0;
{     distances between Paris and other towns }
      D^[Paris,Pau]        :=778.0;
      D^[Paris,Perigueux]  :=556.0;
      D^[Paris,Poitiers]   :=506.0;
      D^[Paris,Reims]      :=142.0;
      D^[Paris,Rennes]     :=348.0;
      D^[Paris,Rouen]      :=139.0;
      D^[Paris,StEtienne]  :=516.0;
      D^[Paris,Strasbourg] :=488.0;
      D^[Paris,Toulouse]   :=705.0;
      D^[Paris,Tours]      :=234.0;
      D^[Paris,Troyes]     :=178.0;
{     distances between Pau and other towns }  
      D^[Pau,Perigueux]  :=261.0;
      D^[Pau,Poitiers]   :=427.0;
      D^[Pau,Reims]      :=920.0;
      D^[Pau,Rennes]     :=631.0;
      D^[Pau,Rouen]      :=848.0;
      D^[Pau,StEtienne]  :=700.0;
      D^[Pau,Strasbourg] :=1171.0;
      D^[Pau,Toulouse]   :=173.0;
      D^[Pau,Tours]      :=544.0;
      D^[Pau,Troyes]     :=959.0;
{     distances between Perigueux and other towns }
      D^[Perigueux,Poitiers]   :=196.0;
      D^[Perigueux,Reims]      :=635.0;
      D^[Perigueux,Rennes]     :=433.0;
      D^[Perigueux,Rouen]      :=617.0;
      D^[Perigueux,StEtienne]  :=504.0;
      D^[Perigueux,Strasbourg] :=807.0;
      D^[Perigueux,Toulouse]   :=248.0;
      D^[Perigueux,Tours]      :=340.0;
      D^[Perigueux,Troyes]     :=514.0;
{     distances between Poitiers and other towns }
      D^[Poitiers,Reims]      :=521.0;
      D^[Poitiers,Rennes]     :=259.0;
      D^[Poitiers,Rouen]      :=424.0;
      D^[Poitiers,StEtienne]  :=431.0;
      D^[Poitiers,Strasbourg] :=741.0;
      D^[Poitiers,Toulouse]   :=436.0;
      D^[Poitiers,Tours]      :=104.0;
      D^[Poitiers,Troyes]     :=444.0;
{     distances between Reims and other towns }  
      D^[Reims,Rennes]     :=483.0;
      D^[Reims,Rouen]      :=230.0;
      D^[Reims,StEtienne]  :=533.0;
      D^[Reims,Strasbourg] :=346.0;
      D^[Reims,Toulouse]   :=840.0;
      D^[Reims,Tours]      :=369.0;
      D^[Reims,Troyes]     :=123.0;
{     distances between Rennes and other towns }
      D^[Rennes,Rouen]      :=300.0;
      D^[Rennes,StEtienne]  :=660.0;
      D^[Rennes,Strasbourg] :=829.0;
      D^[Rennes,Toulouse]   :=676.0;
      D^[Rennes,Tours]      :=236.0;
      D^[Rennes,Troyes]     :=485.0;
{     distances between Rouen and other towns }
      D^[Rouen,StEtienne]  :=644.0;
      D^[Rouen,Strasbourg] :=576.0;
      D^[Rouen,Toulouse]   :=792.0;
      D^[Rouen,Tours]      :=277.0;
      D^[Rouen,Troyes]     :=317.0;
{     distances between Saint-Etienne and other towns }
      D^[StEtienne,Strasbourg] :=536.0;
      D^[StEtienne,Toulouse]   :=527.0;
      D^[StEtienne,Tours]      :=424.0;
      D^[StEtienne,Troyes]     :=399.0;
{     distances between Strasbourg and other towns }
      D^[Strasbourg,Toulouse]  :=1015.0;
      D^[Strasbourg,Tours]     :=715.0;
      D^[Strasbourg,Troyes]    :=365.0;
{     distance between Toulouse et Tours, Troyes }
      D^[Toulouse,Tours] :=535.0;
      D^[Toulouse,Troyes]:=786.0;
{     distance Tours-Troyes }        
      D^[Tours,Troyes]:=300.0;
{     end by symmetry }
      D^[1,1]:=0.0;
      for I:=2 to NV do
      begin
        for J:=1 to I-1 do D^[I,J]:=D^[J,I];
        D^[I,I]:=0.0
      end
    End;                            

{main program}
BEGIN

    New(CH); New(Vs); New(ICO); New(D); New(L); New(C); New(ID);
    writeln;
    writeln('  Distances between towns (km):');
    writeln('      1 : as the crow flies');
    writeln('      2 : by road');
    write('  Your choice ( 1 or 2 ): '); readln(choice);
    writeln;	
    write('  Starting number of iterations: '); readln(nbreiter);
    write('  Maximum number of iterations : '); readln(nmax);
    write('  Increment value of iterations: '); readln(diter);

    Assign(fp1,'villes.lst'); Rewrite(fp1);

    if choice=1 then
    begin
      writeln(fp1);
      writeln(fp1,'        SALESMAN''S  PROBLEM');
      writeln(fp1);
      writeln(fp1,'   NO     TOWN          COORDINATES (KM)');
      writeln(fp1,'                         EAST SOUTH');
      writeln(fp1)
    end;
    {NV = Number of towns}
    Assign(fp,'villes.dat'); Reset(fp);
    Readln(fp,NV);
    for i:=1 to NV do
    begin
      readln(fp, CH^[i], ICO^[i,1], ICO^[i,2]);
      if choice = 1 then
        writeln(fp1, i:4, CH^[i]:16,'   ',ICO^[i,1]:4,' ',ICO^[i,2]:4)
    end;
    Close(fp);

    for i:=1 to NV-1 do
      for j:=i+1 to NV do
      begin
        X:=ICO^[i,1]-ICO^[j,1];
        Y:=ICO^[i,2]-ICO^[j,2];
        D^[i,j]:=sqrt(X*X+Y*Y);
        ID^[i,j]:= Round(D^[i,j]);
        D^[j,i]:=D^[i,j];
        ID^[j,i]:= Round(D^[i,j])
      end;
    if choice = 1 then               {By air}
    begin
      writeln(fp1);
      writeln(fp1,'  MATRIX OF DISTANCES (KM)');
      for i:=1 to NV do
      begin
        for j:=1 to NV do write(fp1,' ',ID^[i,j]);
        writeln(fp1)
      end
    end;
    if choice = 2 then
    begin
      InitDist1;    {By road}
      InitDist2
    end;
{   Main loop }
    DDMIN:=35000.0;
    while nbreiter < nmax do
    begin
      TIME;
{     aff_L(0); }
      ANNEAL(nbreiter, itermin);
      nbreiter := nbreiter + diter
    end;
{   print results }
    writeln(fp1);
    writeln(fp1,'   Shortest Itinerary:');
    writeln(fp1);
    i:=1;
    while i < NV do
    begin
      writeln(fp1,' ',Vs^[L^[i]],' -->  ', Vs^[L^[i+1]],' -->');
      i:=i+2
    end;
    writeln(fp1);
    writeln(fp1,'   NITER= ',itermin,'   DMIN= ', DELTAT:8:2,' KM.');

    writeln;
    writeln('  Results in file villes.lst'); 
    writeln;

    Close(fp1); 
    ReadKey;

    Dispose(CH); Dispose(Vs); Dispose(ICO); Dispose(D); Dispose(L);
    Dispose(C); Dispose(ID);
    DoneWinCrt

END.

{End of file tanneal.pas}