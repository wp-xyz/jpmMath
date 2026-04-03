{*******************************************************************
!*      Purpose: This program computes the first NT zeros of Airy  * 
!*               functions Ai(x) and Ai'(x), and the associated    *
!*               values of Ai(a') and Ai'(a), and the first NT     *
!*               zeros of Airy functions Bi(x) and Bi'(x), and     *
!*               the associated values of Bi(b') and Bi'(b) using  * 
!*               subroutine AIRYZO                                 *
!*      Input :  NT    --- Total number of zeros                   *
!*               KF    --- Function code                           *
!*                         KF=1 for Ai(x) and Ai'(x)               *
!*                         KF=2 for Bi(x) and Bi'(x)               *
!*      Output:  XA(m) --- a, the m-th zero of Ai(x) or            *
!*                         b, the m-th zero of Bi(x)               *
!*               XB(m) --- a', the m-th zero of Ai'(x) or          *
!*                         b', the m-th zero of Bi'(x)             *
!*               XC(m) --- Ai(a') or Bi(b')                        *
!*               XD(m) --- Ai'(a) or Bi'(b)                        *
!*                         ( m --- Serial number of zeros )        * 
!*      Example: NT=5                                              *
!*                                                                 *
!*      m         a            Ai'(a)         a'          Ai(a')   *
!*     ----------------------------------------------------------- *
!*      1    -2.33810741     .70121082   -1.01879297    .53565666  *
!*      2    -4.08794944    -.80311137   -3.24819758   -.41901548  *
!*      3    -5.52055983     .86520403   -4.82009921    .38040647  *
!*      4    -6.78670809    -.91085074   -6.16330736   -.35790794  *
!*      5    -7.94413359     .94733571   -7.37217726    .34230124  *
!*                                                                 *
!*      m         b            Bi'(b)         b'          Bi(b')   *
!*     ----------------------------------------------------------- *
!*      1    -1.17371322     .60195789   -2.29443968   -.45494438  *
!*      2    -3.27109330    -.76031014   -4.07315509    .39652284  *
!*      3    -4.83073784     .83699101   -5.51239573   -.36796916  *
!*      4    -6.16985213    -.88947990   -6.78129445    .34949912  *
!*      5    -7.37676208     .92998364   -7.94017869   -.33602624  *
!*                                                                 *
!* --------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special Func-   *
!*             tions jin.ece.uiuc.edu/routines/routines.html".     *
!*                                                                 *
!*                               TPW Release By J-P Moreau, Paris. *
!*                                       (www.jpmoreau.fr)         *
!******************************************************************}
Program AiryZerosLaz;

uses
  jpmTypes, jpmSpecialFunc;

var
  K, KF, NT: integer;
  XA, XB, XC, XD: TFloatArray;

begin
  WriteLn;
  WriteLn(' KF=1 for Ai(x) and Ai''(x); KF=2 for Bi(x) and Bi''(x)');
  WriteLn(' NT is the number of the zeros.');
  Write(' Please enter NT: ');
  ReadLn(NT);

  Writeln;
  for KF := 1 to 2 do
  begin
    if (KF=1) then
      WriteLn('  m        a         Ai''(a)         a''         Ai(a'')')
    else if (KF=2) then
      WriteLn('  m        b         Bi''(b)         b''         Bi(b'')');
    WriteLn('----------------------------------------------------------');
        
    AiryZeros(NT, KF, XA, XB, XC, XD);
        
    for K := 0 to NT-1 do
    begin
      WriteLn(K:3, ' ', XA[K]:12:8, ' ', XD[K]:12:8, ' ' ,XB[K]:12:8, ' ', XC[K]:12:8)
    end;
    WriteLn;
  end;
  Write('Press ENTER to close...');
  ReadLn;
end.

{end of file AiryZerosLaz.pas}
