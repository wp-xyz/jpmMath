{************************************************}
  {                                                }
  {   Turbo Pascal for Windows                     }
  {   Tips & Techniques Demo Program               }
  {   Copyright (c) 1991 by Borland International  }
  {                                                }
  { ---------------------------------------------- }
  {         Improved version by j-P Moreau, Paris  }
  {                   (www.jpmoreau.fr)            }
  {************************************************}

  {This unit allows to send graphic commands to a printer
   See example of use in program billard.pas}

  UNIT WinPrint;
  {$R PRINTER.RES}

  interface

  USES WinTypes, WinProcs, WOBJECTS, Strings;

  TYPE

  { TComboXferRec }
  { The transfer buffer used FOR the ComboBox in the TPrinterInfo method
  SelectPrinter.  The fields, Strings AND Selection, are set up in the
  TPrinterInfo CONSTRUCTOR Init.  The routine GetCurrentPrinter is used
  to find current printing device which is placed in Selection.  And the
  routine GetPrinterTypes is used to fill out the Strings field.}

  TComboXferRec = RECORD
    Strings: PStrCollection;
    Selection: ARRAY[0..80] OF Char;
  END;

  { TAbortDialog }
  { A descendant OF TDialog used FOR the Abort Dialog seen when printing is
  in progress. The AbortDialog is installed as a data field OF TPrinterInfo
  and is initialized AND displayed in its StartDoc method. The ENDDoc
  method will Close the dialog if necessary.}

  PAbortDialog = ^TAbortDialog;
  TAbortDialog = OBJECT(TDlgWindow)
    PROCEDURE SetUpWindow; virtual;
    PROCEDURE WMCommand(VAR Msg: TMessage);
    virtual wm_First + wm_Command;
  END;
  
  { TPrinterInfo 
  The controlling object for printing.  It is intended that this OBJECT be
  initialized as a data field of a TWindow OR TApplication descENDant. This
  printing object must be used OWL based applications. The data fields are
  not supposed to be used directly but may need TO be accessed in special
  situations.  PrintDC and Error are the two most likely TO be used WITHout
  a specific method call.  The description OF the data fields are as
  follows.

  -AbortDialog holds a pointer to the abort dialog when it valid.  It is
  valid only after a call to the method StartDoc AND before the call TO the
  method EndDoc.
  
  -AbortCallBackProc holds the address of the Abort Dialog's callback
  FUNCTION.  It's definition is found in the function AbortCallBack in the
  implementation section of this unit.
  
  -SelectDialog is a pointer to the dialog used when selecting the current
  printer. To be used when overriding the function of the SelectPrinter
  method.

  -SelectInfo is the transfer record used in SelectDialog.  Holds
  descriptions of all printers available and the currently selected printer.

  -Driver, PrinterType, Port are null terminated strings holding information
  relevant TO the current printer.

  -DriverHandle is a handle to the library of the current printer driver. It
  is setup in Init constructor and is freed in the Done DESTRUCTor.  It is
  used FOR setting up the DeviceMode configuration call.

  -PrintDC is the device control established for printing. It is created by
  the StartDoc method and valid until the EndDoc method call. May be
  accessed directly or by the GetPrinterDC method call.

  -ErrOR holds the results of printer escape calls.  If an errOR occurs, the
  result is placed here.  Is tested to determine if further printing output
  is appropriate.

  -ExtDeviceMode holds the ExtDeviceMode procedure used FOR retrieving,
  installing, and prompting for printing configurations.

  -DeviceModeVar holds the DeviceMode procedure used for prompting the
  user FOR printer configurations.
  }

  PPrinterInfo = ^TPrinterInfo;
  TPrinterInfo = OBJECT
    AbortDialog: PAbortDialog;
    AbortCallBackProc: TFarProc;
    SelectDialog: PDialog;
    SelectInfo: TComboXferRec;
    Driver,PrinterType,Port: PChar;
    DriverHandle: THandle;
    PrintDC: HDC;
    Error,LineHeight,LinesPerPage: Integer;
    id_error:Boolean;
    ExtDeviceMode: TExtDeviceMode;
    DeviceModeVAR: TDeviceMode;
    CurrentLine,RasterCaps: INTEGER;
    CONSTRUCTOR Init;
    DESTRUCTOR Done;
    PROCEDURE SelectPrinter; virtual;
    FUNCTION GetPrinterDC: HDC;
    PROCEDURE InitPrintParams;
    PROCEDURE PrnLine(P:PChar);
    PROCEDURE DeviceMode;
    FUNCTION BitMapCapable: BOOLEAN;
    FUNCTION BandingRequired: BOOLEAN;
    PROCEDURE StartDoc(Name: PChar); virtual;
    PROCEDURE NewFrame; virtual;
    PROCEDURE NextBand(VAR R:TRect); virtual;
    PROCEDURE ENDDoc; virtual;
  END;


  VAR
  PrinterAbort,Printing: Boolean;
  { Holds true when the user has aborted printing. }

  CONST
  LinesAtBottom = 3;     { Nombre de ligne dans la marge du bas }

  implementation

  CONST
  LeftMargin = 5;        { Largeur de la marge gauche }
  LinesAtTop = 3;        { Nombre de ligne dans la marge haute }
  MinimumLines = LinesAtTop + LinesAtBottom + 1;
  id_ComboBox = 101;
  { ID for the ComboBox used for Selecting the current printer }

  VAR
  AbortWindow: HWnd;
  { Window handle for the Abort Dialog.  It is used by the
  AbortCallBackProc.}

  FUNCTION GetItem(VAR S: PChar): PChar;
  { Retrieves comma separated data from a null terminated STRING. It
  returns the first data item AND advances the pointer S TO the next
  data item in the string.}
  VAR
  P: PChar;
  I: Integer;

  BEGIN
    I:=0;
    WHILE (S[I]<>',') AND (S[I]<>#0) DO
    inc(I);
    S[I]:=#0;
    GetMem(P, Strlen(S)+1);
    StrCopy(P,S);
    GetItem:=P;
    IF S[0]<>#0 THEN S:=@S[I+1];
  END;

  PROCEDURE GetPrinterTypes(VAR PrinterTypes: PStrCollection);
  { Retrieves all the device TYPEs from the WIN.INI AND places this
  information into the PStrCollection parameter.}
  VAR
  Buffer, BufferItem: PChar;
  Item: PChar;
  Count, I: Integer;

  BEGIN
    New(PrinterTypes, init(5,1));
    GetMem(Buffer, 1024);
    Count:=GetProfileString('devices', nil, ',,', Buffer, 1024);
    BufferItem:=Buffer;
    I:=0;
    WHILE I<Count DO
    BEGIN
      GetMem(Item, StrLen(BufferItem)+1);
      StrCopy(Item, BufferItem);
      PrinterTypes^.Insert(Item);
      WHILE (BufferItem[i]<>#0) AND (I<Count) DO
      inc(I);
      inc(I);
      IF BufferItem[I]=#0 THEN I:=Count;
      IF I<Count THEN
      BEGIN
        BufferItem:=@BufferItem[I];
        Count:=Count-I;
        I:=0;
      END;
    END;
    FreeMem(Buffer, 1024);
  END;

  PROCEDURE GetCurrentPrinter(VAR Driver, PrinterType, Port: PChar);
  { Retrieves the current printing device information from the WIN.INI
  file.}
  VAR
  ProfileInfo, CurrentItem: PChar;
  BEGIN
    GetMem(ProfileInfo, 80+1);
    GetProfileString('windows', 'device', ',,', ProfileInfo, 80);
    CurrentItem:=ProfileInfo;
    PrinterType:=GetItem(CurrentItem);
    Driver:=GetItem(CurrentItem);
    Port:=GetItem(CurrentItem);
    FreeMem(ProfileInfo, 80+1);
  END;
  
  PROCEDURE GetPrinter(PrinterType: PChar; VAR Driver, Port: PChar);
  { Given a PrinterType string, this PROCEDURE returns the appropriate
  driver and port information.}
  
  VAR
  ProfileInfo, CurrentItem: PChar;
  
  BEGIN
    GetMem(ProfileInfo, 80+1);
    GetProfileString('devices', PrinterType, ',', ProfileInfo, 80);
    CurrentItem:=ProfileInfo;
    Driver:=GetItem(CurrentItem);
    Port:=GetItem(CurrentItem);
  END;
  
  PROCEDURE TAbortDialog.SetUpWindow;
  { Initializes PrinterAbort and AbortWindow. THEN set the focus TO the
  AbortDialog.}
  BEGIN
    PrinterAbort:=false;
    SetFocus(HWindow);
    AbortWindow:=HWindow;
  END;
  
  PROCEDURE TAbortDialog.WMCommand(VAR Msg: TMessage);
  { If any commAND messages occur, a user abort has taken place.  Normally,
  this will include pressing ENTER, ESCAPE, the SPACEBAR  OR clicking the
  mouse on the Abort Dialog's Escape button.}
  BEGIN
    PrinterAbort:=true;
  END;
  
  FUNCTION AbortCallBack(DC: HDC; Code: Integer): Bool; export;
  { While printing is taking place, checks TO see IF PrinterAbort is
  true.  Otherwise messages are passed on.}
  VAR
  Msg: TMsg;
  BEGIN
    WHILE (NOT PrinterAbort) AND PeekMessage(Msg, 0, 0, 0, pm_Remove) DO
    IF NOT IsDialogMessage(AbortWindow, Msg) THEN
    BEGIN
      TranslateMessage(Msg);
      DispatchMessage(Msg);
    END;
    IF PrinterAbort THEN AbortCallBack:=false ELSE AbortCallBack:=true;
  END;
  
  CONSTRUCTOR TPrinterInfo.Init;
  { Gets the current printer information (Type, Driver, & Port) and
  the printer types currently available.  THEN retrieves the
  ExtDeviceMode and DeviceModeVar address from the current printer's
  library.}
  VAR
  I: Integer;
  FullDriverName: PChar;
  P: TFarProc;
  
  BEGIN
    GetCurrentPrinter(Driver, PrinterType, Port);
    FOR I:= 0 TO StrLen(PrinterType) DO
    SelectInfo.Selection[I]:=PrinterType[I];
    GetPrinterTypes(SelectInfo.Strings);

  
    GetMem(FullDriverName, 12+1);
    StrLCat(StrCopy(FullDriverName, Driver), '.DRV', 12);
    DriverHandle:=LoadLibrary(FullDriverName);
    FreeMem(FullDriverName, 12+1);
  
  
    P:=GetProcAddress(DriverHandle, 'ExtDeviceMode');
    ExtDeviceMode:=TExtDeviceMode(P);
    P:=GetProcAddress(DriverHandle, 'DeviceMode');
    DeviceModeVAR:=TDeviceMode(P);
    PrintDC:=0; id_error:=FALSE
  END;
  
  DESTRUCTOR TPrinterInfo.Done;
  { Frees up the library taken in the constructor Init.}
  BEGIN
    FreeLibrary(DriverHandle);
  END;
  
  PROCEDURE TPrinterInfo.SelectPrinter;
  { Displays a Printer Select dialog called PISELECT AND changes the
  current printer information as is done in Init.}
  VAR
  FullDriverName: PChar;
  P: TFarProc;
  ComboBox: PComboBox;
  
  BEGIN
    new(SelectDialog, Init(Application^.MainWindow,
    'PISELECT'));
    New(ComboBox, InitResource(SelectDialog, id_ComboBox, 80));
    SelectDialog^.TransferBuffer:=@SelectInfo;
    IF Application^.ExecDialog(SelectDialog) = id_Ok THEN
    BEGIN
      FreeLibrary(DriverHandle);
      IF PrintDC<>0 THEN DeleteDC(PrintDC);
      FreeMem(PrinterType, StrLen(PrinterType)+1);
      GetMem(PrinterType, StrLen(@SelectInfo.Selection)+1);
      StrCopy(PrinterType, @SelectInfo.Selection);
      FreeMem(Driver, StrLen(Driver)+1);
      FreeMem(Port, StrLen(Port)+1);
      GetPrinter(PrinterType, Driver, Port);
      GetMem(FullDriverName, 12+1);
      StrLCat(StrCopy(FullDriverName, Driver), '.DRV', 12);
      DriverHandle:=LoadLibrary(FullDriverName);
      FreeMem(FullDriverName, 12+1);
      P:=GetProcAddress(DriverHandle, 'ExtDeviceMode');
      ExtDeviceMode:=TExtDeviceMode(P);
      P:=GetProcAddress(DriverHandle, 'DeviceMode');
      DeviceModeVAR:=TDeviceMode(P);
    END;
  END;

  FUNCTION TPrinterInfo.GetPrinterDC: HDC;
  { Retrieves the Device control associated WITH the printer.  May only be
  called after a call to the StartDoc method. }
  BEGIN
    GetPrinterDC:=PrintDC;
  END;
  
  { Initialise les paramètres globaux d'impression }
  PROCEDURE TPrinterInfo.InitPrintParams;
  VAR
  TM: TTextMetric;
  PageWidth, PageHeight: Integer;
  BEGIN
    GetTextMetrics(PrintDc, TM);
    PageWidth := GetDeviceCaps(PrintDc, HorzRes); { Non utilisé }
    PageHeight := GetDeviceCaps(PrintDc, VertRes);
    LineHeight := TM.tmHeight + TM.tmHeight DIV 2;
    IF LineHeight <= 0 THEN
    LineHeight := 10;  { Empêche la division par zéro }
    LinesPerPage := PageHeight DIV LineHeight;
    IF LinesPerPage < MinimumLines THEN
    LinesPerPage := MinimumLines;
    CurrentLine := LinesAtTop
  END;

  PROCEDURE TPrinterInfo.StartDoc(Name: PChar);
  { Called immediately before printing is to begin.  Establishes the
  device control.  Sets up the Abort Dialog. And sEND the STARTDOC
  escape call.}
  BEGIN
    Error:=0;  Printing:=TRUE;
    PrintDC:=CreateDC(Driver, PrinterType, Port, NIL);
    IF LowMemory THEN
    AbortDialog:=Nil
    ELSE
    BEGIN
      new(AbortDialog, Init(Application^.MainWindow, 'PIABORT'));
      AbortDialog^.Create;
    END;
    IF AbortDialog<>Nil THEN
    BEGIN
      AbortCallBackProc:=MakeProcInstance(@AbortCallBack, HInstance);
      Escape(PrintDC, SETABORTPROC, 0, AbortCallBackProc, NIL);
    END;
    RasterCaps:=GetDeviceCaps(PrintDC, WINTYPES.RASTERCAPS);
    Error:=Escape(PrintDC, WINTYPES.STARTDOC, StrLen(Name), Name, NIL);
    InitPrintParams
  END;

  PROCEDURE TPrinterInfo.NewFrame;
  { Sends the NEWFRAME escape call, performs appropriate error
  checking and changes page }
  BEGIN
    CurrentLine := LinesAtTop;    { changement de page   }
    IF id_errOR THEN exit;        { erreur déjà signalée }
    IF Error>=0 THEN                        { cas normal }
    Error:=Escape(PrintDC, WINTYPES.NEWFRAME, 0, NIL, NIL);
    IF Error<0 THEN
    BEGIN
      CASE ErrOR OF
        SP_ERROR: MessageBox(Application^.MainWindow^.HWindow,
        'Erreur générale imprimante', NIL, mb_Ok OR mb_IconStop);
        SP_OUTOFDISK: MessageBox(Application^.MainWindow^.HWindow,
        'Plus d''espace disque pour stocker l''impression', NIL, mb_Ok OR mb_IconStop);
        SP_OUTOFMEMORY: MessageBox(Application^.MainWindow^.HWindow,
        'Plus de mémoire pour stocker l''impression', NIL, mb_Ok OR mb_IconStop);
        SP_USERABORT: MessageBox(Application^.MainWindow^.HWindow,
        'Impression arrêtée par l''utilisateur', NIL, mb_Ok OR mb_IconStop);
        ELSE
        MessageBox(Application^.MainWindow^.HWindow,
        'Impression arrêtée', NIL, mb_OK OR mb_IconStop);
      END;
      id_error:=TRUE;
      AbortDialog^.Destroy;   {fermeture dialogue impression }
      PrinterAbort:=TRUE
    END
  END;

  { Imprime une ligne, dont l'adresse est donnée par P }
  PROCEDURE TPrinterInfo.PrnLine(P: PChar);
  BEGIN
    Inc(CurrentLine);
    TextOut(PrintDc, LeftMargin, CurrentLine * LineHeight, P, StrLen(P));
    {IF CurrentLine >= LinesPerPage - LinesAtBottom THEN
    NewFrame }
  END;

  PROCEDURE TPrinterInfo.NextBand(VAR R:TRect);
  { When Bitmap banding is required, this routine returns the next
  rectangular region to be printed.  This method is not required but
  can speed up printing bitmaps.}
  BEGIN
    IF Error>=0 THEN
    Error:=Escape(PrintDC, WINTYPES.NEXTBAND, 0, NIL, @R);
    IF Error<0 THEN
    CASE ErrOR OF
      SP_ERROR: MessageBox(Application^.MainWindow^.HWindow,
      'General Printer Error', NIL, mb_Ok OR mb_IconStop);
      SP_OUTOFDISK: MessageBox(Application^.MainWindow^.HWindow,
      'No disk space for spooling', NIL, mb_Ok OR mb_IconStop);
      SP_OUTOFMEMORY: MessageBox(Application^.MainWindow^.HWindow,
      'No memory space for spooling', NIL, mb_Ok OR mb_IconStop);
      SP_USERABORT: MessageBox(Application^.MainWindow^.HWindow,
      'Printing Terminated by User', NIL, mb_Ok OR mb_IconStop);
      ELSE
      MessageBox(Application^.MainWindow^.HWindow,
      'Printing Halted', NIL, mb_OK OR mb_IconStop);
    END;
  END;

  PROCEDURE TPrinterInfo.ENDDoc;
  { Sends the ENDDOC escape call AND closes the Abort Dialog IF no errors
  have occurred.}
  BEGIN
    IF Error>=0 THEN
    Error:=Escape(PrintDC, WINTYPES.ENDDOC, 0, NIL, NIL);
    IF Error>=0 THEN
    BEGIN
      DeleteDC(PrintDC);
      IF AbortDialog<>Nil THEN AbortDialog^.CloseWindow;
    END;
    Printing:=FALSE
  END;
  
  PROCEDURE TPrinterInfo.DeviceMode;
  { Calls the printer driver's DeviceMode routine.  Normally displays a
  dialog allowing the user TO change the printer's configuration.}
  BEGIN
    DeviceModeVAR(Application^.MainWindow^.HWindow,
    DriverHandle, PrinterType, Port);
  END;
  
  FUNCTION TPrinterInfo.BitMapCapable: BOOLEAN;
  { Returns true if the current printing device can handle bitmap
  graphics.}
  BEGIN
    BitMapCapable:=(RasterCaps AND RC_BITBLT)<>0;
  END;

  FUNCTION TPrinterInfo.BandingRequired: BOOLEAN;
  { Returns true if banding of bitmap images will enhance printing speed.}
  BEGIN
    BandingRequired:=(RasterCaps AND RC_BANDING)<>0;
  END;

  END.

  { Here are the descriptions of the dialogs PIABORT and PISELECT found in
  the resources file PRINTER.RES

  PIABORT DIALOG DISCARDABLE LOADONCALL PURE MOVEABLE 44, 46, 175, 78
  STYLE WS_POPUP | WS_VISIBLE | WS_CAPTION | 0x80L
  CAPTION "Printing in Progress"
  BEGIN
  CONTROL "Press Escape to Halt Printing" 101, "STATIC", WS_CHILD |
  WS_VISIBLE, 37, 17, 98, 12
  CONTROL "Escape" 102, "BUTTON", WS_CHILD | WS_VISIBLE | WS_TABSTOP,
  73, 49, 40, 13
  END
  
  PISELECT DIALOG DISCARDABLE LOADONCALL PURE MOVEABLE 44, 37, 145, 85
  STYLE WS_POPUP | WS_VISIBLE | WS_CAPTION | 0x80L
  CAPTION "Select Printer"
  BEGIN
  CONTROL "COMBOBOX" 101, "COMBOBOX", WS_CHILD | WS_VISIBLE | WS_VSCROLL |
  0x101L, 26, 11, 84, 43
  CONTROL "Ok" 1, "BUTTON", WS_CHILD | WS_VISIBLE | WS_TABSTOP,
  29, 61, 40, 12
  CONTROL "Cancel" 2, "BUTTON", WS_CHILD | WS_VISIBLE | WS_TABSTOP,
  86, 61, 40, 12
  END
  }

{End of file winprint.pas}
