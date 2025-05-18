{*********************************************
*   Error messages used by program EF1.pas   *
*********************************************}
UNIT Errors;

INTERFACE

Type Error_String = String[41];   {Longest error message}
     Error_Type = Record
                  Error_ID:Byte;
                  Error_Msg:Error_String;
                  End;

Const Run_Errors:Array[1..11] of Error_Type=(
                (Error_ID:$01; Error_Msg:'Floating point overflow'),
                (Error_ID:$02; Error_Msg:'Division by zero'),
                (Error_ID:$03; Error_Msg:'Sqrt of a negative number'),
                (Error_ID:$04; Error_Msg:'Log of a negative number'),
                (Error_ID:$10; Error_Msg:'String length error'),
                (Error_ID:$11; Error_Msg:'Invalid string index'),
                (Error_ID:$90; Error_Msg:'Index out of range'),
                (Error_ID:$91; Error_Msg:'Scaler or subrange out of range'),
                (Error_ID:$92; Error_Msg:'Out of integer range'),
                (Error_ID:$F0; Error_Msg:'Overlay file not found'),
                (Error_ID:$FF; Error_Msg:'Out of memory'));
Const IO_Errors:Array [1..16] of Error_Type=(
                (Error_ID:$01; Error_Msg:'File does not exist'),
                (Error_ID:$02; Error_Msg:'File not open for input'),
                (Error_ID:$03; Error_Msg:'File not open for output'),
                (Error_ID:$04; Error_Msg:'File not open'),
                (Error_ID:$10; Error_Msg:'Error in numeric format'),
                (Error_ID:$20; Error_Msg:'Operation not allowed on a logical device'),
                (Error_ID:$21; Error_Msg:'Not allowed in direct mode'),
                (Error_ID:$22; Error_Msg:'Assign to standard file'),
                (Error_ID:$90; Error_Msg:'Record length mismatch'),
                (Error_ID:$91; Error_Msg:'Seek beyond end of file'),
                (Error_ID:$99; Error_Msg:'Unexpected end of file'),
                (Error_ID:$F0; Error_Msg:'Disk write error'),
                (Error_ID:$F1; Error_Msg:'Directory is full'),
                (Error_ID:$F2; Error_Msg:'File size overflow'),
                (Error_ID:$F3; Error_Msg:'Too many open files'),
                (Error_ID:$FF; Error_Msg:'File disappeared'));
IMPLEMENTATION

END.
