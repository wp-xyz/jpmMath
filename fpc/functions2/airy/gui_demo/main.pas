unit main;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, ExtCtrls, StdCtrls,
  Spin, TAGraph, TASeries, TAFuncSeries, TAChartListbox,
  jpmTypes, jpmSpecialFunc;

type

  { TMainForm }

  TMainForm = class(TForm)
    Chart: TChart;
    AiSeries: TFuncSeries;
    BiSeries: TFuncSeries;
    ChartListbox: TChartListbox;
    DerivAiSeries: TFuncSeries;
    DerivBiSeries: TFuncSeries;
    feXMin: TFloatSpinEdit;
    feXMax: TFloatSpinEdit;
    lblXMin: TLabel;
    lblXMax: TLabel;
    Panel1: TPanel;
    Panel2: TPanel;
    procedure AiSeriesCalculate(const AX: Double; out AY: Double);
    procedure BiSeriesCalculate(const AX: Double; out AY: Double);
    procedure DerivAiSeriesCalculate(const AX: Double; out AY: Double);
    procedure DerivBiSeriesCalculate(const AX: Double; out AY: Double);
    procedure feXMaxChange(Sender: TObject);
    procedure feXMinChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    X, Ai, Bi, Ad, Bd: Float;

  public

  end;

var
  MainForm: TMainForm;

implementation

{$R *.lfm}

{ TMainForm }

procedure TMainForm.AiSeriesCalculate(const AX: Double; out AY: Double);
begin
  if X <> AX then
  begin
    X := AX;
    AiryA(AX, Ai, Bi, Ad, Bd);
  end;
  AY := Ai;
end;

procedure TMainForm.BiSeriesCalculate(const AX: Double; out AY: Double);
begin
  if X <> AX then
  begin
    X := AX;
    AiryA(AX, Ai, Bi, Ad, Bd);
  end;
  AY := Bi;
end;

procedure TMainForm.DerivAiSeriesCalculate(const AX: Double; out AY: Double);
begin
  if X <> AX then
  begin
    X := AX;
    AiryA(AX, Ai, Bi, Ad, Bd);
  end;
  AY := Ad;
end;

procedure TMainForm.DerivBiSeriesCalculate(const AX: Double; out AY: Double);
begin
  if X <> AX then
  begin
    X := AX;
    AiryA(AX, Ai, Bi, Ad, Bd);
  end;
  AY := Bd;
end;

procedure TMainForm.feXMaxChange(Sender: TObject);
begin
  Chart.BottomAxis.Range.Max := feXMax.Value;
  Chart.BottomAxis.Range.UseMax := true;
end;

procedure TMainForm.feXMinChange(Sender: TObject);
begin
  Chart.BottomAxis.Range.Min := feXMin.Value;
  Chart.BottomAxis.Range.UseMin := true;
end;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  Chart.BottomAxis.Range.Min := feXMin.Value;
  Chart.BottomAxis.Range.Max := feXMax.Value;
  Chart.BottomAxis.Range.UseMin := true;
  Chart.BottomAxis.Range.UseMax := true;
end;

end.

