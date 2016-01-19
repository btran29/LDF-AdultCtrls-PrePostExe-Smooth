Sub perimed_data_conversion()
'
' perimed_data_conversion Macro
'
' Keyboard Shortcut: Ctrl+q
'
    ActiveCell.FormulaR1C1 = "=(RC3-R" & ActiveCell.Row & "C3)*86000"
    ActiveCell.Offset(0, -1).Range("A1").Select
    Range(Selection, Selection.End(xlDown)).Select
    rCount = Selection.Rows.Count
    ActiveCell.Offset(0, 1).Select
    Dim aFillRange As String
    Let aFillRange = "A" & 1 & ":" & "A" & rCount
    Selection.AutoFill Destination:=ActiveCell.Range(aFillRange)
    
    ActiveCell.Offset(-1, 0).Range("A1").Select
    ActiveCell.FormulaR1C1 = "time"
    ActiveCell.Offset(0, 1).Range("A1").Select
    ActiveCell.FormulaR1C1 = "temp"
    ActiveCell.Offset(0, 1).Range("A1").Select
    ActiveCell.FormulaR1C1 = "pressure"
    ActiveCell.Offset(0, 1).Range("A1").Select
    ActiveCell.FormulaR1C1 = "perfusion"
    ActiveCell.Offset(0, 1).Range("A1").Select
	ActiveCell.FormulaR1C1 = ActiveWorkbook.FullName
	
    ActiveCell.Offset(1, -6).Range("A1").Select
    Range(Selection, Selection.End(xlDown)).Select
    Application.CutCopyMode = False
    Selection.Copy
    ActiveCell.Offset(0, 3).Range("A1").Select
    ActiveSheet.Paste
    
    ActiveCell.Offset(0, -4).Range("A1").Select
    Range(Selection, Selection.End(xlDown)).Select
    Application.CutCopyMode = False
    Selection.Copy
    ActiveCell.Offset(0, -3).Range("A1").Select
    ActiveCell.Offset(0, 8).Range("A1").Select
    ActiveSheet.Paste
    
    ActiveCell.Offset(0, -3).Range("A1").Select
    Range(Selection, Selection.End(xlDown)).Select
    Application.CutCopyMode = False
    Selection.Copy
    ActiveCell.Offset(0, -5).Range("A1").Select
    ActiveCell.Offset(0, 9).Range("A1").Select
    ActiveSheet.Paste
    
    ActiveCell.Offset(-1, -3).Range("A1").Select
    Range(Selection, Selection.End(xlToRight)).Select
    Range(Selection, Selection.End(xlDown)).Select
    
    Selection.Copy
    Workbooks.Add
    Selection.PasteSpecial Paste:=xlPasteValues, Operation:=xlNone, SkipBlanks _
        :=False, Transpose:=False
End Sub
