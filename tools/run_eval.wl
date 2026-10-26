(* tools/run_eval.wl -- evaluate 7+1.nb in place in Wolfram 15.0.1 and save it with its outputs.

   Usage, from any notebook of a running Wolfram 15.0.1 front end:
       Get["D:/Developer/github/7+1/tools/run_eval.wl"]
   The script closes any open copy of <repository>/7+1.nb (an open window would otherwise be reused with its old
   contents), opens the file from disk, deletes old output and message cells, evaluates every input cell with
   NotebookEvaluate[..., InsertResults -> True], saves and closes the notebook, and writes a short log to
   <repository>/build/eval_log.txt (build/ is not part of the repository). *)
Module[{root, target, norm, nb, t0, t1, msgs, closed = 0},
  root = ParentDirectory[DirectoryName[$InputFileName]];
  target = FileNameJoin[{root, "7+1.nb"}];
  norm[f_String] := ToLowerCase[StringReplace[ExpandFileName[f], "\\" -> "/"]];
  Scan[With[{f = Quiet[NotebookFileName[#]]}, If[StringQ[f] && norm[f] === norm[target], NotebookClose[#]; closed++]] &, Notebooks[]];
  t0 = AbsoluteTime[];
  nb = NotebookOpen[target];
  NotebookDelete[Cells[nb, CellStyle -> {"Output", "Message", "Print", "Echo"}]];
  NotebookEvaluate[nb, InsertResults -> True];
  msgs = Length[Cells[nb, CellStyle -> "Message"]];
  NotebookSave[nb];
  NotebookClose[nb];
  t1 = AbsoluteTime[];
  Quiet[CreateDirectory[FileNameJoin[{root, "build"}]]];
  Export[FileNameJoin[{root, "build", "eval_log.txt"}],
    StringJoin["target: ", target, "\nversion: ", $Version, "\nseconds: ", ToString[N[t1 - t0]],
      "\nmessage cells: ", ToString[msgs], "\nopen copies closed first: ", ToString[closed],
      "\nfinished: ", DateString[], "\n"], "Text"];
  {target, N[t1 - t0], msgs, closed}]
