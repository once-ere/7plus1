(* tools/run_eval.wl -- evaluate 7+1.nb in place in a FRESH kernel of Wolfram 15.0.1 and save it with its outputs.

   Usage, from any notebook of a running Wolfram 15.0.1 front end:
       Get["D:/Developer/github/7+1/tools/run_eval.wl"]
   7+1.nb needs a fresh kernel: the tensor packages xAct (section 12) keep their definitions in the kernel, while a notebook that
   is closed and opened again gets a new private context (CellContext -> Notebook).  The script therefore
     1. closes any open copy of <repository>/7+1.nb and opens the file from disk, and deletes its old output and message cells;
     2. adds the kernel configuration "Eval71" (a local kernel, like "Local") to the front end and makes it the evaluator of the
        notebook, so that the evaluation starts a new kernel;
     3. appends a marker cell that writes build/eval_done.txt and then quits that kernel, evaluates all cells of the notebook
        (SelectionEvaluate) and waits for the marker (at most one hour; otherwise the evaluation is aborted);
     4. deletes the marker cell, gives the notebook the default evaluator again, saves and closes it, removes the kernel
        configuration "Eval71" again, and writes a short log to build/eval_log.txt (build/ is not part of the repository). *)
Module[{root, target, doneF, norm, nb, t0, t1, msgs, outs, closed = 0, ev = "Eval71", evs0, marker, finished},
  root = ParentDirectory[DirectoryName[$InputFileName]];
  target = FileNameJoin[{root, "7+1.nb"}];
  doneF = FileNameJoin[{root, "build", "eval_done.txt"}];
  Quiet[CreateDirectory[FileNameJoin[{root, "build"}]]];
  If[FileExistsQ[doneF], DeleteFile[doneF]];
  norm[f_String] := ToLowerCase[StringReplace[ExpandFileName[f], "\\" -> "/"]];
  Scan[With[{f = Quiet[NotebookFileName[#]]}, If[StringQ[f] && norm[f] === norm[target], NotebookClose[#]; closed++]] &, Notebooks[]];
  evs0 = EvaluatorNames /. Options[$FrontEnd, EvaluatorNames];
  If[! MemberQ[Keys[evs0], ev], SetOptions[$FrontEnd, EvaluatorNames -> Append[evs0, ev -> {"AutoStartOnLaunch" -> False}]]];
  t0 = AbsoluteTime[];
  nb = NotebookOpen[target];
  NotebookDelete[Cells[nb, CellStyle -> {"Output", "Message", "Print", "Echo"}]];
  SetOptions[nb, Evaluator -> ev];
  marker = "Export[" <> ToString[StringReplace[doneF, "\\" -> "/"], InputForm] <> ", DateString[], \"Text\"]; Quit[]";
  SelectionMove[nb, After, Notebook];
  NotebookWrite[nb, Cell[marker, "Input", CellTags -> {"evalDoneMarker7p1"}]];
  SelectionMove[nb, All, Notebook];
  SelectionEvaluate[nb];
  While[! FileExistsQ[doneF] && AbsoluteTime[] - t0 < 3600, Pause[2]];
  finished = FileExistsQ[doneF];
  If[! finished, FrontEndTokenExecute[nb, "EvaluatorAbort"]; Pause[10]];
  Pause[3];
  NotebookDelete[Cells[nb, CellTags -> "evalDoneMarker7p1"]];
  SetOptions[nb, Evaluator -> Inherited];
  msgs = Length[Cells[nb, CellStyle -> "Message"]];
  outs = Length[Cells[nb, CellStyle -> "Output"]];
  NotebookSave[nb];
  NotebookClose[nb];
  SetOptions[$FrontEnd, EvaluatorNames -> evs0];
  t1 = AbsoluteTime[];
  Export[FileNameJoin[{root, "build", "eval_log.txt"}],
    StringJoin["target: ", target, "\nversion: ", $Version, "\nkernel: a fresh kernel (front-end kernel configuration ", ev, ")",
      "\nseconds: ", ToString[N[t1 - t0]], "\nmessage cells: ", ToString[msgs], "\noutput cells: ", ToString[outs],
      "\ncompleted: ", ToString[finished], "\nopen copies closed first: ", ToString[closed], "\nfinished: ", DateString[], "\n"], "Text"];
  {target, N[t1 - t0], msgs, finished}]
