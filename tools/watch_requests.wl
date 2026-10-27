(* tools/watch_requests.wl -- let the build of 7+1 request Wolfram 15.0.1 jobs without keyboard input.

   Load ONCE in a notebook of the running Wolfram 15.0.1 front end
   (C:\Program Files\Wolfram Research\Wolfram\15.0.1\WolframNB.exe):
       Get["D:/Developer/github/7+1/tools/watch_requests.wl"]
   Every 10 seconds a scheduled task reads <repository>/build/request.txt.  A line  "<token> <job> [<file>]"  with a
   token that has not been handled yet queues the job as an ordinary evaluation in the notebook "7+1 runner":
       push           ->  Get["<repository>/tools/git_push.wl"]   (log: build/git_push_log.txt)
       eval           ->  Get["<repository>/tools/run_eval.wl"]   (log: build/eval_log.txt)
       run <file.wl>  ->  Get["<repository>/<file.wl>"]          (a Wolfram Language file inside the repository)
   The handled token is recorded in build/request_done.txt.  Remove the task with  TaskRemove[watchTask7p1]. *)
Module[{root = ParentDirectory[DirectoryName[$InputFileName]], done},
  watchRoot7p1 = StringReplace[root, "\\" -> "/"];
  watchReq7p1 = watchRoot7p1 <> "/build/request.txt";
  watchDone7p1 = watchRoot7p1 <> "/build/request_done.txt";
  done = If[FileExistsQ[watchDone7p1], Quiet[StringSplit[StringTrim[ReadString[watchDone7p1]]]], {}];
  watchLast7p1 = If[ListQ[done] && Length[done] >= 1, First[done], ""]];
watchQueue7p1[file_String] := Module[{nb},
  nb = If[Head[watchNB7p1] === NotebookObject && MemberQ[Notebooks[], watchNB7p1], watchNB7p1,
     watchNB7p1 = CreateDocument[{TextCell["7+1 runner: jobs queued by tools/watch_requests.wl", "Text"]},
        WindowTitle -> "7+1 runner"]];
  SelectionMove[nb, After, Notebook];
  NotebookWrite[nb, Cell[BoxData[RowBox[{"Get", "[", "\"" <> file <> "\"", "]"}]], "Input"], All];
  SelectionEvaluate[nb]];
watchFile7p1[parts_List] := Module[{f},
  Switch[parts[[2]],
   "push", watchRoot7p1 <> "/tools/git_push.wl",
   "eval", watchRoot7p1 <> "/tools/run_eval.wl",
   "run", If[Length[parts] >= 3 && StringEndsQ[parts[[3]], ".wl"] && ! StringContainsQ[parts[[3]], ".."] &&
       FileExistsQ[f = watchRoot7p1 <> "/" <> parts[[3]]], f, None],
   _, None]];
watchCheck7p1[] := Module[{txt, parts, file},
  If[! FileExistsQ[watchReq7p1], Return[Null]];
  txt = Quiet[ReadString[watchReq7p1]];
  If[! StringQ[txt], Return[Null]];
  parts = StringSplit[StringTrim[txt]];
  If[Length[parts] < 2 || parts[[1]] === watchLast7p1, Return[Null]];
  watchLast7p1 = parts[[1]];
  file = watchFile7p1[parts];
  If[file =!= None, watchQueue7p1[file]];
  Export[watchDone7p1, StringRiffle[Join[parts, {If[file === None, "unknown-job", "queued"], DateString[]}], " "], "Text"]];
If[ValueQ[watchTask7p1], Quiet[TaskRemove[watchTask7p1]]];
watchTask7p1 = SessionSubmit[ScheduledTask[watchCheck7p1[], 10]];
{"watching " <> watchReq7p1, "last handled token: " <> watchLast7p1}
