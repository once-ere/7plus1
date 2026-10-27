(* tools/eval_cells_log.wl -- development check of a generated notebook: evaluate every input cell of the notebook named in
   <repository>/build/dev/target.txt (default build/dev/7+1_dev.nb) in a fresh context, expression by expression, with a time
   limit of 300 s per expression, and write a log with the timing, the messages and the result of every check["...", ...] to
   <repository>/build/dev/eval_cells_log.txt.  It stops at the first expression that fails (a syntax error, a message, a check
   that is not True, or a time-out), because later cells depend on earlier ones.  Run it with
       Get["<repository>/tools/eval_cells_log.wl"]
   or through tools/watch_requests.wl with the request line  "<token> run tools/eval_cells_log.wl". *)
Module[{root, tfile, target, nbx, codes, log, put, trunc600, t0, ctx, ncell = 0, nexpr = 0, nmsg = 0, nfail = 0, lim = 300, logfile, stop = False},
  root = ParentDirectory[DirectoryName[$InputFileName]];
  tfile = FileNameJoin[{root, "build", "dev", "target.txt"}];
  target = If[FileExistsQ[tfile], FileNameJoin[{root, StringTrim[ReadString[tfile]]}], FileNameJoin[{root, "build", "dev", "7+1_dev.nb"}]];
  logfile = FileNameJoin[{root, "build", "dev", "eval_cells_log.txt"}];
  log = {"target: " <> target, "version: " <> $Version, "started: " <> DateString[]};
  put[s_String] := (AppendTo[log, s]; Export[logfile, StringRiffle[log, "\n"], "Text"]);
  trunc600[s_String] := StringTake[s, Min[StringLength[s], 600]];
  nbx = Get[target];
  codes = Cases[nbx, Cell[c_String, "Input", ___] :> c, Infinity];
  put["input cells: " <> ToString[Length[codes]]];
  ctx = "dev7p1x" <> ToString[RandomInteger[10^9]] <> "`";
  t0 = AbsoluteTime[];
  Block[{$Context = ctx, $ContextPath = {ctx, "System`"}},
   Do[If[stop, Break[]]; ncell++;
    Module[{held, k = 0},
     held = Quiet[ReadList[StringToStream[code], Hold[Expression]], {General::shdw}];
     If[! ListQ[held] || MemberQ[held, $Failed], put["cell " <> ToString[ncell] <> ": SYNTAX ERROR in " <> ToString[StringTake[code, Min[300, StringLength[code]]], InputForm]]; nfail++; stop = True,
      Do[k++; nexpr++;
       Module[{r, t, msgs, bad = False, mlist = {}, mtext},
        {t, {r, msgs}} = AbsoluteTiming[Block[{$MessageList = {}},
            {Internal`HandlerBlock[{"Message", If[Last[#] === True, AppendTo[mlist, trunc600[ToString[#, InputForm]]]] &},
               Quiet[TimeConstrained[ReleaseHold[h], lim, "TIMEOUT"], {General::shdw}]], $MessageList}]];
        mtext = StringRiffle[mlist, "\n            "];
        If[msgs =!= {}, nmsg += Length[msgs]; bad = True];
        put[StringJoin["cell ", ToString[ncell], ".", ToString[k], "  ", ToString[NumberForm[t, {7, 2}]], " s",
          If[msgs =!= {}, "  MESSAGES: " <> ToString[msgs, InputForm] <> "\n      TEXT: " <> StringTake[mtext, Min[StringLength[mtext], 2000]], ""],
          Which[r === "TIMEOUT", (bad = True; "  TIMEOUT"),
           MatchQ[r, _String -> (True | False)], (If[r[[2]] =!= True, bad = True]; "  " <> ToString[r, InputForm]),
           r === Null, "",
           True, "  -> " <> StringTake[#, Min[StringLength[#], 300]] &[ToString[Short[r, 2], InputForm]]]]];
        If[bad, nfail++; stop = True; put["STOPPED at cell " <> ToString[ncell] <> "." <> ToString[k] <> ": " <> ToString[Short[h, 4], InputForm]]; Break[]]],
       {h, held}]]],
    {code, codes}]];
  put["finished: " <> DateString[] <> "  total seconds: " <> ToString[N[AbsoluteTime[] - t0]]];
  Export[FileNameJoin[{root, "build", "dev", "defined_symbols.txt"}],
    StringRiffle[Sort[StringDrop[#, StringLength[ctx]] & /@ Select[Names[ctx <> "*"],
       ToExpression[#, InputForm, Function[sym, OwnValues[sym] =!= {} || DownValues[sym] =!= {} || SubValues[sym] =!= {}, HoldAll]] &]], "\n"], "Text"];
  put["cells: " <> ToString[ncell] <> "  expressions: " <> ToString[nexpr] <> "  messages: " <> ToString[nmsg] <> "  failures: " <> ToString[nfail] <>
    If[nfail == 0 && nmsg == 0, "  ALL GOOD", ""]];
  {ncell, nexpr, nmsg, nfail}]
