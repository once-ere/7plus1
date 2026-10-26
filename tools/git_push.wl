(* tools/git_push.wl -- push the local repository (the folder that contains tools/) to its remote "origin"
   with the user's own Git for Windows and the credentials it already has stored.

   Usage, from any notebook of a running Wolfram front end:
       Get["D:/Developer/github/7+1/tools/git_push.wl"]
   Nothing is typed or prompted: git runs with GIT_TERMINAL_PROMPT=0 and credential.interactive=never, so a
   missing login makes the push fail with a message instead of opening a dialog.  The output of every git
   command is written to <repository>/build/git_push_log.txt (build/ is not part of the repository). *)
Module[{root, rootF, git, env, run, results, logText},
  root = ParentDirectory[DirectoryName[$InputFileName]];
  rootF = StringReplace[root, "\\" -> "/"];
  git = SelectFirst[{"git", "C:\\Program Files\\Git\\cmd\\git.exe", "C:\\Program Files\\Git\\bin\\git.exe"},
     Quiet[Check[RunProcess[{#, "--version"}, "ExitCode"], 1]] === 0 &, None];
  If[git === None, Return["git was not found on this computer"]];
  env = Join[Association[GetEnvironment[]], <|"GIT_TERMINAL_PROMPT" -> "0", "GCM_INTERACTIVE" -> "Never"|>];
  run[args_List] := Module[{r = RunProcess[Join[{git, "-c", "safe.directory=" <> rootF, "-c", "credential.interactive=never",
        "-C", rootF}, args], ProcessEnvironment -> env]},
     {StringRiffle[args, " "], r["ExitCode"], StringTrim[r["StandardOutput"]], StringTrim[r["StandardError"]]}];
  results = {run[{"--version"}], run[{"remote", "-v"}], run[{"log", "--oneline", "-5"}], run[{"push", "-u", "origin", "main"}],
     run[{"ls-remote", "origin"}]};
  logText = StringRiffle[Table[StringJoin["$ git ", r[[1]], "\nexit code: ", ToString[r[[2]]], "\n", r[[3]], "\n", r[[4]], "\n"], {r, results}], "\n"];
  Quiet[CreateDirectory[FileNameJoin[{root, "build"}]]];
  Export[FileNameJoin[{root, "build", "git_push_log.txt"}], StringJoin["finished: ", DateString[], "\n\n", logText], "Text"];
  Column[Table[Row[{Style["git " <> r[[1]], Bold], "   exit code ", r[[2]]}], {r, results}]]]
