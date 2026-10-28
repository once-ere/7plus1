#!/usr/bin/env python3
"""Build the documentation of 7+1.nb.

    python3 tools/make_docs.py [--reuse-runs | --reuse-refs] [NAME ...]

writes, in the repository root, for each NAME of the list below the files NAME.md, NAME.tex and NAME.pdf:
    README               README_PROVENANCE               README_SYMBOLS                   (the whole notebook)
    DIRAC16COMPLEX       DIRAC16COMPLEX_PROVENANCE       DIRAC16COMPLEX_SYMBOLS           (sections 7-10)
    GRAVITY              GRAVITY_PROVENANCE              GRAVITY_SYMBOLS                  (sections 11-14)
    KOHN_SHAM_GRAVITY    KOHN_SHAM_GRAVITY_PROVENANCE    KOHN_SHAM_GRAVITY_SYMBOLS        (sections 15-16)
    STANDARD_MODEL       STANDARD_MODEL_PROVENANCE       STANDARD_MODEL_SYMBOLS           (sections 17-19)
    COUPLED_SOLUTIONS    COUPLED_SOLUTIONS_PROVENANCE    COUPLED_SOLUTIONS_SYMBOLS        (sections 20-23)
    DARK_SECTOR          DARK_SECTOR_PROVENANCE          DARK_SECTOR_SYMBOLS              (section 24)
    AUXILIARY_TORUS      AUXILIARY_TORUS_PROVENANCE      AUXILIARY_TORUS_SYMBOLS          (the auxiliary torus; section 25 and 14-24)
    WORLD_BRANE          WORLD_BRANE_PROVENANCE          WORLD_BRANE_SYMBOLS              (the world brane; sections 25-26 and 17-24)
    RADION               RADION_PROVENANCE               RADION_SYMBOLS                   (the radion on the 4-torus and the 4-sphere; section 26)
    7plus1_status_2026-10-28   7plus1_status_2026-10-28_PROVENANCE   7plus1_status_2026-10-28_SYMBOLS
                         (the state of the repository at commit 102f399: the contents of the spreadsheet 7plus1_status_2026-10-28.xlsx,
                          which tools/make_status_xlsx.py writes from tools/status_2026-10-28.json)
Every document is self-contained.  The LaTeX files are compiled with latexmk -xelatex inside build/docs (excluded by .gitignore);
only the PDFs are copied back.  Every matrix, table, check list and number is computed by the doc_*.py modules or parsed from the
evaluated 7+1.nb and from the outputs of the independent Python programs, which are run here (or, with --reuse-runs, taken from
build/docs/runs.json of the previous run; with --reuse-refs only tools/verify_7p1_nb.py, which reads the notebook, is run again).
"""
import os, sys, subprocess, shutil, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
from doc_dirac16_common import gather
from doc_stage3_common import gather3
import doc_readme, doc_symbols, doc_provenance
import doc_dirac16, doc_dirac16_symbols, doc_dirac16_provenance
import doc_gravity, doc_gravity_symbols, doc_gravity_provenance
import doc_ks_gravity, doc_ks_gravity_symbols, doc_ks_gravity_provenance
import doc_sm, doc_sm_symbols, doc_sm_provenance
import doc_coupled, doc_coupled_symbols, doc_coupled_provenance
import doc_dark, doc_dark_symbols, doc_dark_provenance
import doc_torus, doc_torus_symbols, doc_torus_provenance
import doc_brane, doc_brane_symbols, doc_brane_provenance
import doc_radion, doc_radion_symbols, doc_radion_provenance
import doc_status, doc_status_symbols, doc_status_provenance
from doc_status_common import NAME as STATUS, present as status_present

PROGRAMS = [('verify_7p1_nb.py', 'parses the evaluated 7+1.nb'), ('ref44.py', 'rebuilds 4+4.nb'), ('ref71.py', 'rebuilds 7+1.nb objects'),
            ('extra71.py', 'further identities'), ('ref_dirac16.py', 'independent checks of sections 7-9'),
            ('ref_gravity.py', 'independent checks of sections 11-14'), ('ref_ks_gravity.py', 'independent checks of sections 15-16'),
            ('ref_sm.py', 'independent checks of sections 17-19'), ('ref_coupled.py', 'independent checks of sections 20-23'),
            ('ref_dark.py', 'independent checks of section 24'), ('ref_torus_brane.py', 'independent checks of section 25'),
            ('ref_radion.py', 'independent checks of section 26')]


def run(cmd):
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return (p.stdout + p.stderr).strip(), p.returncode


def all_runs(rfile):
    """run the verification programs; their outputs go into the provenance documents.  The independent references run first
    and are saved to rfile, because tools/verify_7p1_nb.py compares the notebook with their outputs."""
    runs = []
    for f, what in PROGRAMS[1:]:
        out, rc = run([sys.executable, '-B', 'tools/' + f])
        runs.append(['tools/%s (%s; exit code %d)' % (f, what, rc), out, rc])
    json.dump(runs, open(rfile, 'w', encoding='utf-8'))
    f, what = PROGRAMS[0]
    out, rc = run([sys.executable, '-B', 'tools/' + f])
    runs.insert(0, ['tools/%s (%s; exit code %d)' % (f, what, rc), out, rc])
    return runs


def main(argv):
    reuse = '--reuse-runs' in argv
    reuse_refs = '--reuse-refs' in argv
    only = [a for a in argv if not a.startswith('--')]
    bdir = os.path.join(ROOT, 'build', 'docs'); os.makedirs(bdir, exist_ok=True)
    rfile = os.path.join(bdir, 'runs.json')
    if reuse and os.path.exists(rfile):
        runs = json.load(open(rfile, encoding='utf-8'))
    elif reuse_refs and os.path.exists(rfile):
        # the independent programs do not read the notebook: keep their saved outputs (running those that have none yet), and run
        # verify_7p1_nb.py again
        runs = [r for r in json.load(open(rfile, encoding='utf-8')) if not r[0].startswith('tools/' + PROGRAMS[0][0])]
        have = set(r[0].split(' ')[0] for r in runs)
        for f, what in PROGRAMS[1:]:
            if 'tools/' + f not in have:
                out, rc = run([sys.executable, '-B', 'tools/' + f])
                runs.append(['tools/%s (%s; exit code %d)' % (f, what, rc), out, rc])
        json.dump(runs, open(rfile, 'w', encoding='utf-8'))
        f, what = PROGRAMS[0]
        out, rc = run([sys.executable, '-B', 'tools/' + f])
        runs.insert(0, ['tools/%s (%s; exit code %d)' % (f, what, rc), out, rc])
        json.dump(runs, open(rfile, 'w', encoding='utf-8'))
    else:
        runs = all_runs(rfile)
        json.dump(runs, open(rfile, 'w', encoding='utf-8'))
    D = gather()
    D3 = gather3(runs)
    D['s3'] = D3
    today = D['today']
    byname = {r[0].split(' ')[0]: r for r in runs}
    rc = {k: v[2] for k, v in byname.items()}
    runs_readme = [(t, o) for t, o, c in runs]
    runs16 = [(byname['tools/verify_7p1_nb.py'][0], byname['tools/verify_7p1_nb.py'][1]),
              (byname['tools/ref_dirac16.py'][0], byname['tools/ref_dirac16.py'][1])]
    builders = {'README': lambda: doc_readme.build(D),
                'README_SYMBOLS': lambda: doc_symbols.build(D),
                'README_PROVENANCE': lambda: doc_provenance.build(D, runs_readme),
                'DIRAC16COMPLEX': lambda: doc_dirac16.build(D),
                'DIRAC16COMPLEX_SYMBOLS': lambda: doc_dirac16_symbols.build(D),
                'DIRAC16COMPLEX_PROVENANCE': lambda: doc_dirac16_provenance.build(D, runs16),
                'GRAVITY': lambda: doc_gravity.build(D3, today),
                'GRAVITY_SYMBOLS': lambda: doc_gravity_symbols.build(D3, today),
                'GRAVITY_PROVENANCE': lambda: doc_gravity_provenance.build(D3, today, runs),
                'KOHN_SHAM_GRAVITY': lambda: doc_ks_gravity.build(D3, today),
                'KOHN_SHAM_GRAVITY_SYMBOLS': lambda: doc_ks_gravity_symbols.build(D3, today),
                'KOHN_SHAM_GRAVITY_PROVENANCE': lambda: doc_ks_gravity_provenance.build(D3, today, runs),
                'STANDARD_MODEL': lambda: doc_sm.build(D3, today),
                'STANDARD_MODEL_SYMBOLS': lambda: doc_sm_symbols.build(D3, today),
                'STANDARD_MODEL_PROVENANCE': lambda: doc_sm_provenance.build(D3, today, runs),
                'COUPLED_SOLUTIONS': lambda: doc_coupled.build(D3, today),
                'COUPLED_SOLUTIONS_SYMBOLS': lambda: doc_coupled_symbols.build(D3, today),
                'COUPLED_SOLUTIONS_PROVENANCE': lambda: doc_coupled_provenance.build(D3, today, runs),
                'DARK_SECTOR': lambda: doc_dark.build(D3, today),
                'DARK_SECTOR_SYMBOLS': lambda: doc_dark_symbols.build(D3, today),
                'DARK_SECTOR_PROVENANCE': lambda: doc_dark_provenance.build(D3, today, runs),
                'AUXILIARY_TORUS': lambda: doc_torus.build(D3, today),
                'AUXILIARY_TORUS_SYMBOLS': lambda: doc_torus_symbols.build(D3, today),
                'AUXILIARY_TORUS_PROVENANCE': lambda: doc_torus_provenance.build(D3, today, runs),
                'WORLD_BRANE': lambda: doc_brane.build(D3, today),
                'WORLD_BRANE_SYMBOLS': lambda: doc_brane_symbols.build(D3, today),
                'WORLD_BRANE_PROVENANCE': lambda: doc_brane_provenance.build(D3, today, runs),
                'RADION': lambda: doc_radion.build(D3, today),
                'RADION_SYMBOLS': lambda: doc_radion_symbols.build(D3, today),
                'RADION_PROVENANCE': lambda: doc_radion_provenance.build(D3, today, runs)}
    if status_present():
        builders.update({STATUS: lambda: doc_status.build(today),
                         STATUS + '_SYMBOLS': lambda: doc_status_symbols.build(today),
                         STATUS + '_PROVENANCE': lambda: doc_status_provenance.build(D, today, runs)})
    docs = {k: f() for k, f in builders.items() if not only or k in only}
    report = []
    for name, d in docs.items():
        md, tex = d.to_md(), d.to_tex()
        with open(os.path.join(ROOT, name + '.md'), 'w', encoding='utf-8') as f: f.write(md)
        with open(os.path.join(ROOT, name + '.tex'), 'w', encoding='utf-8') as f: f.write(tex)
        shutil.copyfile(os.path.join(ROOT, name + '.tex'), os.path.join(bdir, name + '.tex'))
        p = subprocess.run(['latexmk', '-xelatex', '-interaction=nonstopmode', '-halt-on-error', name + '.tex'],
                           cwd=bdir, capture_output=True, text=True)
        pdf = os.path.join(bdir, name + '.pdf')
        ok = p.returncode == 0 and os.path.exists(pdf)
        if ok: shutil.copyfile(pdf, os.path.join(ROOT, name + '.pdf'))
        logf = os.path.join(bdir, name + '.log')
        log = open(logf, encoding='utf-8', errors='replace').read() if os.path.exists(logf) else ''
        report.append((name, len(md), len(tex), p.returncode, ok, log.count('Overfull \\hbox'), log.count('Underfull \\hbox'),
                       log.count('\n! '), log.count('invalid in math mode')))
    for r in report:
        print('%-30s md %6d chars  tex %6d chars  latexmk exit %d  pdf %s  overfull %d  underfull %d  errors %d  math-mode warnings %d' % r)
    print('program exit codes:', ', '.join('%s %d' % (k.split('/')[-1], v) for k, v in rc.items()))
    if D3.get('xact_tex_errors'):
        print('xAct LaTeX conversion problems:', D3['xact_tex_errors'])
    print('checks of sections 11-26 (grav, ks, sm, cp, dk, tb, rd):', {k: (len(v), sum(1 for c in v if c[1] == 'True')) for k, v in D3['checks'].items()})
    print('all checks of the notebook:', D3.get('n_checks_total'), 'all passed:', D3.get('all_passed'), 'message cells:', len(D3['messages']))
    return all(r[4] for r in report) and all(v == 0 for v in rc.values())


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1:]) else 1)
