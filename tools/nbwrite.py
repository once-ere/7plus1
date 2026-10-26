"""Minimal writer for Wolfram notebooks (.nb).

Cells are given as (style, content) pairs.  Input cells hold Wolfram Language
source text; text cells hold plain text or a list of (text, style-options)
runs.  The writer escapes every non-ASCII character as a Wolfram named
character (\\[Name]) or \\:xxxx so that the notebook file is pure ASCII, and it
nests the cells into CellGroupData by heading level, exactly as the front end
does when a notebook is typed by hand.
"""
import re


NAMED = {
 'α':'Alpha','β':'Beta','γ':'Gamma','δ':'Delta','ε':'Epsilon','ζ':'Zeta','η':'Eta','θ':'Theta',
 'ι':'Iota','κ':'Kappa','λ':'Lambda','μ':'Mu','ν':'Nu','ξ':'Xi','ο':'Omicron','π':'Pi','ρ':'Rho',
 'σ':'Sigma','τ':'Tau','υ':'Upsilon','φ':'Phi','ϕ':'CurlyPhi','χ':'Chi','ψ':'Psi','ω':'Omega',
 'Α':'CapitalAlpha','Β':'CapitalBeta','Γ':'CapitalGamma','Δ':'CapitalDelta','Ε':'CapitalEpsilon',
 'Ζ':'CapitalZeta','Η':'CapitalEta','Θ':'CapitalTheta','Ι':'CapitalIota','Κ':'CapitalKappa',
 'Λ':'CapitalLambda','Μ':'CapitalMu','Ν':'CapitalNu','Ξ':'CapitalXi','Ο':'CapitalOmicron',
 'Π':'CapitalPi','Ρ':'CapitalRho','Σ':'CapitalSigma','Τ':'CapitalTau','Υ':'CapitalUpsilon',
 'Φ':'CapitalPhi','Χ':'CapitalChi','Ψ':'CapitalPsi','Ω':'CapitalOmega',
 '∂':'PartialD','ℓ':'ScriptL','→':'RightArrow','←':'LeftArrow','↔':'LeftRightArrow',
 '⇒':'DoubleRightArrow','⇔':'DoubleLeftRightArrow','⊗':'CircleTimes','⊕':'CirclePlus',
 '·':'CenterDot','×':'Times','≅':'TildeFullEqual','≠':'NotEqual','≤':'LessEqual','≥':'GreaterEqual',
 '∈':'Element','∉':'NotElement','ℝ':'DoubleStruckCapitalR','ℂ':'DoubleStruckCapitalC',
 'ℍ':'DoubleStruckCapitalH','𝕆':'DoubleStruckCapitalO','∘':'SmallCircle','√':'Sqrt',
 '±':'PlusMinus','∓':'MinusPlus','—':'LongDash','–':'Dash','…':'Ellipsis','⋯':'CenterEllipsis',
 '“':'OpenCurlyDoubleQuote','”':'CloseCurlyDoubleQuote','‘':'OpenCurlyQuote','’':'CloseCurlyQuote',
 '∑':'Sum','∏':'Product','∞':'Infinity','≡':'Congruent','≈':'TildeTilde','∼':'Tilde',
 '⟨':'LeftAngleBracket','⟩':'RightAngleBracket','†':'Dagger','∧':'Wedge','∨':'Vee','¬':'Not',
 '∀':'ForAll','∃':'Exists','∅':'EmptySet','✓':'Checkmark','⚠':'WarningSign','•':'Bullet',
 'ᵀ':'Transpose','§':'Section','°':'Degree','′':'Prime','∖':'Backslash',
}

def esc_char(ch):
    if ch in NAMED: return '\\[' + NAMED[ch] + ']'
    o = ord(ch)
    if o < 0x10000: return '\\:%04x' % o
    return '\\|%06x' % o

_ESC = re.compile(r'\\(\[[A-Za-z]+\]|:[0-9a-fA-F]{4})')

def wl_string(s):
    """Return the Wolfram Language string literal (with quotes) for python str s.
    Named-character escapes already present in s (\\[Name], \\:xxxx) are kept as single escapes,
    so that the front end shows the character itself."""
    out = ['"']
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == '\\':
            m = _ESC.match(s, i)
            if m:
                out.append(m.group(0)); i = m.end(); continue
            out.append('\\\\'); i += 1; continue
        i += 1
        if ch == '"': out.append('\\"')
        elif ch == '\n': out.append('\\n')
        elif ch == '\t': out.append('    ')
        elif ord(ch) < 128: out.append(ch)
        else: out.append(esc_char(ch))
    out.append('"')
    return ''.join(out)

LEVEL = {'Title':0, 'Subtitle':1, 'Chapter':1, 'Section':2, 'Subsection':3, 'Subsubsection':4}

# Characters that the text font of the Wolfram front end shows badly: the non-BMP double-struck O
# (the next character is drawn on top of it), the script l (it looks like a slash, so x = a + l b would
# read like a division) and the combining arrow above (missing glyph).  In text cells the double-struck O
# and the arrow are written as small inline formulas, which use the math fonts,
# Cell[BoxData[FormBox[..., TraditionalForm]]], and the script l in the font Cambria Math.
SPECIAL_RUN = re.compile('(\U0001D546|\u2113|[A-Za-z]\u20d7)')

def inline_formula(tok):
    if tok == '\U0001D546':
        box = '"\\[DoubleStruckCapitalO]"'
    else:                                   # letter + combining arrow above: OverVector
        box = 'OverscriptBox["%s", "\\[RightVector]"]' % tok[0]
    return 'Cell[BoxData[FormBox[%s, TraditionalForm]]]' % box

def text_runs(s, opts=None):
    """WL TextData elements for the plain text s (optionally styled with opts)."""
    o = ', '.join('%s->%s' % (k, v) for k, v in (opts or {}).items())
    parts = []
    for tok in SPECIAL_RUN.split(s):
        if not tok: continue
        if tok == '\u2113':
            parts.append('StyleBox["\\[ScriptL]", FontFamily->"Cambria Math"%s]' % (', ' + o if o else ''))
        elif SPECIAL_RUN.fullmatch(tok): parts.append(inline_formula(tok))
        elif o: parts.append('StyleBox[%s, %s]' % (wl_string(tok), o))
        else: parts.append(wl_string(tok))
    return parts

def text_content(t):
    """t: str  or list of runs: str | (str, dict_of_options)"""
    if isinstance(t, str):
        if not SPECIAL_RUN.search(t):
            return wl_string(t)
        return 'TextData[{' + ', '.join(text_runs(t)) + '}]'
    parts = []
    for r in t:
        if isinstance(r, str):
            parts.extend(text_runs(r))
        else:
            s, opts = r
            parts.extend(text_runs(s, opts))
    return 'TextData[{' + ', '.join(parts) + '}]'

def cell_expr(style, content, opts=None):
    o = ''
    if opts:
        o = ', ' + ', '.join('%s->%s' % (k, v) for k, v in opts.items())
    if style in ('Input', 'Code'):
        return 'Cell[%s, "%s"%s]' % (wl_string(content), style, o)
    return 'Cell[%s, "%s"%s]' % (text_content(content), style, o)

def group(cells):
    """cells: list of (style, content, opts). Build nested CellGroupData by heading level."""
    def build(i, level):
        out = []
        while i < len(cells):
            st = cells[i][0]
            lv = LEVEL.get(st, 99)
            if lv <= level:
                break
            if lv < 99:
                head = cell_expr(*cells[i])
                children, j = build(i + 1, lv)
                out.append('Cell[CellGroupData[{\n' + ',\n'.join([head] + children) + '\n}, Open  ]]')
                i = j
            else:
                out.append(cell_expr(*cells[i]))
                i += 1
        return out, i
    out, _ = build(0, -1)
    return out

def notebook(cells, nb_options=None):
    body = ',\n\n'.join(group(cells))
    o = nb_options or {}
    opts = ',\n'.join('%s->%s' % (k, v) for k, v in o.items())
    head = ('(* Content-type: application/vnd.wolfram.mathematica *)\n\n'
            '(*** Wolfram Notebook File ***)\n(* http://www.wolfram.com/nb *)\n\n')
    return head + 'Notebook[{\n' + body + '\n}' + (',\n' + opts if opts else '') + '\n]\n'
