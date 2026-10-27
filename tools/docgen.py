"""Tiny document model that renders the same content to GitHub Markdown and to XeLaTeX.

Inline markup inside paragraph strings:
    `code`      monospace (Mathematica names, may contain Unicode)
    **bold**    bold
    *italic*    italic
    $math$      inline mathematics (LaTeX; identical in both outputs)
Display mathematics is given as LaTeX source and rendered as a ```math block (Markdown)
or as an equation*/align* environment (LaTeX).

Unicode in the LaTeX output.  The text font (Latin Modern Roman) has no Greek letters and few
mathematical symbols, and the monospace font (DejaVu Sans Mono) lacks a few symbols such as the
double-struck O.  Every character that the font in use cannot show is therefore written as the
equivalent LaTeX mathematics (UNI_MATH below) or, failing that, in DejaVu Serif.  The coverage of the
fonts is read from the font files themselves with fontTools.
"""
import re, os, subprocess

INLINE = re.compile(r'(`[^`]+`|\*\*[^*]+\*\*|\$[^$]+\$|\*[^*]+\*)')

# ----------------------------------------------------------------------------- Unicode -> LaTeX math
UNI_MATH = {
    'α': r'\alpha', 'β': r'\beta', 'γ': r'\gamma', 'δ': r'\delta', 'ε': r'\varepsilon', 'ϵ': r'\epsilon',
    'ζ': r'\zeta', 'η': r'\eta', 'θ': r'\theta', 'ϑ': r'\vartheta', 'ι': r'\iota', 'κ': r'\kappa',
    'λ': r'\lambda', 'μ': r'\mu', 'µ': r'\mu', 'ν': r'\nu', 'ξ': r'\xi', 'ο': r'o', 'π': r'\pi',
    'ρ': r'\rho', 'σ': r'\sigma', 'ς': r'\varsigma', 'τ': r'\tau', 'υ': r'\upsilon', 'φ': r'\varphi',
    'ϕ': r'\phi', 'χ': r'\chi', 'ψ': r'\psi', 'ω': r'\omega',
    'Α': r'\mathrm{A}', 'Β': r'\mathrm{B}', 'Γ': r'\Gamma', 'Δ': r'\Delta', 'Ε': r'\mathrm{E}',
    'Ζ': r'\mathrm{Z}', 'Η': r'\mathrm{H}', 'Θ': r'\Theta', 'Ι': r'\mathrm{I}', 'Κ': r'\mathrm{K}',
    'Λ': r'\Lambda', 'Μ': r'\mathrm{M}', 'Ν': r'\mathrm{N}', 'Ξ': r'\Xi', 'Ο': r'\mathrm{O}',
    'Π': r'\Pi', 'Ρ': r'\mathrm{P}', 'Σ': r'\Sigma', 'Τ': r'\mathrm{T}', 'Υ': r'\Upsilon',
    'Φ': r'\Phi', 'Χ': r'\mathrm{X}', 'Ψ': r'\Psi', 'Ω': r'\Omega',
    'ℝ': r'\mathbb{R}', 'ℂ': r'\mathbb{C}', 'ℍ': r'\mathbb{H}', '𝕆': r'\mathbb{O}', 'ℤ': r'\mathbb{Z}',
    'ℕ': r'\mathbb{N}', 'ℚ': r'\mathbb{Q}', '𝕌': r'\mathbb{U}',
    '⊕': r'\oplus', '⊗': r'\otimes', '→': r'\to', '←': r'\leftarrow', '↔': r'\leftrightarrow',
    '⇒': r'\Rightarrow', '⇔': r'\Leftrightarrow', '↦': r'\mapsto', '≅': r'\cong', '≠': r'\neq',
    '≤': r'\leq', '≥': r'\geq', '∈': r'\in', '∉': r'\notin', '≡': r'\equiv', '≈': r'\approx',
    '∼': r'\sim', '≃': r'\simeq', '∘': r'\circ', '√': r'\surd', '±': r'\pm', '∓': r'\mp',
    '⋯': r'\cdots', '∑': r'\sum', '∏': r'\prod', '∞': r'\infty', '⟨': r'\langle', '⟩': r'\rangle',
    '†': r'\dagger', '∧': r'\wedge', '∨': r'\vee', '¬': r'\neg', '∀': r'\forall', '∃': r'\exists',
    '∅': r'\emptyset', '∂': r'\partial', 'ℓ': r'\ell', '▸': r'\triangleright', '·': r'\cdot',
    '✓': r'\checkmark', 'ᵀ': r'{}^{\mathsf T}', 'ᴴ': r'{}^{\dagger}', '×': r'\times', '−': r'-',
    '′': r"{}'", '″': r"{}''", '∖': r'\setminus', '⊂': r'\subset', '⊆': r'\subseteq', '∩': r'\cap',
    '∪': r'\cup', 'ℏ': r'\hbar', '∇': r'\nabla', '⊥': r'\perp', '∥': r'\parallel', '□': r'\square',
    '■': r'\blacksquare', '•': r'\bullet', '°': r'{}^{\circ}', '≔': r'\coloneqq',
    '∫': r'\int', '∝': r'\propto', '⋮': r'\vdots', '↪': r'\hookrightarrow', '↳': r'\hookrightarrow',
    '⁰': r'{}^{0}', '¹': r'{}^{1}', '²': r'{}^{2}', '³': r'{}^{3}', '⁴': r'{}^{4}', '⁵': r'{}^{5}',
    '⁶': r'{}^{6}', '⁷': r'{}^{7}', '⁸': r'{}^{8}', '⁹': r'{}^{9}', '⁻': r'{}^{-}', '⁺': r'{}^{+}',
    '₀': r'{}_{0}', '₁': r'{}_{1}', '₂': r'{}_{2}', '₃': r'{}_{3}', '₄': r'{}_{4}', '₅': r'{}_{5}',
    '₆': r'{}_{6}', '₇': r'{}_{7}', '₈': r'{}_{8}', '₉': r'{}_{9}',
}
GREEK = set('αβγδεϵζηθϑικλμµνξπρσςτυφϕχψωΓΔΘΛΞΠΣΥΦΨΩ')
# combining accents placed on the preceding character
COMBINING = {'̄': r'\bar', '̅': r'\overline', '̂': r'\hat', '̃': r'\tilde',
             '̇': r'\dot', '̈': r'\ddot', '⃗': r'\vec', '́': r'\acute', '̀': r'\grave'}

def _font_file(query, kpse=None):
    try:
        if kpse:
            p = subprocess.run(['kpsewhich', kpse], capture_output=True, text=True).stdout.strip()
            if p: return p
        return subprocess.run(['fc-match', '-f', '%{file}', query], capture_output=True, text=True).stdout.strip()
    except Exception:
        return ''

def _coverage(path):
    try:
        from fontTools.ttLib import TTFont
        return set(TTFont(path, lazy=True).getBestCmap().keys())
    except Exception:
        return None

COVER_TEXT = _coverage(_font_file('Latin Modern Roman', 'lmroman10-regular.otf'))
COVER_MONO = _coverage(_font_file('DejaVu Sans Mono'))
COVER_UNI = _coverage(_font_file('DejaVu Serif'))

def _covered(ch, cover):
    if ord(ch) < 128: return True
    return cover is not None and ord(ch) in cover

def _math(latex, original):
    """Typeset latex as mathematics; PDF strings (bookmarks) get the original Unicode text instead."""
    return r'\texorpdfstring{\ensuremath{' + latex + '}}{' + original + '}'

def _uni(ch, mode):
    """LaTeX for one non-ASCII character in mode 'text' (roman) or 'mono' (\\texttt, Verbatim).
    In running text Greek letters are always set as mathematics (as in the displayed formulas)."""
    cover = COVER_TEXT if mode == 'text' else COVER_MONO
    if mode == 'text' and ch in GREEK: return _math(UNI_MATH[ch], ch)
    if _covered(ch, cover): return ch
    if ch in UNI_MATH: return _math(UNI_MATH[ch], ch)
    if COVER_UNI is not None and ord(ch) in COVER_UNI: return r'{\unifont ' + ch + '}'
    return ch

SPECIAL = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '#': r'\#', '_': r'\_', '$': r'\$',
           '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}

def tex_escape_text(s, mode='text'):
    """Escape plain text for LaTeX; mode 'text' (roman) or 'mono' (\\texttt)."""
    out = []; i = 0; n = len(s)
    cover = COVER_TEXT if mode == 'text' else COVER_MONO
    while i < n:
        ch = s[i]
        nxt = s[i + 1] if i + 1 < n else ''
        if nxt in COMBINING and not (mode == 'mono' and _covered(ch, cover) and _covered(nxt, cover)):
            base = SPECIAL.get(ch, ch) if ord(ch) < 128 else UNI_MATH.get(ch, ch)
            out.append(_math(COMBINING[nxt] + '{' + base + '}', ch + nxt if ord(ch) >= 128 or ch.isalnum() else ch)); i += 2; continue
        if ch in SPECIAL: out.append(SPECIAL[ch])
        elif ord(ch) < 128: out.append(ch)
        else: out.append(_uni(ch, mode))
        i += 1
    return ''.join(out)

# plain-text (Unicode) rendering of inline mathematics, for PDF bookmarks of headings
PDF_MATH = [(r'\mathbb O', '\U0001D546'), (r'\mathbb H', '\u210d'), (r'\mathbb R', '\u211d'), (r'\mathbb C', '\u2102'),
            (r'\mathbb{O}', '\U0001D546'), (r'\mathbb{H}', '\u210d'), (r'\mathbb{R}', '\u211d'), (r'\mathbb{C}', '\u2102'),
            (r'\bar\tau', '\u03c4\u0304'), (r'\Gamma', '\u0393'), (r'\tau', '\u03c4'), (r'\sigma', '\u03c3'),
            (r'\zeta', '\u03b6'), (r'\eta', '\u03b7'), (r'\ell', '\u2113'), (r'\oplus', '\u2295'), (r'\cong', '\u2245'),
            (r'\mathrm{Re}', 'Re'), (r'\widehat', ''), (r'\not ', ''), (r'\,', ' '), (r'\ ', ' ')]

def pdf_plain(math):
    t = math
    for a, b in PDF_MATH: t = t.replace(a, b)
    t = re.sub(r'\\[A-Za-z]+', '', t)
    t = t.replace('{', '').replace('}', '').replace('^', '').replace('_', '')
    return re.sub(r'\s+', ' ', t).strip()

SCI = re.compile(r'(?<![\w.\\])(-?\d+(?:\.\d+)?\\times10\^\{-?\d+\})')

def sci_math(s):
    r"""A number written as 1.2\times10^{-5} outside inline mathematics is set as mathematics ($1.2\times10^{-5}$)."""
    if '\\times10^{' not in s: return s
    parts = []
    for tok in INLINE.split(s):
        if tok and not (tok.startswith('`') or tok.startswith('$')):
            tok = SCI.sub(r'$\1$', tok)
        parts.append(tok)
    return ''.join(parts)

def top_breaks(latex):
    r"""True if the display latex contains a line break \\ outside braces and outside \begin...\end environments."""
    depth = env = 0; i = 0; n = len(latex)
    while i < n:
        for key, d in ((r'\begin{', 1), (r'\end{', -1)):
            if latex.startswith(key, i):
                env += d; i = latex.index('}', i) + 1
                break
        else:
            c = latex[i]
            if c == '\\':
                if latex.startswith('\\\\', i) and depth == 0 and env == 0: return True
                i += 2; continue
            if c == '{': depth += 1
            elif c == '}': depth -= 1
            i += 1
    return False

def tex_inline(s, pdf=False):
    """Inline markup -> LaTeX.  pdf=True (headings): inline mathematics gets a plain-text twin for the PDF bookmarks."""
    s = sci_math(s)
    out = []
    for tok in INLINE.split(s):
        if not tok: continue
        if tok.startswith('`') and tok.endswith('`') and len(tok) > 1:
            out.append(r'\texttt{' + tex_escape_text(tok[1:-1], 'mono') + '}')
        elif tok.startswith('**') and tok.endswith('**') and len(tok) > 3:
            out.append(r'\textbf{' + tex_inline(tok[2:-2], pdf) + '}')
        elif tok.startswith('$') and tok.endswith('$') and len(tok) > 1:
            out.append(r'\texorpdfstring{' + tok + '}{' + pdf_plain(tok[1:-1]) + '}' if pdf else tok)
        elif tok.startswith('*') and tok.endswith('*') and len(tok) > 1:
            out.append(r'\emph{' + tex_inline(tok[1:-1], pdf) + '}')
        else:
            out.append(tex_escape_text(tok, 'text'))
    return ''.join(out)

def guard(t):
    """A line that follows \\\\ or \\item must not begin with '[' (LaTeX would read an optional argument)."""
    return '{}' + t if t.lstrip().startswith('[') else t

def md_inline(s):
    r"""Inline markup -> GitHub Markdown.  Inline mathematics $...$ is written $`...`$, GitHub's delimiter for
    expressions that contain characters with a meaning in Markdown (backslash escapes such as \, \{ or \\
    would otherwise be eaten before the mathematics is rendered)."""
    s = sci_math(s)
    out = []
    for tok in INLINE.split(s):
        if not tok: continue
        if tok.startswith('`') and tok.endswith('`') and len(tok) > 1:
            out.append(tok)
        elif tok.startswith('**') and tok.endswith('**') and len(tok) > 3:
            out.append('**' + md_inline(tok[2:-2]) + '**')
        elif tok.startswith('$') and tok.endswith('$') and len(tok) > 1:
            out.append('$`' + tok[1:-1] + '`$')
        elif tok.startswith('*') and tok.endswith('*') and len(tok) > 1:
            out.append('*' + md_inline(tok[1:-1]) + '*')
        else:
            out.append(tok)
    return ''.join(out)

def pmat(M, small=False):
    """LaTeX pmatrix of a sympy/list matrix; entries rendered with sympy.latex."""
    import sympy as sp
    M = sp.Matrix(M)
    rows = []
    for i in range(M.rows):
        rows.append(' & '.join(sp.latex(M[i, j]) for j in range(M.cols)))
    if small:
        return r'\left(\begin{smallmatrix}' + r' \\ '.join(rows) + r'\end{smallmatrix}\right)'
    return r'\begin{pmatrix}' + r' \\ '.join(rows) + r'\end{pmatrix}'

VERB_WIDTH = 96          # characters per line of a \small Verbatim block on A4 with 2.2 cm margins

def wrap_lines(text, width=VERB_WIDTH):
    """Wrap long lines; a continuation line starts with '    ...'."""
    lines = []
    for ln in text.rstrip('\n').split('\n'):
        ln = ln.replace('\t', '    ')
        first = True
        while len(ln) > (width if first else width - 7):
            w = width if first else width - 7
            cut = ln.rfind(' ', w // 2, w)
            if cut <= 0: cut = w
            lines.append(ln[:cut] if first else '    ...' + ln[:cut])
            ln = ln[cut:].lstrip(' '); first = False
        lines.append(ln if first else '    ...' + ln)
    return '\n'.join(lines)

def verbatim(text):
    """Verbatim block; long lines wrapped; characters the monospace font lacks are typeset as
    mathematics (fancyvrb's commandchars are then switched on and \\ { } escaped)."""
    body = wrap_lines(text)
    if all(_covered(ch, COVER_MONO) for ch in body):
        return [r'\begin{Verbatim}[fontsize=\small,frame=single,framesep=2mm]', body, r'\end{Verbatim}']
    esc = []
    for ch in body:
        if ch == '\\': esc.append(r'\textbackslash{}')
        elif ch == '{': esc.append(r'\{')
        elif ch == '}': esc.append(r'\}')
        elif ord(ch) < 128: esc.append(ch)
        else: esc.append(_uni(ch, 'mono'))
    return [r'\begin{Verbatim}[fontsize=\small,frame=single,framesep=2mm,commandchars=\\\{\}]', ''.join(esc), r'\end{Verbatim}']

class Doc:
    def __init__(self, title, subtitle='', meta=None, abstract=''):
        self.title = title; self.subtitle = subtitle; self.meta = meta or []; self.abstract = abstract
        self.blocks = []
    def h1(self, t): self.blocks.append(('h1', t))
    def h2(self, t): self.blocks.append(('h2', t))
    def h3(self, t): self.blocks.append(('h3', t))
    def p(self, t): self.blocks.append(('p', t))
    def eq(self, latex): self.blocks.append(('eq', latex))
    def align(self, lines): self.blocks.append(('align', lines))
    def ul(self, items): self.blocks.append(('ul', items))
    def ol(self, items): self.blocks.append(('ol', items))
    def code(self, text, lang='mathematica'): self.blocks.append(('code', text, lang))
    def table(self, header, rows, widths=None, small=False): self.blocks.append(('table', header, rows, widths, small))
    def quote(self, t): self.blocks.append(('quote', t))
    def rule(self): self.blocks.append(('rule',))

    # ------------------------------------------------------------------ Markdown
    def to_md(self):
        L = ['# ' + self.title, '']
        if self.subtitle: L += ['*' + self.subtitle + '*', '']
        for m in self.meta: L += [m + '  ']
        if self.meta: L.append('')
        if self.abstract: L += ['> ' + md_inline(self.abstract), '']
        for b in self.blocks:
            k = b[0]
            if k == 'h1': L += ['## ' + md_inline(b[1]), '']
            elif k == 'h2': L += ['### ' + md_inline(b[1]), '']
            elif k == 'h3': L += ['#### ' + md_inline(b[1]), '']
            elif k == 'p': L += [md_inline(b[1]), '']
            elif k == 'quote': L += ['> ' + md_inline(b[1]), '']
            elif k == 'eq':
                # a display with line breaks \\ is a gathered block (centred lines)
                L += ['```math'] + ([r'\begin{gathered}', b[1].strip(), r'\end{gathered}'] if top_breaks(b[1]) else [b[1].strip()]) + ['```', '']
            elif k == 'align':
                L += ['```math', r'\begin{aligned}'] + [ln + r' \\' if i < len(b[1]) - 1 else ln for i, ln in enumerate(b[1])] + [r'\end{aligned}', '```', '']
            elif k == 'ul': L += ['- ' + md_inline(it) for it in b[1]] + ['']
            elif k == 'ol': L += ['%d. ' % (i + 1) + md_inline(it) for i, it in enumerate(b[1])] + ['']
            elif k == 'code': L += ['```' + b[2], b[1].rstrip('\n'), '```', '']
            elif k == 'rule': L += ['---', '']
            elif k == 'table':
                header, rows = b[1], b[2]
                esc = lambda c: md_inline(str(c)).replace('|', r'\|')
                L.append('| ' + ' | '.join(esc(h) for h in header) + ' |')
                L.append('|' + '|'.join(['---'] * len(header)) + '|')
                for r in rows: L.append('| ' + ' | '.join(esc(c) for c in r) + ' |')
                L.append('')
        return '\n'.join(L) + '\n'

    # ------------------------------------------------------------------ LaTeX
    def to_tex(self):
        T = [r'\documentclass[10pt,a4paper]{article}',
             r'\usepackage[margin=2.2cm]{geometry}',
             r'\usepackage{fontspec}',
             r'\setmonofont{DejaVu Sans Mono}[Scale=0.84]',
             r'\newfontfamily\unifont{DejaVu Serif}[Scale=0.9]',
             r'\usepackage{amsmath,amssymb,mathtools}',
             r'\setcounter{MaxMatrixCols}{20}',
             r'\usepackage{longtable,booktabs,array}',
             r'\usepackage{fancyvrb}',
             r'\usepackage[dvipsnames]{xcolor}',
             r'\usepackage[colorlinks=true,linkcolor=NavyBlue,urlcolor=NavyBlue]{hyperref}',
             r'\usepackage{microtype}',
             r'\usepackage{needspace}',
             r'\usepackage{tocloft}\setlength{\cftbeforesecskip}{0.45em}',
             r'\setlength{\parindent}{0pt}\setlength{\parskip}{0.55em}',
             r'\setlength{\emergencystretch}{8em}',
             r'\renewcommand{\arraystretch}{1.15}',
             r'\allowdisplaybreaks',
             r'\begin{document}',
             r'\title{' + tex_inline(self.title) + (r'\\[0.4em]\large ' + tex_inline(self.subtitle) if self.subtitle else '') + '}',
             r'\author{' + r'\\ '.join(guard(tex_inline(m)) for m in self.meta) + '}',
             r'\date{}', r'\maketitle']
        if self.abstract:
            T += [r'\begin{abstract}', tex_inline(self.abstract), r'\end{abstract}']
        T.append(r'\tableofcontents')
        for idx, b in enumerate(self.blocks):
            k = b[0]
            nxt = self.blocks[idx + 1][0] if idx + 1 < len(self.blocks) else None
            if k == 'h1': T.append(r'\Needspace{12\baselineskip}\section*{' + tex_inline(b[1], True) + r'}\addcontentsline{toc}{section}{' + tex_inline(b[1], True) + '}')
            elif k == 'h2': T.append(r'\Needspace{8\baselineskip}\subsection*{' + tex_inline(b[1], True) + r'}\addcontentsline{toc}{subsection}{' + tex_inline(b[1], True) + '}')
            elif k == 'h3': T.append(r'\Needspace{6\baselineskip}\subsubsection*{' + tex_inline(b[1], True) + '}')
            elif k == 'p':
                # a paragraph that introduces a display stays in the same TeX paragraph (no page break before the display)
                T.append(tex_inline(b[1]) + ('' if nxt in ('eq', 'align') else '\n'))
            elif k == 'quote': T += [r'\begin{quote}', tex_inline(b[1]), r'\end{quote}']
            elif k == 'eq':
                env = 'gather*' if top_breaks(b[1]) else 'equation*'
                T += [r'\begin{%s}' % env, b[1].strip(), r'\end{%s}' % env]
            elif k == 'align':
                T += [r'\begin{align*}'] + [guard(ln) + r' \\' if i < len(b[1]) - 1 else guard(ln) for i, ln in enumerate(b[1])] + [r'\end{align*}']
            elif k in ('ul', 'ol'):
                env = 'itemize' if k == 'ul' else 'enumerate'
                T.append(r'\begin{%s}' % env)
                for it in b[1]: T.append(r'\item ' + guard(tex_inline(it)))
                T.append(r'\end{%s}' % env)
            elif k == 'code':
                T += verbatim(b[1])
            elif k == 'rule': T.append(r'\par\noindent\rule{\linewidth}{0.4pt}\par')
            elif k == 'table':
                header, rows, widths, small = b[1], b[2], b[3], b[4]
                n = len(header)
                if not widths: widths = [1.0 / n] * n
                s = float(sum(widths)); widths = [w / s for w in widths]
                spec = ''.join(r'>{\raggedright\arraybackslash}p{\dimexpr %.4f\linewidth-2\tabcolsep\relax}' % w for w in widths)
                T.append(r'{\small' if not small else r'{\footnotesize')
                T.append(r"\begin{longtable}{" + spec + "}")
                head = ' & '.join(guard(r'\textbf{' + tex_inline(str(h)) + '}') for h in header) + r' \\'
                # the bottom rule belongs to the last foot, so it can never start a page of its own
                T += [r'\toprule', head, r'\midrule\endfirsthead', r'\toprule', head, r'\midrule\endhead',
                      r'\endfoot', r'\bottomrule\endlastfoot']
                for j, r in enumerate(rows):
                    # \\* after the last row: no page break between it and the closing rule
                    T.append(' & '.join(guard(tex_inline(str(c))) for c in r) + (r' \\*' if j == len(rows) - 1 else r' \\'))
                T.append(r'\end{longtable}}')
        T.append(r'\end{document}')
        return '\n'.join(T) + '\n'
