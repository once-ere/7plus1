#!/usr/bin/env python3
"""Lenient parser for Wolfram notebook (.nb) files.

Produces a tree of:
  Expr(head, args)   for  head[args...]
  ListE(items)       for  {items...}
  Str(s)             for  "..." (decoded)
  Num(text)          for  numbers
  Sym(name)          for  symbols
  Op(text)           for  operator tokens
  Seq(items)         for  argument slots that contain operators (e.g. a -> b)
and a box-to-text renderer for reading BoxData of Input/Output cells.
"""
import re, sys

NAMED = {
 'Tau':'τ','Gamma':'γ','Eta':'η','Mu':'μ','Nu':'ν','Alpha':'α','Beta':'β','Delta':'δ',
 'Epsilon':'ε','Sigma':'σ','Psi':'ψ','Phi':'φ','Rho':'ρ','Lambda':'λ','Omega':'ω','Kappa':'κ',
 'Theta':'θ','Chi':'χ','Xi':'ξ','Zeta':'ζ','Iota':'ι','Pi':'π','Upsilon':'υ','Omicron':'ο',
 'CapitalGamma':'Γ','CapitalDelta':'Δ','CapitalOmega':'Ω','CapitalSigma':'Σ','CapitalPhi':'Φ',
 'CapitalPsi':'Ψ','CapitalLambda':'Λ','CapitalTheta':'Θ','CapitalPi':'Π','CapitalXi':'Ξ',
 'CapitalUpsilon':'Υ','CapitalChi':'Χ','CapitalEpsilon':'Ε','CapitalEta':'Η','CapitalTau':'Τ',
 'CapitalMu':'Μ','CapitalNu':'Ν','CapitalAlpha':'Α','CapitalBeta':'Β','CapitalKappa':'Κ',
 'CapitalRho':'Ρ','CapitalIota':'Ι','CapitalZeta':'Ζ','CapitalOmicron':'Ο',
 'CurlyPhi':'ϕ','CurlyEpsilon':'ϵ','CurlyTheta':'ϑ','CurlyRho':'ϱ','CurlyPi':'ϖ','CurlyKappa':'ϰ',
 'PartialD':'∂','Element':'∈','Equal':'==','Rule':'->','RuleDelayed':':>','Times':'×',
 'Transpose':'ᵀ','ConjugateTranspose':'ᴴ','Dagger':'†','IndentingNewLine':'\n','NewLine':'\n',
 'InvisibleSpace':'','InvisibleTimes':' ','InvisibleComma':',','InvisibleApplication':'',
 'UnderBrace':'⏟','OverBrace':'⏞','Infinity':'∞','ImaginaryI':'I','ExponentialE':'E',
 'ImaginaryJ':'I','Function':'↦','LeftDoubleBracket':'[[','RightDoubleBracket':']]',
 'CenterDot':'·','CircleTimes':'⊗','TensorProduct':'⊗','TensorWedge':'∧','Wedge':'∧',
 'Sum':'Σ','Product':'∏','Integral':'∫','DifferentialD':'d','Degree':'°','Cross':'×',
 'LongEqual':'==','NotEqual':'!=','LessEqual':'<=','GreaterEqual':'>=','And':'&&','Or':'||',
 'Not':'!','Implies':'⇒','Equivalent':'⇔','RightArrow':'→','LeftArrow':'←','LeftRightArrow':'↔',
 'DoubleRightArrow':'⇒','Rule':'->','Distributed':'~','Vee':'∨','SmallCircle':'∘',
 'Hacek':'ˇ','Bullet':'•','Ellipsis':'…','CenterEllipsis':'⋯','VerticalEllipsis':'⋮',
 'Placeholder':'□','SelectionPlaceholder':'■','EmptySet':'∅','DotlessI':'ı','Micro':'µ',
 'HBar':'ℏ','ScriptCapitalL':'ℒ','ScriptCapitalD':'𝒟','ScriptCapitalE':'ℰ','ScriptL':'ℓ',
 'Del':'∇','Square':'□','PlusMinus':'±','MinusPlus':'∓','Sqrt':'√','Proportional':'∝',
 'Congruent':'≡','Tilde':'~','TildeTilde':'≈','Prime':'′','DoublePrime':'″',
 'LeftAngleBracket':'⟨','RightAngleBracket':'⟩','LeftBracketingBar':'|','RightBracketingBar':'|',
 'LeftDoubleBracketingBar':'‖','RightDoubleBracketingBar':'‖','DirectedEdge':'→',
 'Conjugate':'*','HermitianConjugate':'†','NonBreakingSpace':' ','VeryThinSpace':' ',
 'ThinSpace':' ','MediumSpace':' ','ThickSpace':' ','SpaceIndicator':' ','AlignmentMarker':'',
 'Continuation':'','LineSeparator':'\n','ParagraphSeparator':'\n\n','Null':'',
 'DoubleStruckCapitalR':'ℝ','DoubleStruckCapitalC':'ℂ','DoubleStruckCapitalH':'ℍ',
 'DoubleStruckCapitalO':'𝕆','DoubleStruckCapitalZ':'ℤ','DoubleStruckCapitalN':'ℕ',
 'DoubleStruckCapitalQ':'ℚ','DoubleStruckOne':'𝟙','GothicCapitalC':'ℭ','ScriptCapitalC':'𝒞',
 'Equilibrium':'⇌','UpArrow':'↑','DownArrow':'↓','Diameter':'⌀','Angle':'∠','Star':'⋆',
 'FivePointedStar':'★','Checkmark':'✓','WarningSign':'⚠','LongDash':'—','Dash':'–',
 'OpenCurlyQuote':'‘','CloseCurlyQuote':'’','OpenCurlyDoubleQuote':'“','CloseCurlyDoubleQuote':'”',
 'Section':'§','Paragraph':'¶','Copyright':'©','RegisteredTrademark':'®','Trademark':'™',
 'Eth':'ð','Aring':'å','ADoubleDot':'ä','ODoubleDot':'ö','UDoubleDot':'ü','EAcute':'é',
 'CapitalADoubleDot':'Ä','CapitalODoubleDot':'Ö','CapitalUDoubleDot':'Ü','SZ':'ß',
 'OverBracket':'⎴','UnderBracket':'⎵','OverParenthesis':'⏜','UnderParenthesis':'⏝',
 'HorizontalLine':'─','VerticalLine':'│','DoubleVerticalBar':'∥','VerticalBar':'|',
 'Integers':'ℤ','Reals':'ℝ','Complexes':'ℂ','Rationals':'ℚ','Colon':':','Backslash':'\\',
 'Divide':'÷','Minus':'−','Piecewise':'{','Union':'∪','Intersection':'∩','Subset':'⊂',
 'SubsetEqual':'⊆','Superset':'⊃','ForAll':'∀','Exists':'∃','NotElement':'∉',
 'Therefore':'∴','Because':'∵','Mho':'℧','Angstrom':'Å','GothicCapitalR':'ℜ','GothicCapitalI':'ℑ',
}

class Expr:
    __slots__=('head','args')
    def __init__(s,h,a): s.head=h; s.args=a
    def __repr__(s): return f"Expr({s.head!r},{s.args!r})"
class ListE:
    __slots__=('items',)
    def __init__(s,i): s.items=i
    def __repr__(s): return f"List({s.items!r})"
class Str:
    __slots__=('s',)
    def __init__(s,x): s.s=x
    def __repr__(s): return f"Str({s.s!r})"
class Num:
    __slots__=('t',)
    def __init__(s,x): s.t=x
    def __repr__(s): return f"Num({s.t})"
class Sym:
    __slots__=('n',)
    def __init__(s,x): s.n=x
    def __repr__(s): return f"Sym({s.n})"
class Op:
    __slots__=('t',)
    def __init__(s,x): s.t=x
    def __repr__(s): return f"Op({s.t})"
class Seq:
    __slots__=('items',)
    def __init__(s,i): s.items=i
    def __repr__(s): return f"Seq({s.items!r})"

# every named character that tools/nbwrite.py may write is decoded as well
try:
    from nbwrite import NAMED as _WRITE_NAMED
    for _ch, _nm in _WRITE_NAMED.items(): NAMED.setdefault(_nm, _ch)
except Exception:
    pass

def decode_named(name):
    return NAMED.get(name, '\\['+name+']')

def decode_string(raw, keep_markers=False):
    """raw: content between the quotes, with escapes."""
    out=[]; i=0; n=len(raw)
    while i<n:
        c=raw[i]
        if c=='\\' and i+1<n:
            d=raw[i+1]
            if d=='\n': i+=2; continue
            if d=='\r' and i+2<n and raw[i+2]=='\n': i+=3; continue
            if d=='"': out.append('"'); i+=2; continue
            if d=='\\': out.append('\\'); i+=2; continue
            if d=='n': out.append('\n'); i+=2; continue
            if d=='t': out.append('\t'); i+=2; continue
            if d=='r': i+=2; continue
            if d=='<' or d=='>':
                if keep_markers: out.append('\\'+d)
                i+=2; continue
            if d=='[':
                j=raw.find(']',i+2)
                if j>0:
                    out.append(decode_named(raw[i+2:j])); i=j+1; continue
            if d==':' and i+5<n+1:
                try:
                    out.append(chr(int(raw[i+2:i+6],16))); i+=6; continue
                except Exception: pass
            if d=='|':
                try:
                    out.append(chr(int(raw[i+2:i+8],16))); i+=8; continue
                except Exception: pass
            if d=='.':
                try:
                    out.append(chr(int(raw[i+2:i+4],16))); i+=4; continue
                except Exception: pass
            if d in '!%^&*()_+=-{}|;\'?/,`~ ':
                # box escape sequences like \! \( \) \* \^ \_ etc. keep minimal
                out.append({'!':'','(':'(',')':')','*':'','^':'^','_':'_','%':'%','&':'&','+':'+','/':'/','@':'@','`':'`',' ':' '}.get(d,d)); i+=2; continue
            out.append(d); i+=2; continue
        out.append(c); i+=1
    return ''.join(out)

TOK_RE = re.compile(r'''
 (?P<ws>\s+)
|(?P<comment>\(\*)
|(?P<str>")
|(?P<num>(?:\d+\^\^[0-9a-zA-Z.]+|\d+\.?\d*(?:`[\d.\-+]*`?)?(?:\*\^[-+]?\d+)?|\.\d+(?:`[\d.\-+]*)?(?:\*\^[-+]?\d+)?))
|(?P<sym>(?:[A-Za-z$]|\\\[[A-Za-z]+\])(?:[A-Za-z0-9$`]|\\\[[A-Za-z]+\])*|`[A-Za-z$][A-Za-z0-9$`]*)
|(?P<slot>\#\#?\d*)
|(?P<op>:>|->|\^:=|:=|===|=!=|==|!=|<=|>=|&&|\|\||//\.|/\.|//@|//|/@|@@@|@@|<>|\+\+|--|\+=|-=|\*=|/=|\^=|;;|::|\.\.\.|\.\.|[-+*/^=<>!&@;.?:'|~_%])
|(?P<open>[\[\{\(])
|(?P<close>[\]\}\)])
|(?P<comma>,)
''', re.X)

class Parser:
    def __init__(s, text):
        s.t=text; s.i=0; s.n=len(text)
    def skip_ws(s):
        t=s.t
        while s.i<s.n:
            c=t[s.i]
            if c in ' \t\r\n': s.i+=1; continue
            if c=='\\' and s.i+1<s.n and t[s.i+1] in '\r\n': s.i+=2; continue
            if t.startswith('(*',s.i):
                depth=1; s.i+=2
                while s.i<s.n and depth>0:
                    if t.startswith('(*',s.i): depth+=1; s.i+=2
                    elif t.startswith('*)',s.i): depth-=1; s.i+=2
                    else: s.i+=1
                continue
            break
    def read_string(s):
        assert s.t[s.i]=='"'
        j=s.i+1; t=s.t
        while True:
            k=t.find('"',j)
            # count preceding backslashes
            b=0; m=k-1
            while t[m]=='\\': b+=1; m-=1
            if b%2==0:
                raw=t[s.i+1:k]; s.i=k+1; return raw
            j=k+1
    def next_token(s):
        s.skip_ws()
        if s.i>=s.n: return ('eof',None)
        if s.t[s.i]=='"':
            raw=s.read_string(); return ('str',raw)
        m=TOK_RE.match(s.t,s.i)
        if not m: raise ValueError(f"bad token at {s.i}: {s.t[s.i:s.i+50]!r}")
        s.i=m.end()
        k=m.lastgroup
        return (k,m.group(k))
    def peek(s):
        save=s.i; tok=s.next_token(); s.i=save; return tok
    def parse_seq(s, closers):
        """parse items until a comma or closer at this level; return node"""
        items=[]
        while True:
            save=s.i
            k,v=s.next_token()
            if k=='eof': break
            if k=='comma' or (k=='close' and v in closers):
                s.i=save; break
            if k=='close':
                raise ValueError(f"unexpected closer {v} at {s.i}")
            if k=='str': node=Str(v)
            elif k=='num': node=Num(v)
            elif k=='sym': node=Sym(v)
            elif k=='slot': node=Sym(v)
            elif k=='op': node=Op(v)
            elif k=='open':
                if v=='{':
                    node=ListE(s.parse_args('}'))
                elif v=='(':
                    inner=s.parse_seq(')'); kk,vv=s.next_token(); assert vv==')', (kk,vv,s.i)
                    node=Expr(Sym('Paren'),[inner])
                else:  # '[' -> part or call on previous
                    if s.t.startswith('[',s.i):  # [[ part
                        s.i+=1
                        args=s.parse_args(']'); kk,vv=s.next_token(); assert vv==']'
                        prev=items.pop() if items else Sym('?')
                        node=Expr(Sym('Part'),[prev]+args)
                    else:
                        args=s.parse_args(']')
                        prev=items.pop() if items else Sym('?')
                        node=Expr(prev,args)
            else:
                raise ValueError((k,v))
            # call syntax: head immediately followed by [
            items.append(node)
            # handle chained calls f[a][b]
            while s.i<s.n and s.t[s.i]=='[' and not s.t.startswith('[[',s.i):
                s.i+=1
                args=s.parse_args(']')
                items[-1]=Expr(items[-1],args)
        if len(items)==1: return items[0]
        if not items: return Seq([])
        return Seq(items)
    def parse_args(s, closer):
        args=[]
        # empty?
        s.skip_ws()
        if s.i<s.n and s.t[s.i]==closer:
            s.i+=1; return args
        while True:
            node=s.parse_seq((closer,))
            args.append(node)
            k,v=s.next_token()
            if k=='comma': continue
            if k=='close' and v==closer: return args
            raise ValueError(f"expected , or {closer} got {k} {v} at {s.i}: {s.t[s.i-40:s.i+40]!r}")

def parse_notebook(path):
    text=open(path,encoding='utf-8',errors='replace').read()
    start=text.find('Notebook[')
    p=Parser(text); p.i=start
    node=p.parse_seq(())
    return node

# ------------------------------------------------------------------ box rendering
def hname(e):
    if isinstance(e,Expr) and isinstance(e.head,Sym): return e.head.n
    return None

def box_to_text(b, ctx='in'):
    """Render box expression to approximate InputForm text."""
    if isinstance(b,Str):
        s=decode_string(b.s)
        return s
    if isinstance(b,(Sym,Num)):
        return b.n if isinstance(b,Sym) else b.t
    if isinstance(b,ListE):
        return ''.join(box_to_text(x,ctx) for x in b.items)
    if isinstance(b,Seq):
        return ' '.join(box_to_text(x,ctx) for x in b.items)
    if isinstance(b,Expr):
        h=hname(b); a=b.args
        R=lambda x: box_to_text(x,ctx)
        if h=='RowBox': return R(a[0])
        if h=='BoxData': return R(a[0]) if a else ''
        if h=='SubscriptBox': return f"Subscript[{R(a[0])}, {R(a[1])}]"
        if h=='SuperscriptBox': return f"({R(a[0])})^({R(a[1])})"
        if h=='SubsuperscriptBox': return f"Subsuperscript[{R(a[0])}, {R(a[1])}, {R(a[2])}]"
        if h=='OverscriptBox':
            if R(a[1]) in ('\\[RightVector]', '\u21c0', '\u2192'): return R(a[0]) + '\u20d7'   # OverVector: x with an arrow
            return f"Overscript[{R(a[0])}, {R(a[1])}]"
        if h=='UnderscriptBox': return f"Underscript[{R(a[0])}, {R(a[1])}]"
        if h=='UnderoverscriptBox': return f"Underoverscript[{R(a[0])}, {R(a[1])}, {R(a[2])}]"
        if h=='FractionBox': return f"(({R(a[0])})/({R(a[1])}))"
        if h=='SqrtBox': return f"Sqrt[{R(a[0])}]"
        if h=='RadicalBox': return f"Surd[{R(a[0])}, {R(a[1])}]"
        if h in ('StyleBox','TagBox','FormBox','AdjustmentBox','FrameBox','ItemBox','PaneBox','PanelBox','TooltipBox','ButtonBox','ActionMenuBox','PaneSelectorBox','Cell','TextData','HighlightedBox','DynamicWrapperBox'):
            return R(a[0]) if a else ''
        if h=='InterpretationBox': return R(a[0])
        if h=='TemplateBox':
            inner=a[0]; tag=R(a[1]) if len(a)>1 else ''
            if isinstance(inner,ListE):
                parts=[R(x) for x in inner.items]
                if tag in ('NotationTemplateTag','NotationPatternTag','NotationMadeBoxesTag'): return ''.join(parts)
                return f"{tag}<" + ', '.join(parts) + ">"
            return R(inner)
        if h=='GridBox':
            rows=a[0]
            if isinstance(rows,ListE):
                rr=[]
                for r in rows.items:
                    if isinstance(r,ListE): rr.append('{'+', '.join(R(c) for c in r.items)+'}')
                    else: rr.append(R(r))
                return '{'+', '.join(rr)+'}'
            return R(rows)
        if h=='GraphicsBox' or h=='Graphics3DBox': return '<<Graphics>>'
        if h=='DynamicBox' or h=='DynamicModuleBox': return '<<Dynamic>>'
        if h=='CheckboxBox' or h=='SliderBox' or h=='InputFieldBox': return '<<control>>'
        # generic
        return f"{h}[" + ', '.join(R(x) for x in a) + "]"
    return str(b)

def cell_info(c):
    """Return (style, content_node, options) for Cell expression."""
    if hname(c)!='Cell': return None
    a=c.args
    content=a[0] if a else None
    style=None; opts=[]
    for x in a[1:]:
        if isinstance(x,Str) and style is None: style=decode_string(x.s)
        else: opts.append(x)
    return style, content, opts

def iter_cells(node, depth=0, out=None):
    if out is None: out=[]
    if hname(node)=='Notebook':
        for c in node.args[0].items: iter_cells(c,depth,out)
        return out
    if hname(node)=='Cell':
        a=node.args
        if a and hname(a[0])=='CellGroupData':
            cg=a[0]
            for c in cg.args[0].items: iter_cells(c,depth+1,out)
            return out
        out.append((depth,node))
    return out

if __name__=='__main__':
    nb=parse_notebook(sys.argv[1])
    cells=iter_cells(nb)
    print(len(cells),'cells')
