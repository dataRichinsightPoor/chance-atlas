# Builds web/methods.html from web/MATH.md so the two never drift. Headings may carry explicit ids: "## Title {#id}".
import re, html
src=open('web/MATH.md').read().split('\n')
out=[]; i=0; para=[]
def inline(t):
    t=html.escape(t,quote=False)
    return re.sub(r'(https?://[^\s<]+)',r'<a href="\1">\1</a>',t)
def flush():
    global para
    if para: out.append('<p>'+inline(' '.join(para))+'</p>'); para=[]
while i<len(src):
    l=src[i]
    if l.startswith('\\['):
        flush(); blk=[]
        while not src[i].startswith('\\]'): blk.append(src[i]); i+=1
        blk.append(src[i]); out.append('<div class="math">'+html.escape('\n'.join(blk),quote=False)+'</div>')
    elif l.startswith('# '): flush(); out.append('<h1>'+inline(l[2:])+'</h1>')
    elif l.startswith('## '):
        flush(); t=l[3:]; m=re.search(r'\s*\{#([a-z0-9-]+)\}\s*$',t)
        slug = m.group(1) if m else re.sub(r'[^a-z0-9]+','-',t.lower()).strip('-')
        if m: t=t[:m.start()]
        out.append(f'<h2 id="{slug}">'+inline(t)+'</h2>')
    elif l.strip()=='': flush()
    else: para.append(l)
    i+=1
flush()
body='\n'.join(out).replace('<h1>Chance Atlas: mathematical contract</h1>','<h1>Where the randomness enters,<br>with the assumptions showing.</h1>',1)
page=f'''<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chance Atlas · Math, assumptions &amp; tests · v0.1.0-alpha</title><link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Instrument+Serif&display=swap" rel="stylesheet"><link rel="stylesheet" href="style.css">
<script>window.MathJax={{tex:{{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}},svg:{{fontCache:'global'}}}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<style>.methods h2{{scroll-margin-top:20px}}</style></head><body>
<header><a href="index.html">← Chance Atlas</a><nav aria-label="Resources"><a href="models.js" download>Download model code</a><a href="MATH.md">Math as Markdown</a></nav></header>
<main class="methods"><span class="badge">DATA-RICH, INSIGHT-POOR · 99 SMALL PROBLEMS · ATLAS · v0.1.0-alpha</span>
{body}
</main></body></html>'''
open('web/methods.html','w').write(page); print('ok', body.count('<h2'))
