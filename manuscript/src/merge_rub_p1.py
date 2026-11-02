# -*- coding: utf-8 -*-
"""Selective merge of the Rubriq (AJE) language edit of the Paper 2 manuscript: keep the edit where it changed only spelling, hyphenation or punctuation (American style);
restore the original paragraph (with American spelling applied) wherever words were added, removed or rewritten, or where symbols, tabs or equation numbering were damaged.
    python merge_rub_p2.py <mine.docx> <rub.docx> <out.docx>"""
import copy
import re
import sys

import docx
from docx.oxml.ns import qn

MINE, RUB, OUT = sys.argv[1:4]
AM = [("favourable", "favorable"), ("Favourable", "Favorable"), ("unfavourable", "unfavorable"), ("Unfavourable", "Unfavorable"), ("annualised", "annualized"), ("Annualised", "Annualized"), ("standardised", "standardized"),
      ("Standardised", "Standardized"), ("standardisation", "standardization"), ("centred", "centered"), ("judgement", "judgment"), ("maximised", "maximized"), ("maximises", "maximizes"), ("minimised", "minimized"),
      ("modelled", "modeled"), ("programmes", "programs"), ("Programmes", "Programs"), ("programme", "program"), ("emphasised", "emphasized"), ("non-dominated", "nondominated"), ("co-benefits", "cobenefits"),
      ("multi-objective", "multiobjective"), ("Multi-objective", "Multiobjective"), ("towards", "toward"), ("favour", "favor"), ("landfill-avoidance", "landfill avoidance"), ("water-quality", "water quality"),
      ("reef-survival", "reef survival"), ("oyster-survival", "oyster survival"), ("decision-maker", "decision maker"), ("analyses an", "analyzes an"), ("licence", "license"), ("behaviour", "behavior")]
BLACK = ["Leave one", "sh, gross", "upper mode", "are drawn", "is understood", "validate the", "a total of"]


def am(t):
    for a, b in AM:
        t = t.replace(a, b)
    return t


def norm(t):
    return re.sub(r"[^a-z0-9]", "", am(t).lower())


def numtok(t):
    t = t.replace("−", "-").replace("–", "-")
    return re.findall(r"[-+]?\d[\d.,]*", t)


def plist(d):
    return [p for p in d.element.body.iter(qn("w:p")) if "".join(x.text or "" for x in p.iter(qn("w:t"))).strip()]


def text(p):
    return "".join(x.text or "" for x in p.iter(qn("w:t")))


mine = docx.Document(MINE); rub = docx.Document(RUB)
pm, pr = plist(mine), plist(rub)
assert len(pm) == len(pr), (len(pm), len(pr))
acc = res = 0; log = []
for i, (a, b) in enumerate(zip(pm, pr)):
    ta, tb = text(a), text(b)
    if ta == tb:
        continue
    force = ta.startswith("A Coupled Chemo-Mechanical")
    ok = (not force) and norm(ta) == norm(tb) and numtok(ta) == numtok(tb) and not any(k in tb for k in BLACK) and (("\t" in ta) == ("\t" in tb)) and ("," not in ta.replace(",", "", 0) or True)
    if ok:
        acc += 1; continue
    new = copy.deepcopy(a)
    for t in new.iter(qn("w:t")):
        if t.text:
            t.text = am(t.text)
    b.addnext(new); b.getparent().remove(b); res += 1; log.append(i)
rub.save(OUT)
print("accepted (spelling/punctuation only):", acc, "restored:", res)
