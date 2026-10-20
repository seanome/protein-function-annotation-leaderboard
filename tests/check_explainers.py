#!/usr/bin/env python3
"""Drive the widgets on explainers.html and the "?" definitions on index.html.

The made-up example is read from the page itself (#toy-data). Every expected number is
computed here in Python, separately from docs/assets/explainers.js:

1. Figure 1: at several cutoffs, the kept and correct counts, precision, recall and F; the
   reported Fmax is the maximum F over cutoffs 0..100, and "Go to the best cutoff" lands on it.
2. Figure 2: AUPRC = sum_k (R_k - R_(k-1)) P_k over the pooled list, for protein A alone and
   with protein B; Fmax of both; the area of the calls taken at a cutoff.
3. Figure 3: IoU for several start and end settings, and the median over correct calls.
4. Figure 4: raw and dedup counts, precision with every call kept, Fmax and AUPRC.
5. Figure 5: cost_usd_est for the default inputs and one change.
6. Keyboard: Tab onto a "?" in Table 1 and in the Figure-free key; its definition shows, ends
   with a link to the explainer, and Escape hides it. A tap elsewhere hides it.
7. No console errors; no sideways page scroll at 400 px.
"""

import argparse
import functools
import http.server
import json
import socketserver
import statistics
import sys
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
TOL = 5e-4  # numbers on the page are shown to 3 decimals


def iou(a, b):
    inter = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    return inter, (a["end"] - a["start"] + 1) + (b["end"] - b["start"] + 1) - inter


def judge(protein, calls):
    taken, out = set(), []
    for c in sorted(calls, key=lambda c: -c["score"]):
        hit = next((i for i, t in enumerate(protein["true"]) if i not in taken and t["family"] == c["family"]
                    and iou(c, t)[0] / iou(c, t)[1] >= 0.5), None)
        if hit is not None:
            taken.add(hit)
        out.append(dict(c, correct=hit is not None, domain=hit))
    return out


def dedup(calls):
    keep = []
    for c in sorted(calls, key=lambda c: -c["score"]):
        if not any(k["family"] == c["family"] and iou(k, c)[0] / iou(k, c)[1] >= 0.5 for k in keep):
            keep.append(c)
    return keep


def at_cutoff(judged, proteins, t):
    precs, recs, kept, correct = [], [], 0, 0
    for calls, p in zip(judged, proteins):
        k = [c for c in calls if c["score"] >= t]
        ok = sum(c["correct"] for c in k)
        kept += len(k); correct += ok
        if k:
            precs.append(ok / len(k))
        recs.append(ok / len(p["true"]))
    P = statistics.fmean(precs) if precs else None
    R = statistics.fmean(recs)
    F = 0.0 if P is None or P + R == 0 else 2 * P * R / (P + R)
    return {"P": P, "R": R, "F": F, "kept": kept, "correct": correct}


def fmax(judged, proteins):
    # Fmax over the call scores themselves (where the kept set changes), not over 0..100,
    # so the page's sweep over whole numbers is checked against a different route.
    scores = sorted({c["score"] for calls in judged for c in calls})
    Fs = {t: at_cutoff(judged, proteins, t)["F"] for t in scores}
    best = max(Fs.values())
    return best, max(t for t, f in Fs.items() if f == best)


def pooled(judged, proteins):
    allc = sorted((c for calls in judged for c in calls), key=lambda c: -c["score"])
    n_true = sum(len(p["true"]) for p in proteins)
    tp, prev, pieces = 0, 0.0, []
    for k, c in enumerate(allc, 1):
        tp += c["correct"]
        R, P = tp / n_true, tp / k
        pieces.append((c["score"], (R - prev) * P))
        prev = R
    return sum(p for _, p in pieces), pieces


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=ROOT / "docs")
    a = ap.parse_args()
    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=str(a.site)))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}/"
    fails, n_checked, errors = [], 0, []

    def close(what, got, want, tol=TOL):
        nonlocal n_checked
        n_checked += 1
        if want is None or got is None:
            if got != want:
                fails.append(f"{what}: page {got}, expected {want}")
        elif abs(float(got) - float(want)) > tol:
            fails.append(f"{what}: page {got}, expected {want:.4f}")

    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        page = browser.new_page(viewport={"width": 1280, "height": 900})
        page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(base + "explainers.html")
        page.wait_for_function("window.__explainers && window.__explainers.ready")
        toy = json.loads(page.eval_on_selector("#toy-data", "e => e.textContent"))
        A, B = toy["proteins"]["A"], toy["proteins"]["B"]
        jA, jB = judge(A, A["calls"]), judge(B, B["calls"])

        def slide(sel, v):
            page.eval_on_selector(sel, "(e, v) => { e.value = v; e.dispatchEvent(new Event('input', {bubbles: true})); }", v)

        def val(sel):
            return page.eval_on_selector(sel, "e => e.dataset.value === undefined ? e.textContent.trim() : e.dataset.value")

        def shown(sel):
            return page.eval_on_selector(sel, "e => e.textContent.trim()")

        # ---- Figure 1 ----
        want_fmax, want_best = fmax([jA], [A])
        for t in (0, 10, 15, 16, 30, 40, 41, 55, 70, 79, 80, 81, 95, 96, 100):
            slide("#f1-cut", t)
            r = at_cutoff([jA], [A], t)
            close(f"fig1 t={t} kept", val("#f1-kept"), r["kept"])
            close(f"fig1 t={t} correct", val("#f1-correct"), r["correct"])
            pv = val("#f1-p")
            close(f"fig1 t={t} precision", float(pv) if pv else None, r["P"])
            close(f"fig1 t={t} recall", val("#f1-r"), r["R"])
            close(f"fig1 t={t} F", val("#f1-f"), r["F"])
            close(f"fig1 t={t} F shown", shown("#f1-f"), r["F"])
            kept_marks = page.eval_on_selector_all("#f1-track [data-state='kept']", "es => es.length")
            close(f"fig1 t={t} kept bars drawn", kept_marks, r["kept"])
        close("fig1 Fmax", val("#f1-fmax"), want_fmax)
        sweep = max(at_cutoff([jA], [A], t)["F"] for t in range(101))
        close("fig1 Fmax is the max over cutoffs 0..100", val("#f1-fmax"), sweep, 1e-12)
        page.click("#f1-go")
        close("fig1 best cutoff button", page.eval_on_selector("#f1-cut", "e => +e.value"), want_best)
        close("fig1 F at best cutoff", val("#f1-f"), want_fmax)

        # ---- Figure 2 ----
        aA, piecesA = pooled([jA], [A])
        aAB, piecesAB = pooled([jA, jB], [A, B])
        fA, fAB = fmax([jA], [A])[0], fmax([jA, jB], [A, B])[0]
        close("fig2 AUPRC A", val("#f2-a-only-auprc"), aA)
        close("fig2 AUPRC A+B", val("#f2-ab-auprc"), aAB)
        close("fig2 Fmax A", val("#f2-a-only-fmax"), fA)
        close("fig2 Fmax A+B", val("#f2-ab-fmax"), fAB)
        if not (fAB > fA and aAB < aA):
            fails.append(f"fig2: the text says adding B raises Fmax and lowers AUPRC; got Fmax {fA:.3f}->{fAB:.3f}, AUPRC {aA:.3f}->{aAB:.3f}")
        for withB, pieces, auc in ((False, piecesA, aA), (True, piecesAB, aAB)):
            page.click(f"#f2-which button[data-b='{int(withB)}']")
            close(f"fig2 B={withB} AUPRC shown", val("#f2-auprc"), auc)
            close(f"fig2 B={withB} rows", page.eval_on_selector_all("#f2-list tbody tr", "es => es.length"), len(pieces))
            for t in (0, 50, 60, 81, 100):
                slide("#f2-cut", t)
                close(f"fig2 B={withB} t={t} area so far", val("#f2-area"), sum(x for s, x in pieces if s >= t))
        page.click("#f2-which button[data-b='0']")

        # ---- Figure 3 ----
        T = toy["iou"]["true"]
        correct_A = [iou(c, A["true"][c["domain"]]) for c in jA if c["correct"]]
        correct_A = [i / u for i, u in correct_A]
        for s, e in ((30, 129), (10, 109), (1, 200), (60, 70), (110, 150), (5, 120), (50, 109), (109, 109)):
            slide("#f3-end", e)
            slide("#f3-start", s)
            slide("#f3-end", e)
            inter, union = iou({"start": s, "end": e}, T)
            close(f"fig3 {s}-{e} shared", val("#f3-inter"), inter)
            close(f"fig3 {s}-{e} either", val("#f3-union"), union)
            close(f"fig3 {s}-{e} IoU", val("#f3-iou"), inter / union)
            vals = correct_A + ([inter / union] if inter / union >= 0.5 else [])
            close(f"fig3 {s}-{e} median", val("#f3-median"), statistics.median(vals))
            close(f"fig3 {s}-{e} n", val("#f3-n"), len(vals))
        # The Methods page's example: 80 / 120 = 0.67.
        slide("#f3-start", 30); slide("#f3-end", 129)
        if shown("#f3-iou") != "0.667":
            fails.append(f"fig3: Methods example shows {shown('#f3-iou')}, not 0.667")

        # ---- Figure 4 ----
        jD = judge(A, dedup(A["calls"]))
        for k, j in (("raw", jA), ("dedup", jD)):
            page.click(f"#f4-which button[data-mode='{k}']")
            close(f"fig4 {k} calls", val(f"#f4-{k}-n"), len(j))
            close(f"fig4 {k} precision all kept", val(f"#f4-{k}-p0"), at_cutoff([j], [A], 0)["P"])
            close(f"fig4 {k} Fmax", val(f"#f4-{k}-fmax"), fmax([j], [A])[0])
            close(f"fig4 {k} AUPRC", val(f"#f4-{k}-auprc"), pooled([j], [A])[0])
            drawn = page.eval_on_selector_all("#f4-track [data-state='kept']", "es => es.length")
            close(f"fig4 {k} bars that count", drawn, len(j))

        # ---- Figure 5 ----
        nums = [float(page.eval_on_selector(f"#{i}", "e => e.value")) for i in ("f5-cpu", "f5-gpu", "f5-rcpu", "f5-rgpu")]
        close("fig5 default cost", val("#f5-out"), nums[0] * nums[2] + nums[1] * nums[3], 1e-12)
        # Drag the core-hours number 40 px right: 10 steps of 0.005 up, and the cost follows.
        sp = page.locator("#fig5 .scrub[data-for='f5-cpu']")
        sp.scroll_into_view_if_needed()
        bb = sp.bounding_box()
        page.mouse.move(bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2)
        page.mouse.down()
        page.mouse.move(bb["x"] + bb["width"] / 2 + 40, bb["y"] + bb["height"] / 2, steps=8)
        page.mouse.up()
        cpu2 = float(page.eval_on_selector("#f5-cpu", "e => e.value"))
        close("fig5 drag core-hours +40 px", cpu2, nums[0] + 10 * 0.005, 1e-9)
        close("fig5 cost after drag", val("#f5-out"), cpu2 * nums[2] + nums[1] * nums[3], 1e-12)
        if sp.inner_text().strip() != str(cpu2):
            fails.append(f"fig5: dragged number shows {sp.inner_text()!r}, box holds {cpu2}")
        # Keyboard: focus the GPU-hours number and press the right arrow twice.
        page.focus("#fig5 .scrub[data-for='f5-gpu']")
        page.keyboard.press("ArrowRight"); page.keyboard.press("ArrowRight")
        close("fig5 arrow keys on GPU-hours", float(page.eval_on_selector("#f5-gpu", "e => e.value")), 0.02, 1e-9)
        page.click("#fig5 details.ex-exact summary")
        page.fill("#f5-cpu", str(nums[0]))
        page.fill("#f5-gpu", "2.5")
        page.fill("#f5-rgpu", "1.006")
        close("fig5 with GPU", val("#f5-out"), nums[0] * nums[2] + 2.5 * 1.006, 1e-12)

        # ---- Figure 6 drag: move the whole call, then one end ----
        slide("#f3-start", 30); slide("#f3-end", 129)
        bar = page.locator("#f3-track [data-drag='move']")
        bar.scroll_into_view_if_needed()
        bb = bar.bounding_box()
        px_per_res = bb["width"] / (129 - 30)  # the bar spans residues 30 to 129
        page.mouse.move(bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2)
        page.mouse.down()
        page.mouse.move(bb["x"] + bb["width"] / 2 + 20 * px_per_res, bb["y"] + bb["height"] / 2, steps=6)
        page.mouse.up()
        a = int(page.eval_on_selector("#f3-start", "e => e.value")); b = int(page.eval_on_selector("#f3-end", "e => e.value"))
        if b - a != 99 or abs(a - 50) > 1:
            fails.append(f"fig6 drag call +20 residues: now {a}-{b}, expected about 50-149")
        inter = max(0, min(b, 109) - max(a, 10) + 1); union = max(b, 109) - min(a, 10) + 1
        close("fig6 IoU after drag", val("#f3-iou"), inter / union)
        h = page.locator("#f3-track [data-drag='end']").bounding_box()
        page.mouse.move(h["x"] + h["width"] / 2, h["y"] + h["height"] / 2)
        page.mouse.down()
        page.mouse.move(h["x"] + h["width"] / 2 - 40 * px_per_res, h["y"] + h["height"] / 2, steps=6)
        page.mouse.up()
        a2 = int(page.eval_on_selector("#f3-start", "e => e.value")); b2 = int(page.eval_on_selector("#f3-end", "e => e.value"))
        if a2 != a or abs(b2 - (b - 40)) > 1:
            fails.append(f"fig6 drag end handle -40: now {a2}-{b2}, expected {a}-{b - 40}")

        # ---- phone width ----
        page.set_viewport_size({"width": 400, "height": 900})
        page.reload()
        page.wait_for_function("window.__explainers && window.__explainers.ready")
        over = page.evaluate("document.documentElement.scrollWidth - window.innerWidth")
        if over > 0:
            fails.append(f"explainers.html at 400 px scrolls sideways by {over} px")
        wide = page.evaluate("[...document.querySelectorAll('.ex-box svg')].filter(s => s.getBoundingClientRect().right > window.innerWidth).length")
        if wide:
            fails.append(f"explainers.html at 400 px: {wide} figure(s) wider than the screen")

        # ---- 6. keyboard and tap on the definitions ----
        page.set_viewport_size({"width": 1280, "height": 900})
        page.goto(base + "index.html")
        page.wait_for_load_state("networkidle")
        heads = page.eval_on_selector_all("#lb thead button.term", "es => es.map(e => e.dataset.term)")
        for t in ("fmax", "auprc", "iou", "n", "cpu", "ram", "cost"):
            if t not in heads:
                fails.append(f"index.html: no definition button for {t} in Table 1")
        page.focus("#lb thead th:nth-child(1) button")  # the Tool column's sort button
        for _ in range(4):  # Tab onto the Fmax "?" (sort Fmax, then its "?")
            page.keyboard.press("Tab")
            if page.evaluate("document.activeElement.dataset.term") == "fmax":
                break
        if page.evaluate("document.activeElement.dataset.term") != "fmax":
            fails.append("keyboard: Tab did not reach the Fmax definition button")
        tip = page.evaluate("""() => { const b = document.activeElement, t = document.getElementById(b.getAttribute('aria-controls'));
            const r = t.getBoundingClientRect();
            return {hidden: t.hidden, w: r.width, h: r.height, desc: b.getAttribute('aria-describedby') === t.id,
                    role: t.getAttribute('role'), link: t.querySelector('a') && t.querySelector('a').getAttribute('href'),
                    last: t.lastElementChild && t.lastElementChild.textContent}; }""")
        if tip["hidden"] or tip["w"] < 50 or tip["h"] < 20:
            fails.append(f"keyboard: Fmax definition did not show on focus: {tip}")
        if not tip["desc"] or tip["role"] != "tooltip" or tip["link"] != "explainers.html#fmax" or tip["last"] != "how it works →":
            fails.append(f"keyboard: Fmax definition is not linked or does not end with its link: {tip}")
        page.keyboard.press("Tab")  # onto the link inside the definition
        if page.evaluate("document.activeElement.getAttribute('href')") != "explainers.html#fmax" or \
                page.evaluate("document.getElementById('tip-fmax-th').hidden"):
            fails.append("keyboard: Tab from the ? did not reach the definition's link with it open")
        page.keyboard.press("Escape")
        if not page.evaluate("document.getElementById('tip-fmax-th').hidden"):
            fails.append("keyboard: Escape did not hide the definition")
        # Tap (click) shows; a click elsewhere hides.
        page.click("button.term[aria-controls='tip-auprc-th']")
        if page.evaluate("document.getElementById('tip-auprc-th').hidden"):
            fails.append("tap: AUPRC definition did not show")
        page.mouse.click(5, 880)
        if not page.evaluate("document.getElementById('tip-auprc-th').hidden"):
            fails.append("tap: AUPRC definition did not hide on a tap elsewhere")
        # Clicking the "?" must not sort the column.
        before = page.evaluate("window.__lb.state().sort")
        page.click("button.term[aria-controls='tip-cost-th']")
        if page.evaluate("window.__lb.state().sort") != before:
            fails.append("clicking the ? of the cost column sorted the table")
        # Every definition link on every page with definitions points at an explainer anchor.
        anchors = None
        for pg in ("index.html", "cost.html"):
            page.goto(base + pg)
            page.wait_for_load_state("networkidle")
            links = page.eval_on_selector_all(".term-tip a", "es => es.map(e => e.getAttribute('href'))")
            if not links:
                fails.append(f"{pg}: no definitions found")
            if anchors is None:
                page.goto(base + "explainers.html")
                anchors = set(page.eval_on_selector_all("[id]", "es => es.map(e => e.id)"))
                page.goto(base + pg)
            for h in links:
                if h.split("#")[1] not in anchors:
                    fails.append(f"{pg}: definition link {h} has no anchor on explainers.html")
        browser.close()
    httpd.shutdown()

    fails += [f"console: {e}" for e in errors]
    print(f"explainers: {n_checked} numbers checked")
    if fails:
        print(f"{len(fails)} FAILURES")
        for f in fails[:60]:
            print("  " + f)
        return 1
    print("explainer checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
