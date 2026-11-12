"""Ask a model for gene cards, a batch of proteins per call. Each call gets prompt.md, the
evidence sheets (never a sequence) and a JSON schema for the answer. Answers are appended to
out/cards.<model>.jsonl; a protein already there is not asked again, so a stopped run resumes.

Two ways to reach a model:

* Claude, through the `claude` command line (Olga's Claude subscription; no API key on this
  machine). --model takes the command line's names (haiku, sonnet, opus); the model that
  answered, as the command line reports it, is stored with each card.

    python gene_cards/run_cards.py gene_cards/out --model haiku --batch 10
    python gene_cards/run_cards.py gene_cards/out --model opus --effort high --only-kind conflict_different_gene --batch 5

* An open-weights model served by vLLM on the same Sherlock GPU node (--server), so the
  evidence never leaves Sherlock; see deepseek_cards.sbatch. --model is the name vLLM serves.
  DeepSeek's R1 models are asked without a system prompt, as DeepSeek recommends, so prompt.md
  goes at the top of the user message; vLLM's deepseek_r1 reasoning parser keeps the model's
  thinking out of the answer, and the answer is held to the JSON schema.

    python gene_cards/run_cards.py gene_cards/out --server http://127.0.0.1:8000 \
        --model DeepSeek-R1-Distill-Qwen-32B --batch 1 --parallel 16
    python gene_cards/run_cards.py gene_cards/out --server http://127.0.0.1:8000 \
        --model DeepSeek-V4-Flash --temperature 1 --top-p 1 --reasoning-effort high --batch 1 --parallel 16
"""
import argparse
import json
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERDICTS = ["existing", "proposed", "same_gene", "family_only", "neither", "not_applicable"]
CARD = {
    "type": "object", "additionalProperties": False,
    "required": ["protein_id", "card_name", "one_line", "summary", "confidence", "confidence_reason", "name_verdict"],
    "properties": {
        "protein_id": {"type": "string"},
        "card_name": {"type": "string"},
        "one_line": {"type": "string"},
        "summary": {"type": "string"},
        "confidence": {"type": "string", "enum": ["high", "medium", "low", "none"]},
        "confidence_reason": {"type": "string"},
        "name_verdict": {"type": "object", "additionalProperties": False,
                         "required": ["choice", "recommended_name", "reason"],
                         "properties": {"choice": {"type": "string", "enum": VERDICTS},
                                        "recommended_name": {"type": "string"},
                                        "reason": {"type": "string"}}},
    },
}
SCHEMA = {"type": "object", "additionalProperties": False, "required": ["cards"],
          "properties": {"cards": {"type": "array", "items": CARD}}}


def ask(batch, model, effort):
    user = "Write a gene card for each of these %d proteins.\n\n" % len(batch) + "\n\n".join(
        "=== %s ===\n%s" % (b["protein_id"], b["evidence"]) for b in batch)
    cmd = ["claude", "-p", "--model", model, "--effort", effort, "--output-format", "json",
           "--system-prompt", (HERE / "prompt.md").read_text(), "--json-schema", json.dumps(SCHEMA),
           "--tools", "", "--strict-mcp-config", "--setting-sources", "", "--no-session-persistence"]
    for attempt in range(3):
        r = subprocess.run(cmd, input=user, capture_output=True, text=True, timeout=1800)
        try:
            res = json.loads(r.stdout)
            out = res.get("structured_output") or json.loads(res["result"])
            cards = {c["protein_id"]: c for c in out["cards"]}
            missing = [b["protein_id"] for b in batch if b["protein_id"] not in cards]
            if missing:
                raise ValueError("no card for %s" % missing)
            used = ",".join(res.get("modelUsage", {})) or model
            return [dict(cards[b["protein_id"]], model=used, effort=effort,
                         cost_usd=res.get("total_cost_usd", 0) / len(batch)) for b in batch]
        except Exception as e:
            print("batch %s attempt %d failed: %s %s" % (batch[0]["protein_id"], attempt + 1, e, (r.stderr or r.stdout)[:300]), file=sys.stderr)
            time.sleep(10)
    return []


def ask_server(batch, server, model, temperature=0.6, top_p=0.95, effort=None):
    """One call to an OpenAI-compatible server (vLLM) on this node."""
    user = (HERE / "prompt.md").read_text() + "\n\nWrite a gene card for each of these %d proteins.\n\n" % len(batch) + \
        "\n\n".join("=== %s ===\n%s" % (b["protein_id"], b["evidence"]) for b in batch)
    body = dict(model=model, messages=[{"role": "user", "content": user}], temperature=temperature, top_p=top_p,
                max_tokens=16384, response_format={"type": "json_schema",
                                                   "json_schema": {"name": "cards", "schema": SCHEMA, "strict": True}})
    if effort:
        body["reasoning_effort"] = effort
    req = urllib.request.Request(server.rstrip("/") + "/v1/chat/completions", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=3600) as r:
                res = json.load(r)
            msg = res["choices"][0]["message"]
            cards = {c["protein_id"]: c for c in json.loads(msg["content"])["cards"]}
            missing = [b["protein_id"] for b in batch if b["protein_id"] not in cards]
            if missing:
                raise ValueError("no card for %s" % missing)
            return [dict(cards[b["protein_id"]], model=res.get("model", model), effort=effort or "",
                         reasoning_tokens=len((msg.get("reasoning_content") or msg.get("reasoning") or "").split()),
                         cost_usd=0) for b in batch]
        except Exception as e:
            print("batch %s attempt %d failed: %s" % (batch[0]["protein_id"], attempt + 1, str(e)[:300]), file=sys.stderr)
            time.sleep(5)
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--model", default="haiku")
    ap.add_argument("--effort", default="medium")
    ap.add_argument("--batch", type=int, default=10)
    ap.add_argument("--parallel", type=int, default=3)
    ap.add_argument("--only-kind", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--server", default=None, help="an OpenAI-compatible server on this node, such as vLLM")
    ap.add_argument("--temperature", type=float, default=0.6, help="--server only; DeepSeek's advice: R1 0.6, V4 1.0")
    ap.add_argument("--top-p", type=float, default=0.95, help="--server only; DeepSeek's advice: R1 0.95, V4 1.0")
    ap.add_argument("--reasoning-effort", default=None, help="--server only; V4: high (Think High) or max")
    a = ap.parse_args()
    ev = [json.loads(l) for l in open(a.out / "evidence.jsonl")]
    if a.only_kind:
        ev = [e for e in ev if e["kind"] in a.only_kind]
    path = a.out / ("cards.%s.jsonl" % a.model)
    done = {json.loads(l)["protein_id"] for l in open(path)} if path.exists() else set()
    todo = [e for e in ev if e["protein_id"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    batches = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
    print(len(todo), "proteins in", len(batches), "calls to", a.model, file=sys.stderr)
    with ThreadPoolExecutor(a.parallel) as pool, open(path, "a") as fh:
        call = (lambda b: ask_server(b, a.server, a.model, a.temperature, a.top_p, a.reasoning_effort)) if a.server else (lambda b: ask(b, a.model, a.effort))
        for cards in pool.map(call, batches):
            for c in cards:
                fh.write(json.dumps(c) + "\n")
            fh.flush()
            print("written", len(cards), file=sys.stderr)


if __name__ == "__main__":
    main()
