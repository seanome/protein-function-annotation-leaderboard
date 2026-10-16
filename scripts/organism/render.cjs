// Run docs/assets/organism.js in Node for every pre-built species, so pre-built pages and pages
// built in the browser come from one piece of code.
//   node scripts/organism/render.cjs <strata.json> <lookups.jsonl> <out.jsonl>
// Each output line: {"answer": ..., "html": ...}.
"use strict";
const fs = require("fs");
const path = require("path");
const O = require(path.resolve(__dirname, "../../docs/assets/organism.js"));
const [strataPath, lookupsPath, outPath] = process.argv.slice(2);
const S = JSON.parse(fs.readFileSync(strataPath, "utf8"));
const out = fs.openSync(outPath, "w");
let n = 0;
for (const line of fs.readFileSync(lookupsPath, "utf8").split("\n")) {
  if (!line.trim()) continue;
  const L = JSON.parse(line);
  const A = O.computeAnswer(L, S);
  fs.writeSync(out, JSON.stringify({ answer: A, html: O.renderAnswer(A, S) }) + "\n");
  n++;
}
fs.closeSync(out);
console.log(`rendered ${n} species`);
