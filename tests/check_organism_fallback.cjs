// No pre-built species has a mean pLDDT in a band with fewer than 200 measured domains, so
// this drives the "fewer than 200" rule with the oyster's lookups and its AlphaFold models
// moved to pLDDT 95-100 (the 90-100 band holds 186 domains).
//   node tests/check_organism_fallback.cjs docs
"use strict";
const fs = require("fs");
const path = require("path");
const site = process.argv[2] || "docs";
const O = require(path.resolve(site, "assets/organism.js"));
const S = JSON.parse(fs.readFileSync(path.join(site, "data/organism/strata.json"), "utf8"));
const page = JSON.parse(fs.readFileSync(path.join(site, "organism/magallana_gigas.html.json"), "utf8"));
const L = JSON.parse(JSON.stringify(page.data.lookups));
for (const k of Object.keys(L.alphafold.bands5)) L.alphafold.bands5[k] = 0;
L.alphafold.bands5.plddt_95_100 = 1000;
const A = O.computeAnswer(L, S);
const html = O.renderAnswer(A, S);
const fails = [];
if (A.stratum.band_asked !== "90-100") fails.push(`band asked ${A.stratum.band_asked}`);
if (!(A.stratum.n_asked < 200)) fails.push(`n asked ${A.stratum.n_asked}`);
if (!A.stratum.fallback || A.stratum.band !== "70-90") fails.push(`fallback ${A.stratum.fallback} band ${A.stratum.band}`);
const notice = html.indexOf('class="notice">Fewer than 200');
if (notice < 0 || notice > html.indexOf('class="verdict"')) fails.push("the fewer-than-200 notice is not shown before the tools");
console.log(fails.length ? "FAIL " + fails.join("; ") : `ok   fallback: asked 90-100 (${A.stratum.n_asked} domains), ranked in ${A.stratum.band}`);
process.exit(fails.length ? 1 : 0);
