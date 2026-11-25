/*
 * CHLI core: the same pipeline as the Python package (src/chli/pipeline.py), in plain
 * JavaScript so that it can run in a browser without sending data anywhere.
 *
 * Deterministic outputs (scores, alpha, item-total correlations, PCA weights and
 * variances, cut-points, labels, correlations) match the Python package to ~1e-12; this
 * is checked by tests/test_js_parity.py. Two outputs differ by design: McDonald's omega
 * is computed only in Python, and the bootstrap CI for alpha uses a different random
 * stream, so it can differ from Python in the third decimal.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.CHLI = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const isNum = (v) => typeof v === "number" && !Number.isNaN(v);
  const TIE_DECIMALS = 10;

  // ---------------------------------------------------------------- basic statistics
  function mean(a) { let s = 0, n = 0; for (const v of a) if (isNum(v)) { s += v; n++; } return n ? s / n : NaN; }
  function sd(a) {
    const m = mean(a); let s = 0, n = 0;
    for (const v of a) if (isNum(v)) { s += (v - m) * (v - m); n++; }
    return n > 1 ? Math.sqrt(s / (n - 1)) : NaN;
  }
  function variance(a) { const s = sd(a); return s * s; }
  function median(a) {
    const v = a.filter(isNum).sort((x, y) => x - y); const n = v.length;
    if (!n) return NaN; return n % 2 ? v[(n - 1) / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
  }
  function zscore(a) { const m = mean(a), s = sd(a); return a.map((v) => (isNum(v) ? (v - m) / s : NaN)); }
  function percentile(a, q) { // numpy default ("linear")
    const v = a.filter(isNum).sort((x, y) => x - y); const h = (v.length - 1) * q / 100;
    const lo = Math.floor(h), hi = Math.ceil(h); return v[lo] + (h - lo) * (v[hi] - v[lo]);
  }
  function round(v, d) { const f = Math.pow(10, d); return Math.round(v * f) / f; }
  function ranks(a) { // average ranks for ties, as in pandas/scipy
    const idx = a.map((v, i) => [round(v, TIE_DECIMALS), i]).sort((x, y) => x[0] - y[0]);
    const r = new Array(a.length); let i = 0;
    while (i < idx.length) {
      let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
      const avg = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = avg; i = j + 1;
    }
    return r;
  }
  function pearson(x, y) {
    const mx = mean(x), my = mean(y); let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) { const a = x[i] - mx, b = y[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
    return sxy / Math.sqrt(sxx * syy);
  }
  function spearman(x, y) {
    const xs = [], ys = [];
    for (let i = 0; i < x.length; i++) if (isNum(x[i]) && isNum(y[i])) { xs.push(x[i]); ys.push(y[i]); }
    return pearson(ranks(xs), ranks(ys));
  }

  // matrices are arrays of columns: X[j][i] = item j, participant i
  function completeRows(X) {
    const n = X[0].length, keep = [];
    for (let i = 0; i < n; i++) if (X.every((c) => isNum(c[i]))) keep.push(i);
    return keep;
  }
  function subset(X, rows) { return X.map((c) => rows.map((i) => c[i])); }
  function rowSum(X) { return X[0].map((_, i) => X.reduce((s, c) => s + c[i], 0)); }
  function rowMean(X) { return rowSum(X).map((v) => v / X.length); }

  function cronbachAlpha(X) {
    const C = subset(X, completeRows(X)), k = C.length;
    if (k < 2 || C[0].length < 3) return NaN;
    const sumVar = C.reduce((s, c) => s + variance(c), 0), tot = variance(rowSum(C));
    return tot > 0 ? (k / (k - 1)) * (1 - sumVar / tot) : NaN;
  }
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function alphaBootstrapCI(X, rand, nBoot) {
    const C = subset(X, completeRows(X)), n = C[0].length, k = C.length, out = [];
    if (n < 3 || k < 2) return null;
    for (let b = 0; b < nBoot; b++) {
      const idx = Array.from({ length: n }, () => Math.floor(rand() * n));
      const S = C.map((c) => idx.map((i) => c[i])); const t = variance(rowSum(S));
      if (t > 0) out.push((k / (k - 1)) * (1 - S.reduce((s, c) => s + variance(c), 0) / t));
    }
    return out.length ? [percentile(out, 2.5), percentile(out, 97.5)] : null;
  }
  function itemTotal(X, names, method) {
    const C = subset(X, completeRows(X)), out = {};
    const R = (a) => a.map((v) => round(v, TIE_DECIMALS));
    C.forEach((c, j) => { const rest = rowSum(C.filter((_, k) => k !== j));
      out[names[j]] = method === "pearson" ? pearson(R(c), R(rest)) : spearman(c, rest); });
    return out;
  }

  // symmetric eigendecomposition (cyclic Jacobi), descending eigenvalues
  function eigSym(A) {
    const n = A.length, a = A.map((r) => r.slice()), V = a.map((_, i) => a.map((_, j) => (i === j ? 1 : 0)));
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
      if (off < 1e-30) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    const order = a.map((r, i) => [r[i], i]).sort((x, y) => y[0] - x[0]);
    return { values: order.map((o) => o[0]), vectors: order.map((o) => V.map((r) => r[o[1]])) };
  }
  function pca(X) { // X: columns, no missing values
    const k = X.length, n = X[0].length, Xc = X.map((c) => { const m = mean(c); return c.map((v) => v - m); });
    const C = Xc.map((a) => Xc.map((b) => a.reduce((s, v, i) => s + v * b[i], 0) / (n - 1)));
    const { values, vectors } = eigSym(C); const tot = values.reduce((s, v) => s + v, 0);
    const scores = vectors.map((v) => Xc[0].map((_, i) => v.reduce((s, w, j) => s + w * Xc[j][i], 0)));
    return { ratio: values.map((v) => v / tot), components: vectors, scores };
  }

  // ---------------------------------------------------------------- reading numbers
  // Same rules as src/chli/dataio.py: strict numbers, decimal comma only on its own,
  // usual missing-value codes are silent, anything else is missing *and reported*.
  const MISSING = new Set(["", "na", "n/a", "nan", "null", "none", "#n/a", "#na", "<na>", "-nan", ".", "-"]);
  const NUM = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/, DEC_COMMA = /^[+-]?\d+,\d+$/;
  const THOUSANDS_COMMA = /^[+-]?[1-9]\d{0,2},\d{3}$/, DEC_POINT = /^[+-]?(\d+\.\d*|\.\d+)$/;
  // Same text as dataio.number_format_warning in the Python package.
  function numberFormatWarning(column, values) {
    const texts = values.filter((v) => typeof v === "string").map((v) => v.trim());
    const thousands = texts.filter((t) => THOUSANDS_COMMA.test(t)), commaCells = texts.filter((t) => DEC_COMMA.test(t));
    // ambiguous only if every comma value could be a thousands separator (0,982 or 1,15 shows a decimal comma)
    if (thousands.length && thousands.length === commaCells.length) return `Column '${column}': ${thousands.length} cell(s) such as ${thousands.slice(0, 3).map((t) => JSON.stringify(t)).join(", ")} were read as decimal numbers; if the comma is a thousands separator, the file needs to be corrected`;
    if (texts.some((t) => DEC_COMMA.test(t)) && texts.some((t) => DEC_POINT.test(t))) return `Column '${column}': some cells use a decimal comma and others a decimal point; the values should be checked`;
    return null;
  }
  function parseNumber(v) { // -> [value, unreadable]
    if (v === null || v === undefined) return [NaN, false];
    if (typeof v === "boolean") return [NaN, true];
    if (typeof v === "number") return Number.isNaN(v) ? [NaN, false] : Number.isFinite(v) ? [v, false] : [NaN, true];
    const t = String(v).trim();
    if (MISSING.has(t.toLowerCase())) return [NaN, false];
    let f;
    if (NUM.test(t)) f = Number(t); else if (DEC_COMMA.test(t)) f = Number(t.replace(",", ".")); else return [NaN, true];
    return Number.isFinite(f) ? [f, false] : [NaN, true];
  }

  // ---------------------------------------------------------------- configuration
  const RESERVED = new Set(["CHLI", "CHLI_equal", "CHLI_pca2", "phenotype_tertile", "phenotype_median"]);
  const ID_RE = /^[A-Za-z][A-Za-z0-9_]{0,31}$/;
  const DEFAULTS = { screening_method: "spearman", screening_threshold: 0.2, max_missing_reflective: 0.3,
    pc1_min_variance: 0.2, bootstrap_n: 2000, seed: 20261020 };
  // Checks a configuration the same way as src/chli/config.py and returns a clean copy.
  function validateConfig(cfg) {
    if (!cfg || !Array.isArray(cfg.domains)) throw new Error("The configuration needs a list of 'domains'");
    if (cfg.domains.length < 2) throw new Error("The index needs at least two domains");
    const ids = new Set();
    const domains = cfg.domains.map((d) => {
      const id = d && d.id;
      if (typeof id !== "string" || !ID_RE.test(id)) throw new Error(`Domain id ${JSON.stringify(id)}: a short name with letters, digits or '_' that starts with a letter is required (for example D1)`);
      if (RESERVED.has(id) || id.startsWith("n_items_observed_")) throw new Error(`Domain id '${id}' is reserved for an output column; another id is needed`);
      if (ids.has(id)) throw new Error("Domain ids must be unique"); ids.add(id);
      const type = String(d.type ?? "").trim().toLowerCase();
      if (type !== "reflective" && type !== "formative") throw new Error(`Domain '${id}': type must be 'reflective' or 'formative', got ${JSON.stringify(d.type)}`);
      const items = (Array.isArray(d.items) ? d.items : []).map((it) => {
        const name = it && it.name;
        if (typeof name !== "string" || !name.trim()) throw new Error("Every item needs a non-empty 'name'");
        if (typeof it.column !== "string" || !it.column.trim()) throw new Error(`Item '${name}': 'column' must be the name of a column in the data`);
        const pv = it.polarity === undefined ? 1 : it.polarity;
        const p = typeof pv === "boolean" ? NaN : Number(typeof pv === "string" && !pv.trim() ? NaN : pv);
        if (p !== 1 && p !== -1) throw new Error(`Item '${name}': direction (polarity) must be 1 (higher is better) or -1 (higher is worse), got ${JSON.stringify(pv)}`);
        const transform = it.transform === undefined || it.transform === "" ? null : it.transform;
        if (transform !== null && transform !== "neg_abs_distance") throw new Error(`Item '${name}': unknown transform ${JSON.stringify(transform)} (use none or 'neg_abs_distance')`);
        let target = null;
        if (transform) {
          target = typeof it.target === "number" ? it.target : (typeof it.target === "string" && it.target.trim() ? Number(it.target) : NaN);
          if (!Number.isFinite(target)) throw new Error(`Item '${name}': transform 'neg_abs_distance' needs a numeric 'target'`);
        }
        return { name, column: it.column, polarity: p, transform, target };
      });
      if (items.length < 2) throw new Error(`Domain '${id}': needs at least two items`);
      const names = items.map((i) => i.name), dup = [...new Set(names.filter((nm, i) => names.indexOf(nm) !== i))].sort();
      if (dup.length) throw new Error(`Domain '${id}': repeated item names: ${dup.join(", ")}`);
      return { id, label: d.label || id, type, items };
    });
    const raw = cfg.settings || {}, unknown = Object.keys(raw).filter((k) => !(k in DEFAULTS)).sort();
    if (unknown.length) throw new Error(`Unknown settings: ${unknown.join(", ")} (allowed: ${Object.keys(DEFAULTS).sort().join(", ")})`);
    const st = Object.assign({}, DEFAULTS, raw);
    if (st.screening_method !== "spearman" && st.screening_method !== "pearson") throw new Error("screening_method must be 'spearman' or 'pearson'");
    for (const [k, lo, hi] of [["screening_threshold", -1, 1], ["max_missing_reflective", 0, 1], ["pc1_min_variance", 0, 1]])
      if (typeof st[k] !== "number" || !(st[k] >= lo && st[k] <= hi)) throw new Error(`${k} must be a number between ${lo} and ${hi}, got ${JSON.stringify(st[k])}`);
    for (const k of ["bootstrap_n", "seed"])
      if (!Number.isInteger(st[k]) || st[k] < 0) throw new Error(`${k} must be a whole number of at least 0, got ${JSON.stringify(st[k])}`);
    return { name: cfg.name || "CHLI", description: cfg.description || "", domains, settings: st };
  }

  // ---------------------------------------------------------------- the pipeline
  const MIN_COMPLETE = 3, FEW_COMPLETE = 30;
  function rawItem(rows, it) {
    const bad = [];
    let x = rows.map((r) => { const [f, u] = parseNumber(r[it.column]); if (u) bad.push(String(r[it.column])); return f; });
    if (it.transform === "neg_abs_distance") x = x.map((v) => -Math.abs(v - it.target));
    return { x: x.map((v) => it.polarity * v), bad };
  }

  function build(allRows, cfgIn, idColumn) {
    const cfg = validateConfig(cfgIn), st = cfg.settings;
    const cols = allRows.length ? new Set(Object.keys(allRows[0])) : new Set();
    for (const dom of cfg.domains) for (const it of dom.items)
      if (!cols.has(it.column)) throw new Error(`Column '${it.column}' (item '${it.name}') is not in the data`);
    if (idColumn !== undefined && idColumn !== null && idColumn !== "") {
      if (!cols.has(idColumn)) throw new Error(`Identifier column '${idColumn}' is not in the data`);
      if (RESERVED.has(idColumn) || cfg.domains.some((d) => d.id === idColumn) || idColumn.startsWith("n_items_observed_"))
        throw new Error(`The identifier column cannot be called '${idColumn}'; it needs another name`);
    }

    const raw = {}, unreadable = {};
    for (const dom of cfg.domains) {
      raw[dom.id] = dom.items.map((it) => {
        const { x, bad } = rawItem(allRows, it);
        if (bad.length) unreadable[it.column] = { count: bad.length, examples: bad.slice(0, 3) };
        return x;
      });
    }
    const keepIdx = allRows.map((_, i) => i).filter((i) => cfg.domains.some((d) => raw[d.id].some((x) => isNum(x[i]))));
    if (keepIdx.length < MIN_COMPLETE) throw new Error(`Only ${keepIdx.length} rows have any usable value; the index needs more participants`);
    const n = keepIdx.length, rand = mulberry32(st.seed);
    const rep = { config: cfg.name, n, n_rows_in_file: allRows.length, n_rows_without_data: allRows.length - n,
      settings: st, unreadable_cells: unreadable, warnings: [], domains: {} };
    for (const [col, u] of Object.entries(unreadable))
      rep.warnings.push(`Column '${col}': ${u.count} cell(s) could not be read as numbers and were treated as missing (e.g. ${u.examples.map((e) => JSON.stringify(e)).join(", ")})`);
    const seenCols = new Set();
    for (const dom of cfg.domains) for (const it of dom.items) {
      if (seenCols.has(it.column)) continue; seenCols.add(it.column);
      const w = numberFormatWarning(it.column, allRows.map((r) => r[it.column])); if (w) rep.warnings.push(w);
    }
    if (rep.n_rows_without_data) rep.warnings.push(`${rep.n_rows_without_data} row(s) have no usable value and were left out`);
    const finalX = {}, finalNames = {};

    for (const dom of cfg.domains) {
      let names = dom.items.map((it) => it.name), X = raw[dom.id].map((x) => keepIdx.map((i) => x[i]));
      const miss = {}, frac = {};
      names.forEach((nm, j) => { frac[nm] = X[j].filter((v) => !isNum(v)).length / n; miss[nm] = round(100 * frac[nm], 1); });
      const d = { label: dom.label, type: dom.type, items_considered: names.slice(), missing_pct: miss,
        excluded_unusable: [], excluded_missing: [], excluded_screening: [], item_total: {}, warnings: [] };
      // no variation is judged on rounded values (same rule as the Python package)
      const usable = X.map((c) => c.filter(isNum).length >= MIN_COMPLETE && new Set(c.filter(isNum).map((v) => round(v, TIE_DECIMALS))).size >= 2);
      d.excluded_unusable = names.filter((_, j) => !usable[j]);
      X = X.filter((_, j) => usable[j]).map(zscore); names = names.filter((_, j) => usable[j]);
      if (dom.type === "reflective") {
        const keep = names.map((nm) => frac[nm] <= st.max_missing_reflective);
        d.excluded_missing = names.filter((_, j) => !keep[j]); X = X.filter((_, j) => keep[j]); names = names.filter((_, j) => keep[j]);
      } else {
        d.warnings = names.filter((nm) => frac[nm] > st.max_missing_reflective).map((nm) => `${nm}: ${miss[nm]}% missing (formative domain, kept)`);
      }
      if (X.length < 2) throw new Error(`Domain ${dom.id}: fewer than two usable items (without variation: ${d.excluded_unusable.join(", ") || "none"}; too much missing: ${d.excluded_missing.join(", ") || "none"})`);
      if (dom.type === "reflective") {
        const ncc = completeRows(X).length;
        if (ncc < MIN_COMPLETE) throw new Error(`Domain ${dom.id}: only ${ncc} participant(s) have every item, so the reliability screening cannot be done. The column mapping should be checked, or the item with the most missing values removed.`);
        d.alpha_before_screening = cronbachAlpha(X);
        d.item_total = itemTotal(X, names, st.screening_method);
        const keep = names.map((nm) => d.item_total[nm] >= st.screening_threshold);
        d.excluded_screening = names.filter((_, j) => !keep[j]); X = X.filter((_, j) => keep[j]); names = names.filter((_, j) => keep[j]);
        if (X.length < 2) throw new Error(`Domain ${dom.id}: fewer than two items left after the reliability screening`);
      } else {
        d.alpha_before_screening = cronbachAlpha(X);
      }
      d.items_retained = names.slice(); d.n_complete_cases = completeRows(X).length;
      if (d.n_complete_cases < FEW_COMPLETE) d.warnings.push(`only ${d.n_complete_cases} complete cases: reliability estimates are unstable`);
      finalX[dom.id] = X; finalNames[dom.id] = names; rep.domains[dom.id] = d;
    }
    for (const dom of cfg.domains) {
      const d = rep.domains[dom.id], X = finalX[dom.id];
      d.alpha = cronbachAlpha(X);
      d.alpha_ci95 = st.bootstrap_n ? alphaBootstrapCI(X, rand, st.bootstrap_n) : null;
      d.omega = null; // computed in the Python package only
    }

    const S = {}, nobs = {};
    for (const dom of cfg.domains) {
      const d = rep.domains[dom.id], X = finalX[dom.id], names = finalNames[dom.id];
      nobs[dom.id] = X[0].map((_, i) => X.reduce((s, c) => s + (isNum(c[i]) ? 1 : 0), 0));
      const Xi = X.map((c) => { const m = median(c); return c.map((v) => (isNum(v) ? v : m)); });
      if (dom.type === "reflective") {
        const p = pca(Xi); let s = p.scores[0], comp = p.components[0];
        if (pearson(s, rowMean(Xi)) < 0) { s = s.map((v) => -v); comp = comp.map((v) => -v); }
        S[dom.id] = zscore(s); d.aggregation = "PC1"; d.pc1_variance = p.ratio[0];
        d.pc1_meets_threshold = p.ratio[0] >= st.pc1_min_variance;
        d.weights = Object.fromEntries(names.map((nm, j) => [nm, comp[j]]).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])));
      } else {
        S[dom.id] = zscore(rowMean(Xi)); d.aggregation = "mean of z-scores (formative)"; d.pc1_variance = null;
        d.weights = Object.fromEntries(names.map((nm) => [nm, 1 / names.length]));
        d.descriptive_pca_variance = pca(Xi).ratio;
      }
      d.participants_without_any_item = nobs[dom.id].filter((v) => v === 0).length;
    }

    const ids = cfg.domains.map((d) => d.id), CH = ids.map((k) => S[k]);
    const equal = rowMean(CH), chliZ = zscore(equal);
    const p2 = pca(CH); let L = p2.components[0], pca2 = p2.scores[0];
    if (L.reduce((s, v) => s + v, 0) < 0) { L = L.map((v) => -v); pca2 = pca2.map((v) => -v); }
    const p33 = percentile(equal, 100 / 3), p67 = percentile(equal, 200 / 3), med = median(equal);
    const R = (v) => round(v, TIE_DECIMALS), r33 = R(p33), r67 = R(p67), rmed = R(med);
    const tert = equal.map((v) => (R(v) >= r67 ? "Healthy" : R(v) <= r33 ? "Non-healthy" : "Intermediate"));
    const binl = equal.map((v) => (R(v) >= rmed ? "Healthy" : "Non-healthy"));
    if (r33 === r67) rep.warnings.push("The two tertile cut-points are equal (many tied index values); tertile labels are not informative in this cohort");
    const count = (a) => a.reduce((o, v) => ((o[v] = (o[v] || 0) + 1), o), {});

    rep.second_level_pca = { pc1_variance: p2.ratio[0], loadings: Object.fromEntries(ids.map((k, j) => [k, L[j]])),
      spearman_with_CHLI_equal: spearman(equal, pca2) };
    rep.inter_domain_spearman = {};
    ids.forEach((a, i) => ids.slice(i + 1).forEach((b) => { rep.inter_domain_spearman[`${a}-${b}`] = spearman(S[a], S[b]); }));
    rep.leave_one_domain_out = {};
    for (const k of ids) {
      const rest = zscore(rowMean(ids.filter((x) => x !== k).map((x) => S[x])));
      rep.leave_one_domain_out[k] = { rho_with_full_index: spearman(chliZ, rest), rho_with_excluded_domain: spearman(rest, S[k]) };
    }
    rep.CHLI_equal = { mean: mean(equal), sd: sd(equal), median: med, min: Math.min(...equal), max: Math.max(...equal) };
    rep.phenotypes = { tertile_cutpoints: { P33: p33, P67: p67 }, tertile_counts: count(tert), median_cutpoint: med, median_counts: count(binl) };

    const pos = new Map(keepIdx.map((r, i) => [r, i]));
    const scores = allRows.map((row, r) => {
      const o = {}, i = pos.get(r), has = i !== undefined;
      if (idColumn) o[idColumn] = row[idColumn];
      for (const k of ids) o[k] = has ? S[k][i] : NaN;
      Object.assign(o, has ? { CHLI_equal: equal[i], CHLI: chliZ[i], CHLI_pca2: pca2[i], phenotype_tertile: tert[i], phenotype_median: binl[i] }
        : { CHLI_equal: NaN, CHLI: NaN, CHLI_pca2: NaN, phenotype_tertile: "", phenotype_median: "" });
      for (const k of ids) o[`n_items_observed_${k}`] = has ? nobs[k][i] : 0;
      return o;
    });
    return { scores, report: rep };
  }

  // ---------------------------------------------------------------- CSV
  // Same rules as src/chli/dataio.py (read_csv_text).
  function parseCSV(text) {
    text = text.replace(/^﻿/, "");
    if (!text.trim()) throw new Error("The file is empty");
    const first = text.split(/\r\n|\r|\n/)[0];
    const delim = [",", ";", "\t"].map((d) => [d, first.split(d).length]).sort((a, b) => b[1] - a[1])[0][0];
    // Same rules as Python's csv module with strict=True: a quotation mark opens a quoted
    // value only as the first character of a cell, elsewhere it is an ordinary character.
    const rows = []; let row = [], field = "", q = false, fresh = true, closed = false, line = 1;
    for (let i = 0; i < text.length; i++) {
      const ch = text[i];
      if (q) {
        if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
        else if (ch === '"') { q = false; closed = true; }
        else { if (ch === "\n" || (ch === "\r" && text[i + 1] !== "\n")) line++; field += ch; }
      } else if (ch === delim) { row.push(field); field = ""; fresh = true; closed = false; }
      else if (ch === "\n" || ch === "\r") { if (ch === "\r" && text[i + 1] === "\n") i++; row.push(field); rows.push(row); row = []; field = ""; fresh = true; closed = false; line++; }
      else if (closed) throw new Error(`Line ${line}: after a closing quotation mark (") a separator or the end of the line is expected.`);
      else if (ch === '"' && fresh) { q = true; fresh = false; }
      else { field += ch; fresh = false; }
    }
    if (q) throw new Error('The file ends inside a value that opens with a quotation mark (") and is never closed. A quotation mark may be missing or there may be one too many.');
    if (field !== "" || row.length || closed) { row.push(field); rows.push(row); }
    const head = rows.shift().map((h) => h.trim()), named = head.filter((h) => h);
    const dup = [...new Set(named.filter((h, i) => named.indexOf(h) !== i))].sort();
    if (dup.length) throw new Error(`Repeated column names: ${dup.join(", ")}. Each column needs a unique name.`);
    const out = []; let blank = 0;
    rows.forEach((r, k) => {
      if (r.every((v) => !v.trim())) { if (r.length > 1) blank++; return; }
      if (r.length > head.length && r.slice(head.length).some((v) => v.trim()))
        throw new Error(`Row ${k + 2} has ${r.length} cells but the header has ${head.length}. The row may have an extra separator or an unquoted comma inside a value.`);
      out.push(Object.fromEntries(head.map((h, j) => [h, r[j] === undefined ? "" : r[j]])));
    });
    return { columns: head, rows: out, delimiter: delim, blank_rows: blank };
  }
  // Same rules as dataio.decode_bytes in the Python package: no guessing of the encoding.
  function decodeBytes(bytes, encoding) {
    const enc = String(encoding || "utf-8").toLowerCase();
    let b = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
    if (enc === "utf-8" || enc === "utf8") {
      if (b[0] === 0xef && b[1] === 0xbb && b[2] === 0xbf) b = b.subarray(3);
      try { return new TextDecoder("utf-8", { fatal: true }).decode(b); }
      catch (e) {
        let pos = 0, line = 1;
        while (pos < b.length) {   // find the first invalid byte to report it
          const c = b[pos]; let len = c < 0x80 ? 1 : c >= 0xc2 && c < 0xe0 ? 2 : c >= 0xe0 && c < 0xf0 ? 3 : c >= 0xf0 && c < 0xf5 ? 4 : 0;
          let ok = len > 0 && pos + len <= b.length;
          if (ok) { try { new TextDecoder("utf-8", { fatal: true }).decode(b.subarray(pos, pos + len)); } catch (e2) { ok = false; } }
          if (!ok) break;
          if (c === 0x0a) line++; pos += len;
        }
        const hex = pos < b.length ? b[pos].toString(16).toUpperCase().padStart(2, "0") : "??";
        throw new Error(`The file is not UTF-8 text (byte 0x${hex} on line ${line}). It can be saved again as 'CSV UTF-8', or the encoding can be set to windows-1252.`);
      }
    }
    if (enc === "windows-1252" || enc === "cp1252") return new TextDecoder("windows-1252").decode(b);
    throw new Error(`Unknown encoding '${encoding}' (use utf-8 or windows-1252)`);
  }

  function toCSV(objs) {
    if (!objs.length) return "";
    const cols = Object.keys(objs[0]);
    const esc = (v) => { const s = typeof v === "number" ? (Number.isFinite(v) ? String(v) : "") : String(v ?? ""); return /[",\n;]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
    return [cols.join(","), ...objs.map((o) => cols.map((c) => esc(o[c])).join(","))].join("\n") + "\n";
  }

  return { build, validateConfig, parseCSV, parseNumber, decodeBytes, numberFormatWarning, toCSV, DEFAULTS, stats: { mean, sd, median, zscore, percentile, spearman, cronbachAlpha, itemTotal, pca, eigSym } };
});
