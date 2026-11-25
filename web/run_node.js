// Run the browser pipeline from the command line (used by tests/test_js_parity.py).
// node web/run_node.js <data.csv> <config.json> [utf-8|windows-1252] [id column]  -> prints {scores, report} as JSON
const fs = require("fs");
const CHLI = require("./chli-core.js");
const [, , dataPath, cfgPath, encoding, idColumn] = process.argv;
try {
  const { rows } = CHLI.parseCSV(CHLI.decodeBytes(fs.readFileSync(dataPath), encoding || "utf-8"));
  const cfg = JSON.parse(fs.readFileSync(cfgPath, "utf8"));
  process.stdout.write(JSON.stringify(CHLI.build(rows, cfg, idColumn || undefined)));
} catch (e) {
  process.stdout.write(JSON.stringify({ error: e.message }));
}
