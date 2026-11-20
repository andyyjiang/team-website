// Local trial: runs the login service on http://localhost:8787 with a fake
// GitHub that starts from this repository's data/curation.json and
// data/feedback.json (nothing is sent to the real GitHub).
//
//   npm run dev        then serve the site:  python3 -m http.server 8000  (in the repository folder)
//   open http://localhost:8000/mpstruc.html?login-service=http://localhost:8787#/curators
//   setup code: open-sesame-42
//
// http://localhost:8788/calls lists what reached the fake GitHub; POST /reset-limits clears sign-in failure counts.
import { createServer } from "node:http";
import { readFileSync, existsSync } from "node:fs";
import { startWorker } from "./harness.mjs";

const root = new URL("../../", import.meta.url);
const seed = {};
for (const f of ["data/curation.json", "data/feedback.json"]) {
  const p = new URL(f, root);
  if (existsSync(p)) seed[f] = readFileSync(p, "utf8");
}
const port = Number(process.env.PORT || 8787);
const w = await startWorker({ seed, port });
createServer(async (req, res) => {
  res.setHeader("Content-Type", "application/json");
  if (req.url === "/calls") return res.end(JSON.stringify(w.github.calls));
  if (req.url === "/reset-limits" && req.method === "POST") {
    const kv = await w.kv();
    const { keys } = await kv.list({ prefix: "fail:" });
    await Promise.all(keys.map((k) => kv.delete(k.name)));
    return res.end(JSON.stringify({ cleared: keys.length }));
  }
  res.statusCode = 404;
  res.end("{}");
}).listen(port + 1, "127.0.0.1");
console.log(`login service: ${w.base}  (control: http://127.0.0.1:${port + 1})`);
