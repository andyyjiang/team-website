// A tiny stand-in for the parts of the GitHub REST API the login service uses,
// for tests and local trials. Every request is recorded in `calls`.
import { createHash } from "node:crypto";

export function fakeGitHub(seed = {}, { token = "bot-token" } = {}) {
  const files = new Map(); // path -> { text, sha }
  const blobs = new Map(); // sha -> text
  const sha = (t) => createHash("sha1").update(t).digest("hex");
  const put = (path, text) => { const s = sha(text + path); files.set(path, { text, sha: s }); blobs.set(s, text); return s; };
  for (const [p, t] of Object.entries(seed)) put(p, t);
  const calls = [];
  const reply = (status, data) => new Response(JSON.stringify(data), { status, headers: { "Content-Type": "application/json" } });

  async function handler(request) {
    const url = new URL(request.url);
    const call = { method: request.method, host: url.host, path: url.pathname, search: url.search,
      authorization: request.headers.get("Authorization"), userAgent: request.headers.get("User-Agent") };
    if (request.method === "PUT") call.body = JSON.parse(await request.text());
    calls.push(call);
    if (url.host !== "api.github.com") return reply(599, { message: "unexpected host" });
    if (call.authorization !== `Bearer ${token}`) return reply(401, { message: "Bad credentials" });
    let m = /^\/repos\/andyyjiang\/team-website\/contents\/(.+)$/.exec(url.pathname);
    if (m && request.method === "GET") {
      const f = files.get(m[1]);
      if (!f) return reply(404, { message: "Not Found" });
      return reply(200, { path: m[1], sha: f.sha, encoding: "base64", content: Buffer.from(f.text).toString("base64") });
    }
    if (m && request.method === "PUT") {
      const f = files.get(m[1]);
      if ((f ? f.sha : undefined) !== call.body.sha) return reply(409, { message: `${m[1]} does not match ${call.body.sha}` });
      const s = put(m[1], Buffer.from(call.body.content, "base64").toString("utf8"));
      return reply(f ? 200 : 201, { content: { path: m[1], sha: s }, commit: { message: call.body.message } });
    }
    m = /^\/repos\/andyyjiang\/team-website\/git\/blobs\/([0-9a-f]{40})$/.exec(url.pathname);
    if (m && blobs.has(m[1])) return reply(200, { sha: m[1], encoding: "base64", content: Buffer.from(blobs.get(m[1])).toString("base64") });
    return reply(404, { message: "Not Found" });
  }
  return { handler, calls, files, put };
}
