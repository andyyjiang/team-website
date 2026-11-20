// Runs worker.js in Miniflare (the Workers runtime) with an in-memory D1 database and a
// fake api.github.com, for the tests and for local trials of the site.
import { Miniflare } from "miniflare";
import { fileURLToPath } from "node:url";
import { fakeGitHub } from "./fake-github.mjs";

export const SECRETS = { GITHUB_TOKEN: "bot-token", SESSION_SECRET: "test-session-secret-0123456789abcdef", SETUP_CODE: "open-sesame-42" };

export async function startWorker({ seed = {}, port, secrets = SECRETS, githubAccepts = secrets.GITHUB_TOKEN, d1 = true, vars = {} } = {}) {
  const github = fakeGitHub(seed, { token: githubAccepts });
  const mf = new Miniflare({
    modules: true,
    scriptPath: fileURLToPath(new URL("../worker.js", import.meta.url)),
    compatibilityDate: "2024-11-01",
    ...(d1 ? { d1Databases: ["DB"] } : {}),
    bindings: { ...secrets, ALLOWED_ORIGIN: "https://andyyjiang.github.io", ...vars },
    outboundService: github.handler,
    ...(port ? { host: "127.0.0.1", port } : {}),
  });
  const base = port ? (await mf.ready).toString().replace(/\/$/, "") : "http://worker.test";
  // call(path, { method, body, token, origin, ip })
  async function call(path, { method = "GET", body, token, origin = "https://andyyjiang.github.io", ip = "203.0.113.1", raw } = {}) {
    const headers = { "CF-Connecting-IP": ip };
    if (origin) headers.Origin = origin;
    if (token) headers.Authorization = `Bearer ${token}`;
    if (body !== undefined) headers["Content-Type"] = "application/json";
    const res = await mf.dispatchFetch(base + path, { method, headers, body: raw ?? (body === undefined ? undefined : JSON.stringify(body)) });
    const text = await res.text();
    let data = null;
    try { data = JSON.parse(text); } catch { data = text; }
    return { status: res.status, data, headers: res.headers };
  }
  // sql(query, ...params): rows, straight from the Worker's database.
  async function sql(query, ...params) {
    const db = await mf.getD1Database("DB");
    return (await db.prepare(query).bind(...params).all()).results;
  }
  return { mf, call, github, base, sql };
}
