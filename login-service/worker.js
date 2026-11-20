// mpstruc login service: username and password sign-in for curators.
//
// A Cloudflare Worker (one file, no build step). It keeps the curator
// accounts in a KV namespace (binding USERS), checks passwords, hands out
// signed session tokens, and forwards a short allow-list of GitHub API calls
// with ONE bot token, so curators do not need GitHub accounts.
//
// Secrets (Worker settings > Variables and Secrets, type "Secret"):
//   GITHUB_TOKEN    fine-grained GitHub token: Contents read/write on andyyjiang/team-website only
//   SESSION_SECRET  a long random string that signs session tokens
//   SETUP_CODE      a one-off code that creates the first administrator
// Binding: KV namespace as USERS.
// Optional variable: ALLOWED_ORIGIN (default https://andyyjiang.github.io).
//
// Passwords and tokens are never logged or returned.

const REPO = "andyyjiang/team-website";
const BRANCH = "main";
const GITHUB_API = "https://api.github.com";
const DEFAULT_ORIGIN = "https://andyyjiang.github.io";
// Cloudflare's WebCrypto allows at most 100,000 PBKDF2 iterations.
const PBKDF2_ITERATIONS = 100000;
const SHORT_SESSION = 12 * 3600;        // seconds
const LONG_SESSION = 30 * 24 * 3600;    // "Keep me signed in"
const MAX_BODY = 2 * 1024 * 1024;       // bytes
const MIN_PASSWORD = 10, MAX_PASSWORD = 200;
const FAIL_LIMIT = 10, FAIL_WINDOW = 15 * 60; // failures per username and per IP, seconds
const USERNAME_RE = /^[a-z0-9._-]{3,32}$/;
const ROLES = ["admin", "curator"];

// The only GitHub calls the site may make through this service.
const FILES = ["data/curation.json", "data/feedback.json"];
function allowedGitHubCall(method, path, search) {
  for (const f of FILES) {
    const p = `/repos/${REPO}/contents/${f}`;
    if (path === p && method === "GET" && search === `?ref=${BRANCH}`) return "read";
    if (path === p && method === "PUT" && search === "") return "write";
  }
  if (method === "GET" && search === "" && new RegExp(`^/repos/${REPO}/git/blobs/[0-9a-f]{40}$`).test(path)) return "read";
  return null;
}

class HttpError extends Error {
  constructor(status, message, extra = {}) { super(message); this.status = status; this.extra = extra; }
}

// ---- CORS ----
function allowedOrigin(origin, env) {
  if (!origin) return false;
  const site = (env.ALLOWED_ORIGIN || DEFAULT_ORIGIN).replace(/\/+$/, "");
  if (origin === site) return true;
  return /^http:\/\/(localhost|127\.0\.0\.1)(:\d{1,5})?$/.test(origin);
}
function corsHeaders(origin) {
  return {
    "Access-Control-Allow-Origin": origin,
    "Access-Control-Allow-Methods": "GET, POST, PUT, DELETE, OPTIONS",
    "Access-Control-Allow-Headers": "Authorization, Content-Type, Accept, X-GitHub-Api-Version",
    "Access-Control-Max-Age": "600",
    Vary: "Origin",
  };
}
function json(data, status, origin) {
  return new Response(JSON.stringify(data), {
    status,
    headers: { "Content-Type": "application/json; charset=utf-8", "Cache-Control": "no-store", ...(origin ? corsHeaders(origin) : {}) },
  });
}

// ---- encoding helpers ----
const enc = new TextEncoder();
function b64url(bytes) {
  let bin = "";
  for (const b of bytes) bin += String.fromCharCode(b);
  return btoa(bin).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}
function fromB64url(s) {
  const b = atob(s.replace(/-/g, "+").replace(/_/g, "/") + "===".slice((s.length + 3) % 4));
  return Uint8Array.from(b, (c) => c.charCodeAt(0));
}
function equalBytes(a, b) { // constant time for equal lengths
  if (a.length !== b.length) return false;
  let diff = 0;
  for (let i = 0; i < a.length; i++) diff |= a[i] ^ b[i];
  return diff === 0;
}

// ---- passwords ----
export async function hashPassword(password, saltB64, iterations = PBKDF2_ITERATIONS) {
  const salt = saltB64 ? fromB64url(saltB64) : crypto.getRandomValues(new Uint8Array(16));
  const key = await crypto.subtle.importKey("raw", enc.encode(password), "PBKDF2", false, ["deriveBits"]);
  const bits = await crypto.subtle.deriveBits({ name: "PBKDF2", hash: "SHA-256", salt, iterations }, key, 256);
  return { hash: b64url(new Uint8Array(bits)), salt: b64url(salt), iterations };
}
async function checkPassword(user, password) {
  const { hash } = await hashPassword(password, user.salt, user.iterations);
  return equalBytes(fromB64url(hash), fromB64url(user.hash));
}
// Used when the username does not exist, so the answer takes as long either way.
const DUMMY = { salt: "AAAAAAAAAAAAAAAAAAAAAA", iterations: PBKDF2_ITERATIONS, hash: "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" };

function checkNewPassword(pw) {
  if (typeof pw !== "string" || pw.length < MIN_PASSWORD) throw new HttpError(400, `The password must be at least ${MIN_PASSWORD} characters long.`);
  if (pw.length > MAX_PASSWORD) throw new HttpError(400, `The password can be at most ${MAX_PASSWORD} characters long.`);
}

// ---- sessions ----
async function hmacKey(env) {
  if (!env.SESSION_SECRET || env.SESSION_SECRET.length < 16) throw new HttpError(500, "The login service is not set up yet (SESSION_SECRET is missing).");
  return crypto.subtle.importKey("raw", enc.encode(env.SESSION_SECRET), { name: "HMAC", hash: "SHA-256" }, false, ["sign"]);
}
async function sign(env, payload) {
  const body = b64url(enc.encode(JSON.stringify(payload)));
  const sig = new Uint8Array(await crypto.subtle.sign("HMAC", await hmacKey(env), enc.encode(body)));
  return `${body}.${b64url(sig)}`;
}
async function readToken(env, token) {
  const m = /^([A-Za-z0-9_-]{10,2000})\.([A-Za-z0-9_-]{43})$/.exec(token || "");
  if (!m) return null;
  const want = new Uint8Array(await crypto.subtle.sign("HMAC", await hmacKey(env), enc.encode(m[1])));
  if (!equalBytes(want, fromB64url(m[2]))) return null;
  let p;
  try { p = JSON.parse(new TextDecoder().decode(fromB64url(m[1]))); } catch { return null; }
  return p && typeof p.u === "string" && typeof p.exp === "number" ? { ...p, sig: m[2] } : null;
}
const now = () => Math.floor(Date.now() / 1000);
const nonce = () => b64url(crypto.getRandomValues(new Uint8Array(9))); // makes every token distinct, so one can be signed out alone

const SIGNED_OUT = "Your sign-in has ended. Please sign in again.";
// The signed-in user for this request, checked against the account as stored now:
// a disabled account, a changed password or a reset ends every older session.
async function session(request, env) {
  const auth = request.headers.get("Authorization") || "";
  const token = auth.startsWith("Bearer ") ? auth.slice(7).trim() : "";
  const p = await readToken(env, token);
  if (!p) throw new HttpError(401, SIGNED_OUT);
  if (p.exp <= now()) throw new HttpError(401, "Your sign-in has expired. Please sign in again.");
  if (await env.USERS.get(`revoked:${p.sig}`)) throw new HttpError(401, SIGNED_OUT);
  const user = await getUser(env, p.u);
  if (!user) throw new HttpError(401, "This account no longer exists. Ask an administrator.");
  if (user.disabled) throw new HttpError(401, "This account has been disabled. Ask an administrator to enable it again.");
  if ((user.sessionVersion || 0) !== p.v) throw new HttpError(401, "Your password was changed or reset, so you have been signed out. Please sign in again.");
  return { user, payload: p };
}
async function requireAdmin(request, env) {
  const s = await session(request, env);
  if (s.user.role !== "admin") throw new HttpError(403, "Only administrators can manage curators.");
  return s;
}

// ---- users ----
const userKey = (name) => `user:${name}`;
async function getUser(env, name) {
  if (typeof name !== "string" || !USERNAME_RE.test(name)) return null;
  const raw = await env.USERS.get(userKey(name));
  return raw ? JSON.parse(raw) : null;
}
const putUser = (env, u) => env.USERS.put(userKey(u.username), JSON.stringify(u));
async function allUsers(env) {
  const out = [];
  let cursor;
  do {
    const page = await env.USERS.list({ prefix: "user:", cursor });
    for (const k of page.keys) {
      const raw = await env.USERS.get(k.name);
      if (raw) out.push(JSON.parse(raw));
    }
    cursor = page.list_complete ? undefined : page.cursor;
  } while (cursor);
  return out.sort((a, b) => a.username.localeCompare(b.username));
}
const publicUser = (u) => ({ username: u.username, displayName: u.displayName, role: u.role, disabled: Boolean(u.disabled),
  createdAt: u.createdAt, createdBy: u.createdBy });
function cleanName(v) {
  const s = typeof v === "string" ? v.trim().toLowerCase() : "";
  if (!USERNAME_RE.test(s)) throw new HttpError(400, "A username is 3 to 32 characters: lowercase letters, digits, dot, dash or underscore.");
  return s;
}
function cleanDisplay(v, fallback) {
  const s = typeof v === "string" ? v.replace(/\s+/g, " ").trim() : "";
  if (s.length > 80) throw new HttpError(400, "The display name can be at most 80 characters.");
  return s || fallback;
}
async function newUser(env, { username, displayName, password, role }, createdBy) {
  const name = cleanName(username);
  if (!ROLES.includes(role)) throw new HttpError(400, "The role must be admin or curator.");
  checkNewPassword(password);
  if (await getUser(env, name)) throw new HttpError(409, `There is already an account called ${name}.`);
  const u = { username: name, displayName: cleanDisplay(displayName, name), role, ...(await hashPassword(password)),
    createdAt: new Date().toISOString(), createdBy, disabled: false, sessionVersion: 0 };
  await putUser(env, u);
  return u;
}
async function activeAdmins(env) { return (await allUsers(env)).filter((u) => u.role === "admin" && !u.disabled); }

// ---- rate limit (failed sign-ins) ----
async function failures(env, key) {
  const raw = await env.USERS.get(key);
  if (!raw) return { n: 0, reset: 0 };
  const r = JSON.parse(raw);
  return r.reset > now() ? r : { n: 0, reset: 0 };
}
async function addFailure(env, key) {
  const r = await failures(env, key);
  const reset = r.reset || now() + FAIL_WINDOW;
  await env.USERS.put(key, JSON.stringify({ n: r.n + 1, reset }), { expirationTtl: Math.max(60, reset - now()) });
}
function tooMany(reset) {
  const mins = Math.max(1, Math.ceil((reset - now()) / 60));
  return new HttpError(429, `Too many failed sign-in attempts. Try again in ${mins} ${mins === 1 ? "minute" : "minutes"}.`, { retryMinutes: mins });
}

// ---- request body ----
async function body(request) {
  const len = Number(request.headers.get("Content-Length") || 0);
  if (len > MAX_BODY) throw new HttpError(413, "This change is too large to save (over 2 MB).");
  const text = await request.text();
  if (enc.encode(text).length > MAX_BODY) throw new HttpError(413, "This change is too large to save (over 2 MB).");
  try {
    const v = JSON.parse(text || "{}");
    if (!v || typeof v !== "object" || Array.isArray(v)) throw new Error();
    return v;
  } catch { throw new HttpError(400, "The request was not valid JSON."); }
}

// ---- handlers ----
async function login(request, env) {
  const b = await body(request);
  const name = typeof b.username === "string" ? b.username.trim().toLowerCase() : "";
  const password = typeof b.password === "string" ? b.password : "";
  if (!name || !password) throw new HttpError(400, "Enter your username and password.");
  const ip = request.headers.get("CF-Connecting-IP") || "unknown";
  const uKey = `fail:user:${name.slice(0, 64)}`, ipKey = `fail:ip:${ip}`;
  const [fu, fi] = await Promise.all([failures(env, uKey), failures(env, ipKey)]);
  if (fu.n >= FAIL_LIMIT) throw tooMany(fu.reset);
  if (fi.n >= FAIL_LIMIT) throw tooMany(fi.reset);
  const user = await getUser(env, name);
  const ok = await checkPassword(user || DUMMY, password.slice(0, MAX_PASSWORD));
  if (!user || !ok) {
    await Promise.all([addFailure(env, uKey), addFailure(env, ipKey)]);
    throw new HttpError(401, "Wrong username or password.");
  }
  if (user.disabled) throw new HttpError(403, "This account has been disabled. Ask an administrator to enable it again.");
  await env.USERS.delete(uKey);
  const exp = now() + (b.keep === true ? LONG_SESSION : SHORT_SESSION);
  const token = await sign(env, { u: user.username, r: user.role, v: user.sessionVersion || 0, exp, n: nonce() });
  return { token, expires: new Date(exp * 1000).toISOString(), user: publicUser(user) };
}

async function logout(request, env) {
  try {
    const { payload } = await session(request, env);
    const ttl = payload.exp - now();
    if (ttl > 0) await env.USERS.put(`revoked:${payload.sig}`, "1", { expirationTtl: Math.max(60, ttl) });
  } catch { /* already signed out: nothing to do */ }
  return { ok: true };
}

async function setup(request, env) {
  if ((await activeAdmins(env)).length) throw new HttpError(403, "Setup is already done. Sign in with your username and password.");
  const b = await body(request);
  const ip = request.headers.get("CF-Connecting-IP") || "unknown";
  const ipKey = `fail:ip:${ip}`;
  const fi = await failures(env, ipKey);
  if (fi.n >= FAIL_LIMIT) throw tooMany(fi.reset);
  const code = typeof b.setupCode === "string" ? b.setupCode.trim() : "";
  if (!env.SETUP_CODE) throw new HttpError(500, "The login service is not set up yet (SETUP_CODE is missing).");
  const digest = async (s) => new Uint8Array(await crypto.subtle.digest("SHA-256", enc.encode(s)));
  if (!code || !equalBytes(await digest(code), await digest(String(env.SETUP_CODE).trim()))) {
    await addFailure(env, ipKey);
    throw new HttpError(401, "That setup code is not right.");
  }
  const u = await newUser(env, { ...b, role: "admin" }, "setup");
  return { ok: true, user: publicUser(u) };
}

async function changeOwnPassword(request, env) {
  const { user } = await session(request, env);
  const b = await body(request);
  if (!(await checkPassword(user, String(b.oldPassword || "").slice(0, MAX_PASSWORD)))) throw new HttpError(400, "Your current password is not right.");
  checkNewPassword(b.newPassword);
  Object.assign(user, await hashPassword(b.newPassword), { sessionVersion: (user.sessionVersion || 0) + 1 });
  await putUser(env, user);
  // A fresh session for this browser; every other one has ended.
  const p = (await readToken(env, (request.headers.get("Authorization") || "").slice(7).trim()));
  const token = await sign(env, { u: user.username, r: user.role, v: user.sessionVersion, exp: p.exp, n: nonce() });
  return { ok: true, token, user: publicUser(user) };
}

async function manageUser(request, env, name, action) {
  const { user: me } = await requireAdmin(request, env);
  const u = await getUser(env, name);
  if (!u) throw new HttpError(404, `There is no account called ${name}.`);
  const lastAdmin = async () => u.role === "admin" && !u.disabled && (await activeAdmins(env)).length <= 1;
  if (action === "reset") {
    const b = await body(request);
    checkNewPassword(b.password);
    Object.assign(u, await hashPassword(b.password), { sessionVersion: (u.sessionVersion || 0) + 1 });
  } else if (action === "disable") {
    if (u.username === me.username) throw new HttpError(400, "You cannot disable your own account.");
    if (await lastAdmin()) throw new HttpError(400, "This is the last administrator; it cannot be disabled.");
    u.disabled = true;
    u.sessionVersion = (u.sessionVersion || 0) + 1;
  } else if (action === "enable") {
    u.disabled = false;
  } else if (action === "delete") {
    if (u.username === me.username) throw new HttpError(400, "You cannot delete your own account.");
    if (await lastAdmin()) throw new HttpError(400, "This is the last administrator; it cannot be deleted.");
    await env.USERS.delete(userKey(u.username));
    return { ok: true };
  }
  await putUser(env, u);
  return { ok: true, user: publicUser(u) };
}

// GitHub proxy: /gh/<GitHub API path>, allow-list only.
async function proxy(request, env, url) {
  const path = url.pathname.slice(3);
  const kind = allowedGitHubCall(request.method, path, url.search);
  if (!kind) throw new HttpError(403, "The login service does not allow this request.");
  const { user } = await session(request, env);
  if (!ROLES.includes(user.role)) throw new HttpError(403, "This account cannot save changes.");
  if (!env.GITHUB_TOKEN) throw new HttpError(500, "The login service is not set up yet (GITHUB_TOKEN is missing).");
  let payload;
  if (kind === "write") {
    const b = await body(request);
    if (typeof b.content !== "string" || typeof b.message !== "string" || (b.sha !== undefined && typeof b.sha !== "string")
        || (b.branch !== undefined && b.branch !== BRANCH)) throw new HttpError(400, "This change is not in the expected form.");
    // The commit is made by the bot; the curator is named in the message.
    const rest = b.message.replace(/^Curator [^:\n]{1,80}:\s*/, "").replace(/[\r\n]+/g, " ").slice(0, 900);
    payload = { message: `Curator ${user.username}: ${rest}`, content: b.content, branch: BRANCH, ...(b.sha ? { sha: b.sha } : {}) };
  }
  let r;
  try {
    r = await fetch(GITHUB_API + path + url.search, {
      method: request.method,
      headers: { Accept: "application/vnd.github+json", Authorization: `Bearer ${env.GITHUB_TOKEN}`, "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "mpstruc-login-service", ...(payload ? { "Content-Type": "application/json" } : {}) },
      body: payload ? JSON.stringify(payload) : undefined,
    });
  } catch { throw new HttpError(502, "The login service could not reach GitHub. Try again in a minute."); }
  const text = await r.text();
  // A refused bot token is the site owner's problem, not the curator's: never answer 401 for it.
  if (r.status === 401) throw new HttpError(502, "GitHub refused the login service's token (it may have expired). Ask the site owner to renew GITHUB_TOKEN.");
  if (r.status === 403 && !/rate limit/i.test(text)) {
    throw new HttpError(502, "GitHub refused the login service's token for this repository. Ask the site owner to check that GITHUB_TOKEN has Contents: Read and write on team-website.");
  }
  return new Response(text, { status: r.status, headers: { "Content-Type": "application/json; charset=utf-8", "Cache-Control": "no-store" } });
}

async function route(request, env, url) {
  const p = url.pathname.replace(/\/+$/, "") || "/";
  const m = request.method;
  if (p.startsWith("/gh/")) return proxy(request, env, url);
  if (p === "/status" && m === "GET") return { ok: true, needsSetup: (await activeAdmins(env)).length === 0 };
  if (p === "/login" && m === "POST") return login(request, env);
  if (p === "/logout" && m === "POST") return logout(request, env);
  if (p === "/me" && m === "GET") { const { user, payload } = await session(request, env); return { user: publicUser(user), expires: new Date(payload.exp * 1000).toISOString() }; }
  if (p === "/setup" && m === "POST") return setup(request, env);
  if (p === "/password" && m === "POST") return changeOwnPassword(request, env);
  if (p === "/users" && m === "GET") { await requireAdmin(request, env); return { users: (await allUsers(env)).map(publicUser) }; }
  if (p === "/users" && m === "POST") {
    const { user: me } = await requireAdmin(request, env);
    const b = await body(request);
    return { ok: true, user: publicUser(await newUser(env, b, me.username)) };
  }
  const um = /^\/users\/([^/]+)(?:\/(reset|disable|enable))?$/.exec(p);
  if (um) {
    const name = decodeURIComponent(um[1]).toLowerCase();
    if (um[2] && m === "POST") return manageUser(request, env, name, um[2]);
    if (!um[2] && m === "DELETE") return manageUser(request, env, name, "delete");
  }
  throw new HttpError(404, "Not found.");
}

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const origin = request.headers.get("Origin") || "";
    // A plain visit (no Origin), e.g. the owner opening the Worker's address to check it runs.
    if (!origin && request.method === "GET" && url.pathname === "/") {
      return new Response("The mpstruc login service is running.\n", { headers: { "Content-Type": "text/plain; charset=utf-8" } });
    }
    if (!allowedOrigin(origin, env)) return json({ error: "This service only answers the mpstruc site." }, 403);
    if (request.method === "OPTIONS") return new Response(null, { status: 204, headers: corsHeaders(origin) });
    try {
      const out = await route(request, env, url);
      return out instanceof Response ? withCors(out, origin) : json(out, 200, origin);
    } catch (e) {
      if (e instanceof HttpError) return json({ error: e.message, ...e.extra }, e.status, origin);
      console.log("login service error:", e && e.name); // never the request, which may hold a password
      return json({ error: "The login service had a problem. Try again in a minute." }, 500, origin);
    }
  },
};
function withCors(res, origin) {
  const h = new Headers(res.headers);
  for (const [k, v] of Object.entries(corsHeaders(origin))) h.set(k, v);
  return new Response(res.body, { status: res.status, headers: h });
}
