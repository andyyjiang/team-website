// mpstruc login service: username and password sign-in for curators.
//
// A Cloudflare Worker (one file, no build step). It keeps the curator
// accounts in a D1 database (binding DB), checks passwords, hands out signed
// session tokens, and forwards a short allow-list of GitHub API calls with ONE
// bot token, so curators do not need GitHub accounts.
//
// D1 rather than Workers KV: D1 is strongly consistent and updates rows
// atomically, so sign-in attempt counters cannot be raced, a disabled account
// or a signed-out session ends at once everywhere, and the free plan allows
// 100,000 row writes a day (KV: 1,000).
//
// Secrets (Worker settings > Variables and Secrets, type "Secret"):
//   GITHUB_TOKEN    fine-grained GitHub token: Contents read/write on andyyjiang/team-website only
//   SESSION_SECRET  a random string of at least 32 characters that signs session tokens
//   SETUP_CODE      a one-off code that creates the first administrator
// Binding: D1 database as DB (the Worker creates its tables itself).
// Optional variables: ALLOWED_ORIGIN (default https://andyyjiang.github.io),
//   PBKDF2_ITERATIONS (default 20000; up to 100000 on the Workers Paid plan).
//
// Passwords and tokens are never logged or returned.

const REPO = "andyyjiang/team-website";
const BRANCH = "main";
const GITHUB_API = "https://api.github.com";
const DEFAULT_ORIGIN = "https://andyyjiang.github.io";
// PBKDF2 rounds for new password hashes. The Workers Free plan allows 10 ms of
// CPU per request; 100,000 rounds take about 15 ms, 20,000 about 3 ms (and
// changing a password hashes twice). Each account stores its own round count,
// and a sign-in re-hashes an account below the current setting, so raising
// PBKDF2_ITERATIONS later (Workers Paid) strengthens accounts as people sign in.
// Cloudflare's WebCrypto allows at most 100,000.
const DEFAULT_ITERATIONS = 20000, MAX_ITERATIONS = 100000;
const iterationsFor = (env) => {
  const n = Number.parseInt(env.PBKDF2_ITERATIONS || "", 10);
  return Number.isFinite(n) ? Math.min(MAX_ITERATIONS, Math.max(DEFAULT_ITERATIONS, n)) : DEFAULT_ITERATIONS;
};
const SHORT_SESSION = 12 * 3600;        // seconds
const LONG_SESSION = 30 * 24 * 3600;    // "Keep me signed in"
const MAX_BODY = 64 * 1024;             // bytes, for sign-in and account requests
// Saves send the whole data file base64-encoded (a third larger). data/curation.json
// is about 0.85 MB in November 2026 and grows by about 60 KB a week; 40 MB is far
// beyond that and below Cloudflare's 100 MB request limit.
const MAX_SAVE = 40 * 1024 * 1024;
const MIN_PASSWORD = 10, MAX_PASSWORD = 200;
const MIN_SECRET = 32;
// Sign-in attempts are counted BEFORE the password is checked, atomically (a D1
// upsert), and a successful sign-in gives its attempt back. Per 15 minutes:
//   per username and address  10   (one person mistyping, or one attacker on one account)
//   per address               30   (one attacker trying many usernames; an office may share an address)
//   per username, new devices 20   (an attacker spread over many addresses)
//   per username, known device 10  (a browser that has signed in to this account before;
//                                   it is not held up by other people's failures, so
//                                   strangers cannot lock a curator out)
// IPv6 addresses are counted per /64, the block one connection usually gets.
const WINDOW = 15 * 60;
const LIMITS = { userIp: 10, ip: 30, user: 20, device: 10 };
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
export async function hashPassword(password, saltB64, iterations = DEFAULT_ITERATIONS) {
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
const dummy = (env) => ({ salt: "AAAAAAAAAAAAAAAAAAAAAA", iterations: iterationsFor(env), hash: "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" });

function checkNewPassword(pw) {
  if (typeof pw !== "string" || pw.length < MIN_PASSWORD) throw new HttpError(400, `The password must be at least ${MIN_PASSWORD} characters long.`);
  if (pw.length > MAX_PASSWORD) throw new HttpError(400, `The password can be at most ${MAX_PASSWORD} characters long.`);
}

// ---- set-up checks ----
// What is missing, in words the site owner can act on. GET / lists all of it;
// the site's requests stop at the first problem that affects them.
function setupProblems(env) {
  const out = [];
  if (!env.DB || typeof env.DB.prepare !== "function") {
    out.push("The database is not connected. In the Worker, open Settings > Bindings, add a D1 database binding with the variable name DB (capitals) and the database mpstruc-users, then deploy.");
  }
  if (!env.SESSION_SECRET) out.push("The secret SESSION_SECRET is missing. Add it under Settings > Variables and Secrets (type Secret), then deploy.");
  else if (String(env.SESSION_SECRET).length < MIN_SECRET) out.push(`The secret SESSION_SECRET is too short. Use at least ${MIN_SECRET} characters.`);
  if (!env.GITHUB_TOKEN) out.push("The secret GITHUB_TOKEN is missing. Add it under Settings > Variables and Secrets (type Secret), then deploy.");
  return out;
}
function requireSetUp(env) {
  const p = setupProblems(env).filter((x) => !/GITHUB_TOKEN/.test(x));
  if (p.length) throw new HttpError(503, `The sign-in service is not set up completely. Ask the site owner to open the service's address, which lists what is missing.`, { setupProblem: true });
}

// ---- database ----
let schemaReady = false; // per Worker instance; the statements are safe to repeat
async function db(env) {
  requireSetUp(env);
  if (!schemaReady) {
    await env.DB.batch([
      env.DB.prepare("CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, data TEXT NOT NULL)"),
      env.DB.prepare("CREATE TABLE IF NOT EXISTS revoked (sig TEXT PRIMARY KEY, until INTEGER NOT NULL)"),
      env.DB.prepare("CREATE TABLE IF NOT EXISTS attempts (key TEXT PRIMARY KEY, n INTEGER NOT NULL, reset INTEGER NOT NULL)"),
    ]);
    schemaReady = true;
  }
  return env.DB;
}

// ---- sessions ----
async function hmacKey(env) {
  requireSetUp(env);
  return crypto.subtle.importKey("raw", enc.encode(env.SESSION_SECRET), { name: "HMAC", hash: "SHA-256" }, false, ["sign"]);
}
async function sign(env, payload) {
  const body = b64url(enc.encode(JSON.stringify(payload)));
  const sig = new Uint8Array(await crypto.subtle.sign("HMAC", await hmacKey(env), enc.encode(body)));
  return `${body}.${b64url(sig)}`;
}
// A signed token's payload (session or device token), or null.
async function readSigned(env, token) {
  const m = /^([A-Za-z0-9_-]{10,2000})\.([A-Za-z0-9_-]{43})$/.exec(typeof token === "string" ? token : "");
  if (!m) return null;
  const want = new Uint8Array(await crypto.subtle.sign("HMAC", await hmacKey(env), enc.encode(m[1])));
  if (!equalBytes(want, fromB64url(m[2]))) return null;
  try {
    const p = JSON.parse(new TextDecoder().decode(fromB64url(m[1])));
    return p && typeof p === "object" ? { ...p, sig: m[2] } : null;
  } catch { return null; }
}
async function readToken(env, token) {
  const p = await readSigned(env, token);
  return p && typeof p.u === "string" && typeof p.a === "string" && typeof p.exp === "number" ? p : null;
}
const now = () => Math.floor(Date.now() / 1000);
const nonce = () => b64url(crypto.getRandomValues(new Uint8Array(9))); // makes every token distinct, so one can be signed out alone
const bearer = (request) => { const a = request.headers.get("Authorization") || ""; return a.startsWith("Bearer ") ? a.slice(7).trim() : ""; };

const SIGNED_OUT = "Your sign-in has ended. Please sign in again.";
// The signed-in user for this request, checked against the account as stored now:
// a deleted, disabled or re-created account, a changed password or a reset ends
// every older session.
async function session(request, env) {
  const p = await readToken(env, bearer(request));
  if (!p) throw new HttpError(401, SIGNED_OUT);
  if (p.exp <= now()) throw new HttpError(401, "Your sign-in has expired. Please sign in again.");
  const d = await db(env);
  const [u, r] = await d.batch([
    d.prepare("SELECT data FROM users WHERE username = ?1").bind(p.u),
    d.prepare("SELECT 1 AS hit FROM revoked WHERE sig = ?1").bind(p.sig),
  ]);
  if (r.results.length) throw new HttpError(401, SIGNED_OUT);
  const user = u.results.length ? JSON.parse(u.results[0].data) : null;
  // The account id is random per account, so a token for a deleted account never
  // works for a new account that is given the same username.
  if (!user || user.id !== p.a) throw new HttpError(401, "This account no longer exists. Ask an administrator.");
  if (user.disabled) throw new HttpError(401, "This account has been disabled. Ask an administrator to enable it again.");
  if ((user.sessionVersion || 0) !== p.v) throw new HttpError(401, "Your password was changed or reset, so you have been signed out. Please sign in again.");
  return { user, payload: p };
}
async function requireAdmin(request, env) {
  const s = await session(request, env);
  if (s.user.role !== "admin") throw new HttpError(403, "Only administrators can manage curators.");
  return s;
}
const sessionToken = (env, user, exp) => sign(env, { u: user.username, a: user.id, r: user.role, v: user.sessionVersion || 0, exp, n: nonce() });

// ---- users ----
async function getUser(env, name) {
  if (typeof name !== "string" || !USERNAME_RE.test(name)) return null;
  const row = await (await db(env)).prepare("SELECT data FROM users WHERE username = ?1").bind(name).first();
  return row ? JSON.parse(row.data) : null;
}
async function putUser(env, u) {
  await (await db(env)).prepare("INSERT INTO users (username, data) VALUES (?1, ?2) ON CONFLICT (username) DO UPDATE SET data = excluded.data")
    .bind(u.username, JSON.stringify(u)).run();
}
async function allUsers(env) {
  const { results } = await (await db(env)).prepare("SELECT data FROM users ORDER BY username").all();
  return results.map((r) => JSON.parse(r.data));
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
  const u = { username: name, id: nonce(), displayName: cleanDisplay(displayName, name), role, ...(await hashPassword(password, null, iterationsFor(env))),
    createdAt: new Date().toISOString(), createdBy, disabled: false, sessionVersion: 0 };
  // DO NOTHING on a clash: two admins adding the same name at once cannot overwrite each other.
  const res = await (await db(env)).prepare("INSERT INTO users (username, data) VALUES (?1, ?2) ON CONFLICT (username) DO NOTHING")
    .bind(u.username, JSON.stringify(u)).run();
  if (!res.meta || !res.meta.changes) throw new HttpError(409, `There is already an account called ${name}.`);
  return u;
}
async function activeAdmins(env) { return (await allUsers(env)).filter((u) => u.role === "admin" && !u.disabled); }

// ---- attempt limits (sign-in and setup) ----
// IPv4 as is; IPv6 by its first four groups (the /64 one connection usually has).
export function addressKey(ip) {
  const s = String(ip || "unknown").toLowerCase();
  if (!s.includes(":")) return s;
  const [head, tail] = s.split("::");
  const h = head ? head.split(":") : [];
  const t = tail === undefined ? [] : tail ? tail.split(":") : [];
  const groups = tail === undefined ? h : [...h, ...Array(Math.max(0, 8 - h.length - t.length)).fill("0"), ...t];
  return `${groups.slice(0, 4).map((g) => g.replace(/^0+(?=.)/, "")).join(":")}::/64`;
}
// Counts one attempt against each [key, limit], atomically. Returns the time
// the earliest blocking window ends, or 0 when every count is within its limit.
// A count already at its limit is not written again, so a flood of refused
// attempts costs no database writes.
const HIT = `INSERT INTO attempts (key, n, reset) VALUES (?1, 1, ?2 + ?3)
  ON CONFLICT (key) DO UPDATE SET
    n = CASE WHEN reset <= ?2 THEN 1 ELSE n + 1 END,
    reset = CASE WHEN reset <= ?2 THEN ?2 + ?3 ELSE reset END
  WHERE reset <= ?2 OR n < ?4
  RETURNING n`;
async function countAttempt(env, pairs) {
  const d = await db(env), t = now();
  const res = await d.batch([
    ...pairs.map(([key, limit]) => d.prepare(HIT).bind(key, t, WINDOW, limit)),
    d.prepare(`SELECT key, reset FROM attempts WHERE key IN (${pairs.map((_, i) => `?${i + 1}`).join(", ")})`).bind(...pairs.map(([k]) => k)),
  ]);
  const resets = new Map(res[pairs.length].results.map((r) => [r.key, r.reset]));
  const blocked = pairs.filter((_, i) => !res[i].results.length).map(([k]) => resets.get(k) || t + WINDOW);
  // Now and then, clear out windows and revocations that have ended.
  if (Math.random() < 0.05) {
    await d.batch([d.prepare("DELETE FROM attempts WHERE reset <= ?1").bind(t), d.prepare("DELETE FROM revoked WHERE until <= ?1").bind(t)]);
  }
  return blocked.length ? Math.max(...blocked) : 0;
}
// A successful sign-in gives its attempt back.
async function giveBack(env, clearKey, keys) {
  const d = await db(env);
  await d.batch([
    d.prepare("DELETE FROM attempts WHERE key = ?1").bind(clearKey),
    ...keys.map((k) => d.prepare("UPDATE attempts SET n = n - 1 WHERE key = ?1 AND n > 0").bind(k)),
  ]);
}
function tooMany(reset) {
  const mins = Math.max(1, Math.ceil((reset - now()) / 60));
  return new HttpError(429, `Too many failed sign-in attempts. Try again in ${mins} ${mins === 1 ? "minute" : "minutes"}.`, { retryMinutes: mins });
}

// ---- request body ----
async function body(request, max = MAX_BODY) {
  const tooLarge = () => new HttpError(413, max === MAX_SAVE ? `This change is too large to save (over ${max / 1024 / 1024} MB).` : "This request is too large.");
  const len = Number(request.headers.get("Content-Length") || 0);
  if (len > max) throw tooLarge();
  const text = await request.text();
  if (text.length > max) throw tooLarge();
  try {
    const v = JSON.parse(text || "{}");
    if (!v || typeof v !== "object" || Array.isArray(v)) throw new Error();
    return v;
  } catch { throw new HttpError(400, "The request was not valid JSON."); }
}

// ---- handlers ----
async function login(request, env) {
  const b = await body(request);
  const name = typeof b.username === "string" ? b.username.trim().toLowerCase().slice(0, 64) : "";
  const password = typeof b.password === "string" ? b.password : "";
  if (!name || !password) throw new HttpError(400, "Enter your username and password.");
  const ip = addressKey(request.headers.get("CF-Connecting-IP"));
  // A device token: this browser signed in to this account before.
  const dev = await readSigned(env, b.device);
  const known = dev && dev.d === name;
  const userIp = `ui:${name}|${ip}`, ipKey = `ip:${ip}`, broad = known ? `dev:${dev.sig}` : `u:${name}`;
  const blocked = await countAttempt(env, [[ipKey, LIMITS.ip], [userIp, LIMITS.userIp], [broad, known ? LIMITS.device : LIMITS.user]]);
  if (blocked) throw tooMany(blocked);
  const user = await getUser(env, name);
  const ok = await checkPassword(user || dummy(env), password.slice(0, MAX_PASSWORD));
  if (!user || !ok) throw new HttpError(401, "Wrong username or password.");
  if (user.disabled) throw new HttpError(403, "This account has been disabled. Ask an administrator to enable it again.");
  await giveBack(env, userIp, [ipKey, broad]);
  if ((user.iterations || 0) < iterationsFor(env)) { // stronger hashing was switched on: upgrade this account
    Object.assign(user, await hashPassword(password.slice(0, MAX_PASSWORD), null, iterationsFor(env)));
    await putUser(env, user);
  }
  const exp = now() + (b.keep === true ? LONG_SESSION : SHORT_SESSION);
  const token = await sessionToken(env, user, exp);
  const device = known ? b.device : await sign(env, { d: user.username, n: nonce() });
  return { token, device, expires: new Date(exp * 1000).toISOString(), user: publicUser(user) };
}

async function logout(request, env) {
  let s = null;
  try { s = await session(request, env); } catch { /* already signed out: nothing to do */ }
  if (s && s.payload.exp > now()) {
    await (await db(env)).prepare("INSERT INTO revoked (sig, until) VALUES (?1, ?2) ON CONFLICT (sig) DO NOTHING").bind(s.payload.sig, s.payload.exp).run();
  }
  return { ok: true };
}

async function setup(request, env) {
  if ((await activeAdmins(env)).length) throw new HttpError(403, "Setup is already done. Sign in with your username and password.");
  if (!env.SETUP_CODE) throw new HttpError(503, "The sign-in service is not set up completely (the secret SETUP_CODE is missing). Ask the site owner.");
  const b = await body(request);
  const ip = addressKey(request.headers.get("CF-Connecting-IP"));
  const blocked = await countAttempt(env, [[`setup:${ip}`, LIMITS.userIp]]);
  if (blocked) throw tooMany(blocked);
  const code = typeof b.setupCode === "string" ? b.setupCode.trim() : "";
  const digest = async (s) => new Uint8Array(await crypto.subtle.digest("SHA-256", enc.encode(s)));
  if (!code || !equalBytes(await digest(code), await digest(String(env.SETUP_CODE).trim()))) throw new HttpError(401, "That setup code is not right.");
  const u = await newUser(env, { ...b, role: "admin" }, "setup");
  return { ok: true, user: publicUser(u) };
}

async function changeOwnPassword(request, env) {
  const { user, payload } = await session(request, env);
  const b = await body(request);
  checkNewPassword(b.newPassword);
  if (!(await checkPassword(user, String(b.oldPassword || "").slice(0, MAX_PASSWORD)))) throw new HttpError(400, "Your current password is not right.");
  Object.assign(user, await hashPassword(b.newPassword, null, iterationsFor(env)), { sessionVersion: (user.sessionVersion || 0) + 1 });
  await putUser(env, user);
  // A fresh session for this browser; every other one has ended.
  return { ok: true, token: await sessionToken(env, user, payload.exp), user: publicUser(user) };
}

async function manageUser(request, env, name, action) {
  const { user: me } = await requireAdmin(request, env);
  const u = await getUser(env, name);
  if (!u) throw new HttpError(404, `There is no account called ${name}.`);
  const lastAdmin = async () => u.role === "admin" && !u.disabled && (await activeAdmins(env)).length <= 1;
  if (action === "reset") {
    const b = await body(request);
    checkNewPassword(b.password);
    Object.assign(u, await hashPassword(b.password, null, iterationsFor(env)), { sessionVersion: (u.sessionVersion || 0) + 1 });
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
    await (await db(env)).prepare("DELETE FROM users WHERE username = ?1").bind(u.username).run();
    return { ok: true };
  }
  await putUser(env, u);
  return { ok: true, user: publicUser(u) };
}

// GitHub proxy: /gh/<GitHub API path>, allow-list only.
const B64_RE = /^[A-Za-z0-9+/]*={0,2}$/, SHA_RE = /^[0-9a-f]{40}$|^[0-9a-f]{64}$/;
async function proxy(request, env, url) {
  const path = url.pathname.slice(3);
  const kind = allowedGitHubCall(request.method, path, url.search);
  if (!kind) throw new HttpError(403, "The login service does not allow this request.");
  const { user } = await session(request, env);
  if (!ROLES.includes(user.role)) throw new HttpError(403, "This account cannot save changes.");
  if (!env.GITHUB_TOKEN) throw new HttpError(503, "The sign-in service is not set up completely (the secret GITHUB_TOKEN is missing). Ask the site owner.");
  let payload;
  if (kind === "write") {
    const b = await body(request, MAX_SAVE);
    if (typeof b.content !== "string" || !B64_RE.test(b.content) || typeof b.message !== "string"
        || (b.sha !== undefined && !(typeof b.sha === "string" && SHA_RE.test(b.sha)))
        || (b.branch !== undefined && b.branch !== BRANCH)) throw new HttpError(400, "This change is not in the expected form.");
    // The commit is made by the bot; the curator is named in the message.
    const rest = b.message.replace(/^Curator [^:\n]{1,80}:\s*/, "").replace(/[\r\n]+/g, " ").slice(0, 900);
    // Built by hand so the large content is not re-encoded (that costs CPU time,
    // which the free plan limits). It is safe to splice in: it is plain base64.
    payload = `{"message":${JSON.stringify(`Curator ${user.username}: ${rest}`)},"branch":"${BRANCH}",${b.sha ? `"sha":"${b.sha}",` : ""}"content":"${b.content}"}`;
  }
  let r;
  try {
    r = await fetch(GITHUB_API + path + url.search, {
      method: request.method,
      headers: { Accept: "application/vnd.github+json", Authorization: `Bearer ${env.GITHUB_TOKEN}`, "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "mpstruc-login-service", ...(payload ? { "Content-Type": "application/json" } : {}) },
      body: payload,
    });
  } catch { throw new HttpError(502, "The login service could not reach GitHub. Try again in a minute."); }
  // A refused bot token is the site owner's problem, not the curator's: never answer 401 for it.
  if (r.status === 401) throw new HttpError(502, "GitHub refused the login service's token (it may have expired). Ask the site owner to renew GITHUB_TOKEN.");
  if (r.status === 403) {
    const text = await r.text();
    if (!/rate limit/i.test(text)) {
      throw new HttpError(502, "GitHub refused the login service's token for this repository. Ask the site owner to check that GITHUB_TOKEN has Contents: Read and write on team-website.");
    }
    return new Response(text, { status: 403, headers: { "Content-Type": "application/json; charset=utf-8", "Cache-Control": "no-store" } });
  }
  // Passed through as a stream: the data files are large, and decoding them would cost CPU time.
  return new Response(r.body, { status: r.status, headers: { "Content-Type": "application/json; charset=utf-8", "Cache-Control": "no-store" } });
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

// The page the site owner sees when opening the Worker's address: whether every
// part of the set-up is in place, and what to fix if not. The GitHub check is
// remembered for a minute, so reloading this page does not use up GitHub's limits.
let githubCheck = { at: 0, token: "", text: "" };
async function checkGitHub(env) {
  if (githubCheck.token === env.GITHUB_TOKEN && Date.now() - githubCheck.at < 60000) return githubCheck.text;
  let text;
  try {
    const r = await fetch(`${GITHUB_API}/repos/${REPO}`, { headers: { Accept: "application/vnd.github+json",
      Authorization: `Bearer ${env.GITHUB_TOKEN}`, "X-GitHub-Api-Version": "2022-11-28", "User-Agent": "mpstruc-login-service" } });
    await r.body?.cancel();
    text = r.status === 200 ? ""
      : r.status === 401 ? "GitHub does not accept GITHUB_TOKEN: it is mistyped or has expired. Make a new token (setup step 1) and paste it into GITHUB_TOKEN."
      : r.status === 404 || r.status === 403 ? `GITHUB_TOKEN cannot see ${REPO}. When making the token, choose 'Only select repositories' and pick team-website.`
      : `GitHub answered ${r.status} when checking GITHUB_TOKEN. Try again in a few minutes.`;
  } catch { text = "Could not reach GitHub to check GITHUB_TOKEN. Try again in a few minutes."; }
  githubCheck = { at: Date.now(), token: env.GITHUB_TOKEN, text };
  return text;
}
async function healthPage(env) {
  const problems = setupProblems(env);
  let needsSetup = null;
  if (!problems.some((x) => /database|SESSION_SECRET/.test(x))) {
    try { needsSetup = (await activeAdmins(env)).length === 0; } catch (e) {
      problems.push(`The database could not be used (${e && e.message ? String(e.message).slice(0, 200) : "unknown error"}). Check the DB binding points to the D1 database mpstruc-users.`);
    }
  }
  if (env.GITHUB_TOKEN) { const g = await checkGitHub(env); if (g) problems.push(g); }
  if (needsSetup && !env.SETUP_CODE) problems.push("The secret SETUP_CODE is missing. It is needed to create the first administrator.");
  const text = problems.length
    ? `NOT READY: the mpstruc login service needs attention.\n\n${problems.map((x) => `- ${x}`).join("\n")}\n\nAfter fixing, reload this page.\n`
    : `READY: the mpstruc login service is running.\n- Database: connected.\n- GITHUB_TOKEN: GitHub accepts it and it can see ${REPO}.\n- ${needsSetup ? "No administrator yet: create one on the Curators tab with SETUP_CODE." : "Administrators: set up."}\n`;
  return new Response(text, { status: problems.length ? 503 : 200, headers: { "Content-Type": "text/plain; charset=utf-8", "Cache-Control": "no-store" } });
}

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const origin = request.headers.get("Origin") || "";
    // A plain visit (no Origin), e.g. the owner opening the Worker's address to check it.
    if (!origin && request.method === "GET" && url.pathname === "/") return healthPage(env);
    if (!allowedOrigin(origin, env)) return json({ error: "This service only answers the mpstruc site." }, 403);
    if (request.method === "OPTIONS") return new Response(null, { status: 204, headers: corsHeaders(origin) });
    try {
      const out = await route(request, env, url);
      return out instanceof Response ? withCors(out, origin) : json(out, 200, origin);
    } catch (e) {
      if (e instanceof HttpError) return json({ error: e.message, ...e.extra }, e.status, origin);
      console.log("login service error:", e && e.name, e && e.message ? String(e.message).slice(0, 200) : ""); // never the request, which may hold a password
      return json({ error: "The login service had a problem. Try again in a minute. If it keeps happening, ask the site owner to open the service's address and check its logs in Cloudflare." }, 500, origin);
    }
  },
};
function withCors(res, origin) {
  const h = new Headers(res.headers);
  for (const [k, v] of Object.entries(corsHeaders(origin))) h.set(k, v);
  return new Response(res.body, { status: res.status, headers: h });
}
