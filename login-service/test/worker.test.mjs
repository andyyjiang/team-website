import { test, before, after } from "node:test";
import assert from "node:assert/strict";
import { createHmac } from "node:crypto";
import { startWorker, SECRETS } from "./harness.mjs";

const REPO = "/repos/andyyjiang/team-website";
const CUR = `${REPO}/contents/data/curation.json`, FB = `${REPO}/contents/data/feedback.json`;
let w;
const ADMIN = { username: "Owner", displayName: "Site Owner", password: "admin-password-1" };
let adminToken;

before(async () => {
  w = await startWorker({ seed: { "data/curation.json": '{"version":1,"entries":{},"log":[]}\n', "data/feedback.json": '{"version":1,"items":[]}\n' } });
});
after(async () => { await w.mf.dispose(); });

const b64url = (b) => Buffer.from(b).toString("base64url");
function forge(payload, secret = SECRETS.SESSION_SECRET) {
  const body = b64url(JSON.stringify(payload));
  return `${body}.${createHmac("sha256", secret).update(body).digest("base64url")}`;
}
let ipN = 10;
const freshIp = () => `198.51.100.${ipN++}`;

test("status reports that setup is needed, and /setup checks the code", async () => {
  let r = await w.call("/status");
  assert.equal(r.status, 200);
  assert.equal(r.data.needsSetup, true);
  r = await w.call("/setup", { method: "POST", body: { ...ADMIN, setupCode: "wrong" }, ip: freshIp() });
  assert.equal(r.status, 401);
  assert.match(r.data.error, /setup code/);
  r = await w.call("/setup", { method: "POST", body: { ...ADMIN, password: "short", setupCode: SECRETS.SETUP_CODE }, ip: freshIp() });
  assert.equal(r.status, 400);
  assert.match(r.data.error, /at least 10/);
  r = await w.call("/setup", { method: "POST", body: { ...ADMIN, setupCode: SECRETS.SETUP_CODE }, ip: freshIp() });
  assert.equal(r.status, 200, JSON.stringify(r.data));
  assert.deepEqual([r.data.user.username, r.data.user.role, r.data.user.displayName], ["owner", "admin", "Site Owner"]);
});

test("setup is locked once an administrator exists", async () => {
  const r = await w.call("/setup", { method: "POST", body: { username: "intruder", password: "intruder-password", setupCode: SECRETS.SETUP_CODE } });
  assert.equal(r.status, 403);
  assert.equal((await w.call("/status")).data.needsSetup, false);
});

test("passwords are stored as salted PBKDF2-SHA256, never in clear", async () => {
  const raw = JSON.parse(await (await w.kv()).get("user:owner"));
  assert.equal(raw.iterations, 100000);
  assert.equal(raw.salt.length, 22); // 16 random bytes
  assert.equal(Buffer.from(raw.hash, "base64url").length, 32);
  assert.ok(!JSON.stringify(raw).includes(ADMIN.password));
  const { hashPassword } = await import("../worker.js");
  const again = await hashPassword(ADMIN.password, raw.salt, raw.iterations);
  assert.equal(again.hash, raw.hash);
  const other = await hashPassword(ADMIN.password);
  assert.notEqual(other.salt, raw.salt);
  assert.notEqual(other.hash, raw.hash);
});

test("login succeeds with the right password (any case of username) and fails otherwise", async () => {
  let r = await w.call("/login", { method: "POST", body: { username: "OWNER ", password: ADMIN.password } });
  assert.equal(r.status, 200);
  assert.equal(r.data.user.username, "owner");
  assert.equal(r.data.user.role, "admin");
  assert.ok(!("hash" in r.data.user) && !("salt" in r.data.user));
  adminToken = r.data.token;
  const exp = JSON.parse(Buffer.from(adminToken.split(".")[0], "base64url")).exp;
  assert.ok(Math.abs(exp - Date.now() / 1000 - 12 * 3600) < 60, "12 hours by default");
  r = await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password, keep: true } });
  const exp2 = JSON.parse(Buffer.from(r.data.token.split(".")[0], "base64url")).exp;
  assert.ok(Math.abs(exp2 - Date.now() / 1000 - 30 * 86400) < 60, "30 days when kept");
  r = await w.call("/login", { method: "POST", body: { username: "owner", password: "not-the-password" }, ip: freshIp() });
  assert.equal(r.status, 401);
  assert.equal(r.data.error, "Wrong username or password.");
  r = await w.call("/login", { method: "POST", body: { username: "nobody", password: "not-the-password" }, ip: freshIp() });
  assert.equal(r.status, 401);
  assert.equal(r.data.error, "Wrong username or password.");
  r = await w.call("/me", { token: adminToken });
  assert.equal(r.status, 200);
  assert.equal(r.data.user.username, "owner");
});

test("sessions: tampered, expired and logged-out tokens are refused", async () => {
  const [body] = adminToken.split(".");
  let r = await w.call("/me", { token: `${body}.${"A".repeat(43)}` });
  assert.equal(r.status, 401);
  r = await w.call("/me", { token: forge({ u: "owner", r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) + 600 }, "some-other-secret-value") });
  assert.equal(r.status, 401);
  r = await w.call("/me", { token: forge({ u: "owner", r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) - 5 }) });
  assert.equal(r.status, 401);
  assert.match(r.data.error, /expired/);
  r = await w.call("/me", { token: forge({ u: "owner", r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) + 600 }) });
  assert.equal(r.status, 200, "a correctly signed, unexpired token works");
  const t = (await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password } })).data.token;
  assert.equal((await w.call("/logout", { method: "POST", token: t })).status, 200);
  assert.equal((await w.call("/me", { token: t })).status, 401);
  assert.equal((await w.call("/me", { token: adminToken })).status, 200, "other sessions continue");
  assert.equal((await w.call("/logout", { method: "POST" })).status, 200, "logout without a session is fine");
});

let curatorToken;
test("admins manage curators; curators cannot", async () => {
  let r = await w.call("/users", { method: "POST", token: adminToken, body: { username: "Ana.Lee", displayName: "Ana Lee", password: "temporary-pass-1", role: "curator" } });
  assert.equal(r.status, 200, JSON.stringify(r.data));
  assert.equal(r.data.user.username, "ana.lee");
  assert.equal(r.data.user.createdBy, "owner");
  r = await w.call("/users", { method: "POST", token: adminToken, body: { username: "ana.lee", displayName: "x", password: "temporary-pass-1", role: "curator" } });
  assert.equal(r.status, 409);
  r = await w.call("/users", { method: "POST", token: adminToken, body: { username: "a b", password: "temporary-pass-1", role: "curator" } });
  assert.equal(r.status, 400);
  r = await w.call("/users", { method: "POST", token: adminToken, body: { username: "bob", password: "short", role: "curator" } });
  assert.equal(r.status, 400);
  r = await w.call("/users", { method: "POST", token: adminToken, body: { username: "bob", password: "long-enough-1", role: "boss" } });
  assert.equal(r.status, 400);
  r = await w.call("/users", { token: adminToken });
  assert.deepEqual(r.data.users.map((u) => [u.username, u.role, u.disabled]), [["ana.lee", "curator", false], ["owner", "admin", false]]);
  assert.ok(r.data.users.every((u) => !u.hash && !u.salt));

  curatorToken = (await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "temporary-pass-1" } })).data.token;
  assert.ok(curatorToken);
  for (const [path, method, body] of [["/users", "GET"], ["/users", "POST", { username: "eve", password: "eve-password-1", role: "admin" }],
    ["/users/owner/disable", "POST"], ["/users/owner/reset", "POST", { password: "hijacked-pass" }], ["/users/owner", "DELETE"]]) {
    r = await w.call(path, { method, body, token: curatorToken });
    assert.equal(r.status, 403, `${method} ${path}`);
  }
  r = await w.call("/users/owner/disable", { method: "POST", token: adminToken });
  assert.equal(r.status, 400, "cannot disable yourself");
  r = await w.call("/users/owner", { method: "DELETE", token: adminToken });
  assert.equal(r.status, 400, "cannot delete yourself");
});

test("changing your own password needs the old one and ends other sessions", async () => {
  let r = await w.call("/password", { method: "POST", token: curatorToken, body: { oldPassword: "wrong-old-pass", newPassword: "brand-new-pass-1" } });
  assert.equal(r.status, 400);
  r = await w.call("/password", { method: "POST", token: curatorToken, body: { oldPassword: "temporary-pass-1", newPassword: "tiny" } });
  assert.equal(r.status, 400);
  const other = (await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "temporary-pass-1" } })).data.token;
  r = await w.call("/password", { method: "POST", token: curatorToken, body: { oldPassword: "temporary-pass-1", newPassword: "brand-new-pass-1" } });
  assert.equal(r.status, 200);
  assert.equal((await w.call("/me", { token: other })).status, 401, "other session revoked");
  assert.equal((await w.call("/me", { token: curatorToken })).status, 401, "the old token too");
  curatorToken = r.data.token;
  assert.equal((await w.call("/me", { token: curatorToken })).status, 200, "the new token works");
  assert.equal((await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "temporary-pass-1" }, ip: freshIp() })).status, 401);
  assert.equal((await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "brand-new-pass-1" } })).status, 200);
});

test("reset, disable, enable and delete by an admin", async () => {
  let r = await w.call("/users/ana.lee/reset", { method: "POST", token: adminToken, body: { password: "reset-by-admin-1" } });
  assert.equal(r.status, 200);
  assert.equal((await w.call("/me", { token: curatorToken })).status, 401, "reset revokes sessions");
  curatorToken = (await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "reset-by-admin-1" } })).data.token;
  r = await w.call("/users/ana.lee/disable", { method: "POST", token: adminToken });
  assert.equal(r.status, 200);
  assert.equal(r.data.user.disabled, true);
  r = await w.call("/me", { token: curatorToken });
  assert.equal(r.status, 401);
  assert.match(r.data.error, /disabled/);
  r = await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "reset-by-admin-1" } });
  assert.equal(r.status, 403);
  assert.match(r.data.error, /disabled/);
  r = await w.call("/users/ana.lee/enable", { method: "POST", token: adminToken });
  assert.equal(r.data.user.disabled, false);
  curatorToken = (await w.call("/login", { method: "POST", body: { username: "ana.lee", password: "reset-by-admin-1" } })).data.token;
  assert.ok(curatorToken);
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "temp", password: "temp-password-1", role: "curator" } });
  r = await w.call("/users/temp", { method: "DELETE", token: adminToken });
  assert.equal(r.status, 200);
  assert.equal(await (await w.kv()).get("user:temp"), null);
  r = await w.call("/users/temp", { method: "DELETE", token: adminToken });
  assert.equal(r.status, 404);
});

test("rate limit: 10 failures per username, then 429 even with the right password", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "victim", password: "victim-password-1", role: "curator" } });
  for (let i = 0; i < 10; i++) {
    const r = await w.call("/login", { method: "POST", body: { username: "victim", password: `guess-${i}-xxxxx` }, ip: freshIp() });
    assert.equal(r.status, 401);
  }
  const r = await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip: freshIp() });
  assert.equal(r.status, 429);
  assert.match(r.data.error, /^Too many failed sign-in attempts\. Try again in 15 minutes\.$/);
  assert.equal(r.data.retryMinutes, 15);
  // When the window has passed, the count starts again.
  await (await w.kv()).put("fail:user:victim", JSON.stringify({ n: 10, reset: Math.floor(Date.now() / 1000) - 1 }));
  assert.equal((await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip: freshIp() })).status, 200);
});

test("rate limit: 10 failures per IP address across usernames", async () => {
  const ip = "192.0.2.77";
  for (let i = 0; i < 10; i++) await w.call("/login", { method: "POST", body: { username: `user${i}`, password: "whatever-pass" }, ip });
  const r = await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password }, ip });
  assert.equal(r.status, 429);
  assert.equal((await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password }, ip: freshIp() })).status, 200);
});

test("CORS: only the site and localhost", async () => {
  let r = await w.call("/status", { origin: "https://evil.example" });
  assert.equal(r.status, 403);
  assert.equal(r.headers.get("Access-Control-Allow-Origin"), null);
  r = await w.call("/status", { origin: "https://andyyjiang.github.io.evil.example" });
  assert.equal(r.status, 403);
  r = await w.call("/status", { origin: null });
  assert.equal(r.status, 403);
  r = await w.call("/status", { origin: "https://andyyjiang.github.io" });
  assert.equal(r.headers.get("Access-Control-Allow-Origin"), "https://andyyjiang.github.io");
  r = await w.call("/status", { origin: "http://localhost:8000" });
  assert.equal(r.headers.get("Access-Control-Allow-Origin"), "http://localhost:8000");
  r = await w.call("/login", { method: "OPTIONS", origin: "http://127.0.0.1:5500" });
  assert.equal(r.status, 204);
  assert.match(r.headers.get("Access-Control-Allow-Headers"), /Authorization/);
  r = await w.call("/login", { method: "OPTIONS", origin: "http://localhost.evil.example" });
  assert.equal(r.status, 403);
  r = await w.call("/", { origin: null });
  assert.equal(r.status, 200, "a plain visit says the service runs");
});

test("GitHub proxy: reads go through with the bot token", async () => {
  const before = w.github.calls.length;
  let r = await w.call(`/gh${CUR}?ref=main`, { token: curatorToken });
  assert.equal(r.status, 200);
  assert.equal(Buffer.from(r.data.content, "base64").toString(), '{"version":1,"entries":{},"log":[]}\n');
  const c = w.github.calls[before];
  assert.equal(c.host, "api.github.com");
  assert.equal(c.authorization, "Bearer bot-token");
  assert.equal(c.userAgent, "mpstruc-login-service");
  r = await w.call(`/gh${REPO}/git/blobs/${r.data.sha}`, { token: curatorToken });
  assert.equal(r.status, 200);
  r = await w.call(`/gh${FB}?ref=main`, { token: curatorToken });
  assert.equal(r.status, 200);
  r = await w.call(`/gh${CUR}?ref=main`);
  assert.equal(r.status, 401, "a session is needed");
});

test("GitHub proxy: everything outside the allow-list is refused before GitHub is called", async () => {
  const before = w.github.calls.length;
  const bad = [
    ["GET", "/user"], ["GET", `${REPO}`], ["GET", `${REPO}/contents/js/main.js?ref=main`], ["PUT", `${REPO}/contents/js/main.js`],
    ["GET", `${CUR}`], ["GET", `${CUR}?ref=other`], ["GET", `${CUR}?ref=main&x=1`], ["PUT", `${CUR}?ref=main`], ["DELETE", CUR],
    ["POST", CUR], ["GET", `${REPO}/contents/data/../js/main.js?ref=main`], ["GET", `${REPO}/git/blobs/xyz`],
    ["GET", "/repos/someone/else/contents/data/curation.json?ref=main"], ["PATCH", `${REPO}/git/refs/heads/main`],
  ];
  for (const [method, path] of bad) {
    const r = await w.call(`/gh${path}`, { method, token: adminToken, body: method === "GET" || method === "DELETE" ? undefined : { message: "x", content: "" } });
    assert.equal(r.status, 403, `${method} ${path} -> ${r.status}`);
  }
  assert.equal(w.github.calls.length, before);
  const r = await w.call(`/gh${CUR}?ref=main`, { token: curatorToken, origin: "https://evil.example" });
  assert.equal(r.status, 403, "wrong origin");
  assert.equal(w.github.calls.length, before);
});

test("GitHub proxy: a PUT carries the bot token and names the signed-in curator", async () => {
  const sha = (await w.call(`/gh${CUR}?ref=main`, { token: curatorToken })).data.sha;
  const content = Buffer.from('{"version":1,"entries":{"1ABC":{"status":"rejected"}},"log":[]}\n').toString("base64");
  const before = w.github.calls.length;
  let r = await w.call(`/gh${CUR}`, { method: "PUT", token: curatorToken,
    body: { message: "Curator someone-else: reject 1ABC (not a membrane protein)\nCo-authored-by: x", content, branch: "main", sha,
      committer: { name: "Fake", email: "fake@example.com" }, author: { name: "Fake", email: "fake@example.com" } } });
  assert.equal(r.status, 200, JSON.stringify(r.data));
  const c = w.github.calls[before];
  assert.equal(c.method, "PUT");
  assert.equal(c.authorization, "Bearer bot-token");
  assert.equal(c.body.message, "Curator ana.lee: reject 1ABC (not a membrane protein) Co-authored-by: x");
  assert.deepEqual(Object.keys(c.body).sort(), ["branch", "content", "message", "sha"], "no author or committer: the bot makes the commit");
  assert.equal(c.body.branch, "main");
  // A stale sha comes back as 409 so the site re-reads and retries.
  r = await w.call(`/gh${CUR}`, { method: "PUT", token: curatorToken, body: { message: "Curator ana.lee: x", content, branch: "main", sha } });
  assert.equal(r.status, 409);
  // Another branch is refused.
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, body: { message: "x", content, branch: "gh-pages" } });
  assert.equal(r.status, 400);
  // Bodies over 2 MB are refused.
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, raw: JSON.stringify({ message: "x", content: "A".repeat(2 * 1024 * 1024 + 10) }) });
  assert.equal(r.status, 413);
  // A disabled curator cannot write.
  await w.call("/users/ana.lee/disable", { method: "POST", token: adminToken });
  const n = w.github.calls.length;
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, body: { message: "x", content } });
  assert.equal(r.status, 401);
  assert.match(r.data.error, /disabled/);
  assert.equal(w.github.calls.length, n);
  await w.call("/users/ana.lee/enable", { method: "POST", token: adminToken });
});

test("a refused bot token is reported as a service problem, not a sign-out", async () => {
  const w2 = await startWorker({ secrets: { ...SECRETS, GITHUB_TOKEN: "expired-token" }, githubAccepts: "bot-token" });
  try {
    await w2.call("/setup", { method: "POST", body: { ...ADMIN, setupCode: SECRETS.SETUP_CODE } });
    const t = (await w2.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password } })).data.token;
    const r = await w2.call(`/gh${CUR}?ref=main`, { token: t });
    assert.equal(r.status, 502);
    assert.match(r.data.error, /site owner/);
  } finally { await w2.mf.dispose(); }
});

test("the last administrator cannot be removed; a second admin can do it", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "second", password: "second-admin-1", role: "admin" } });
  const t2 = (await w.call("/login", { method: "POST", body: { username: "second", password: "second-admin-1" } })).data.token;
  let r = await w.call("/users/owner/disable", { method: "POST", token: t2 });
  assert.equal(r.status, 200);
  r = await w.call("/users/second/disable", { method: "POST", token: t2 });
  assert.equal(r.status, 400);
  await w.call("/users/owner/enable", { method: "POST", token: t2 });
});
