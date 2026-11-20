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
  const raw = JSON.parse((await w.sql("SELECT data FROM users WHERE username = 'owner'"))[0].data);
  assert.equal(raw.iterations, 20000, "fits the free plan's CPU limit");
  assert.equal(typeof raw.id, "string");
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
  const a = JSON.parse((await w.sql("SELECT data FROM users WHERE username = 'owner'"))[0].data).id;
  let r = await w.call("/me", { token: `${body}.${"A".repeat(43)}` });
  assert.equal(r.status, 401);
  r = await w.call("/me", { token: forge({ u: "owner", a, r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) + 600 }, "some-other-secret-value") });
  assert.equal(r.status, 401);
  r = await w.call("/me", { token: forge({ u: "owner", a, r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) - 5 }) });
  assert.equal(r.status, 401);
  assert.match(r.data.error, /expired/);
  r = await w.call("/me", { token: forge({ u: "owner", a, r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) + 600 }) });
  assert.equal(r.status, 200, "a correctly signed, unexpired token works");
  r = await w.call("/me", { token: forge({ u: "owner", r: "admin", v: 0, exp: Math.floor(Date.now() / 1000) + 600 }) });
  assert.equal(r.status, 401, "a token without the account id is refused");
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
  assert.equal((await w.sql("SELECT * FROM users WHERE username = 'temp'")).length, 0);
  r = await w.call("/users/temp", { method: "DELETE", token: adminToken });
  assert.equal(r.status, 404);
});

test("rate limit: 10 attempts per username and address, then 429 even with the right password", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "victim", password: "victim-password-1", role: "curator" } });
  const ip = "192.0.2.10";
  for (let i = 0; i < 10; i++) {
    const r = await w.call("/login", { method: "POST", body: { username: "victim", password: `guess-${i}-xxxxx` }, ip });
    assert.equal(r.status, 401);
  }
  const r = await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip });
  assert.equal(r.status, 429);
  assert.match(r.data.error, /^Too many failed sign-in attempts\. Try again in 15 minutes\.$/);
  assert.equal(r.data.retryMinutes, 15);
  // Someone else (another address) is not held up by it.
  assert.equal((await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip: freshIp() })).status, 200);
  // When the window has passed, the count starts again.
  await w.sql("UPDATE attempts SET reset = ?1 WHERE key = ?2", Math.floor(Date.now() / 1000) - 1, `ui:victim|${ip}`);
  assert.equal((await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip })).status, 200);
});

test("rate limit: parallel guesses are counted before the password is checked", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "racer", password: "racer-password-1", role: "curator" } });
  const ip = "192.0.2.20";
  const guesses = Array.from({ length: 60 }, (_, i) => w.call("/login", { method: "POST", body: { username: "racer", password: `guess-${i}-xxxxx` }, ip }));
  const right = w.call("/login", { method: "POST", body: { username: "racer", password: "racer-password-1" }, ip });
  const all = await Promise.all([...guesses, right]);
  const checked = all.filter((r) => r.status !== 429).length;
  assert.equal(checked, 10, `only 10 of 61 parallel attempts were checked (got ${checked})`);
  const row = (await w.sql("SELECT n FROM attempts WHERE key = ?1", `ui:racer|${ip}`))[0];
  assert.ok(row.n <= 10, `counter ${row.n}`);
});

test("rate limit: spread over many addresses, an account pauses for new devices but not for the curator's own browser", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "target", password: "target-password-1", role: "curator" } });
  // The curator signed in on this browser before: it holds a device token.
  const first = await w.call("/login", { method: "POST", body: { username: "target", password: "target-password-1" }, ip: freshIp() });
  assert.equal(first.status, 200);
  const device = first.data.device;
  assert.ok(device);
  // Strangers on 20 addresses make 20 wrong guesses.
  for (let i = 0; i < 20; i++) {
    assert.equal((await w.call("/login", { method: "POST", body: { username: "target", password: `guess-${i}-xxxxx` }, ip: `203.0.113.${100 + i}` })).status, 401);
  }
  // Another stranger, even with the right password, is paused.
  const r = await w.call("/login", { method: "POST", body: { username: "target", password: "target-password-1" }, ip: freshIp() });
  assert.equal(r.status, 429);
  // The curator, from a new address but on their own browser, gets in.
  const mine = await w.call("/login", { method: "POST", body: { username: "target", password: "target-password-1", device }, ip: freshIp() });
  assert.equal(mine.status, 200, JSON.stringify(mine.data));
  assert.equal(mine.data.device, device, "the device token is kept");
  // A device token for another account does not help.
  const other = await w.call("/login", { method: "POST", body: { username: "victim", password: "victim-password-1" }, ip: freshIp() });
  const r2 = await w.call("/login", { method: "POST", body: { username: "target", password: "target-password-1", device: other.data.device }, ip: freshIp() });
  assert.equal(r2.status, 429);
});

test("rate limit: 30 attempts per address across usernames; IPv6 counts per /64", async () => {
  const ip = "192.0.2.77";
  for (let i = 0; i < 30; i++) await w.call("/login", { method: "POST", body: { username: `user${i}`, password: "whatever-pass" }, ip });
  const r = await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password }, ip });
  assert.equal(r.status, 429);
  assert.equal((await w.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password }, ip: freshIp() })).status, 200);
  const { addressKey } = await import("../worker.js");
  assert.equal(addressKey("2001:db8:a:b:1:2:3:4"), "2001:db8:a:b::/64");
  assert.equal(addressKey("2001:0db8:000a:000b::9"), "2001:db8:a:b::/64");
  assert.equal(addressKey("2001:db8::1"), "2001:db8:0:0::/64");
  assert.equal(addressKey("198.51.100.4"), "198.51.100.4");
});

test("refused attempts past the limit write nothing to the database", async () => {
  const ip = "192.0.2.88";
  for (let i = 0; i < 12; i++) await w.call("/login", { method: "POST", body: { username: "victim", password: "nope-nope-nope" }, ip });
  const before = (await w.sql("SELECT n, reset FROM attempts WHERE key = ?1", `ui:victim|${ip}`))[0];
  assert.equal(before.n, 10);
  await w.call("/login", { method: "POST", body: { username: "victim", password: "nope-nope-nope" }, ip });
  assert.deepEqual((await w.sql("SELECT n, reset FROM attempts WHERE key = ?1", `ui:victim|${ip}`))[0], before);
});

test("a deleted account's sessions never come back, even if the name is reused", async () => {
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "gone", password: "gone-password-1", role: "curator" } });
  const old = (await w.call("/login", { method: "POST", body: { username: "gone", password: "gone-password-1", keep: true }, ip: freshIp() })).data.token;
  assert.equal((await w.call("/me", { token: old })).status, 200);
  assert.equal((await w.call("/users/gone", { method: "DELETE", token: adminToken })).status, 200);
  assert.equal((await w.call("/me", { token: old })).status, 401);
  await w.call("/users", { method: "POST", token: adminToken, body: { username: "gone", password: "another-pass-9", role: "admin" } });
  assert.equal((await w.call("/me", { token: old })).status, 401, "old token does not work for the new account");
  assert.equal((await w.call("/users", { token: old })).status, 401);
  assert.equal((await w.call("/users/gone", { method: "DELETE", token: adminToken })).status, 200);
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
  assert.equal(r.status, 200, "a plain visit says the service is ready");
  assert.match(r.data, /^READY/);
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
  // Content that is not plain base64 (for example, extra JSON fields smuggled in) is refused.
  const n0 = w.github.calls.length;
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, body: { message: "x", content: 'AAAA","author":{"name":"x","email":"x@example.com"},"x":"' } });
  assert.equal(r.status, 400);
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, body: { message: "x", content, sha: 'abc","author":"x' } });
  assert.equal(r.status, 400);
  assert.equal(w.github.calls.length, n0);
  // A large file (3 MB, about 3.5 times today's curation.json) saves.
  const fbSha = (await w.call(`/gh${FB}?ref=main`, { token: curatorToken })).data.sha;
  const big = Buffer.from(JSON.stringify({ version: 1, items: [], pad: "x".repeat(3 * 1024 * 1024) })).toString("base64");
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, body: { message: "big", content: big, branch: "main", sha: fbSha } });
  assert.equal(r.status, 200, JSON.stringify(r.data).slice(0, 200));
  assert.equal(w.github.calls.at(-1).body.content, big);
  // Bodies over 40 MB are refused.
  r = await w.call(`/gh${FB}`, { method: "PUT", token: curatorToken, raw: JSON.stringify({ message: "x", content: "A".repeat(40 * 1024 * 1024 + 10) }) });
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

test("the address page and /status say what is missing from the set-up", async () => {
  const w3 = await startWorker({ secrets: {}, d1: false });
  try {
    let r = await w3.call("/", { origin: null });
    assert.equal(r.status, 503);
    assert.match(r.data, /^NOT READY/);
    assert.match(r.data, /D1 database binding with the variable name DB/);
    assert.match(r.data, /SESSION_SECRET is missing/);
    assert.match(r.data, /GITHUB_TOKEN is missing/);
    r = await w3.call("/status");
    assert.equal(r.status, 503);
    assert.match(r.data.error, /not set up completely/);
    r = await w3.call("/login", { method: "POST", body: { username: "owner", password: "whatever-pass" } });
    assert.equal(r.status, 503);
  } finally { await w3.mf.dispose(); }
  const w4 = await startWorker({ secrets: { ...SECRETS, SESSION_SECRET: "short-secret", SETUP_CODE: "" }, githubAccepts: "other-token" });
  try {
    const r = await w4.call("/", { origin: null });
    assert.equal(r.status, 503);
    assert.match(r.data, /SESSION_SECRET is too short/);
    assert.match(r.data, /GitHub does not accept GITHUB_TOKEN/);
  } finally { await w4.mf.dispose(); }
  const w5 = await startWorker({ secrets: { ...SECRETS, SETUP_CODE: "" } });
  try {
    const r = await w5.call("/", { origin: null });
    assert.match(r.data, /SETUP_CODE is missing/);
  } finally { await w5.mf.dispose(); }
});

test("raising PBKDF2_ITERATIONS upgrades an account when it signs in", async () => {
  const w6 = await startWorker();
  try {
    await w6.call("/setup", { method: "POST", body: { ...ADMIN, setupCode: SECRETS.SETUP_CODE } });
    await w6.mf.setOptions({ modules: true, scriptPath: new URL("../worker.js", import.meta.url).pathname, compatibilityDate: "2024-11-01",
      d1Databases: ["DB"], bindings: { ...SECRETS, PBKDF2_ITERATIONS: "50000" }, outboundService: w6.github.handler });
    assert.equal((await w6.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password } })).status, 200);
    const raw = JSON.parse((await w6.sql("SELECT data FROM users WHERE username = 'owner'"))[0].data);
    assert.equal(raw.iterations, 50000);
    assert.equal((await w6.call("/login", { method: "POST", body: { username: "owner", password: ADMIN.password } })).status, 200);
  } finally { await w6.mf.dispose(); }
});
