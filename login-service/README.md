# Login service: username and password sign-in for curators

Curators normally sign in on the Curators tab with a personal GitHub token.
This small Cloudflare Worker lets them use a **username and password**
instead. The Worker keeps one GitHub token, checks passwords, and makes only
the GitHub calls the site needs: reading and saving `data/curation.json` and
`data/feedback.json`. It runs on Cloudflare's free plan.

Until the Worker's address is in `js/main.js`, the site shows only the GitHub-token sign-in, as before.

## Setup checklist (about 20 minutes, once)

### 1. Make the GitHub token for the Worker
1. Signed in to GitHub as **andyyjiang**, open <https://github.com/settings/personal-access-tokens/new> (a *fine-grained* token).
2. **Token name:** `mpstruc login service`.
3. **Expiration:** choose *Custom* and a date about a year away (or the longest the page offers). Put a reminder in your calendar a week before that date. If you choose *No expiration*, no reminder is needed.
4. **Resource owner:** andyyjiang.
5. **Repository access:** *Only select repositories*, then pick **team-website**.
6. **Permissions:** add the **Contents** repository permission and set it to **Read and write**. Add nothing else. (GitHub adds *Metadata: Read-only* by itself; that is expected.)
7. Click **Generate token** and copy it (it starts with `github_pat_`). Keep the tab open for step 5.

### 2. Make a free Cloudflare account
Sign up at <https://dash.cloudflare.com/sign-up> and confirm your email address.

### 3. Create the Worker (in the dashboard, no command line)
1. In the Cloudflare dashboard, open **Workers & Pages** (under *Compute*), click **Create**, then **Start with Hello World!**.
2. Name it `mpstruc-login`, then click **Deploy**.
3. Click **Edit code**. Delete everything in the editor, paste in the whole of [`worker.js`](worker.js), and click **Deploy**.
4. Copy the Worker's address shown at the top, for example `https://mpstruc-login.YOURNAME.workers.dev`.

### 4. Create the database and connect it
1. In the dashboard, open **Storage & Databases > D1 SQL Database**, click **Create** (or *Create database*), name it `mpstruc-users`, and click **Create**. You do not need to add any tables; the Worker does that itself.
2. Go back to the `mpstruc-login` Worker: **Settings > Bindings > Add > D1 database**. **Variable name:** `DB` (capitals). **D1 database:** `mpstruc-users`. Click **Deploy** (or **Save**).

### 5. Add the three secrets
In the Worker, under **Settings**, open **Variables and Secrets** and click **Add** three times. Each time, set **Type** to **Secret**:

| Variable name    | Value |
|------------------|-------|
| `GITHUB_TOKEN`   | the `github_pat_…` token from step 1 |
| `SESSION_SECRET` | a random string of at least 32 characters (for example, mash the keyboard) |
| `SETUP_CODE`     | a one-off code that you choose, for example `blue-heron-7392`. You need it once, in step 8 |

Click **Deploy** after the last one.

### 6. Check the Worker
Open the Worker's address in a browser. The first word should be **READY**. If it says **NOT READY**, the page lists exactly what is missing or wrong (for example, a mistyped binding name or token). Fix it, click **Deploy** in Cloudflare, and reload the page until it says READY.

### 7. Tell the site the Worker's address
Ask Claude: *"Set LOGIN_SERVICE_URL to https://mpstruc-login.YOURNAME.workers.dev"*. Or do it yourself on GitHub:
open `js/main.js`, click the pencil, find this line

```js
  const LOGIN_SERVICE_URL = "";
```

and paste the address between the quotes (no `/` at the end):

```js
  const LOGIN_SERVICE_URL = "https://mpstruc-login.YOURNAME.workers.dev";
```

Then commit the change (through a pull request, as usual). The site updates a minute or two after it reaches `main`.

### 8. Create the first administrator and the curators
Open <https://andyyjiang.github.io/team-website/mpstruc.html#/curators>. You should see **Create the first administrator**.
Enter the setup code from step 5, choose your username and a password (at least 10 characters), and click **Create administrator**.
After that, the setup form is gone for good. Under **Manage curators**, add each curator (username, display name, temporary password, role) and send them their username and temporary password privately. They can change it under **Change my password**.

### 9. Final check: one real save
Still signed in with your username and password, make one decision you were going to make anyway (or leave a piece of feedback) and check that it appears in the **Activity log**.
If instead the page says *The sign-in service could not be reached* while the Worker's address still says READY, open the Worker in Cloudflare, click **Logs** (or *Observability*), and try again. If the log shows **Exceeded CPU time limit**, the free plan is too small: switch the account to **Workers Paid** ($5 a month, under *Workers & Pages > Plans*). This has not happened in testing, but it has not been tried on a real Cloudflare account yet.

## Day to day
- **Someone forgot their password:** an administrator clicks **Reset password** next to their name and sends them the new one.
- **Someone leaves:** **Disable** (they can be enabled again) or **Delete**. Either way they are signed out at once. A deleted account's old sign-ins never work again, even if someone is later given the same username.
- **The GitHub token expires:** saving fails with a message asking the site owner to renew it. Make a new token (step 1), then, in the Worker, **Settings > Variables and Secrets**, edit `GITHUB_TOKEN`, paste it, and deploy. The Worker's address then says READY again.
- **"Too many failed sign-in attempts":** sign-in pauses for up to 15 minutes. A curator who has signed in on their own browser before is not held up by other people's failed attempts on their account; on a new browser they may have to wait.
- **All administrators are locked out:** in Cloudflare, open **Storage & Databases > D1 SQL Database > mpstruc-users > Console** and run `DELETE FROM users WHERE username = 'their-username';` for each administrator. When no administrator is left, the setup form comes back (with the same `SETUP_CODE`).
- **Usernames are public.** Decisions in `data/curation.json` record who made them, and that file is on GitHub Pages. Choose usernames that differ from curators' GitHub logins, so the two kinds of sign-in are never confused in the Curators table. Passwords are what protects the accounts: at least 10 characters, never reused from elsewhere.

## Free plan limits
Cloudflare's free plan allows 100,000 requests a day to the Worker, 10 ms of CPU time per request, and in D1 100,000 row writes and 5 million row reads a day. The site uses a tiny fraction of that. If a limit is ever reached (for example, someone floods the sign-in page), sign-in and saving stop until the next day (UTC) and the Worker's address says what failed. GitHub-token sign-in keeps working regardless.

## How it works
- Accounts live in the D1 database. Passwords are stored only as salted PBKDF2-SHA256 hashes (20,000 rounds, sized for the free plan's CPU limit; each account records its own round count). On Workers Paid you can add a plain-text variable `PBKDF2_ITERATIONS` = `100000`: each account is re-hashed with the new count the next time it signs in.
- Signing in gives a signed session token that lasts 12 hours, or 30 days with *Keep me signed in*. Disabling, deleting, resetting a password or changing your own password ends that person's other sessions, at once (D1 is consistent everywhere, unlike Workers KV).
- Every sign-in attempt is counted before the password is checked, in one atomic database step, so sending guesses in parallel does not get around the limits. Per 15 minutes: 10 attempts per username from one internet address; 30 from one address in all; 20 per username from browsers that have not signed in to it before; 10 per username from a browser that has (it carries a *device token*). Successful sign-ins do not count. IPv6 addresses count per /64.
- The Worker answers only `https://andyyjiang.github.io` (and `http://localhost` for testing). Its GitHub proxy (`/gh/…`) allows only reading and saving the two data files (and reading large files as blobs), on `main`. Saves are commits made with the Worker's token, and the message starts with `Curator <username>:`. Decisions in the record carry the username (`by`) and display name (`byName`). Saves of up to 40 MB are accepted (today's `curation.json` is under 1 MB).
- Endpoints: `GET /` (set-up check, for people), `POST /login`, `GET /me`, `POST /logout`, `POST /password`, `GET /status`, `POST /setup`. For administrators: `GET /users`, `POST /users`, `POST /users/<name>/reset`, `POST /users/<name>/disable`, `POST /users/<name>/enable`, `DELETE /users/<name>`.

## Deploying from the command line instead (optional)
```sh
cd login-service
npx wrangler login
npx wrangler d1 create mpstruc-users        # paste the database_id it prints into wrangler.toml
npx wrangler secret put GITHUB_TOKEN
npx wrangler secret put SESSION_SECRET
npx wrangler secret put SETUP_CODE
npx wrangler deploy
```

## Tests and local trial (for developers)
```sh
cd login-service
npm install
npm test        # runs worker.js in Miniflare (the Workers runtime) against a fake GitHub
npm run dev     # the service on http://localhost:8787, with a fake GitHub seeded from data/
```
With `npm run dev` running, serve the repository (`python3 -m http.server 8000`) and open
`http://localhost:8000/mpstruc.html?login-service=http://localhost:8787#/curators` (setup code `open-sesame-42`).
The `login-service` parameter works only when the page itself is opened from localhost.
Miniflare does not enforce Cloudflare's CPU limits, so the final check in step 9 matters.
