# Login service: username and password sign-in for curators

Curators normally sign in on the Curators tab with a personal GitHub token.
This small Cloudflare Worker lets them use a **username and password**
instead. The Worker keeps one GitHub token, checks passwords, and makes only
the GitHub calls the site needs: reading and saving `data/curation.json` and
`data/feedback.json`. It runs on Cloudflare's free plan.

Until the Worker's address is in `js/main.js`, the site shows only the GitHub-token sign-in, as before.

## Setup checklist (about 20 minutes, once)

### 1. Make the GitHub token for the Worker
1. Signed in to GitHub as **andyyjiang**, open <https://github.com/settings/personal-access-tokens/new>.
2. **Token name:** `mpstruc login service`. **Expiration:** 1 year, or the longest offered. Put a reminder in your calendar to renew it.
3. **Resource owner:** andyyjiang.
4. **Repository access:** *Only select repositories*, then pick **team-website**.
5. **Permissions > Repository permissions > Contents:** *Read and write*. Leave everything else as it is. GitHub adds *Metadata: Read-only* by itself.
6. Click **Generate token** and copy it (it starts with `github_pat_`). Keep it open in a tab for step 4.

### 2. Make a free Cloudflare account
Sign up at <https://dash.cloudflare.com/sign-up> and confirm your email address.

### 3. Create the Worker (in the dashboard, no command line)
1. In the Cloudflare dashboard, open **Workers & Pages** (under *Compute*), click **Create**, then **Start with Hello World!**.
2. Name it `mpstruc-login`, then click **Deploy**.
3. Click **Edit code**. Delete everything in the editor, paste in the whole of [`worker.js`](worker.js), and click **Deploy**.
4. Copy the Worker's address shown at the top, for example `https://mpstruc-login.YOURNAME.workers.dev`.

### 4. Add storage and the three secrets
1. Open **Storage & Databases > KV** (Workers KV), click **Create** (or *Create instance*), name it `mpstruc-users`, and click **Create**.
2. Go back to the `mpstruc-login` Worker: **Settings > Bindings > Add > KV namespace**. **Variable name:** `USERS` (capitals). **KV namespace:** `mpstruc-users`. Click **Deploy** (or **Save**).
3. Still under **Settings**, open **Variables and Secrets** and click **Add** three times. Each time, set **Type** to **Secret**:

   | Variable name    | Value |
   |------------------|-------|
   | `GITHUB_TOKEN`   | the `github_pat_…` token from step 1 |
   | `SESSION_SECRET` | a long random string: at least 32 characters (for example, mash the keyboard) |
   | `SETUP_CODE`     | a one-off code that you choose, for example `blue-heron-7392`. You need it once, in step 6 |

   Click **Deploy** after the last one.
4. To check it runs, open the Worker's address in a browser. It should say *The mpstruc login service is running.*

### 5. Tell the site the Worker's address
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

### 6. Create the first administrator
Open <https://andyyjiang.github.io/team-website/mpstruc.html#/curators>. You should see **Create the first administrator**.
Enter the setup code from step 4, choose your username and a password (at least 10 characters), and click **Create administrator**.
After that, the setup form is gone for good. Under **Manage curators**, add each curator (username, display name, temporary password, role) and send them their username and temporary password privately. They can change it under **Change my password**.

## Day to day
- **Someone forgot their password:** an administrator clicks **Reset password** next to their name and sends them the new one.
- **Someone leaves:** **Disable** (they can be enabled again) or **Delete**. Either way they are signed out at once.
- **The GitHub token expires:** saving fails with a message asking the site owner to renew it. Make a new token (step 1), then, in the Worker, **Settings > Variables and Secrets**, edit `GITHUB_TOKEN`, paste it, and deploy.
- **All administrators are locked out:** in Cloudflare, open the `mpstruc-users` KV namespace and delete the `user:…` entries of the administrators. When no administrator is left, the setup form comes back (with the same `SETUP_CODE`).

## How it works
- Accounts live in the KV namespace (`user:<username>`). Passwords are stored only as salted PBKDF2-SHA256 hashes (100,000 rounds, the most Cloudflare allows).
- Signing in gives a signed session token that lasts 12 hours, or 30 days with *Keep me signed in*. Disabling an account, resetting its password, or changing your own password ends that person's other sessions.
- After 10 failed sign-ins in 15 minutes for one username, or from one internet address, sign-in pauses for the rest of those 15 minutes.
- The Worker answers only `https://andyyjiang.github.io` (and `http://localhost` for testing). Its GitHub proxy (`/gh/…`) allows only reading and saving the two data files (and reading large files as blobs), on `main`. Saves are commits made with the Worker's token, and the message starts with `Curator <username>:`. Decisions in the record carry the username (`by`) and display name (`byName`).
- Endpoints: `POST /login`, `GET /me`, `POST /logout`, `POST /password`, `GET /status`, `POST /setup`. For administrators: `GET /users`, `POST /users`, `POST /users/<name>/reset`, `POST /users/<name>/disable`, `POST /users/<name>/enable`, `DELETE /users/<name>`.

## Deploying from the command line instead (optional)
```sh
cd login-service
npx wrangler login
npx wrangler kv namespace create USERS      # paste the id it prints into wrangler.toml
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
