# Project guide for Claude

A website for our club/business, built by two people who each use their own Claude.
Follow these rules so both of our sessions produce consistent work.

## Stack
- Plain HTML, CSS, and vanilla JavaScript. No frameworks, no build step, no npm.
- Hosted on GitHub Pages from the `main` branch, root folder.

## Structure
- `index.html`: home page. New pages go at the root as `about.html`, `events.html`, etc.
- `css/style.css`: all styles. Colors and fonts live as variables in `:root` at the top; use them instead of hard-coded values.
- `js/main.js`: all scripts.
- `images/`: images. Keep files under ~500 KB and use descriptive lowercase names (`team-photo.jpg`).

## Conventions
- Every page uses the same `<header>` and `<footer>` as `index.html`, so copy them exactly when adding a page and update the nav on every page.
- Layout must work on phones (test at ~375px wide) and in dark mode.
- Every `<img>` needs meaningful `alt` text.
- Use relative links (`css/style.css`, not `/css/style.css`) so the site works on GitHub Pages.

## Running locally
Open `index.html` in a browser, or run `python3 -m http.server` in this folder and visit http://localhost:8000.

## Workflow
- Never commit directly to `main`. Make a branch per change (`add-events-page`, `fix-mobile-nav`), then open a pull request.
- Keep each pull request to one feature so it's easy for the other person to review.
- Split work by page or section to avoid editing the same lines at the same time.

## Sharing session files
Cloud sessions are private and temporary, so anything not pushed to GitHub is lost and invisible to the rest of the team.
- Save files the user uploads to `shared/uploads/` and files you generate (reports, exports, data) to `shared/outputs/`, named `YYYY-MM-DD-topic.ext`. Code and site pages go in their normal place in the project.
- Before ending work in a session, commit everything on the current branch and push it, then tell the user the branch name so teammates can find it.
- This repository is public. If an upload looks private (personal details, unpublished data, credentials), ask the user before committing it.
- If you publish an Artifact page, add its link to `shared/outputs/artifacts.md` with a one-line description.
