#!/usr/bin/env node
// Read the other group's demo (https://mpstruc-viewer-demo.vercel.app/) the way a
// visitor does: open the rendered page in a headless browser, read the visible
// group cards, and click "Load 50 more groups" until we are past the weeks we
// want. Nothing else is touched (no API calls, no keys, no bundle reading).
//
//   node read_page.js [--weeks 8] [--all] [--out groups.jsonl] [--wait 1500]
//   node read_page.js --parse-text page.txt     (parse a saved page text; no browser)
//
// Output: one JSON object per group (newest first), then {"_meta": {...}}.
// Merge into the snapshot with: python -m mpstruc.vercel_import page groups.jsonl
'use strict';

const fs = require('fs');

const URL = 'https://mpstruc-viewer-demo.vercel.app/';
process.stdout.on('error', (e) => { if (e.code === 'EPIPE') process.exit(0); throw e; });

function parseArgs(argv) {
  const a = { weeks: 8, all: false, out: null, wait: 1500, url: URL, parseText: null, maxClicks: 400 };
  for (let i = 0; i < argv.length; i++) {
    const k = argv[i];
    if (k === '--weeks') a.weeks = Number(argv[++i]);
    else if (k === '--all') a.all = true;
    else if (k === '--out') a.out = argv[++i];
    else if (k === '--wait') a.wait = Number(argv[++i]);
    else if (k === '--url') a.url = argv[++i];
    else if (k === '--max-clicks') a.maxClicks = Number(argv[++i]);
    else if (k === '--parse-text') a.parseText = argv[++i];
    else if (k === '-h' || k === '--help') { console.log(fs.readFileSync(__filename, 'utf8').split('\n').slice(1, 11).join('\n')); process.exit(0); }
    else throw new Error('unknown argument ' + k);
  }
  if (!(a.weeks > 0)) throw new Error('--weeks must be a positive number');
  return a;
}

// ---- text parsing (mirrors mpstruc/vercel_import.py parse_group) ----
const FIELD = /^(Title|Release Date\(s\)|Number of TM segment|Resolution|UniProt ID\(s\)|Classification\(s\)|Membrane Protein Verdict):\s*(.*)$/;
const PRED = /^(.+? model on [^:]+):\s*(.*?)\s+[—-]\s+([\d.]+)%$/;
const PDB = /^[0-9][A-Z0-9]{3}$|^PDB_[0-9A-Z]{8}$/;
const blank = (v) => { v = (v || '').trim(); return v || null; };
const list = (v) => v.split(/,\s*/).filter(Boolean);

function parseGroup(text) {
  let lines = text.split('\n').map((s) => s.trim()).filter(Boolean);
  // Lines above the label (a "#1" badge, say) are not part of the group.
  const at = lines.findIndex((s) => s.startsWith('PubMed Group:'));
  if (at < 0) return null;
  lines = lines.slice(at);
  const g = { group: blank(lines[0].split(':').slice(1).join(':')), journal: null, paperTitle: null, doi: null,
    title: null, pdbIds: [], releaseDates: [], tm: null, resolution: null, uniprot: [], classifications: [],
    predictions: [], verdict: null };
  let i = 1;
  let head = [];
  while (i < lines.length && !lines[i].startsWith('Title:') && lines[i] !== 'PDB IDs:') head.push(lines[i++]);
  if (head[0] === 'Journal:') {
    head = head.slice(1);
    if (head.length) { g.journal = head[0] === 'N/A' ? null : head[0]; head = head.slice(1); }
  }
  const urls = head.filter((h) => h.startsWith('http'));
  g.doi = urls[0] || null;
  g.paperTitle = head.filter((h) => !h.startsWith('http')).join(' ') || null;
  while (i < lines.length) {
    const ln = lines[i];
    if (ln === 'PDB IDs:') {
      i++;
      while (i < lines.length && PDB.test(lines[i].toUpperCase())) g.pdbIds.push(lines[i++].toUpperCase());
      continue;
    }
    const m = FIELD.exec(ln);
    if (m) {
      const key = m[1]; const val = m[2].trim();
      if (key === 'Title') g.title = val || null;
      else if (key === 'Release Date(s)') g.releaseDates = val.match(/\d{4}-\d{2}-\d{2}/g) || [];
      else if (key === 'Number of TM segment') g.tm = val || null;
      else if (key === 'Resolution') g.resolution = val || null;
      else if (key === 'UniProt ID(s)') g.uniprot = list(val).filter((u) => u !== 'N/A');
      else if (key === 'Classification(s)') g.classifications = list(val);
      else g.verdict = val || null;
    } else {
      const p = PRED.exec(ln);
      if (p) g.predictions.push({ model: p[1], label: p[2], pct: Number(p[3]) });
    }
    i++;
  }
  return g;
}

function parsePageText(text) {
  return text.split(/^(?=\s*PubMed Group:)/m).map(parseGroup).filter(Boolean);
}

// ---- browser ----
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// Visible text of each group card from index `from` on. A card is the largest
// ancestor of a "PubMed Group:" label that holds no other such label.
function cardTexts(from) {
  const labels = [...document.querySelectorAll('strong, b, span, p, div')]
    .filter((e) => e.childElementCount === 0 && e.textContent.trim() === 'PubMed Group:');
  const count = (el) => (el.innerText.match(/PubMed Group:/g) || []).length;
  return labels.slice(from).map((el) => {
    let card = el;
    while (card.parentElement && count(card.parentElement) === 1) card = card.parentElement;
    return card.innerText;
  });
}

async function main() {
  const a = parseArgs(process.argv.slice(2));
  const out = a.out ? fs.createWriteStream(a.out) : process.stdout;
  const emit = (o) => out.write(JSON.stringify(o) + '\n');

  if (a.parseText) {
    const groups = parsePageText(fs.readFileSync(a.parseText, 'utf8'));
    groups.forEach(emit);
    emit({ _meta: { source: a.parseText, groups: groups.length } });
    if (a.out) out.end();
    return;
  }

  const cutoff = new Date(Date.now() - a.weeks * 7 * 86400e3).toISOString().slice(0, 10);
  const { chromium } = require('playwright');
  const browser = await chromium.launch();
  const page = await browser.newPage();
  const readAt = new Date().toISOString();
  let total = null; let clicks = 0; let stopReason = 'end of list';
  const groups = [];
  // Cards read so far (labels seen), counted apart from the groups parsed, so a card that
  // does not parse is not read again, and a layout we cannot parse stops the reader early
  // instead of clicking through the whole list.
  let seen = 0; let unparsed = 0; let failure = null;
  try {
    await page.goto(a.url, { waitUntil: 'domcontentloaded', timeout: 90000 });
    await page.getByText('PubMed Group:').first().waitFor({ timeout: 90000 });
    await sleep(a.wait);
    for (;;) {
      const texts = await page.evaluate(cardTexts, seen);
      seen += texts.length;
      for (const t of texts) { const g = parseGroup(t); if (g && g.pdbIds.length) { groups.push(g); emit(g); } else unparsed++; }
      if (seen && !groups.length) { failure = `none of the ${seen} group cards on the page could be read (has the layout changed?)`; break; }
      if (unparsed > Math.max(2, seen / 10)) { failure = `${unparsed} of ${seen} group cards could not be read (has the layout changed?)`; break; }
      const btn = page.getByRole('button', { name: /Load \d+ more groups/ });
      const label = (await btn.count()) ? await btn.first().innerText() : null;
      const m = label && /\((\d+)\s*\/\s*(\d+)\)/.exec(label);
      if (m) total = Number(m[2]);
      process.stderr.write(`groups ${groups.length}${unparsed ? ` (${unparsed} unreadable)` : ''}${total ? '/' + total : ''}\n`);
      if (!a.all && groups.length) {
        // Stop once the last group read was released entirely before the cutoff:
        // the list is newest first, so everything after it is older still.
        const last = groups[groups.length - 1];
        if (last.releaseDates.length && last.releaseDates.every((d) => d < cutoff)) { stopReason = `past ${a.weeks} weeks (before ${cutoff})`; break; }
      }
      if (!label) break;
      if (++clicks > a.maxClicks) { stopReason = 'max clicks'; break; }
      const before = groups.length;
      await btn.first().click();
      // Wait politely for the next cards to render.
      await page.waitForFunction((n) => (document.body.innerText.match(/PubMed Group:/g) || []).length > n,
        before, { timeout: 60000, polling: 500 }).catch(() => {});
      await sleep(a.wait);
    }
  } finally {
    await browser.close();
  }
  if (!failure && !groups.length) failure = 'no group cards were found on the page';
  emit({ _meta: { source: a.url, readAt, weeks: a.all ? null : a.weeks, cutoff: a.all ? null : cutoff,
    groups: groups.length, cards: seen, unreadable: unparsed, pdbs: groups.reduce((s, g) => s + g.pdbIds.length, 0),
    total, clicks, stopReason: failure ? `failed: ${failure}` : stopReason } });
  if (a.out) out.end();
  if (failure) {
    console.error(`read_page.js: ${failure}`);
    process.exitCode = 2;
  }
}

module.exports = { parseGroup, parsePageText };
if (require.main === module) main().catch((e) => { console.error(e); process.exit(1); });
