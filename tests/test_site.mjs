import test from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';

// All scientific values below are test-only schema fixtures, never exported data.
const appURL = new URL('../site/app.mjs', import.meta.url);
const app = existsSync(appURL) ? await import(appURL) : {};
const figure = (src = 'plots/test-only.svg') => ({src, title: 'Fixture spectrum', caption: 'Test-only caption.', alt: 'Test-only spectral axes.'});
const passive = (kind = 'whistle', range_m = 100, surface = true) => ({
  id: `test-${kind}-${range_m}-${surface}`, kind, range_m, surface,
  title: 'Test-only passive case', summary: 'In-memory fixture, not a computed experiment.',
  metrics: [{label: 'Fixture pressure', value: 0.125, unit: 'Pa'}], figures: [figure()],
  downloads: [{href: 'data/test-only.json', label: 'Fixture metadata'}],
});
const active = (range_m = 150, snr_db = -18) => ({
  id: `test-active-${range_m}-${snr_db}`, range_m, snr_db,
  title: 'Test-only active case', summary: 'In-memory fixture.', metrics: [], figures: [figure()], downloads: [],
});
const catalog = () => ({
  schema_version: 1, title: 'Test-only OceanSignatures', version: 'test', provenance_note: 'In-memory test fixture only.',
  passive_cases: [passive()], active_cases: [active()], reference: null, gallery: [],
  verification: {status: 'not evaluated', summary: 'Test-only status.', checks: []},
  downloads: [], sources: [{id: 'test-source', title: 'Test reference', url: 'https://example.org/science'}],
});

test('catalog boundary accepts the frozen schema without manufacturing values', () => {
  assert.equal(typeof app.validateCatalog, 'function', 'validateCatalog is implemented');
  const input = catalog();
  assert.equal(app.validateCatalog(input), input);
  assert.equal(input.reference, null);
  assert.equal(input.verification.status, 'not evaluated');
});

test('catalog boundary rejects malformed children, ambiguous identities and unsafe assets', () => {
  const bad = [
    c => { c.schema_version = 2; },
    c => { delete c.provenance_note; },
    c => { c.passive_cases[0].range_m = NaN; },
    c => { c.passive_cases[0].surface = 'false'; },
    c => { c.passive_cases[0].kind = 'unknown'; },
    c => { c.passive_cases[0].metrics[0].value = Infinity; },
    c => { c.passive_cases[0].metrics[0].value = {}; },
    c => { delete c.passive_cases[0].figures[0].alt; },
    c => { c.passive_cases.push({...c.passive_cases[0], id: 'other-id'}); },
    c => { c.active_cases[0].id = c.passive_cases[0].id; },
    c => { c.active_cases[0].snr_db = '-18'; },
    c => { c.reference = {title: 'Incomplete'}; },
    c => { c.verification.checks = [{name: 'Missing detail', status: 'pass'}]; },
    c => { c.sources[0].url = 'javascript:alert(1)'; },
    c => { c.downloads = [{href: '/outside-prefix.json', label: 'Unsafe'}]; },
    c => { c.gallery = [figure('https://elsewhere.example/image.png')]; },
  ];
  for (const mutate of bad) {
    const input = catalog(); mutate(input);
    assert.throws(() => app.validateCatalog(input), /catalog/i, mutate.toString());
  }
});

test('asset paths remain relative to the GitHub Pages project prefix', () => {
  assert.equal(typeof app.safeAssetPath, 'function');
  for (const path of ['plots/pressure.svg', './data/case-100.json', 'plots/a%20b.png']) {
    assert.equal(app.safeAssetPath(path), path);
    assert.ok(new URL(path, 'https://example.org/oceansignatures/').pathname.startsWith('/oceansignatures/'));
  }
  for (const path of ['', '/plot.svg', '//host/plot.svg', '../secret', 'plots/../../x', 'plots/%2e%2e/x', 'plots/%252e%252e/x', 'plots\\x', 'https://host/x', 'data:x', 'plots/a.svg#x', 'plots/a.svg?x=1', 'a\nx']) {
    assert.throws(() => app.safeAssetPath(path), /catalog/i, path);
  }
});

test('reference and verification fields survive without invented calibration or counts', () => {
  const c = catalog();
  c.reference = {title: 'Test recording', summary: 'Test clip', calibration_status: 'uncalibrated digital units', source_url: 'https://example.org/recording', figures: [figure()], metrics: [{label: 'Test duration', value: '2.0', unit: 's'}], downloads: []};
  c.verification.checks.push({name: 'Test fixture check', status: 'pending', detail: 'Not measured.'});
  assert.equal(app.validateCatalog(c), c);
  assert.equal(c.reference.calibration_status, 'uncalibrated digital units');
});

test('every discrete selection retrieves its exact case without nearest-case fallback', () => {
  assert.equal(typeof app.selectCase, 'function');
  const c = catalog();
  c.passive_cases = ['whistle', 'clicks', 'vessel'].flatMap(kind => [25, 100, 400].flatMap(range => [false, true].map(surface => passive(kind, range, surface))));
  c.active_cases = [75, 150, 300].flatMap(range => [-18, -6].map(snr => active(range, snr)));
  app.validateCatalog(c);
  const selected = new Set();
  for (const item of c.passive_cases) {
    assert.equal(app.selectCase(c, 'passive', {kind: item.kind, range_m: item.range_m, surface: item.surface}), item);
    selected.add(item.id);
  }
  for (const item of c.active_cases) {
    assert.equal(app.selectCase(c, 'active', {range_m: item.range_m, snr_db: item.snr_db}), item);
    selected.add(item.id);
  }
  assert.equal(selected.size, 24);
  assert.equal(app.selectCase(c, 'passive', {kind: 'whistle', range_m: 99, surface: false}), null);
  assert.equal(app.selectCase(c, 'passive', {kind: 'whistle', range_m: '100', surface: true}), null);
  assert.equal(app.selectCase(null, 'active', {range_m: 150, snr_db: -18}), null);
  assert.throws(() => app.selectCase(c, 'typo', {}), /family/i);
});

test('rendered selection includes the original case identity, units, caption and download', () => {
  assert.equal(typeof app.renderCase, 'function');
  const item = passive();
  item.title = '<script>not markup</script>';
  item.metrics[0].label = 'Pressure < RMS';
  const html = app.renderCase(item, 'passive');
  assert.ok(html.includes(`data-case-id="${item.id}"`));
  for (const part of ['0.125', 'Pa', 'Fixture metadata', 'data/test-only.json', 'Test-only caption.', 'Test-only spectral axes.', '&lt;script&gt;not markup&lt;/script&gt;', 'Pressure &lt; RMS']) assert.ok(html.includes(part), part);
  assert.ok(!html.includes('<script>'));
  const unavailable = app.renderCase(null, 'passive');
  assert.match(unavailable, /No computed case matches/);
  assert.doesNotMatch(unavailable, /<img|data-case-id|0\.125/);
});

test('failed image markup removes image and open-plot link, retains honest context and retry', () => {
  assert.equal(typeof app.renderFigure, 'function');
  const fig = figure();
  const ready = app.renderFigure(fig);
  assert.match(ready, /loading="lazy"/);
  assert.match(ready, /data-figure-src="plots\/test-only.svg"/);
  const failed = app.renderFigure(fig, {failed: true});
  assert.doesNotMatch(failed, /<img|data-open-plot/);
  assert.match(failed, /Plot unavailable/);
  assert.match(failed, /data-retry-figure/);
  assert.match(failed, /Test-only caption/);
});

test('catalog loading clears previous results and exposes recoverable fetch/schema failures', async () => {
  assert.equal(typeof app.createCatalogLoader, 'function');
  const states = [];
  let next = catalog();
  let httpOK = true;
  const loader = app.createCatalogLoader({
    fetcher: async (url, options) => {
      assert.ok(String(url).endsWith('/site/data/catalog.json'));
      assert.equal(options.cache, 'no-store');
      return {ok: httpOK, status: 503, json: async () => next};
    }, onState: state => states.push(state),
  });
  await loader.load();
  assert.deepEqual(states.map(s => s.phase), ['loading', 'ready']);
  assert.equal(states[0].catalog, null);
  assert.equal(states[1].catalog, next);
  httpOK = false;
  await loader.load();
  assert.equal(states.at(-2).catalog, null);
  assert.equal(states.at(-1).catalog, null);
  assert.equal(states.at(-1).phase, 'error');
  assert.match(states.at(-1).error, /503/);
  httpOK = true; next = {schema_version: 999};
  await loader.load();
  assert.equal(states.at(-1).phase, 'error');
  assert.equal(states.at(-1).catalog, null);
  next = catalog(); await loader.load();
  assert.equal(states.at(-1).phase, 'ready');
  loader.dispose();
});

test('late responses cannot replace a newer selection catalogue', async () => {
  assert.equal(typeof app.createCatalogLoader, 'function');
  const states = [];
  let resolveFirst;
  let n = 0;
  const newer = catalog(); newer.version = 'newer';
  const loader = app.createCatalogLoader({fetcher: () => ++n === 1 ? new Promise(resolve => {resolveFirst = resolve;}) : Promise.resolve({ok: true, json: async () => newer}), onState: state => states.push(state)});
  const first = loader.load();
  await loader.load();
  resolveFirst({ok: true, json: async () => catalog()});
  await first;
  assert.equal(states.at(-1).catalog.version, 'newer');
  assert.equal(states.filter(s => s.phase === 'ready').length, 1);
  loader.dispose();
});

test('a stalled request has a visible timeout instead of indefinite loading', async () => {
  assert.equal(typeof app.createCatalogLoader, 'function');
  const states = [];
  const loader = app.createCatalogLoader({fetcher: () => new Promise(() => {}), onState: s => states.push(s), timeoutMs: 10});
  await loader.load();
  assert.equal(states.at(-1).phase, 'error');
  assert.match(states.at(-1).error, /timed out/i);
  loader.dispose();
});

test('reference and evidence renderers preserve provenance and do not imply validation', () => {
  assert.equal(typeof app.renderReference, 'function');
  assert.equal(typeof app.renderVerification, 'function');
  assert.equal(typeof app.renderSources, 'function');
  assert.match(app.renderReference(null), /No reference recording is available/);
  assert.doesNotMatch(app.renderReference(null), /<img/);
  const reference = {title: 'Test recording', summary: 'Test clip', calibration_status: 'Uncalibrated digital units', source_url: 'https://example.org/recording', figures: [figure()], metrics: [], downloads: []};
  const html = app.renderReference(reference);
  assert.match(html, /Uncalibrated digital units/);
  assert.match(html, /https:\/\/example.org\/recording/);
  assert.match(html, /not absolute pressure/i);
  const c = catalog();
  c.verification.checks = [{name: '<check>', status: 'pending', detail: 'No result yet'}];
  const verification = app.renderVerification(c.verification);
  assert.match(verification, /not evaluated/);
  assert.match(verification, /&lt;check&gt;/);
  assert.doesNotMatch(verification, /all tests passed|validated system/i);
  assert.match(app.renderSources(c.sources), /https:\/\/example.org\/science/);
});

test('all data-owned panels lose stale plots, values and source links on reload or failure', () => {
  assert.equal(typeof app.renderCatalogView, 'function');
  const c = catalog(); c.gallery = [figure()];
  const selections = {passive: {kind: 'whistle', range_m: 100, surface: true}, active: {range_m: 150, snr_db: -18}};
  const ready = app.renderCatalogView({phase: 'ready', catalog: c}, selections);
  assert.match(ready.passive, /data-case-id/);
  assert.match(ready.gallery, /<img/);
  for (const phase of ['loading', 'error']) {
    const view = app.renderCatalogView({phase, catalog: c, error: '<network>'}, selections);
    for (const key of ['passive', 'active', 'reference', 'gallery', 'verification', 'downloads', 'sources', 'provenance']) {
      assert.doesNotMatch(view[key], /<img|0\.125|test-only.svg|example.org|data-case-id/);
    }
    assert.equal(view.disabled, true);
    assert.ok(view.status.includes(phase === 'loading' ? 'Loading' : '&lt;network&gt;'));
  }
});

test('page shell provides labelled discrete controls, learning anchors and local assets', () => {
  const url = new URL('../site/index.html', import.meta.url);
  assert.ok(existsSync(url), 'scientific page shell exists');
  const html = readFileSync(url, 'utf8');
  for (const id of ['main', 'learn', 'passive', 'active', 'spectra', 'reference', 'evidence', 'sources', 'catalog-status', 'catalog-retry', 'passive-kind', 'passive-range', 'passive-surface', 'active-range', 'active-snr', 'passive-results', 'active-results']) assert.ok(html.includes(`id="${id}"`), id);
  for (const id of ['passive-kind', 'passive-range', 'passive-surface', 'active-range', 'active-snr']) assert.ok(html.includes(`for="${id}"`), `label ${id}`);
  assert.match(html, /<math/);
  assert.match(html, /source.*channel.*receiver/is);
  assert.match(html, /precomputed/);
  assert.match(html, /verification/i);
  assert.match(html, /validation/i);
  assert.match(html, /<noscript>/);
  assert.match(html, /type="module" src="\.\/app.mjs"/);
  assert.doesNotMatch(html, /(?:src|href)="https?:[^"\s]+\.(?:js|css)|cdn\./i);
  assert.equal(typeof app.mount, 'function', 'browser entry point is wired');
  const cssURL = new URL('../site/styles.css', import.meta.url);
  assert.ok(existsSync(cssURL), 'responsive style sheet exists');
  const css = readFileSync(cssURL, 'utf8');
  assert.match(css, /min-height:\s*44px/);
  assert.match(css, /prefers-reduced-motion/);
  assert.match(css, /:focus-visible/);
  assert.match(css, /@media/);
});

test('technical explanations distinguish controlled SNR, receiver delay and recording calibration', () => {
  const html = readFileSync(new URL('../site/index.html', import.meta.url), 'utf8');
  for (const phrase of ['fixed-SNR range', 'known simulation noise', 'finite emission', 'receiver group delay', 'recording-pipeline']) assert.ok(html.includes(phrase), phrase);
});

test('tablet geometry stays full width and masthead links retain 44px targets', () => {
  const css = readFileSync(new URL('../site/styles.css', import.meta.url), 'utf8');
  const tablet = css.split('@media (max-width: 820px)')[1].split('@media (max-width: 560px)')[0];
  assert.match(tablet, /\.hero\s*\{[^}]*grid-template-columns:\s*1fr\s*;/);
  assert.match(css, /\.masthead a\s*\{[^}]*min-width:\s*44px/);
});

test('metric display rounds numerical noise without changing exported values', () => {
  const metrics = [{label: 'SNR', value: 29.860244829523342, unit: 'dB'},
    {label: 'Range', value: 149.9375, unit: 'm'},
    {label: 'Status', value: 'Not calibrated', unit: ''}];
  const original = JSON.stringify(metrics);
  const html = app.renderMetrics(metrics);
  assert.match(html, /data-raw-value="29\.860244829523342">29\.860245/);
  assert.match(html, />149\.9375/);
  assert.match(html, /Not calibrated/);
  assert.equal(JSON.stringify(metrics), original);
});

test('responsive headings retain word boundaries at authored line breaks', () => {
  const css = readFileSync(new URL('../site/styles.css', import.meta.url), 'utf8');
  assert.doesNotMatch(css, /h2\s+br\s*\{[^}]*display:\s*none/,
    'hiding a br joins adjacent words and overflows the 320px layout');
});
