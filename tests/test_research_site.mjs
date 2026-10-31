import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, readFileSync} from 'node:fs';

const moduleURL = new URL('../site/research.mjs', import.meta.url);
const research = existsSync(moduleURL) ? await import(moduleURL) : {};
// Explicitly in-memory schema fixtures, never written to public result paths.
const study = () => ({
  id: 'test-passage', title: 'Test-only passage', summary: 'A schema fixture, not a result.',
  evidence_kind: 'predicted_acoustic_with_synthetic_observation',
  metrics: [{label: 'Test-only residual', value: 6.000000000000001, unit: 'dB'}],
  figures: [{src: 'data/research/test.png', title: 'Test-only figure', caption: 'Fixture caption', alt: 'Fixture axes'}],
  downloads: [{href: 'data/research/test.json', label: 'Fixture data'}],
  limitations: ['This fixture establishes no scientific result.'],
});
const catalog = () => ({schema_version: 1, project_version: '0.1.0', studies: [study()]});

test('research catalog keeps project identity and evidence classes explicit', () => {
  assert.equal(typeof research.validateResearch, 'function', 'research schema boundary is missing');
  const input = catalog();
  assert.equal(research.validateResearch(input), input);
  assert.equal(input.project_version, '0.1.0');
});

test('research catalog rejects ambiguous, unsafe and unsupported computed results', () => {
  for (const mutate of [
    c => {c.project_version = '2.0';},
    c => {c.schema_version = 0;},
    c => {c.studies = [];},
    c => {c.studies.push(study());},
    c => {c.studies[0].evidence_kind = 'field validated';},
    c => {c.studies[0].metrics[0].value = Infinity;},
    c => {c.studies[0].figures[0].src = '../outside.png';},
    c => {c.studies[0].downloads[0].href = 'https://elsewhere.invalid/data';},
    c => {c.studies[0].limitations = [];},
    c => {delete c.studies[0].figures[0].alt;},
  ]) {
    const input = catalog(); mutate(input);
    assert.throws(() => research.validateResearch(input), /research|catalog/i, mutate.toString());
  }
});

test('research rendering preserves exact identity, raw values, units and boundaries', () => {
  assert.equal(typeof research.renderStudy, 'function', 'research rendering is missing');
  const input = study(); input.title = '<script>fixture</script>';
  const output = research.renderStudy(input);
  assert.match(output, /data-study-id="test-passage"/);
  assert.match(output, /data-raw-value="6.000000000000001"/);
  assert.match(output, /&lt;script&gt;fixture&lt;\/script&gt;/);
  assert.doesNotMatch(output, /<script>/);
  assert.match(output, /synthetic observation/i);
  assert.match(output, /This fixture establishes no scientific result/);
  input.evidence_kind = 'native_hydrodynamic';
  assert.match(research.renderStudy(input), /not radiated acoustic pressure/i);
  assert.match(research.renderStudy(null), /unavailable/i);
});

test('research loading reports real ready/error states without a fallback dataset', async () => {
  assert.equal(typeof research.createResearchLoader, 'function', 'research loader is missing');
  const states = [];
  let broken = false;
  const loader = research.createResearchLoader({
    onState: state => states.push(state),
    fetcher: async url => ({ok: !broken, status: broken ? 404 : 200, url: String(url), json: async () => catalog()}),
  });
  await loader.load();
  assert.deepEqual(states.map(state => state.phase), ['loading', 'ready']);
  assert.equal(states.at(-1).catalog.studies[0].id, 'test-passage');
  broken = true;
  await loader.load();
  assert.deepEqual(states.slice(-2).map(state => state.phase), ['loading', 'error']);
  assert.equal(states.at(-1).catalog, undefined);
});

test('research loader rejects redirected identities', async () => {
  const states = [];
  await research.createResearchLoader({
    onState: state => states.push(state),
    fetcher: async () => ({ok: true, url: 'https://elsewhere.invalid/other.json', json: async () => catalog()}),
  }).load();
  assert.equal(states.at(-1).phase, 'error');
});

test('a superseded research request cannot replace a newer result', async () => {
  const states = [];
  let resolveFirst;
  let calls = 0;
  const loader = research.createResearchLoader({onState: state => states.push(state), fetcher: url => {
    const value = catalog(); value.studies[0].id = ++calls === 1 ? 'old' : 'new';
    const response = {ok: true, url: String(url), json: async () => value};
    return calls === 1 ? new Promise(resolve => {resolveFirst = () => resolve(response);}) : Promise.resolve(response);
  }});
  const first = loader.load();
  await loader.load();
  resolveFirst(); await first;
  assert.equal(states.at(-1).catalog.studies[0].id, 'new');
  assert.equal(states.filter(state => state.phase === 'ready').length, 1);
});

test('research requests have a real abort deadline', async () => {
  const states = [];
  await research.createResearchLoader({timeoutMs: 2, onState: state => states.push(state),
    fetcher: async (_url, {signal}) => new Promise((_resolve, reject) => {
      signal.addEventListener('abort', () => reject(new Error('aborted')));
    }),
  }).load();
  assert.equal(states.at(-1).phase, 'error');
  assert.match(states.at(-1).message, /timed out/i);
});

test('mounted research controls select exact studies without touching legacy results', async () => {
  assert.equal(typeof research.mountResearch, 'function', 'research mounting is missing');
  const elements = Object.fromEntries(['research-state', 'research-status', 'research-study', 'research-results', 'research-retry'].map(id => [id, {
    value: '', innerHTML: '', textContent: '', disabled: true, dataset: {}, events: {},
    setAttribute(name, value) {this[name] = value;},
    addEventListener(name, handler) {this.events[name] = handler;},
  }]));
  const data = catalog();
  data.studies.push({...study(), id: 'test-cfd', title: 'Native fixture', evidence_kind: 'native_hydrodynamic'});
  const mount = research.mountResearch({getElementById: id => elements[id]}, {
    fetcher: async url => ({ok: true, url: String(url), json: async () => data}),
  });
  await mount.initial;
  assert.equal(elements['research-state'].dataset.state, 'ready');
  assert.equal(elements['research-study'].disabled, false);
  elements['research-study'].value = 'test-cfd';
  elements['research-study'].events.change();
  assert.match(elements['research-results'].innerHTML, /data-study-id="test-cfd"/);
  assert.doesNotMatch(elements['research-results'].innerHTML, /data-study-id="test-passage"/);
});

test('existing page mounts the additive workflow with explicit evidence boundaries', () => {
  const html = readFileSync(new URL('../site/index.html', import.meta.url), 'utf8');
  assert.match(html, /src="\.\/research.mjs"/, 'research entry point is missing');
  for (const id of ['research', 'research-state', 'research-status', 'research-study', 'research-results', 'research-retry', 'passive-form', 'active-form']) {
    assert.match(html, new RegExp(`id="${id}"`));
  }
  assert.match(html, /sensor-coupling model/);
});

test('native field movies use declared assets and never autoplay', () => {
  const input = catalog();
  input.studies[0].evidence_kind = 'native_hydrodynamic';
  input.studies[0].video = {src: 'data/research/test-wake.mp4', poster: 'data/research/test-wake.png', title: 'Test-only movie', caption: 'Synthetic schema fixture, not CFD.'};
  research.validateResearch(input);
  const html = research.renderStudy(input.studies[0]);
  assert.match(html, /<video\b[^>]*controls/, 'native field video rendering is missing');
  assert.match(html, /data\/research\/test-wake.mp4/);
  assert.doesNotMatch(html, /\sautoplay(?:\s|=|>)/);
  input.studies[0].video.src = 'https://untrusted.invalid/test.mp4';
  assert.throws(() => research.validateResearch(input), /research|catalog/i);
});
