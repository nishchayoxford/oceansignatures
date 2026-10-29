// OceanSignatures: a read-only adapter for the versioned, computed catalogue.
// Imports are DOM-free; the browser entry point is guarded at the end of this file.
const KINDS = ['whistle', 'clicks', 'vessel'];
const fail = message => { throw new Error(`Invalid catalog: ${message}`); };
const object = (value, where) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) fail(`${where} must be an object`);
};
const text = (value, where, empty = false) => {
  if (typeof value !== 'string' || (!empty && !value.trim())) fail(`${where} must be text`);
};
const list = (value, where, check) => {
  if (!Array.isArray(value)) fail(`${where} must be a list`);
  value.forEach((item, index) => check(item, `${where}[${index}]`));
};
const number = (value, where, positive = false) => {
  if (typeof value !== 'number' || !Number.isFinite(value) || (positive && value <= 0)) fail(`${where} must be a finite ${positive ? 'positive ' : ''}number`);
};

export function safeAssetPath(path) {
  text(path, 'asset path');
  let decoded = path;
  for (let depth = 0; depth < 6; depth++) {
    if (decoded !== decoded.trim() || /[\\\u0000-\u001f\u007f:#?]/u.test(decoded) || decoded.startsWith('/') || decoded.split('/').includes('..')) fail('asset path must be site-relative without traversal');
    let next;
    try { next = decodeURIComponent(decoded); } catch { fail('asset path has invalid encoding'); }
    if (next === decoded) return path;
    decoded = next;
  }
  fail('asset path is excessively encoded');
}

function sourceURL(value) {
  text(value, 'source URL');
  let url;
  try { url = new URL(value); } catch { fail('source URL must be an absolute HTTP(S) URL'); }
  if (!['https:', 'http:'].includes(url.protocol) || url.username || url.password) fail('source URL must be HTTP(S) without credentials');
  return value;
}

function figureCheck(value, where) {
  object(value, where);
  safeAssetPath(value.src);
  for (const key of ['title', 'caption', 'alt']) text(value[key], `${where}.${key}`);
  if ('section' in value) text(value.section, `${where}.section`);
}
function downloadCheck(value, where) {
  object(value, where); safeAssetPath(value.href); text(value.label, `${where}.label`);
}
function resultCheck(value, where) {
  object(value, where);
  for (const key of ['title', 'summary']) text(value[key], `${where}.${key}`);
  list(value.figures, `${where}.figures`, figureCheck);
  list(value.downloads, `${where}.downloads`, downloadCheck);
  list(value.metrics, `${where}.metrics`, (metric, at) => {
    object(metric, at); text(metric.label, `${at}.label`); text(metric.unit, `${at}.unit`, true);
    if (typeof metric.value === 'number') number(metric.value, `${at}.value`);
    else text(metric.value, `${at}.value`);
  });
}

export function validateCatalog(catalog) {
  object(catalog, 'root');
  if (catalog.schema_version !== 1) fail('schema_version must be 1');
  for (const key of ['title', 'version', 'provenance_note']) text(catalog[key], key);
  const ids = new Set();
  for (const family of ['passive', 'active']) {
    const combinations = new Set();
    list(catalog[`${family}_cases`], `${family}_cases`, (item, where) => {
      resultCheck(item, where); text(item.id, `${where}.id`);
      if (ids.has(item.id)) fail(`duplicate case id ${item.id}`);
      ids.add(item.id); number(item.range_m, `${where}.range_m`, true);
      if (family === 'passive') {
        if (!KINDS.includes(item.kind)) fail(`${where}.kind is unsupported`);
        if (typeof item.surface !== 'boolean') fail(`${where}.surface must be boolean`);
      } else number(item.snr_db, `${where}.snr_db`);
      const key = family === 'passive' ? JSON.stringify([item.kind, item.range_m, item.surface]) : JSON.stringify([item.range_m, item.snr_db]);
      if (combinations.has(key)) fail(`ambiguous ${family} case selection`);
      combinations.add(key);
    });
  }
  if (catalog.reference !== null) {
    resultCheck(catalog.reference, 'reference');
    text(catalog.reference.calibration_status, 'reference.calibration_status');
    sourceURL(catalog.reference.source_url);
  }
  list(catalog.gallery, 'gallery', figureCheck);
  list(catalog.downloads, 'downloads', downloadCheck);
  list(catalog.sources, 'sources', (item, where) => {
    object(item, where); text(item.id, `${where}.id`); text(item.title, `${where}.title`); sourceURL(item.url);
  });
  object(catalog.verification, 'verification');
  text(catalog.verification.status, 'verification.status'); text(catalog.verification.summary, 'verification.summary');
  list(catalog.verification.checks, 'verification.checks', (item, where) => {
    object(item, where);
    for (const key of ['name', 'status', 'detail']) text(item[key], `${where}.${key}`);
  });
  return catalog;
}

export function selectCase(catalog, family, selection) {
  if (!['passive', 'active'].includes(family)) throw new Error('Unknown experiment family');
  if (!catalog) return null;
  return catalog[`${family}_cases`].find(item => item.range_m === selection.range_m && (family === 'passive'
    ? item.kind === selection.kind && item.surface === selection.surface
    : item.snr_db === selection.snr_db)) ?? null;
}

export function escapeHTML(value) {
  return String(value).replace(/[&<>"']/g, char => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[char]));
}
const esc = escapeHTML;

export function renderDownloads(downloads) {
  if (!downloads.length) return '<p class="quiet">No downloads are listed for this result.</p>';
  return `<ul class="download-list">${downloads.map(item => `<li><a href="${esc(safeAssetPath(item.href))}" download>${esc(item.label)} <span aria-hidden="true">↓</span></a></li>`).join('')}</ul>`;
}

export function renderMetrics(metrics) {
  if (!metrics.length) return '';
  return `<dl class="metrics">${metrics.map(item => {
    const display = typeof item.value === 'number' ? Number(item.value.toPrecision(8)).toString() : item.value;
    return `<div><dt>${esc(item.label)}</dt><dd data-raw-value="${esc(item.value)}">${esc(display)}${item.unit ? ` <span class="unit">${esc(item.unit)}</span>` : ''}</dd></div>`;
  }).join('')}</dl>`;
}

export function renderFigure(figure, {failed = false} = {}) {
  const src = esc(safeAssetPath(figure.src));
  return `<figure class="result-figure${failed ? ' figure-failed' : ''}" data-figure-src="${src}">
    <div class="figure-heading"><h4>${esc(figure.title)}</h4>${failed ? '' : `<a data-open-plot href="${src}" target="_blank" rel="noopener">Open full size <span aria-hidden="true">↗</span><span class="sr-only">: ${esc(figure.title)} (new tab)</span></a>`}</div>
    <div class="figure-media">${failed
      ? `<div class="image-error" role="status"><strong>Plot unavailable.</strong><p>The image could not be loaded. No substitute plot is shown.</p><button type="button" data-retry-figure>Retry this plot</button></div>`
      : `<img src="${src}" alt="${esc(figure.alt)}" loading="lazy" decoding="async">`}</div>
    <figcaption>${esc(figure.caption)}</figcaption>
  </figure>`;
}

export function renderCase(item, family) {
  if (!item) return '<div class="result-message" data-state="missing"><h3>No computed case matches this selection.</h3><p>No neighbouring case has been substituted. Choose another combination or reload the catalogue.</p></div>';
  const conditions = family === 'passive'
    ? `${item.kind} · ${item.range_m} m horizontal range · ${item.surface ? 'direct + surface reflection' : 'direct path only'}`
    : `${item.range_m} m point-target range · ${item.snr_db} dB input SNR`;
  return `<article class="computed-case" data-case-id="${esc(item.id)}" data-family="${esc(family)}">
    <header class="result-heading"><p class="eyebrow">Computed experiment / ${esc(conditions)}</p><h3>${esc(item.title)}</h3><p>${esc(item.summary)}</p><p class="case-identity">Case <code>${esc(item.id)}</code></p></header>
    ${renderMetrics(item.metrics)}
    <div class="figure-stack">${item.figures.length ? item.figures.map(figure => renderFigure(figure)).join('') : '<p class="result-message">No figures are listed for this case.</p>'}</div>
    <details class="case-downloads"><summary>Case data &amp; metadata</summary>${renderDownloads(item.downloads)}<p class="quiet">Use the exported metadata for exact sample rate, analysis window, signal normalization and numerical settings. The controls load saved results; they do not run a simulation.</p></details>
  </article>`;
}

export function renderReference(reference) {
  if (!reference) return '<div class="result-message" data-state="unavailable"><h3>No reference recording is available.</h3><p>The synthetic lessons remain useful, but no observed recording has been substituted or invented.</p></div>';
  return `<article class="reference-result"><h3>${esc(reference.title)}</h3><p>${esc(reference.summary)}</p>
    <p class="calibration"><strong>Calibration status:</strong> ${esc(reference.calibration_status)}</p>
    <p>Without a documented calibration, recording amplitudes are normalized digital units, not absolute pressure. Compare time–frequency structure, not absolute source level.</p>
    <a class="text-link" href="${esc(sourceURL(reference.source_url))}" target="_blank" rel="noopener noreferrer">Recording source &amp; provenance <span aria-hidden="true">↗</span><span class="sr-only"> (new tab)</span></a>
    ${renderMetrics(reference.metrics)}${reference.figures.map(figure => renderFigure(figure)).join('')}
    <details><summary>Reference data &amp; metadata</summary>${renderDownloads(reference.downloads)}</details></article>`;
}

export function renderVerification(verification) {
  return `<p class="verification-status"><strong>Reported status:</strong> ${esc(verification.status)}</p><p>${esc(verification.summary)}</p>
    ${verification.checks.length ? `<ul class="check-list">${verification.checks.map(check => `<li><div><strong>${esc(check.name)}</strong><span class="check-status">${esc(check.status)}</span></div><p>${esc(check.detail)}</p></li>`).join('')}</ul>` : '<p class="quiet">No individual verification checks are listed.</p>'}`;
}

export function renderSources(sources) {
  if (!sources.length) return '<p class="result-message">No source references are listed in this catalogue.</p>';
  return `<ol class="source-list">${sources.map(source => `<li><span class="source-id">${esc(source.id)}</span><a href="${esc(sourceURL(source.url))}" target="_blank" rel="noopener noreferrer">${esc(source.title)} <span aria-hidden="true">↗</span><span class="sr-only"> (new tab)</span></a><span class="source-url">${esc(source.url)}</span></li>`).join('')}</ol>`;
}

export function createCatalogLoader({fetcher = globalThis.fetch, onState, timeoutMs = 15000} = {}) {
  let generation = 0;
  let controller;
  const url = new URL('data/catalog.json', import.meta.url);
  return {
    async load() {
      const request = ++generation;
      controller?.abort();
      const ownController = new AbortController();
      controller = ownController;
      onState({phase: 'loading', catalog: null, error: null});
      let timer;
      try {
        const timeout = new Promise((_, reject) => {
          timer = setTimeout(() => {
            reject(new Error('Catalogue request timed out. Check your connection and retry.'));
            ownController.abort();
          }, timeoutMs);
        });
        const read = async () => {
          const response = await fetcher(url, {signal: ownController.signal, cache: 'no-store'});
          if (!response.ok) throw new Error(`Catalogue request failed (HTTP ${response.status}).`);
          return validateCatalog(await response.json());
        };
        const catalog = await Promise.race([read(), timeout]);
        if (request === generation) onState({phase: 'ready', catalog, error: null});
      } catch (error) {
        if (request === generation) onState({phase: 'error', catalog: null, error: error instanceof Error ? error.message : 'Catalogue could not be loaded.'});
      } finally { clearTimeout(timer); }
    },
    dispose() { ++generation; controller?.abort(); },
  };
}

export function renderCatalogView(state, selections) {
  const keys = ['passive', 'active', 'reference', 'gallery', 'verification', 'downloads', 'sources', 'provenance'];
  if (state.phase !== 'ready') {
    const loading = state.phase === 'loading';
    const message = loading ? 'Loading the computed catalogue…' : 'Computed results are unavailable. Use “Retry catalogue” above; no previous results are retained.';
    return {...Object.fromEntries(keys.map(key => [key, `<p class="result-message" data-state="${loading ? 'loading' : 'error'}">${message}</p>`])), disabled: true,
      status: loading ? 'Loading computed experiments…' : `<strong>Catalogue unavailable.</strong> ${esc(state.error || 'The catalogue could not be read.')} The lessons below remain readable.`,
    };
  }
  const c = state.catalog;
  return {
    disabled: false,
    status: `<strong>${esc(c.title)}</strong> · catalogue ${esc(c.version)} · precomputed results loaded`,
    passive: renderCase(selectCase(c, 'passive', selections.passive), 'passive'),
    active: renderCase(selectCase(c, 'active', selections.active), 'active'),
    reference: renderReference(c.reference),
    gallery: c.gallery.length ? c.gallery.map(figure => `${figure.section ? `<p class="eyebrow gallery-section">${esc(figure.section)}</p>` : ''}${renderFigure(figure)}`).join('') : '<p class="quiet">No supplementary plates are listed in this catalogue.</p>',
    verification: renderVerification(c.verification), downloads: renderDownloads(c.downloads),
    sources: renderSources(c.sources), provenance: `<p>${esc(c.provenance_note)}</p>`,
  };
}

export function mount(root = document) {
  const byId = id => root.getElementById(id);
  let currentCatalog = null;
  let retryNumber = 0;
  const selections = () => ({
    passive: {kind: byId('passive-kind').value, range_m: Number(byId('passive-range').value), surface: byId('passive-surface').value === 'true'},
    active: {range_m: Number(byId('active-range').value), snr_db: Number(byId('active-snr').value)},
  });
  const targets = {passive: 'passive-results', active: 'active-results', reference: 'reference-results', gallery: 'gallery-results', verification: 'verification-results', downloads: 'download-results', sources: 'source-results', provenance: 'provenance-results'};
  const announceCase = family => {
    const item = selectCase(currentCatalog, family, selections()[family]);
    byId(`${family}-announcement`).textContent = item ? `Showing ${item.title}. Case ${item.id}.` : 'No computed case matches this selection.';
  };
  const loader = createCatalogLoader({onState: state => {
    currentCatalog = state.phase === 'ready' ? state.catalog : null;
    const view = renderCatalogView(state, selections());
    for (const [key, id] of Object.entries(targets)) {
      byId(id).innerHTML = view[key];
      byId(id).setAttribute('aria-busy', String(state.phase === 'loading'));
    }
    byId('catalog-status').innerHTML = view.status;
    byId('catalog-state').dataset.state = state.phase;
    byId('catalog-retry').disabled = state.phase === 'loading';
    byId('catalog-retry').textContent = state.phase === 'error' ? 'Retry catalogue' : 'Reload catalogue';
    for (const family of ['passive', 'active']) {
      byId(`${family}-controls`).disabled = view.disabled;
      byId(`${family}-announcement`).textContent = state.phase === 'ready' ? '' : (state.phase === 'loading' ? 'Loading results.' : 'Results unavailable.');
    }
  }});
  for (const family of ['passive', 'active']) {
    byId(`${family}-form`).addEventListener('submit', event => event.preventDefault());
    byId(`${family}-controls`).addEventListener('change', () => {
      // Only replace this experiment, never the other selection or the reference.
      byId(`${family}-results`).innerHTML = renderCase(selectCase(currentCatalog, family, selections()[family]), family);
      announceCase(family);
    });
  }
  byId('catalog-retry').addEventListener('click', () => loader.load());

  const figureFromNode = node => ({
    src: node.dataset.figureSrc,
    title: node.querySelector('h4').textContent,
    caption: node.querySelector('figcaption').textContent,
    alt: node.querySelector('img')?.alt || node.dataset.originalAlt || 'Computed plot; see caption and full-size figure.',
  });
  root.addEventListener('error', event => {
    if (event.target.tagName !== 'IMG' || !event.target.isConnected) return;
    const node = event.target.closest('[data-figure-src]');
    if (!node) return;
    const figure = figureFromNode(node);
    const holder = root.createElement('div');
    holder.innerHTML = renderFigure(figure, {failed: true});
    const replacement = holder.firstElementChild;
    replacement.dataset.originalAlt = figure.alt;
    node.replaceWith(replacement);
  }, true);
  root.addEventListener('click', event => {
    const retry = event.target.closest?.('[data-retry-figure]');
    if (!retry) return;
    const node = retry.closest('[data-figure-src]');
    const figure = figureFromNode(node);
    const holder = root.createElement('div');
    holder.innerHTML = renderFigure(figure);
    const replacement = holder.firstElementChild;
    // Retry only this validated local asset; defeat a cached image failure.
    replacement.querySelector('img').src = `${safeAssetPath(figure.src)}?retry=${++retryNumber}`;
    node.replaceWith(replacement);
    replacement.querySelector('[data-open-plot]').focus({preventScroll: true});
  });
  loader.load();
  return loader;
}

if (typeof document !== 'undefined' && document.getElementById('catalog-state')) mount(document);
