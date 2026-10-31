// Additive read-only research results; the original experiment adapter is unchanged.
import {safeAssetPath, escapeHTML as esc, renderMetrics, renderFigure, renderDownloads} from './app.mjs';

const fail = message => {throw new Error(`Invalid research catalog: ${message}`);};
const text = value => {
  if (typeof value !== 'string' || !value.trim()) fail('nonempty text required');
};
const list = (value, check) => {
  if (!Array.isArray(value) || !value.length) fail('nonempty list required');
  value.forEach(check);
};

export function validateResearch(catalog) {
  if (!catalog || catalog.schema_version !== 1 || catalog.project_version !== '0.1.0') fail('project/schema identity');
  const ids = new Set();
  list(catalog.studies, item => {
    if (!item || !/^[a-z0-9-]+$/.test(item.id) || typeof item.id !== 'string' || ids.has(item.id)) fail('study identity');
    ids.add(item.id);
    text(item.title); text(item.summary);
    if (!['predicted_acoustic_with_synthetic_observation', 'native_hydrodynamic'].includes(item.evidence_kind)) fail('evidence kind');
    list(item.metrics, metric => {
      text(metric.label);
      if (typeof metric.unit !== 'string' || typeof metric.value !== 'number' || !Number.isFinite(metric.value)) fail('metric value/unit');
    });
    list(item.figures, figure => {
      safeAssetPath(figure.src);
      text(figure.title); text(figure.caption); text(figure.alt);
    });
    list(item.downloads, download => {safeAssetPath(download.href); text(download.label);});
    list(item.limitations, text);
    if ('video' in item) {
      if (!item.video || item.evidence_kind !== 'native_hydrodynamic') fail('native video definition');
      safeAssetPath(item.video.src); safeAssetPath(item.video.poster);
      text(item.video.title); text(item.video.caption);
    }
  });
  return catalog;
}

export function renderStudy(study) {
  if (!study) return '<p class="result-message">Research result unavailable. No substitute result is shown.</p>';
  const evidence = study.evidence_kind === 'native_hydrodynamic'
    ? 'Native OpenFOAM / local hydrodynamic pressure — not radiated acoustic pressure'
    : 'Empirical acoustic prediction / controlled synthetic observation — not field validation';
  return `<article class="computed-research" data-study-id="${esc(study.id)}">
    <header class="result-heading"><p class="eyebrow">${esc(evidence)}</p><h3>${esc(study.title)}</h3><p>${esc(study.summary)}</p></header>
    ${renderMetrics(study.metrics)}
    <div class="research-boundary"><h4>Read these results within their scope</h4><ul>${study.limitations.map(item => `<li>${esc(item)}</li>`).join('')}</ul></div>
    ${study.video ? `<figure class="result-figure research-movie"><div class="figure-heading"><h4>${esc(study.video.title)}</h4><a href="${esc(safeAssetPath(study.video.src))}" download>Download native-field movie ↓</a></div><video controls playsinline preload="metadata" poster="${esc(safeAssetPath(study.video.poster))}" src="${esc(safeAssetPath(study.video.src))}" aria-label="${esc(study.video.title)}">Your browser cannot display this video. Use the movie download.</video><figcaption>${esc(study.video.caption)}</figcaption></figure>` : ''}
    ${study.figures.map(figure => renderFigure(figure)).join('')}
    <details class="research-downloads"><summary>Numerical data, native evidence &amp; reproduction</summary>${renderDownloads(study.downloads)}</details>
  </article>`;
}

export function createResearchLoader({fetcher = globalThis.fetch, onState, timeoutMs = 15000} = {}) {
  const url = new URL('./data/research/catalog.json', import.meta.url);
  let generation = 0;
  let pending;
  return {
    async load() {
      const request = ++generation;
      pending?.abort();
      const controller = new AbortController();
      pending = controller;
      const timer = setTimeout(() => controller.abort(), timeoutMs);
      onState({phase: 'loading'});
      try {
        const response = await fetcher(url, {cache: 'no-cache', signal: controller.signal});
        if (!response.ok || response.url !== url.href) {
          await response.body?.cancel();
          throw new Error(!response.ok ? `HTTP ${response.status}` : 'Unexpected research response identity');
        }
        const catalog = validateResearch(await response.json());
        if (request === generation) onState({phase: 'ready', catalog});
      } catch (error) {
        if (request === generation) onState({phase: 'error', message: controller.signal.aborted ? 'Research request timed out.' : String(error.message)});
      } finally {
        clearTimeout(timer);
      }
    },
  };
}

export function mountResearch(root = document, {fetcher = globalThis.fetch} = {}) {
  const get = id => root.getElementById(id);
  const select = get('research-study');
  const results = get('research-results');
  const status = get('research-status');
  const retry = get('research-retry');
  let catalog = null;
  const show = () => {
    const study = catalog?.studies.find(item => item.id === select.value);
    results.innerHTML = renderStudy(study);
    status.textContent = study ? `Showing ${study.title}. Stored results; no live simulation.` : 'No research result matches this selection.';
  };
  const loader = createResearchLoader({fetcher, onState: state => {
    const previous = select.value;
    catalog = state.phase === 'ready' ? state.catalog : null;
    get('research-state').dataset.state = state.phase;
    results.setAttribute('aria-busy', String(state.phase === 'loading'));
    select.disabled = state.phase !== 'ready';
    retry.disabled = state.phase === 'loading';
    if (catalog) {
      select.innerHTML = catalog.studies.map(item => `<option value="${esc(item.id)}">${esc(item.title)}</option>`).join('');
      select.value = catalog.studies.some(item => item.id === previous) ? previous : catalog.studies[0].id;
      show();
    } else {
      select.innerHTML = '<option value="">Results not loaded</option>';
      results.innerHTML = renderStudy(null);
      status.textContent = state.phase === 'loading' ? 'Loading computed research…' : `Research unavailable: ${state.message}. The original experiments remain independent.`;
    }
  }});
  select.addEventListener('change', show);
  retry.addEventListener('click', () => loader.load());
  return {...loader, initial: loader.load()};
}

if (typeof document !== 'undefined' && document.getElementById('research-state')) mountResearch(document);
