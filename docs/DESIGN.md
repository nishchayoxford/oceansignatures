---
version: alpha
name: OceanSignatures
description: An editorial field guide to underwater sound.
colors:
  primary: "#172E33"
  secondary: "#52676A"
  tertiary: "#1C665E"
  neutral: "#F5F4EE"
  plate: "#FFFFFF"
  deep: "#0E252B"
  sea: "#85D6C7"
  on-deep: "#F4F5ED"
  dim-deep: "#B6C8C4"
  line: "#C6D0C9"
  accent-soft: "#E4EEE7"
  danger: "#963D29"
typography:
  h1:
    fontFamily: Georgia
    fontSize: 6.25rem
    fontWeight: 400
    lineHeight: 0.99
    letterSpacing: "-0.035em"
  h2:
    fontFamily: Georgia
    fontSize: 3.6rem
    fontWeight: 400
    lineHeight: 1.09
    letterSpacing: "-0.035em"
  body-md:
    fontFamily: Segoe UI
    fontSize: 17px
    fontWeight: 400
    lineHeight: 1.65
  annotation:
    fontFamily: Consolas
    fontSize: 0.7rem
    lineHeight: 1.7
    letterSpacing: 0.08em
rounded:
  sm: 2px
spacing:
  sm: 8px
  md: 16px
  lg: 24px
  xl: 40px
  section: 96px
components:
  body:
    backgroundColor: "{colors.neutral}"
    textColor: "{colors.primary}"
  caption:
    backgroundColor: "{colors.plate}"
    textColor: "{colors.secondary}"
  button-primary:
    backgroundColor: "{colors.sea}"
    textColor: "{colors.deep}"
    rounded: "{rounded.sm}"
    padding: 12px
  button-secondary:
    backgroundColor: "{colors.neutral}"
    textColor: "{colors.tertiary}"
    rounded: "{rounded.sm}"
    padding: 12px
  hero:
    backgroundColor: "{colors.deep}"
    textColor: "{colors.on-deep}"
  hero-note:
    backgroundColor: "{colors.deep}"
    textColor: "{colors.dim-deep}"
  control-panel:
    backgroundColor: "{colors.accent-soft}"
    textColor: "{colors.primary}"
  rule:
    backgroundColor: "{colors.line}"
  error:
    backgroundColor: "{colors.neutral}"
    textColor: "{colors.danger}"
---

## Overview

**Primary surface: Decide / Learn. Secondary surface: Explore.**

An original editorial ocean-lab field guide, not a generic SaaS dashboard. A restrained midnight opening asks one physical question; off-white technical plates then explain the experiment before presenting its controls. A conceptual source/channel/receiver schematic is explicitly labelled as a schematic, never a computed field. No fake metrics, testimonial, institution, company, interview or field-performance claim fills missing data.

The learning sequence is source → channel → receiver; passive listening → spectral interpretation → active ranging → observations → evidence. Beginners get a concrete question beside each result. An informed scientific reader can open the underlying equations, source conventions, SNR caveats and reproducible metadata without losing the selected experiment.

## Colors

- Midnight `#0E252B` anchors the masthead, introductory plate and footer.
- Paper `#F5F4EE` supports sustained reading; white plots preserve the producer's original figure colours.
- Ink `#172E33` is primary text. Muted text `#52676A` is used for captions, never low-contrast decorative greys.
- Teal `#1C665E` denotes links and controls; light sea `#85D6C7` is the readable accent on midnight.
- Quiet green technical inserts distinguish spectral interpretation and evidence from the experiment sequence.
- Failure uses explicit words and a warm text colour, not colour alone. Verification status is copied literally from the catalogue, not converted into an invented success badge.

Colour contrast and responsive-layout checks are separate from scientific verification.

## Typography

Georgia is the intentionally editorial headline family. The body uses locally available Aptos / Segoe UI / Helvetica Neue / Arial; annotations use local Consolas / Liberation Mono. There are no downloaded fonts, CDN requests, frameworks, analytics or external runtime dependencies.

Body text is 17 px at desktop and 16 px at mobile with generous line spacing. Plot captions are secondary prose, not tiny metadata. Native MathML provides fractions, superscripts, subscripts, sums and square roots. Each equation scrolls locally and has a keyboard-focusable named region. Delimiters are scoped to their operands so unrelated expressions do not stretch parentheses.

## Layout

- Reading width: a 1404 px maximum container, with 72 px desktop gutters.
- Hero: asymmetric narrative and conceptual geometry, not a centered feature-card composition.
- Experiments: a narrow lesson column beside a wider, uninterrupted figure column. All catalogue figures retain natural aspect ratio and offer a full-size link.
- At 820 px the layout becomes a single reading column. At 560 px all experiment selectors stack. There is no horizontal page scrolling intended; equations alone can scroll.
- Quiet horizontal rules organize methods, metrics and sources. Observed values are not presented as oversized achievements.
- The masthead provides Learn / Listen / Locate / Evidence. A visible keyboard skip link bypasses it.
- Reduced motion is explicitly supported; the site has no continuous animation or loading spinner.

## Elevation & Depth

No glass, blur, gradient or shadow system. A plot is a flat technical plate with a small heading and an explanatory caption. Hierarchy comes from typography, whitespace and rules.

## Shapes

Controls use restrained 2 px corners. Native selects, details and semantic links remain identifiable. Controls are at least 44 px high; selects are 48 px. No icon-topped tile grids or floating decorative badges.

## Components

### Static shell

`site/index.html` remains an intelligible scientific guide when JavaScript or data loading fails. It contains the physical assumptions, equations and no-JavaScript message. The inline SVG depicts geometry only and is labelled not to scale. The scientific source register, reference provenance, numerical verification status, result metrics and plots all come from the catalogue.

### Read-only data adapter

`site/app.mjs` implements DOM-free validation, exact selection and render helpers. Its `validateCatalog` function defines schema version 1; child fields are not inferred from undocumented producer arrays. Relative assets resolve correctly under `/oceansignatures/`, not the domain root.

All expected control tuples are discrete:

- Passive: whistle / clicks / vessel × 25 / 100 / 400 m × direct / surface.
- Active: 75 / 150 / 300 m × −18 / −6 dB.

The default passive selection is whistle, 100 m, direct only. The default active selection is 150 m, −18 dB. The adapter does not invent a nearest case when a requested tuple is absent. Uniqueness checks reject ambiguous case IDs and tuples. Catalogue content is HTML-escaped; local asset paths reject remote URLs, root-relative paths, traversal, unsafe schemes and encoded traversal. Public reference URLs allow HTTP(S) only.

### Failure and recovery

Loading and failed catalogue requests remove every data-owned panel's previous results immediately. Controls disable during unavailable states; the static explanation remains. Network, HTTP, JSON, schema and timeout failures are visible, with a retry control. A request-generation guard prevents late responses overwriting newer catalogue state.

An image failure removes the image and its open-plot link, retains the current caption and shows a per-figure retry. A controlled retry query bypasses cached image failures without changing catalogue metadata. Case changes replace the result subtree; detached image events cannot alter a new selection. Nothing generates a substitute plot.

### Scientific scope

Source recipes, emission-window RMS, receiver group delay, SNR support, correlation normalization and recording-pipeline details are linked to the implementation documentation. A fixed-SNR range comparison explicitly does not represent fixed ambient noise. A matched-filter estimate is not a detection probability. Reference recordings remain digital-unit observations unless documented calibration supports a pressure conversion.

## Do's and Don'ts

- Do preserve the relationship between case ID, chosen parameters, caption, metric units and downloads.
- Do keep source normalization, spectral density and playback normalization distinct.
- Do copy verification status literally and separate it from physical validation.
- Do place important plot interpretation beside the figures and exact assumptions in disclosure sections.
- Do not add a frontend solver, fake production data, synthesized statistics or fabricated verification totals.
- Do not normalize plotted arrays in the browser, relabel an old image after a selection change, or hide a failed result behind a plausible placeholder.
- Do not make analytical recipes sound like validated biological or machinery models.

## Checks and retained evidence

Command: `node --test tests/test_site.mjs`.

Fixtures are constructed in memory inside the test module. No test catalogue or plot is written into production data paths. These interface tests do not substitute for exercising the real computed catalogue in a browser.

The test suite covers schema integrity, typed exact selection across all planned tuples, escaping and safe paths, metrics/captions/downloads, unavailable cases, reference absence and calibration copy, stale-state clearing, HTTP/schema recovery, request races, timeout, image-error markup, static accessibility hooks and key scientific caveats. Node syntax is checked separately.

Real-data figure inspection and browser acceptance are separate checks. A structural/unit pass is not visual verification, scientific validation or download-integrity evidence.

## Proposed browser acceptance selectors

1. Serve the real `site/` at a project prefix. Wait for `#catalog-state[data-state="ready"]`; verify there are no console exceptions and `#catalog-retry` is enabled.
2. Select all exact passive combinations with `#passive-kind`, `#passive-range`, `#passive-surface`. Assert `#passive-results [data-case-id]` equals the catalogue lookup and its `img[src]`, metric units and download `href` values match that exact case. Active equivalents: `#active-range`, `#active-snr`, `#active-results [data-case-id]`.
3. Use real keyboard focus and arrow keys on native selects; assert a changed case ID, not merely a focused widget. Confirm passive changes do not rebuild `#active-results` or `#reference-results`.
4. Await each visible image's `complete && naturalWidth > 0` before screenshotting. Inspect `.result-figure`, `.figure-heading`, `.metrics`, `figcaption` at 1440 px, 768 px, 390 px and 320 px. Use `[data-open-plot]` to inspect original scientific axes.
5. Expand `#passive-theory > summary`, `#active-theory > summary`, `#verification-details > summary` and `#downloads-details > summary` by keyboard. Inspect `.equation math` fractions and local scrolling; no page-level horizontal overflow.
6. Inject HTTP failure / invalid schema / missing case at the in-memory fetch boundary only. Click `#catalog-retry`; expect `#catalog-state[data-state="error"]`, disabled `#passive-controls` / `#active-controls`, and zero stale `img` / `[data-case-id]` inside data-owned panels. Restore the real response and retry.
7. Fail a selected image request. Assert `[data-figure-src].figure-failed` contains no `img` or `[data-open-plot]`, has its current caption and `[data-retry-figure]`. Restore the asset, click retry and await a decoded replacement image. Repeat across a rapid case switch to catch stale events.
8. Verify `#reference-results` against `reference`, including the null branch; source URLs in `#source-results`; provenance text in `#provenance-results`; literal check statuses in `#verification-results`.
9. Exercise `.case-downloads`, `#download-results a[download]` and `[data-open-plot]` with actual returned bytes under `/oceansignatures/`. Do not infer successful downloads from link existence.
10. Check reduced-motion preference, focus visibility, 44 px controls, readable long captions and no-JavaScript fallback. The interface should not need any third-party network request.
