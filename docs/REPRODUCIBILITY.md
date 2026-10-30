# Reproduce OceanSignatures

The project remains **0.1.0**. These routes separate rebuilding synthetic
acoustics, reanalysing saved native CFD and launching a new native simulation.
Only the first two are needed to inspect or reproduce the demonstration.

Run commands from the repository root unless stated otherwise. Use new output
directories: differing or incomplete existing research products are deliberately
not overwritten. Keep every part of the native archive beside its `package.json`.

## 1. Environment and software tests

```sh
uv sync --frozen
uv run pytest -q
node --test tests/test_site.mjs tests/test_research_site.mjs
uv run ruff check src tests scripts
uv run ruff format --check src tests scripts
```

The default interpreter is Python 3.11; numerical dependencies are in `uv.lock`.
The CI frontend runtime is Node.js 24. FFmpeg is needed only to regenerate the
native movie; use `--no-movie` on the CFD builder when it is unavailable, and do
not expect a movie-less export to match the complete export manifest.

Pytest does **not** run OpenFOAM. Runner tests use explicitly test-only surrogate
cases, recording stubs and small Python subprocesses; those are software tests,
not native flow results. Real flow evidence is the separately retained OpenFOAM
case data and its externally pinned completion receipts.

## 2. View the complete existing site

```sh
python -m http.server 8000 --bind 127.0.0.1 --directory site
```

Open `http://127.0.0.1:8000/`. No online solver, CDN or hardware is required.
The Research selector loads the real acoustic and native-flow bundles already
under `site/data/research/`. Original passive/active controls and downloads remain
separate and unchanged.

## 3. Recompute the original teaching experiments

```sh
uv run python scripts/build_site.py --output build/teaching-replay
```

This produces the original generated data/plots. It is **not** the research
builder and does not calculate a vessel passage or execute CFD. The website
shell and its research assets are not replaced by this command.

## 4. Recompute the acoustic research experiment

```sh
uv run python scripts/build_research.py \
  --output build/acoustic-replay --seed 1729 --gain-db 6
```

This actually evaluates the empirical spectra, synthesizes the declared WAV
chunks, calibrates and processes them, computes comparisons and sensitivity,
and renders the six figures. It does not download weather, AIS or recordings.
The chosen UTC epoch is synthetic. See [RESEARCH.md](RESEARCH.md) for the full
scenario, estimator, missing-data policy and exported array definitions.

Repeat the exact command to exercise protected deterministic replay. An
identical rerun is accepted; a different seed, implementation or numerical-stack
provenance needs another destination. The manifest pins every payload except
itself, and the CLI prints the manifest digest. Source-code identities and
software versions are retained in `metadata.json`.

Check numerical results and provenance separately. A changed source hash is a
real provenance difference even when WAV bytes and plotted values are unchanged;
do not replace a historical manifest with newly calculated hashes to make a
comparison pass. Bitwise rendering identity is conditional on the software/font
stack, not promised across arbitrary operating systems.

## 5. Replay the saved native CFD without OpenFOAM

```sh
uv run python scripts/build_cfd_results.py \
  --baseline site/data/research/cfd/native/baseline \
  --sensitivity site/data/research/cfd/native/sensitivity \
  --output build/cfd-replay
```

The public reader consumes the compressed native subset directly, verifies its
pins, calculates pressure/force statistics, and renders the figures and movie.
It does **not** invoke a mesh generator, solver or reconstruction utility.
No original private case directory is required. The packaged `replay/` source
also supplies a standalone route described in [CFD_RESULTS.md](CFD_RESULTS.md).

This replay checks 45 retained original pins per case and the anchored original
completion receipt. The original full-case audit checked 206 pins per case.
Omitted full-volume and unused plane files are enumerated, not silently treated
as verified by the subset. Consequently `verification.json` and its manifest
legitimately differ between a full-native-input export and a subset replay.
Compare the same verification scope for byte-identical replay, and numerical
payloads separately when comparing those two scopes.

A fresh destination is required for each CFD replay. Repeating a CFD export onto
an existing destination is refused, unlike the verified-identical acoustic
rerun. Native data, logs, dictionaries and completion identities remain unchanged.

## 6. Reassemble the authenticated research site bundle

This route copies the existing scientific bundles into a new directory
and recreates their catalogue. It is an **assembly check**, not a scientific
recalculation. Obtain expected hashes from the retained catalogue, not by
silently trusting whatever payloads happen to be present:

```sh
ACOUSTIC_SHA=$(uv run python -c 'import json; print(json.load(open("site/data/research/catalog.json"))["source_manifests"]["acoustic"])')
CFD_SHA=$(uv run python -c 'import json; print(json.load(open("site/data/research/catalog.json"))["source_manifests"]["cfd"])')
uv run python scripts/build_research_site.py \
  --acoustic site/data/research/acoustic \
  --cfd site/data/research/cfd \
  --acoustic-sha256 "$ACOUSTIC_SHA" \
  --cfd-sha256 "$CFD_SHA" \
  --output build/research-assembly
```

The assembler checks each expected manifest, complete membership, payload hashes,
scientific evidence kinds and safe paths before staging the output. Its 85-file
family contains both complete scientific bundles, the catalogue and a master
manifest. Exact reruns are accepted; differing existing output is refused.

A newly computed acoustic or CFD bundle is a **new candidate**, not an automatic
replacement for the served assets. Validate it and record its accepted manifest
identity before assembly or promotion. Keep earlier accepted and failed evidence;
do not overwrite source data or the previous results to conceal a difference.

## 7. New native computation is a separate operation

The genuine retained cases use OpenFOAM **v2206**, `pimpleFoam`, a two-dimensional
laminar Re100 cylinder and matched 0.02 s probe/force clocks. The setup,
boundaries, numerical choices and opt-in commands are in [CFD.md](CFD.md).

Generating dictionaries is not executing a solver. A native run requires the
explicit `--execute` option, the correct local OpenFOAM environment, a fresh
case and the campaign's runtime/storage/locking gates. Do not use new execution
to reproduce figures that can already be replayed from saved evidence. A new
run must earn its own completion and provenance record; it cannot inherit the
identity or physical conclusions of the existing pair.

## What the checks establish

- Unit and independent-oracle tests check numerical definitions and supported
  boundary behavior.
- Synthetic gain recovery checks calibration, spectral processing and exact
  comparison alignment, not environmental source-model accuracy.
- Saved-native replay checks traceable data reduction and deterministic products,
  not a second independent physical realization.
- The two native timesteps provide sensitivity evidence, not spatial convergence,
  an extrapolated convergence order or a validated instrument self-noise model.
- Website checks establish that the inspected controls and downloads expose the
  intended assets. They do not upgrade synthetic observations into measurements.

No authorized calibrated field recording is bundled. Optional external recording
acquisition is outside these commands and remains subject to rights, calibration
and deployment-metadata checks.
