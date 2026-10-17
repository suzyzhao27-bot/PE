# Reproducibility Repair and B10 Closeout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore byte-exact provenance checks, remove the dielectric figure annotation overlap, and close the locally verifiable portion of P1/B10 without changing scientific equations or reported model values.

**Architecture:** Protect the one registered byte-copy with a path-specific Git attribute, make figure layout a small testable function, then rerun the existing scientific test and calculation gates. Research conclusions are written only to the existing formal board, reproducibility record, and P1 stage summary.

**Tech Stack:** Git attributes, Python 3.14, unittest, NumPy, SciPy, Matplotlib.

---

### Task 1: Protect the registered PB source bytes

**Files:**
- Create: `.gitattributes`
- Verify: `src/long_tube_diagnostic/test_diagnostic.py`
- Restore bytes: `src/long_tube_diagnostic/pb_v1_readonly_copy.py`

- [x] **Step 1: Run the existing hash identity test and record the expected failure**

Run from `src/long_tube_diagnostic/`: `py -3.14 -m unittest test_diagnostic.py -v`
Expected: exactly `test_registered_copy_hash_identity` fails because the working copy SHA256 differs from `3c681f6408af6e133586b1a6f8b2b0f474ce4aa014ce8c1fb7db070406c96cd2`.

- [x] **Step 2: Add the narrow Git attribute**

Add exactly:

```gitattributes
src/long_tube_diagnostic/pb_v1_readonly_copy.py text eol=lf
```

- [x] **Step 3: Restore the tracked file from the Git index under the new attribute**

Temporarily stage `.gitattributes`, run `git checkout-index --force -- src/long_tube_diagnostic/pb_v1_readonly_copy.py`, verify its SHA256, then immediately unstage `.gitattributes`. The explicit `text eol=lf` attribute is required because both `checkout-index` and `git archive` otherwise apply this checkout's CRLF conversion.
Expected: `git ls-files --eol` reports `i/lf w/lf attr/text eol=lf` and the SHA256 equals the registered value.

- [x] **Step 4: Run the long-tube tests**

Run from `src/long_tube_diagnostic/`: `py -3.14 -m unittest test_diagnostic.py -v`
Expected: 7 tests pass.

### Task 2: Make dielectric plot layout testable and non-overlapping

**Files:**
- Modify: `src/surface_pb/plot_surface.py`
- Create: `src/surface_pb/test_plot_surface.py`

- [x] **Step 1: Add a failing layout test**

The test imports `build_dielectric_control_figure`, draws the canvas, obtains the display-space bounding boxes for the x-axis label and note, and asserts a positive vertical gap:

```python
figure, axis, note = build_dielectric_control_figure()
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
label_box = axis.xaxis.label.get_window_extent(renderer)
note_box = note.get_window_extent(renderer)
self.assertGreater(label_box.y0 - note_box.y1, 4.0)
```

- [x] **Step 2: Verify the new test fails before implementation**

Run from `src/surface_pb/`: `py -3.14 -m unittest test_plot_surface.py -v`
Expected: import failure because `build_dielectric_control_figure` does not yet exist.

- [x] **Step 3: Extract the minimal figure builder and reserve bottom space**

Implement `build_dielectric_control_figure()` in `plot_surface.py`. It computes the same `z`, `dg`, `lb`, and `eta`, creates the same curve and labels, sets an explicit bottom margin, places the note below the axis label, and returns `(figure, axis, note, data)` without saving. Keep the analytical formula and numerical inputs unchanged.

- [x] **Step 4: Verify the focused plot test and surface physics suite**

Run from `src/surface_pb/`: `py -3.14 -m unittest test_plot_surface.py test_surface.py -v`
Expected: 10 tests pass.

- [x] **Step 5: Regenerate and inspect the two surface figures**

Run from `src/surface_pb/`: `py -3.14 run_surface.py` then `py -3.14 plot_surface.py`.
Expected: two PNG files and the analytical control CSV are regenerated; `dielectric_symmetric_control.png` has no annotation/x-label overlap.

### Task 3: Re-run all scientific gates and the diagnostic entry point

**Files:**
- Generated ignored outputs only under `src/**/outputs/` and `src/qa/`

- [x] **Step 1: Run all four test suites plus the new layout test**

Run the v2, long-tube, surface, finite-PNP, and layout suites with Python 3.14.
Expected: 43 tests pass, with no failures or errors (including the 6 v2 pipeline tests).

- [x] **Step 2: Run the long-tube diagnostic**

Run from `src/long_tube_diagnostic/`: `py -3.14 run_diagnostic.py`.
Expected: completion without the source-hash exception and a manifest recording `copy_hash_matches_registered=true`.

- [x] **Step 3: Confirm primary numerical invariants remain unchanged**

Compare current 32-grid five-row values against the established baseline: flux ratio `1.0275731111`, `t_Na=0.6528317904`, static mean sigma `-0.0047360405 C/m2`, open-circuit mean sigma `-0.0077769197 C/m2`, and liquid potential difference `-10.43103665 mV`.

### Task 4: Close the locally verifiable part of B10

**Files:**
- Modify: `项目进展计划与任务看板.md`
- Modify: `复现记录_20261015.md`
- Create: `阶段总结/20261015_P1_B10_zeta源数据审计.md`

- [x] **Step 1: Record source-identity and uncertainty audit results**

Document 36 observations, 9 groups, 3 source hashes, four unique ramps per group, exact reproduction of means and both SD conventions, absence of blank/duplicate source identities, and independent membrane count `unknown`.

- [x] **Step 2: State the IEP limitation**

Record that 4.1381, 4.0532, and 3.9964 are two-point linear interpolations between approximately pH 3 and 7, not direct dense-scan measurements.

- [x] **Step 3: Update card status without overstating completion**

Mark B10 complete for the available files; keep B11 and B12 blocked by Origin exports and instrument method/cell/background-electrolyte records. Keep P1 overall `进行中`.

- [x] **Step 4: Record rerun environment and verification evidence**

Add Python, NumPy, SciPy, and Matplotlib versions, test counts, hashes, regenerated output locations, and unchanged scientific interpretation to the reproducibility record.

### Task 5: Final consistency check

**Files:**
- Review all modified tracked files and `.planning` records

- [x] **Step 1: Review the diff for unrelated changes**

Run: `git diff --check` and `git status --short`.
Expected: only the approved Git attribute, plot/test repair, existing formal research records, P1 summary, and the pre-existing `.gitignore` change appear.

- [x] **Step 2: Verify the project folder remains organized**

Confirm root files are formal project entry points or reproducibility controls; implementation/process notes remain under `.planning`; generated data remain under module output directories.
