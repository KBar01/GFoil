# Production-readiness cleanup — report (readability only, zero behaviour change)

## Verification (the hard success criterion)

- **Regression golden 10/10 bit-identical** (CL/CD/CM/OASPL **and** all three AD
  gradient arrays, rel_err 0.000e+00) after **every** commit.
- **Full bench cold-start sweep** (1381 cases) re-run after the changes and
  diffed against the pre-cleanup `post-A` reference: **identical converged set
  (1208 converged), bit-identical CL/CD/CM/OASPL, 0 diffs.**
- Every commit's diff was verified to contain **only comments and the one token
  rename** — zero re-expressed arithmetic, zero non-comment line removals (except
  the intentional dead-comment deletions).

## Per-area summary

### Stage 1 — comments & documentation
- **`NOMENCLATURE.md`** (new): one glossary for the terse domain-standard symbols
  — the `[th, ds, sa, ue]` state vector, the `_U`/`_x` derivative convention,
  Hk/Hs/Hss/Ret/cf/cDi/cteq/cttr/damp/Us/uq/upw, compressibility terms — kept
  rather than renamed (they match Drela/XFOIL and Fidkowski 2021), with refs.
- **`get_funcs.hpp`** (was 2.7 % comments, the densest file): full file header
  (pipeline role, the `_U` convention, CoDi-tape warning, glossary pointer)
  replacing the stale optimisation-TODO; 1–4 line docstrings on every closure
  (`get_cp/uk/Mach2/H/Hw/Hk/Hss/de/Ret/cf/cfxt/Hs/Us/uq/cDi*+total/cDixt/upw/cteq/damp/cttr`).
- **File-header blocks** added to: `coupled.cpp`, `update_state.cpp`,
  `update_transition.cpp`, `build_global_sys.cpp`, `solve_glob.cpp`,
  `init_BL.cpp`, `residuals_shared.hpp`, `extract_BL_TE.hpp` — each stating the
  function in the solve pipeline and the CoDi-tape caveat.
- **Comment typos fixed**: "sqaured"→"squared", "parametrs"→"parameters",
  "calcuolating"→removed, "fofr"→"for".
- **Corrected a wrong comment**: the "Amplification rate dn/dxi" label sat above
  `get_de`, which actually returns the BL thickness δ; relabelled, and the
  amplification rate identified as `get_damp`.
- **Cruft removed**: the stale optimisation-musing blocks (get_funcs.hpp top,
  init_BL.cpp) folded into concise perf `NOTE:`s; 3 commented-out
  `//check_for_nans(...)` debug lines deleted.

### Stage 2 — internal naming
- Renamed the opaque `aux1`/`aux2` → **`wgap1`/`wgap2`** across all 5 files that
  use them (`residuals_shared.hpp`, `residuals.h`, `build_global_sys.cpp`,
  `init_BL.cpp`, `srcAD/.../main_func.hpp`). They are the trailing-edge/wake gap
  subtracted from δ* in `residual_station` (`ds = U[1] - wgap`); the new name
  states the quantity. Pure token substitution; golden + full-sweep unchanged.
- Domain-standard shorthand (θ/δ*/Hk/…) deliberately **left alone** (glossary
  covers them) per the convention that matches the literature.

### Stage 3 — intra-file structure
- Navigable section banners added to `get_funcs.hpp`: Compressibility &
  thermodynamics · Shape factors & integral thicknesses · Skin friction · Energy
  shape factor & slip · Dissipation coefficient · Discretisation upwinding ·
  Shear-stress closure & e^N transition · Inviscid edge velocity. (Pure comments;
  no function bodies moved — moving statements inside active CoDi functions is
  unsafe and was not done.)

### Stage 4 — cross-file organisation
- **Proposed, NOT executed** (await review): none rose to the bar of "clearly
  mis-housed". The sparse-assembly helpers (`equate_block_inplace_sparse`,
  `findColumnIndices`, `addColumnValues`) are tightly coupled to `glob.R_V_*`
  and reasonably live in `build_global_sys.cpp`; no duplicated helpers found.
- `Faddeeva.{cpp,hh}` left untouched (third-party).

## Latent issues identified but deliberately NOT acted on (for separate review)

1. **`init_BL.cpp` `#ifndef USE_CODIPACK` branch** (lines ~37–78): a second
   `solve_linear_system` using a direct `colPivHouseholderQr` solve. If no
   non-CODIPACK build target exists, this branch is dead. Behaviour-neutral to
   remove, but it is `#ifdef`-guarded build-variant code — confirm the build
   matrix before deleting.
2. **`extract_BL_TE.hpp` finite-difference step** (`h = 1e-6`) carries an
   existing `TODO: verify step is correct (convergence)` — a potential accuracy
   concern in a derivative approximation feeding the acoustic inputs. Not a
   readability item; flagged for a numerical review (changing it moves numbers).
3. **`newAmiet.hpp` denominator clipping** (several `TODO: replace clipping with
   limiting forms` near G_c/G_d/G_e and D±2k = 0): a robustness/accuracy item in
   the Amiet radiation integral; out of scope here.
4. **`residuals_shared.hpp` repeated-closure recompute** and the analogous note
   in `init_BL.cpp`/`get_funcs.hpp`: a performance opportunity (cache Hk/Ret/Mach2
   and the fixed start-node state across inner-Newton steps), explicitly left
   as-is and marked with `NOTE:`.

None of these are behaviour bugs in the current functional state; they are
documentation-flagged for a future, separately-reviewed change.

## Commits (all "readability only, golden 10/10 bit-identical")
1. `NOMENCLATURE.md + document get_funcs.hpp`
2. `file-header blocks to core solver/acoustic files`
3. `rename aux1/aux2 -> wgap1/wgap2`
4. `section banners in get_funcs.hpp + drop stale init_BL TODO`
5. `remove dead commented debug lines + typo fix in init_BL.cpp`
