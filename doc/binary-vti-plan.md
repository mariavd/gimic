# Binary `.vti` output and `vti2csv` — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** GIMIC writes its image-data `.vti` files in binary by default,
`vti_format=ascii` restores the text format, and `vti2csv` turns either
into a CSV table.

**Architecture:** The two image-data writers in `vtkplot.f90` branch on
`settings%vti_ascii`; the ASCII branch is today's code, the binary branch
writes the XML header and one raw appended block through stream access.
`vti2csv` is a standalone standard-library Python script configured and
installed like `gimic-profile`.

**Tech Stack:** Fortran 2003 (stream I/O), Python 3 standard library,
CMake/CTest, runtest.

**Spec:** `doc/binary-vti-design.md`

## Global Constraints

- File names of all outputs are unchanged.
- No new build or run-time dependency (no numpy in `vti2csv`).
- ASCII output stays byte-for-byte what it is today.
- Binary: `version="1.0"`, `header_type="UInt64"`, Float32, raw appended
  data, byte order detected at run time.
- `.vtu` writers are not touched.

## Review Focus

- A 2D grid (one point along an axis) in binary: file must be readable and
  agree with the ASCII one. Test: Task 3, second case.
- `vti_format` with an unknown value: the run stops with a message naming
  the keyword. Test: Task 3.
- `vti2csv` given a file that is not image data (a `.vtu`, a missing file):
  a one-line error and non-zero exit, no traceback. Test: Task 3.
- `vti2csv` on an ASCII file gives the same table as on the binary one.
  Test: Task 3.
- Open-shell runs write tagged files through the same writer; covered by
  `open-shell/3d` in ASCII, and the writer has no tag-dependent code.

---

### Task 1: Keyword and binary writer

**Files:**
- Modify: `src/gimic.in` (Advanced section), `src/libgimic/settings.f90`,
  `src/fgimic/gimic.F90`, `src/fgimic/vtkplot.f90`

**Interfaces:**
- Produces: `settings%vti_ascii` (logical, default `.false.`), set from
  `Advanced.vti_format` (`'binary'` default, `'ascii'`; anything else is a
  fatal error). Signatures of `write_vtk_imagedata` and
  `write_vtk_vector_imagedata` are unchanged.
- Binary layout: point array `scalars` (1 component) or `vectors`
  (3 components), x running fastest; for vector files a cell array `foo`
  holding the per-cell norm. Each array in the appended block is an
  8-byte length followed by Float32 values.

- [ ] Add `adv.add_kw('vti_format', 'STR', ('binary',))`, the settings
  field, and the `getkw` call with validation.
- [ ] Add the binary branch to both writers; keep the ASCII code as the
  other branch without edits to its statements.
- [ ] Build; run `benzene/3d` by hand in both formats; check the ASCII
  files are identical to the references' format and the binary files are
  smaller.
- [ ] Commit.

### Task 2: `vti2csv`

**Files:**
- Create: `tools/vti2csv.py.in`
- Modify: `tools/CMakeLists.txt` (configure to `vti2csv`, install, and
  wherever `gimic-profile` is copied into the build `bin`)

**Interfaces:**
- `vti2csv FILE.vti [-o OUT]`; default output `FILE.csv`, `-o -` is
  standard output. Columns `x,y,z` then `jx,jy,jz` (3 components) or
  `value` (1 component); with several point arrays the array name is
  prefixed. Reads `format="appended"` raw data (UInt32 or UInt64 headers,
  either byte order, Float32 or Float64) and `format="ascii"`. Exit status
  1 with a one-line message on any unreadable input.

- [ ] Write the script; convert the two `benzene/3d` outputs of Task 1
  and compare by eye.
- [ ] Commit.

### Task 3: Tests

**Files:**
- Create: `test/vti-binary/test`
- Modify: `test/CMakeLists.txt`; the inputs of `benzene/2d`, `3d`,
  `vectors`, `3d-keyword-magnet`, `2d-keyword-magnet`, `divj`,
  `tropicity-cdens`, `open-shell/3d`, `open-shell/divj` (add
  `vti_format=ascii` to `Advanced`); the input template in
  `test/divj-consistency/test`.

- [ ] New test, same command-line interface as `divj-consistency`
  (`--binary-dir --work-dir --case`). For `benzene/3d` and `benzene/2d`:
  run with the default and with `vti_format=ascii`; assert the default
  file contains `format="appended"` and is smaller; convert every `.vti`
  of both runs with `vti2csv`; assert equal row counts and coordinates and
  values equal within `2e-6*|v| + 1e-12`. Assert a run with
  `vti_format=foo` fails and mentions `vti_format`. Assert `vti2csv` on a
  missing file and on a non-`.vti` file exits 1 without a traceback.
- [ ] Add `vti_format=ascii` to the reference tests; run the whole suite.
- [ ] Commit.

### Task 4: Documentation

**Files:**
- Modify: `doc/input.rst` (Advanced keyword), `doc/output.rst` (format
  note, remove `jvec*.txt` entries), `doc/scripts.rst` (`vti2csv`),
  `doc/visualization.rst`, `doc/testing.rst` if it describes the `.vti`
  comparison.

- [ ] Write, commit, push.
