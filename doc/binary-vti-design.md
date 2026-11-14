# Binary `.vti` output and a `vti2csv` converter — design

## Why

Large `jvec.vti` files take many minutes to load in ParaView because every
number is written as text. The files should load quickly without the user
having to configure anything. Whoever needs the explicit numbers should
still be able to get them.

## Goals

- All image-data files (`jvec*`, `jmod*`, `acid`, `divj`, `tropicity`,
  `jvec_dia/para/undef`) are written in binary by default.
- One input keyword brings the old text format back.
- A command-line tool turns any `.vti` GIMIC writes into a CSV table.

## Non-goals

- The unstructured `.vtu` files stay as they are.
- No compression.
- No CSV output from GIMIC itself; the converter covers that.

## Anti-goals

- File names must not change: the ParaView state files and the
  documentation refer to them.
- No new build dependency, neither for GIMIC nor for the converter.

## Approaches

1. **Binary by default, ASCII behind a keyword, CSV through a separate
   converter** — chosen. The fast format needs no setup, existing tests
   keep their references, and GIMIC writes no extra files.
2. Binary only, tests decode the files before comparing. Covers the
   default path everywhere, but every `.vti` test and reference has to be
   redone and there is no way back to text for a user.
3. Binary as an opt-in. Rejected: the beginner who follows the
   documentation would still get the slow files.

## Design

### Writer (`src/fgimic/vtkplot.f90`)

- `write_vtk_imagedata` and `write_vtk_vector_imagedata` gain a binary
  path: the XML header as now, the data arrays declared with
  `format="appended"` and an offset, then one
  `<AppendedData encoding="raw">` block written with stream access.
- Numbers are stored as Float32. That is about seven digits, more than
  the six of the text format, at 12 bytes per vector instead of about 43.
- The file header declares `version="1.0" header_type="UInt64"`, so that
  an array larger than 4 GB does not overflow its length field.
- Byte order is detected at run time and declared in the header.
- The per-cell norm block of `jvec.vti` is kept, also in binary.
- The ASCII path is the present code, unchanged, so its output stays
  byte-for-byte what it is today.

### Keyword

`vti_format` in the `Advanced` section, `binary` (default) or `ascii`.
Declared in `src/gimic.in`, read in `src/fgimic/gimic.F90` into
`settings`. Any other value stops the run with a clear message.

### Converter (`tools/vti2csv.py.in`, installed as `vti2csv`)

- Python standard library only, like `gimic-profile`.
- Reads both the binary and the ASCII files.
- Writes one row per grid point: `x,y,z` in bohr (from the origin and
  spacing in the header) followed by `jx,jy,jz` for a vector file or
  `value` for a scalar file, with a header line.
- `vti2csv jvec.vti` writes `jvec.csv`; `-o FILE` names the output, `-o -`
  writes to standard output.
- Point data only; the per-cell norm is not exported.

### Tests

- The tests that compare `.vti` files against references set
  `vti_format=ascii` in their inputs. References stay untouched.
- `test/divj-consistency` also sets `vti_format=ascii`, since it parses
  `jvec.vti` as text.
- One new test runs the benzene 3D case in both formats, converts both
  with `vti2csv` and checks that they agree to the six digits of the text
  format. This covers the default path and the converter.
- `test/threads` compares files bitwise and needs no change.

### Documentation

`doc/output.rst`, `doc/input.rst`, `doc/scripts.rst`,
`doc/visualization.rst`: the binary default, the keyword, the converter.
The entries for `jvec.txt`, which the code no longer writes, are removed.

## Left to the builder

- Internal layout of the writer (shared helper for the header and the
  appended block).
- The exact number format of the CSV columns.
