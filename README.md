# commensal

Gut commensal analysis: amplicon sequences, consortium dynamics and host
receptor kinetics, in Rust with Python bindings.

## Overview

The work splits into a sequence half, which determines which organisms a 16S
sample contains, and a kinetic half, which determines what they produce and
how the host responds. The two halves meet at flagellin: a defined consortium
has known flagellin sequences, and whether a given one activates TLR5 decides
what the host detects.

## Modules

The design specifies six modules. Three are built; three are planned.

- `amplicon` (built): reads to unique sequences. Dereplication collapses
  reads to unique sequences with counts; bimera detection flags PCR chimeras
  by testing whether a sequence reads as the left of one abundant parent and
  the right of another. Denoising, with an error model and banded alignment,
  is planned for a later phase.
- `community` (planned): taxonomy, diversity and the flagellome.
- `consortium` (built): colonisation and butyrate production for a defined
  gut consortium, modelled as a chemostat: substrate arrives, a strain grows
  on it by Monod kinetics, and butyrate comes off the fermentation while
  transit dilutes everything and the epithelium absorbs butyrate.
- `inference` (planned): distance, ordination, permanova and differential
  abundance.
- `phenotype` (planned): dose response, ED50 and power.
- `receptor` (built): binding against activation for flagellin at TLR5. A
  flagellin can bind tightly while leaving the receptor silent; the module
  splits an observed activation ratio into the part a measured 1:1 binding
  lifetime accounts for and the part the 2:2 assembly must supply.

## Measured comparison with dada2 1.38.0

One sample of PRJNA1091036, 54 058 reads, on one core:

| | commensal | dada2 1.38.0 |
|---|---|---|
| unique sequences | 25 989 | 25 989 |
| bimeras called | 2143 | 2192 |
| dereplication | 0.009 s | 1.451 s |
| bimera detection | 0.181 s | 50.877 s |
| total | 0.255 s | 52.328 s |

**Dereplication is exact.** The unique set and every abundance match, so that
step is a drop-in.

**Bimera detection agrees on 98.5% of sequences**, recovering 90% of the
reference's calls and adding 167 of its own. The residual is alignment: the
reference compares a parent to a query with a banded Needleman-Wunsch, which
absorbs an indel anywhere in the sequence, while commensal compares them at
every offset within `max_shift`, which absorbs one only at the ends. Closing
that gap requires implementing the alignment, and until it closes this stays a
fast near-agreement rather than a drop-in.

## Running the binary

```bash
cargo build --release --bin commensal
```

```bash
target/release/commensal amplicon reads_1.fastq.gz [reads_2.fastq.gz ...] --out table.tsv
```

It accepts one or more FASTQ files, plain or gzipped, and reports the read,
unique sequence and bimera counts along with the time each step takes. With
`--out`, it writes a table of sequence, abundance and bimera flag to the given
path.

## Running the Python package

The `commensal-py` crate builds Python bindings with PyO3 and maturin. From
`crates/commensal-py`, build a wheel and install it into a uv-managed
environment:

```bash
cd crates/commensal-py
uv venv
uv run --with maturin maturin build --release
uv pip install --python .venv/bin/python target/wheels/commensal-*.whl
```

```python
from commensal import amplicon, consortium, receptor

uniques = amplicon.dereplicate(reads)
flags = amplicon.flag_bimeras(uniques)
state = consortium.steady_state(input=0.4)
result = receptor.decompose(observed_fold=94.0)
```
