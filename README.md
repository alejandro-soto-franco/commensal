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

### Bimera detection options

Named as dada2 names them, so a run can be matched directly against it.

| option | default | meaning |
|---|---|---|
| `min_fold_parent_over_abundance` | 2 | a parent must be this many times as abundant as the query |
| `min_parent_abundance` | 8 | a parent must be at least this abundant outright |
| `max_shift` | 16 | how far a parent may be offset against the query |

## Binding against activation for flagellin at TLR5

Against the D0 swap of Clasen et al. 2023, which makes `RhFlaB` a hundredfold
more stimulatory while leaving D1 binding untouched, `receptor` decomposes the
observed ratio as:

| contribution | factor |
|---|---|
| observed | 100x |
| the measured 1:1 lifetime supplies | 2.88x |
| the 2:2 complex must supply | 34.7x |

The 1:1 kinetics account for 23% of the gap on a log scale. Bell et al. 2024
propose that "rapid dissociation of the silent flagellin may limit the ability
for TLR5 to trigger a downstream immune response"; that effect is real and
minor. The rest sits in the lifetime of the assembled complex.

### Case against kinetic proofreading

If activation ran through `n` sequential races against the same 1:1
dissociation, as the T cell receptor does (McKeithan 1995), `n = 4.35` would
account for the whole ratio with no help from D0. The same exponent applied to
the `StFliC-PIM` mutant, whose dissociation is 57 times faster, predicts it to
be 4.4e7-fold less active. The PIM constructs are used as an activity
threshold in the classification, so they are measurably active, and a spread
of that size is absent. A single race survives; proofreading does not.

### Predictions

**Receptor density.** Raising TLR5 surface density collapses the agonist to
silent ratio from about 90-fold at 0.5 receptors per square micrometre to
under 3-fold at 200, because more encounters fit inside a shorter-lived
complex. Titrating TLR5 expression in a reporter line tests this directly.

**Dose response.** Both partners in the assembly encounter are occupied
receptors, so the assembly rate goes as the square of occupancy. The model
gives a log-log slope of 2.04 below saturation, falling as occupancy
saturates. A slope near two is the signature of the 2:2 architecture and
distinguishes it from anything that assembles on a single receptor.

**What to measure.** For adaptor recruitment at 0.1 per second, the agonist's
2:2 complex has to last about 10 s and the silent one about 0.15 s. That is a
single-molecule measurement, and it is the number the model most wants.

## Running the binary

```bash
cargo build --release --bin commensal
```

```bash
target/release/commensal amplicon reads_1.fastq.gz [reads_2.fastq.gz ...] --out table.tsv
target/release/commensal consortium
target/release/commensal receptor
```

`amplicon` accepts one or more FASTQ files, plain or gzipped, and reports the
read, unique sequence and bimera counts along with the time each step takes.
With `--out`, it writes a table of sequence, abundance and bimera flag to the
given path. `consortium` and `receptor` take no arguments and print the
tables above.

## Running the Python package

The `commensal-py` crate builds Python bindings with PyO3 and maturin. Build a
wheel and install it into a uv-managed environment. Never `maturin develop`:
it installs into whatever `VIRTUAL_ENV` points at, and uv adopts and prunes an
exported `VIRTUAL_ENV` for its own project environments.

```bash
uv venv --python 3.12 /home/venvs/commensal-py
cd ~/commensal
uv run --with 'maturin>=1.4,<2' maturin build --release -m crates/commensal-py/Cargo.toml
uv pip install --python /home/venvs/commensal-py/bin/python target/wheels/commensal-*.whl
```

```python
from commensal import amplicon, consortium, receptor

uniques = amplicon.dereplicate(reads)
flags = amplicon.flag_bimeras(uniques)
state = consortium.steady_state(input=0.4)
result = receptor.decompose(observed_fold=94.0)
```

## Sources

- `callahan-2016-dad-hig-sam`, DADA2, the reference implementation
  dereplication and bimera detection are measured against.
- `hesser-2024-syn-ana-cac`, the synbiotic paper whose methods name `dada2`
  v1.14 and whose reads the dada2 comparison uses, and the consortium and
  cohort the `consortium` module models.
- `clasen-2023-sil-rec-fla`, silent recognition of flagellins from human gut
  commensal bacteria by TLR5, for the D0 swap the `receptor` decomposition
  explains.
- `bell-2024-kin-str-com`, kinetic and structure-based comparisons of silent
  and stimulatory flagellin interactions with TLR5, for the measured
  constructs.
- `forstneric-2017-rol-cte-d0`, the role of the C-terminal D0 domain, for the
  architecture that assigns D0 the job of crosslinking the two receptors and
  prolonging the assembled complex rather than deciding whether it forms.
- `yoon-2012-str-bas-tlr`, the structural basis of TLR5-flagellin
  recognition, for the 2:2 assembly the model is built on.
- `mckeithan-1995-kin-pro-tce`, kinetic proofreading, for the T cell receptor
  architecture the case against kinetic proofreading rules out for TLR5.
- `kemter-2023-com-bac-sig`, commensal signalling through TLR5 and AhR.
