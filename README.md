# RatGenePy

**English** | [简体中文](README.zh-CN.md)

Python implementation of **RatGene** — growth-to-production ratio-based gene
and reaction modification algorithms for constraint-based metabolic networks.

RatGene searches for gene/reaction deletion (and optional addition) strategies
that couple target metabolite production to cell growth. It supports GPR-based
gene rules, a two-step strategy-size reduction, model integration, and
worst-case/best-case target-rate analysis.

This package is a Python port of the MATLAB implementation in
[`RatGene-main/`](https://github.com/Ma-Yier/RatGene). It supports **SCIP** through
[`PySCIPOpt`](https://github.com/scipopt/PySCIPOpt) as the default solver, and
can also use **IBM ILOG CPLEX** or **Gurobi** when their Python APIs are
installed. The original MATLAB implementation uses CPLEX.

---

## Features

- Gene deletion and gene deletion/addition strategies
- Reaction deletion and reaction deletion/addition strategies
- GPR (gene–protein–reaction) rule parsing and Boolean constraints
- Target exchange detection or automatic introduction
- MILP solution-pool based candidate verification
- Pluggable LP/MILP backends: SCIP (default), CPLEX, and Gurobi
- Two-step reduction of gene strategy size
- Integration of a core model and an addition model
- Worst-case / best-case target production analysis
- Loading BiGG-style `.mat` models without MATLAB
- Python API, command-line interface, and backward-compatible `RatGene.py`
  module

---

## Requirements

- Python >= 3.9
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- At least one LP/MILP backend:
  - [PySCIPOpt](https://github.com/scipopt/PySCIPOpt) (installed by default)
  - [CPLEX Python API](https://pypi.org/project/cplex/)
  - [Gurobi Python API](https://www.gurobi.com/) (`gurobipy`)

SCIP is the default backend. Recent PySCIPOpt wheels usually bundle SCIP; if
your platform does not, install SCIP separately. CPLEX and Gurobi require
their own installations/licenses (the CPLEX PyPI package provides the
Community Edition runtimes).

Tested with Python 3.12, NumPy 2.1.3, SciPy 1.14.1, and PySCIPOpt 5.2.1.

---

## Installation

### Standard installation

From the project root:

```bash
python setup.py install
```

### Recommended: pip installation

```bash
pip install .
```

### Development installation

```bash
pip install -e .
```

The development install makes the package importable while editing the source
tree.

### Optional solver backends

SCIP is installed automatically. To use Gurobi, install `gurobipy` and make
sure a valid Gurobi license is available:

```bash
pip install gurobipy
# or, with the optional extra:
pip install ".[gurobi]"
```

To use CPLEX, install its Python API. The PyPI package provides the CPLEX
Community Edition runtimes, while a full IBM ILOG CPLEX Optimization Studio
installation can also be used:

```bash
pip install cplex
# or, with the optional extra:
pip install ".[cplex]"
```

### Verify the installation

```bash
python -c "import ratgenepy; print(ratgenepy.__version__)"
ratgene --help
python -m ratgenepy --help
```

The `ratgene` console command is installed by `setup.py`.

---

## Quick start

### Python API

```python
from ratgenepy import Model, RatGene

# Load a BiGG-style .mat model.
model = Model.load_mat("data/iMM904.mat")

# target_met must match an entry in model.mets.
target_met = "glycerol_c"  # replace with a metabolite ID from model.mets

# The biomass reaction is auto-detected from the model objective.
# Pass biomass="BIOMASS" to override it explicitly.
x_target, knockouts = RatGene(
    model,
    target_met,
    carbon="EX_glc__D_e",
    oxygen="EX_o2_e",
)

print(f"Target rate: {x_target:.4f} mmol/gDW/h")
print(f"Number of knockouts: {len(knockouts)}")
print(knockouts)
```

Select a solver backend with `backend="scip"` (default), `backend="cplex"`, or
`backend="gurobi"`:

```python
from ratgenepy import available_solvers

print(available_solvers())  # e.g. ["scip", "gurobi"]

x_target, knockouts = RatGene(
    model,
    target_met,
    biomass="BIOMASS",
    backend="gurobi",
)
```

The same backend is used for the ratio search, candidate verification, and
strategy reduction.

For a deletion/addition problem, provide a second model:

```python
from ratgenepy import Model, RatGene

core = Model.load_mat("data/iMM904.mat")
addition = Model.load_mat("data/iJR904.mat")

x_target, deletions, additions = RatGene(
    core,
    "formate_c",
    biomass="BIOMASS",
    addition=addition,
)

print("Deletions:", deletions)
print("Additions:", additions)
```

### Reaction-level strategies

```python
x_target, knockouts = RatGene(
    model,
    target_met,
    type=0,        # 0 = reaction, 1 = gene
    size=0,        # reduction is only available for gene strategies
)
```

### Command line

```bash
ratgene path/to/model.mat TARGET_MET_ID \
    --biomass BIOMASS \
    --carbon EX_glc__D_e \
    --oxygen EX_o2_e
```

Useful options:

```bash
ratgene model.mat TARGET_MET_ID \
    --type gene \
    --size true \
    --maxLoop 1000 \
    --pool 10 \
    --timeLimit 3600

ratgene model.mat TARGET_MET_ID \
    --type reaction \
    --size false \
    --backend cplex

ratgene core.mat TARGET_MET_ID \
    --biomass BIOMASS \
    --addition addition.mat \
    --backend gurobi
```

The CLI also accepts:

```text
--type gene|reaction|1|0
-s, --size true|false|True|False|1|0
-t, --timeLimit SECONDS
-p, --pool N
--backend scip|cplex|gurobi   (alias: --solver)
--integerTolerance VALUE
```

The original single-file entry point remains available:

```bash
python RatGene.py path/to/model.mat TARGET_MET_ID
```

and from Python:

```python
import RatGene as rg

rg.RatGene(model, target_met)
```

---

## `RatGene` parameters

```python
RatGene(
    model,
    targetMet,
    biomass=None,
    carbon="EX_glc__D_e",
    oxygen="EX_o2_e",
    LBbiomass=0.05,
    LBcarbon=-15.0,
    LBoxygen=-15.0,
    maxLoop=1000,
    timeLimit=float("inf"),
    pool=10,
    type=1,
    size=1,
    addition=None,
    backend="scip",
    integerTolerance=None,
)
```

| Parameter | Default | Description |
|---|---|---|
| `model` | required | A `Model` object or a MATLAB struct-like model |
| `targetMet` | required | Target metabolite ID, must be in `model.mets` |
| `biomass` | `None` | Biomass reaction ID; if omitted, the objective reaction is used |
| `carbon` | `"EX_glc__D_e"` | Carbon source exchange reaction |
| `oxygen` | `"EX_o2_e"` | Oxygen source exchange reaction |
| `LBbiomass` | `0.05` | Lower bound of the biomass reaction |
| `LBcarbon` | `-15.0` | Lower bound of the carbon source exchange reaction |
| `LBoxygen` | `-15.0` | Lower bound of the oxygen source exchange reaction |
| `maxLoop` | `1000` | Maximum iterations of the ratio-based search |
| `timeLimit` | `inf` | Maximum computation time in seconds |
| `pool` | `10` | Number of MILP solutions considered for verification |
| `type` | `1` | `1` = gene strategy, `0` = reaction strategy |
| `size` | `1` | `1` = reduce gene strategy size, `0` = no reduction |
| `addition` | `None` | Additional model for a deletion/addition problem |
| `backend` | `"scip"` | LP/MILP backend: `"scip"`, `"cplex"`, or `"gurobi"` |
| `integerTolerance` | `None` | MILP integrality tolerance; `None` uses each backend's minimum |

### Output

- `type=1`, no addition:
  `(x_target, knockouts)`
- `type=1`, with addition:
  `(x_target, deletions, additions)`
- `type=0`, no addition:
  `(x_target, knockouts)`
- `type=0`, with addition:
  `(x_target, deletions, additions)`

`x_target` is the target exchange reaction rate after applying the strategy.
Gene strategy vectors use `1` for an active/present gene and `0` for a deleted
gene; addition entries use `1` for an added gene.

> **Note:** The strategy-size reduction (`size=1`) is only available for gene
> strategies. Reaction strategies are returned without reduction.

---

## Model loading

`Model.load_mat` reads a `.mat` file saved by MATLAB in BiGG/COBRA format:

```python
from ratgenepy import Model

model = Model.load_mat("data/iMM904.mat")
print(model)
```

If the `.mat` file contains multiple variables, specify the model variable:

```python
model = Model.load_mat("data/models.mat", model_name="model")
```

A model can also be built from a MATLAB struct loaded by SciPy:

```python
import scipy.io as sio
from ratgenepy import Model

data = sio.loadmat("data/iMM904.mat", struct_as_record=False, squeeze_me=True)
model = Model.from_mat(data["model"])
```

### `Model` fields

The `Model` container mirrors the fields used by the MATLAB code:

```text
S, b, c, lb, ub,
mets, metNames, metFormulas, metCharge,
rxns, rxnNames, subSystems,
grRules, genes, rev, rxnGeneMat, description
```

Useful methods:

```python
model.copy()                 # deep copy
Model.from_mat(struct)       # build from a MATLAB struct / dict
Model.load_mat(filename)     # load from a .mat file
```

---

## Low-level API

The main algorithm is also available as individual functions. Indices are
**0-based** in Python, unlike the 1-based MATLAB code.

### Model integration and exchange

```python
from ratgenepy import uniteGeneModel, introExchange

united = uniteGeneModel(core_model, addition_model)
new_model, id_target, TMPR = introExchange(
    united, id_biomass, [id_carbon, id_oxygen], id_target_met,
    backend="scip",
)
```

- `uniteGeneModel(model1, model2)` merges two models and their gene/reaction
  relationships.
- `introExchange(...)` finds an existing `EX_<met>` or `DM_<met>` reaction, or
  appends a new target exchange reaction named `EX_target`.
- `TMPR` is the theoretical maximum production rate of the target.

### Ratio-based search

```python
from ratgenepy import ratioGene, ratioMethod

x_target, alpha, knockout = ratioGene(
    model,
    id_biomass,
    id_target,
    TMGR,
    max_loop,
    gap,
    numMultiStrat,
    timeLimit,
    backend="scip",   # or "cplex" / "gurobi"
)

x_target, alpha, knockout = ratioMethod(
    model,
    id_biomass,
    id_target,
    max_loop,
    gap,
    timeLimit,
    backend="scip",
)
```

- `ratioGene` handles gene-level strategies.
- `ratioMethod` handles reaction-level strategies.
- `alpha` is the final target/biomass ratio used by the search.
- `knockout` is the strategy vector.

For gene strategies the typical gap is:

```python
gap = TMPR / (max_loop * LBbiomass)
```

### Strategy reduction

```python
from ratgenepy import rmRedundancy

pre, after, reduced = rmRedundancy(
    model,
    knockout,
    id_target,
    id_biomass,
    split,
    backend="scip",   # or "cplex" / "gurobi"
)
```

- Pure deletion problem: `split=[0, 0]`; returns integer deletion counts.
- Deletion/addition problem:
  `split=[len(core.rxns), len(core.genes)]`; returns arrays
  `[deletions, additions]` for `pre` and `after`.

### Bound analysis

```python
from ratgenepy import boundAnalysis

worst, best = boundAnalysis(
    model, knockout, id_biomass, id_target, backend="scip"
)
```

Returns the worst-case and best-case target reaction rates for a gene
strategy.

### GPR and verification helpers

```python
from ratgenepy import (
    maxParseGPR,
    locateGene,
    constructMatrix,
    verifyRatioGene,
    verifyGeneKnock,
    verifyKnock,
)
```

- `maxParseGPR(model)` parses GPR rules into Boolean relationships.
- `constructMatrix(model)` builds the matrices used by `ratioGene`.
- `verifyRatioGene(model, gene_values)` converts a gene strategy into a
  reaction-knockout vector.
- `verifyGeneKnock(...)` and `verifyKnock(...)` validate gene and reaction
  strategies.

### Solver interface

`ratgenepy.solvers` provides a common interface for SCIP, CPLEX, and Gurobi:

```python
from ratgenepy import (
    available_solvers,
    cplex_solver,
    gurobi_solver,
    optimizeCbModel,
    pool_solutions,
    scip_solver,
    solver,
)

print(available_solvers())  # backends importable in this environment

handle, x = solver(
    objfunc,
    backend="gurobi",  # "scip" (default), "cplex", or "gurobi"
    Aeq=Aeq,
    beq=beq,
    A=A,
    b=b,
    lb=lb,
    ub=ub,
    opt="milp",        # "lp" or "milp"
    vtype="B" * n,
    pool=10,           # keep up to 10 MILP candidates
    timeLimit=60,
    integerTolerance=None,  # None = backend minimum; 1e-9 also valid
)

# Backend-specific entry points are also available:
handle, x = cplex_solver(objfunc, ...)
handle, x = gurobi_solver(objfunc, ...)
handle, x = scip_solver(objfunc, ...)

result = optimizeCbModel(model, backend="gurobi")
print(result.stat, result.f, result.x)

candidates = pool_solutions(handle, limit=10)
```

`solver` always minimizes the linear objective and returns `(handle, x)`,
where `handle` is the backend model object and `x` is the best solution vector
or `None` if the problem is infeasible. The `pool` argument configures the
CPLEX/Gurobi solution pool; use `pool_solutions(handle, limit)` to retrieve
candidate solutions in best-objective-first order.

For MILP problems, `integerTolerance` controls the integrality tolerance. When
it is omitted, each backend is set to its minimum supported value:

| Backend | Parameter | Minimum |
|---|---|---|
| SCIP | `numerics/feastol` (closest equivalent) | `1e-17` |
| CPLEX | `mip.tolerances.integrality` | `0.0` |
| Gurobi | `IntFeasTol` | `1e-9` |

Pass `integerTolerance=1e-9` explicitly to use the MATLAB
`OPTIONS.IntegerTolerance` value on all three backends. Very tight tolerances
can increase runtime and expose numerical difficulties.

Installing `gurobipy` or the CPLEX Python API is optional. If a selected
backend is unavailable, `solver` raises an informative `ImportError`.

### Compatibility wrappers

The previous Python draft's reaction-level API is preserved:

```python
from ratgenepy import RatMethodRxn, verifyRxn

x_target, knockout, alpha = RatMethodRxn(
    model, id_biomass, id_target, max_loop, gap, time_limit
)
x_target, deleted_reactions = verifyRxn(
    model, knockout, id_biomass, id_target, LBbiomass
)
```

Both wrappers accept an optional `backend="scip"|"cplex"|"gurobi"` argument.

---

## Toy example

The repository includes the same toy models used by MATLAB `test3.m`:

```text
RatGene-main/toyData/toyCore.mat
RatGene-main/toyData/toyEdge.mat
```

```python
import numpy as np

from ratgenepy import Model, uniteGeneModel, ratioGene, rmRedundancy

core = Model.load_mat("RatGene-main/toyData/toyCore.mat")
edge = Model.load_mat("RatGene-main/toyData/toyEdge.mat")

# Integrate the core and edge models.
model = uniteGeneModel(core, edge)

id_biomass = int(np.where(model.c != 0)[0][0])
id_target = model.rxns.index("R11")
model.lb[id_biomass] = 0.05

# Ratio-based gene strategy search.
TMGR = 2.0
max_loop = 1000
gap = TMGR / model.lb[id_biomass] / max_loop

x_target, alpha, knockout = ratioGene(
    model,
    id_biomass,
    id_target,
    TMGR,
    max_loop,
    gap,
    10,
    1000,
)

gene_split = len(core.genes)
rxn_split = len(core.rxns)

deletions = [
    model.genes[i] for i in np.where(knockout[:gene_split] == 0)[0]
]
additions = [
    model.genes[gene_split + i]
    for i in np.where(knockout[gene_split:] == 1)[0]
]

print("Target rate:", x_target)
print("Deleted genes:", deletions)
print("Added genes:", additions)

# Reduce the strategy size.
pre, after, reduced = rmRedundancy(
    model,
    knockout,
    id_target,
    id_biomass,
    [rxn_split, gene_split],
)

reduced_deletions = [
    model.genes[i] for i in np.where(reduced[:gene_split] == 0)[0]
]
reduced_additions = [
    model.genes[gene_split + i]
    for i in np.where(reduced[gene_split:] == 1)[0]
]

print("Reduction:", pre, "->", after)
print("Reduced deletions:", reduced_deletions)
print("Reduced additions:", reduced_additions)
```

The MATLAB/CPLEX reference result printed by `test3.m` is:

```text
Deleted genes:
  gene1, gene2, gene3, gene4, gene6, gene7

Added genes:
  gene8, gene9

Reduced deletions: 6 -> 5
Reduced additions: 2 -> 2

Reduced deletion genes:
  gene1, gene2, gene4, gene6, gene7

Reduced addition genes:
  gene8, gene9
```

> **Solver note:** RatGenePy uses SCIP by default, while the MATLAB version
> uses CPLEX. You can select CPLEX or Gurobi with `backend="cplex"` or
> `backend="gurobi"`. MILP solution pools are solver-dependent, so the selected
> strategy can differ slightly between implementations even when both are
> valid. The reduction procedure itself follows the MATLAB logic.

---

## Algorithm overview

1. **Target exchange setup** — `introExchange` locates or creates the target
   exchange reaction and computes the theoretical maximum production rate.
2. **GPR matrix construction** — `maxParseGPR` and `constructMatrix` convert
   gene–protein–reaction rules into linear Boolean constraints.
3. **Ratio-based MILP search** — `ratioGene` or `ratioMethod` iteratively
   increases the target/biomass ratio and solves MILPs until a verified
   strategy is found.
4. **Verification** — candidate strategies are checked by maximizing biomass
   and confirming target production (`verifyGeneKnock`, `verifyKnock`).
5. **Optional strategy reduction** — `rmRedundancy` removes unnecessary gene
   deletions/additions while preserving the reaction-knockout pattern.
6. **Optional bound analysis** — `boundAnalysis` reports worst-case and
   best-case target rates for the final strategy.

---

## Package layout

```text
ratgenepy/
├── __init__.py         # public API and version
├── __main__.py         # python -m ratgenepy
├── _utils.py           # MATLAB conversion, matrix and indexing helpers
├── model.py            # Model container and .mat loading
├── solvers/
│   ├── __init__.py     # backend dispatch, solver, pool_solutions
│   ├── _common.py      # shared LP/MILP argument normalization
│   ├── scip_backend.py  # SCIP backend
│   ├── cplex_backend.py # CPLEX backend
│   └── gurobi_backend.py # Gurobi backend
├── scip.py             # backward-compatible SCIP interface
├── gpr.py              # GPR parsing, constructMatrix, verifyRatioGene
├── verification.py     # verifyGeneKnock, verifyKnock
├── exchange.py         # introExchange
├── candidates.py       # candidate deletion/addition sets
├── ratio.py            # ratioGene, ratioMethod, compatibility wrappers
├── reduction.py        # rmRedundancy and strategy reduction
├── integrate.py        # findDiffGene, uniteGeneModel
├── analysis.py         # boundAnalysis
├── core.py             # RatGene main function
└── cli.py              # command-line interface
```

`RatGene.py` at the repository root is a compatibility shim that re-exports the
package and preserves the original single-file import style.

---

## MATLAB compatibility notes

- Python uses **0-based** indices. Convert MATLAB indices with `-1`.
- Function names and argument order intentionally follow the MATLAB files.
- `ratioGene` follows `ratioGene.m`; `ratioMethod` follows `ratioMethod.m`.
- `rmRedundancy` follows the two-step reduction in `rmRedundancy.m`.
- For reaction strategies, `rmRedundancy` mirrors the MATLAB early return
  `(-1, -1, [inf])`; `RatGene` does not call it for reaction strategies.
- The MILP solution pool may produce a different valid strategy depending on
  the selected backend (SCIP, CPLEX, or Gurobi).
- No MATLAB installation is required. SCIP is the default backend; CPLEX and
  Gurobi are optional when their Python APIs are available.

---

## Citation

If you use RatGene in your work, please cite:

> Ma Y, Tamura T. **RatGene: Gene deletion-addition algorithms using growth to
> production ratio for growth-coupled production in constraint-based metabolic
> networks.** IEEE Transactions on Computational Biology and Bioinformatics,
> 2025.

---

## License

MIT License, as declared in `setup.py`.
