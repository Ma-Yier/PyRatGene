# RatGenePy

[English](README.md) | **简体中文**

**RatGene** 的 Python 实现：面向约束型代谢网络的、基于生长-生产比率的基因与反应改造算法。

RatGene 用于搜索能够将目标代谢物生产与细胞生长偶联起来的基因/反应删除（以及可选的添加）策略。它支持基于 GPR 规则的基因逻辑、两步式策略规模缩减、模型整合，以及目标产率的最坏/最好情况分析。

本包是 [`RatGene-main/`](RatGene-main/) 中 MATLAB 实现的 Python 移植版本。默认使用 **SCIP**（通过 [`PySCIPOpt`](https://github.com/scipopt/PySCIPOpt)）；如果安装了相应的 Python API，也可以使用 **IBM ILOG CPLEX** 或 **Gurobi**。原 MATLAB 实现使用 CPLEX。

---

## 功能特性

- 基因删除策略与基因删除/添加策略
- 反应删除策略与反应删除/添加策略
- GPR（基因-蛋白-反应）规则解析与布尔约束
- 目标交换反应检测或自动添加
- 基于 MILP 解池的候选策略验证
- 可切换的 LP/MILP 求解后端：SCIP（默认）、CPLEX 和 Gurobi
- 两步式基因策略规模缩减
- 核心模型与添加模型的整合
- 目标生产的最坏/最好情况分析
- 无需 MATLAB 即可读取 BiGG 风格的 `.mat` 模型
- 提供 Python API、命令行接口，以及向后兼容的 `RatGene.py` 模块

---

## 环境依赖

- Python >= 3.9
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- 至少一个 LP/MILP 求解后端：
  - [PySCIPOpt](https://github.com/scipopt/PySCIPOpt)（默认安装）
  - [CPLEX Python API](https://pypi.org/project/cplex/)
  - [Gurobi Python API](https://www.gurobi.com/)（`gurobipy`）

默认求解后端为 SCIP。较新的 PySCIPOpt wheel 通常已经内置 SCIP；如果所在平台没有内置，请单独安装 SCIP。CPLEX 和 Gurobi 需要各自的产品安装和许可证（CPLEX 的 PyPI 包提供 Community Edition 运行时）。

已在 Python 3.12、NumPy 2.1.3、SciPy 1.14.1 和 PySCIPOpt 5.2.1 环境下测试。

---

## 安装

### 标准安装

在项目根目录执行：

```bash
python setup.py install
```

### 推荐：使用 pip 安装

```bash
pip install .
```

### 开发模式安装

```bash
pip install -e .
```

开发模式安装后，修改源码即可直接生效，无需反复安装。

### 可选求解后端

SCIP 会随包自动安装。如需使用 Gurobi，请安装 `gurobipy`，并确保 Gurobi 许可证可用：

```bash
pip install gurobipy
# 或使用可选依赖：
pip install ".[gurobi]"
```

如需使用 CPLEX，请安装其 Python API。PyPI 包提供 CPLEX Community Edition 运行时；也可以使用完整的 IBM ILOG CPLEX Optimization Studio 安装：

```bash
pip install cplex
# 或使用可选依赖：
pip install ".[cplex]"
```

### 验证安装

```bash
python -c "import ratgenepy; print(ratgenepy.__version__)"
ratgene --help
python -m ratgenepy --help
```

`ratgene` 命令行工具由 `setup.py` 自动安装。

---

## 快速开始

### Python API

```python
from ratgenepy import Model, RatGene

# 读取 BiGG 风格的 .mat 模型。
model = Model.load_mat("data/iMM904.mat")

# target_met 必须是 model.mets 中的代谢物 ID。
target_met = "glycerol_c"  # 请替换为 model.mets 中的实际 ID

# 生物量反应会根据模型目标函数自动识别。
# 如需显式指定，可传入 biomass="BIOMASS"。
x_target, knockouts = RatGene(
    model,
    target_met,
    carbon="EX_glc__D_e",
    oxygen="EX_o2_e",
)

print(f"目标产率: {x_target:.4f} mmol/gDW/h")
print(f"敲除数量: {len(knockouts)}")
print(knockouts)
```

可以通过 `backend="scip"`（默认）、`backend="cplex"` 或 `backend="gurobi"`
选择求解后端：

```python
from ratgenepy import available_solvers

print(available_solvers())  # 例如 ["scip", "gurobi"]

x_target, knockouts = RatGene(
    model,
    target_met,
    biomass="BIOMASS",
    backend="gurobi",
)
```

比率搜索、候选验证和策略缩减都会使用同一个求解后端。

如果是删除/添加问题，需要再提供一个用于添加信息的模型：

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

print("删除:", deletions)
print("添加:", additions)
```

### 反应层级策略

```python
x_target, knockouts = RatGene(
    model,
    target_met,
    type=0,        # 0 = 反应策略，1 = 基因策略
    size=0,        # 规模缩减仅适用于基因策略
)
```

### 命令行

```bash
ratgene path/to/model.mat TARGET_MET_ID \
    --biomass BIOMASS \
    --carbon EX_glc__D_e \
    --oxygen EX_o2_e
```

常用选项：

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

命令行还支持：

```text
--type gene|reaction|1|0
-s, --size true|false|True|False|1|0
-t, --timeLimit SECONDS
-p, --pool N
--backend scip|cplex|gurobi   （别名：--solver）
```

原有的单文件入口仍然可用：

```bash
python RatGene.py path/to/model.mat TARGET_MET_ID
```

在 Python 中也可以保持原有写法：

```python
import RatGene as rg

rg.RatGene(model, target_met)
```

---

## `RatGene` 参数

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
)
```

| 参数 | 默认值 | 说明 |
|---|---|---|
| `model` | 必填 | `Model` 对象，或 MATLAB struct 形式的模型 |
| `targetMet` | 必填 | 目标代谢物 ID，必须存在于 `model.mets` 中 |
| `biomass` | `None` | 生物量反应 ID；省略时使用模型目标反应 |
| `carbon` | `"EX_glc__D_e"` | 碳源交换反应 |
| `oxygen` | `"EX_o2_e"` | 氧源交换反应 |
| `LBbiomass` | `0.05` | 生物量反应下界 |
| `LBcarbon` | `-15.0` | 碳源交换反应下界 |
| `LBoxygen` | `-15.0` | 氧源交换反应下界 |
| `maxLoop` | `1000` | 比率搜索的最大迭代次数 |
| `timeLimit` | `inf` | 最大计算时间（秒） |
| `pool` | `10` | 用于验证的 MILP 解池大小 |
| `type` | `1` | `1` = 基因策略，`0` = 反应策略 |
| `size` | `1` | `1` = 缩减基因策略规模，`0` = 不缩减 |
| `addition` | `None` | 删除/添加问题中提供添加信息的模型 |
| `backend` | `"scip"` | LP/MILP 求解后端：`"scip"`、`"cplex"` 或 `"gurobi"` |

### 返回值

- `type=1`，无添加模型：
  `(x_target, knockouts)`
- `type=1`，有添加模型：
  `(x_target, deletions, additions)`
- `type=0`，无添加模型：
  `(x_target, knockouts)`
- `type=0`，有添加模型：
  `(x_target, deletions, additions)`

`x_target` 表示应用策略后的目标交换反应速率。基因策略向量中，`1` 表示基因保留/活跃，`0` 表示基因被删除；添加策略中，`1` 表示基因被添加。

> **注意：** 策略规模缩减（`size=1`）仅适用于基因策略。反应策略会直接返回，不执行缩减。

---

## 模型加载

`Model.load_mat` 可读取 MATLAB 保存的 BiGG/COBRA 格式 `.mat` 文件：

```python
from ratgenepy import Model

model = Model.load_mat("data/iMM904.mat")
print(model)
```

如果 `.mat` 文件中包含多个变量，可以显式指定模型变量名：

```python
model = Model.load_mat("data/models.mat", model_name="model")
```

也可以通过 SciPy 读取 MATLAB struct，再转换为 `Model`：

```python
import scipy.io as sio
from ratgenepy import Model

data = sio.loadmat("data/iMM904.mat", struct_as_record=False, squeeze_me=True)
model = Model.from_mat(data["model"])
```

### `Model` 字段

`Model` 容器与 MATLAB 代码使用的字段保持一致：

```text
S, b, c, lb, ub,
mets, metNames, metFormulas, metCharge,
rxns, rxnNames, subSystems,
grRules, genes, rev, rxnGeneMat, description
```

常用方法：

```python
model.copy()                 # 深拷贝
Model.from_mat(struct)       # 从 MATLAB struct / dict 构建
Model.load_mat(filename)     # 从 .mat 文件读取
```

---

## 底层 API

主算法也可以按步骤单独调用。需要注意：Python 使用 **0 基索引**，而 MATLAB 使用 1 基索引。

### 模型整合与交换反应

```python
from ratgenepy import uniteGeneModel, introExchange

united = uniteGeneModel(core_model, addition_model)
new_model, id_target, TMPR = introExchange(
    united, id_biomass, [id_carbon, id_oxygen], id_target_met,
    backend="scip",
)
```

- `uniteGeneModel(model1, model2)`：整合两个模型及其基因/反应关系。
- `introExchange(...)`：查找已有的 `EX_<met>` 或 `DM_<met>` 反应；如果不存在，则添加名为 `EX_target` 的目标交换反应。
- `TMPR`：目标代谢物的理论最大生产速率。

### 比率搜索

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
    backend="scip",   # 或 "cplex" / "gurobi"
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

- `ratioGene`：处理基因层级策略。
- `ratioMethod`：处理反应层级策略。
- `alpha`：搜索结束时使用的目标/生物量比率。
- `knockout`：策略向量。

对于基因策略，常用的 `gap` 计算方式为：

```python
gap = TMPR / (max_loop * LBbiomass)
```

### 策略规模缩减

```python
from ratgenepy import rmRedundancy

pre, after, reduced = rmRedundancy(
    model,
    knockout,
    id_target,
    id_biomass,
    split,
    backend="scip",   # 或 "cplex" / "gurobi"
)
```

- 纯删除问题：`split=[0, 0]`，返回删除数量的整数值。
- 删除/添加问题：
  `split=[len(core.rxns), len(core.genes)]`，`pre` 和 `after` 返回
  `[删除数, 添加数]` 数组。

### 边界分析

```python
from ratgenepy import boundAnalysis

worst, best = boundAnalysis(
    model, knockout, id_biomass, id_target, backend="scip"
)
```

返回指定基因策略下目标反应速率的最坏情况和最好情况。

### GPR 与验证辅助函数

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

- `maxParseGPR(model)`：将 GPR 规则解析为布尔关系。
- `constructMatrix(model)`：构建 `ratioGene` 所需的矩阵。
- `verifyRatioGene(model, gene_values)`：将基因策略转换为反应敲除向量。
- `verifyGeneKnock(...)`、`verifyKnock(...)`：验证基因策略和反应策略。

### 求解器接口

`ratgenepy.solvers` 为 SCIP、CPLEX 和 Gurobi 提供统一接口：

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

print(available_solvers())  # 当前环境中可导入的后端

handle, x = solver(
    objfunc,
    backend="gurobi",  # "scip"（默认）、"cplex" 或 "gurobi"
    Aeq=Aeq,
    beq=beq,
    A=A,
    b=b,
    lb=lb,
    ub=ub,
    opt="milp",        # "lp" 或 "milp"
    vtype="B" * n,
    pool=10,           # 保留最多 10 个 MILP 候选解
    timeLimit=60,
)

# 也可以直接调用特定后端：
handle, x = cplex_solver(objfunc, ...)
handle, x = gurobi_solver(objfunc, ...)
handle, x = scip_solver(objfunc, ...)

result = optimizeCbModel(model, backend="gurobi")
print(result.stat, result.f, result.x)

candidates = pool_solutions(handle, limit=10)
```

`solver` 始终最小化线性目标函数，并返回 `(handle, x)`。其中 `handle` 是后端模型对象，`x` 是最优解向量；如果问题不可行则为 `None`。`pool` 参数用于配置 CPLEX/Gurobi 解池；可以通过 `pool_solutions(handle, limit)` 按目标值从优到劣获取候选解。

`gurobipy` 和 CPLEX Python API 都是可选依赖。如果选择的后端不可用，`solver` 会抛出带有说明的 `ImportError`。

### 兼容包装函数

保留了此前 Python 草稿中的反应层级 API：

```python
from ratgenepy import RatMethodRxn, verifyRxn

x_target, knockout, alpha = RatMethodRxn(
    model, id_biomass, id_target, max_loop, gap, time_limit
)
x_target, deleted_reactions = verifyRxn(
    model, knockout, id_biomass, id_target, LBbiomass
)
```

这两个兼容包装函数也支持可选的 `backend="scip"|"cplex"|"gurobi"` 参数。

---

## Toy 示例

仓库中包含与 MATLAB `test3.m` 相同的 toy 模型：

```text
RatGene-main/toyData/toyCore.mat
RatGene-main/toyData/toyEdge.mat
```

```python
import numpy as np

from ratgenepy import Model, uniteGeneModel, ratioGene, rmRedundancy

core = Model.load_mat("RatGene-main/toyData/toyCore.mat")
edge = Model.load_mat("RatGene-main/toyData/toyEdge.mat")

# 整合核心模型与边缘模型。
model = uniteGeneModel(core, edge)

id_biomass = int(np.where(model.c != 0)[0][0])
id_target = model.rxns.index("R11")
model.lb[id_biomass] = 0.05

# 基于比率的基因策略搜索。
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

print("目标产率:", x_target)
print("删除基因:", deletions)
print("添加基因:", additions)

# 缩减策略规模。
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

print("缩减情况:", pre, "->", after)
print("缩减后的删除基因:", reduced_deletions)
print("缩减后的添加基因:", reduced_additions)
```

MATLAB/CPLEX 在 `test3.m` 中的参考输出为：

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

> **求解器说明：** RatGenePy 默认使用 SCIP，而 MATLAB 版本使用 CPLEX。可以通过 `backend="cplex"` 或 `backend="gurobi"` 选择 CPLEX 或 Gurobi。MILP 解池依赖于具体求解器，因此不同实现可能返回不同但都有效的策略。策略缩减过程本身遵循 MATLAB 的逻辑。

---

## 算法流程

1. **目标交换反应设置**：`introExchange` 查找或创建目标交换反应，并计算理论最大生产速率。
2. **GPR 矩阵构建**：`maxParseGPR` 和 `constructMatrix` 将基因-蛋白-反应规则转换为线性布尔约束。
3. **基于比率的 MILP 搜索**：`ratioGene` 或 `ratioMethod` 逐步提高目标/生物量比率并求解 MILP，直到找到通过验证的策略。
4. **验证**：通过最大化生物量并确认目标生产来验证候选策略（`verifyGeneKnock`、`verifyKnock`）。
5. **可选策略缩减**：`rmRedundancy` 在保持反应敲除模式不变的前提下，移除不必要的基因删除/添加。
6. **可选边界分析**：`boundAnalysis` 给出最终策略下目标产率的最坏和最好情况。

---

## 包结构

```text
ratgenepy/
├── __init__.py         # 公共 API 与版本号
├── __main__.py         # python -m ratgenepy
├── _utils.py           # MATLAB 数据转换、矩阵与索引辅助函数
├── model.py            # Model 容器与 .mat 模型加载
├── solvers/
│   ├── __init__.py     # 后端分发、solver、pool_solutions
│   ├── _common.py      # LP/MILP 参数标准化
│   ├── scip_backend.py  # SCIP 后端
│   ├── cplex_backend.py # CPLEX 后端
│   └── gurobi_backend.py # Gurobi 后端
├── scip.py             # 向后兼容的 SCIP 接口
├── gpr.py              # GPR 解析、constructMatrix、verifyRatioGene
├── verification.py     # verifyGeneKnock、verifyKnock
├── exchange.py         # introExchange
├── candidates.py       # 候选删除/添加集合
├── ratio.py            # ratioGene、ratioMethod、兼容包装函数
├── reduction.py        # rmRedundancy 与策略缩减
├── integrate.py        # findDiffGene、uniteGeneModel
├── analysis.py         # boundAnalysis
├── core.py             # RatGene 主函数
└── cli.py              # 命令行接口
```

项目根目录的 `RatGene.py` 是兼容层，会重新导出包内函数，并保留原来的单文件导入方式。

---

## MATLAB 兼容性说明

- Python 使用 **0 基索引**；MATLAB 索引需要减 1。
- 函数名和参数顺序尽量与 MATLAB 文件保持一致。
- `ratioGene` 对应 `ratioGene.m`；`ratioMethod` 对应 `ratioMethod.m`。
- `rmRedundancy` 对应 `rmRedundancy.m` 中的两步缩减流程。
- 对于反应策略，`rmRedundancy` 保持 MATLAB 的提前返回行为
  `(-1, -1, [inf])`；`RatGene` 不会对反应策略调用它。
- MILP 解池可能因所选后端（SCIP、CPLEX 或 Gurobi）不同而产生不同但有效的策略。
- Python 包不依赖 MATLAB。默认使用 SCIP；如果安装了相应的 Python API，也可以使用 CPLEX 或 Gurobi。

---

## 引用

如果您在研究中使用了 RatGene，请引用：

> Ma Y, Tamura T. **RatGene: Gene deletion-addition algorithms using growth to
> production ratio for growth-coupled production in constraint-based metabolic
> networks.** IEEE Transactions on Computational Biology and Bioinformatics,
> 2025.

---

## 许可证

MIT 许可证，详见 `setup.py`。
