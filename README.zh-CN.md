[English](README.md) | **简体中文**

# ChemWeave

一个基于 RDKit 反应模板的 Python 反应网络生成包。输入初始物种、反应规则和碳数范围，输出物种表、反应表、计量关系和规则来源。规则与生成算法分开维护，便于随文献整理逐步扩展。

当前版本是网络生成工具。内置四条规则用于检验总体结构变换：烷烃脱氢、烯烃双键迁移、受限的端烯烃偶联和烯烃裂解。**这些示例没有覆盖完整芳构化网络，也不是经过验证的 Ga/H-ZSM-5 基元机理。**

## 安装与测试

需要 Python 3.11 或更新版本，以及 RDKit。可以使用已有环境，例如 Conda 的 `base` 或 `ml`。在项目的上级目录运行：

```bash
conda activate base
cd chemweave
python -m pip install -e .
python -m pytest -q
```

pytest 作为可选开发依赖声明在 `.[test]` 中。需要一起安装时，可运行 `python -m pip install -e ".[test]"`；已有 pytest 的环境可以使用上面的命令。

不安装时也可以在项目目录直接运行：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/propane.toml --output results/propane
```

## Python 接口

```python
from chemweave import Species, Generator, GenerationConfig, load_rules, Network

network = Generator(
    rules=load_rules(),
    config=GenerationConfig(max_carbon=6, max_species=200, max_reactions=1000),
).generate(seeds=[Species.from_smiles("CCC", name="propane")])

print(network.summary())
network.export_json("results/network.json")
network.export_csv("results")
network.export_dot("results/network.dot")

restored = Network.from_json("results/network.json")
```

`load_rules()` 读取内置示例；`load_rules("my_rules.json")` 读取自己的规则文件。`Species` 使用 `from_smiles` 或 `from_mol` 创建，以统一结构、分子式和唯一标识。

在 base 环境（Python 3.12、RDKit 2026.03.6）运行 C6 示例得到 22 个物种、45 条反应和 45 个对称归并通道。这个数字是运行示例，不作为化学正确性的判据。

## 自定义规则

规则使用 JSON 列表存储。以下例子明确表示丙烷脱氢的总体变换，包含 H2：

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[C;H3:1](-[H:4])-[C;H2:2](-[H:5])-[C;H3:3]>>[C:1]=[C:2]-[C:3].[H:4]-[H:5]",
    "family": "dehydrogenation",
    "kind": "overall",
    "version": "1",
    "source": "在此填写文献或规则依据",
    "notes": "总体变换，未展开催化循环"
  }
]
```

也可以用 `ReactionRule(...)` 在 Python 中直接定义。主要字段如下：

| 字段 | 含义 |
| --- | --- |
| `id` | 规则唯一名称；不同机理即使总反应相同，也可以使用不同规则 ID |
| `smarts` | RDKit Reaction SMARTS，原子映射编号用于定义结构变化 |
| `family` | 反应族标签 |
| `kind` | 步骤性质，默认 `overall`；经机理审查后可标记为 `elementary` |
| `version`、`source`、`notes` | 规则版本、依据和假设 |
| `exclude_smarts` | 任一反应物命中其中一个模式时，不应用此规则 |
| `reactant_phases`、`product_phases` | 按模板顺序指定各物种相态，默认全部为 `gas` |
| `reactant_sites`、`product_sites` | 按模板顺序指定各物种单个位点标签，默认全部为 `null` |

默认不自动反转规则。逆反应需要显式提供模板；某些双键迁移模板自身可以生成正反两个方向。规则被声明为 `elementary` 也不等于程序已经证明了它的机理。

## 网络如何扩展

1. 标准化种子结构，以结构、相态和位点标签确定物种身份。
2. 每个新物种只与每条规则的各个反应物位置匹配一次。
3. 将新加入候选池的物种与已有候选物种组合，执行尚未处理的有序反应物组合。允许同一物种占据两个位置，以表示自偶联。
4. 用显式氢执行模板，恢复来源映射号，检查原子一一对应（包括元素和同位素）、电荷和位点标签数守恒。
5. 按 `max_carbon` 筛选，并按规则 ID 和两侧计量合并重复反应。每种非等价映射只在 `channels` 中保存一个完整代表；等价执行仅累加其 `match_count`，不保留重复记录。
6. 把新物种加入队列，直到队列耗尽或达到物种/反应预算。

检查放在输入读取和候选反应生成处；正常枚举不重复校验整张网络。分子描述和反应模板在创建时处理一次。当前枚举不估算能垒或通量，也不会仅因某个物种未在实验中检出就删除它。

`max_carbon` 对每个物种生效，种子也必须处于这个范围内。它定义本次网络的范围，不表示范围外路径在实际体系中不重要。

## 表面状态与 Single-event 的边界

可以创建 `Species.from_smiles("C=CC", phase="surface", site="BAS")`。位点标签参与物种去重，反应规则可以限定输入和输出位点，并检查每类位点的数目守恒。

这是**单个位点标签的抽象表示**，适合相同位点上有明确分子结构的变换。当前未提供空位、质子交换、双位点物种以及 Ga 催化循环的专门表示或内置规则。真实表面结构及其氢、电荷需要在进一步扩展机理时一致定义；给气相分子增加标签本身不能建立真实吸附机理。

抽象位点标签不改变原子编号。比较反应通道时，相态和位点标签也参与等价性判断。

## 原子映射与对称性

物种身份由不带映射号的标准 SMILES、相态和位点确定。`species.mapped_smiles` 单独保存显式氢表示。已有且不冲突的映射号保留；缺失编号的原子（包括补出的氢）从已有最大编号之后依次分配。输入内部重复的编号会被消除冲突。`species.input_maps` 保存 `(分配后的编号, 输入编号)` 对，输入编号为零表示原先未提供标签。

```python
seed = Species.from_smiles("[CH3:101][CH2:205][CH3:309]")
print(seed.smiles)             # CCC: used for species identity
print(seed.mapped_smiles)      # Retains carbon maps 101, 205, 309; includes mapped H
print(dict(seed.input_maps))   # Assigned labels mapped back to supplied labels
print(seed.symmetry_classes)  # (map_number, symmetry_class) pairs
```

多分子反应中，每份反应物独立编号。与前一份反应物冲突的编号从全部预留编号之后重新分配。每个匹配记录的 `atom_mapping` 保存反应物位置、来源原子索引、物种内部编号、用户输入编号及本次反应编号。因此，同种物种的两份分子可以参与反应，而不会共享原子 ID。结构相同的种子输入记录在 `network.seed_inputs` 中；网络扩展对每种物种使用一个代表。

模板在显式氢分子上执行。涉及氢转移的规则必须明确移动带映射的氢原子及其化学键。若生成无法追溯的新原子、丢失来源原子或在产物中重复使用编号，候选会以 `atom_mapping` 原因被拒绝；电荷和位点守恒另行检查。旧的隐式氢模板可能需要改写，本版本不会自动猜测氢的对应关系。内置第 2 版规则提供总体反应的原子对应，并不代表已经验证了真实的基元氢转移机理。

`reaction.channels` 对每种非等价反应映射只保存一份完整代表。记录包括 `mapped_reaction_smiles`、两侧带映射的 SMILES 列表、`atom_mapping` 来源编号表、产物原子来源、`symmetry_key` 和 `match_count`。即使计量关系合并了两份相同物种，记录仍区分反应物实例顺序。遇到重复时，仅增加已有通道的 `match_count`，丢弃该次完整映射；不再单独保存 `reaction.matches` 或重复匹配索引列表。

例如，丙烷脱氢的 24 个模板匹配只保存为一个通道，其 `match_count=24`。摘要中的 `template_matches` 和 CSV 同名列为这些计数的总和，不是实际存储的记录数。

通道通过反应前后完整图的标准化表示归并。图中的分子节点保留分组、相态和位点；原子标签保留元素、同位素、电荷、自由基及已指定的原子 CIP 标签；化学键标签保留两侧键级和已定义的立体信息。这样不会仅凭原子各自属于同一对称类别，就认定多位点反应等价；化学上相同的反应物实例则允许整体交换。

基础原子类别通过 `CanonicalRankAtoms(breakTies=False, includeAtomMaps=False)` 计算，不枚举分子自匹配。显式氢新戊烷的类别大小为 1、4、12。真实映射号始终唯一，不会被对称类别标签覆盖。已定义的立体化学会保留；程序不穷举新生成的立体异构体。

这些是单步反应映射。编号属于物种代表或某次反应实例，并非全网唯一的原子身份。多步来源追踪需要连接各步映射，并将产物原子与物种库代表的原子对应起来，不能仅凭编号数值直接连接。通道数及 `match_count` 都不是 SEMK 简并度；物理路径计数和动力学参数共享需要单独处理。本版本不计算速率、转化率或选择性。

新的 JSON 导出使用第 3 版格式，保留精简后的通道代表、计数、编号、对称类别和种子输入。读取第 2 版文件时，会将原始匹配和通道索引转换为精简表示。第 1 版文件仍可作为结构网络读取，其中不完整的旧来源记录不会保留为映射通道；不能从旧的不带编号结构中恢复来源映射。

## 命令行与输出

安装后可运行：

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

可使用 `examples/mapped_propane.toml` 从部分编号的种子生成同样的网络。

自定义配置中的 `rules = "my_rules.json"` 相对于配置文件所在目录解析。输出目录相对于运行命令的目录解析。

| 文件 | 内容 |
| --- | --- |
| `network.json` | 完整物种、反应、通道映射代表与计数、规则、配置及 RDKit 版本，可重新载入 |
| `species.csv` | 物种 ID、名称、不带映射的结构、分子式、电荷、相态、位点及带映射的 SMILES |
| `reactions.csv` | 反应 ID、规则、反应族、步骤性质、两侧计量、模板匹配数和通道数 |
| `channels.csv` | 每个对称等价通道的一份代表性 mapped reaction 及其 `match_count` |
| `stoichiometry.csv` | 每条反应每一侧的物种及正整数计量数；生成计量矩阵时反应物取负 |
| `network.dot` | Graphviz 图；每条反应有独立节点，保留多反应物和多产物关系 |

DOT 导出不依赖 Graphviz；如已安装 Graphviz，可用 `dot -Tsvg results/propane/network.dot -o results/propane/network.svg` 绘图。

`termination="exhausted"` 表示在所选规则和碳数范围内没有待处理组合。`max_species` 或 `max_reactions` 表示预算中止，此时结果是部分网络，仍会保存。命令行正常完成退出码为 0，预算中止为 2，输入错误由 Python/RDKit 报告。

## 后续扩展

先用带来源和适用条件的规则逐步补齐裂解、支化、环化和芳烃生成路径，再比较提高碳数上限后的变化。随后根据选定的表面机理加入空位和催化剂状态表示，最后接入 Single-event 参数组织和动力学求解。

## 参考

- [RDKit 化学反应接口](https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions)
- [RDKit Reaction SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)
- [Python package 打包指南](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

