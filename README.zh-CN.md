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

在 base 环境（Python 3.12、RDKit 2026.03.6）运行 C6 示例得到 22 个物种、45 条反应。这个数字是运行示例，不作为化学正确性的判据。

## 自定义规则

规则使用 JSON 列表存储。以下例子明确表示丙烷脱氢的总体变换，包含 H2：

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[CH3:1][CH2:2][CH3:3]>>[CH2:1]=[CH:2][CH3:3].[H][H]",
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
4. 标准化产物，检查价态、元素（包括氢）、电荷和位点标签数守恒。
5. 按 `max_carbon` 筛选，并按规则 ID、反应物计量、产物计量合并重复反应。等价模板匹配记录在 `matches` 中。
6. 把新物种加入队列，直到队列耗尽或达到物种/反应预算。

检查放在输入读取和候选反应生成处；正常枚举不重复校验整张网络。分子描述和反应模板在创建时处理一次。当前枚举不估算能垒或通量，也不会仅因某个物种未在实验中检出就删除它。

`max_carbon` 对每个物种生效，种子也必须处于这个范围内。它定义本次网络的范围，不表示范围外路径在实际体系中不重要。

## 表面状态与 Single-event 的边界

可以创建 `Species.from_smiles("C=CC", phase="surface", site="BAS")`。位点标签参与物种去重，反应规则可以限定输入和输出位点，并检查每类位点的数目守恒。

这是**单个位点标签的抽象表示**，适合相同位点上有明确分子结构的变换。当前未提供空位、质子交换、双位点物种以及 Ga 催化循环的专门表示或内置规则。真实表面结构及其氢、电荷需要在进一步扩展机理时一致定义；给气相分子增加标签本身不能建立真实吸附机理。

每条反应保留规则和模板匹配记录，记录包含反应物顺序、产物顺序、带原子索引标签的反应物 SMILES 及 RDKit 原始产物原子来源。索引 SMILES 的映射值为输入原子索引加一；来源中的 `react_atom_idx` 为零起始索引。H2 等模板新生成的原子可能没有输入来源。

这些是结构生成的追溯信息，**`len(reaction.matches)` 不是 Single-event 简并度**。转动对称数、等价路径的物理计数、过渡态和动力学参数共享需要后续独立实现。本版本不计算速率、转化率或选择性。

## 命令行与输出

安装后可运行：

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

自定义配置中的 `rules = "my_rules.json"` 相对于配置文件所在目录解析。输出目录相对于运行命令的目录解析。

| 文件 | 内容 |
| --- | --- |
| `network.json` | 完整物种、反应、模板匹配、规则、配置及 RDKit 版本，可重新载入 |
| `species.csv` | 物种 ID、名称、结构、分子式、电荷、相态和位点 |
| `reactions.csv` | 反应 ID、规则、反应族、步骤性质、两侧计量和模板匹配数 |
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

