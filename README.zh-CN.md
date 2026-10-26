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

在 macOS/Linux 中，不安装时也可以在项目目录直接运行：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/propane.toml --output results/propane
```

## 网页界面

Windows、macOS 和 Linux 使用相同的 Python 启动入口。在包含 `pyproject.toml` 的项目目录打开终端，使用 Python 3.11 或更新版本的环境。首次安装网页可选依赖，随后启动本地工作台：

```bash
python -m pip install -e ".[web]"
python -m chemweave web
```

安装完成且对应环境已激活时，也可以使用等价命令 `chemweave web`。使用网页期间保持终端打开。浏览器会自动打开 `http://127.0.0.1:8765`；如果没有自动打开，可手动输入该地址。网页界面使用英文。输入初始物种 SMILES 和可选名称，勾选规则，设置碳数、深度后点击“Generate”（生成）。可添加目标物种来记录实验产物。“Stop when all targets are found”（找到所有目标后停止）默认开启；寻找替代生成途径时请关闭；展开各物种的“Phase and site”（相态与位点）可以设置对应状态。

中间图使用分子结构节点和独立的反应节点，边上显示非 1 的计量系数。点击节点可以查看物种或反应详情、切换映射通道、显示该反应实际保存的 atom-map 编号。可以搜索、按深度筛选或查看选中节点附近。超过 80 个物种时，默认显示到第 2 层，可通过深度选项展开。筛选后只显示包含全部反应物和共生成物的完整反应。运行期间更新统计信息，结束后绘制网络。

“Stop”（停止）在候选反应之间请求停止；当前正在执行的 RDKit 底层匹配调用需要返回后，停止请求才能生效。已接受的反应记录保持完整，结果标记为 `termination="user_stopped"`，可以继续保存；匹配计数和替代路径可能尚未枚举完整。

“Save project”（保存项目）下载包含输入、所选规则、设置和网络的 JSON 文件。“Open”（打开）支持项目文件及现有 `network.json`，包含受支持的旧版格式。“Export”（导出）支持网络 JSON、CSV 表格及 DOT 的 ZIP 包，或当前网络视图的 PNG 图片。下载位置遵循浏览器的正常设置。

每个本地服务在内存中维护一个工作区，多个浏览器标签页共享它。刷新会恢复当前工作区；保存项目后可以跨服务重启恢复。“New”（新建）在确认后清空工作区。服务只监听本机回环地址；可在终端按 Ctrl+C 关闭。

如果希望使用独立环境，在项目目录中创建一次：

```text
python -m venv .venv-web
```

然后直接用该环境的 Python 安装和启动，无需先执行环境激活命令。

**Windows（PowerShell 或命令提示符）：**

```text
.\.venv-web\Scripts\python.exe -m pip install -e ".[web]"
.\.venv-web\Scripts\python.exe -m chemweave web
```

**macOS / Linux：**

```bash
.venv-web/bin/python -m pip install -e ".[web]"
.venv-web/bin/python -m chemweave web
```

以后使用时，只需运行对应环境的启动命令，不必重复创建环境和安装。如果系统中的 Python 命令是 `python3`，创建环境时使用 `python3`。已有 Conda 环境也可以直接使用本节开头的通用安装和启动命令。如果要复用创建虚拟环境时所在 Python 环境中的 RDKit 等包，可在 `venv` 命令中加入 `--system-site-packages`。

需要更换端口时运行 `python -m chemweave web --port 8766`；加上 `--no-browser` 表示不自动打开浏览器。使用独立环境时，将这里的 `python` 替换为上面对应系统的 Python 路径。在终端按 Ctrl+C 可以关闭服务。

后台使用 [FastAPI](https://fastapi.tiangolo.com/) 和 Uvicorn；分子 SVG 由 RDKit 绘制。[Cytoscape.js](https://js.cytoscape.org/) 3.34.4 及 MIT 许可证已随包提供，运行时不连接 CDN，不需要 Node.js 或前端构建。Python 接口和命令行生成仍可在未安装网页依赖时使用。

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

新的 JSON 导出使用第 4 版格式，保留精简后的通道代表、计数、编号、对称类别、种子输入、物种及反应深度和目标状态。第 3 版文件仍可读取；旧文件缺失的深度保持未知。读取第 2 版文件时，会将原始匹配和通道索引转换为精简表示。第 1 版文件仍可作为结构网络读取，其中不完整的旧来源记录不会保留为映射通道；不能从旧的不带编号结构中恢复来源映射。

## 深度与目标物种停止条件

设置 `GenerationConfig(max_depth=2)`，允许从种子出发最多经过两步连续反应。种子的深度为 0，反应深度定义为：

```text
反应深度 = 1 + 所有反应物深度的最大值
新产物深度 = 反应深度
```

每个物种保留其最小可达深度。FIFO 队列按深度递增发现物种，多种种子或双分子反应也适用。循环反应不会增加已登记物种的深度。对于双分子反应，深度衡量最长依赖链，不将两条并行制备支路的反应次数相加。生成已有浅层物种的反应，仍记录自身的反应深度。

达到深度上限的产物及其生成反应、映射都会保留，但不再从这些物种继续扩展。`max_depth=0` 只保留种子；不设置该参数则不限制深度。如果网络提前穷尽，则不必达到深度上限。使用内置规则时，丙烷深度为 0，丙烯和 H2 为 1，丙烯自偶联生成的 C6 烯烃为 2。可运行 `examples/depth_limited.toml`。

指定具体的目标物种，在它们**全部出现**时停止：

```python
network = Generator(
    load_rules(),
    GenerationConfig(max_carbon=6, max_depth=4),
).generate(
    [Species.from_smiles("CCC")],
    targets=[Species.from_smiles("C=CC"), Species.from_smiles("C=C")],
)
print(network.summary())
print(network.species_depths)  # 物种 ID -> 最小深度
print(network.targets)        # 目标结构及已发现的物种 ID
```

目标沿用物种库的身份判定：标准化结构、相态和位点相同，忽略 atom-map 编号和名称。重复目标只计一次；已在种子中的目标立即计为找到。指定目标不会将它加入网络。这是结构可达性条件，不判断实验浓度或检测限。`max_species` 仍限制网络中存储的物种总数。

TOML 使用与 `[[seeds]]` 相同字段的 `[[targets]]` 条目：

```toml
[[seeds]]
smiles = "CCC"

[[targets]]
smiles = "C=CC"

[[targets]]
smiles = "C=C"

[generation]
max_carbon = 6
max_depth = 4
```

`examples/target_species.toml` 提供带标签和数量预算的完整示例。两个条件可以同时设置：达到深度边界时停止；如果目标提前全部出现，则提前停止。摘要中提供 `depth_reached`、`targets_total`、`targets_found` 和 `missing_targets`。`network.targets` 每个目标的 `species_id` 为找到的物种 ID，未找到时为 `None`。

目标达成后，程序保存当前接受的完整反应、全部共生成物和满足条件的通道映射，随即返回。尚未遍历的匹配、通道和替代路径不会继续枚举，因此匹配计数可能不完整。这样的部分网络受种子和规则的枚举顺序影响，不表示已收集所有最短路径。

## 目标产物生成途径

生成或打开网络后，在左侧 **Formation pathways**（生成途径）中选择产物，设置 **Max reactions per route**（每条途径的最大反应数）和 **Max pathways**（最多返回的途径数），点击 **Analyze formation**。也可以在物种详情中直接点击 **Analyze formation**。为收集替代途径，生成网络前应关闭 **Stop when all targets are found**。分析只使用已有网络，不会继续扩展此前提前停止的网络。

选择单条途径，可以按一个可执行顺序查看反应；选择 **Union of found pathways** 可以对比所有已找到途径的并集。每个反应保留完整反应物、共生成物、计量系数和映射通道。并集用于比较，不代表其中任意组合都构成一条有效途径。途径视图不受网络深度显示过滤器限制，以保留完整反应；**Reset view** 返回完整网络。分析期间可以点击 **Stop**，已找到的途径会保留。**Save project** 保存最近的分析结果，**Export pathways JSON** 将分析结果与原网络一同导出。

```python
from chemweave import GenerationConfig, Generator, Species, find_pathways, load_rules

network = Generator(load_rules(), GenerationConfig(
    max_carbon=6, max_depth=4, stop_when_targets_found=False,
)).generate([Species.from_smiles("CCC")], targets=[Species.from_smiles("C=C")])

target_id = network.targets[0]["species_id"]
if target_id is not None:
    result = find_pathways(network, target_id, max_reactions=8, max_pathways=100)
    for route in result["pathways"]:
        print(route["reactions"])  # 按一个有效执行顺序排列的反应 ID
    print(result["termination"], result["network_search_complete"])
```

使用内置规则，这个例子会得到两条各含四个反应的乙烯生成途径。`examples/formation_pathways.toml` 可以通过命令行生成同一网络。分析保存的网络时，将下面的 `S00008` 替换为文件中目标产物的 ID：

```bash
python -m chemweave generate --config examples/formation_pathways.toml --output results
python -m chemweave pathways --network results/network.json --target S00008 --max-reactions 8 --max-pathways 100 --output results/pathways.json
```

`stop_when_targets_found = false` 放在 TOML 的 `[generation]` 段中；省略时保留原先的目标提前停止行为。`pathways` 命令在途径枚举触及数量上限或原网络提前停止时返回退出码 2，仍保存结果。枚举完成但没有匹配途径时返回 0，结果列表为空。

分析器先反向收集候选生成反应，再从种子出发正向枚举反应集合。只有全部反应物都有来源时，反应才可执行；仅当删除集合中的任意一个反应都会使目标不可达时，才保留该途径。这会排除多余制备支路和产物回收绕路，同时保留长短不同的替代途径。共享制备反应只计一次，共生成物可以供其他分支使用。如果图中的环同时产生了必需的共生成物，不会仅因结构成环而被删除；排除的是依赖尚无种子来源的物质来启动的循环。相同反应集合的不同独立分支执行顺序只算一条途径。若目标已经属于初始种子，返回一条零反应途径。

参数与结果的含义：

- `max_reactions` 统计**所有分支**中的不同反应。例如 `A → B`、`A → C`、`B + C → P` 包含三个反应，尽管依赖深度为二。网络记录的最小深度不代表每条替代途径的长度上限。
- 结果判断物种是否有来源，不计算有限投料、重复反应事件、净产率或动力学。例如自偶联要求该物种有来源，但不核算库存中是否还剩两个分子；计量系数仍保留供查看。
- 映射通道保留在原始反应记录中，不枚举通道组合，也不为整条途径重新拼接连续的 atom mapping。
- `termination="exhausted"` 表示在**当前保存的网络和途径反应数上限内**完成枚举。`max_pathways` 表示触及数量上限，未确认枚举完整；`user_stopped` 表示用户停止。搜索依照保存的反应顺序，不是最短途径或动力学排序。
- 只有原网络因 `exhausted` 或 `max_depth` 停止时，`network_search_complete` 才为真，且仅相对于其规则与配置范围。目标提前停止、物种/反应数量上限和用户停止均可能遗漏替代途径，即使途径枚举本身完成也一样。网页分别提示网络搜索和途径枚举的完整性。

## 命令行与输出

安装后可运行：

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

可使用 `examples/mapped_propane.toml` 从部分编号的种子生成同样的网络。

自定义配置中的 `rules = "my_rules.json"` 相对于配置文件所在目录解析。输出目录相对于运行命令的目录解析。

| 文件 | 内容 |
| --- | --- |
| `network.json` | 完整物种、反应、深度、目标状态、通道映射代表与计数、规则、配置及 RDKit 版本，可重新载入 |
| `species.csv` | 物种 ID、名称、不带映射的结构、分子式、电荷、相态、位点、带映射的 SMILES 及最小深度 |
| `reactions.csv` | 反应 ID、规则、反应族、步骤性质、两侧计量、模板匹配数、通道数及反应深度 |
| `channels.csv` | 每个对称等价通道的一份代表性 mapped reaction 及其 `match_count` |
| `stoichiometry.csv` | 每条反应每一侧的物种及正整数计量数；生成计量矩阵时反应物取负 |
| `network.dot` | Graphviz 图；每条反应有独立节点，保留多反应物和多产物关系 |

DOT 导出不依赖 Graphviz；如已安装 Graphviz，可用 `dot -Tsvg results/propane/network.dot -o results/propane/network.svg` 绘图。

保存结果中的 `termination` 表示停止原因：

| 值 | 含义 |
| --- | --- |
| `exhausted` | 所选规则和碳数范围内没有待处理组合；目标仍可能缺失 |
| `max_depth` | 达到深度边界，包含该深度的产物 |
| `targets_reached` | 所有目标均已出现，也可能已经包含在种子中 |
| `user_stopped` | 调用者或网页请求停止 |
| `max_species`、`max_reactions` | 数量预算不允许加入下一条完整反应，保存结构一致的部分网络 |

网络穷尽、达到深度边界或找到全部目标时，命令行在**没有缺失目标**的前提下返回退出码 0。数量预算中止、用户停止或仍缺失目标时返回 2，结果仍会保存。输入错误由 Python/RDKit 报告。

## Python 中的进度与用户停止

`Generator.generate` 新增可选的 `progress(network)` 和 `stop_requested()` 回调。完成种子登记和每次完整的模板/反应物组合执行后报告进度；进度回调访问当前网络，应尽快返回。停止回调返回真时，生成器返回结构一致的部分网络，并标记 `termination="user_stopped"`。已有调用无需修改。`Network.to_dict()`、`Network.from_dict()` 与文件持久化共用 JSON 表示。

## 后续扩展

先用带来源和适用条件的规则逐步补齐裂解、支化、环化和芳烃生成路径，再比较提高碳数上限后的变化。随后根据选定的表面机理加入空位和催化剂状态表示，最后接入 Single-event 参数组织和动力学求解。

## 参考

- [RDKit 化学反应接口](https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions)
- [RDKit Reaction SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)
- [Python package 打包指南](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

