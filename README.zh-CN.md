[English](README.md) | **简体中文**

# ChemWeave

给定初始物种和反应模板，ChemWeave 帮你探索反应网络如何展开。你可以在本地网页或 Python 中生成网络、查看原子映射、追踪产物生成路径并保存结果。它也能自动计算 **single-event kinetics 使用的反应路径因子 `n_e`**，支持气相反应和抽象催化位点模型。

## 目录

- [安装与测试](#installation)
- [网页界面](#web-interface)
- [Python 接口](#python-api)
- [自定义规则](#custom-templates)
- [网络如何扩展](#network-expansion)
- [表面状态与 Single-event 的边界](#catalyst-sites)
- [原子映射与对称性](#atom-mapping)
- [深度与目标物种停止条件](#search-limits)
- [目标产物生成途径](#formation-pathways)
- [命令行与输出](#cli-and-outputs)
- [Single-event kinetics：计算 n_e](#single-event-kinetics)
  - [自动计算 n_e](#automatic-ne)
  - [气相与催化模板共同运行的示例](#combined-example)
  - [独立催化示例：H₂ 在金属位点上解离吸附](#catalytic-example)
  - [原有对称数比验证示例](#symmetry-example)
- [Python 中的进度与用户停止](#progress-and-stop)
- [后续扩展](#development)
- [参考](#references)

<a id="installation"></a>

## 安装与测试

需要 Python 3.11 或更新版本，以及 RDKit。可以使用已有环境，例如 Conda 的 `base` 或 `ml`。在项目的上级目录运行：

```bash
conda activate base
cd chemweave
python -m pip install -e .
python -m pytest -q
```

pytest 作为可选开发依赖声明在 `.[test]` 中。需要一起安装时，可运行 `python -m pip install -e ".[test]"`；已有 pytest 的环境可以使用上面的命令。

在 macOS/Linux 中，已有 Python 和 RDKit 环境时，可以直接从项目源码运行，无需另行安装 ChemWeave：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/propane.toml --output results/propane
```

<a id="web-interface"></a>

## 网页界面

Windows、macOS 和 Linux 使用相同的 Python 启动入口。在包含 `pyproject.toml` 的项目目录打开终端，使用 Python 3.11 或更新版本的环境。首次安装网页可选依赖，随后启动本地工作台：

```bash
python -m pip install -e ".[web]"
python -m chemweave web
```

安装完成且对应环境已激活时，也可以使用等价命令 `chemweave web`。使用网页期间保持终端打开。浏览器会自动打开 `http://127.0.0.1:8765`；如果没有自动打开，可手动输入该地址。网页界面使用英文。输入初始物种 SMILES 和可选名称，勾选规则，设置碳数、深度后点击“Generate”（生成）。可添加目标物种来记录实验产物。“Stop when all targets are found”（找到所有目标后停止）默认开启；寻找替代生成途径时请关闭；展开各物种的“Phase and catalyst sites”（相态与催化位点）可以设置对应状态。

中间图使用分子结构节点和独立的反应节点，边上显示非 1 的计量系数。点击节点可以查看物种或反应详情、切换映射通道、显示该反应实际保存的 atom-map 编号。可以搜索、按深度筛选或查看选中节点附近。超过 80 个物种时，默认显示到第 2 层，可通过深度选项展开。筛选后只显示包含全部反应物和共生成物的完整反应。运行期间更新统计信息，结束后绘制网络。

“Stop”（停止）在候选反应之间请求停止；当前正在执行的 RDKit 底层匹配调用需要返回后，停止请求才能生效。已接受的反应记录保持完整，结果标记为 `termination="user_stopped"`，可以继续保存；匹配计数和替代路径可能尚未枚举完整。

“Save project”（保存项目）下载包含输入、所选规则、设置和网络的 JSON 文件。“Open”（打开）支持项目文件及现有 `network.json`，包含受支持的旧版格式。“Export”（导出）支持网络 JSON、CSV 表格及 DOT 的 ZIP 包，或当前网络视图的 PNG 图片。下载位置遵循浏览器的正常设置。

每个本地服务在内存中维护一个工作区，多个浏览器标签页共享它。刷新会恢复当前工作区；保存项目后可以跨服务重启恢复。**New** 在确认后创建空白工作区，不预填种子或模板。**Load example** 将当前工作区替换为气相＋催化示例输入，但不自动运行；**Generate** 只根据当前输入生成网络。服务只监听本机回环地址；可在终端按 Ctrl+C 关闭。

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

<a id="python-api"></a>

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

<a id="custom-templates"></a>

## 自定义规则

规则使用 JSON 列表存储。以下例子表示含 H2 的映射丙烷脱氢步骤。软件按用户提供的基元步骤处理，不负责证明其基元性：

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[C;H3:1](-[H:4])-[C;H2:2](-[H:5])-[C;H3:3]>>[C:1]=[C:2]-[C:3].[H:4]-[H:5]",
    "family": "dehydrogenation",
    "kind": "elementary",
    "version": "1",
    "source": "在此填写文献或规则依据",
    "notes": "映射步骤，未展开催化循环"
  }
]
```

也可以用 `ReactionRule(...)` 在 Python 中直接定义。主要字段如下：

| 字段 | 含义 |
| --- | --- |
| `id` | 规则唯一名称；不同机理即使总反应相同，也可以使用不同规则 ID |
| `smarts` | RDKit Reaction SMARTS，原子映射编号用于定义结构变化 |
| `family` | 反应族标签 |
| `kind` | 可选历史元数据，默认 `elementary`；所有输入步骤按基元反应处理，此字段不限制因子计算 |
| `version`、`source`、`notes` | 规则版本、依据和假设 |
| `exclude_smarts` | 任一反应物命中其中一个模式时，不应用此规则 |
| `reactant_phases`、`product_phases` | 按模板顺序指定各物种相态，默认全部为 `gas` |
| `reactant_sites`、`product_sites` | 按模板顺序指定旧式物种级位点标签，默认全部为 `null` |
| `site_constraints` | 按模板映射号约束显式催化节点属性；节点查询使用 `[#0]` |

默认不自动反转规则。逆反应需要显式提供模板；某些双键迁移模板自身可以生成正反两个方向。规则被声明为 `elementary` 也不等于程序已经证明了它的机理。

<a id="network-expansion"></a>

## 网络如何扩展

1. 标准化种子结构，以结构、相态、旧式位点标签及显式催化节点属性确定物种身份。
2. 每个新物种只与每条规则的各个反应物位置匹配一次。
3. 将新加入候选池的物种与已有候选物种组合，执行尚未处理的有序反应物组合。允许同一物种占据两个位置，以表示自偶联。
4. 用显式氢执行模板，恢复来源映射号，检查原子一一对应（包括元素和同位素）、电荷、旧式位点标签数，以及显式催化节点的守恒和容量。
5. 按 `max_carbon` 筛选，并按规则 ID 和两侧计量合并重复反应。每种非等价映射只在 `channels` 中保存一个完整代表；等价执行仅累加其 `match_count`，不保留重复记录。
6. 把新物种加入队列，直到队列耗尽或达到物种/反应预算。

检查放在输入读取和候选反应生成处；正常枚举不重复校验整张网络。分子描述和反应模板在创建时处理一次。当前枚举不估算能垒或通量，也不会仅因某个物种未在实验中检出就删除它。

`max_carbon` 对每个物种生效，种子也必须处于这个范围内。它定义本次网络的范围，不表示范围外路径在实际体系中不重要。

<a id="catalyst-sites"></a>

## 表面状态与 Single-event 的边界

物种可以包含显式映射的催化节点。物种 SMILES 中的 `[*:1]` 表示抽象位点；反应 SMARTS 必须用 `[#0:1]`，避免与匹配任意原子的通配符混淆。节点属性单独提供：

```python
metal = Species.from_smiles("[*:1]", phase="surface",
    site_atoms={1: {"kind": "metal", "type": "Fe_model", "capacity": 1}})
acid = Species.from_smiles("[*:1]O", phase="surface",
    site_atoms={1: {"kind": "bronsted", "type": "ZSM5_model", "capacity": 1}})
```

抽象节点可为空位，或通过键级 1–4 的一条键连接一个原子；省略属性时默认为 kind `surface`、type `X`、capacity 1。多个节点表示多位点占用。Brønsted 酸模型包含真实 O–H 键和可转移的映射氢；质子转移后，框架 O− 与碳正离子作为一个表面关联物种保存。属性、吸附连接、同位素和电荷参与身份、匹配、映射、守恒与保存。晶格邻接和吸附取向不在这一抽象图模型的范围内。详细约定及文档案例见[计算与验证报告](docs/path_degeneracy.md)。

旧式 `Species.from_smiles("C=CC", phase="surface", site="BAS")` 仍可用于通用带标签网络，但物种级标签不足以自动计算催化路径因子。网页展开 **Phase and catalyst sites** 后填写节点属性 JSON。保存时需要 JSON 或配套 `site_atoms`；仅靠紧凑 SMILES 不能恢复催化位点类型。

<a id="atom-mapping"></a>

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

模板在显式氢分子上执行。涉及氢转移的规则必须明确移动带映射的氢原子及其化学键。若生成无法追溯的新原子、丢失来源原子或在产物中重复使用编号，候选会以 `atom_mapping` 原因被拒绝；电荷和位点守恒另行检查。旧的隐式氢模板可能需要改写，本版本不会自动猜测氢的对应关系。内置第 2 版规则提供完整原子对应。与用户模板一样，基元性是输入假设，不是软件验证结果。

`reaction.channels` 对每种非等价反应映射只保存一份完整代表。记录包括 `mapped_reaction_smiles`、两侧带映射的 SMILES 列表、`atom_mapping` 来源编号表、产物原子来源、`symmetry_key` 和 `match_count`。即使计量关系合并了两份相同物种，记录仍区分反应物实例顺序。遇到重复时，仅增加已有通道的 `match_count`，丢弃该次完整映射；不再单独保存 `reaction.matches` 或重复匹配索引列表。

例如，丙烷脱氢的 24 个模板匹配只保存为一个通道，其 `match_count=24`。摘要中的 `template_matches` 和 CSV 同名列为这些计数的总和，不是实际存储的记录数。

通道通过反应前后完整图的标准化表示归并。图中的分子节点保留分组、相态和位点；原子标签保留元素、同位素、电荷、自由基及已指定的原子 CIP 标签；化学键标签保留两侧键级和已定义的立体信息。这样不会仅凭原子各自属于同一对称类别，就认定多位点反应等价；化学上相同的反应物实例则允许整体交换。

基础原子类别通过 `CanonicalRankAtoms(breakTies=False, includeAtomMaps=False)` 计算，不枚举分子自匹配。显式氢新戊烷的类别大小为 1、4、12。真实映射号始终唯一，不会被对称类别标签覆盖。已定义的立体化学会保留；程序不穷举新生成的立体异构体。

这些是单步反应映射。编号属于物种代表或某次反应实例，并非全网唯一的原子身份。多步来源追踪需要连接各步映射，并将产物原子与物种库代表的原子对应起来，不能仅凭编号数值直接连接。通道数及 `match_count` 都不是 SEMK 简并度；独立因子和动力学参数共享由可选的 single-event 标注模块处理。本版本不拟合参数或求解速率、转化率和选择性。

新的 JSON 导出使用第 6 版格式，保留精简后的通道代表、计数、编号、对称类别、种子输入、物种及反应深度和目标状态。第 3、4、5 版文件仍可读取；旧文件缺失的深度保持未知。读取第 2 版文件时，会将原始匹配和通道索引转换为精简表示。第 1 版文件仍可作为结构网络读取，其中不完整的旧来源记录不会保留为映射通道；不能从旧的不带编号结构中恢复来源映射。

<a id="search-limits"></a>

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

<a id="formation-pathways"></a>

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

<a id="cli-and-outputs"></a>

## 命令行与输出

安装后可运行：

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

可使用 `examples/mapped_propane.toml` 从部分编号的种子生成同样的网络。

自定义配置中的 `rules = "my_rules.json"` 相对于配置文件所在目录解析。输出目录相对于运行命令的目录解析。

| 文件 | 内容 |
| --- | --- |
| `network.json` | 完整物种、反应、深度、目标状态、通道映射代表与计数、可选的 single-event 标注与共享参数、规则、配置及 RDKit 版本，可重新载入 |
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

<a id="single-event-kinetics"></a>

## Single-event kinetics：计算 n_e

SEMK 是 ChemWeave 的一种用途。`Generator` 仍然生成通用的规则反应网络，不带动力学假设。可以在生成后调用 `annotate_single_events(network, ...)`，在 CLI 配置中加入 `[single_event]`，或在英文网页的 **Single-event annotation → Annotate after generation** 中启用。省略该选项时网络不带标注。所有输入的映射步骤（包括默认五条示例）均可进入路径因子计算。旧项目中的 `kind="overall"` 历史标记不再阻止计算；对旧项目重新应用标注即可替换原来的跳过结果，无需重新生成网络。

`single_event.py` 是与反应族无关的框架，只负责中心提取、调用显式选择的计算器、接受导入因子及管理共享参数。默认不选择计算器：`annotate_single_events(network)` 只提取反应中心，不假定任何反应族或赋予因子。自动结构路径计数位于 `path_degeneracy.py`，可选示例模板位于 `data/path_families.json`，计数算法不读取这些模板。原有甲烷对称数示例仍单独保留在 `single_event_calculators.py`。

替代计算方法可通过 `register_single_event_calculator(name, callback)` 注册，再用 `annotate_single_events(network, calculator=name)` 选择，无需修改核心模块。回调接收关键字参数 `network`、`reaction`、`channel`、`center`；不适用时返回 `None`，适用时返回字典，包含 `convention`、`source`，以及 `factor` 或 `sigma_reactants`/`sigma_ts`。可选的 `references`、`kinetic_subclass`、`parameter_class` 也会保留，但不是因子计算前提。具体化学条件及相态/位点检查属于计算器；通用层只检查和保存统一的因子/参数记录。导入记录优先于相应通道的计算结果。在新进程中重放自定义计算器的保存配置前，需要重新注册回调；可执行 Python 代码不会被序列化。

<a id="automatic-ne"></a>

### 自动计算 n_e

选择 `calculator="mapped_event_v1"`（`path_degeneracy_v1` 保留为兼容名称）。计算器直接从完整映射基元反应提取断成键、原子状态变化和产物物种划分，在反应物图的对称变换下枚举等价事件。新增模板不需要第二套 recipe，也不需要内置反应族名称。保存方法为 **`reaction_path_degeneracy`**，算法为 **`mapped_event_orbit_v1`**；它是明确计数约定下的结构路径因子，不是三维转动对称数计算。

```python
from chemweave import GenerationConfig, Generator, ReactionRule, Species, annotate_single_events

rule = ReactionRule(id="user_CH_homolysis",
                    smarts="[C;X4:1]-[#1:2]>>[C;H0:1].[#1;H0:2]")
network = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("C")])
annotate_single_events(network, calculator="mapped_event_v1")
assert network.reactions[0].channels[0]["single_event"]["n_e"] == 4
```

纯函数 `calculate_path_factor(mapped_reactants, mapped_products, ...)` 接受完整映射、显式氢的物种 SMILES 列表。每个列表元素代表一个物种实例，包括作为整体保存的不连通关联物种。网络生成器会根据映射模板自动提供这些完整映射；可选相态/位点状态及显式节点属性用于保留催化环境差异。参数类别、共享参数、`rate_basis`、`A_units` 均不是前提；纯因子计算不创建参数类别或速率表达式。

**计数约定：**枚举时保留不同反应物实例的身份；同一组编号原子上重复生成的事件先去重，得到 `N = unique_labeled_paths`，此时尚未商去相同反应物的交换。令 `D = product(m_i!)`，其中 `m_i` 是相同反应物物种的实例数，只修正一次：**`n_e = N / D`**。例如 H2 + 2 X → 2 XH 为 `2 / 2 = 1`，逆反应为 `1 / 2 = 0.5`。`match_count`、产物相同和已保存通道数均不能代替因子。保存的中间计数明确记录这一次修正。使用者负责在自己的动力学表达式中一致使用该因子，避免再次施加相同反应物除数。

运行未列入内置反应族的新模板示例：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/custom_ne.toml --output results/custom-ne
```

网页展开 **Custom templates (JSON)**，输入模板列表，点击 **Use these templates**；填写种子 SMILES（必要时指定相态/位点属性），勾选 **Annotate after generation**，保留配置 `{"calculator":"mapped_event_v1"}`，再点击 **Generate**。点击反应即可查看因子、状态、事件变化及路径数/除数。**Apply to current network** 也可标注已有映射网络。气相/金属/酸位示例只是快捷入口，不限制用户输入。

算法不设 C/H、无环或反应族白名单。验证覆盖均裂、夺氢、歧化/复合、金属位点吸附/脱附/氢化，以及酸位质子转移、氢负离子迁移和 β 断裂，另检查含氧及环状输入。计算对象是输入的电子结构图，不自动枚举共振结构集合。已指定立体化学目前返回 `unsupported_algorithm`；映射不完整或缺少催化节点返回 `missing_information`；输入矛盾返回 `invalid_input`；计算超限返回 `incomplete`。这些状态均附原因且不填因子，仍可手动导入。抽象位点不隐含晶格、吸附取向或三维简并度。

[实现与验证报告](docs/path_degeneracy.md)列出计数依据、11 个文档对应或明确标记为派生的案例、人工枚举、独立穷举交叉检查及固定 RMG 源码引用。**本地没有运行 RMG**。本路线不采用 2023 年加氢裂化论文的自由态抵消假设，也不声称与 Vandewiele 转动对称数约定普遍相等：甲基—乙基歧化的路径因子为 3，而该文对称数约定为 6。

<a id="combined-example"></a>

### 气相与催化模板共同运行的示例

点击 **Load example** 才会同时载入五条模板，种子为丙烷和金属空位，此演示默认启用可选因子标注。新建工作区保持空白。丙烷脱氢在深度 1 产生 H₂，随后在深度 2 参与 H₂ + 2 M → 2 MH。催化吸附模板与气相模板一起勾选，也可单独取消。旧项目保留其原有模板和种子；点击 **Load example** 载入组合示例，再点击 **Generate** 运行。

`examples/combined_gas_catalytic.toml` 可运行同一组合：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/combined_gas_catalytic.toml --output results/combined-gas-catalytic
```

这些示例步骤通过生成的 H₂ 连接；这里没有增加金属参与丙烷脱氢的机理。

<a id="catalytic-example"></a>

### 独立催化示例：H₂ 在金属位点上解离吸附

与自定义气相模板并列，`examples/custom_catalytic_rules.json` 提供自定义催化模板，不注册反应族，也不另写 recipe；`examples/custom_catalytic_ne.toml` 配置氢气和显式金属空位：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/custom_catalytic_ne.toml --output results/custom-catalytic-ne
```

来源：Lozano-Blanco 等（2011），p.427 Table 2 (1)，H₂ + 2 M → 2 MH。节点 `[*:1]` 的相态为 `surface`，属性为 `kind="metal", type="M", capacity=1`。只需输入一种空位物种，生成器会将其用于两个独立、相同的反应物实例。模板以 `[#0]` 匹配位点，断开 H1–H2，形成 H1–M3 和 H2–M4。不同编号配对只有 `(H1–M3, H2–M4)`、`(H1–M4, H2–M3)` 两种；同时交换氢和位点只会重复已有配对，不增加路径。因此 `N=2`，`D=2!`，**n_e=1**。位点及其属性在吸附物中保留，不附加晶格邻接或取向倍数。

英文网页点击 **Single-event annotation → H₂ adsorption on metal → Generate**，可以载入相同反应的快捷示例。使用独立自定义示例时，导入上述规则 JSON 并填写种子，或直接 **Open** 命令行生成的 `results/custom-catalytic-ne/network.json`。自定义模板使用独立 ID，未填写反应族名称，仍通过通用算法计算。反应详情和 CSV 均记录两条编号路径及除数 2。

<a id="symmetry-example"></a>

### 原有对称数比验证示例

原有对称数验证族是正向气相 **H/D 原子 + 甲烷(H/D) → 双氢分子 + 甲基自由基**。基元规则单独放在 `src/chemweave/data/single_event_rules.json`；网页的 **Load methane validation example** 可以载入该例，不改变默认规则。对应 CLI 示例：

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/single_event.toml --output results/single-event
```

计算器 `methane_h_abstraction_v1` 采用四面体甲烷反应物和沿 C–H/D–H/D 轴的过渡态鞍点，另有三个旁观配体。**以正旋转对称数定义** `n_SE = product(sigma_reactants) / sigma_TS`，不采用完整图自同构数。原子反应物的对称数为 1。该约定不包含产物对称数、电子/核自旋权重、隧穿和再穿越，也不从 atom mapping 推断这些效应；该验证族没有内转子。

| 甲烷 | 转移同位素 | 甲烷 sigma | TS sigma | n_SE |
| --- | --- | --- | --- | --- |
| CH4 | H | 12 | 3 | 4 |
| CH3D | H / D | 3 | 1 / 3 | 3 / 1 |
| CH2D2 | H / D | 2 | 1 / 1 | 2 / 2 |
| CHD3 | H / D | 3 | 3 / 1 | 1 / 3 |
| CD4 | D | 12 | 3 | 4 |

这些 H/D 结果来自上述理想几何及保持同位素的正旋转：Td 有 12 个、C3v 有 3 个、C2v 有 2 个、Cs 有 1 个。[NIST 对称数说明](https://cccbdb.nist.gov/thermox.asp) 介绍了相应约定；[NIST CCCBDB 第 8 号反应](https://cccbdb.nist.gov/ts2x.asp?ri=8) 收录 H + CH4 基元反应；[Hu 等，J. Phys. Chem. A 2006](https://web.stanford.edu/group/Zarelab/publinks/764.pdf) 描述了其 C3v 鞍点。同位素变体验证的是这一对称数约定，不是实验测得的动力学同位素效应。包中未提供拟合或从头算速率常数。

每个映射通道新增 `single_event.reaction_center`：变化的化学键（含键级/立体信息）、变化的原子电荷/自由基/立体状态，以及参与变化的**本次反应局部映射号**。甲烷对称数计算器检查完整的断成键/自由基变化及气相无位点条件，不依赖规则名称。请统一普通氢的标签：甲烷混用未指定同位素 H 与显式 1H 时不计算因子，因为 RDKit 会将其拆成不同通道。动力学子类包括碳同位素、转移同位素、夺氢原子同位素和旁观氢同位素组成，因此不同同位素环境不会自动共享拟合参数。不在旧甲烷计算器适用范围内的输入仍保留中心，并明确标为 `unsupported`，不赋因子；历史步骤分类标签不限制通用路径计算器。没有映射通道的反应列在 `network.single_event["unmapped_reactions"]` 中。

**`match_count` 和通道数均不参与 single-event 因子的计算。** 测试故意改变模板，使匹配次数从 4 变为 12，独立计算的因子仍为 4；提前停止导致的计数不完整也不影响对称数比。这一旧对称数计算器只覆盖正向甲烷族；独立路径计算器的范围和归一化见上文，不能把它理解为对称数比算法的通用扩展。

```python
from chemweave import Network, annotate_single_events

network = Network.from_json("results/single-event/network.json")
annotate_single_events(network, calculator="methane_h_abstraction_v1")
channel = network.reactions[0].channels[0]
print(channel["single_event"]["factor"])  # 4.0，来自 12 / 3
class_id = channel["single_event"]["parameter_class"]
print(network.single_event["parameters"][class_id])
# A、n、Ea 初始为 None；status 为 "unfitted"。
network.export_json("results/single-event/annotated.json")
network.export_csv("results/single-event")
```

旧计算器和显式导入的参数记录仍保留可选共享参数支持：通道可以引用 `parameter_class`，参数仅在 `network.single_event["parameters"]` 中保存一份。只有存在参数类别时才保留原有 modified Arrhenius 表达式。通用路径计算器不创建这些记录，也不要求速率表达式、单位或拟合值；本模块不执行拟合或速率求值。

也可导入有限正因子，或通过反应物/TS 对称数比指定因子。导入需提供 `convention` 和 `source`；子类、共享参数类、`rate_basis`、`A_units` 均可省略。如果同时提供因子和对称数，两者必须一致。以下旧甲烷输入展示可选共享参数记录，**不是实验拟合数据**：

```json
{
  "calculator": null,
  "factors": [{
    "reaction_id": "R00001", "channel": 1,
    "kinetic_subclass": "methane-H-abstraction",
    "parameter_class": "methane-H-fit",
    "sigma_reactants": [12, 1], "sigma_ts": 3,
    "convention": "Forward axial saddle; proper rotations; n_SE = 12/3",
    "source": "Independent symmetry assignment under the documented methane convention"
  }],
  "parameters": {
    "methane-H-fit": {"A": null, "n": null, "Ea": null, "status": "unfitted", "source": "Awaiting experimental fitting"}
  }
}
```

保存为 `single-event.json` 后，可对已有网络标注：

```bash
python -m chemweave annotate --network results/single-event/network.json \
  --annotations single-event.json --output results/annotated
```

这些键同样适用于 `annotate_single_events(network, **options)` 和网页 **Calculation options and imported factors (JSON)** 编辑区。**Apply to current network** 可以更新已载入的网络，无需重新生成。通道索引从 1 开始，只指向该保存网络中的记录；重新生成后应重新核对选择器。同一参数类不能混用不同分子数、速率基准或单位。原有甲烷计算器分配参数类别，通用路径计算器不分配。每次应用配置会替换旧标注和参数记录，因此需要保留的导入因子及参数值也应包含在配置中。`calculator: null`（省略时也是这一默认行为）关闭自动因子计算，仍提取反应中心。网页验证按钮会显式填入甲烷计算器；空配置对象不选择反应族。`annotate` 命令未提供配置文件时也只提取中心。

网络 JSON 第 6 版及网页项目文件保留中心、子类/共享参数类引用、独立因子、约定、来源、共享参数记录和可复现的标注配置。第 1–5 版仍可读取；缺少完整映射的旧记录在获得可用映射前不会被推断赋因子。仅在应用标注后，CSV 导出增加 **`single_event.csv`**（每个映射通道一条，包括不支持状态及中心）与 **`kinetic_parameters.csv`**（每个共享参数类一条）；`species.csv` 追加 `site_atoms`，标注导出追加 `n_e`、方法、算法、反应族、状态/原因、提取的事件、动力学通道、编号路径数及相同反应物除数，原有列保留。英文网页的反应详情可切换通道，查看中心、独立因子、对称数/过渡态约定、动力学子类、共享参数类以及已提供或待拟合的参数值。Save/Open 与 JSON/CSV 下载保留相同数据。

<a id="progress-and-stop"></a>

## Python 中的进度与用户停止

`Generator.generate` 新增可选的 `progress(network)` 和 `stop_requested()` 回调。完成种子登记和每次完整的模板/反应物组合执行后报告进度；进度回调访问当前网络，应尽快返回。停止回调返回真时，生成器返回结构一致的部分网络，并标记 `termination="user_stopped"`。已有调用无需修改。`Network.to_dict()`、`Network.from_dict()` 与文件持久化共用 JSON 表示。

<a id="development"></a>

## 后续扩展

后续可改进模板编写、不同规则集合的比较、物种生成路径探索，以及大型网络的生成性能。可选路径因子算法可继续扩展立体化学与共振处理。实验拟合和动力学求解可作为独立应用使用导出的网络。

<a id="references"></a>

## 参考

- [RDKit 化学反应接口](https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions)
- [RDKit Reaction SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)
- [Python package 打包指南](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

