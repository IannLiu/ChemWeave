> **2026-11-09 更新：下文为历史交接记录，不代表当前实现边界。** 现已实现无需额外 recipe 或反应族白名单的 `mapped_event_v1`，参数类别及单位不再是纯因子计算前提。当前接口、计数约定、验证和限制以 [docs/path_degeneracy.md](docs/path_degeneracy.md) 及中英文 README 为准。

# Single-event 自动计算研究与实施交接

> 更新：下文保留的是此前研究阶段的历史快照。用户随后选择了 RMG 风格反应路径简并度路线；现已实现 `path_degeneracy_v1`、显式金属/Brønsted 位点和 schema 6。当前能力、文献来源、计数约定及验证结果以 [docs/path_degeneracy.md](docs/path_degeneracy.md) 和中英文 README 为准。旧文中的“尚未实现”“等待建议”以及工作区/测试数量不再代表当前状态。后续已取得并阅读本地 2023 年加氢裂化正文及 2011 年 Fischer–Tropsch 机理论文。

更新日期：2026-11-08（Asia/Shanghai）。项目：`/Users/liuyan/Desktop/aromatization/chemweave`。

## 1. 下一轮如何继续

用户准备先提出建议，再开始相应实现。请先阅读本文及当前代码，听取并整合用户建议，再确定第一批修改。本文记录的是研究结论和建议方案，不代表用户已逐项批准全部架构、依赖或催化假设。

此前研究阶段明确要求不修改代码；本次仅新增本交接文档。通用自动对称性引擎、通用过渡态构造和催化环境模型尚未实现。不要将建议方案描述为现有能力。

用户要求：

- 保持 ChemWeave 为通用规则网络生成器；SEMK 只是可选用途之一。
- 复用 atom mapping 提取反应中心，关联动力学子类，独立计算或导入 single-event 因子，关联按类别共享、以后由实验拟合的参数。
- **绝不能把 `match_count` 或当前通道数当作 single-event 数 `n_e`。**
- 通用核心不能包含针对某个具体分子的查表或条件分支；反应族专用的化学与 TS 规则应与核心分离。
- 目标包括催化反应。不能把气相结果直接用于表面；必须区分通用对称性算法与催化环境模型。
- 支持保存、导出和英文网页；同步中英文 README。
- 避免过度防御、过度拆分大量小函数。可用本地 base 或 ml 环境测试。
- 本轮没有要求提交 Git，也没有授权覆盖现有未提交工作。

新对话可用的开场说明：

> 请先阅读 `/Users/liuyan/Desktop/aromatization/chemweave/SINGLE_EVENT_HANDOFF.md`，核对当前代码状态。我们要继续开发可选的自动 single-event 数计算模块，目标包括催化反应。我会先提出建议，请结合建议确定本轮范围后再实施。保留现有未提交修改；不要把 `match_count` 当作 `n_e`，不要默认自由态结果适用于催化环境。

## 2. 工作区快照与既有实现

交接时 HEAD：`664e7d8 add single-evnet module`。

新增本文之前，`git status --short` 为：

```text
 M README.md
 M README.zh-CN.md
 M src/chemweave/__init__.py
 M src/chemweave/single_event.py
 M src/chemweave/static/index.html
 M src/chemweave/web.py
 M tests/test_single_event.py
 M tests/test_web.py
?? src/chemweave/single_event_calculators.py
```

这些修改是先前将甲烷专用逻辑移出通用核心的工作，必须保留。下一轮重新检查状态，因为用户可能继续修改文件。

现有实现：

| 文件 | 现有职责 |
| --- | --- |
| `src/chemweave/generator.py` | 通用网络枚举、物种和反应去重、映射通道归并、停止条件 |
| `src/chemweave/rules.py` | 执行显式 H 的 RDKit 反应模板，恢复来源映射，检查原子、元素/同位素、总电荷和 site 标签数守恒 |
| `src/chemweave/mapping.py` | 每次反应的唯一局部编号；构造前后联合图的 canonical `symmetry_key`；不计算自同构群阶 |
| `src/chemweave/models.py` | Species/Reaction/Network；JSON schema 5；兼容旧 schema 1–4；CSV 等导出 |
| `src/chemweave/single_event.py` | 通用反应中心提取、计算器注册/调度、因子导入、共享参数登记 |
| `src/chemweave/single_event_calculators.py` | 目前唯一生产计算器：气相 H/D 原子从甲烷 H/D 同位素体夺氢的专用验证逻辑 |
| `src/chemweave/__init__.py` | 注册 `methane_h_abstraction_v1`；注册不等于默认选择 |
| `src/chemweave/data/single_event_rules.json` | 独立验证规则，未改变默认通用演示规则 |
| `examples/single_event.toml` | 显式选择甲烷计算器的示例 |
| `src/chemweave/cli.py`、`web.py`、`static/*` | 可选标注、保存打开、JSON/CSV 导出、英文界面展示 |
| `tests/test_single_event.py` | 甲烷 H/D、映射重编号、与匹配次数无关、导入/共享参数/往返保存、适用范围等测试 |

核心调用关系：网络生成和 mapping → 已生成的网络 → 可选 single-event 标注。标注不影响网络生成规则。

当前 API：

```python
annotate_single_events(network, calculator=None, factors=(), parameters=None)
register_single_event_calculator(name, callback)
```

- 未选择计算器时仅提取反应中心，不自动赋因子。
- callback 接收 `network, reaction, channel, center` 关键字参数；不适用时返回 `None`。
- 返回记录包含 `kinetic_subclass, parameter_class, convention, source`，以及独立 `factor` 或 `sigma_reactants/sigma_ts`。
- 当前仅保存规定的字段及可选 `references`；新增审计数据需要修改核心，不能假设任意 callback 字段会自动保留。
- 已有因子状态为 `calculated/imported/unsupported`。参数 A/n/Ea 默认 `None`，并不拟合或求解动力学。
- 调用标注会替换此前标注和参数登记。保存计算器名称不等于保存 Python 可执行代码；自定义计算器重放时需重新注册。
- `single_event.csv` 和 `kinetic_parameters.csv` 仅在有标注时添加。

甲烷计算器明确限制气相且无 site，排除反向反应和甲烷中混用未指定同位素 H 与显式 1H 的情况。它使用已知对称数约定，不是通用 SIGMA 实现。

## 3. 已取得并阅读的文献

### Vandewiele 正文与 SI

- 正文：`/Users/liuyan/Desktop/aromatization/J Comput Chem - 2014 - Vandewiele - Symmetry calculation for molecules and transition states.pdf`
- SI：`/Users/liuyan/Downloads/jcc23788-sup-0001-suppinfo.doc`
- 文献：Vandewiele et al., *Symmetry calculation for molecules and transition states*, J. Comput. Chem. **2015**, 36, 181–192，DOI `10.1002/jcc.23788`。2014 为在线发表年。
- 正文和 SI 都已读取。此前出版社下载失败的阻塞已由用户提供本地 SI 解除。
- 本文页码均为正文印刷页码；PDF 页码 = 印刷页码 − 180。

正文主要核实结果：

1. p.182 式（2）：`sigma = sigma_ext * product(sigma_int_i)`。输出包括内部转动贡献，不是单一刚性构象的外部点群对称数。
2. p.183 式（7）：`sigma = |Aut(G)| / n_d`。图节点/边属性区分原子与键；自同构群阶是起点，不是最终物理对称数。
3. pp.183–184 式（8）、（9）：`n_d = 2**n_label`，`n_label = n_pot - n_real`。这是该文采用的 label-stereocenter 计数方法，存在适用边界，不是任意结构的普适公式。
4. label-stereoisomers 指把相同取代基人为区分后形成的可区分构型；**不能把现有唯一 atom map 编号直接作为对称性图颜色**。
5. p.185 表 3：甲烷 sigma=12，氘代甲烷=3，乙烷=18，新戊烷=972。后两者包含内部转动。
6. pp.184、186：孤对电子与自由基采取不同几何约定。未配对电子不直接作为四面体中心的第四个取代基；部分自由基采用最对称构型。
7. pp.186–187：需要处理共振杂化和芳香键等价，不能任取一个共振式导致对称性丢失。
8. p.187：部分柔性环结果可能是时间平均意义上的上限，而非最低能构象值。
9. p.189 Scheme 11：立方烷的耦合 label 中心及 SF6 等特殊几何暴露了简单计数的限制。正文 SF6 段落同时列出 720、预期 sigma=24 和 n_d=20，按式（7）前两者应对应 n_d=30；属于需保留说明的数字不一致，不能照抄为真值。

### SI 数据审计（已实际执行）

SI 仅含 `Section S1: Validation data set`，列为 `InChI [expected, calculated value]`。没有额外算法源代码、TS 图集、反应 n_e 表、催化剂/位点/受限转动模型。

- 594 条记录，594 个不同 InChI。
- 527 条 expected=calculated，占 88.7205%。
- 67 条不一致。
- `max(expected, calculated)/min(expected, calculated) <= 2` 有 570 条，占 95.9596%。
- 13 条 InChI 含 `/i` 同位素层；2 条含 `/t` 四面体立体层；23 条含 `/b` 双键立体层。这些只是编码字段统计，不等于覆盖了所有相关化学情况。
- base 的 RDKit 2026.03.6 用 `Chem.MolFromInchi(..., removeHs=False)` 读取时，590 条返回结构，4 条返回 None：`InChI=1S/CHCl/c1-2/h1H`、`InChI=1S/CHBr/c1-2/h1H`、`InChI=1S/CCl2/c2-1-3`、`InChI=1S/CClF/c2-1-3`。
- 成功解析不代表电子态、共振、立体信息和输入约定已经全部正确恢复。
- 上述为 SI 数据及导入审计，**不是 ChemWeave 新算法的测试成绩**。

重要的 SI 数值：

| 物种 | expected | calculated |
| --- | ---: | ---: |
| 甲烷 | 12 | 12 |
| 乙烷 | 18 | 18 |
| 丙烷 | 18 | 18 |
| 异丁烷 | 81 | 81 |
| 新戊烷 | 972 | 972 |
| 异丙基自由基 | 9 | 18 |
| 叔丁基自由基 | 81 | 162 |
| S8 | 8 | 16 |
| 立方烷 | 24 | 1 |
| SF6 | 24 | 360 |

应同时保留 expected 和 calculated。不能把 SIGMA 输出全部视为正确答案，也不能在未分析约定前假定所有 expected 都能直接用于我们的模型。尤其正文展示叔丁基自由基 calculated=162，而 SI 的 expected=81，说明自由基几何约定需要明确。

### 2023 年加氢裂化论文

Wang et al., *Single event kinetic modeling for paraffin hydrocracking over an industrial Ni–W silica–aluminum catalyst*, React. Chem. Eng. 2023, 8, 123–136，DOI `10.1039/D2RE00286H`。

来源：<https://pubs.rsc.org/en/content/articlehtml/2023/re/d2re00286h>。相关 §2.6 正文已通过网页搜索返回内容核实；没有保存本地全文 PDF。

- §2.6 末段明确承认对称数按自由态计算；作者假设反应物和过渡态结合活性位时损失相同自由度，因此 n_e 不变。
- 这是特定动力学模型采用的抵消假设，不是对任意表面反应成立的定理或验证。
- Figures 7–8 给出碳正离子 β-裂解的部分键表示（形成双键的边为 1.5，断裂边为 0.5）；它仍不等于显式表面 TS 模型。
- 不能因为两篇文献都使用“global symmetry”名称，就忽略物种/对映体合并和统计约定的差异。

## 4. 从 sigma 到 n_e 的边界

Vandewiele p.182 式（5）：`n_e = product(sigma_reactants) / sigma_TS`，配合 `k = n_e * k_single_event`。

该式要求已确定过渡态及通道统计约定。论文没有给出由任意 atom-mapped 反应自动恢复所有正确 TS 的完整算法。

- pp.184–185 Scheme 4：H 加成丙烯的 TS 中，变化中的 C1–C2 键仍具有刚性约束；不能把过渡键一律当自由单键。
- p.187 Scheme 8：H 原子从甲烷夺氢，含一条断裂键和一条形成键。
- pp.187–188 Scheme 9：甲基加成丙烯存在两个不可叠合、等能且通向同一产物的 TS，需要在单一 TS 的对称比之外合并贡献。
- Scheme 10：若不同 TS 通向分别保存的不同对映产物，应分别处理，不能再盲目乘 2。
- 作者明确将这些附加修正的自动分配留在本文范围之外。

建议区分：单个 TS 的 `symmetry_ratio`、TS 集合/贡献、最终动力学通道的 `n_e`。只有路径可以共享同一单事件速率系数时才直接合并贡献；不等能路径分别处理。不要重复计算既显式枚举又作为倍数纳入的路径。也不要无依据将 n_e 四舍五入为整数。

## 5. 催化环境结论

Vandewiele 正文和 SI 没有显式催化剂结构、吸附锚点、位点占用/邻位关系、表面法向或孔道约束模型，也没有一般受阻转子势能与配分函数计算。

可复用：带属性图、自同构计算、立体修正思想、可审计结果结构。

必须补充：环境允许的置换/运动、吸附连接、取向、位点等价关系、统计基准，以及反应族表面 TS。给分子加 `phase="surface"`、`site="BAS"` 不等于已经具备这些信息。

“受阻转动”不自动意味着 sigma=1；若势能仍保留等价周期极小值，相关对称性可能保留。要区分对称计数与转动势能对熵/速率的影响。

用于检查自由态抵消假设的分析表达式（不是文献原编号公式）：

```text
f_R  = sigma_R_cat  / sigma_R_free
f_TS = sigma_TS_cat / sigma_TS_free
n_e_cat = n_e_free * f_R / f_TS
```

前提是相同通道、位点归一化和其他统计约定。`sigma_R` 指整个反应物状态的对称因子，独立自由分子时才直接写成分子对称数乘积。沿用自由态值要求 `f_R == f_TS`；自由度数量损失相同，不足以证明对称因子按同一倍数变化。还须检查等能路径数量是否改变。

当前判断：

- 气相 H 原子夺氢验证不能证明表面 C–H 活化或表面 H 反应也可沿用同一 n_e。
- 酸位碳正离子 β-裂解可作为催化扩展的首个候选，但抵消假设尚未针对用户体系验证。
- 用户实际催化剂和位点结构尚未选定，不能默认等同于 2023 年论文中的催化体系。
- 如果采用抵消近似，应显式标为 `Adsorption cancellation assumed`，与显式环境计算结果分开。
- 位点数量、覆盖度和吸附构型计数不得无说明地塞进分子 n_e 或重复计入速率式。

## 6. 现有代码的关键差距

1. `reaction_symmetry_key` 是前后联合图的 canonical hash，不是自同构群阶、物理对称数或 n_e。
2. `Species.symmetry_classes` 是等价原子类别，不含所需的完整群和物理修正。
3. `reaction_center` 提供成断键和原子状态变化，但不提供 TS 几何、刚性、转子、吸附连接或 TS 集合。
4. 当前中心记录 CIP 标签变化；CIP 改变可能只因取代基优先级改变，不必然表示构型翻转。后续需比较映射邻接顺序下的构型关系。
5. 当前 `single_event.py` 强制 `factor == product(sigma_R)/sigma_TS`。若支持额外 TS 合并贡献，需要扩展数据与检查，不能通过篡改 sigma_TS 隐藏修正。
6. 当前一个 mapping channel 不等于唯一动力学 TS。建议在其下记录多个 TS/动力学贡献，保留来源关系。
7. 吸附构型若应是不同物种，需在物种去重之前进入身份定义；标注后无法恢复已经合并的物理状态。通用状态表示可以扩展，生成器无需理解 SEMK 公式。
8. 当前 site 仅为抽象单标签，无完整空位、多位占用和催化状态表示。
9. 当前 A 的默认单位由分子数推出气相浓度单位；催化速率的覆盖度、位点或催化剂质量基准需要显式定义。
10. 动力学参数共享应考虑环境和同位素适用性，不应仅凭有机反应中心相同跨催化环境共享。

## 7. 建议架构与实施顺序（待结合用户建议确定）

调用链：

```text
通用网络生成 + atom mapping
  → 可选 single-event 标注
  → 映射反应中心
  → 反应族 TS 构造 + 环境/统计约定
  → 通用对称性算法
  → 单个 TS 对称比及动力学通道归并
  → n_e、子类、共享参数
  → 保存 / 导出 / 英文 UI
```

建议少量职责完整的模块：

- 拟新增 `src/chemweave/symmetry.py`：通用增强图、群阶、受支持立体修正及计算明细；不识别具体分子名称来返回答案。
- 修改 `single_event_calculators.py`：各族适用性、TS 构造和分类；先集中放置，不需要大量插件文件。
- 后续拟新增 `catalytic_environment.py`：吸附和位点约束、环境允许的操作、统计约定；由用户目标决定第一批实际范围。
- 修改 `single_event.py`：保留调度职责，增加结构化计算依据、TS 集合和最终因子语义。

分阶段建议：

| 阶段 | 内容与涉及文件 | 验收 |
| --- | --- | --- |
| A | 整理 `tests/data/symmetry_reference.json` | 保留 SI 原始两列、来源、结构导入问题和约定；对异常分类，不伪造通过 |
| B | `symmetry.py`、`tests/test_symmetry.py` | 先通过结构/几何明确的子集；逐项比较群阶、label 修正和 sigma |
| C | TS 族构造、`single_event.py`、`tests/test_single_event.py`、验证规则和示例 | 新分子自动计算、通道归并可追溯、不依赖 match_count |
| D | 催化环境、必要的通用状态扩展、环境测试 | β-裂解或用户选定族中，分别验证自由态、抵消假设和显式环境结果 |
| E | `models.py`、`cli.py`、`web.py`、`static/app.js/index.html/style.css`、两份 README、CLI/web 测试 | 保存重放、CSV、英文 UI、旧 schema 兼容、可选功能不干扰生成 |

阶段 B–C 的首个候选族：H 原子从无环饱和烃夺氢。以正文 Scheme 8 为起点，扩展到乙烷、丙烷、异丁烷、新戊烷等未硬编码结构。这些是算法验证选择，不意味着忽略催化目标，也未由用户最终确认。

在所选自由态模型和通道约定下，乙烷夺氢的 6、丙烷端位/中间的 6/2 是待独立验证的扩展预期，不是 SI 反应参考值。正丁烷仲位夺氢可检验 TS 立体构型和同产物 TS 合并，不能仅用数 H 的方法完成。

验证至少包括：

- 论文分子 sigma 与正文 TS 示例，分别验证。
- 未硬编码分子、同位素和立体变化。
- 原子重编号、等价 SMILES、SMARTS 增加旁观原子、提前停止/修改 match_count 不改变物理结果。
- 同一约定下正逆向一致性；反向及自由基状态未验证时不声称支持。
- 环境对称性确实抵消和不抵消的受控案例。
- 不支持结构的明确结果；不以错误数值满足“总有输出”。
- 参数拟合好坏不能单独证明 n_e 正确：指前因子可能吸收常数倍误差。

建议新增审计字段覆盖：算法/版本、原子与键规范化约定、各物种和 TS 的群阶/label 修正/sigma、TS 规则版本、环境、通道合并约定、symmetry_ratio、最终 n_e、参考来源和验证状态。现有 `factor` 可兼容保留为最终使用的因子。

第一版不承诺任意表面/孔道、金属配位几何、任意受阻转子势能积分、耦合多环立体、完整共振体系、所有未定义立体异构体枚举、覆盖度相关或不等能路径合并。

## 8. 环境与验证记录

- 本地 base：`/Users/liuyan/miniconda3/bin/python`。研究时实际读取版本为 Python 3.12.2、RDKit 2026.03.6、NetworkX 3.5。
- 项目 `.venv-web/bin/python` 先前可运行测试。ml 环境先前含 RDKit 2025.03.5；下一轮重新核实。
- `pyproject.toml` 当前核心依赖仅 RDKit；NetworkX 不是项目已声明依赖。
- 建议初期用 NetworkX 对小图精确枚举/交叉检查；高对称大图需评估直接求群阶的后端，不能截断枚举后返回错误计数。依赖选择未最终确定。
- 历史记录：初始实现曾完整测试 84 passed；分离甲烷逻辑后，`tests/test_single_event.py tests/test_web.py tests/test_cli.py` 曾运行 30 passed。**这些是前序记录，不是本次研究阶段重新运行的结果，也不能声称当前完整测试为 86 passed。**
- 本轮只进行了文献读取、SI 统计和 RDKit InChI 导入审计，没有实现或测试新的对称性算法。
- `.doc` 可直接只读提取：`textutil -convert txt -stdout /Users/liuyan/Downloads/jcc23788-sup-0001-suppinfo.doc`。
- base 中有 PyMuPDF (`fitz`)，用于读取/渲染本地正文。不要依赖 `/tmp` 中上轮临时文件仍然存在。
- 没有创建新的安装包，也没有部署或启动需要保留的服务。

