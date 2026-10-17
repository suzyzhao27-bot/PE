# 超薄 PE 纤维膜的界面荷电、离子选择性与盐输运

本项目由同一位研究者通过两台电脑和两个 Codex 账号维护。`dev-zhao`、`dev-yi` 是设备开发分支，`main` 是统一的权威共享分支。`原始资料/` 保留原始交付包、文献、实验快照、历史版本和既有结果，按根目录 `.gitignore` 规则不进入后续 Git 管理。

2026-10-17 导师将当前任务明确为基于已有实验的纯理论计算，并交付 v5.0.0 主模型。当前路线见[正式看板第 6 节](项目进展计划与任务看板.md)、[导师反馈评阅](阶段总结/20261017_导师反馈与v5主模型接入.md)和[本工作区复现记录](复现记录_20261017_v5.md)；[v5 仓库工作副本及复算入口](src/effective_network_model/README.md)已接入。旧 O/S/W 等新增实验草案保留追溯，不作为当前计算关卡。

## 项目目录说明

```text
.planning/                 按任务归档的规划、发现、实施方案和进度记录
src/
  v2_training/             标准库训练、情景扫描与测试
  long_tube_diagnostic/    长管 PB-NP 受限诊断
  surface_pb/              纤维外水域 PB 与荷电调节
  finite_pnp/              有限二维 PNP 与五排纤维计算
  data_ingest/             Origin 只读解码、盐输运长表和 QC
data/processed/            审阅后的轻量导出、审计清单和 QC 图
paper/review-traces/       结果到声明等独立审查轨迹
阶段总结/                  各阶段和工作包的正式总结
原始资料/                  原始交付和数据，只读参考，不纳入 Git
AGENTS.md                  两台设备共同遵守的长期工作规范
项目进展计划与任务看板.md   唯一的正式任务状态来源
项目说明与行动清单.md      研究背景、边界与待办
复现记录_*.md              环境、运行和验证证据
资料分类索引.md            资料地图
```

| 位置 | 作用与代表内容 | 后续归档规则 |
|---|---|---|
| `src/` | 科研代码、测试、计算和数据接入脚本 | 进入现有对应模块；新增模块目录前先说明用途并获得确认 |
| `data/processed/` | 审阅后的轻量 CSV/JSON、审计清单和 QC 图 | 只放可追溯的处理结果；不得替代或覆盖 `原始资料/` |
| `src/**/outputs/`、`src/outputs/` | 可由代码重新生成的运行输出 | 作为本地运行产物，不作为正式科研记录的唯一副本 |
| `paper/review-traces/` | 结果到声明等独立审查过程与判定 | 同类审查继续沿用现有层级，不在根目录另建副本 |
| `阶段总结/` | P0、P1、P2 等阶段或工作包的正式总结 | 每项结论保留一份权威总结；草案在文件名中明确版本 |
| `.planning/` | 需要持续规划或跨设备交接的任务过程记录 | 使用 `YYYY-MM-DD-任务名称/`；同一任务持续更新原目录，完成后作为历史保留 |
| `原始资料/` | 原始交付、文献、实验快照、历史版本和整理记录 | 只读保存；`03_发送用ZIP/`、`90_历史与工作底稿/` 等历史区不接收日常新产出 |
| 根目录 Markdown | 正式入口、看板、行动清单、资料索引和复现记录 | 不放临时笔记；任务状态只以正式看板为准 |

`.venv/`、`.mplconfig/`、`.test_tmp/`、`__pycache__/` 以及模块输出目录属于本机环境、缓存或可再生成产物，不作为双设备科研记忆。根目录已有 `原始资料.zip` 是历史交付归档，不应继续在根目录堆放新的运行压缩包。

数据接入脚本放在 `src/data_ingest/`；它们不改写 `原始资料/`。同类科研产出应集中进入上述既有位置，不为每次运行、测试或分析另建文件夹，也不擅自移动、重命名或覆盖已有成果。

## 代码来源与复制映射

| 共享副本 | 原始位置 |
|---|---|
| `src/v2_training/models.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/models.py` |
| `src/v2_training/pipeline.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/pipeline.py` |
| `src/v2_training/run_demo.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/run_demo.py` |
| `src/v2_training/run_scenario.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/run_scenario.py` |
| `src/v2_training/test_models.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/test_models.py` |
| `src/v2_training/tests_pipeline.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/03_学生计算_v1训练/tests_pipeline.py` |
| `src/long_tube_diagnostic/long_tube_np.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/02_长管PB_NP诊断/long_tube_np.py` |
| `src/long_tube_diagnostic/pb_v1_readonly_copy.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/02_长管PB_NP诊断/pb_v1_readonly_copy.py` |
| `src/long_tube_diagnostic/run_diagnostic.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/02_长管PB_NP诊断/run_diagnostic.py` |
| `src/long_tube_diagnostic/test_diagnostic.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/02_长管PB_NP诊断/test_diagnostic.py` |
| `src/long_tube_diagnostic/plot_charge_flux.py` | `原始资料/赵珈绮/01_研究框架与资料_v2/02_长管PB_NP诊断/plot_charge_flux.py` |
| `src/surface_pb/fibre_pb.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/surface/fibre_pb.py` |
| `src/surface_pb/run_surface.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/surface/run_surface.py` |
| `src/surface_pb/test_surface.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/surface/test_surface.py` |
| `src/surface_pb/plot_surface.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/surface/plot_surface.py` |
| `src/finite_pnp/solver.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/solver.py` |
| `src/finite_pnp/run.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/run.py` |
| `src/finite_pnp/run_multifiber.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/run_multifiber.py` |
| `src/finite_pnp/refine_multifiber.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/refine_multifiber.py` |
| `src/finite_pnp/tests.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/tests.py` |
| `src/finite_pnp/plot_fields.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/plot_fields.py` |
| `src/finite_pnp/plot_multifiber_fields.py` | `原始资料/赵珈绮/02_荷电与选择性实算_v3/finite_pnp/plot_multifiber_fields.py` |

还复制了三个运行所需的小型文件：

- `src/v2_training/scenario_template.json`：`run_scenario.py` 的示例输入。
- `src/v2_training/DATA_TEMPLATE.csv`：真实实验数据输入模板，不含原始实验数据。
- `src/long_tube_diagnostic/pb_source_identity.json`：PB 只读副本的身份登记；其中的旧绝对路径仅为溯源文本，计算代码不读取该旧位置。

## 模块说明与运行

### 运行环境（Windows）

已在 Python 3.14.7 下验证。建议在项目根目录创建隔离环境；Windows 请使用 `py -3.14`，而非可能指向 Microsoft Store 占位程序的 `python`。

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

之后，所有需要科学计算库的命令均用 `..\..\.venv\Scripts\python.exe` 执行。绘图前，可在相应模块目录将 Matplotlib 缓存置于本地工作目录，避免默认缓存目录权限警告：

```powershell
New-Item -ItemType Directory -Force .mplconfig | Out-Null
$env:MPLCONFIGDIR = (Resolve-Path .mplconfig).Path
```

### `src/v2_training`

这是标准库训练模块，提供圆柱 PB、Donnan/中性对照、储液池曲线、参数情景和训练/冻结演示。

```powershell
cd src/v2_training
python -m unittest -v test_models tests_pipeline
python run_demo.py
python run_scenario.py --config scenario_template.json --out my_scenario_001
python pipeline.py demo
```

依赖：Python 标准库。`run_demo.py` 和 `pipeline.py` 默认生成 synthetic 输出。若要处理真实数据，需由使用者提供符合 `DATA_TEMPLATE.csv` 字段的轻量 CSV；真实数据来源、样品身份和实验边界均待确认。

### `src/long_tube_diagnostic`

这是受限长圆管 PB-NP 诊断，用于说明固定表面电荷不必然降低共同盐通量。它不是有限长度真实膜模型。

```powershell
cd src/long_tube_diagnostic
python test_diagnostic.py
python run_diagnostic.py
python plot_charge_flux.py  # 可选，需要已有 matplotlib
```

依赖：核心测试与计算只用 Python 标准库；绘图需要已有 matplotlib。`run_diagnostic.py` 会在本目录创建 `outputs/`。

### `src/surface_pb`

这是纤维外水域的非线性 PB 与自洽荷电调节模块，包含固定电荷、H/OH 或 Na/Cl 吸附、酸性位点等假设边界。

```powershell
cd src/surface_pb
python test_surface.py
python run_surface.py
python plot_surface.py  # 可选，需要已有 matplotlib
```

依赖：NumPy、SciPy；绘图额外需要 matplotlib。`run_surface.py` 会生成 `src/outputs/`。模型位点密度、结合参数和化学机制均为待实验约束的假设，不得表述为实测 PE 性质。

### `src/finite_pnp`

这是有限二维 Poisson-Nernst-Planck 主模型，计算实体纤维、间隙、入口和外水域的 Na+/Cl- 输运、理论开路电位及五排纤维结果。

```powershell
cd src/finite_pnp
python tests.py
python run.py
python run_multifiber.py
python refine_multifiber.py
python plot_fields.py             # 可选，需要已有 matplotlib 和 run.py 输出
python plot_multifiber_fields.py  # 可选，需要已有 matplotlib 和 refine_multifiber.py 输出
```

依赖：NumPy、SciPy；绘图额外需要 matplotlib。`run.py` 首次创建本目录 `outputs/`；因此请先执行 `run.py`，再运行 `run_multifiber.py` 或 `refine_multifiber.py`。后两者与绘图脚本都依赖此前生成的 JSON/CSV/NPZ 输出，但不依赖原始论文或实验数据。

### `src/data_ingest`

该模块通过已安装的 Origin COM 接口读取 OPJU 工作表。源工程只读处理，不调用保存；每个文件关闭后会复算 SHA256。随后将盐输运源表整理为带证据标记的长表，并生成仅用于转录/QC 的描述性汇总图。

```powershell
$src = '原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/source_snapshots'
$inputs = Get-ChildItem -LiteralPath $src -Filter '*.opju' | Select-Object -ExpandProperty FullName
.\src\data_ingest\extract_origin_opju.ps1 -InputPath $inputs -OutputRoot '.\data\processed\origin_exports'
.\.venv\Scripts\python.exe .\src\data_ingest\normalize_salt_transport.py
.\.venv\Scripts\python.exe .\src\data_ingest\analyze_salt_transport.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_origin_ingest.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_mass_balance.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_osw_protocol.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_independent_readouts.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_project_readiness.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_conductivity_calibration.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_osw_randomization.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_salt_transport_runs.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_process_salt_transport_runs.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_blank_registry.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_quality_thresholds.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_liquid_compatibility.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_liquid_compatibility_thresholds.py
.\.venv\Scripts\python.exe .\src\data_ingest\test_evaluate_liquid_compatibility.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_experiment_records.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_osw_protocol.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_independent_readouts.py
.\.venv\Scripts\python.exe .\src\data_ingest\assess_project_readiness.py
.\.venv\Scripts\python.exe .\src\data_ingest\conductivity_calibration.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_salt_transport_runs.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_blank_registry.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_quality_thresholds.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_liquid_compatibility.py
.\.venv\Scripts\python.exe .\src\data_ingest\validate_liquid_compatibility_thresholds.py
.\.venv\Scripts\python.exe .\src\data_ingest\evaluate_liquid_compatibility.py
```

当前长表不具备机制拟合资格：缺失电导-浓度校准、温度、有效面积、两池体积、取样/补液、润湿历史和独立样品映射。`salt_transport_origin_qc.png` 不能解释为浓度、通量或渗透率结果。

`data/processed/salt_transport/salt_transport_record_tracker.csv` 为 17 条曲线逐条登记这些外部记录。未知项必须写 `unknown`；校验通过只表示字段结构有效，不自动授权机制拟合。

阶段 2 的空表、校验和端到端运行顺序集中见 `阶段总结/P2_开测前执行清单_v0.md`。该清单是操作草案，不代表协议已冻结或实验已获准开始。

## 未复制的脚本

以下脚本留在 `原始资料/`，因为它们服务于原始资料溯源、文献获取、数据审计或制报环境，而不是当前共享计算入口：

- `04_文献_v1阅读包` 的文献获取/构建脚本；它们涉及来源资产和版权资料。
- `05_真实数据_v1核查/audit_sources.py` 与 `verify_original_conversion.py`；它们依赖原始 XLS、OPJU、PDF、DOCX 快照。
- v2/v3 的 `build_*_report.py`、`research/` 脚本与验收脚本；它们依赖制报资源、原机来源或研究数据登记。

若需要将其中某一类脚本纳入协作，先确认其原始数据是否可公开、是否需要去标识化，以及其输入是否应当在仓库外通过配置路径提供。

## 路径、数据与双设备协作注意事项

- 当前复制的主计算脚本未发现会打开原机绝对路径的代码。原始 README 中的 macOS 解释器路径仅是历史运行记录。
- 代码使用同目录导入，例如 `from solver import ...` 和 `from models import ...`。请在对应模块目录内运行，或保持脚本与同目录模块一起移动。
- 已有 `.gitignore` 会忽略 `*.npz`、原始 Office/PDF/Origin 文件、压缩包和 `原始资料/`。新运行产生的 CSV/JSON/PNG 是否进入版本库，须在代码审阅时按“可复现、小型、非敏感”原则逐项决定。
- 已提供 `requirements.txt`；在 Windows/Python 3.14.7 环境中已验证 NumPy 2.5.4、SciPy 1.18.1 与 Matplotlib 3.11.2。虚拟环境 `.venv/` 会被 Git 忽略。
- 复制后代码与原始副本分别演化。修改时应记录对应来源版本和物理假设变化；除既有溯源流程或明确完整性风险外，不额外重复执行文件内容哈希检查。
- 两台设备共同遵守根目录 `AGENTS.md`。切换设备前更新正式看板和当前任务的 `.planning/` 记录；新设备先安全获取最新 `main`，再在自己的设备分支继续工作。
