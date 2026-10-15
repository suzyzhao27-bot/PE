# 超薄 PE 纤维膜的界面荷电、离子选择性与盐输运

本目录用于双人协作维护代码、说明文档和轻量配置。`原始资料/` 保留原始交付包、文献、实验快照、历史版本和既有结果，按根目录 `.gitignore` 规则不进入后续 Git 管理。

## 目录结构

```text
src/
  v2_training/             标准库训练、情景扫描与测试
  long_tube_diagnostic/    长管 PB-NP 受限诊断
  surface_pb/              纤维外水域 PB 与荷电调节
  finite_pnp/              有限二维 PNP 与五排纤维计算
原始资料/                  原始交付和数据，只读参考，不纳入 Git
项目说明与行动清单.md      研究背景、边界与待办
资料分类索引.md            资料地图
```

本次仅创建 `src/`。未创建独立 `scripts/`，因为现有入口脚本依赖同目录中的模块；保持同目录可避免改变其导入方式。

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

## 未复制的脚本

以下脚本留在 `原始资料/`，因为它们服务于原始资料溯源、文献获取、数据审计或制报环境，而不是当前共享计算入口：

- `04_文献_v1阅读包` 的文献获取/构建脚本；它们涉及来源资产和版权资料。
- `05_真实数据_v1核查/audit_sources.py` 与 `verify_original_conversion.py`；它们依赖原始 XLS、OPJU、PDF、DOCX 快照。
- v2/v3 的 `build_*_report.py`、`research/` 脚本与验收脚本；它们依赖制报资源、原机来源或研究数据登记。

若需要将其中某一类脚本纳入协作，先确认其原始数据是否可公开、是否需要去标识化，以及其输入是否应当在仓库外通过配置路径提供。

## 路径、数据与协作注意事项

- 当前复制的主计算脚本未发现会打开原机绝对路径的代码。原始 README 中的 macOS 解释器路径仅是历史运行记录。
- 代码使用同目录导入，例如 `from solver import ...` 和 `from models import ...`。请在对应模块目录内运行，或保持脚本与同目录模块一起移动。
- 已有 `.gitignore` 会忽略 `*.npz`、原始 Office/PDF/Origin 文件、压缩包和 `原始资料/`。新运行产生的 CSV/JSON/PNG 是否进入版本库，须在代码审阅时按“可复现、小型、非敏感”原则逐项决定。
- 没有新增 `requirements.txt`、虚拟环境或依赖安装。Python、NumPy、SciPy、matplotlib 的具体版本兼容范围仍为**待确认**。
- 复制后代码与原始副本分别演化。每次共同修改前，应在提交说明中记录对应原始版本/文件哈希和物理假设变化。
