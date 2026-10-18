# 导师 v5.0.0 主模型的仓库工作副本

本目录是 `原始资料/赵珈绮/08_主模型交接_v5.0.0/02_程序/` 中核心、共享拟合和有限 PNP 程序及其必要轻量输入的**独立工作副本**。保留导师交付的原有目录结构与源码内容；未复制交付包预生成的结果。`原始资料/` 始终只读。`core/inputs/` 与 `fit/inputs/` 的同表副本是上游两个独立入口所需，不能当作两组独立实验。

在仓库根目录安装完毕的 Python 环境中运行以下命令；脚本不安装依赖、不写入仓库，要求每次指定一个**尚不存在、且位于仓库外**的非同步输出目录：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\src\effective_network_model\run_reproduction.ps1 -OutputDirectory 'D:\CodexWorkspaces\pe-v5-run-001'
```

默认使用仓库 `.venv/Scripts/python.exe`；其他已有环境可用 `-PythonExecutable` 指定。运行依次执行 core 12 项、六端点计算、共享拟合、PNP 新增 7 项与旧版 10 项检查，以及六个代表空间工况。输出目录中保留复制的程序和相应结果；原工作副本不会就地运行。单独运行 `pnp/run_energy_examples.py` 可重现交付中的示范能场，但不属于这条最小主模型复算路径。

学生先完成不依赖 NumPy/SciPy 的 core 本人复算，可在上述命令末尾添加 `-CoreOnly`，并为这次运行指定新的仓库外输出目录。此模式只复制 `core/`，运行 `test_model.py` 和 `run_student.py` 两步；`run_record.json` 标记 `mode=core_only`。具备已安装 NumPy/SciPy 的环境后，再用新的输出目录运行不带 `-CoreOnly` 的完整六步入口。两个目录分别保留，不能把本机成功记录当作学生本人操作。

模型结论与数值边界见根目录 `复现记录_20261017_v5.md`、`阶段总结/20261017_v5一页推导与正反结果.md`。此处参数是六个旧厚度端点的共享解释性拟合；`ell` 非 Debye 长，局部半径关系为条件假设，空间工况不是独立实验预测。此入口的通过不代表学生已在本人设备完成复算。

固定浓度比的虚拟盐度对照使用 `pnp/salinity_sensitivity.py`，从仓库根目录运行；命令、冻结参数、同 `σ` 零场分母与 24/32 格结果见 `复现记录_20261017_v5.md`。输出路径必须位于已有目录。该对照只改变预设浴液浓度，不是旧电导记录的浓度重建或真实 PE 盐度响应。

每次运行在指定输出目录生成 `run_record.json`，记录模式、实际解释器路径、版本、已执行命令与退出码；失败时也保留已执行步骤及失败状态。学生在本人设备运行后，先核对 `core/student_results/run_summary.json`、`core/student_results/shared_model_comparison.csv`。完整模式另核对 `fit/results/fit_results.json`、`pnp/energy_test_report.json`、`pnp/test_report.json` 和 `pnp/outputs/root_fitted_energy_spatial_check.json`。运行记录是可复核的执行轨迹，学生仍须本人说明操作、结果及不一致项，不能仅凭他人机器生成的记录替代本人复算。无需把仓库外整套结果复制进原始资料目录。

本工作区新增的 `time_shape_check.py` 不属于导师原始交付。用 `python -B .\src\effective_network_model\time_shape_check.py` 可只读检查六条旧曲线在共同时间标签 `20/40/60/100 min` 的归一化形状。常 `P`、充分混合有限池且 `theta≥0` 的曲线不低于 `t/120`；负差额只表明与该条件模型不完全一致，不提供误差显著性或迟滞机制。80/85 min 歧义行不参与检查。
同一入口还输出 `late_over_early_60=[κ(120)−κ(60)]/[κ(60)−κ(0)]`，仅描述前后等长时间窗的电导增量比，不等于渗透率比。数值与限制见 `阶段总结/20261017_旧时间曲线内部检验.md`。

若要检查六条旧曲线能否共享一个非负时间开启过程且只差固定强弱系数，仍用此入口加 `--common-gate-output data/processed/v5_common_gate_order.csv`。它使用 20/40/60/100 min 四个共同时间点，写出 60 项端点排序与中途归一化排序比较；不拟合 `g(t)`。成立需六列初始池状态、体积、面积和读数比例可比，具体推导及结果见同一时间形状阶段总结。

更少依赖跨列幅值的检验加 `--crossing-output data/processed/v5_normalized_curve_crossings.csv`。它只比较每列各自归一化后的曲线是否在四个共同时间点改变高低顺序；共同非负 `g(t)` 与每膜固定系数的有限池形状不允许曲线交叉。输出 15 对厚度的逐时差值及交叉差幅，限制见同一阶段总结。

本工作区另增 `boundary_refit_check.py`，用已有 NumPy/SciPy 对预设的六个稳定 Robin 边界分别重拟合六端点，对 `r=0/0.75` 做固定边界留一，并在嵌套留一中每折只用五点选择边界。命令：`python -B .\src\effective_network_model\boundary_refit_check.py`。只打印结果，不改动原始资料、导师程序或参考参数；发现与限制见 `阶段总结/20261017_共享边界离散重拟合.md`。

两个旧厚度比例的定向约束检验使用 `ratio_constraint_check.py`，以同一个 `r=0, a=25 nm` 场求共享 `ws、ell`，强制匹配 350/200 与 2600/700；700/350 及 65/120 保持未用于求解。从根目录运行 `.\.venv\Scripts\python.exe -B .\src\effective_network_model\ratio_constraint_check.py --output .\data\processed\v5_ratio_constraint_check.csv`。输出是诊断性相对量，不替换导师全六点参考拟合；见 `阶段总结/20261017_共享边界离散重拟合.md`。

可在同一命令后加 `--tradeoff-output .\data\processed\v5_ratio_tradeoff.csv`，比较原六点参考、两比例条件解和三相邻比例折中解的比值误差。三比例折中仍只用共享 `ws、ell`，输出 12 行含未用于选择比例目标的 65/120 对照；多初值局部优化不是全域不可行证明。

最厚标签的参数辨识作用由 `parameter_sensitivity.py` 在**冻结六点参考解**处评估。运行 `.\.venv\Scripts\python.exe -B .\src\effective_network_model\parameter_sensitivity.py --output .\data\processed\v5_parameter_sensitivity.csv`。输出全六点与逐点留出的七行，列出消去共同幅值后的两列 log 参数灵敏度奇异值及条件数。需要现有 NumPy；结果仅是局部几何诊断，不是五点重拟合或统计置信区间。

同一命令加 `--leverage-output .\data\processed\v5_heldout_leverage.csv`，得到六行“沿其余五点最弱参数方向，留出点相对变化 / 五点内部均方根变化”。该比值只在冻结参考解附近有局部意义，不等于留一预测误差倍率。

图3分母口径的独立替代图由 `plot_pnp_baselines.py` 生成。先在仓库外 PNP 工作副本运行 `run_energy_examples.py` 和 `run_fitted_field_check.py`；再把所得 `representative_energy_results.json`、`root_fitted_energy_spatial_check.json` 分别传给脚本的 `--representative-json`、`--fitted-json`，用 `--output` 指定已有输出目录内的 PNG 路径。本工作区审核图为 `data/processed/figures/pnp_same_sigma_baseline.png`。图中普通 24 格工况与另一次强场转移工况均按**同 `σ` 零能场**归一化；它是虚拟几何的数值比较，不替代导师原图，也不是实验通量。

图4的可追溯替代图由 `plot_report_limits.py` 绘制，输入为导师交付包 `fit/results/leave_one_out.csv` 与本目录 `core/inputs/Book1__Sheet1.csv`，输出为现有 `data/processed/figures/v5_report_limits.png`。命令示例（从仓库根目录运行）：

```powershell
python -B .\src\effective_network_model\plot_report_limits.py --loo-csv '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\fit\results\leave_one_out.csv' --book1-csv '.\src\effective_network_model\core\inputs\Book1__Sheet1.csv' --output '.\data\processed\figures\v5_report_limits.png'
```

左图是六个旧端点的逐点留一预测/观测，右图仅用共同确认的时间标签；灰线为条件模型下界而非拟合。此图不提供独立新样品预测、统计误差或动态机制结论。

曲线交叉的补充图沿用此脚本的 `--layout crossings`，不用留一 CSV：

```powershell
$env:MPLCONFIGDIR = (Resolve-Path -LiteralPath '.venv').Path
.\.venv\Scripts\python.exe -B .\src\effective_network_model\plot_report_limits.py --layout crossings --book1-csv .\src\effective_network_model\core\inputs\Book1__Sheet1.csv --output .\data\processed\figures\v5_time_curve_crossings.png
```

左幅画六条各自按 120 min 增量归一化的旧曲线；右幅画 120 nm 减 2600 nm 的完成比例，展示 20–40 min 高低次序反转。只使用 A/C 共同确认的 0/20/40/60/100/120 min 行；连线仅帮助阅读，不能解释为新增观测。它是图四时间限制的补充证据，不替换原图四或已审阅的 `v5_report_limits.png`。

图1的可追溯替代图由 `plot_report_endpoints.py` 绘制，仅使用导师交付包 `fit/results/endpoint_fit.csv` 与 `thickness_curve.csv`，输出为现有 `data/processed/figures/v5_report_endpoints.png`。从仓库根目录运行：

```powershell
python -B .\src\effective_network_model\plot_report_endpoints.py --endpoint-csv '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\fit\results\endpoint_fit.csv' --curve-csv '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\fit\results\thickness_curve.csv' --output '.\data\processed\figures\v5_report_endpoints.png'
```

黑点为旧 0–120 min 电导增量分别除以 2600 nm 旧端点增量；蓝线及空心点为参考模型预测分别除以其 2600 nm 预测值。两个系列各自以 2600 nm 为 1，均为无量纲。跨厚度解释为相对盐转移仍须共同校准、测量几何及体积历史可比。模型使用全部六点训练；图仅覆盖旧标签区间，不绘制导师原图中身份未确认的紫色曲线。

图二主题的**闭合规则补充图**由 `plot_report_structure.py` 默认版式绘制。输入为导师交付包 `core/student_results/thickness_scan.csv`、`local_confinement_scan.csv` 和本目录 `core/inputs/model_parameters.json`；固定 `ws` 的入口长度对照直接调用同一 `theory_model.log_resistance_nm`。从仓库根目录运行：

```powershell
python -B .\src\effective_network_model\plot_report_structure.py --thickness-csv '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\core\student_results\thickness_scan.csv' --local-csv '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\core\student_results\local_confinement_scan.csv' --parameters-json '.\src\effective_network_model\core\inputs\model_parameters.json' --output '.\data\processed\figures\v5_report_structure.png'
```

图中全部是冻结参数下的模型相对量；`a` 与 `ws(a)` 未经湿态孔道实测证实。

若需与导师 v5 原图二的**左右比较对象一致**，沿用此脚本而不覆盖上述闭合规则对照图：

```powershell
$env:MPLCONFIGDIR = (Resolve-Path -LiteralPath '.venv').Path
.\.venv\Scripts\python.exe -B .\src\effective_network_model\plot_report_structure.py --layout mentor-v5 --parameters-json .\src\effective_network_model\core\inputs\model_parameters.json --output .\data\processed\figures\v5_report_h_a.png
```

新图左侧固定 `a=25 nm` 绘出导师相同四个厚度的 `w(z/h)`；右侧冻结参考参数并按 `ws(a)=ws(25)×25/a` 绘出三种 `a=20/25/30 nm` 随 `h` 的相对 `P`，共同分母为 `P(2600 nm,25 nm)`。它是导师图二主题的仓库复算图，不是导师原 PNG 的逐像素复制；原件保留在 `原始资料/`，`v5_report_structure.png` 仍保留闭合规则对照用途。

**补充 PNP 对照（不是导师 v5 图二）**由 `plot_charge_entry_comparison.py` 绘制，使用导师冻结的 `pnp/outputs/representative_energy_results.json` 中四个匹配的 24 格工况，输出 `data/processed/figures/v5_charge_entry_comparison.png`。从仓库根目录运行：

```powershell
$env:MPLCONFIGDIR = (Resolve-Path -LiteralPath '.venv').Path
.\.venv\Scripts\python.exe -B .\src\effective_network_model\plot_charge_entry_comparison.py --results-json '.\原始资料\赵珈绮\08_主模型交接_v5.0.0\02_程序\pnp\outputs\representative_energy_results.json' --output '.\data\processed\figures\v5_charge_entry_comparison.png'
```

左图为开路电位响应，右图为共同盐通量相对中性零能场通量；四种情景使用相同五排虚拟几何及 30/1 mM NaCl。开路电位并非实测离子选择性系数，开路下两个离子跨膜通量近似相等。导师 v5 的正式图二为 `原始资料/赵珈绮/08_主模型交接_v5.0.0/04_结果图/图2_厚度与局部限域的分工.png`；`v5_report_structure.png` 是其问题的条件扫描补充图，构图和比较切片与导师原图不同。

图二厚度机制的局部斜率分解由 `thickness_balance.py` 完成，直接读取本目录冻结的 `core/inputs/model_parameters.json`，不重拟合、不改导师图。运行：

```powershell
.\.venv\Scripts\python.exe -B .\src\effective_network_model\thickness_balance.py --output .\data\processed\v5_thickness_balance.csv
```

若要量化图二预设 `ws(a)=ws_ref×25/a` 对 200/350 nm 反常排序的要求，可在同一命令加 `--local-scale-threshold-output .\data\processed\v5_local_scale_threshold.csv`。此诊断固定 200 nm 的 `a=25 nm` 与全部参考参数，仅反求 350 nm 达到等传盐及旧端点条件比例所需的 `a`；不是按膜调参后的新模型，更不是实测孔径。

CSV 对旧六个厚度标签给出膜内阻力占比、路径长度项、双面场减弱项及其和 `dlnP/dlnh`。脚本以既有 `thickness_log_slope` 和数值差分分别核对六点斜率；每个厚度只用冻结的 `r=0, a=25 nm, ws_ref, ell`。正斜率是该位置的条件模型趋势，不是旧实验观测导数。通俗解释和四张导师原图的坐标说明见 `阶段总结/20261018_v5科研通俗说明.md`。

`build_report_preview.py` 将现有 `阶段总结/20261017_v5学生报告工作稿.md` 排成 `阶段总结/20261017_v5学生报告审阅预览.pdf`，用于检查正文、表格与四图的阅读顺序。此可选入口需要 ReportLab、Pillow 和 Windows `SimHei` 字体，不属于 core/fit/PNP 科研复算依赖；在无这些排版依赖的设备上直接阅读 Markdown 即可。PDF 页脚明确标记内部审阅稿，不能替代学校模板、学生本人稿或正式提交。命令示例：

仓库既有 `.gitignore` 忽略 `*.pdf`，本机预览 PDF 保留在上述位置；跨设备可由受版本控制的 Markdown、四图和本脚本重新生成。不要为该临时排版件解除现有 PDF 忽略规则。

```powershell
python -B .\src\effective_network_model\build_report_preview.py --source '.\阶段总结\20261017_v5学生报告工作稿.md' --output '.\阶段总结\20261017_v5学生报告审阅预览.pdf'
```
