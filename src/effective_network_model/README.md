# 导师 v5.0.0 主模型的仓库工作副本

本目录是 `原始资料/赵珈绮/08_主模型交接_v5.0.0/02_程序/` 中核心、共享拟合和有限 PNP 程序及其必要轻量输入的**独立工作副本**。保留导师交付的原有目录结构与源码内容；未复制交付包预生成的结果。`原始资料/` 始终只读。`core/inputs/` 与 `fit/inputs/` 的同表副本是上游两个独立入口所需，不能当作两组独立实验。

在仓库根目录安装完毕的 Python 环境中运行以下命令；脚本不安装依赖、不写入仓库，要求每次指定一个**尚不存在、且位于仓库外**的非同步输出目录：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\src\effective_network_model\run_reproduction.ps1 -OutputDirectory 'D:\CodexWorkspaces\pe-v5-run-001'
```

默认使用仓库 `.venv/Scripts/python.exe`；其他已有环境可用 `-PythonExecutable` 指定。运行依次执行 core 12 项、六端点计算、共享拟合、PNP 新增 7 项与旧版 10 项检查，以及六个代表空间工况。输出目录中保留复制的程序和相应结果；原工作副本不会就地运行。单独运行 `pnp/run_energy_examples.py` 可重现交付中的示范能场，但不属于这条最小主模型复算路径。

模型结论与数值边界见根目录 `复现记录_20261017_v5.md`、`阶段总结/20261017_v5一页推导与正反结果.md`。此处参数是六个旧厚度端点的共享解释性拟合；`ell` 非 Debye 长，局部半径关系为条件假设，空间工况不是独立实验预测。此入口的通过不代表学生已在本人设备完成复算。

每次运行在指定输出目录生成 `run_record.json`，记录实际解释器路径、版本、六步命令与退出码；失败时也保留已执行步骤及失败状态。学生在本人设备运行后，应连同 `core/student_results/run_summary.json`、`core/student_results/shared_model_comparison.csv`、`fit/results/fit_results.json`、`pnp/energy_test_report.json`、`pnp/test_report.json` 和 `pnp/outputs/root_fitted_energy_spatial_check.json` 一起核对。运行记录是可复核的执行轨迹，学生仍须本人说明操作、结果及不一致项，不能仅凭他人机器生成的记录替代本人复算。无需把仓库外整套结果复制进原始资料目录。

本工作区新增的 `time_shape_check.py` 不属于导师原始交付。用 `python -B .\src\effective_network_model\time_shape_check.py` 可只读检查六条旧曲线在共同时间标签 `20/40/60/100 min` 的归一化形状。常 `P`、充分混合有限池且 `theta≥0` 的曲线不低于 `t/120`；负差额只表明与该条件模型不完全一致，不提供误差显著性或迟滞机制。80/85 min 歧义行不参与检查。
同一入口还输出 `late_over_early_60=[κ(120)−κ(60)]/[κ(60)−κ(0)]`，仅描述前后等长时间窗的电导增量比，不等于渗透率比。数值与限制见 `阶段总结/20261017_旧时间曲线内部检验.md`。

本工作区另增 `boundary_refit_check.py`，用已有 NumPy/SciPy 对预设的六个稳定 Robin 边界分别重拟合六端点，对 `r=0/0.75` 做固定边界留一，并在嵌套留一中每折只用五点选择边界。命令：`python -B .\src\effective_network_model\boundary_refit_check.py`。只打印结果，不改动原始资料、导师程序或参考参数；发现与限制见 `阶段总结/20261017_共享边界离散重拟合.md`。

图3分母口径的独立替代图由 `plot_pnp_baselines.py` 生成。先在仓库外 PNP 工作副本运行 `run_energy_examples.py` 和 `run_fitted_field_check.py`；再把所得 `representative_energy_results.json`、`root_fitted_energy_spatial_check.json` 分别传给脚本的 `--representative-json`、`--fitted-json`，用 `--output` 指定已有输出目录内的 PNG 路径。本工作区审核图为 `data/processed/figures/pnp_same_sigma_baseline.png`。图中普通 24 格工况与另一次强场转移工况均按**同 `σ` 零能场**归一化；它是虚拟几何的数值比较，不替代导师原图，也不是实验通量。
