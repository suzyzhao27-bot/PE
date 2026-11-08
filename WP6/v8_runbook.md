# v8.0.1 真实运行清单与复现边界

执行设备：本机Windows，2026-10-26。工作目录为项目根目录；实际解释器`.\.venv\Scripts\python.exe`，Python 3.14.7、NumPy 2.5.4、SciPy 1.18.1、Matplotlib 3.11.2。系统Matplotlib默认缓存目录权限不足；制图脚本显式使用项目已有`.mplconfig`，六图生成退出0。无新依赖安装、无外部实验、无新PNP或MD。

## 依赖顺序与真实命令

```powershell
.\.venv\Scripts\python.exe src\data_ingest\v8_step1_units.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_prescreen.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_verify.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_compare.py tables
.\.venv\Scripts\python.exe src\effective_network_model\v8_compare.py compare
.\.venv\Scripts\python.exe src\effective_network_model\v8_validate_seeds.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_validate_domain.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_validate_physics.py fixed
.\.venv\Scripts\python.exe src\effective_network_model\v8_validate_physics.py water
.\.venv\Scripts\python.exe src\effective_network_model\v8_validate_transfer.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_decide.py
.\.venv\Scripts\python.exe src\effective_network_model\v8_paper_figures.py
```

这些命令本轮实际退出0。第3步后的局部核接口显式化`r_s`后，第2步预筛与8项技术检查定向重跑，输出与门槛结论未变；制图中将图3纵轴限制为非负后定向重绘。主稳态表和优化格未因这些改动重算。所有模型计算由Codex本机执行，未记作研究者本人复算或导师独立验收。

## 一次小规模回放

首选第3步：运行`v8_verify.py`，它不依赖大规模拟合表，约数秒，输出8行`WP3/v8_技术检查.csv`和单图`v8_技术探针.json`。如需检查完整数字，按上面顺序运行；`v8_compare.py compare`要求先有稳态表，后续验证要求先有24行选定参数。大规模运行仅在可控本机目录，原16号输入只读。

## 实际计数与关键状态

预筛2352行但**不解盐网络**；主稳态2592次、独立种子1296次、24×24域240次、十情景384次，合计4512次盐网络线性解，低于冻结5616上限。另144次条件水力网络估计，正好等于冻结上限。Ncal从Scal β0复用；预筛排除21个q不进行Tcal稳态拟合；γ两情景复用稳态。主格点80520次实算损失，30744次合法跳过状态；L装置549格不增加材料稳态解。8项新增技术检查8/8通过；主Tcal科学判定`SCIENTIFIC_NOT_SUPPORTED`。本轮没有结果后扩格、删点或更改阈值。

本项目原有`.planning/2026-10-15-pe-research/progress.md`保存逐阶段检查点；正式看板第9节保存四态。完整输入来源、已有来源哈希与复制路径仍由16号`source_manifest.csv`登记。`WP4/v8_全部逐点预测.csv`每点保存冻结配置哈希，`WP6/v8_actual_spec.json`与16号冻结配置相同。本轮不另对全部输出增加内容哈希。17号正式包已获研究者同意创建；包内第1步量纲及8项技术检查代表回放退出0，实际记录见`17_学生物理约束接续成果_v8.0.1/checks/v8_replay_log.md`。正式包未重复运行完整矩阵。
