# 导师 v5.0.0 主模型的仓库工作副本

本目录是 `原始资料/赵珈绮/08_主模型交接_v5.0.0/02_程序/` 中核心、共享拟合和有限 PNP 程序及其必要轻量输入的**独立工作副本**。保留导师交付的原有目录结构与源码内容；未复制交付包预生成的结果。`原始资料/` 始终只读。`core/inputs/` 与 `fit/inputs/` 的同表副本是上游两个独立入口所需，不能当作两组独立实验。

在仓库根目录安装完毕的 Python 环境中运行以下命令；脚本不安装依赖、不写入仓库，要求每次指定一个**尚不存在、且位于仓库外**的非同步输出目录：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\src\effective_network_model\run_reproduction.ps1 -OutputDirectory 'D:\CodexWorkspaces\pe-v5-run-001'
```

默认使用仓库 `.venv/Scripts/python.exe`；其他已有环境可用 `-PythonExecutable` 指定。运行依次执行 core 12 项、六端点计算、共享拟合、PNP 新增 7 项与旧版 10 项检查，以及六个代表空间工况。输出目录中保留复制的程序和相应结果；原工作副本不会就地运行。单独运行 `pnp/run_energy_examples.py` 可重现交付中的示范能场，但不属于这条最小主模型复算路径。

模型结论与数值边界见根目录 `复现记录_20261017_v5.md`、`阶段总结/20261017_v5一页推导与正反结果.md`。此处参数是六个旧厚度端点的共享解释性拟合；`ell` 非 Debye 长，局部半径关系为条件假设，空间工况不是独立实验预测。此入口的通过不代表学生已在本人设备完成复算。
