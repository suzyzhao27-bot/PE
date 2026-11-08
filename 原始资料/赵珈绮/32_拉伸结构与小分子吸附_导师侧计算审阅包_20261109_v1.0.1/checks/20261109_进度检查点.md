# 表面吸附任务进度

## 最近检查点（2026-11-08，几何代理与导师包收尾验收）

- 已完成：结构约束表、位点密度比和 `-RT ln(R)` 标尺；`face_charge.py` 九组逐列复现 PASS，最大差 `1.4210854715202004e-14`。
- 已生成：`data/processed/surface_adsorption/` 下三张 CSV、复算输出、三张 PNG/SVG 图、状态审计、干态输入合同、结构模板、湿态启动门、复算边界记录，以及自由膜/背衬膜几何外层代理表和代理合同。
- 已知阻塞：没有三组原子表面结构、辛酸/PE 力场和 GROMACS/WHAM/MBAR；不能生成 W(z)、Kdry、Rdry 或 Rwet。
- 远端核查：GitHub `main` 公共克隆同样没有结构坐标或分子模拟输入，不能从仓库补齐阻塞项。
- 代理判断：182/357/604 的层数为 11/6/4，自由膜外层周长占比分别为 0.18182/0.33333/0.50000；该结果只作为结构情景代理，不代表原子可达面积或吸附容量。
- 交付：旧 `v1.0.0` 包保留；本轮已按同一格式生成 `原始资料/赵珈绮/32_拉伸结构与小分子吸附_导师侧计算审阅包_20261109_v1.0.1/`，manifest 47 项、缺失 0 项，纳入图3、几何代理表、合同、说明和源码。
- 未完成：取得输入后运行 MD；MD 完成并通过收敛审计后才能生成 `W(z)`、`K_dry`、`R_dry` 或 `R_wet`。

## 验收证据

- `py -3.14 -m unittest -v src.test_surface_adsorption_baseline src.test_surface_geometry_proxy src.test_build_adsorption_package`：9/9 PASS。
- `py -3.14 -m py_compile src/mechanism_surface_adsorption_baseline.py src/surface_geometry_proxy.py src/build_adsorption_package.py`：退出码 0。
- 新包验收：manifest 与实际文件均为 47 项，清单覆盖、JSON 解析、SVG 解析和 `00/01/02/03 + checks/figure_data/figures/source_inputs/src` 结构均通过。
- `py -3.14 src/mechanism_surface_adsorption_baseline.py`：退出码 0；`face_charge_reproduction=PASS`；无字体警告。
- 图件已用 `view_image` 目视检查；未见裁切、遮挡或虚构吸附曲线。
