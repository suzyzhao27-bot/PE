# v6 WP1–WP6 复现包说明

本目录包含 `v6_reproduction.zip`、`inventory.csv`、`environment.json`、`requirements.txt`。压缩包保留项目相对路径，可解压到**新的空目录**后从根目录运行；本地工作区的 `原始资料/`、既有结果和历史版本不需覆盖。清单列出相对路径、文件角色与字节数，未额外计算内容哈希；源 OPJU 的既有 SHA256 只在 WP2 运行记录中保留用于既定溯源。

## 环境与输入

记录环境为 Windows、Python 3.14.7、NumPy 2.5.4、SciPy 1.18.1；Matplotlib 版本见 `environment.json`。新环境按 `requirements.txt` 安装相同版本。输入包括已审阅的 `WP1/03_数据与参数清单.csv`、其上游整理表、Book1 CSV 工作副本和原始 OPJU 只读副本；原始读数不在脚本中重标定。`WP1/02_原始文献比较表.csv` 标明文献实际阅读层级。压缩包中的原始数据副本仅为交接，不改变项目内原件。

## 从冻结输入重建关键证据

在**解压目录的根目录**依次运行以下命令（PowerShell 可将 `python` 换成本机虚拟环境解释器）。每一行是独立命令，不建议在已有科研工作区直接重跑以免覆盖原交付表。

```text
python src/data_ingest/build_wp1_manifest.py
python WP2/src/validate.py
python src/effective_network_model/v6_scale_diagnostics.py
python src/effective_network_model/v6_regime_map.py
python src/effective_network_model/v6_plot_regime.py
python src/effective_network_model/v6_fit_wp4.py M0 --max-nfev 150
python src/effective_network_model/v6_fit_wp4.py M1 --max-nfev 150
python src/effective_network_model/v6_fit_wp4.py M2 --max-nfev 150
python src/effective_network_model/v6_report_wp4.py
python src/effective_network_model/v6_wp5_holdouts.py M0 --max-nfev 120
python src/effective_network_model/v6_wp5_holdouts.py M1 --max-nfev 120
python src/effective_network_model/v6_wp5_holdouts.py M2 --max-nfev 120
python src/effective_network_model/v6_wp5_report.py
python src/effective_network_model/v6_wp5_identifiability.py
python src/effective_network_model/v6_wp5_sensitivity.py --max-nfev 120
python src/effective_network_model/v6_wp5_worst_cases.py
python src/effective_network_model/v6_wp6_figures.py
```

WP5 的 M2 留出搜索可中断后用脚本 `--resume` 接续，再用 `--rebuild-only` 从完整搜索记录重建逐膜表，避免重复行；完整语法见脚本帮助。耗时的 WP3 网格、M2 留出、剖面和 44 情景搜索均已在包内保存原输出和全部起点记录，复现者可先用已保存表核对、再按需重算。最终报表采用 BDF 轨迹，快速特征/矩阵指数只是优化阶段的数值求值；近上界/失败起点保留在 JSON。`v6_wp6_figures.py` 只读已保存表，不重新拟合。

## 可核对的目标与停止边界

- WP2：46 PASS、4 OBSERVED、0 FAIL；OBSERVED 是辨识观察。
- WP4：M0/M1/M2 训练损失 BDF 约 `1.209716/0.667126/0.552434`；M2 主五十格不满足每点 1% 网格差，百格精化另见 `WP4/12_M2百格精化.json`。
- WP5：E_time 约 `0.884752/0.586806/0.780691`，E_h 约 `1.655928/1.501394/1.028292`；2600 nm 留出终点误差倍数约 `38.45/40.52/1.881`。三项代价剖面点未成功的状态应保持可见。
- 六张图位于 `WP6/figures/`。报告只允许据 `WP6/03_论文主张与证据对应表.csv` 的范围作结论；留出为旧数据内部检验，投稿判断仍待导师。

这些数字是给复现者的回归锚点，不代表另行验证过压缩包全部计算。当前实际执行者为本机 Codex；学生本人讲解、独立运行及导师第三次审阅未被本包代替。
