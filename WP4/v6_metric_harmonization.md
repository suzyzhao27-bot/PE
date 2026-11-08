# v6 与 v7 留出误差的同定义核对（2026-10-23）

v6 的正式表把六膜的逐膜 asinh RMSE 作算术平均；v7 的 `E_time` 和 `E_h` 是六膜等权均方根：

\[
E_{\mathrm{v7}}=\sqrt{\frac{1}{6}\sum_{h=1}^{6}E_h^2},\qquad
E_h=\sqrt{\frac{1}{n_h}\sum_{t\in\mathrm{test}(h)}
[\operatorname{asinh}(\hat y_{ht}/1)-\operatorname{asinh}(y_{ht}/1)]^2}.
\]

两者不能直接并列作细微优劣。脚本 `src/effective_network_model/v7_audit_v6_metrics.py` 仅读取14号包内冻结的 v6 WP5 逐点 CSV，逐膜复算并与 v6 存档 JSON 的六个逐膜值和算术平均逐一断言一致；没有重拟合旧模型。复算结果在 `WP4/v6_metric_harmonization.csv`：

| 历史模型 | v6 原 E_time / E_h（算术平均） | 按 v7 定义 E_time / E_h（均方根） |
|---|---:|---:|
| M0 | 0.884752 / 1.655928 | 0.977634 / 2.276692 |
| M1 | 0.586806 / 1.501394 | 0.620206 / 2.206858 |
| M2 | 0.780691 / 1.028292 | 0.878072 / 1.125299 |

这只统一**误差汇总定义**；旧模型的物理假设、几何和技术状态未被 v7 重新验证。尤其旧 M2 有数值预测失败及局部格点误差记录，保持“历史技术条件性”，不能据其较低 `E_h` 认证新孔网或细微模型优劣。正式模型比较仍需 NS/ND 六折全部结束后并列原单位残差、最差膜、2600 nm 端点比及所有失败起点。

运行命令：`.\.venv\Scripts\python.exe src\effective_network_model\v7_audit_v6_metrics.py`，退出码 0；每个模型的时间/整膜逐点数分别为 18/36，六膜逐项断言均通过。数据是原 Book1 六膜的回顾性留出，没有独立新实验。
