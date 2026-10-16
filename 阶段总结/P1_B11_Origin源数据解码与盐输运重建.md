# P1/B11 Origin 源数据解码与盐输运重建

**工作包：** B11 - 目标 OPJU 导出为带样品身份的 CSV  
**执行日期：** 2026-10-16  
**状态：** 进行中（表格解码完成；实验记录闭合未完成）  
**原始资料权限：** 只读；全部源文件提取前后 SHA256 一致

## 1. 本轮目标与必要性

本轮不再从图片取点，而是直接读取 Origin 工程中的工作表、列标签、单位、注释和数值。该步骤用于关闭 G01“盐时间-电导原表”中的表格缺口，并判断 G11“长时曲线身份”能否由源工程本身解释。电导-浓度校准、有效面积、储液池体积、取样/补液、温度、润湿历史和独立膜身份仍须由实验记录提供，不能由源表结构推断。

## 2. 输入与只读保证

读取 `source_snapshots/` 下 5 个唯一 OPJU：

- `AM_nanofiltration_origin.opju`
- `Polymer_KPFM_origin.opju`
- `Polymer_salt_2h.opju`
- `Polymer_salt_dry_6h.opju`
- `Polymer_zeta_origin.opju`

提取使用独立的 Origin COM 会话；脚本没有调用 `Save`。每个工程在会话退出并释放文件句柄后重新计算 SHA256。5/5 个工程的提取前后哈希一致，完整记录见 `data/processed/origin_exports/manifest.json`。

## 3. 实际产出

| 产出 | 内容 | 用途 |
|---|---|---|
| `src/data_ingest/extract_origin_opju.ps1` | 枚举工程内全部工作表，导出 CSV、列元数据和哈希清单 | 可重复的 OPJU 解码 |
| `data/processed/origin_exports/` | 32 张工作表的 CSV/metadata，加项目级 manifest，共 70 个文件 | 保留源表结构与派生拟合表 |
| `src/data_ingest/normalize_salt_transport.py` | 对盐输运原始表建立显式列映射和证据标记 | 生成长表，不混入拟合表 |
| `data/processed/salt_transport/salt_transport_tidy.csv` | 17 条曲线、167 个有效时间-电导点 | 描述性 QC 输入 |
| `data/processed/salt_transport/normalization_audit.json` | 每条曲线的样品、时间列、单位证据与拟合禁入原因 | 不确定性台账 |
| `data/processed/salt_transport/exact_duplicate_curves.csv` | 3 对完全重复曲线 | 防止把复制曲线当独立重复 |
| `data/processed/salt_transport/salt_transport_curve_summary.csv` | 曲线范围、起终值、描述性 OLS、单调性检查 | QC，不用于机制拟合 |
| `data/processed/figures/salt_transport_origin_qc.png` | 2 h 与 dry/unwetted 源表重建图 | 视觉核查，不代表浓度、通量或渗透率 |

## 4. 工程级审计结果

| 工程 | 工作表 | 判定 |
|---|---:|---|
| AM nanofiltration | 2 | 文献/本工作表观量对比，不是时间-电导原表 |
| Polymer KPFM | 1 | 4 点表格，可追溯但本轮不把干态电势换算为液相电荷 |
| Polymer salt 2 h | 2 | 恢复 14 条曲线；显式或表内继承的时间/电导列和样品注释可用 |
| Polymer salt dry 6 h | 25 | 仅“未润湿体系数据”为原始表；其余 24 张为线性/多项式拟合及拟合曲线，已排除出原始观测长表 |
| Polymer zeta | 2 | 恢复 draw-ratio 汇总和 3 组 pH × 3 draw-ratio 表；method/cell 仍未出现 |

## 5. 盐输运关键发现

1. 共建立 17 条曲线、167 个非空时间-电导点：2 h 工程 14 条，dry/unwetted 工程 3 条。
2. `Polymer_salt_dry_6h.opju` 中最长时间为 `2831 min = 47.18 h`。因此文件名中的 `6h` 不能解释为表内全部曲线的实际最长时长；源工程本身确认了长时数据存在，但它与公开 Fig.15c 各图线的一一映射仍需记录确认。
3. Book1 与 Book2 的 `65 nm`、`700 nm`、`2600 nm` 三条曲线在所有时间和电导数值上完全相同，属于复制关系，不得计为独立实验重复。
4. 2 h 表中 `Time (min)` 与 `Ionic conductivity (μS/cm)` 对部分列为 Origin 显式标签；其余列按同表相邻/共享时间列继承。dry/unwetted 表仅有 `700/350/200 nm` 注释，时间和电导单位来自伴随盐工程继承，必须在原实验记录中确认。
5. 当前图和 OLS 仅用于检查转录、时间范围、数量级与异常；未进行电导到浓度、通量或渗透率换算。

## 6. 质量检查

- OPJU 源文件完整性：5/5 提取前后 SHA256 一致。
- 数据接入测试：4/4 通过，包括哈希、点数、长时文件名冲突和重复曲线检测。
- 原有科学计算基线：43/43 通过（v2 16、长管 7、surface 10、finite PNP 10）。v2 首次在受限沙箱内因 Python 临时目录 ACL 失败；在正常临时目录权限下重跑 16/16 通过，未修改科学代码。
- 图形人工检查：QC 图可读，使用对数纵轴；标题明确声明不是浓度、通量或渗透率。

## 7. 未闭合项与准入决定

以下信息仍缺失：电导-浓度校准与温度补偿、blank、有效膜面积与密封、两池体积、取样/补液、搅拌、曲线级润湿历史、母膜/裁片/曲线身份和 dry 表单位确认。

因此：

- **B11 从“阻塞”转为“进行中”**：OPJU 数值表已可重复解码，G01 的表格部分已关闭。
- **阶段 1 仍未通过**：缺失的实验记录使这些曲线仍不具备机制拟合资格。
- **允许继续做**：描述性 QC、重复/异常检查、实验记录追索和 B20 协议草案。
- **禁止继续做**：把电导直接当浓度，计算通量/渗透率，或用旧曲线反演 `sigma`、润湿分数与孔网参数。
