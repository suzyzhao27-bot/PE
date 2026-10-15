# Findings

## 2026-10-15 上下文恢复

- 项目主题：超薄 PE 纤维膜的界面荷电、离子选择性与盐输运。
- 共享计算代码分为四条链：`v2_training`、`long_tube_diagnostic`、`surface_pb`、`finite_pnp`。
- README 明确：长管 PB-NP 仅是受限诊断，不是有限长度真实膜；二维有限 PNP 是主模型；表面位点密度和结合参数仍是假设，不能当作实测 PE 性质。
- 当前工作区根目录不是 Git 仓库，`PE/` 是仓库；读取 Git 状态时触发 Windows dubious ownership，需要命令级安全目录处理。
- 两个用户提供的 thread ID 为 `01a13ebd-2f48-7940-a7f1-43494bcdd454` 与 `01a13a17-91d3-72a1-af66-83fc6f9b3909`。应用索引未找到，需查本机会话存档。
- 本机 `sessions/` 与 `archived_sessions/` 的文件名检索没有找到两个目标 ID；当前 `.codex` 文本检索只命中本次聊天对这两个链接的引用。
- 应用数据目录检索仅在 Codex 桌面日志与 Sentry scope 中命中这两个 ID，尚未发现目标聊天的消息正文或本地实体；这与“同伴账户私有链接、当前账户不可读”一致，仍需用日志上下文作最后排除。
- P0 已通过：v2 基础模型 10/10、长管诊断 7/7、外纤维 PB 9/9、有限二维 PNP 10/10；`tests_pipeline.py` 的 6 项因受控终端临时目录 ACL 未执行成功，不能写作全部通过。
- 主五排模型在假设 `Gamma=0.05 nm^-2`、`K_Na=0`、`K_Cl=10000 M^-1` 下给出 `t_Na=0.65283`、共同盐通量比 `1.02757`、理论液相电位 `-10.431 mV`；它支持“选择性偏好可形成”，不支持“负电必然挡盐”。
- 原项目路线规定下一阶段是 P1 源数据闭合与不确定性台账。现阶段允许推进 B10（36 ramp/9 zeta 点身份与不确定性表）；B11/B12 依赖 Origin 导出和仪器记录，仍是外部条件。
- 原看板 `HOLD` 是历史暂停状态。本轮用户已明确要求继续运行，视为恢复项目执行，但不能越过各阶段的科学准入条件。
- 原始资料中已确认存在正式 v2/v3 核心 PDF、学生交接说明、数据核查报告、`observations.csv`、`zeta_summary.csv`、PB/PNP README 及相应历史副本。阅读时以 `01_研究框架与资料_v2` 和 `02_荷电与选择性实算_v3` 为正式版本，历史副本仅用于冲突追溯。
- 原始资料共 1133 个文件、约 670.8 MB；正式 v2/v3 之外的大量历史副本不应作为日常入口。
- `gaps.csv` 登记 13 个缺口，其中 G01-G08、G10、G13 等关键项仍未闭合；这阻止旧盐曲线机制拟合，但不阻止 zeta 描述性审计。
- zeta 审计为 36 条 ramp、9 个汇总点、504 个非空原格；原表 SD 精确对应四次 ramp 的总体 SD（ddof=0），独立膜数量仍为 `unknown`。
- 三个 IEP 仅由 pH 约 3 与 7 的均值线性插值获得：4.1381、4.0532、3.9964；不能表述为直接密集扫描或高精度拟合。
- 当前仓库没有 `.venv` 和既有运行输出。系统有 Python 3.14.3；Codex 捆绑 Python 缺 SciPy，不能用于 v3。
- 系统 Python 3.14.3 已有 NumPy 2.4.4、SciPy 1.17.1、Matplotlib 3.10.8。
- 当前复跑：`v2_training` 16/16、`surface_pb` 9/9、`finite_pnp` 10/10 通过；长管模块的 6 个物理/数值测试通过，唯一失败是字节身份哈希。
- 长管哈希失败的根因已完整证明：当前文件 14605 字节，含 339 个 CRLF；将 CRLF 在内存中恢复为 LF 后为 14266 字节，SHA256 恰为登记值 `3c681f...`。原始正式副本与历史副本也均为该哈希。Git 报告索引 LF、工作树 CRLF，且 `core.autocrlf=true`。
- 同样的换行转换也发生在 `src/v2_training/models.py`；因此问题属于 Windows checkout 字节规范，不是模型内容漂移。
- 主计算已在 Python 3.14.3 / NumPy 2.4.4 / SciPy 1.17.1 下完成：v2 demo（24 PB、54 flux 情景）、synthetic scenario、surface PB（100 情景）、有限二维 PNP（10 transport、12 G 情景）、五排初版与 16/24/32 细化均成功。
- 二维主模型最大缩放残差 `1.4634e-10`、最大相对净电流 `1.6409e-15`、最大离子通量守恒误差 `4.0741e-14`。
- 五排 32 cells/fiber 结果：`J_chemical/J_neutral=1.0275731111`、`t_Na=0.6528317904`、静态平均 `sigma=-0.0047360405 C/m2`、开路平均 `sigma=-0.0077769197 C/m2`、`Delta V=-10.43103665 mV`；与项目历史记录显示值一致。
- 绘图脚本均成功：surface 2 张 PNG、有限单层场图 1 张、五排场图 2 张。
- 视觉检查：`surface_regulation_and_selectivity.png`、`finite_PNP_fields.png` 和两张五排场图无明显裁切或互相遮挡；`dielectric_symmetric_control.png` 的底部说明文字与 x 轴标题明显重叠，需要判定是历史脚本缺陷还是环境回归。
- 当前仓库无 `.gitattributes`；因此 `core.autocrlf=true` 会把受 SHA256 登记约束的只读 PB 副本检出为 CRLF。最小修复应只对该字节身份文件声明 `-text`，避免把全仓库换行策略一次性改写。
- 介电对照图由 `src/surface_pb/plot_surface.py` 单行生成：`constrained_layout=True` 只管理坐标轴，不为 `fig.text(..., y=.01)` 预留稳定底边，因此底注与 x 轴标题重叠。修复应改为明确底部留白和位置，并加可自动检查的版式验收。
- 修复后 `dielectric_symmetric_control.png` 的底注与 x 轴标题有清晰可见间隔；自动测试确认显示坐标中的垂直间距大于 4 px。`surface_regulation_and_selectivity.png` 保持三联图内容完整，无标题、图例或坐标文字裁切。
- 修复后完整自动化测试为 43/43：v2 16、长管 7、surface 10（含新版式回归）、finite PNP 10。长管主入口已通过来源哈希门，清单中的只读副本哈希为 `3c681f...` 且 `copy_hash_matches_registered=true`。
- 三份 B10 正式只读输入在核查前后 SHA256 完全一致；独立内存核查为 36 observations、9 groups、3 hashes、0 blanks、0 duplicate source rows，每组 4 个唯一 ramp，均值/SD 与汇总表仅有浮点舍入级差异。
- PB 副本的 `git hash-object` 与 `HEAD` blob 均为 `a8a54ab...`；暂时的 Git `M` 状态不是文件内容差异，需在属性规则进入索引或刷新索引后复核。
- 布局重构后 `dielectric_image_control.csv` 与正式 v3 输出 SHA256 同为 `13f9976...`。`surface_scenarios.csv` 和 `surface_radial_profiles.csv` 的行数、列名、文本字段与正式输出一致；前者最大绝对数值差 `7.84e-10` 出现在约 229 的分配比（约 `3.4e-12` 相对），后者最大相对差 `2.54e-13`。场景表表面最大相对差 1 仅是 Gauss error 的 `0` 对 `1.52e-13`。

## 待验证

- B11 所需的 Origin 时间-电导表、校准、面积、体积、补液和润湿记录尚未取得。
- B12 所需的 zeta method、cell、背景电解质、温度和安装记录尚未取得。
- 独立膜数量仍为 `unknown`，不能从四次 ramp 反推生物学/材料学独立样本数。

## 2026-10-15 最终复核

- 新一轮完整自动化测试为 43/43：v2 16、长管 7、surface 10、finite PNP 10。
- `run_diagnostic.py` 再次成功运行；长管结果仍明确限定为假想长管、局部径向平衡诊断，不是有限长度真实 PE 验证。
- 五排 32 网格五个关键值与既有基线逐位一致：通量比 `1.0275731111161168`、`t_Na=0.6528317904475466`、静态平均 `sigma=-0.004736040481486473 C/m2`、开路平均 `sigma=-0.007776919660934164 C/m2`、液相电位差 `-10.431036649323906 mV`。
- 三份正式 B10 输入保持只读且 SHA256 与登记值一致：`118773dc...`、`2835ac75...`、`7e86a62b...`。
