# PE 科研项目续跑计划

> 历史任务计划。2026-10-17 导师 v5.0.0 将当前路线改为纯理论计算；下方 Phase 5 与末尾 Next Step 的实验接入安排不再是当前执行指令。现行任务以根目录 `项目进展计划与任务看板.md` 第 6 节为准；断点从同目录 `progress.md` 顶部“最近检查点”恢复。保留原计划供追溯，不覆盖历史。

## 目标

综合两段前期聊天、本项目全部已有说明与计算资产，确认可复现基线和未决科学问题；在不虚构实验参数的前提下，选择证据最充分的下一阶段计算，产出可复现结果、验证记录和阶段结论。

## 阶段

### Phase 1：恢复上下文与建立基线

- Status: complete
- 定位并读取两段指定聊天。
- 阅读项目说明、行动清单、任务看板、复现记录、P0 阶段总结和资料索引。
- 盘点代码、既有输出、原始资料与版本状态。

### Phase 2：复现与诊断

- Status: complete
- 运行现有自动化测试和代表性计算入口。
- 核对数值收敛、守恒、边界条件、单位和输出可读性。
- 将失败或不可靠环节定位到最小责任层。

### Phase 3：确定并实现下一研究增量

- Status: complete
- 根据前期聊天和现有证据选择一个可检验、无需臆造实验参数的下一步。
- 先定义验收量与回归测试，再做最小实现。

### Phase 4：运行、分析与交付

- Status: complete
- 运行新增计算，生成版本化结果。
- 对照基线解释结果，明确哪些结论受支持、哪些仍待实验约束。
- 更新项目看板、复现记录和阶段总结。

### Phase 5：外部输入接入与协议准备

- Status: in_progress
- 完成最终测试、主入口、来源哈希、JSON 和工作树一致性检查。
- 保持 P1 整体为进行中，不把局部完成扩大为材料机制验证。
- B11 已完成 Origin 数值表解码；逐曲线 tracker 等待校准、面积、体积、补液、润湿和样品身份记录。
- B12 等待 zeta method、cell、背景电解质、温度和安装记录。

### Phase 6：GitHub 同步

- Status: complete
- 将本轮代码修复、正式研究记录、声明审查和必要规划记录提交到 `dev-yi`。
- 排除 `原始资料/` 与可再生成的模型输出。
- 推送后核对本地与 `origin/dev-yi` 的提交一致性。

## 决策约束

- 原始资料只读；不覆盖既有原始数据或历史结果。
- 不把假设位点密度、结合常数或表面化学机制表述为实测 PE 性质。
- 每个结论必须能追溯到代码、配置和输出。
- 先复现，后扩展；数值验证优先于增加模型复杂度。
- 旧看板中的 `HOLD` 由本轮用户“继续跑该科研项目”的明确指令解除；这不解除看板对科学结论与外部实验输入的准入限制。

## Errors Encountered

| Error | Attempt | Resolution |
|---|---:|---|
| `read_thread` 的 `turnLimit=50` 超过上限 10 | 1 | 改为 10 |
| 两个目标 thread ID 未被当前 `local` 主机索引 | 2 | 改查本机会话存档与项目文件 |
| `list_threads` / `list_archived_threads` 的 `limit=100` 超过上限 50 | 1 | 改为 50 |
| `git -C PE status` 因 dubious ownership 被拒绝 | 1 | 后续用命令级 `safe.directory` 配置读取，不修改全局配置 |
| 目标聊天在本机返回 `thread not loaded` / `no rollout found` | 3 | 认定当前账户无法读取；以同伴留下的正式文件为证据源 |
| 项目 `.venv` 不存在 | 1 | 检查系统 Python 3.14 的现有科学库，不自动安装 |
| Codex 捆绑 Python 缺少 SciPy | 1 | 不用于 v3；改查系统 Python 3.14 |
| 长管诊断的只读 PB 副本哈希测试失败 | 1 | 根因确认：Git `core.autocrlf=true` 将 339 个 LF 转为 CRLF；未改科学内容，待设计最小跨平台修复 |
| 介电对照图底部注释与 x 轴标题重叠 | 1 | 已在视觉检查中确认；待对照原图和绘图代码定位版式来源 |
| 当前 PowerShell 无 `python` 命令别名 | 1 | 统一改用 Windows 启动器 `py -3.14`，不重复失败入口 |
| 首次 `checkout-index` 仍输出 CRLF | 1 | 未提交属性尚未参与索引过滤；改为单命令 `-c core.autocrlf=false` 直接恢复索引 LF，不改全局配置 |
| 命令级关闭 `core.autocrlf` 后仍输出 CRLF | 2 | `checkout-index` 的属性来源仍是索引；第三次改为临时暂存 `.gitattributes`、恢复副本、立即取消暂存 |
| 临时暂存属性后 `checkout-index` 仍输出 CRLF | 3 | 停止使用该机制；改由 `git archive` 导出提交中的原始 blob 并原样解包，校验后删除临时归档 |
| `git archive` 也导出 14605 字节 CRLF | 4 | 提交 blob 为 14266 字节；将路径规则从 `-text` 改为显式 `text eol=lf`，确保 Windows 工作树固定 LF |
| 误读 `src/surface_pb/outputs/run_manifest.json` | 1 | 该模块按现有约定输出到共享 `src/outputs/`；改读正确路径，不重复失败入口 |
| 假定 surface 清单名为 `run_manifest.json` | 2 | 源码确认实际名称为 `src/outputs/surface_manifest.json`；后续以写出代码定位文件名，不再猜测 |
| JSON 汇总命令在 `foreach` 后直接接管道导致 PowerShell `empty pipe element` | 1 | 改用显式 `$results` 数组后再格式化，不重复原命令 |
| 五排 refinement JSON 被误假定含顶层 `results` 键 | 1 | 先读取实际顶层键，再按真实结构提取 32 网格记录 |
| 首次提交因当前仓库未配置 Git 作者身份而失败 | 1 | 采用最近主线提交的作者 `suzyzhao27-bot <suzyzhao27@gmail.com>`，仅写入仓库本地配置 |

## Next Step

接收并只读核验 B11/B12 外部实验资料，按 tracker 逐条登记可定位证据；同时继续完善不依赖未知参数的协议与质量门。在资料到位前，不启动依赖真实电导曲线或 zeta 测量条件的参数反演与机制判定。
