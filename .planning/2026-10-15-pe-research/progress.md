# Progress

## 2026-10-15

- 已读取 `using-superpowers`、`scientific-toolkit-skill` 和 `planning-with-files` 工作流说明。
- 已尝试通过应用读取两段目标聊天；当前主机索引未命中。
- 已读取项目 README 和顶层文件清单，确认四条计算链及其用途。
- 当前阶段：Phase 1（上下文恢复与基线建立）。
- 已读取六份核心项目文档；确认 P0 通过、P1/B10 为当前可推进卡片，B11/B12 受外部实验记录阻塞。
- 已确认 Git 分支 `dev-yi` 相对 `origin/dev-yi` ahead 2；本轮仅新增三份规划记录，未改动既有研究代码。
- 已在原始资料中定位正式交接、源数据审计和 v2/v3 报告；已确认目标聊天未出现在本机会话文件索引中。
- 已将三份过程记录迁入隐藏的 `.planning/2026-10-15-pe-research/`，根目录恢复整洁，并将 `.planning/` 加入 Git 忽略。
- 已完成 Phase 1：目标聊天不可读的限制已确认，正式文件与原始审计成为当前证据基线。
- 已开始 Phase 2：发现项目 `.venv` 缺失、Codex 捆绑 Python 缺 SciPy，待检查系统 Python 3.14。
- 已确认系统 Python 3.14 科学环境并复跑 45 项测试：44 通过，1 项字节身份哈希失败；根因是 Git 自动 CRLF 转换，6 个长管物理/数值断言本身均通过。
- 已完成 v2 demo/scenario、surface PB、有限二维 PNP、五排初版与 16/24/32 网格细化，并生成 5 张图；关键结果复现历史显示值。
- 已恢复本轮上下文并核对仓库：仅 `.gitignore` 有未提交改动，`.planning/` 已被忽略；根目录没有新增过程文件。
- 已定位两个待修复点：`pb_v1_readonly_copy.py` 的工作树 CRLF 导致登记哈希失败；`plot_surface.py` 的底部 `fig.text` 与 x 轴标题重叠。
- 当前等待最小修复设计确认，尚未改动科学代码、登记哈希或历史结果。
- 实施入口第一次调用 `python` 失败，因为当前 PowerShell 未注册该别名；后续统一使用 `py -3.14`。
- 添加精确 `-text` 属性后，首次普通 `checkout-index` 仍受索引阶段的 `core.autocrlf=true` 影响；已改用命令级关闭转换恢复该工作副本。
- 命令级 `core.autocrlf=false` 仍未改变 `checkout-index` 输出，证明其属性解析依赖索引；改为临时暂存属性文件后再恢复，并在完成后取消暂存。
- 临时暂存属性仍不能改变本机 `checkout-index` 的 CRLF 输出；已停止重复该路径，改用 `git archive` 从提交对象导出原始字节。
- `git cat-file -s` 证明提交 blob 为 14266 字节，而 `git archive` 条目为 14605 字节；因此将属性策略修正为精确 `text eol=lf`，不再依赖 `-text`。
- Task 1 完成：只读 PB 工作副本为 `i/lf w/lf attr/text eol=lf`，SHA256 为登记值 `3c681f...`，长管 7/7 测试通过；临时归档已删除。
- Task 2 完成：新增布局回归测试先以缺少构图函数失败，提取函数并预留底边后 10/10 surface 相关测试通过；两张 surface 图已重新生成并视觉复核。
- Task 3 测试门完成：v2 16/16、长管 7/7、surface 10/10、finite PNP 10/10，共 43/43 通过。
- 长管诊断主入口成功；`run_manifest.json` 登记 `copy_hash_matches_registered=true`，Python 3.14.3，且明确保持“受限长管诊断、非真实 PE 验证”的范围。
- Task 3 完成：五排 32 网格五项主数值与既有基线一致，未因可复现性和绘图修复改变科学结论。
- Task 4 完成：B10 正式总结、项目看板和复现记录已更新；原始输入读取前后 SHA256 一致，B11/B12 继续阻塞。
- 建立五类变量证据台账前读取模型清单时误用 surface 模块内输出路径；已确认正确位置为 `src/outputs/`。
- 进一步确认 surface 清单实际名为 `surface_manifest.json`；目录与生成结果本身正常。
- 与正式 v3 输出比较：介电控制 CSV SHA256 完全一致；surface 场景/径向 CSV 行列与文本字段一致，数值差仅为浮点求解噪声。
- 最终复核再次通过 43/43：v2 16、长管 7、surface 10、finite PNP 10。
- 长管诊断主入口再次成功；五排 32 网格五个主数值与既有基线逐位一致。
- 三份正式 B10 原始输入仅作读取，完整 SHA256 仍为 `118773dc16e45567dc68169c84b48115280a5bc2222cd0044aa8c993ba6a3272`、`2835ac752e755bb3c59ceeb48661f286bdb695ca512b4d8dc0c78f68cbcbb064`、`7e86a62bb33b37eeb4e51bf1c3835f1403464c3b19c28d9897b5167248b2d314`。
- 首次 JSON 汇总命令因 PowerShell 在 `foreach` 后直接接管道而解析失败；这是验证命令语法错误，已改用显式结果数组重试。
- JSON 重试通过：长管 `run_manifest.json`、五排 refinement JSON 和声明审查 `verdict.json` 均可解析；清单确认 `copy_hash_matches_registered=true`。
- Task 5 完成：`git diff --check` 返回 0（仅有既有换行转换提示），工作树文件均属于已批准的可复现修复、正式研究记录和 P1 交付；根目录只保留正式入口与控制文件，过程记录集中在 `.planning/`。
- 最终结构化提取五排数值时误用顶层 `results` 键；JSON 文件本身解析正常，改为读取实际键名后完成提取。
- 按实际 `grid_study` 键提取 32 网格记录成功；五个主数值与既有基线严格相等，`exact_baseline_match=True`。
- 用户要求将新完成工作同步至 GitHub；已解除 `.planning/` 的忽略规则，使四份必要规划记录进入版本控制。
- 同步前重新运行四组完整测试，结果仍为 43/43；三份正式 B10 原始输入 SHA256 保持登记值，待提交列表不含 `原始资料/` 与生成输出。
- 首次 `git commit` 因仓库未配置作者身份而失败，未产生提交；已从最近主线提交确认作者身份，准备仅设置仓库本地配置后重试。

### 下一步

- 提交并推送 `dev-yi`，随后核对 `HEAD` 与 `origin/dev-yi`。
