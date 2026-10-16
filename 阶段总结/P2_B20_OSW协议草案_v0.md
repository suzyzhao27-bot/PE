# P2/B20 O/S/W 最小来源试点协议草案 v0

**日期：** 2026-10-16  
**状态：** 草案，不可直接开测，不代表 B20 已冻结  
**目标：** 用两张独立母膜的同膜随机分片，检验接触/洗脱历史是否可重复改变电动响应、离子选择性和湿态通路；不预设效应方向，不把洗不掉等同“本征 PE”。

## 1. 草案依据

本草案整理自以下只读材料：

- `原始资料/赵珈绮/01_研究框架与资料_v2/01_研究报告/第二课题_完整研究报告_v2.txt`
- `原始资料/赵珈绮/01_研究框架与资料_v2/06_独立审查与新近邻/最小来源干预链_适用范围与验收.txt`
- `原始资料/赵珈绮/01_研究框架与资料_v2/08_学生执行与导师关卡/02_导师首周须闭合的协议.txt`
- `项目进展计划与任务看板.md`
- `阶段总结/P1_B11_Origin源数据解码与盐输运重建.md`

## 2. 最小设计

### 2.1 实验单位与层级

| 层级 | 最小要求 | 统计/解释边界 |
|---|---|---|
| 母膜 | 2 张独立制备母膜，记为 M1、M2 | 母膜才是制备重复的最低候选单位 |
| 处理分片 | 每张母膜随机分到 O、S、W | 同母膜比较用于控制片间制备差异 |
| 伴随片 | 按仪器需求从同一母膜/邻近位置取得 | 必须记录与输运片的对应关系；不能默认可互换 |
| 技术读数 | 同片重复读数或压力 ramp | 估计测量重复性，不增加独立膜数量 |

### 2.2 三组处理

| 组别 | 操作定义 | 允许解释 | 禁止解释 |
|---|---|---|---|
| O | 原始状态；只接受开测所必需且预先冻结的共同平衡步骤 | 当前批次基线 | 不自动代表“纯 PE”或未接触任何加工助剂 |
| S | 与 W 等时、同温、同搅拌、同初始协议液体，但不换新液 | 接触/老化历史对照 | S 也可能发生抽提，不能称完美假洗 |
| W | 与 S 使用同种兼容协议液体，在冻结的节点换新液 | 相对 S 的“加强抽提/交换历史”效应 | W−S 不能唯一指认具体分子或电子起源 |

所有组在处理后进入同一终末测量液、同一平衡温度和时间。处理液种类、面积/体积比、总时长、换液次数/时间、搅拌方式和终末 NaCl 平衡必须由导师在开测前填写；本草案不补默认值。

## 3. 随机化与身份登记

1. 为每张母膜绘制位置图，固定坐标、机器方向/横纵向和边缘区。
2. 裁片后立即赋不可重用 ID：`批次-母膜-位置-用途-处理-序号`。
3. 在不知道后续测量结果时生成 O/S/W 分配；随机种子、执行人和分配时间写入记录。
4. 仪器测量顺序另行随机化或平衡；记录夹持、正反面、首测时间和顺序。
5. 同一片的顺序测量只算时间序列；须有同等老化/顺序的平行片，不能替代母膜重复。

## 4. 开测前必须冻结的记录字段

### 4.1 样品与处理

`batch_id, parent_membrane_id, coupon_id, position_x, position_y, machine_direction, intended_use, treatment_group, protocol_liquid_batch, initial_liquid_volume, membrane_exposed_area, area_volume_ratio, start_time, end_time, temperature, agitation, change_count, change_times, final_measurement_liquid, equilibration_start, equilibration_end, operator`

### 4.2 盐输运装置

`run_id, coupon_id, fixture_id, exposed_diameter, effective_area, seal_material, seal_lot, side_left_initial_volume, side_right_initial_volume, side_left_initial_solution, side_right_initial_solution, actual_mass_or_volume_preparation, bath_temperature, stirring_left, stirring_right, sampling_time, sampling_volume, returned_or_replenished, replenishment_solution, conductivity_raw, conductivity_temperature, calibration_id, blank_id, endpoint_Na_method, endpoint_Cl_method, pre_run_mass, post_run_mass, leakage_flag`

### 4.3 电动/电学与通路伴随读数

`measurement_id, coupon_id_or_companion_id, parent_membrane_id, instrument, method_file, cell, mounting_orientation, electrolyte, concentration_definition, pH_actual, temperature, electrode_or_salt_bridge, sequence_index, raw_file, calibration_id, operator, exclusion_flag, exclusion_reason`

任何字段若暂时无法得到，写 `unknown`；不得留空后由模型反推。

## 5. 必须随主样运行的空白与核对

| 项目 | 目的 | 最低记录 |
|---|---|---|
| 无膜全过程空白 | 限定装置、夹具、液体和操作带来的响应 | 与主样相同液体、时间、温度、搅拌、取样/补液 |
| 不透膜空白 | 检查密封、旁路、外液层与泄漏 | 不透材料、有效面积、压紧/密封方式 |
| 液体/容器空白 | 检查浸出、蒸发和容器污染 | 容器批次、液体批次、质量/体积前后 |
| 电极/盐桥/液接空白 | 分离电位测量基线 | 对称液、交换左右两侧、漂移与偏置 |
| 搅拌对照 | 判断边界层和混合敏感性 | 转速/方式、启动时间和左右一致性 |
| 少量独立 Na/Cl 或中性示踪 | 核对电导是否主要代表目标盐及湿态通路 | 方法、检出限、自由溶液/吸附空白 |

无膜空白只能给装置背景上限，不能唯一反演有膜运行时的外层厚度。

## 6. 校准与质量闭合

开测前冻结：

- 电导标准液浓度定义、原始电导、温度、blank、拟合式、残差和适用区间；
- 左右池初始与每次取样后的实际体积；返还/补液的溶液身份与体积；
- 有效面积的图纸/照片/实测依据和密封压紧记录；
- 独立端点 Na/Cl 或示踪的回收率；
- 体积与盐量闭合公式、容许误差和超界处置。

在这些字段冻结前，电导曲线只能保持原始单位，不转换为浓度、通量或渗透率。

## 7. 坏样与排除规则框架

阈值须用空白和仪器分辨率确定，不能在看到 O/S/W 差异后修改。至少预先定义：

1. 可见破损、折皱、密封越界或旁路泄漏；
2. 体积/盐量闭合失败；
3. 温度或搅拌偏离冻结范围；
4. 校准超出适用区间或校准漂移超界；
5. 电极/盐桥/液接空白超界；
6. 终末平衡未达到预定稳定判据；
7. 原始文件、样品 ID 或处理日志缺失；
8. 伴随片与输运片来源映射断裂。

排除决定必须保留原始行、理由、时间和决定人；不得删除失败数据痕迹。

## 8. 试点结果决策树

- 若 W−S 在 M1、M2 方向一致，且通路/空白/质量闭合合格：进入“回加是否值得触发”的导师关卡；仍只称加强抽提历史效应。
- 若 O、S、W 的变化与湿态通路或装置空白同步：优先归入润湿/结构/装置竞争解释，不用 `sigma` 补偿。
- 若母膜间方向不一致：扩大原因审计或停止来源结论；技术重复不能替代新母膜。
- 若处理液无可检出化学信息：可以报告操作性可去除贡献，但不能指认 AO、凡士林、氧化物或其他具体来源。
- 若效应小于漂移与模型误差包络：结论为当前条件下不可分辨，不升级为“无效”或“本征”。

## 9. 冻结前的明确阻塞项

以下项目全部签字/登记前，本文件保持 v0 草案，B20 不标记完成：

- 当前膜批次来源、配方、抽提、存放和既往润湿历史；
- 兼容协议液体及兼容性/形貌/通路验证；
- 面积/体积比、时间、温度、搅拌、换液和终末平衡；
- 可用仪器、cell、方法页、分辨率与样品尺寸；
- 校准、空白、坏样阈值、盲法/随机化执行人；
- 母膜、伴随片和输运片的实际数量与映射；
- 数据文件命名、备份、原始记录签字和不可修改副本位置。

## 10. 可执行记录接口

`src/data_ingest/osw_sample_registry_template.csv` 固定了母膜、裁片、位置、O/S/W 分组、随机种子、处理条件、终末平衡、输入等级和来源记录字段。模板只含表头，不含虚构样品。

登记表同时要求协议液体身份/批次、兼容性记录 ID、裁片长宽/单位和实际仪器或夹具 ID。缺任一项的记录不能标记 ready，防止分组完成但样品无法装夹或处理液破坏膜结构。

`src/data_ingest/liquid_compatibility_template.csv` 和 `validate_liquid_compatibility.py` 记录具体液体/批次、暴露条件、裁片、质量、尺寸、视觉损伤、润湿观察、原始文件和批准决定。`compatible` 只适用于已登记条件，不可跨液体、批次或暴露条件外推。

`src/data_ingest/validate_osw_protocol.py` 检查字段结构、裁片 ID 唯一性、两张母膜、每张母膜的 O/S/W 覆盖、`unknown`、排除理由和 ready 状态。结构完整仍不自动授权开测；兼容性、仪器、空白和签字审查继续保留。

`src/data_ingest/osw_randomization.py` 对每张母膜的真实裁片 ID 使用研究者登记的随机种子生成可复现、近似平衡的 O/S/W 分配。函数不读取测量结果，不自动生成样品 ID；少于 3 个裁片、空 ID 或重复 ID 直接失败。若裁片数不是 3 的整数倍，各组数量最多相差 1。
