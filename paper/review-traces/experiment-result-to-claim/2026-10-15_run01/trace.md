# Experiment Result-to-Claim Review Trace

## Prompt

对 PE 项目做只读 RESULT-TO-CLAIM 独立审查，不得修改文件。读取项目说明、复现记录、B10/B13 正式总结以及长管、surface PB、finite PNP 输出清单。评估：

- C1 给定候选荷电机制可产生阳离子选择性。
- C2 负电/阳离子选择性必然降低共同 NaCl 通量。
- C3 现有证据已证明真实 PE 的荷电来源或真实 `sigma`。
- C4 当前数值结果可用于设计下一轮判别实验。

每项返回 `claim_supported`、`confidence`、支持/不支持内容、缺失证据、建议改写和下一实验。模型输出与实测必须分开。

## Raw Reviewer Response

### C1

- `claim_supported: yes`，仅限计算存在性/机制可行性。
- `confidence: high`（对所求模型）；真实 PE 外推为 low。
- 两套独立候选化学都产生阳离子偏好，且数值收敛/守恒通过；不证明真实 PE 位点、参数、`sigma` 或实验迁移数。
- 缺 zeta method/cell/背景盐/温度、湿态通路、独立化学约束和跨膜选择性测量。

### C2

- `claim_supported: no`，`confidence: high`。
- 长管 `J/J0=1.0560138` 和五排 `J/J0=1.0275731` 直接否定“必然降低”。
- 结果不能推出负电总会增大真实通量；仍缺同样品的联合测量、空白和质量守恒。

### C3

- `claim_supported: no`，`confidence: high`。
- B10 只证明 zeta 身份和统计口径可复算；真实 cell、背景液、滑移面、温度、安装和独立膜身份未闭合。
- 干态 KPFM 不等于液相 zeta；候选化学参数未经实测，不能确认来源或真实 `sigma`。

### C4

- `claim_supported: partial`，`confidence: high`。
- 当前模型足以设计联合判别逻辑和控制项，但不能冻结真实效应量、参数范围、样本量或阈值。
- 先闭合 B11/B12，再预注册两母膜 O/S/W、联合电动/传输/通路测量和留出验证。

## Review Metadata

- `paper_claim_audit: unavailable`
- `review_status: pending external review`
- Reviewer mode: independent read-only agent
- Files modified by reviewer: none
