# v7 G 静态极限模型适用图执行记录

权威范围为v7任务包04与只读configs/frozen_spec.json：单一NS诊断候选；H/S各12个65–2600 nm对数等距厚度、β=0–6的9点、主种子7001–7008，共1728个稳态组合。每点不重拟合；α取既有全B候选主值，但稳态P不依赖池容量α。图读及文字锚点保留原身份，未添加组合扫描。

实际命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_g_regime.py compute
    .\.venv\Scripts\python.exe src\effective_network_model\v7_g_regime.py plot

两命令退出0；第二条在图中文字明确标出所有网格8/8种子贯通后定向重绘一次，没有重算1728个组合。生成192张几何图（每个系列12厚度×8种子），最大体积相对误差1.96×10⁻¹⁶。1728行(series,h,β,seed)均唯一；所有图至少存在贯通节点。定向四角图的P与既有fixed_bath_p一致；全网格稳态能量平衡最大相对误差1.58405×10⁻¹²。P总范围9.35831×10⁻¹²–1.78548×10⁻³ m/s。

固定池边界下P=J/A。双端口access耗散比例定义为Σ_port J_e²/[4D_b a_e]除以Σ_all J_e Δu_e，只算既有access阻力，不把全部端口通道阻力混入。八种子单格平均的范围：H图P=1.19056×10⁻¹¹–1.67016×10⁻³ m/s、access比例7.74458×10⁻⁹–0.219082；S图P=1.85979×10⁻¹⁰–1.62975×10⁻³ m/s、access比例1.40525×10⁻⁷–0.253618。两系列12个厚度上，平均P均随β的九档增大而严格下降；全网格8/8种子有贯通节点，未见连接断裂区间。高β时access耗散比例很小，主要电阻在内通道势垒；这些均是当前代理孔网的数值行为，不是材料特定机制证明。

图源为data/processed/figures/v7_G_static_regime.csv；几何审计为WP5/v7_G_geometry_audit.csv，摘要为WP5/v7_G_summary.json。两张图分别为data/processed/figures/v7_G_H_static_regime.svg/.png与v7_G_S_static_regime.svg/.png。颜色为log10 P，白色等值线为双端口access耗散比例，虚线为现有结构锚点厚度；标题/图内说明8/8种子贯通。由于图未揭示额外材料特定结论，按04将其作为拟议补充材料图，不作为新增投稿核心证据；现有2600 nm和L晚期失配仍需在论文结论中保留。
