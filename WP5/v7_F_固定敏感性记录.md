# v7 F 固定敏感性执行记录

按v7任务包04和configs/frozen_spec.json的18项清单执行；每项只改列明量，先用主参数前向，再对单一NS候选作B全数据三起点重拟合。每起点实际目标函数调用上限200。所有结果均是已见全B数据的敏感性，不能称留出。

## S18：仅ND的初始f0=0.5，不适用于单一NS诊断候选

04的D2判定曾条件启用ND作为比较模型；D3后因ND整膜三折多起点不稳定，第5步按冻结规则选NS为单一诊断候选，见WP4/dynamic_branch_decision.md及阶段总结/20261023_v7_第4步模型比较.md。S18的f0仅属于ND动态边状态，NS无此变量，因此本项标为**不适用，未运行、无数值结果**。不把“ND已条件启用”误写为“ND未启用”，也不为S18另换候选或增加拟合任务。F清单处置17项计算、1项条件不适用。

## S17：每个外端口加固定L_b=10 μm扩散层阻力

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S17 compare

仅对外端口在原通道与access阻力之外再串联L_b/(D_b A_e)，L_b=1×10⁻⁵ m，内部边与容量不变。定向核对单图72个端口的阻力增量，内部边保持原值。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.158646047547807、最大盐量相对误差5.07731×10⁻¹³。

三起点88/71/69次均solver_success，但损失分别为4.158645982107509、3.6403468494633837、4.158645983875481。最低分支仅一个起点到达，另两起点并未落在其冻结1%范围内，状态**local_optimum_unstable**；成功不等于稳定，不追加冻结规则以外的起点或延长。最低分支变换参数(log10α,log10Q,β)=(1.1143828607113798×10⁻¹⁵,2.7634197482364913,4.5050316666543757×10⁻¹⁷)，log10α和β都贴下界。42点预测，最大盐量相对误差5.07731×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数0.859045–0.999356倍，最低分支重拟合后13.716422–8800.121110倍。200/350 nm的120 min响应比，原主结果1.637303、S17固定主参数1.637690、最低分支1.048424、原读数1.170732；不能挑这个更接近的比值作稳健支持。2600 nm的120 min最低分支预测23.505369 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S17_comparison.json，42点见WP5/v7_F_S17_point_predictions.csv。

其余输出为WP5/v7_F_S17_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有正式命令退出0；S17技术守恒通过但多起点判据**不通过**，E2整体仍不记PASS，论文未就绪。

## S16：膜内预载u=0.25并从供给池扣盐

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S16 compare

按02将所有膜节点初态u_i=0.25，接收池u_R=0，供给池u_L=1−0.25Σ_mem M_i/V_L，使初始总盐量仍等于V_L；若u_L<0即标无效，不作裁剪。正式运行前对六厚度、八种子、β=0/6、α=10⁵ m⁻¹的96个边界组合检查，最低u_L=0.97985，均非负且初始盐量守恒。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.158770836919058、最大盐量相对误差4.42×10⁻¹³。

三起点70/100/69次均solver_success，最佳重拟合L=4.158770838287824，另两起点在冻结1%内，multi_start_stable。该重拟合L较主参数前向高1.37×10⁻⁹，仅为数值量级，不能称改善。最佳变换参数(log10α,log10Q,β)=(4.999999956062299,1.5520281206958773,5.318423367584815)，log10α仍贴上界。42点预测，最大盐量相对误差7.42303×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数均为1，因固定池稳态不依赖初态；重拟合后0.999914–0.999942倍。200/350 nm的120 min响应比，原主结果1.637303、S16固定主参数1.637293、S16重拟合1.637367、原读数1.170732。2600 nm的120 min固定主参数预测17.780120、重拟合预测17.780146 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S16_comparison.json，42点见WP5/v7_F_S16_point_predictions.csv。

其余输出为WP5/v7_F_S16_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有正式命令退出0；S16自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现16/18，其余2项未运行。

## S15：供给池与接收池体积比r_V=2

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S15 compare

按02定义α=A/V_R，V_L=r_V V_R；仅令r_V从1变为2。定向检查确认供给池容量加倍，接收池及膜容量、孔网与全部边导通不变；主初态仍u_L=1、u_R=0，因此总初始盐量随V_L加倍，封闭两池的后续盐量守恒。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.227776605151786、最大盐量相对误差2.38×10⁻¹³。

三起点77/111/77次均solver_success，最佳重拟合L=4.164214581571358，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999855287,1.4250561806187303,5.165214540544507)，log10α仍贴上界。42点预测，最大盐量相对误差4.81107×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数均为1，因P是固定两池浓度边界下的几何/边导通结果；重拟合后1.249061–1.389397倍。200/350 nm的120 min响应比，原主结果1.637303、S15固定主参数1.921656、S15重拟合1.596537、原读数1.170732。2600 nm的120 min固定主参数预测23.569933、重拟合预测17.690135 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S15_comparison.json，42点见WP5/v7_F_S15_point_predictions.csv。

其余输出为WP5/v7_F_S15_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有正式命令退出0；S15自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现15/18，其余3项未运行。

## S14：固定D_b乘1.3

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S14 compare

仅将固定D_b乘预定1.3，D_Na、D_Cl和孔网不变；D不参与重拟合。定向检查确认容量、面积不变，所有内部边和端口导通系数均乘1.3。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.170120765233697、最大盐量相对误差6.01×10⁻¹³。

三起点82/133/56次均solver_success，最佳重拟合L=4.154260948074372，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999992879,1.5536255236387573,5.45683887876358)，log10α仍贴上界。42点预测，最大盐量相对误差7.90977×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数均为1.3倍（浮点偏差约5×10⁻¹⁴），重拟合后0.966137–1.063725倍。200/350 nm的120 min响应比，原主结果1.637303、S14固定主参数1.446418、S14重拟合1.673313、原读数1.170732。2600 nm的120 min重拟合预测17.854964 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S14_comparison.json，42点见WP5/v7_F_S14_point_predictions.csv。

其余输出为WP5/v7_F_S14_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有命令退出0；S14自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现14/18，其余4项未运行。

## S13：固定D_b乘0.7

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S13 compare

仅将D_b=2D_NaD_Cl/(D_Na+D_Cl)乘预定0.7，D_Na、D_Cl和孔网不变；D不参与重拟合。定向检查确认容量、面积不变，所有内部边和端口导通系数均乘0.7。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.2032670478192、最大盐量相对误差2.26×10⁻¹³。

三起点93/132/65次均solver_success，最佳重拟合L=4.166028818038497，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999999397,1.5495504097623345,5.12848934409787)，log10α仍贴上界。42点预测，最大盐量相对误差2.81838×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数均为0.7倍（浮点偏差约3×10⁻¹⁴），重拟合后0.922303–1.052461倍，反映α/Q/β的补偿而非D被重新确定。200/350 nm的120 min响应比，原主结果1.637303、S13固定主参数2.003929、S13重拟合1.586977、原读数1.170732。2600 nm的120 min重拟合预测17.662914 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S13_comparison.json，42点见WP5/v7_F_S13_point_predictions.csv。

其余输出为WP5/v7_F_S13_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有命令退出0；S13自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现13/18，其余5项未运行。

## S12：格距主值下层数改为round，0.5向上取整

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S12 compare

格距仍为主值s=d_f+d_p；仅将层数由max(1,ceil(h/s))改为max(1,floor(h/s+0.5))，半数向上。定向检查含精确0.5边界，六厚度的(s_nm,n_z)为65:(78.97962,1)、120:(78.97962,2)、200:(78.97962,3)、350:(67.12452,5)、700:(75.11260,9)、2600:(95.51265,27)。48张S图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.171602915374056、最大盐量相对误差5.98628×10⁻¹³。

三起点78/112/60次均solver_success，最佳重拟合L=4.170581878880215，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999884281,1.5513077043228425,5.350458267328336)，log10α仍贴上界。42点预测，最大盐量相对误差5.98628×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数1.000000–1.088404倍，重拟合后0.945487–1.016153倍。200/350 nm的120 min响应比，原主结果1.637303、S12固定主参数1.508437、S12重拟合1.563194、原读数1.170732。2600 nm的120 min重拟合预测17.743959 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S12_comparison.json，42点见WP5/v7_F_S12_point_predictions.csv。

其余输出为WP5/v7_F_S12_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有命令退出0；S12自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现12/18，其余6项未运行。

## S11：格距取0.75(d_f+d_p)

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S11 compare

正式运行前的定向几何检查在65 nm首次失败：原截断对数正态求根固定区间不能夹住所需均值，ValueError: f(a) and f(b) must have different signs；120/200/700/2600 nm同因失败，350 nm成功。仅将数值求根区间改为按符号扩展，未改分布均值、σ、截断规则或科学阈值。修复后六厚度均可生成；旧S10单图的孔径定向回归一致。S11按s=0.75(d_f+d_p)和n_z=max(1,ceil(h/s))重算格距、层数、面积及体积预算。六厚度的(s_nm,n_z)为65:(59.2347,2)、120:(59.2347,3)、200:(59.2347,4)、350:(50.3434,7)、700:(56.3344,13)、2600:(71.6345,37)。

48张S图最大体积相对误差2.66×10⁻¹⁶。主参数前向L=4.2627182646121895、最大盐量相对误差4.12430×10⁻¹³。三起点68/92/44次均solver_success，最佳重拟合L=4.068417892476238，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999999999,1.5845309973292483,5.8925222603259275)，log10α仍贴上界。42点预测，最大盐量相对误差4.12430×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数0.178846–2.953554倍，重拟合后0.061000–1.223608倍。200/350 nm的120 min响应比，原主结果1.637303、S11固定主参数1.129712、S11重拟合2.210009、原读数1.170732；全B重拟合使这项局部比值偏离读数更远。2600 nm的120 min重拟合预测19.194220 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S11_comparison.json，42点见WP5/v7_F_S11_point_predictions.csv。

其余输出为WP5/v7_F_S11_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。正式矩阵命令均退出0；S11自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现11/18，其余7项未运行。

## S10：三个图读锚点孔径与面积率联合取高端

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S10 compare

冻结锚点81/182/357的(dp_nm,φ_A)分别由(66,0.310)/(52,0.286)/(43,0.292)改为其图读区间高端(67,0.315)/(53,0.291)/(44,0.297)；604文字锚点仍为(52.4,0.345)，df不变。48张S图最大体积相对误差1.49×10⁻¹⁶。主参数前向L=4.1640289518187314、最大盐量相对误差5.36779×10⁻¹³。

三起点78/124/52次均solver_success，最佳重拟合L=4.15070724053149，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999755226,1.5563091501456667,5.449178463255001)，log10α仍贴上界。42点预测，最大盐量相对误差5.36779×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数1.000000–1.272261倍，重拟合后0.795795–0.969900倍。200/350 nm的120 min响应比，原主结果1.637303、S10固定主参数1.452600、S10重拟合1.664448、原读数1.170732。2600 nm的120 min重拟合预测17.946650 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S10_comparison.json，42点见WP5/v7_F_S10_point_predictions.csv。

其余输出为WP5/v7_F_S10_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有命令退出0；S10自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现10/18，其余8项未运行。

## S09：三个图读锚点孔径与面积率联合取低端

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S09 compare

冻结锚点81/182/357的(dp_nm,φ_A)分别由(66,0.310)/(52,0.286)/(43,0.292)改为其图读区间低端(65,0.305)/(51,0.281)/(42,0.287)；604文字锚点仍为(52.4,0.345)，df不变。48张S图最大体积相对误差1.54×10⁻¹⁶。主参数前向L=4.179185356797148、最大盐量相对误差5.03080×10⁻¹³。

三起点86/121/69次均solver_success，最佳重拟合L=4.1640262648570125，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999999969,1.5494170605834863,5.195271837598393)，log10α仍贴上界。42点预测，最大盐量相对误差5.03080×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数0.778393–1.000000倍，重拟合后1.016714–1.240363倍。200/350 nm的120 min响应比，原主结果1.637303、S09固定主参数1.888704、S09重拟合1.624002、原读数1.170732。2600 nm的120 min重拟合预测17.681251 μS/cm，对原读数1363.25，主要失配未解。逐膜比较见WP5/v7_F_S09_comparison.json，42点见WP5/v7_F_S09_point_predictions.csv。

其余输出为WP5/v7_F_S09_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。所有命令退出0；S09自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现9/18，其余9项未运行。

## S08：锚点范围外改为ln h线性外延

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S08 compare

仅替换S系列的锚点范围外外推；锚点内数值保持原规则，65/120/200 nm用最薄两锚点，2600 nm用最厚两锚点，均按ln h线性外延。四个外延厚度的(df,dp,φ_A)依次为65:(37.276025,77.385513,0.485876)、120:(32.086458,65.263312,0.417527)、200:(27.762622,55.163350,0.360581)、2600:(31.400628,71.316782,0.319114)，长度单位nm。其余配置不变。48张S图最大体积相对误差1.78×10⁻¹⁶。主参数前向L=4.180387907641261、最大盐量相对误差3.62×10⁻¹²。

三起点61/83/48次均solver_success，最佳重拟合L=4.180018785594484，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999999905,1.549518542263316,5.335936817038651)，log10α贴上界。42点预测，最大盐量相对误差3.62446×10⁻¹²、最小状态0。

六膜稳态P相对原主值：固定主参数1.000000–12.313629倍，重拟合后0.963046–12.094135倍。200/350 nm的120 min响应比，原主结果1.637303、S08固定主参数1.637340、S08重拟合1.670504、原读数1.170732。2600 nm的120 min重拟合预测17.700694 μS/cm，对原读数1363.25；稳态P可大幅改变而此时刻厚膜失配不变。逐膜比较见WP5/v7_F_S08_comparison.json，42点见WP5/v7_F_S08_point_predictions.csv。

其余输出为WP5/v7_F_S08_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。S08自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现8/18，其余10项未运行。

## S07：孔体积分数系数κ_V从主值1改为0.25

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 forward
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S07 compare

首次prepare已写出48张S图和审计表，随后Windows沙箱进程池报PermissionError（WinError 5），命令退出1，未生成主参数前向结果。为保留已完成图谱并记录失败，添加只读取既有图谱的forward续接入口；获执行权限后forward退出0。其余命令均退出0。没有更换种子、阈值或矩阵。

48张图最大体积相对误差1.82×10⁻¹⁶。主参数前向L=5.09615365040953，最大盐量相对误差4.41×10⁻¹⁴。三起点分别调用56/74/44次，均solver_success；最佳重拟合L=4.17403597156931，其余两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999998031,1.5476734905938678,4.594008920438267)，log10α贴上界。42点预测，重拟合最大盐量相对误差2.88524×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数0.089088–0.330243倍，重拟合后0.362785–1.299341倍。200/350 nm的120 min响应比，原主结果1.637303、S07固定主参数6.449580、S07重拟合1.538109、原读数1.170732，重拟合仍未达到原读数。2600 nm的120 min重拟合预测17.630709 μS/cm，对原读数1363.25，主要失配持续。逐膜P和120 min比较见WP5/v7_F_S07_comparison.json，42点见WP5/v7_F_S07_point_predictions.csv。

其余输出为WP5/v7_F_S07_networks.jsonl.gz、geometry_audit.csv、main_forward.json、fit_runs.jsonl及report.json。S07自身达到冻结多起点条件，E2整体仍不记PASS，论文未就绪。F现7/18，其余11项未运行。

## S06：孔体log宽度σ_ln从主值0.35改为0.70

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S06 compare

仅将几何配置σ_ln改为0.70，同一六h与八种子重建48张S图；最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.200768176236509；三起点实际调用88/155/78次，均solver_success，最佳重拟合L=4.186460102939421，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999992982212,1.5452991102540021,5.446540381300229)，log10α仍贴上界1%邻域。前向及重拟合最大盐量相对误差8.24755×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数为1.173290–1.528657倍，重拟合后为0.972574–1.219526倍。200/350 nm的120 min响应比，原主结果1.637303、S06固定主参数1.391993、S06重拟合1.557971、原读数1.170732，排序未反转。2600 nm的120 min重拟合预测17.499482 μS/cm，对原读数1363.25，严重失配未解除。逐膜P和120 min比较见WP5/v7_F_S06_comparison.json，42点及所有时间见WP5/v7_F_S06_point_predictions.csv。

输出路径另含WP5/v7_F_S06_networks.jsonl.gz、WP5/v7_F_S06_geometry_audit.csv、WP5/v7_F_S06_main_forward.json、WP5/v7_F_S06_fit_runs.jsonl及WP5/v7_F_S06_report.json。S06自身达到冻结多起点条件，不改变E2整体不记PASS或论文未就绪。F现6/18，其余12项未运行。

## S05：孔体log宽度σ_ln从主值0.35改为0

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S05 compare

仅将几何配置σ_ln改为0，同一六h与八种子重建48张S图；最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.211549554329296；三起点实际调用60/131/40次，均solver_success，最佳重拟合L=4.085298077710736，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999999939,1.576204228090184,5.749909967482574)，log10α仍贴上界1%邻域。前向及重拟合最大盐量相对误差5.82794×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数为0.161634–2.219181倍，重拟合后为0.070939–1.101083倍。200/350 nm的120 min响应比，原主结果1.637303、S05固定主参数1.197742、S05重拟合2.024678、原读数1.170732，排序未反转且重拟合后的比值偏离读数更远。2600 nm的120 min重拟合预测18.822598 μS/cm，对原读数1363.25，严重失配未解除。逐膜P和120 min比较见WP5/v7_F_S05_comparison.json，42点及所有时间见WP5/v7_F_S05_point_predictions.csv。

输出路径另含WP5/v7_F_S05_networks.jsonl.gz、WP5/v7_F_S05_geometry_audit.csv、WP5/v7_F_S05_main_forward.json、WP5/v7_F_S05_fit_runs.jsonl及WP5/v7_F_S05_report.json。S05自身达到冻结多起点条件，不改变E2整体不记PASS或论文未就绪。F现5/18，其余13项未运行。

## S04：轴向边概率p_z从主值0.65改为1

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S04 compare

仅将几何配置p_z改为1，同一六h与八种子重建48张S图；最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.308731020148614；三起点实际调用69/160/85次，均solver_success，最佳重拟合L=4.1383748221612615，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999993953,1.5604669597682927,5.862381823358179)，log10α仍贴上界1%邻域。前向及重拟合最大盐量相对误差5.77399×10⁻¹³、最小状态0。

六膜稳态P相对原主值：固定主参数为1.000000–2.660634倍，重拟合后为0.455741–1.000619倍。65 nm只有一层，没有轴向边，固定主参数P比恰为1。200/350 nm的120 min响应比，原主结果1.637303、S04固定主参数1.096974、S04重拟合1.716830、原读数1.170732，排序未反转。2600 nm的120 min重拟合预测18.135574 μS/cm，对原读数1363.25，严重失配未解除。逐膜P和120 min比较见WP5/v7_F_S04_comparison.json，42点及所有时间见WP5/v7_F_S04_point_predictions.csv。

输出路径另含WP5/v7_F_S04_networks.jsonl.gz、WP5/v7_F_S04_geometry_audit.csv、WP5/v7_F_S04_main_forward.json、WP5/v7_F_S04_fit_runs.jsonl及WP5/v7_F_S04_report.json。S04自身达到冻结多起点条件，不改变E2整体不记PASS或论文未就绪。F现4/18，其余14项未运行。

## S03：轴向边概率p_z从主值0.65改为0.35

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S03 compare

仅将几何配置p_z改为0.35，同一六h与八种子重建48张S图；最大体积相对误差1.82×10⁻¹⁶。主参数前向L=4.640249788857839；三起点实际调用80/115/69次，均solver_success，最佳重拟合L=4.15735227411142，另两起点在冻结1%内，multi_start_stable。最佳变换参数(log10α,log10Q,β)=(4.999999999954938,1.553204525491384,4.72998116885453)，log10α仍贴上界1%邻域。前向及重拟合最大盐量相对误差1.22106×10⁻¹²、最小状态0。

六膜稳态P相对原主值：固定主参数为0.259910–1.000000倍，重拟合后为0.905546–2.356785倍。65 nm只有一层，没有轴向边，固定主参数P比恰为1；其余厚度受到p_z改动。200/350 nm的120 min响应比，原主结果1.637303、S03固定主参数4.054951、S03重拟合1.607562、原读数1.170732，排序未反转。2600 nm的120 min重拟合预测17.799242 μS/cm，对原读数1363.25，严重失配未解除。逐膜P和120 min比较见WP5/v7_F_S03_comparison.json，42点及所有时间见WP5/v7_F_S03_point_predictions.csv。

输出路径另含WP5/v7_F_S03_networks.jsonl.gz、WP5/v7_F_S03_geometry_audit.csv、WP5/v7_F_S03_main_forward.json、WP5/v7_F_S03_fit_runs.jsonl及WP5/v7_F_S03_report.json。S03自身达到冻结多起点条件，不改变E2整体不记PASS或论文未就绪。F现3/18，其余15项未运行。

## S02：喉比q从主值0.5改为1

复用S01代码入口并用冻结case标识隔离输出。实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 extend 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 extend 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py --case S02 compare

48张S孔网最大体积相对误差1.82×10⁻¹⁶，主参数前向L=4.611770655025605。原三起点均用尽200次实际调用而记max_eval，不能记成功；此时report因没有成功起点拒绝生成，失败已保留。按04唯一加时规则，选择未收敛中L最低的起点2（4.297967114115245）和起点0（4.297979289553569），分别加22/18次，累计222/218次后均solver_success。起点1仍为max_eval。最佳重拟合L=4.2979582953312，另一个成功起点的L在冻结1%内，最终multi_start_stable；总实际调用640次。最佳变换参数(log10α,log10Q,β)=(2.9802490743233365,1.5025282193847382,5.999999999644984)，β贴上界1%邻域。前向及重拟合最大盐量相对误差3.07867×10⁻¹¹、最小状态0。

六膜稳态P相对原主值：固定主参数为63.0690–574.0950倍，重拟合后仍为32.3721–236.8768倍。200/350 nm的120 min响应比，原主结果1.637303、S02固定主参数0.999990857（几乎相等，名义上翻转）、S02重拟合1.100981、原读数1.170732；近等号不能当成稳健反向排序。2600 nm的120 min S02重拟合预测15.769911 μS/cm，对原读数1363.25，既有严重失配未解除。逐膜P和120 min比较见WP5/v7_F_S02_comparison.json，42点及所有时间见WP5/v7_F_S02_point_predictions.csv。

输出路径还包括WP5/v7_F_S02_networks.jsonl.gz、WP5/v7_F_S02_geometry_audit.csv、WP5/v7_F_S02_main_forward.json、WP5/v7_F_S02_fit_runs.jsonl、WP5/v7_F_S02_fit_extensions.jsonl和WP5/v7_F_S02_report.json。E2整体技术状态仍不记PASS；S02自身达到冻结多起点条件。F现2/18，其余16项未运行。

## S01：喉比q从主值0.5改为0.25

实际运行命令：

    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py prepare
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py fit 0
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py fit 1
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py fit 2
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py report
    .\.venv\Scripts\python.exe src\effective_network_model\v7_f_s01.py compare

重建S几何六h×八预设种子共48张图；最大体积相对误差1.82×10⁻¹⁶。主参数前向L=14.746003009291458；三起点实际目标调用111/82/98次，均solver_success，最佳重拟合L=4.201007759347661，另外两个在冻结1%相对L内，multi_start_stable。原主结果L=4.1587531270409706。最佳变换参数(log10α,log10Q,β)=(4.999999994249411,1.5412667527120534,2.285054899203925)，log10α仍在上界1%邻域。前向及重拟合最大盐量相对误差3.08697×10⁻¹³、最小状态0。

固定主参数时，六膜稳态P相对原主结果的比值范围3.39020×10⁻⁶–1.39144×10⁻⁴；重拟合补偿后范围0.745065–1.253681。200/350 nm的120 min响应比，原主结果1.637303、S01固定主参数450.388259、S01重拟合1.438853、原读数1.170732；特定排序未反转，但响应幅度与厚度趋势对q和参数补偿敏感。2600 nm的120 min预测，原主结果17.780108 μS/cm、S01固定主参数0.00361249、S01重拟合17.159187，原读数1363.25，原有严重失配未解决。六h的P、120 min响应及42点预测分别在WP5/v7_F_S01_comparison.json、WP5/v7_F_S01_report.json和WP5/v7_F_S01_point_predictions.csv。

数据身份：B原读数来自v7任务包inputs/salt，单位μS/cm；P为模型稳态有效通量/面积，单位m/s，并非实验实测P。结果仅说明S01这个单因素情景，F其余17项未运行。S01代码入口为src/effective_network_model/v7_f_s01.py；图表和日志还见WP5/v7_F_S01_networks.jsonl.gz、WP5/v7_F_S01_geometry_audit.csv、WP5/v7_F_S01_main_forward.json、WP5/v7_F_S01_fit_runs.jsonl。

技术核对：脚本py_compile通过；48图、3起点、42点、守恒、L重算、1%稳定性均定向核对通过。一次py_compile命令因PowerShell路径多写反斜杠而退出1，随后用正确路径重跑退出0；未改变计算。拟合日志中“best L”打印字段曾使用当前L，已修正为记录的最佳L；日志文字问题不影响保存的优化值，不重跑已通过的数值矩阵。
