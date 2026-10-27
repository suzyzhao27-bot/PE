# v7 冻结矩阵批量入口

这是**从干净输入工作副本**重跑已有矩阵的命令清单，不是本次代表复现的实测日志。正式包已含历史输出；不要在正式包上直接执行拟合命令，部分入口会拒绝重复起点，部分会更新归档。先按`README_reproduce.md`完成约25秒的代表前向回放。完整矩阵原执行证据在WP2/WP4/WP5，逐项失败和合法加时在`WP5/v7_H_冻结矩阵闭合与主张审阅.md`。下述命令只列本版既有任务，不增加搜索。

以下以 PowerShell、项目根目录和已有 `python` 环境为例，`$py`指现成解释器，`$m='src/effective_network_model'`，`$p='src/finite_pnp'`。每个起点、任务和报告顺序执行；命令失败时保留输出并按04的限定加时规则处理，不能跳过、换种子或改阈值。

```powershell
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true  # PowerShell 7
$py = '.\.venv\Scripts\python.exe'
$m = 'src/effective_network_model'
$p = 'src/finite_pnp'
& $py "$m/v7_network_geometry.py"
& $py "$m/v7_straight_chain.py"
& $py "$m/v7_network_transport.py"
& $py "$m/v7_shared_chain_validation.py"
& $py "$m/v7_thermo_capacity.py"
foreach ($model in @('N0','NS')) { & $py "$m/v7_t6_formal.py" $model }
& $py "$m/v7_c1_main.py"
& $py "$m/v7_c1_summary.py"
& $py "$m/v7_c1_ablations.py"
& $py "$m/v7_c2_domain.py"
& $py "$m/v7_c2_refinement.py"
& $py "$p/v7_c3_main.py"
& $py "$p/v7_c3_compare.py"
& $py "$p/v7_c3_refinement.py"
```

N0/NS的完整B、时间、六次整膜及L早期条件校准：

```powershell
foreach ($model in @('N0','NS')) {
  foreach ($start in 0..11) { & $py "$m/v7_d1_fit_start.py" $model $start }
  & $py "$m/v7_d1_report_fullB.py" $model
  foreach ($start in 0..3) { & $py "$m/v7_d1_time_fit_start.py" $model $start }
  & $py "$m/v7_d1_report_time.py" $model
  foreach ($h in @(65,120,200,350,700,2600)) {
    foreach ($start in 0..3) { & $py "$m/v7_d1_film_fit_start.py" $model $h $start }
  }
  if ($model -eq 'NS') {
    & $py "$m/v7_d1_film_fit_extend.py" NS 350 2
    & $py "$m/v7_d1_film_fit_extend.py" NS 350 0
  }
  & $py "$m/v7_d1_report_film.py" $model
  foreach ($start in 0..3) { & $py "$m/v7_d1_l_fit_start.py" $model $start }
  & $py "$m/v7_d1_report_l.py" $model
}
& $py "$m/v7_d2_decision.py"
```

NS的350 nm整膜原四起点均`max_eval`；既有记录确定起点2/0为损失最低两项，各作一次限定加时后再执行`report_film`。原起点记录保留。ND仅在D2两条件满足时启用，本轮已有`WP4/v7_D2_decision.json`记录为真：

```powershell
& $py "$m/v7_nd_t6_formal.py"
& $py "$m/v7_nd_t8_formal.py"
foreach ($start in 0..11) { & $py "$m/v7_nd_fit_start.py" $start }
& $py "$m/v7_nd_report_fullb.py"
foreach ($start in 0..3) { & $py "$m/v7_nd_time_fit_start.py" $start }
& $py "$m/v7_nd_report_time.py"
foreach ($h in @(65,120,200,350,700,2600)) {
  foreach ($start in 0..3) { & $py "$m/v7_nd_film_fit_start.py" $h $start }
}
& $py "$m/v7_nd_report_film.py"
foreach ($start in 0..3) { & $py "$m/v7_nd_l_fit_start.py" $start }
& $py "$m/v7_nd_report_l.py"
& $py "$m/v7_compile_wp4.py"
```

第5步的NS单一诊断候选：

```powershell
foreach ($parameter in 0..2) {
  foreach ($position in 0..8) {
    foreach ($start in 0..2) { & $py "$m/v7_e1_profile_start.py" $parameter $position $start }
    if ($parameter -eq 1 -and $position -eq 0) {
      & $py "$m/v7_e1_profile_extend.py" 1 0 0
    }
    & $py "$m/v7_e1_report_point.py" $parameter $position
  }
}
& $py "$m/v7_e1_compile_alpha.py"
& $py "$m/v7_e1_compile_remaining.py" 1
& $py "$m/v7_e1_compile_remaining.py" 2
& $py "$m/v7_e1_equivalent_b_curves.py"
foreach ($rank in 1..3) {
  foreach ($start in 0..3) { & $py "$m/v7_e1_equivalent_l_start.py" $rank $start }
}
& $py "$m/v7_e1_equivalent_l_report.py"
& $py "$m/v7_e1_equivalent_figures.py"
& $py "$m/v7_e2_prepare.py"
foreach ($truth in @('fullB_best','equivalent_1')) {
  foreach ($noise in @(0.0,0.01)) {
    foreach ($sampling in @('sparse_actual','dense_1min')) {
      foreach ($set in 0..1) {
        foreach ($start in 0..3) { & $py "$m/v7_e2_inverse_start.py" $truth $noise $sampling $set $start }
        & $py "$m/v7_e2_report_task.py" $truth $noise $sampling $set
      }
    }
  }
}
& $py "$m/v7_e3_prepare.py"
foreach ($eta in @(0.01,0.05)) {
  foreach ($seed in 8101..8106) {
    foreach ($start in 0..3) { & $py "$m/v7_e3_fit_start.py" $eta $seed $start }
    & $py "$m/v7_e3_report_task.py" $eta $seed
  }
}
foreach ($case in 1..17) {
  $id = 'S{0:D2}' -f $case
  & $py "$m/v7_f_s01.py" --case $id prepare
  foreach ($start in 0..2) { & $py "$m/v7_f_s01.py" --case $id fit $start }
  if ($id -eq 'S02') {
    & $py "$m/v7_f_s01.py" --case S02 extend 2
    & $py "$m/v7_f_s01.py" --case S02 extend 0
  }
  & $py "$m/v7_f_s01.py" --case $id report
  & $py "$m/v7_f_s01.py" --case $id compare
}
& $py "$m/v7_g_regime.py" compute
& $py "$m/v7_g_regime.py" plot
```

E1的Q首点、F/S02按既有原起点状态和04规则只加时一次；若干净重跑出现不同状态，应停下核查，不盲目继续固定索引。F/S18只属于ND，不进入NS诊断候选。作图入口为`v7_plot_full_b.py`、`v7_plot_time_holdout.py`、`v7_plot_film_holdout.py`、`v7_plot_l_transfer.py`和`v7_paper_figures.py`，只读已存结果。以上批量清单未在正式包逐条重跑；已实际重算的四项代表结果见`tests/representative_reproduction.json`。因此批量入口的存在不把全矩阵标为第二次独立复现。
