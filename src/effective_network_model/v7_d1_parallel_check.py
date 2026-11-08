"""Check process-parallel D1 forward evaluation against the serial reference."""

from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves, graphs, init_forward_worker, predict, predict_parallel
from v7_network_transport import SPEC_PATH


def main():
    spec=json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    data,nets=curves(),graphs()
    saved=json.loads((OUT/"v7_d1_fit_fullB_runs.jsonl").read_text(encoding="utf-8").splitlines()[0])
    x=saved["best_transformed"]
    alpha,q=10**x[0],10**x[1]
    start=time.perf_counter()
    serial,_=predict(data,nets,spec,alpha,0.0,q)
    serial_s=time.perf_counter()-start
    start=time.perf_counter()
    with ProcessPoolExecutor(max_workers=6,initializer=init_forward_worker) as pool:
        parallel,_=predict_parallel(pool,alpha,0.0,q)
        first_s=time.perf_counter()-start
        warm_start=time.perf_counter()
        parallel_warm,_=predict_parallel(pool,alpha,0.0,q)
        warm_s=time.perf_counter()-warm_start
    parallel_s=time.perf_counter()-start
    max_diff=max(float(np.max(abs(serial[h]-parallel[h]))) for h in H_VALUES)
    max_diff=max(max_diff,max(float(np.max(abs(serial[h]-parallel_warm[h]))) for h in H_VALUES))
    evidence={"alpha_m_inv":alpha,"Q_original_unit":q,"beta":0.0,
              "workers":6,"serial_seconds":serial_s,"parallel_first_seconds":first_s,
              "parallel_warm_seconds":warm_s,"parallel_total_with_shutdown_seconds":parallel_s,
              "max_absolute_prediction_difference":max_diff,
              "status":"PASS" if max_diff<=1e-12 else "FAIL"}
    (OUT/"v7_d1_parallel_check.json").write_text(json.dumps(evidence,indent=2)+"\n",encoding="utf-8")
    print(evidence)
    assert evidence["status"]=="PASS"


if __name__=="__main__":
    main()
