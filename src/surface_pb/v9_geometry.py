"""Exact square-fiber layer geometry for the v9 equilibrium and observation model."""

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"


@dataclass(frozen=True)
class Layer:
    layer_id: int
    x_start_periodic_nm: float
    z_bottom_nm: float
    z_top_nm: float
    d_nm: float
    wet_perimeter_nm: float
    is_outer: bool
    is_core: bool


@dataclass(frozen=True)
class Geometry:
    kind: str
    h_nm: float
    d_nm: float
    g_nm: float
    N: int
    period_nm: float
    vertical_spacing_nm: float | None
    aligned: bool
    wet_area_m2_per_m: float
    projected_area_m2_per_m: float
    layers: tuple[Layer, ...]


def make_geometry(kind, d_nm, g_nm, N, h_nm, aligned=False):
    if kind not in ("free", "backed") or d_nm <= 0 or g_nm <= 0 or N < 1:
        raise ValueError("kind, width, gap or layer count outside geometry contract")
    period = d_nm + g_nm
    if N == 1:
        if not math.isclose(h_nm, d_nm, rel_tol=0, abs_tol=1e-8):
            raise ValueError("a single square-fiber layer has h=d")
        spacing = None
    else:
        spacing = (h_nm - d_nm) / (N - 1)
        if spacing <= d_nm:
            raise ValueError("square fiber layers overlap vertically")
    z_origin = -h_nm / 2 if kind == "free" else -h_nm
    layers = []
    for i in range(N):
        x_start = (g_nm / 2 + (0 if aligned else (i % 2) * period / 2)) % period
        z_bottom = z_origin + (spacing or 0) * i
        is_outer = i == 0 or i == N - 1
        is_core = N > 1 and (i == (N - 1) // 2 or i == N // 2)
        # At the neutral backing only the bottom face of the first solid
        # fiber touches support. All other faces contact the modelled water.
        perimeter = d_nm * (3 if kind == "backed" and i == 0 else 4)
        layers.append(Layer(i + 1, x_start, z_bottom, z_bottom + d_nm,
                            d_nm, perimeter, is_outer, is_core))
    assert math.isclose(layers[-1].z_top_nm - z_origin, h_nm, abs_tol=1e-8)
    wet_area = sum(layer.wet_perimeter_nm for layer in layers) * 1e-9
    expected = (4 * N - (kind == "backed")) * d_nm * 1e-9
    assert math.isclose(wet_area, expected, rel_tol=1e-14)
    return Geometry(kind, h_nm, d_nm, g_nm, N, period, spacing, aligned,
                    wet_area, period * 1e-9, tuple(layers))


def family_geometry(spec, draw, kind, aligned=False):
    given = spec["geometry"]["family"][str(draw)]
    d_nm, g_nm, h_nm = given["d_nm"], given["g_nm"], given["h_nm"]
    N = 1 + math.floor((h_nm - d_nm) / (d_nm + g_nm) + 0.5)
    if N != given["N"]:
        raise AssertionError(f"family layer rule mismatch for {draw}")
    return make_geometry(kind, d_nm, g_nm, N, h_nm, aligned)


def pure_geometry(spec, N, kind, aligned=False):
    given = spec["geometry"]["pure_h"]
    if N not in (1, 3, 4, 6, 11, 20, 40):
        raise ValueError("N outside frozen pure-thickness or thick-end grid")
    d_nm, g_nm = given["d_nm"], given["g_nm"]
    h_nm = d_nm + (N - 1) * (d_nm + g_nm)
    return make_geometry(kind, d_nm, g_nm, N, h_nm, aligned)


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    geometries = [(f"family_{draw}_{kind}", family_geometry(spec, draw, kind))
                  for draw in (182, 357, 604) for kind in ("free", "backed")]
    geometries += [(f"pure_N{N}_{kind}", pure_geometry(spec, N, kind))
                   for N in (1, 3, 4, 6, 11, 20, 40) for kind in ("free", "backed")]
    records = []
    for case_id, geometry in geometries:
        assert len(geometry.layers) == geometry.N
        for layer in geometry.layers:
            row = {"case_id": case_id, "kind": geometry.kind,
                   "h_nm": geometry.h_nm, "N": geometry.N,
                   "period_nm": geometry.period_nm,
                   "vertical_spacing_nm": geometry.vertical_spacing_nm or "NA",
                   "wet_area_m2_per_m": geometry.wet_area_m2_per_m,
                   "projected_area_m2_per_m": geometry.projected_area_m2_per_m}
            row.update(asdict(layer))
            records.append(row)
    OUT = ROOT / "WP3"
    OUT.mkdir(exist_ok=True)
    with (OUT / "v9_几何层表.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    summary = {"status": "PASS_GEOMETRY_ONLY", "PB_solver_run": False,
               "geometry_cases": len(geometries), "layer_rows": len(records),
               "family_N": {str(draw): family_geometry(spec, draw, "free").N
                            for draw in (182, 357, 604)},
               "free_N1_outer_count": sum(layer.is_outer for layer in pure_geometry(spec, 1, "free").layers),
               "free_N1_core_count": sum(layer.is_core for layer in pure_geometry(spec, 1, "free").layers),
               "supported_bottom_wet_face": False}
    (OUT / "v9_几何接口核对.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
