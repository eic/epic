#!/usr/bin/env python3

import argparse
import ast
import csv
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple


UNIT_MM = {
    "mm": 1.0,
    "cm": 10.0,
    "m": 1000.0,
    "um": 0.001,
}

SAFE_FUNCS = {
    "min": min,
    "max": max,
    "atan": math.atan,
    "exp": math.exp,
    "sqrt": math.sqrt,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
}


class DiskSpec:
    def __init__(self, detector: str, disk: str, rmin_mm: float, rmax_mm: float,
                 layer_length_mm: float, module_thickness_mm: float) -> None:
        self.detector = detector
        self.disk = disk
        self.rmin_mm = rmin_mm
        self.rmax_mm = rmax_mm
        self.layer_length_mm = layer_length_mm
        self.module_thickness_mm = module_thickness_mm


class ExpressionResolver:
    def __init__(self) -> None:
        self.raw = {}  # type: Dict[str, str]
        self.cache = {}  # type: Dict[str, float]

    def add_constants_from_xml(self, xml_path: Path) -> None:
        root = ET.parse(xml_path).getroot()
        for const in root.findall(".//constant"):
            name = const.get("name")
            value = const.get("value")
            if name and value:
                self.raw[name] = value

    def resolve(self, name: str) -> float:
        if name in self.cache:
            return self.cache[name]
        if name in UNIT_MM:
            return UNIT_MM[name]
        if name == "pi":
            return math.pi
        if name not in self.raw:
            raise KeyError(f"unknown constant '{name}'")
        value = self.eval_expr(self.raw[name])
        self.cache[name] = value
        return value

    def eval_expr(self, expr: str) -> float:
        parsed = ast.parse(expr, mode="eval")
        return float(self._eval_node(parsed.body))

    def _eval_node(self, node: ast.AST) -> float:
        if isinstance(node, ast.Num):
            return float(node.n)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            if node.id in SAFE_FUNCS:
                return SAFE_FUNCS[node.id]
            return self.resolve(node.id)
        if isinstance(node, ast.BinOp):
            left = self._eval_node(node.left)
            right = self._eval_node(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                return left**right
            raise ValueError(f"unsupported operator in expression: {ast.dump(node)}")
        if isinstance(node, ast.UnaryOp):
            value = self._eval_node(node.operand)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                return -value
            raise ValueError(f"unsupported unary operator: {ast.dump(node)}")
        if isinstance(node, ast.Call):
            func = self._eval_node(node.func)
            args = [self._eval_node(arg) for arg in node.args]
            return float(func(*args))
        raise ValueError(f"unsupported expression node: {ast.dump(node)}")


def default_definitions_xml(geometry_xml: Path) -> Path:
    tracking_dir = geometry_xml.parent
    candidates = [
        tracking_dir / "definitions_craterlake.xml",
        tracking_dir / "definitions_craterlake_tile.xml",
        tracking_dir / "definitions.xml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("unable to infer definitions XML; pass --definitions-xml")


def parse_module_thickness(detector: ET.Element, resolver: ExpressionResolver) -> float:
    thickness = 0.0
    module = detector.find("module")
    if module is None:
        raise ValueError(f"detector '{detector.get('name')}' has no <module>")
    for component in module.findall("module_component"):
        value = component.get("thickness")
        if not value:
            continue
        thickness += resolver.eval_expr(value)
    return thickness


def parse_disks(geometry_xml: Path, definitions_xml: Path, disk_filter: List[str]) -> List[DiskSpec]:
    resolver = ExpressionResolver()
    resolver.add_constants_from_xml(definitions_xml)
    resolver.add_constants_from_xml(geometry_xml)

    root = ET.parse(geometry_xml).getroot()
    disks = []  # type: List[DiskSpec]
    for detector in root.findall(".//detector"):
        if detector.get("type") != "epic_TileEndcapTracker":
            continue
        detector_name = detector.get("name", "")
        module_thickness = parse_module_thickness(detector, resolver)
        for layer in detector.findall("layer"):
            disk_name = layer.get("name") or f"{detector_name}_disk{layer.get('id', '')}"
            if disk_filter and disk_name not in disk_filter:
                continue
            envelope = layer.find("envelope")
            if envelope is None:
                continue
            disks.append(
                DiskSpec(
                    detector=detector_name,
                    disk=disk_name,
                    rmin_mm=resolver.eval_expr(envelope.get("rmin", "0")),
                    rmax_mm=resolver.eval_expr(envelope.get("rmax", "0")),
                    layer_length_mm=resolver.eval_expr(envelope.get("length", "0")),
                    module_thickness_mm=module_thickness,
                )
            )
    return disks


def corners_inside_annulus(x_min: float, y_min: float, x_size: float, y_size: float,
                           rmin_mm: float, rmax_mm: float) -> bool:
    corners = (
        (x_min, y_min),
        (x_min + x_size, y_min),
        (x_min, y_min + y_size),
        (x_min + x_size, y_min + y_size),
    )
    for x_val, y_val in corners:
        radius = math.hypot(x_val, y_val)
        if radius < rmin_mm or radius > rmax_mm:
            return False
    return True


def z_inside_layer(dz_mm: float, module_thickness_mm: float, layer_length_mm: float) -> bool:
    half_module = module_thickness_mm / 2.0
    half_layer = layer_length_mm / 2.0
    return dz_mm - half_module >= -half_layer and dz_mm + half_module <= half_layer


def generate_tiles_for_disk(disk: DiskSpec, tile_x_mm: float, tile_y_mm: float,
                            x0_mm: Optional[float], y0_mm: Optional[float]) -> List[Tuple[float, float]]:
    seed_x = x0_mm if x0_mm is not None else -tile_x_mm / 2.0
    seed_y = y0_mm if y0_mm is not None else disk.rmin_mm

    min_ix = math.floor((-disk.rmax_mm - seed_x) / tile_x_mm) - 1
    max_ix = math.ceil((disk.rmax_mm - seed_x) / tile_x_mm) + 1
    min_iy = math.floor((-disk.rmax_mm - seed_y) / tile_y_mm) - 1
    max_iy = math.ceil((disk.rmax_mm - seed_y) / tile_y_mm) + 1

    placements = []  # type: List[Tuple[float, float]]
    for iy in range(min_iy, max_iy + 1):
        y_min = seed_y + iy * tile_y_mm
        for ix in range(min_ix, max_ix + 1):
            x_min = seed_x + ix * tile_x_mm
            if corners_inside_annulus(x_min, y_min, tile_x_mm, tile_y_mm, disk.rmin_mm, disk.rmax_mm):
                placements.append((x_min, y_min))
    placements.sort(key=lambda item: (item[1], item[0]))
    return placements


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate CSV tile placements for epic_TileEndcapTracker disks."
    )
    parser.add_argument(
        "--geometry-xml",
        default="disk_layout/epic/compact/tracking/silicon_disks_outer_tiles.xml",
        help="Tile geometry XML with epic_TileEndcapTracker detectors.",
    )
    parser.add_argument(
        "--definitions-xml",
        default=None,
        help="Definitions XML used by the geometry XML. Defaults to definitions_craterlake.xml next to the geometry file.",
    )
    parser.add_argument(
        "--output",
        default="disk_layout/epic/calibrations/tracking/endcap_tiles.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--disk",
        action="append",
        default=[],
        help="Disk name to generate. Repeat to select multiple disks. Default: all tiled disks.",
    )
    parser.add_argument("--tile-x-mm", type=float, default=130.0, help="Tile width in mm.")
    parser.add_argument("--tile-y-mm", type=float, default=30.0, help="Tile height in mm.")
    parser.add_argument(
        "--x0-mm",
        type=float,
        default=None,
        help="Seed tile lower-left x in mm. Default: centered, so x0=-tile_x/2.",
    )
    parser.add_argument(
        "--y0-mm",
        type=float,
        default=None,
        help="Seed tile lower-left y in mm. Default: tile lower edge touches the disk inner radius.",
    )
    parser.add_argument("--dz-mm", type=float, default=0.0, help="Per-tile dz to write to CSV.")
    parser.add_argument(
        "--facing",
        choices=["+z", "-z"],
        default="+z",
        help="Module facing to write to CSV.",
    )
    parser.add_argument(
        "--enabled",
        default="true",
        choices=["true", "false"],
        help="Enabled flag to write to CSV.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    geometry_xml = Path(args.geometry_xml)
    definitions_xml = Path(args.definitions_xml) if args.definitions_xml else default_definitions_xml(geometry_xml)
    output_path = Path(args.output)

    disks = parse_disks(geometry_xml, definitions_xml, args.disk)
    if not disks:
        print("No matching epic_TileEndcapTracker disks found.", file=sys.stderr)
        return 1

    rows = []  # type: List[List[str]]
    for disk in disks:
        if not z_inside_layer(args.dz_mm, disk.module_thickness_mm, disk.layer_length_mm):
            print(
                f"Skipping {disk.disk}: dz={args.dz_mm:.3f} mm exceeds layer thickness "
                f"{disk.layer_length_mm:.3f} mm for module thickness {disk.module_thickness_mm:.3f} mm.",
                file=sys.stderr,
            )
            continue

        placements = generate_tiles_for_disk(disk, args.tile_x_mm, args.tile_y_mm, args.x0_mm, args.y0_mm)
        if not placements:
            print(f"Skipping {disk.disk}: no valid tile placements found.", file=sys.stderr)
            continue

        for x_min, y_min in placements:
            rows.append(
                [
                    disk.disk,
                    f"{x_min:.3f}",
                    f"{y_min:.3f}",
                    f"{args.tile_x_mm:.3f}",
                    f"{args.tile_y_mm:.3f}",
                    f"{args.dz_mm:.3f}",
                    args.facing,
                    args.enabled,
                    "generated",
                ]
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["disk", "x_min_mm", "y_min_mm", "x_size_mm", "y_size_mm", "dz_mm", "facing", "enabled", "comment"]
        )
        writer.writerows(rows)

    print(f"Wrote {len(rows)} tiles to {output_path}")
    for disk in disks:
        count = sum(1 for row in rows if row[0] == disk.disk)
        if count:
            print(
                f"{disk.disk}: {count} tiles, r=[{disk.rmin_mm:.3f}, {disk.rmax_mm:.3f}] mm, "
                f"seed=({args.x0_mm if args.x0_mm is not None else -args.tile_x_mm / 2.0:.3f}, "
                f"{args.y0_mm if args.y0_mm is not None else disk.rmin_mm:.3f}) mm"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
