#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 Wouter Deconinck
"""Map the on-disk storage taken up by the collections in a podio RNTuple file.

The breakdown comes from ROOT's RNTuple treemap
(ROOT::Experimental::CreateTreeMapFromRNTuple), with the podio relation,
vector-member and subset-index fields grouped under the collection they belong
to. If a reference file is given, every collection and field is compared to it.

Writes to the output directory:
  storage.html         interactive treemap and per-collection table
  storage_treemap.svg  ROOT's own treemap rendering, colored by column type
  storage_summary.txt  one-line summary for the PR comment, e.g. "size Δ: +10.2 kB"

Kept in sync with src/scripts/storage_map.py in https://github.com/eic/EICrecon.
"""
import argparse
import html
import json
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare("#include <ROOT/RNTupleTreeMap.hxx>")

# Depth of the treemap below the collections: field, then member.
MAX_DEPTH = 2


def fmt_bytes(n, sign=False):
    for unit, scale in (("MB", 1e6), ("kB", 1e3)):
        if abs(n) >= scale:
            return f"{n / scale:{'+' if sign else ''}.1f} {unit}"
    return f"{n:{'+' if sign and n else ''}d} B"


def fmt_rel(new, ref):
    if ref is None:
        return "new"
    if new is None:
        return "removed"
    if ref == 0:
        return "" if new == 0 else "new"
    return f"{(new - ref) / ref:+.1%}"


def collection_of(field, collections):
    """Collection that a top-level podio field belongs to."""
    if not field.startswith("_"):
        return field.removesuffix("_objIdx")
    # Relations and vector members: _<collection>_<name>
    owners = [c for c in collections if field.startswith(f"_{c}_")]
    return max(owners, key=len) if owners else field


def subfields(nodes, idx):
    """Field (not column) children of a treemap node, looking through the
    '_0' item fields that RNTuple inserts for vectors and arrays."""
    node = nodes[idx]
    for c in range(node.fChildrenIdx, node.fChildrenIdx + node.fNChildren):
        child = nodes[c]
        if str(child.fType):
            continue
        if str(child.fName) == "_0":
            yield from subfields(nodes, c)
        else:
            yield c


def flatten(nodes, idx, path, depth, sizes):
    sizes[path] = nodes[idx].fSize
    if depth > 0:
        for c in subfields(nodes, idx):
            flatten(nodes, c, path + (str(nodes[c].fName),), depth - 1, sizes)


def load(path, ntuple):
    """Compressed sizes keyed by (collection, field, member, ...) paths."""
    insp = ROOT.Experimental.RNTupleInspector.Create(ntuple, path)
    treemap = ROOT.Experimental.CreateTreeMapFromRNTuple(insp)
    nodes = treemap.fNodes
    top = list(subfields(nodes, 0))
    names = [str(nodes[i].fName) for i in top]
    collections = {n.removesuffix("_objIdx") for n in names if not n.startswith("_")}

    sizes = {(): nodes[0].fSize}
    uncompressed = {}
    for idx, field in zip(top, names):
        coll = collection_of(field, collections)
        sizes[(coll,)] = sizes.get((coll,), 0) + nodes[idx].fSize
        uncompressed[coll] = (
            uncompressed.get(coll, 0)
            + insp.GetFieldTreeInspector(field).GetUncompressedSize()
        )
        flatten(nodes, idx, (coll, field), MAX_DEPTH - 1, sizes)

    return {
        "file": os.path.basename(path),
        "file_size": os.path.getsize(path),
        "entries": insp.GetDescriptor().GetNEntries(),
        "compressed": insp.GetCompressedSize(),
        "uncompressed": insp.GetUncompressedSize(),
        "sizes": sizes,
        "collection_uncompressed": uncompressed,
        "treemap": treemap,
    }


def write_root_svg(treemap, filename):
    ROOT.SetOwnership(treemap, False)
    canvas = ROOT.TCanvas("storage", "storage", 1600, 1000)
    treemap.Draw()
    canvas.SaveAs(filename)


def write_html(new, ref, ntuple, filename):
    sizes = new["sizes"]
    ref_sizes = ref["sizes"] if ref else {}

    # Nodes as parallel arrays with parent indices; labels, hover text and
    # colors are built in the browser to keep the page small.
    paths = sorted(sizes)
    index = {path: i for i, path in enumerate(paths)}
    nodes = {
        "label": [path[-1] if path else ntuple for path in paths],
        "parent": [index[path[:-1]] if path else -1 for path in paths],
        "size": [sizes[path] for path in paths],
        "ref": [ref_sizes.get(path) for path in paths] if ref else None,
    }

    # Per-collection table, including collections that were removed
    collections = sorted({p[0] for p in sizes if len(p) == 1} | {p[0] for p in ref_sizes if len(p) == 1})
    rows = []
    for coll in collections:
        size = sizes.get((coll,))
        ref_size = ref_sizes.get((coll,)) if ref else None
        unc = new["collection_uncompressed"].get(coll)
        cells = [
            (coll, coll),
            (size or 0, fmt_bytes(size) if size is not None else "–"),
            (size or 0, f"{(size or 0) / new['compressed']:.2%}"),
            ((size or 0) / max(new["entries"], 1), fmt_bytes(round((size or 0) / max(new["entries"], 1)))),
            (unc / size if unc and size else 0, f"{unc / size:.1f}" if unc and size else "–"),
        ]
        if ref:
            delta = (size or 0) - (ref_size or 0)
            cells += [
                (ref_size or 0, fmt_bytes(ref_size) if ref_size is not None else "–"),
                (delta, fmt_bytes(delta, sign=True) if delta else ""),
                (delta / ref_size if ref_size else ("Infinity" if delta else 0), fmt_rel(size, ref_size) if delta else ""),
            ]
        cls = ""
        if ref and size != ref_size:
            cls = ' class="grew"' if (size or 0) > (ref_size or 0) else ' class="shrank"'
        rows.append(
            f"<tr{cls}>"
            + "".join(
                f'<td data-sort="{html.escape(str(key))}">{html.escape(str(text))}</td>'
                for key, text in cells
            )
            + "</tr>"
        )
    headers = ["collection", "size", "share", "per event", "compression"]
    if ref:
        headers += ["reference", "Δ", "Δ (%)"]

    def summary_row(label, key, fmt=fmt_bytes):
        cells = f"<th>{label}</th><td>{fmt(new[key])}</td>"
        if ref:
            delta = new[key] - ref[key]
            cells += f"<td>{fmt(ref[key])}</td><td>{fmt(delta, sign=True)}</td>"
        return f"<tr>{cells}</tr>"

    summary = "".join([
        "<tr><th></th><th>this file</th>" + ("<th>reference</th><th>Δ</th>" if ref else "") + "</tr>",
        summary_row("file size", "file_size"),
        summary_row(f"'{html.escape(ntuple)}' compressed", "compressed"),
        summary_row(f"'{html.escape(ntuple)}' uncompressed", "uncompressed"),
        summary_row("entries", "entries", lambda n, sign=False: f"{n:{'+' if sign and n else ''}d}"),
    ])
    if ref and ref["entries"] != new["entries"]:
        summary += '<tr><td colspan="4" class="grew">Entry counts differ: sizes are not directly comparable.</td></tr>'

    title = f"Storage map of {new['file']}"
    reference = f"Compared to reference {html.escape(ref['file'])}. " if ref else "No reference file to compare to. "
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 1em 2em; color: #222; background: #fff; }}
  h1 {{ overflow-wrap: anywhere; }}
  .scroll {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; font-size: 0.9em; }}
  th, td {{ padding: 2px 8px; text-align: right; border-bottom: 1px solid #eee; }}
  th:first-child, td:first-child {{ text-align: left; }}
  #collections th {{ cursor: pointer; position: sticky; top: 0; background: #f6f6f6; }}
  tr.grew td {{ background: #fde0dd; }}
  tr.shrank td {{ background: #deebf7; }}
  #treemap {{ height: 80vh; }}
</style>
</head>
<body>
<h1>{html.escape(title)}</h1>
<p>Compressed on-disk size of the collections in the '{html.escape(ntuple)}' RNTuple, from ROOT's
RNTuple treemap with podio relation, vector-member and subset fields grouped under their collection.
{reference}Click a box to zoom in; click the top bar to zoom out.
See also <a href="storage_treemap.svg">ROOT's treemap rendering</a>, colored by column type.</p>
<div class="scroll"><table>{summary}</table></div>
<div id="treemap"></div>
<h2>Collections</h2>
<p>Compression is uncompressed / compressed size. Click a column header to sort.</p>
<div class="scroll"><table id="collections">
<thead><tr>{''.join(f'<th>{h}</th>' for h in headers)}</tr></thead>
<tbody>
{chr(10).join(rows)}
</tbody>
</table></div>
<script>
const nodes = {json.dumps(nodes, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")};
const fmt = b => Math.abs(b) >= 1e6 ? (b / 1e6).toFixed(1) + " MB"
  : Math.abs(b) >= 1e3 ? (b / 1e3).toFixed(1) + " kB" : b + " B";
const signed = (v, text) => (v > 0 ? "+" : "") + text;
const paths = [];
nodes.label.forEach((label, i) => {{
  const p = nodes.parent[i];
  paths.push(p < 0 ? label : paths[p] + "/" + label);
}});
const trace = {{
  type: "treemap",
  ids: nodes.label.map((_, i) => String(i)),
  parents: nodes.parent.map(p => p < 0 ? "" : String(p)),
  labels: nodes.label,
  values: nodes.size,
  branchvalues: "total",
  maxdepth: 3,
  texttemplate: "%{{label}}<br>%{{customdata[0]}}",
  hovertemplate: "<b>%{{customdata[1]}}</b><br>%{{customdata[0]}} (%{{percentRoot:.2%}} of total)"
    + (nodes.ref ? "<br>reference: %{{customdata[2]}}<br>Δ: %{{customdata[3]}}" : "") + "<extra></extra>",
}};
trace.customdata = nodes.size.map((s, i) => {{
  const r = nodes.ref ? nodes.ref[i] : null;
  return [fmt(s), paths[i], r === null ? "–" : fmt(r),
    r === null ? "new" : signed(s - r, fmt(s - r)) + (r ? " (" + signed(s - r, (100 * (s - r) / r).toFixed(1)) + "%)" : "")];
}});
if (nodes.ref) {{
  // Relative change on a log scale: 1% -> 1/3, 10% -> 2/3, 100% -> 1
  const scale = (s, r) => r === null || (r === 0 && s) ? 1 : r === 0 ? 0
    : Math.sign(s - r) * Math.min(Math.log10(1 + 100 * Math.abs((s - r) / r)) / 3, 1);
  trace.marker = {{
    colors: nodes.size.map((s, i) => scale(s, nodes.ref[i])),
    colorscale: "RdBu", cmin: -1, cmax: 1, cmid: 0, showscale: true,
    colorbar: {{
      title: {{text: "size Δ"}},
      tickvals: [-1, -2 / 3, -1 / 3, 0, 1 / 3, 2 / 3, 1],
      ticktext: ["≤−100%", "−10%", "−1%", "0", "+1%", "+10%", "≥+100%"],
    }},
  }};
}}
Plotly.newPlot("treemap", [trace], {{margin: {{t: 30, l: 0, r: 0, b: 0}}}}, {{responsive: true}});
const table = document.getElementById("collections");
table.querySelectorAll("th").forEach((th, col) => th.addEventListener("click", () => {{
  const desc = th.dataset.desc !== "true";
  th.dataset.desc = desc;
  const key = tr => {{
    const v = tr.children[col].dataset.sort;
    return col === 0 ? v : parseFloat(v);
  }};
  const rows = [...table.tBodies[0].rows].sort((a, b) =>
    (key(a) > key(b) ? 1 : key(a) < key(b) ? -1 : 0) * (desc ? -1 : 1));
  table.tBodies[0].append(...rows);
}}));
</script>
</body>
</html>
"""
    with open(filename, "w") as f:
        f.write(page)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", help="podio RNTuple file")
    ap.add_argument("reference", nargs="?", help="podio RNTuple file to compare to")
    ap.add_argument("-o", "--output-dir", default=".", help="directory to write to")
    ap.add_argument("--ntuple", default="events", help="RNTuple to map (default: %(default)s)")
    args = ap.parse_args()

    new = load(args.file, args.ntuple)
    ref = load(args.reference, args.ntuple) if args.reference else None

    os.makedirs(args.output_dir, exist_ok=True)
    write_html(new, ref, args.ntuple, os.path.join(args.output_dir, "storage.html"))
    write_root_svg(new["treemap"], os.path.join(args.output_dir, "storage_treemap.svg"))

    if ref:
        summary = f"size Δ: {fmt_bytes(new['compressed'] - ref['compressed'], sign=True)}"
    else:
        summary = f"size: {fmt_bytes(new['compressed'])}"
    with open(os.path.join(args.output_dir, "storage_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
