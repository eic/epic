#!/usr/bin/env python3
"""Combine SURVEY (positions) and TWISS (strengths) TFS files to update XML beamline definitions.

Usage:
  python update_beamline_from_survey_and_twiss.py survey.tfs twiss.tfs definitions.xml [--dry-run]

Logic:
  - SURVEY file supplies physical positions (X,Y,Z, THETA) referenced to IP6
  - TWISS file supplies magnetic lengths (L), bend angles (ANGLE), and integrated gradients (K1L)
  - Center Z in EPIC frame defined as -(Z - Z_IP) for backwards-going electron beam
  - Center X defined relative to IP (X - X_IP)
  - Length from TWISS 'L'
  - For RBEND elements: BendAngle = ANGLE
  - For QUADRUPOLE elements: K1L = integrated quadrupole strength; we also compute K1 = K1L / L

Adds (if missing) new constants:
  <prefix>_K1L, <prefix>_K1 for quadrupoles
  <prefix>_BendAngle for dipoles (RBEND)

Existing constants updated when present:
  <prefix>_CenterPosition, <prefix>_Length, <prefix>_CenterZ/<prefix>_CenterX variants

"""
import re
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def parse_tfs_file(tfs_path: Path) -> pd.DataFrame:
    with open(tfs_path, 'r') as f:
        content = f.read()
    is_survey = 'SURVEY' in content[:1000]
    header_match = re.search(r'\* NAME\s+KEYWORD.*', content)
    if not header_match:
        raise ValueError(f"Header not found in {tfs_path}")
    header_line = header_match.group(0)
    columns = header_line[1:].split()
    pattern = r'^\s*"([^"]+)"\s+"([^"]+)"(.*)$'
    rows = []
    for line in content.split('\n'):
        m = re.match(pattern, line)
        if not m:
            continue
        name, keyword, rest = m.group(1), m.group(2), m.group(3).strip()
        values = rest.split()
        row = {'NAME': name, 'KEYWORD': keyword}
        for i, col in enumerate(columns[2:]):
            if i < len(values):
                try:
                    row[col] = float(values[i])
                except ValueError:
                    row[col] = np.nan
            else:
                row[col] = np.nan
        row['_is_survey'] = is_survey
        rows.append(row)
    df = pd.DataFrame(rows)

    # Create unique element_key by adding occurrence index for duplicates
    df['element_key'] = df['NAME'] + '_' + df['KEYWORD']
    df['occurrence'] = df.groupby('element_key').cumcount()

    # For duplicated keys, append occurrence number (e.g., D4ER_6_RBEND_0, D4ER_6_RBEND_1, ...)
    mask = df['element_key'].duplicated(keep=False)
    df.loc[mask, 'element_key'] = df.loc[mask, 'element_key'] + '_' + df.loc[mask, 'occurrence'].astype(str)

    df.set_index('element_key', inplace=True)
    return df


def format_value(value: float, unit: str) -> str:
    # Consistent formatting
    if abs(value) < 0.001:
        base = f"{value:.6f}"
    elif abs(value) < 1:
        base = f"{value:.4f}"
    else:
        base = f"{value:.3f}"
    # Normalize unit names
    if unit == 'T/m':
        unit = 'tesla/m'
    elif unit == 'T':
        unit = 'tesla'
    return f"{base}*{unit}" if unit else base


def compute_epic_positions(survey_df: pd.DataFrame, ip_name: str = 'IP6') -> pd.DataFrame:
    """Compute EPIC-frame positions.

    IMPORTANT: In SURVEY files, Z represents the DOWNSTREAM END of each element.
    For upstream (backwards) magnets, we need to compute the center position:
      CenterZ = Z_end + L/2  (moves upstream toward center since Z is negative)

    EPIC convention: upstream magnets have negative Z positions.

    CenterZ (EPIC) = (Z_survey_end + L/2) - Z_IP
    CenterX (EPIC) = X_survey - X_IP
    """
    survey_df = survey_df.copy()
    ip_candidates = [idx for idx in survey_df.index if ip_name in idx]
    if not ip_candidates:
        raise ValueError(f"IP reference {ip_name} not found in survey file")
    ip_row = survey_df.loc[ip_candidates[0]]
    ip_x = ip_row.get('X', 0.0)
    ip_z = ip_row.get('Z', 0.0)
    ip_theta = ip_row.get('THETA', 0.0)

    # Z in SURVEY is the downstream end position
    # Center = end + L/2 (going upstream for backwards beam)
    # Compute center in local frame
    dx = survey_df['X'] - ip_x - (survey_df['L']/2.0) * np.sin(survey_df['THETA'])
    dz = survey_df['Z'] - ip_z - (survey_df['L']/2.0) * np.cos(survey_df['THETA'])
    # Rotate global vector into EPIC/IP frame: use only IP angle
    c = np.cos(ip_theta)
    s = np.sin(ip_theta)
    survey_df['epic_center_x'] = - dx * c + dz * s
    survey_df['epic_center_z'] = - dx * s - dz * c
    # THETA is the beam angle at this element
    # Convert to EPIC frame: subtract IP theta to get relative angle
    survey_df['epic_theta'] = survey_df['THETA'] - ip_theta
    return survey_df


def merge_element_data(survey_df: pd.DataFrame, twiss_df: pd.DataFrame, mapping: dict, reference_energy_GeV: float = 18.0) -> dict:
    """Merge survey and twiss data, computing reference field strengths.

    Args:
        reference_energy_GeV: Reference electron beam energy in GeV for computing field strengths
                             Default 18 GeV (nominal IP6 electron energy)
    """
    merged = {}
    for tfs_name, prefix in mapping.items():
        # Check multiple possible keyword types
        # For indexed entries like D4ER_6_RBEND_0, the tfs_name IS the full key
        if tfs_name in survey_df.index:
            # Direct match - this is an indexed entry from the mapping
            survey_key_opts = [tfs_name]
            # Try exact match in twiss first, then try without index
            twiss_key_opts = [tfs_name]
            base_parts = tfs_name.split('_')
            if base_parts[-1].isdigit():
                # Also try without occurrence index
                base_tfs = '_'.join(base_parts[:-1])
                twiss_key_opts.append(base_tfs)
        else:
            # Not in index - try appending keyword types
            survey_key_opts = [
                f"{tfs_name}_QUADRUPOLE",
                f"{tfs_name}_RBEND",
                f"{tfs_name}_CRABCAVITY",
                f"{tfs_name}_CRABCAVITY_0",  # Handle indexed CRABCAVITY
                f"{tfs_name}_DRIFT",
                f"{tfs_name}_MARKER"
            ]
            twiss_key_opts = [
                f"{tfs_name}_QUADRUPOLE",
                f"{tfs_name}_RBEND",
                f"{tfs_name}_CRABCAVITY"
            ]

        survey_row = None
        twiss_row = None
        for k in survey_key_opts:
            if k in survey_df.index:
                survey_row = survey_df.loc[k]
                break
        for k in twiss_key_opts:
            if k in twiss_df.index:
                twiss_row = twiss_df.loc[k]
                break

        if survey_row is None:
            continue  # Need at least survey data for positions

        data = {}
        # Positions from survey
        data['CenterZ'] = survey_row['epic_center_z']
        data['CenterX'] = survey_row['epic_center_x']
        # Theta from survey
        theta = survey_row.get('epic_theta', np.nan)
        if not np.isnan(theta):
            data['Theta'] = theta
        # Length from TWISS L (fallback to survey L if missing)
        if twiss_row is not None:
            L = twiss_row.get('L', survey_row.get('L', np.nan))
        else:
            L = survey_row.get('L', np.nan)
        if not np.isnan(L):
            data['Length'] = L

        # Strengths - only if TWISS data available
        if twiss_row is not None:
            # Get keyword from whichever source has it
            keyword = ''
            if hasattr(survey_row, 'get'):
                keyword = survey_row.get('KEYWORD', '')
            if not keyword and hasattr(twiss_row, 'get'):
                keyword = twiss_row.get('KEYWORD', '')

            if 'QUADRUPOLE' in keyword:
                K1L = twiss_row.get('K1L', np.nan)
                if not np.isnan(K1L):
                    data['K1L'] = K1L
                    # Note: K1, Gradient_Ref will be calculated in XML from K1L and Length
            elif 'RBEND' in keyword:
                angle = twiss_row.get('ANGLE', survey_row.get('ANGLE', np.nan))
                if not np.isnan(angle):
                    data['BendAngle'] = angle
                    # Note: B_Ref will be calculated in XML from BendAngle and Length
            # CRABCAVITY doesn't need K1L or BendAngle - just position and length

        merged[prefix] = data
    return merged


def update_xml(xml_path: Path, merged_data: dict, dry_run: bool = False,
               constants_path: Path = None, reference_energy_GeV: float = 18.0) -> tuple[int, int]:
    """Update geometry XML (positions/angles only) and constants XML (magnet fundamental params).

    Fundamental parameters (Length, K1L, BendAngle) are written to beamline_constants.xml.
    Reference field strengths (Gradient_Ref, B_Ref) and per-GeV normalization are calculated
    transparently in XML using formulas.
    Only position/orientation (CenterX, CenterZ, Theta) go to definitions.xml.

    Returns:
        (geometry_updates, constants_updates) tuple
    """
    xml_text = xml_path.read_text()
    geometry_updates = 0
    constants_updates = 0

    # Prepare constants XML content if needed
    constants_lines = []
    if constants_path:
        if constants_path.exists():
            constants_text = constants_path.read_text()
        else:
            # Create initial file with reference energy and Brho
            constants_text = f'''<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->
<!-- Generated from TWISS/SURVEY via update_beamline_from_survey_and_twiss.py -->
<!-- Reference field strengths calculated transparently in XML from fundamental parameters -->

  <define>

    <comment>
      Reference beam energy from TWISS file (auto-detected)
      This is the energy at which the reference field strengths are calculated
      Magnetic rigidity: Brho = E[GeV] / 0.2998 [T·m]
    </comment>
    <constant name="FarBackwardMagnets_ReferenceEnergy" value="{reference_energy_GeV:.12f}*GeV"/>
    <constant name="FarBackwardMagnets_Brho" value="FarBackwardMagnets_ReferenceEnergy / 0.2998 / GeV * m"/>

    <comment>
      Magnet physical parameters from TWISS/SURVEY files
      Length, K1L, BendAngle: fundamental parameters from TWISS
      Reference field strengths calculated transparently in XML:
        Gradient_Ref = K1L * Brho / Length  (for quadrupoles)
        B_Ref = Brho * BendAngle / Length   (for dipoles)
    </comment>

  </define>
'''

    for prefix, params in merged_data.items():
        for suffix, value in params.items():
            # Determine unit
            if suffix.lower().endswith('angle') or suffix.lower() == 'theta':
                unit = 'rad'
            elif suffix.lower() == 'k1l':
                unit = '1/m'
            else:
                unit = 'm'

            const_name = f"{prefix}_{suffix}"

            # Only write fundamental magnet parameters (Length, K1L, BendAngle) to constants file
            # K1, Gradient_Ref, B_Ref are calculated in XML from these fundamental parameters
            if suffix in ['Length', 'K1L', 'BendAngle'] and constants_path:
                formatted = format_value(value, unit)

                pattern = re.compile(rf'(<constant\s+name="{const_name}"\s+value=")[^"]*(")')
                if constants_path.exists():
                    m = pattern.search(constants_text)
                    if m:
                        old_value = constants_text[m.start(1)+len(m.group(1)):m.end(2)-len(m.group(2))]
                        if dry_run:
                            print(f"[constants] Would update {const_name}: {old_value} -> {formatted}")
                        else:
                            constants_text = constants_text[:m.start(1)+len(m.group(1))] + formatted + constants_text[m.end(2)-len(m.group(2)):]
                            print(f"[constants] Updated {const_name}: {old_value} -> {formatted}")
                        constants_updates += 1
                    else:
                        # Add new constant before </define>
                        new_line = f'    <constant name="{const_name}" value="{formatted}"/>\n'
                        insert_pos = constants_text.rfind('</define>')
                        if insert_pos > 0:
                            if dry_run:
                                print(f"[constants] Would add {const_name} = {formatted}")
                            else:
                                constants_text = constants_text[:insert_pos] + new_line + '\n' + constants_text[insert_pos:]
                                print(f"[constants] Added {const_name} = {formatted}")
                            constants_updates += 1
                else:
                    constants_lines.append(f'    <constant name="{const_name}" value="{formatted}"/>')
                    constants_updates += 1
            else:
                # Position/orientation parameters go to definitions.xml
                formatted = format_value(value, unit)
                pattern = re.compile(rf'(<constant\s+name="{const_name}"\s+value=")[^"]*(")')
                m = pattern.search(xml_text)
                if m:
                    old_value = xml_text[m.start(1)+len(m.group(1)):m.end(2)-len(m.group(2))]
                    if dry_run:
                        print(f"[geometry] Would update {const_name}: {old_value} -> {formatted}")
                    else:
                        xml_text = xml_text[:m.start(1)+len(m.group(1))] + formatted + xml_text[m.end(2)-len(m.group(2)) :]
                        print(f"[geometry] Updated {const_name}: {old_value} -> {formatted}")
                    geometry_updates += 1
                else:
                    # Skip constants that don't exist in the XML (they should be pre-defined)
                    if suffix not in ['Gradient', 'B']:  # Don't warn for field strengths
                        print(f"  Warning: {const_name} not found in geometry XML, skipping")

    if not dry_run:
        if geometry_updates:
            xml_path.write_text(xml_text)
        if constants_path and (constants_updates or constants_lines):
            if not constants_path.exists() and constants_lines:
                # Create new file with collected lines
                insert_pos = constants_text.rfind('</define>')
                constants_text = constants_text[:insert_pos] + '\n'.join(constants_lines) + '\n\n' + constants_text[insert_pos:]
            constants_path.write_text(constants_text)

    return geometry_updates, constants_updates
def main():
    ap = argparse.ArgumentParser(description='Merge SURVEY and TWISS TFS files to update XML.')
    ap.add_argument('survey_file', type=Path, nargs='?',
                    default=Path('calibrations/electron_storage_ring/esr-survey-tilt_20250723.tfs'),
                    help='Path to SURVEY TFS file (default: calibrations/electron_storage_ring/esr-survey-tilt_20250723.tfs)')
    ap.add_argument('twiss_file', type=Path, nargs='?',
                    default=Path('calibrations/electron_storage_ring/esr-ip6-18GeV.tfs'),
                    help='Path to TWISS TFS file (default: calibrations/electron_storage_ring/esr-ip6-18GeV.tfs)')
    ap.add_argument('xml_file', type=Path, nargs='?',
                    default=Path('compact/far_backward/definitions.xml'),
                    help='Path to XML definitions file (default: compact/far_backward/definitions.xml)')
    ap.add_argument('--dry-run', action='store_true',
                    help='Show what would be updated without modifying files')
    ap.add_argument('--ip-ref', default='IP6',
                    help='Name of interaction point in TFS file (default: IP6)')
    ap.add_argument('--constants-file', type=Path,
                    default=Path('compact/fields/beamline_constants.xml'),
                    help='Path to output constants file for field strengths (default: compact/fields/beamline_constants.xml)')
    ap.add_argument('--reference-energy', type=float, default=None,
                    help='Reference electron beam energy in GeV for field strength calculation (default: auto-detect from TWISS file)')
    args = ap.parse_args()

    if not args.survey_file.exists() or not args.twiss_file.exists() or not args.xml_file.exists():
        raise SystemExit('One or more input files do not exist.')

    mapping = {
        'Q1ER_6': 'Q1eR',
        'Q2ER_6': 'Q2eR',
        'D2AER_6': 'B2AeR',
        'D2BER_6': 'B2BeR',
        'SQ3ER_6': 'SQ3eR',
        'Q3ER_6': 'Q3eR',
        'SQ4ER_6': 'SQ4eR',
        'Q4ER_6': 'Q4eR',
        'SQ5ER_6': 'SQ5eR',
        'Q5ER_6': 'Q5eR',
        'D3ER_6': 'B3eR',
        'SQ6ER_6': 'SQ6eR',
        'SQ7ER_6': 'SQ7eR',
        'SQ8ER_6': 'SQ8eR',
        'Q8ER_6': 'Q8eR',
        'RF_CRAB': 'RF_Crab',
        'Q9ER_6': 'Q9eR',
        'SQ10ER_6': 'SQ10eR',
        'Q10ER_6': 'Q10eR',
        'Q11ER_6': 'Q11eR',
        'Q12ER_6': 'Q12eR',
        'Q13ER_6': 'Q13eR',
    }

    print('Parsing SURVEY file...')
    survey_df = parse_tfs_file(args.survey_file)
    survey_df = compute_epic_positions(survey_df, ip_name=args.ip_ref)
    print(f"  SURVEY elements: {len(survey_df)}")

    # Handle D4ER_6 which appears 5 times - create indexed mapping by occurrence
    # The parse_tfs_file adds _0, _1, _2, etc. to duplicated element_keys
    d4er_mapping = {}
    for key in survey_df.index:
        if key.startswith('D4ER_6_RBEND_'):
            occurrence = int(key.split('_')[-1])
            d4er_mapping[key] = f'D4eR_{occurrence + 1}'  # XML uses 1-based indexing
            s_pos = survey_df.loc[key, 'S']
            print(f"  Mapping {key} (S={s_pos:.3f}) -> D4eR_{occurrence + 1}")

    # Combine mappings
    full_mapping = {**mapping, **d4er_mapping}

    print('Parsing TWISS file...')
    twiss_df = parse_tfs_file(args.twiss_file)
    print(f"  TWISS elements: {len(twiss_df)}")

    # Auto-detect beam energy from TWISS header
    detected_energy = None
    try:
        with open(args.twiss_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('@ ENERGY'):
                    parts = line.split()
                    if len(parts) >= 4:
                        detected_energy = float(parts[3])
                        print(f"  Detected beam energy from @ ENERGY: {detected_energy} GeV")
                        break
                elif line.startswith('@ E '):
                    parts = line.split()
                    if len(parts) >= 4:
                        detected_energy = float(parts[3])
                        print(f"  Detected beam energy from @ E: {detected_energy} GeV")
                        break
                elif line.startswith('@ PC'):
                    parts = line.split()
                    if len(parts) >= 4:
                        pc = float(parts[3])
                        # For ultra-relativistic particles, E ≈ pc
                        detected_energy = pc
                        print(f"  Detected beam energy from @ PC: {detected_energy} GeV (assuming E≈pc)")
                        break
    except Exception as e:
        print(f"  Warning: Could not auto-detect beam energy: {e}")

    # Use detected energy if available, otherwise fall back to argument or default
    if detected_energy is not None:
        reference_energy = detected_energy
        if args.reference_energy is not None:
            print(f"  Warning: Overriding auto-detected energy {detected_energy} GeV with --reference-energy {args.reference_energy} GeV")
            reference_energy = args.reference_energy
    else:
        reference_energy = args.reference_energy if args.reference_energy is not None else 18.0
        if reference_energy == 18.0:
            print(f"  Warning: Using default reference energy {reference_energy} GeV (could not auto-detect from TWISS file)")

    print(f'Merging element data (reference energy: {reference_energy} GeV)...')
    merged = merge_element_data(survey_df, twiss_df, full_mapping, reference_energy_GeV=reference_energy)
    for prefix, data in merged.items():
        print(f"  {prefix}: keys={list(data.keys())}")

    print('Updating XML files...')
    geo_updates, const_updates = update_xml(args.xml_file, merged,
                                            dry_run=args.dry_run,
                                            constants_path=args.constants_file,
                                            reference_energy_GeV=reference_energy)
    print(f"Done. Geometry: {geo_updates} constants {'would be' if args.dry_run else 'were'} updated.")
    print(f"      Field strengths: {const_updates} constants {'would be' if args.dry_run else 'were'} updated/added.")


if __name__ == '__main__':
    main()
