"""QE data readers, orbital-resolved band analysis, and Wannier result checks."""
from __future__ import annotations
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HA_EV = 27.211386245988
BOHR_ANG = 0.529177210903
NUMBER = r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eEdD][-+]?\d+)?'


def numbers(text):
    return [float(x.replace('D', 'E').replace('d', 'e')) for x in re.findall(NUMBER, text or '')]


def read_qe_xml(path, spin_component='none', fermi_ev=None):
    tree = ET.parse(path)
    for node in tree.iter():
        node.tag = node.tag.split('}')[-1]
    out = tree.find('.//output')
    if out is None:
        raise ValueError('Expected a modern QE data-file-schema.xml with an output block')
    bs = out.find('band_structure')
    if bs is None:
        raise ValueError('QE XML has no band_structure output')
    spinors = bs.findtext('noncolin', 'false').strip().lower() in ('true', '.true.')
    lsda = bs.findtext('lsda', 'false').strip().lower() in ('true', '.true.')
    if lsda and spin_component not in ('up', 'down'):
        raise ValueError('Collinear magnetic data need spin_component up or down; run each channel separately')
    if not lsda and spin_component != 'none':
        raise ValueError('spin_component must be none for scalar or spinor calculations')
    bands, kcart = [], []
    nup = int(bs.findtext('nbnd_up', '0'))
    ndn = int(bs.findtext('nbnd_dw', str(nup)))
    for ks in bs.findall('ks_energies'):
        row = numbers(ks.findtext('eigenvalues'))
        if lsda:
            if len(row) != nup + ndn:
                raise ValueError('Unsupported LSDA eigenvalue layout in QE XML')
            row = row[:nup] if spin_component == 'up' else row[nup:]
        bands.append(row)
        kcart.append(numbers(ks.findtext('k_point')))
    eig = np.asarray(bands, dtype=float) * HA_EV
    if eig.ndim != 2 or not eig.size or not np.isfinite(eig).all():
        raise ValueError('Missing or invalid QE eigenvalues')
    fermi = numbers(bs.findtext('fermi_energy'))
    if not fermi:
        fermi = numbers(bs.findtext('two_fermi_energies'))
        if fermi:
            fermi = [fermi[1 if spin_component == 'down' else 0]]
    if fermi_ev is None and fermi:
        fermi_ev = fermi[0] * HA_EV
    if fermi_ev is None:
        # A fixed-occupation insulator may not have a Fermi-energy field.
        occ = []
        for ks in bs.findall('ks_energies'):
            row = numbers(ks.findtext('occupations'))
            if lsda:
                row = row[:nup] if spin_component == 'up' else row[nup:]
            occ.append(row)
        occ = np.asarray(occ)
        if occ.shape == eig.shape and np.any(occ > 1e-5):
            vbm = float(eig[occ > 1e-5].max())
            empty = eig[occ <= 1e-5]
            cbm = float(empty.min()) if empty.size else vbm
            fermi_ev = (vbm + cbm) / 2 if cbm > vbm else vbm
        else:
            raise ValueError('No Fermi energy: provide fermi_ev in the config')
    structure = out.find('atomic_structure')
    if structure is None:
        raise ValueError('QE XML has no atomic_structure')
    cell = np.asarray([numbers(structure.findtext('cell/a' + str(i))) for i in (1, 2, 3)])
    if cell.shape != (3, 3) or abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError('Invalid QE cell')
    atoms = [{'species': a.attrib['name'], 'cart_bohr': numbers(a.text)}
             for a in structure.findall('atomic_positions/atom')]
    for a in atoms:
        a['fractional'] = np.linalg.solve(cell.T, a['cart_bohr']).tolist()
    return {'eigenvalues': eig, 'fermi_ev': float(fermi_ev), 'spinors': spinors,
            'lsda': lsda, 'atoms': atoms, 'cell_ang': cell * BOHR_ANG,
            'kcart': np.asarray(kcart), 'nbnd': eig.shape[1], 'nk': eig.shape[0]}


def parse_projwfc(path, nk, nbnd, lsda=False, spin_component='none'):
    """Read documented stdout projections; QE truncates terms below 0.001."""
    states = {}
    energies, projections = [], []
    current_k = -1
    current_b = None
    with Path(path).open(errors='replace') as f:
        for line in f:
            m = re.search(r'state\s*#\s*(\d+):\s*atom\s*(\d+)\s*\(([^)]+)\),\s*wfc\s*(\d+)\s*\(l=\s*(\d+)', line)
            if m:
                idx, atom, species, wfc, ell = m.groups()
                states[int(idx)] = {'index': int(idx), 'atom': int(atom), 'species': species.strip(),
                                    'wfc': int(wfc), 'l': int(ell)}
            if re.match(r'\s*k\s*=', line):
                current_k += 1
                energies.append({})
                projections.append({})
                current_b = None
            m = re.search(r'====\s*e\(\s*(\d+)\)\s*=\s*(' + NUMBER + r')', line)
            if m and current_k >= 0:
                current_b = int(m[1]) - 1
                energies[-1][current_b] = numbers(m[2])[0]
                projections[-1][current_b] = {}
            if current_k >= 0 and current_b is not None:
                for weight, idx in re.findall(r'(' + NUMBER + r')\s*\*\s*\[#\s*(\d+)\]', line):
                    projections[-1][current_b][int(idx)] = numbers(weight)[0]
    if not states:
        raise ValueError('projwfc output has no atomic-state table; are pseudo-atomic orbitals available?')
    if sorted(states) != list(range(1, len(states) + 1)):
        raise ValueError('Non-contiguous projection indices')
    if lsda:
        if len(energies) != 2 * nk:
            raise ValueError('Expected both collinear spin channels in projwfc output')
        start = nk if spin_component == 'down' else 0
        energies, projections = energies[start:start + nk], projections[start:start + nk]
    if len(energies) != nk or any(len(e) != nbnd for e in energies):
        raise ValueError('projwfc k-point/band counts do not match the saved NSCF calculation')
    weights = np.zeros((nk, nbnd, len(states)))
    for k, bands in enumerate(projections):
        for band, terms in bands.items():
            for idx, weight in terms.items():
                weights[k, band, idx - 1] = weight
    return [states[i] for i in sorted(states)], weights, np.asarray([[e[b] for b in range(nbnd)] for e in energies])


def analyze(data, proj_output, output_dir, target_relative=(-2.0, 2.0), orbital_coverage=(0.7, 0.9, 0.99)):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    states, weights, printed_eig = parse_projwfc(proj_output, data['nk'], data['nbnd'], data['lsda'], data.get('spin_component', 'none'))
    eig, ef = data['eigenvalues'], data['fermi_ev']
    if np.max(abs(printed_eig - eig)) > 0.02:
        raise ValueError('Projected-band energies differ from QE XML: projection output may be stale')
    lo, hi = ef + target_relative[0], ef + target_relative[1]
    mask = (eig >= lo) & (eig <= hi)
    if not mask.any():
        raise ValueError('No calculated bands intersect target_relative_ev; change the energy range')
    # Group all m and spin components and radial channels of the same atom/l.
    group_map = {}
    for state in states:
        key = f"atom{state['atom']}_l{state['l']}"
        g = group_map.setdefault(key, {'id': key, 'atom': state['atom'], 'species': state['species'],
                                       'l': state['l'], 'indices': [], 'wfcs': []})
        g['indices'].append(state['index'])
        if state['wfc'] not in g['wfcs']:
            g['wfcs'].append(state['wfc'])
    groups = list(group_map.values())
    for g in groups:
        g['weight'] = float(weights[:, :, np.array(g['indices']) - 1][mask].sum())
        g['manual_count'] = (2 * g['l'] + 1) * (2 if data['spinors'] else 1)
        g['atomic_count'] = len(g['indices'])
    groups.sort(key=lambda g: g['weight'], reverse=True)
    total = sum(g['weight'] for g in groups)
    if total <= 0:
        raise ValueError('No usable projected weight near the target; inspect projwfc/pseudopotentials')
    candidates, seen = [], set()
    for threshold in orbital_coverage:
        chosen, accum = [], 0.0
        for g in groups:
            chosen.append(g['id']); accum += g['weight']
            if accum / total >= threshold:
                break
        key = tuple(chosen)
        if key not in seen:
            candidates.append({'id': 'orbitals_' + str(len(candidates) + 1), 'groups': chosen,
                               'projected_weight_coverage': accum / total})
            seen.add(key)
    # Always offer one extra shell where possible to repair rank-deficient guesses.
    if len(candidates[-1]['groups']) < len(groups):
        chosen = candidates[-1]['groups'] + [g['id'] for g in groups if g['id'] not in candidates[-1]['groups']][:1]
        candidates.append({'id': 'orbitals_' + str(len(candidates) + 1), 'groups': chosen,
                           'projected_weight_coverage': sum(g['weight'] for g in groups if g['id'] in chosen) / total})
    crossing = np.where(mask.any(axis=0))[0]
    first, last = int(crossing.min()), int(crossing.max())
    gap_lower = float(np.min(eig[:, first] - eig[:, first - 1])) if first else None
    gap_upper = float(np.min(eig[:, last + 1] - eig[:, last])) if last < data['nbnd'] - 1 else None
    isolated = {'first_band': first + 1, 'last_band': last + 1, 'num_wann': last - first + 1,
                'lower_direct_gap_ev': gap_lower, 'upper_direct_gap_ev': gap_upper,
                'supported': (first == 0 or gap_lower > 0.05) and gap_upper is not None and gap_upper > 0.05}
    # Tutorial 32 motivates erfc fitting; rounded stdout weights suffice for an initial guess,
    # never for the normalized AMN projectability-disentanglement test.
    fit = {'available': False, 'reason': 'Insufficient energy/projectability variation'}
    projectability = np.clip(weights.sum(axis=2).ravel(), 0, 1)
    x = eig.ravel()
    if len(x) >= 12 and np.ptp(projectability) > .15:
        from scipy.optimize import curve_fit
        from scipy.special import erfc
        import warnings
        def model(energy, mu, sigma): return .5 * erfc((energy-mu)/sigma)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                params, _ = curve_fit(model, x, projectability,
                    p0=[ef+2, 2.], bounds=([float(x.min())-10, .05], [float(x.max())+10, 20.]), maxfev=5000)
            mu, sigma = map(float, params)
            rmse = float(np.sqrt(np.mean((model(x,mu,sigma)-projectability)**2)))
            fit = {'available': rmse < .25, 'mu_fit_ev': mu, 'sigma_ev': sigma,
                   'mu_scdm_ev': mu-3*sigma, 'rmse': rmse,
                   'note': 'Initial-guess fit to truncated stdout weights; validate against target coverage'}
        except (ValueError, RuntimeError, Warning):
            fit = {'available': False, 'reason': 'Fit failed; use target-centered heuristic'}
    summary = {'fermi_ev': ef, 'energy_reference': 'absolute eV; plotting and target_relative_ev use E-Ef',
               'target_absolute_ev': [lo, hi], 'num_bands': data['nbnd'], 'num_kpoints': data['nk'],
               'spinors': data['spinors'], 'groups': groups, 'orbital_candidates': candidates,
               'isolated_candidate': isolated,
               'maximum_target_bands_at_one_k': int(mask.sum(axis=1).max()),
               'scdm_projectability_fit': fit,
               'average_projectability_in_target': float(weights.sum(axis=2)[mask].mean()),
               'top_band_below_target_top_at_some_k': bool(np.any(eig[:, -1] < hi + 1.0)),
               'projection_precision': 'QE stdout weights: 0.001 precision, omitted terms below 0.001; not used for normalized-projectability disentanglement'}
    (output_dir/'fatbands_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (output_dir/'fatbands.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['k_index', 'band_index', 'energy_ev', 'energy_minus_fermi_ev'] + [g['id'] for g in groups])
        group_weights = [weights[:, :, np.array(g['indices']) - 1].sum(axis=2) for g in groups]
        for k in range(data['nk']):
            for b in range(data['nbnd']):
                writer.writerow([k + 1, b + 1, eig[k, b], eig[k, b] - ef] + [w[k, b] for w in group_weights])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    count = min(8, len(groups))
    fig, axes = plt.subplots(count, 1, figsize=(10, max(3, count * 2.4)), squeeze=False, sharex=True)
    for i, ax in enumerate(axes[:, 0]):
        g = groups[i]
        ax.scatter(np.repeat(np.arange(data['nk']), data['nbnd']), (eig - ef).ravel(),
                   s=2 + 50 * np.clip(group_weights[i].ravel(), 0, 1), alpha=0.55)
        ax.axhline(0, color='black', lw=0.5)
        ax.set_ylim(target_relative[0] - 2, target_relative[1] + 2)
        ax.set_ylabel('E - Ef (eV)')
        ax.set_title(f"{g['species']} atom {g['atom']}, l={g['l']}")
    axes[-1, 0].set_xlabel('Uniform-mesh k-point index (not a symmetry-line band path)')
    fig.tight_layout(); fig.savefig(output_dir/'fatbands.png', dpi=140); plt.close(fig)
    return summary


def window_counts(eig, plan):
    keep = np.ones(eig.shape[1], dtype=bool)
    keep[np.asarray(plan.get('exclude_bands', []), dtype=int) - 1] = False
    e = eig[:, keep]
    outer = ((e >= plan['outer'][0]) & (e <= plan['outer'][1])).sum(axis=1)
    frozen = ((e >= plan['frozen'][0]) & (e <= plan['frozen'][1])).sum(axis=1) if plan.get('frozen') else np.zeros(len(e), int)
    return {'minimum_outer_states': int(outer.min()), 'maximum_frozen_states': int(frozen.max()),
            'retained_bands': int(keep.sum())}


def localization_result(folder, seed, plan):
    folder = Path(folder)
    wout = folder/(seed + '.wout')
    text = wout.read_text(errors='replace') if wout.exists() else ''
    err = folder/(seed + '.werr')
    errors = err.read_text(errors='replace') if err.exists() else ''
    combined = text + '\n' + errors
    final = text.rsplit('Final State', 1)[-1] if 'Final State' in text else ''
    spreads = [numbers(m)[0] for m in re.findall(r'WF centre and spread\s+\d+\s*\([^)]*\)\s*(' + NUMBER + ')', final)]
    localized = bool(re.search(r'Wannieri[sz]ation convergence criteria satisfied', text, re.I))
    disentangled = plan['num_wann'] == plan['num_bands'] or bool(re.search(r'Disentanglement convergence criteria satisfied', text, re.I))
    good_spreads = len(spreads) == plan['num_wann'] and all(np.isfinite(x) and x >= -1e-8 for x in spreads)
    artifacts = all((folder/(seed + ext)).is_file() and (folder/(seed + ext)).stat().st_size > 0 for ext in ('.chk', '_hr.dat'))
    success = localized and disentangled and good_spreads and artifacts and not errors.strip()
    low = combined.lower()
    if 'inner window' in low and ('more' in low or 'only' in low or 'too many' in low): kind = 'frozen_window'
    elif 'fewer states' in low or 'not enough states' in low or 'outer energy window' in low: kind = 'outer_window'
    elif any(x in low for x in ('zgesvd', 'dis_project', 'singular', 'projection is zero', 'not enough projections', 'too many projections')): kind = 'projection_rank'
    elif any(x in low for x in ('unrecognised', 'unrecognized', 'unknown keyword', 'not implemented', 'not supported')): kind = 'unsupported_method'
    elif any(x in low for x in ('error opening', 'cannot open', 'no such file', 'permission denied', 'k-point mismatch', 'kpoint mismatch', 'b1 is not', 'b1 vector')): kind = 'infrastructure'
    elif not disentangled: kind = 'disentanglement'
    else: kind = 'localization'
    return {'success': bool(success), 'failure_kind': None if success else kind,
            'localization_converged': localized, 'disentanglement_converged': disentangled,
            'spreads_ang2': spreads, 'maximum_spread_ang2': max(spreads) if spreads else None,
            'required_artifacts_present': artifacts, 'log_tail': combined[-16000:]}


def interpolation_check(hr_file, kpoints, reference, plan, target):
    """Training-mesh diagnostic only; never equated with an independent band-path test."""
    from scipy.optimize import linear_sum_assignment
    with Path(hr_file).open() as f:
        next(f); nw = int(next(f)); nr = int(next(f))
        degeneracies = []
        while len(degeneracies) < nr:
            degeneracies.extend(map(int, next(f).split()))
        mats = {}; ir_order = []
        for line in f:
            cols = line.split()
            if not cols: continue
            key = tuple(map(int, cols[:3]))
            if key not in mats:
                mats[key] = np.zeros((nw, nw), complex); ir_order.append(key)
            mats[key][int(cols[3])-1, int(cols[4])-1] = float(cols[5]) + 1j*float(cols[6])
    if nw != plan['num_wann'] or len(ir_order) != nr or any(x <= 0 for x in degeneracies):
        raise ValueError('Invalid Hamiltonian dimensions or degeneracies')
    errors = []; sampled = 0
    for k, ref in zip(kpoints, reference):
        hk = sum(mats[key] * np.exp(2j*np.pi*np.dot(k, key)) / degeneracies[i] for i, key in enumerate(ir_order))
        interp = np.linalg.eigvalsh((hk + hk.conj().T)/2)
        selected = ref[(ref >= target[0]) & (ref <= target[1])]
        sampled += len(selected)
        if not len(selected): continue
        rows, cols = linear_sum_assignment(abs(selected[:, None] - interp[None, :]))
        errors.extend((selected[rows]-interp[cols]).tolist())
    return {'description': 'Uniform training-mesh energy matching; independent path validation is still recommended',
            'matched_states': len(errors), 'target_states': sampled,
            'rmse_ev': float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
            'max_error_ev': float(np.max(np.abs(errors))) if errors else None}
