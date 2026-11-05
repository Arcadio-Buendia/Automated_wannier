#!/usr/bin/env python3
"""Jev-directed QE -> orbital analysis -> Wannier90 localization workflow."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
import numpy as np
from extensions.run_logging import timestamp, log
from wannier_analysis import read_qe_xml, analyze, window_counts, localization_result, interpolation_check

HOME = Path(__file__).resolve().parents[1]


def dump(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def clean_input(text):
    # Strip Fortran comments only outside quoted strings; preserve paths with !.
    return re.sub(r"'[^']*'|\"[^\"]*\"|![^\n]*", lambda m: '' if m[0].startswith('!') else m[0], text)


def namelist(text, name):
    # Slash terminators must be outside strings (absolute QE paths contain /).
    m = re.search(r'&' + re.escape(name) + r'\b', text, re.I)
    if not m: return None
    quoted = None
    for i in range(m.end(), len(text)):
        c = text[i]
        if quoted:
            if c == quoted: quoted = None
        elif c in "'\"": quoted = c
        elif c == '/': return (m.start(), m.end(), i, i + 1)
    raise ValueError('Unterminated namelist: ' + name)


def get_qe(text, name, key, default=None):
    block = namelist(text, name)
    if block is None: return default
    match = re.search(r'(?<!\w)' + re.escape(key) + r"\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)", text[block[1]:block[2]], re.I)
    return match[1].strip("'\"") if match else default


def set_qe(text, name, values):
    block = namelist(text, name)
    if block is None:
        text = '&' + name + '\n/\n' + text
        block = namelist(text, name)
    body = text[block[1]:block[2]]
    for key, value in values.items():
        pattern = r'(?<!\w)' + re.escape(key) + r"\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)\s*,?"
        body = re.sub(pattern, '', body, flags=re.I)
        body += f'\n  {key} = {value},\n'
    return text[:block[1]] + body + text[block[2]:]


def quote(value):
    if "'" in str(value) or '\n' in str(value):
        raise ValueError('QE paths/labels containing quotes or newlines are unsupported')
    return "'" + str(value) + "'"


def mesh_points(mesh):
    return np.array([[i/mesh[0], j/mesh[1], k/mesh[2]] for i in range(mesh[0]) for j in range(mesh[1]) for k in range(mesh[2])])


def make_nscf(original, scratch, prefix, pseudo, nbnd, points, structure=None):
    text = clean_input(original)
    for unused in ('IONS', 'CELL'):
        block=namelist(text,unused)
        if block: text=text[:block[0]]+text[block[3]:]
    text = set_qe(text, 'CONTROL', {'calculation': "'nscf'", 'restart_mode': "'from_scratch'",
                                   'prefix': quote(prefix), 'outdir': quote(scratch), 'pseudo_dir': quote(pseudo), 'disk_io': "'low'"})
    text = set_qe(text, 'SYSTEM', {'nbnd': str(nbnd), 'nosym': '.true.', 'noinv': '.true.'})
    text = set_qe(text, 'ELECTRONS', {'startingpot': "'file'", 'startingwfc': "'atomic+random'"})
    # Delete existing K_POINTS card, preserving following cards.
    lines = text.splitlines(keepends=True)
    start = next((i for i, line in enumerate(lines) if re.match(r'^\s*K_POINTS\b', line, re.I)), None)
    if start is not None:
        mode = lines[start].lower()
        if 'gamma' in mode: end = start + 1
        elif 'automatic' in mode: end = start + 2
        else:
            row = start + 1
            while row < len(lines) and not lines[row].strip(): row += 1
            n = int(lines[row].split()[0]); end = row + n + 1
        lines[start:end] = []
    if structure is not None:
        # The saved final geometry takes precedence over the initial relax/vc-relax input.
        for card, count in [('ATOMIC_POSITIONS', len(structure['atoms'])), ('CELL_PARAMETERS', 3)]:
            start_card = next((i for i, line in enumerate(lines) if re.match(r'^\s*'+card+r'\b',line,re.I)),None)
            if start_card is not None:
                end = start_card+1; rows = 0
                while end<len(lines) and rows<count:
                    if lines[end].strip(): rows+=1
                    end+=1
                del lines[start_card:end]
        text=''.join(lines)
        block=namelist(text,'SYSTEM')
        body=text[block[1]:block[2]]
        body=re.sub(r'(?<![\w])(?:celldm\s*\(\s*\d+\s*\)|cosAB|cosAC|cosBC|A|B|C)\s*=\s*[^,\s]+\s*,?', '', body, flags=re.I)
        text=text[:block[1]]+body+text[block[2]:]
        text=set_qe(text,'SYSTEM',{'ibrav':'0'})
        text += '\nCELL_PARAMETERS angstrom\n' + ''.join(' '.join(f'{v:.14f}' for v in row)+'\n' for row in structure['cell_ang'])
        text += '\nATOMIC_POSITIONS crystal\n' + ''.join(a['species']+' '+' '.join(f'{v:.14f}' for v in a['fractional'])+'\n' for a in structure['atoms'])
        lines=text.splitlines(keepends=True)
    text = ''.join(lines).rstrip() + '\nK_POINTS crystal\n' + str(len(points)) + '\n'
    return text + ''.join(' '.join(f'{v:.14f}' for v in p) + f' {1/len(points):.14f}\n' for p in points)


def saved_mesh_points(xml, data, expected_mesh=None):
    """Validate a complete tensor-product mesh while retaining the QE point order."""
    import xml.etree.ElementTree as ET
    tree=ET.parse(xml)
    for node in tree.iter(): node.tag=node.tag.split('}')[-1]
    alat=float(tree.find('.//output/atomic_structure').attrib['alat'])
    points=data['kcart'] @ (data['cell_ang']/0.529177210903).T / alat
    wrapped=np.mod(np.round(np.mod(points,1.),8),1.)
    axes=[np.unique(wrapped[:,i]) for i in range(3)]
    mesh=[len(axis) for axis in axes]
    if math.prod(mesh)!=len(points) or len(np.unique(wrapped,axis=0))!=len(points):
        raise ValueError('Existing NSCF must contain a complete uniform mesh, not a band path or symmetry-reduced mesh')
    for axis in axes:
        spacing=np.diff(np.r_[axis,axis[0]+1.])
        if not np.allclose(spacing,1./len(axis),atol=3e-7,rtol=0):
            raise ValueError('Existing NSCF k points are not uniformly spaced')
    if expected_mesh is not None and list(expected_mesh)!=mesh:
        raise ValueError(f'Existing NSCF mesh {mesh} differs from requested mesh {expected_mesh}; provide its full matching mesh')
    return mesh,points


class Jev:
    def __init__(self, config, rules, work):
        self.config, self.rules, self.work = config, rules, work
        self.key = os.environ.get('TYPESAFE_API_KEY', '')
        envfile = config['_base']/'.env'
        if not self.key and envfile.is_file():
            for line in envfile.read_text().splitlines():
                if line.startswith('TYPESAFE_API_KEY='):
                    self.key = line.split('=', 1)[1].strip().strip("'\"")
        if not self.key:
            raise ValueError('Set TYPESAFE_API_KEY in environment or .env beside the config')

    def choose(self, state, candidates, index):
        if not candidates: return None
        question = copy.deepcopy(self.rules['questions']['decision'])
        question['criteria'] = {c['id']: c['reason'] for c in candidates}
        question['criteria']['stop'] = 'No offered action is appropriate, or a prerequisite needs user intervention.'
        log(f'Requesting Jev decision {index}')
        payload = {'model': self.config.get('model', 'jev-latest'),
                   'state': dict(state, requested_at=timestamp(), rules=self.rules['decision_rules'], candidates=candidates),
                   'questions': {'decision': question}}
        dump(self.work/f'decision_{index:03d}_request.json', payload)
        req = urllib.request.Request(self.config.get('typesafe_url', 'https://api.typesafe.ai/v1/systemone'),
                                     data=json.dumps(payload, allow_nan=False).encode(),
                                     headers={'Authorization': 'Bearer ' + self.key, 'Content-Type': 'application/json'})
        response = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=120) as stream:
                    response = json.load(stream)
                break
            except urllib.error.HTTPError as exc:
                if exc.code not in (429, 500, 502, 503, 504) or attempt == 2:
                    raise RuntimeError(f'TypeSafe HTTP {exc.code}; no decision applied') from exc
            except (urllib.error.URLError, TimeoutError) as exc:
                if attempt == 2: raise RuntimeError('TypeSafe unavailable; no decision applied') from exc
            log(f'Jev request failed temporarily; retrying in {2 ** attempt} seconds')
            time.sleep(2 ** attempt)
        dump(self.work/f'decision_{index:03d}_response.json', response)
        choice = response.get('answers', {}).get('decision', {}).get('choice')
        log(f'Jev decision {index}: {choice}')
        if choice == 'stop': return None
        chosen = next((c for c in candidates if c['id'] == choice), None)
        if chosen is None:
            raise ValueError('Jev returned an unknown candidate; no arbitrary action will be executed')
        return chosen


def orbital_sets(summary):
    groups = {g['id']: g for g in summary['groups']}
    return [(c, [groups[k] for k in c['groups']]) for c in summary['orbital_candidates']]


def prepare_plan(method, selection, summary, data, defaults):
    # QE's SCDM implementation chooses pivots at the explicit Gamma point.
    if method.startswith('scdm') and not np.any(np.all(np.abs(data['kcart']) < 1e-8, axis=1)):
        return None
    selected = next(groups for c, groups in orbital_sets(summary) if c['id'] == selection)
    nw = sum(g['manual_count'] if method == 'manual' else g['atomic_count'] for g in selected)
    if method.startswith('scdm'): nw = sum(g['manual_count'] for g in selected)
    ef = data['fermi_ev']; target = summary['target_absolute_ev']
    outer = [target[0] - defaults['outer_padding_ev'], target[1] + defaults['outer_padding_ev']]
    low = [b+1 for b in range(data['nbnd']) if np.max(data['eigenvalues'][:, b]) < outer[0]]
    exclude = low
    if method == 'bloch_phases':
        iso = summary['isolated_candidate']
        if not iso['supported']: return None
        nw = iso['num_wann']
        exclude = [i for i in range(1, data['nbnd']+1) if i < iso['first_band'] or i > iso['last_band']]
    nb = data['nbnd'] - len(exclude)
    if nb < nw or nw <= 0: return None
    plan = {'method': method, 'orbital_set': selection, 'num_wann': nw, 'num_bands': nb,
            'nscf_nbnd': data['nbnd'], 'exclude_bands': exclude, 'outer': outer, 'frozen': list(target),
            'scdm_entanglement': 'isolated' if nw == nb else ('gaussian' if method == 'scdm_gaussian' else 'erfc'),
            'scdm_mu': target[1] if method != 'scdm_gaussian' else sum(target)/2,
            'scdm_sigma': max(0.5, (target[1]-target[0])/2),
            'num_iter': defaults['num_iter'], 'dis_num_iter': defaults['dis_num_iter'],
            'trial_step': 2.0, 'dis_mix_ratio': 0.5, 'zona': 1.0,
            'search_shells': defaults.get('search_shells',36),
            'conv_tol': defaults['conv_tol'], 'conv_window': 5,
            'dis_conv_tol': defaults['dis_conv_tol'], 'projectability': False}
    fit = summary.get('scdm_projectability_fit', {})
    if method == 'scdm_erfc' and fit.get('available'):
        # Avoid a fitted cutoff that removes the user's target manifold.
        plan['scdm_mu'] = max(target[1], fit['mu_scdm_ev'])
        plan['scdm_sigma'] = min(8., max(.1, fit['sigma_ev']))
    if nw == nb:
        plan['outer'] = [float(data['eigenvalues'].min()) - 1, float(data['eigenvalues'].max()) + 1]
        plan['frozen'] = None
    else:
        repair_windows(plan, data['eigenvalues'], ef)
    return plan


def repair_windows(plan, eig, ef):
    if plan['num_wann'] == plan['num_bands']: return
    # Ensure the outer window has enough states at every point, without guessing from Gamma alone.
    if window_counts(eig, plan)['minimum_outer_states'] < plan['num_wann']:
        plan['outer'] = [float(eig.min()) - 0.1, float(eig.max()) + 0.1]
    if plan.get('frozen'):
        plan['frozen'] = [max(plan['outer'][0], plan['frozen'][0]), min(plan['outer'][1], plan['frozen'][1])]
        for _ in range(60):
            if window_counts(eig, plan)['maximum_frozen_states'] <= plan['num_wann']: break
            lo, hi = plan['frozen']; center = min(hi, max(lo, ef))
            plan['frozen'] = [center + (lo-center)*0.8, center+(hi-center)*0.8]
        if window_counts(eig, plan)['maximum_frozen_states'] > plan['num_wann'] or plan['frozen'][1] - plan['frozen'][0] < 1e-5:
            plan['frozen'] = None


def validate_plan(plan, data):
    if plan['num_bands'] != data['nbnd'] - len(set(plan['exclude_bands'])):
        raise ValueError('num_bands/exclude_bands mismatch')
    if not 0 < plan['num_wann'] <= plan['num_bands']:
        raise ValueError('Invalid Wannier/band dimensions')
    counts = window_counts(data['eigenvalues'], plan)
    if plan['num_wann'] < plan['num_bands']:
        if counts['minimum_outer_states'] < plan['num_wann'] or counts['maximum_frozen_states'] > plan['num_wann']:
            raise ValueError('Window constraints fail at one or more k-points')
        if plan['frozen'] and not plan['outer'][0] <= plan['frozen'][0] <= plan['frozen'][1] <= plan['outer'][1]:
            raise ValueError('Frozen window outside outer window')
    if plan['method'] == 'bloch_phases' and plan['num_wann'] != plan['num_bands']:
        raise ValueError('Bloch phases cannot be used with disentanglement')
    return counts


def generate_inputs(folder, seed, plan, data, summary, config, scratch, prefix, points):
    validate_plan(plan, data)
    selected = next(groups for c, groups in orbital_sets(summary) if c['id'] == plan['orbital_set'])
    lines = [f"num_wann = {plan['num_wann']}", f"num_bands = {plan['num_bands']}",
             f"spinors = {'.true.' if data['spinors'] else '.false.'}", 'write_hr = .true.', 'write_xyz = .true.',
             f"fermi_energy = {data['fermi_ev']:.12g}", 'mp_grid = ' + ' '.join(map(str, config['mesh']))]
    if plan['exclude_bands']: lines.append('exclude_bands = ' + ','.join(map(str, plan['exclude_bands'])))
    for key in ('num_iter', 'dis_num_iter', 'conv_tol', 'conv_window', 'dis_conv_tol', 'dis_mix_ratio', 'search_shells'):
        lines.append(f'{key} = {plan[key]}')
    if 'fixed_step' in plan: lines.append(f"fixed_step = {plan['fixed_step']}")
    else: lines.append(f"trial_step = {plan['trial_step']}")
    if plan['num_bands'] > plan['num_wann']:
        lines.extend([f"dis_win_min = {plan['outer'][0]:.12g}", f"dis_win_max = {plan['outer'][1]:.12g}"])
        if plan['frozen']:
            lines += [f"dis_froz_min = {plan['frozen'][0]:.12g}", f"dis_froz_max = {plan['frozen'][1]:.12g}"]
    if plan.get('projectability'):
        lines += ['dis_froz_proj = .true.', f"dis_proj_min = {plan.get('projectability_min', .01)}", f"dis_proj_max = {plan.get('projectability_max', .95)}"]
    pp = {'prefix': quote(prefix), 'outdir': quote(scratch), 'seedname': quote(seed),
          'spin_component': quote(config.get('spin_component', 'none')), 'write_mmn': '.true.', 'write_amn': '.true.'}
    if plan['method'] == 'manual':
        lines += ['guiding_centres = .true.', 'begin projections']
        for group in selected:
            coord = data['atoms'][group['atom']-1]['fractional']
            site = 'f=' + ','.join(f'{v:.12g}' for v in coord)
            suffix = '(u,d)' if data['spinors'] else ''
            lines.append(f"{site}:l={group['l']}:zona={plan['zona']}{suffix}")
        lines.append('end projections')
    elif plan['method'] == 'bloch_phases':
        lines.append('use_bloch_phases = .true.'); pp['write_amn'] = '.false.'
    else:
        lines.append('auto_projections = .true.')
        if plan['method'] == 'atomic':
            pp.update(atom_proj='.true.', atom_proj_ortho='.true.')
            keep = {i for g in selected for i in g['indices']}
            excluded = [i for i in range(1, sum(g['atomic_count'] for g in summary['groups'])+1) if i not in keep]
            for index, value in enumerate(excluded, 1): pp[f'atom_proj_exclude({index})'] = str(value)
        else:
            pp.update(scdm_proj='.true.', scdm_entanglement=quote(plan['scdm_entanglement']))
            if plan['scdm_entanglement'] != 'isolated':
                pp.update(scdm_mu=str(plan['scdm_mu']), scdm_sigma=str(plan['scdm_sigma']))
    lines += ['begin unit_cell_cart', 'ang']
    lines += [' '.join(f'{v:.14f}' for v in row) for row in data['cell_ang']]
    lines += ['end unit_cell_cart', 'begin atoms_frac']
    lines += [a['species'] + ' ' + ' '.join(f'{v:.14f}' for v in a['fractional']) for a in data['atoms']]
    lines += ['end atoms_frac', 'begin kpoints']
    lines += [' '.join(f'{v:.14f}' for v in row) for row in points]
    lines += ['end kpoints']
    (folder/(seed+'.win')).write_text('\n'.join(lines)+'\n')
    (folder/(seed+'.pw2wan')).write_text('&INPUTPP\n' + ''.join(f'  {k} = {v},\n' for k,v in pp.items()) + '/\n')
    dump(folder/'plan.json', dict(plan, window_counts=window_counts(data['eigenvalues'], plan)))


def matrix_signature(plan):
    return fingerprint({k: plan.get(k) for k in ('method','orbital_set','num_wann','num_bands','nscf_nbnd','exclude_bands','scdm_mu','scdm_sigma','scdm_entanglement','zona','search_shells')})


def recovery_candidates(plan, result, summary, data, rules, tried, method_history):
    local = []
    limits = rules['limits']
    method_records = [x for x in method_history if x.get('method') == plan['method']]
    def already(action, limit=None):
        limit = limits.get('repairs_per_action_per_method', 1) if limit is None else limit
        return sum(x.get('action') == action for x in method_records) >= limit
    def add(action, reason, **updates):
        p = copy.deepcopy(plan); p.update(updates)
        if fingerprint(p) in tried or already(action): return
        try: validate_plan(p, data)
        except ValueError: return
        local.append({'id': action, 'reason': reason, 'plan': p})
    kind = result['failure_kind']
    if kind == 'infrastructure': return []
    if kind == 'neighbor_shells':
        if plan.get('search_shells',36)<200:
            add('expand_neighbor_search', 'Increase neighbor-shell search for elongated cells without bypassing B1 consistency.', search_shells=200)
        return local
    if kind == 'unsupported_method':
        if plan.get('projectability'):
            add('disable_projectability', 'Keep pseudo-atomic initialization; this build cannot use projectability disentanglement.', projectability=False)
        else:
            local = []
    else:
        if plan.get('projectability'):
            if kind == 'frozen_window' and plan.get('projectability_max', .95) < .99:
                add('raise_projectability_freeze', 'Freeze fewer high-projectability states before changing orbitals.', projectability_max=.99)
            if kind in ('frozen_window','outer_window','disentanglement'):
                add('disable_projectability', 'Keep atomic projectors and return to energy-window disentanglement.', projectability=False)
        if plan['num_wann'] < plan['num_bands']:
            eig = data['eigenvalues']; full = [float(eig.min())-.1, float(eig.max())+.1]
            if plan['outer'] != full and kind in ('outer_window','disentanglement','projection_rank','localization'):
                add('widen_outer', 'Keep method and orbitals; include all computed candidate states before changing approach.', outer=full)
            if plan['frozen'] and kind in ('frozen_window','disentanglement','projection_rank','localization'):
                lo, hi = plan['frozen']; center = (lo+hi)/2
                if hi-lo > 0.2:
                    add('shrink_frozen', 'Keep method; halve frozen range to release overconstrained states.', frozen=[center+(lo-center)/2, center+(hi-center)/2])
                add('release_frozen', 'Keep method and target reporting; remove frozen constraint to test subspace flexibility.', frozen=None)
            if kind in ('disentanglement','projection_rank') and plan['dis_mix_ratio'] > .125:
                add('damp_disentanglement', 'Stabilize subspace optimization before switching initialization.', dis_mix_ratio=plan['dis_mix_ratio']/2)
        if kind in ('projection_rank','disentanglement','localization'):
            for c, groups in orbital_sets(summary):
                fresh = prepare_plan(plan['method'], c['id'], summary, data, rules['defaults'])
                if fresh is not None and c['id'] != plan['orbital_set'] and plan['method'] != 'bloch_phases' and not already('repair_orbitals', limits.get('orbital_repairs_per_method', 2)):
                    # Keep numerical learning when repairing the orbital manifold.
                    for key in ('num_iter','dis_num_iter','trial_step','dis_mix_ratio'):
                        fresh[key] = plan[key]
                    if fingerprint(fresh) not in tried:
                        local.append({'id':'repair_'+c['id'], 'family':'repair_orbitals', 'reason':'Same initialization; change projected shell set and corresponding num_wann from fat-band weights.', 'plan':fresh})
            if plan['method'] == 'manual':
                for zona in (.5, 2.):
                    if zona != plan['zona']: add('radial_'+str(zona), 'Same orbitals/method; change trial-orbital radial extent.', zona=zona)
            if plan['method'].startswith('scdm') and plan['scdm_entanglement'] != 'isolated':
                if plan['scdm_sigma'] < limits['scdm_sigma_max_ev']:
                    add('tune_scdm_sigma', 'Keep SCDM; broaden its selection function to include missing hybridized states.', scdm_sigma=min(limits['scdm_sigma_max_ev'], 2*plan['scdm_sigma']))
                if plan['scdm_mu'] < summary['target_absolute_ev'][1] + 2:
                    add('tune_scdm_mu', 'Keep SCDM; raise selection center by 1 eV.', scdm_mu=plan['scdm_mu']+1)
        if kind in ('localization','disentanglement'):
            key = 'dis_num_iter' if kind == 'disentanglement' else 'num_iter'
            if plan[key] < limits['max_iterations']:
                add('extend_'+key, 'Increase iteration budget while preserving projection method.', **{key:min(limits['max_iterations'], 2*plan[key])})
            if kind == 'localization':
                if plan['trial_step'] > .25 and 'fixed_step' not in plan:
                    add('reduce_trial_step', 'Improve line search while preserving the chosen manifold.', trial_step=plan['trial_step']/2)
                if 'fixed_step' not in plan: add('use_fixed_step', 'Try a fixed minimizer step after line-search difficulty.', fixed_step=.1)
        if (kind in ('outer_window','projection_rank','disentanglement') or summary['top_band_below_target_top_at_some_k']) and data['nbnd'] < limits['max_nscf_bands']:
            request = min(limits['max_nscf_bands'], max(data['nbnd']+8, math.ceil(data['nbnd']*1.25)))
            tag = f"increase_bands_{request}"
            if not already('increase_bands'):
                local.append({'id':tag, 'reason':'Keep initialization; rerun NSCF/projections with more unoccupied bands.', 'increase_bands':request, 'family':'increase_bands'})
        if (plan['method']=='atomic' and rules['capabilities'].get('projectability_disentanglement')
                and not plan['projectability'] and plan['num_bands']>plan['num_wann']):
            add('enable_projectability', 'Use normalized pseudo-atomic AMN for projectability disentanglement; requires configured build support.', projectability=True, frozen=None)
    if local: return local
    # No method switches are offered while an applicable untried local repair remains.
    switches = []
    used_methods = {x.get('method') for x in method_history}
    disabled = {x.get('method') for x in method_history if x.get('failure_kind')=='unsupported_method'}
    for method in rules['methods']:
        if not rules['methods'][method]['enabled'] or method == plan['method'] or method in used_methods or method in disabled: continue
        for c, _ in orbital_sets(summary):
            fresh = prepare_plan(method, c['id'], summary, data, rules['defaults'])
            if fresh and fingerprint(fresh) not in tried:
                switches.append({'id':'switch_'+method+'_'+c['id'], 'reason':rules['methods'][method]['when'], 'plan':fresh})
            if method == 'bloch_phases': break
    return switches


class Workflow:
    def __init__(self, config_path):
        path = Path(config_path).resolve()
        self.config = json.loads(path.read_text())
        if not isinstance(self.config, dict):
            raise ValueError(f'Configuration {path} must contain a JSON object')
        self.config['_base'] = path.parent
        self.base = path.parent
        self.rules = json.loads((HOME/'knowledge/wannier_rules.json').read_text())
        if self.config.get('rules_file'):
            self.rules = json.loads((self.base/self.config['rules_file']).read_text())
        if self.config.get('projectability_disentanglement'):
            self.rules['capabilities']['projectability_disentanglement'] = True
        self.work = (self.base/self.config.get('work_dir', 'wannier_run')).resolve()
        scratch_setting=self.config.get('scratch_dir','./scratch')
        if not isinstance(scratch_setting,str) or not scratch_setting.strip():
            raise ValueError('scratch_dir must be a nonempty directory path')
        self.scratch_root=(self.base/Path(scratch_setting).expanduser()).resolve()
        self.seed = self.config.get('seed', 'wannier')
        if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', self.seed): raise ValueError('Invalid seed label')
        self.explicit_mesh = 'mesh' in self.config
        reuse=self.config.get('existing_nscf')
        if isinstance(reuse,str): reuse={'folder':reuse}
        if reuse is not None and (not isinstance(reuse,dict) or not reuse):
            raise ValueError('existing_nscf must be a folder path or an object with folder/input/output/save_dir')
        self.reuse = reuse
        if self.reuse is None:
            required=('qe_input','dft_output')
            missing=[key for key in required if not isinstance(self.config.get(key),str) or not self.config[key].strip()]
            if missing:
                raise ValueError(
                    f'Configuration {path} is missing: {", ".join(missing)}. '
                    'To reuse completed NSCF, add existing_nscf with folder/input/output '
                    '(see wannier_config.reuse.example.json). '
                    'To generate new NSCF, set both qe_input and dft_output. '
                    'Use --config to select the intended configuration file.')
        self.mesh = self.config.get('mesh', [4,4,4])
        if len(self.mesh)!=3 or any(type(x) is not int or x<1 for x in self.mesh): raise ValueError('mesh must contain three positive integers')
        target=self.config.get('target_relative_ev', self.rules['defaults']['target_relative_ev'])
        if len(target)!=2 or not all(isinstance(x,(int,float)) and math.isfinite(x) for x in target) or target[0]>=target[1]:
            raise ValueError('target_relative_ev must be two finite increasing energies')
        budget=self.config.get('max_attempts',self.rules['limits']['max_attempts'])
        if type(budget) is not int or budget<1: raise ValueError('max_attempts must be positive')
        self.config['mesh'] = self.mesh
        self.points = mesh_points(self.mesh)
        self.history = []; self.tried = set(); self.cache = {}
        self.owns_lock = False; self.initialized = False
        self.spin = self.config.get('spin_component', 'none')

    def run_command(self, name, cwd, stdin=None, args=(), output=None):
        defaults = {'pw':['pw.x'], 'projwfc':['projwfc.x'], 'pw2wannier':['pw2wannier90.x'], 'wannier':['wannier90.x']}
        command = self.config.get('commands', {}).get(name, defaults[name])
        if not isinstance(command, list) or not command or any(not isinstance(a,str) for a in command):
            raise ValueError('Each commands entry must be a nonempty JSON argument array')
        output_path = Path(cwd)/(output or name+'.out')
        started = time.monotonic()
        event = {'stage': name, 'cwd': str(cwd), 'output': str(output_path), 'started_at': timestamp()}
        def record_event():
            with (self.work/'stage_events.jsonl').open('a') as stream:
                stream.write(json.dumps(event, allow_nan=False)+'\n')
        record_event()
        log(f'Starting {name}; output: {output_path}')
        infile = None
        try:
            infile = Path(stdin).open('rb') if stdin else subprocess.DEVNULL
            with output_path.open('wb') as out:
                completed = subprocess.run(command + list(args), cwd=cwd, stdin=infile, stdout=out, stderr=subprocess.STDOUT,
                                           timeout=self.config.get('stage_timeout_seconds', 86400), check=False)
            event['returncode'] = completed.returncode
            return completed.returncode
        except Exception as exc:
            event['error'] = type(exc).__name__+': '+str(exc)
            raise
        finally:
            if hasattr(infile,'close'): infile.close()
            event['finished_at'] = timestamp()
            event['elapsed_seconds'] = round(time.monotonic()-started, 3)
            record_event()
            outcome = ('exit '+str(event['returncode'])) if 'returncode' in event else event.get('error', 'interrupted')
            log(f'Finished {name}: {outcome}; elapsed {event["elapsed_seconds"]:.3f}s')

    def prepare(self):
        self.work.mkdir(parents=True, exist_ok=True)
        self.lock = self.work/'.runner.lock'
        try: self.lock.mkdir()
        except FileExistsError: raise ValueError('Runner lock exists; check for an active process before removing ' + str(self.lock))
        self.owns_lock = True
        (self.lock/'pid').write_text(str(os.getpid()))
        if (self.work/'history.json').exists():
            raise ValueError('work_dir already contains a run. Choose a new work_dir to preserve its provenance.')
        if self.reuse:
            self.input_base=(self.base/self.reuse.get('folder','.')).resolve()
            input_path=(self.input_base/self.reuse.get('input','nscf.in')).resolve()
            self.output_path=(self.input_base/self.reuse.get('output','nscf.out')).resolve()
        else:
            self.input_base=self.base
            input_path=(self.base/self.config['qe_input']).resolve()
            self.output_path=(self.base/self.config['dft_output']).resolve()
        self.original=input_path.read_text()
        clean=clean_input(self.original)
        self.prefix=get_qe(clean,'CONTROL','prefix','pwscf')
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',self.prefix): raise ValueError('Invalid QE prefix')
        if self.reuse and get_qe(clean,'CONTROL','calculation','scf').lower()!='nscf':
            raise ValueError('existing_nscf.input must be the corresponding NSCF input (calculation=nscf)')
        outdir=get_qe(clean,'CONTROL','outdir',os.environ.get('ESPRESSO_TMPDIR','./'))
        source_option=(self.reuse if self.reuse else self.config).get('save_dir',str(Path(outdir)/(self.prefix+'.save')))
        source_save=(self.input_base/source_option).resolve()
        xml=source_save/'data-file-schema.xml'
        data=read_qe_xml(xml,self.spin,self.config.get('fermi_ev'))
        self.fermi=data['fermi_ev'];self.structure=data
        completed_output=self.output_path.read_text(errors='replace')
        if 'JOB DONE' not in completed_output or 'convergence NOT achieved' in completed_output:
            raise ValueError('Starting output does not show a completed converged QE calculation: '+str(self.output_path))
        if self.reuse:
            expected=self.mesh if self.explicit_mesh else None
            # An automatic input grid gives an independent check against symmetry reduction.
            automatic=re.search(r'^\s*K_POINTS\s*[({]?\s*automatic[^\n]*\n\s*(\d+)\s+(\d+)\s+(\d+)',clean,re.I|re.M)
            if expected is None and automatic: expected=[int(x) for x in automatic.groups()]
            self.mesh,self.points=saved_mesh_points(xml,data,expected)
            self.config['mesh']=self.mesh
            if not any(f.is_file() for f in source_save.rglob('wfc*')):
                raise ValueError('Existing NSCF has no saved wavefunctions inside its .save directory. Retain the complete wfc*.dat/HDF5 files; a text output alone cannot be reused.')
            requested=self.config.get('initial_nbnd',data['nbnd'])
            if requested!=data['nbnd']:
                raise ValueError('Remove initial_nbnd when reusing NSCF; its saved band count is used and Jev can request expansion later')
        if self.scratch_root==source_save or source_save in self.scratch_root.parents:
            raise ValueError('scratch_dir must be outside the source .save directory')
        self.scratch_root.mkdir(parents=True,exist_ok=True)
        # A per-run directory avoids prefix collisions when several jobs use /tmp.
        self.scratch=Path(tempfile.mkdtemp(prefix='wannier-',dir=self.scratch_root))
        log('Scratch directory: '+str(self.scratch))
        log('Copying saved calculation into this run directory.')
        destination = self.scratch/(self.prefix+'.save')
        if destination.exists(): raise ValueError('Refusing to overwrite existing scratch/save data')
        # Copy, never hard-link, QE save files: NSCF can overwrite them.
        shutil.copytree(source_save, destination, symlinks=False)
        self.xml = destination/'data-file-schema.xml'
        self.pseudo = (self.input_base/get_qe(clean, 'CONTROL', 'pseudo_dir', './')).resolve()
        self.nbnd = data['nbnd'] if self.reuse else self.config.get('initial_nbnd', data['nbnd']+self.rules['defaults']['extra_bands'])
        if type(self.nbnd) is not int or not 1 <= self.nbnd <= self.rules['limits']['max_nscf_bands']:
            raise ValueError('initial_nbnd must be a positive integer within max_nscf_bands')
        self.metadata = {'scratch_root':str(self.scratch_root),'scratch_dir':str(self.scratch),'source_save':str(source_save), 'qe_input_sha256':hashlib.sha256(self.original.encode()).hexdigest(),
                         'dft_output':str(self.output_path), 'input_path':str(input_path), 'reused_nscf':bool(self.reuse), 'config':{k:v for k,v in self.config.items() if not k.startswith('_')},
                         'rules_sha256':hashlib.sha256(json.dumps(self.rules,sort_keys=True).encode()).hexdigest()}
        dump(self.work/'provenance.json', self.metadata)
        dump(self.work/'history.json', [])
        self.initialized = True

    def bands_and_analysis(self, nbnd, reuse=False):
        folder = self.work/f'bands_{nbnd:04d}'
        if folder.exists(): raise ValueError('Band stage already exists: ' + str(folder))
        folder.mkdir()
        nscf=folder/'nscf.in'
        if reuse:
            log(f'Reusing existing NSCF ({nbnd} bands); skipping pw.x')
            nscf.write_text(self.original)
            shutil.copy2(self.output_path,folder/'pw.out')
            dump(folder/'reuse.json',{'source_save':self.metadata['source_save'],'source_output':str(self.output_path),'pw_executed':False})
        else:
            nscf.write_text(make_nscf(self.original,self.scratch,self.prefix,self.pseudo,nbnd,self.points,self.structure))
            log(f'Running NSCF ({nbnd} bands, {len(self.points)} k points)')
            code=self.run_command('pw',folder,stdin=nscf)
            output=(folder/'pw.out').read_text(errors='replace')
            if code or 'JOB DONE' not in output or 'convergence NOT achieved' in output:
                raise RuntimeError('NSCF failed; see '+str(folder/'pw.out'))
        data=read_qe_xml(self.xml,self.spin,self.fermi);data['spin_component']=self.spin
        if data['nbnd']!=nbnd or data['nk']!=len(self.points): raise ValueError('NSCF mesh/band count does not match request')
        _,actual_points=saved_mesh_points(self.xml,data,self.mesh)
        delta=actual_points-self.points
        if not np.allclose(delta,np.round(delta),atol=1e-6):
            raise ValueError('NSCF k-point coordinates/order differ from requested full mesh')
        # Preserve the exact saved coordinates/order for the Wannier interface.
        self.points=actual_points
        # Fat bands use the printed per-band weights, not an energy-broadened DOS.
        # Both flags must be false to bypass QE's separate partialdos MPI phase.
        projection = folder/'projwfc.in'
        projection.write_text('&PROJWFC\n prefix = '+quote(self.prefix)+',\n outdir = '+quote(self.scratch)+",\n lsym = .false.,\n kresolveddos = .false.,\n filpdos = 'projected',\n/\n")
        if self.run_command('projwfc', folder, stdin=projection): raise RuntimeError('projwfc failed; see '+str(folder/'projwfc.out'))
        summary = analyze(data, folder/'projwfc.out', folder, tuple(self.config.get('target_relative_ev', self.rules['defaults']['target_relative_ev'])), self.rules['defaults']['orbital_coverage'])
        self.nbnd, self.data, self.summary = nbnd, data, summary
        log('Fat-band data and plots: '+str(folder))
        dump(self.work/'latest_analysis.json', {'folder':str(folder),'summary':summary})

    def initial_candidates(self):
        options=[]
        for method, rule in self.rules['methods'].items():
            if not rule['enabled']: continue
            for selection,_ in orbital_sets(self.summary):
                plan=prepare_plan(method,selection['id'],self.summary,self.data,self.rules['defaults'])
                if plan:
                    options.append({'id':'start_'+method+'_'+selection['id'],'reason':rule['when'],'plan':plan})
                if method=='bloch_phases': break
        # Jev can request a larger initial band pool before choosing any method.
        smallest = min(sum(g['manual_count'] for g in groups) for _,groups in orbital_sets(self.summary))
        if (not options or self.summary['top_band_below_target_top_at_some_k']) and self.nbnd < self.rules['limits']['max_nscf_bands']:
            more = min(self.rules['limits']['max_nscf_bands'],max(self.nbnd+8,smallest+8))
            options.append({'id':f'initial_bands_{more}', 'reason':'Compute more bands before committing to an orbital manifold; the current pool may truncate the target.', 'increase_bands':more, 'family':'increase_bands'})
        if not options: raise ValueError('No valid initialization within the band-count limit; inspect projection states')
        return options

    def attempt(self, plan, index):
        folder=self.work/f'attempt_{index:03d}'; folder.mkdir()
        generate_inputs(folder,self.seed,plan,self.data,self.summary,self.config,self.scratch,self.prefix,self.points)
        signature=matrix_signature(plan)
        previous=self.cache.get(signature)
        stage='wannier'
        if previous:
            for ext in ('.nnkp','.mmn','.eig','.amn'):
                source=previous/(self.seed+ext)
                if source.exists():
                    try: os.link(source,folder/source.name)
                    except OSError: shutil.copy2(source,folder/source.name)
        else:
            code=self.run_command('wannier',folder,args=['-pp',self.seed],output='preprocess.out')
            if code or not (folder/(self.seed+'.nnkp')).exists():
                log=(folder/'preprocess.out').read_text(errors='replace')
                for ext in ('.werr','.wout'):
                    path=folder/(self.seed+ext)
                    if path.exists(): log+='\n'+path.read_text(errors='replace')
                low=log.lower()
                kind='neighbor_shells' if any(t in low for t in ('b-vectors','b1 condition','not enough shells','unable to satisfy b1')) else 'infrastructure'
                return {'success':False,'failure_kind':kind,'stage':'preprocess','log_tail':log[-16000:]},folder
            stage='pw2wannier'
            code=self.run_command('pw2wannier',folder,stdin=folder/(self.seed+'.pw2wan'))
            required=['.mmn','.eig']+([] if plan['method']=='bloch_phases' else ['.amn'])
            if code or any(not (folder/(self.seed+ext)).is_file() for ext in required):
                log=(folder/'pw2wannier.out').read_text(errors='replace')
                low=log.lower()
                kind='unsupported_method' if any(t in low for t in ('unrecognized','unrecognised','reading inputpp','not implemented','unknown variable')) else ('projection_rank' if any(t in low for t in ('projection','scdm','rank')) else 'infrastructure')
                return {'success':False,'failure_kind':kind,'stage':stage,'log_tail':log[-16000:]},folder
            for ext in ('.mmn', '.amn'):
                if ext == '.amn' and plan['method']=='bloch_phases': continue
                with (folder/(self.seed+ext)).open() as matrix:
                    next(matrix); header=next(matrix).split()
                    # SCDM AMN headers append floating-point mu/sigma after dimensions.
                    dims=list(map(int,header[:3]))
                valid = len(dims)==3 and dims[:2]==[plan['num_bands'],self.data['nk']]
                if ext=='.amn': valid=valid and dims[2]==plan['num_wann']
                if not valid:
                    return {'success':False,'failure_kind':'projection_rank','stage':'matrix_validation',
                            'log_tail':f'{ext} dimensions {dims} do not match requested bands/kpoints/orbitals'},folder
            self.cache[signature]=folder
        stage='wannier'
        code=self.run_command('wannier',folder,args=[self.seed])
        result=localization_result(folder,self.seed,plan)
        if not result['log_tail'].strip():
            result['log_tail']=(folder/'wannier.out').read_text(errors='replace')[-16000:]
        result['stage']=stage; result['returncode']=code
        if code: result['success']=False
        result['window_counts']=window_counts(self.data['eigenvalues'],plan)
        if result['success']:
            try:
                result['training_mesh_comparison']=interpolation_check(folder/(self.seed+'_hr.dat'),self.points,self.data['eigenvalues'],plan,self.summary['target_absolute_ev'])
            except (ValueError, OSError, StopIteration, IndexError) as exc:
                result['training_mesh_comparison']={'unavailable':str(exc)}
        dump(folder/'result.json',result)
        return result,folder

    def run(self, analyze_only=False):
        try:
            self.prepare()
            self.bands_and_analysis(self.nbnd,reuse=bool(self.reuse))
            if analyze_only:
                dump(self.work/'status.json',{'status':'analysis_complete','localized':False})
                return 0
            jev=Jev(self.config,self.rules,self.work)
            candidates=self.initial_candidates(); last=None; plan=None
            for i in range(1,self.config.get('max_attempts',self.rules['limits']['max_attempts'])+1):
                state={'fatbands':self.summary,'history':[{k:v for k,v in h.items() if k!='log_tail'} for h in self.history],'last_result':last,
                       'current_plan':plan,'stage':'initial_selection' if plan is None else 'recovery',
                       'energy_reference':'All decision/window/SCDM energies are absolute eV; Ef is provided.',
                       'method_switch_policy':'Only local repairs are offered until applicable untried repairs have been exhausted.'}
                decision=jev.choose(state,candidates,i)
                if decision is None:
                    dump(self.work/'status.json',{'status':'stopped','localized':False,'reason':'No valid action or Jev selected stop'})
                    return 2
                if 'increase_bands' in decision:
                    self.history.append({'timestamp':timestamp(),'action':decision.get('family', decision['id']),'candidate_id':decision['id'],'method':plan['method'] if plan else None,'requested_nbnd':decision['increase_bands']})
                    old_groups = next(c['groups'] for c in self.summary['orbital_candidates'] if c['id']==plan['orbital_set']) if plan else None
                    self.bands_and_analysis(decision['increase_bands'])
                    if plan is None:
                        candidates=self.initial_candidates()
                    else:
                        matched = next((c['id'] for c in self.summary['orbital_candidates'] if set(c['groups'])==set(old_groups)), self.summary['orbital_candidates'][0]['id'])
                        refreshed=prepare_plan(plan['method'],matched,self.summary,self.data,self.rules['defaults'])
                        if refreshed is None:
                            candidates=self.initial_candidates(); plan=None
                        else:
                            plan=refreshed
                            candidates=[{'id':'run_expanded_bands','reason':'Continue the same method with the expanded NSCF band pool.','plan':plan}]
                    dump(self.work/'history.json',self.history)
                    continue
                plan=decision['plan']; validate_plan(plan,self.data)
                self.tried.add(fingerprint(plan))
                log(f"Attempt {i}: {plan['method']}, {plan['num_wann']} WFs / {plan['num_bands']} bands")
                started_at=timestamp()
                result,folder=self.attempt(plan,i); last=result
                dump(folder/'result.json', result)
                self.history.append({'started_at':started_at,'finished_at':timestamp(),'attempt':i,'action':decision.get('family', decision['id']),'candidate_id':decision['id'],'method':plan['method'],'plan':plan,'folder':str(folder),**result})
                dump(self.work/'history.json',self.history)
                if result['success']:
                    dump(self.work/'status.json',{'status':'localized','localized':True,'attempt':i,'folder':str(folder),'plan':plan,
                                                'note':'Spread convergence reached; no independent-path accuracy requirement was imposed.'})
                    log('Localization reached: '+str(folder))
                    return 0
                candidates=recovery_candidates(plan,result,self.summary,self.data,self.rules,self.tried,self.history)
            dump(self.work/'status.json',{'status':'attempt_limit','localized':False})
            return 2
        except Exception as exc:
            if self.initialized:
                dump(self.work/'status.json',{'status':'error','localized':False,'error':str(exc)})
            raise
        finally:
            if self.owns_lock and self.lock.exists():
                (self.lock/'pid').unlink(); self.lock.rmdir()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default='wannier_config.json')
    parser.add_argument('--analyze-only',action='store_true',help='Run NSCF/projwfc and Python fat-band analysis without a Jev request')
    args=parser.parse_args()
    try: return Workflow(args.config).run(args.analyze_only)
    except (ValueError,RuntimeError,OSError,subprocess.TimeoutExpired,KeyError) as exc:
        log('ERROR: '+str(exc),file=sys.stderr); return 1

if __name__=='__main__': sys.exit(main())
