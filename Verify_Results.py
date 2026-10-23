#!/usr/bin/env python3
"""Read-only evidence replay and current-manuscript table checks; no fitting."""
import csv
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CLINICAL = ROOT / 'results/session/referee_2026-10-22/clinical'

def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def lse(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(x-maximum) for x in values))

def lc(n, k):
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

source = (ROOT / 'main.tex').read_text()
cells = read(ROOT / 'reference/liming_revision/data/application_fixed_margin_updates.csv')
archived = read(ROOT / 'simulations/results/application_weight_eprocess_paths.csv')
lookup = {(r['analysis'], r['weight_code'], float(r['look_time'])): float(r['log_e_value']) for r in archived}
maximum_error = 0.0
checked = 0
paths = {}
for analysis in sorted({c['analysis'] for c in cells}):
    selected = sorted((c for c in cells if c['analysis'] == analysis), key=lambda c: int(c['update_id']))
    coefficients = [-math.log(2), -math.log(1.5), math.log(1.5), math.log(2)] if analysis.startswith('cav') else list(map(math.log, [.5, .67, .8]))
    for weight in ['LR', 'GB', 'PP', 'FH', 'GT']:
        current = [0.0] * len(coefficients)
        for index, cell in enumerate(selected):
            n0, n1, m, y = (int(cell[k]) for k in ['n0', 'n1', 'd', 'x'])
            support = list(range(max(0, m-n0), min(m, n1)+1))
            null = [lc(n1, u)+lc(n0, m-u)-lc(n0+n1, m) for u in support]
            for component, eta in enumerate(coefficients):
                z = eta * float(cell['weight_'+weight])
                current[component] += z*y-lse([p+z*u for p, u in zip(null, support)])
            time = float(cell['look_time'])
            if index+1 == len(selected) or float(selected[index+1]['look_time']) != time:
                value = lse(current)-math.log(len(current))
                key = (analysis, weight, time)
                paths[key] = value
                maximum_error = max(maximum_error, abs(value-lookup[key]))
                checked += 1
    for time in sorted({float(c['look_time']) for c in selected}):
        value = lse([paths[(analysis, 'LR', time)], paths[(analysis, 'GT', time)]])-math.log(2)
        key = (analysis, 'LR--GT', time)
        maximum_error = max(maximum_error, abs(value-lookup[key]))
        paths[key] = value
        checked += 1
assert checked == len(archived) == 1140
assert maximum_error < 1e-10, maximum_error

def table(label):
    position = source.index('\\label{'+label+'}')
    begin, env = max((source.rfind('\\begin{'+env+'}', 0, position), env) for env in ['table', 'longtable'])
    end = source.index('\\end{'+env+'}', position)
    return source[begin:end]

def normalized(text):
    return re.sub(r'\s+', '', text)

def contains(label, row):
    assert normalized(row) in normalized(table(label)), (label, row)

current_rows = 0
for row in read(CLINICAL / 'comparison_counts.csv'):
    app = 'CAV' if row['dataset'] == 'cav' else 'PsA'
    width = f"{float(row['q']):.2f}" if app == 'CAV' else str(int(float(row['q'])))
    unit = 'years' if app == 'CAV' else 'months'
    expected = f"{app} & {width} {unit} & {row['intervals']} & {row['looks']} & {row['informative_comparisons']}/{row['comparisons']} & {row['informative_events']}/{row['events']}"
    contains('tab:clinical-information', expected)
    current_rows += 1

for row in read(CLINICAL / 'coefficient_weight_ranges.csv'):
    app = 'CAV' if row['analysis'].startswith('cav') else 'PsA'
    q = next(c['q'] for c in cells if c['analysis'] == row['analysis'])
    width = f"{float(q):.2f}" if app == 'CAV' else str(int(float(q)))
    all_range = f"$[{float(row['eta_weight_min']):.3f},{float(row['eta_weight_max']):.3f}]$"
    informative = f"$[{float(row['informative_eta_weight_min']):.3f},{float(row['informative_eta_weight_max']):.3f}]$"
    contains('tab:clinical-coefficient-ranges', f"{app} & {width} & {row['weight']} & {all_range} & {informative}")
    current_rows += 1

contributions = read(CLINICAL / 'psa_reference_table_contributions.csv')
for row in [r for r in contributions if r['weight'] == 'LR']:
    gt = next(r for r in contributions if r['weight'] == 'GT' and r['update_id'] == row['update_id'])
    fmt = lambda x: '0' if abs(float(x)) < 1e-10 else ('1' if abs(float(x)-1) < 1e-10 else f'{float(x):.3f}')
    window = f"$({int(float(row['left']))},{int(float(row['right']))}]$"
    effusion = 'Yes' if row['effusion'] == 'high effusion' else 'No'
    expected = f"{window} & {effusion} & {row['events0']}/{row['n0']} & {row['events1']}/{row['n1']} & {fmt(row['variance'])} & {fmt(row['mixture_multiplier'])} & {fmt(gt['mixture_multiplier'])}"
    contains('tab:psa-eight', expected)
    current_rows += 1

stress = read(CLINICAL / 'stress/summary.csv')
mechanisms = {'common_grid': 'Common grid', 'independent_jitter': 'Independent jitter', 'group_differential_jitter': 'Differential jitter'}
for row in [r for r in stress if r['weight'] == 'LR']:
    gt = next(r for r in stress if r['application'] == row['application'] and r['mechanism'] == row['mechanism'] and r['weight'] == 'GT')
    show = lambda r: f"${float(r['rate']):.4f}\\ [{float(r['wilson_lower']):.4f},{float(r['wilson_upper']):.4f}]$"
    spread = '0' if float(row['median_max_probability_range']) == 0 else f"{float(row['median_max_probability_range']):.4f}"
    contains('tab:clinical-stress', f"{row['application']} & {mechanisms[row['mechanism']]} & {show(row)} & {show(gt)} & {spread}")
    current_rows += 1

summary = read(ROOT / 'results/session/todo_production_2026-10-15/summary.csv')
null = [r for r in summary if r['truth'] == 'null' and r['method'] not in ['fixed_score_LR', 'Fay_Shaw_Sun_pclt']]
assert len(null) == 450 and max(float(r['rejection_probability']) for r in null) == .0486

result = {
    'passed': True,
    'source_sha256': hashlib.sha256((ROOT/'main.tex').read_bytes()).hexdigest(),
    'clinical_paths_independently_replayed': checked,
    'maximum_absolute_log_error_against_archived': maximum_error,
    'current_source_new_clinical_table_rows_matched': current_rows,
    'fixed_design_null_eprocess_comparisons': len(null),
    'maximum_fixed_design_null_crossing_probability': .0486,
    'scope': 'Standard-library hypergeometric replay from archived fixed margins, current clinical table rows, saved fixed-design summary; no simulation or patient preprocessing rerun.',
}
print(json.dumps(result, indent=2))
