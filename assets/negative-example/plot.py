#!/usr/bin/env python3
"""Render the source-bound common-parameter adaptation comparison."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[2]
SOURCE_PINS = {
    'theory_20261005/same_chi_static_control_20261005/same_chi_control_result01.json': None,
    'theory_20261005/bulk_surplus_prediction/finite16_source_numbers.json': None,
}


def load_sources():
    manifest = json.loads((HERE / 'source_pins.json').read_text())
    values = {}
    for name in SOURCE_PINS:
        path = RESEARCH / name
        raw = path.read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        if actual != manifest[name]['sha256']:
            raise ArithmeticError('Frozen figure source changed: ' + name)
        values[name] = json.loads(raw)
    macro, finite = (values[name] for name in SOURCE_PINS)
    assert macro['lambda_transfer'] == finite['parameters']['lambda_transfer'] == .70
    assert macro['mu'] == finite['parameters']['mu'] == -18.10
    assert macro['chi'] == 512 and macro['dimensions'] == [8, 6, 4]
    return macro, finite


def waterfall(ax, parts, labels, colors, limits):
    base = 0.0
    for x, (value, color) in enumerate(zip(parts, colors)):
        ax.bar(x, abs(value), bottom=min(base, base + value), width=.67,
               color=color, edgecolor='none')
        end = base + value
        if x < len(parts) - 1:
            ax.plot([x+.34, x+.66], [end, end], color='#88939C', lw=.65)
        offset = .058 * (limits[1]-limits[0])
        ax.text(x, max(base, end)+offset if value >= 0 else min(base, end)-offset,
                f'{value:+.3f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=7)
        base = end
    x = len(parts)
    ax.bar(x, abs(base), bottom=min(0, base), width=.67, color='#234953', edgecolor='none')
    ax.text(x, min(0, base)-.058*(limits[1]-limits[0]), f'{base:+.3f}',
            ha='center', va='top', fontsize=7, fontweight='bold')
    ax.axhline(0, color='#A9B0B4', lw=.6, zorder=0)
    ax.set_xlim(-.58, x+.58)
    ax.set_ylim(*limits)
    ax.set_xticks(range(x+1), labels + ['Net'])
    ax.tick_params(axis='x', length=0, pad=7)
    ax.set_ylabel('Energy difference')


def main(a):
    macro, finite = load_sources()
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6.7,
        'ytick.labelsize': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
        'legend.frameon': False, 'savefig.dpi': 600,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.086614173228346, 3.267716535433071),
                            gridspec_kw={'left': .075, 'right': .98, 'bottom': .23, 'top': .71, 'wspace': .42})
    finite_cost = finite['same_state_N12']
    waterfall(axes[0], [finite_cost['geometry_distribution_cost'],
              finite_cost['genuine_conditional_orbital_excitation_cost'], finite_cost['full_state_K']],
              ['Pattern\ncost', 'Internal\ncost', 'Transfer'],
              ['#B6BEC2', '#B89B6D', '#438891'], (-1.63, 1.61))
    actual_cost = macro['mobile_minus_static']
    waterfall(axes[1], [actual_cost['static_debt'], actual_cost['movement']],
              ['Static\nresponse', 'Transfer'], ['#B89B6D', '#438891'], (-6.81, 8.37))
    axes[0].set_yticks([-1, 0, 1])
    axes[1].set_yticks([-5, 0, 5])
    ax = axes[2]
    for y, key, color, label in [(1, 'static_endpoint_candidate_window', '#AAB5B9', 'Static'),
                                (0, 'mobile_endpoint_candidate_window', '#438891', 'Mobile')]:
        lo, hi = macro[key]
        ax.hlines(y, lo, hi, color=color, lw=5.8, zorder=3)
        ax.plot([lo, hi], [y, y], marker='|', ls='none', ms=10, mew=1, color=color, zorder=4)
    ax.axvline(macro['mu'], ymax=.80, color='#C36E4E', lw=.9, linestyle=(0,(3,2)), zorder=1)
    ax.set_xlim(-18.225, -17.915)
    ax.set_ylim(-.65, 1.65)
    ax.set_yticks([0, 1], ['Mobile', 'Static'])
    ax.set_xticks([-18.2, -18.1, -18.0, -17.9])
    ax.set_xlabel('Chemical potential µ', labelpad=7)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.text(-18.1, 1.40, '−18.10', color='#A95138', ha='center', va='center', fontsize=7)
    titles = ['Complete conditional spectrum', 'Matched orbital control', 'Endpoint-candidate intervals']
    subtitles = ['16 sites: all occupation patterns', '96 sites: both states χ = 512', '96 sites: fixed measured states']
    for label, ax, title, subtitle in zip('abc', axes, titles, subtitles):
        ax.set_label(label)
        ax.text(0, 1.13, label, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')
        ax.text(.12, 1.13, title, transform=ax.transAxes, fontsize=7.3, va='bottom')
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=6.6, va='bottom', color='#59666C')
    fig.text(.075, .91, 'Internal adaptation pays for collective motion', fontsize=10, weight='bold', color='#234953')
    fig.text(.075, .835, 'Common physical parameters: λ = 0.70, µ = −18.10', fontsize=8, color='#59666C')
    fig.text(.075, .075, 'Deterministic energies; candidate intervals in c compare only vacuum, full activity and each displayed state.',
             fontsize=6.4, color='#59666C')
    fig.canvas.draw()
    sys.path.insert(0, str(a.audit_dir))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    require_matplotlib_panel_alignment(fig, json_out=HERE/'figure.alignment.json',
            overlay_svg=HERE/'figure.alignment.svg', tolerance_pt=1.5,
            gutter_tolerance_pt=1.5, strict=True, require_panel_labels=True)
    fig.savefig(HERE/'figure.pdf')
    fig.savefig(HERE/'figure.svg')
    fig.savefig(HERE/'figure.png', dpi=600)
    plt.close(fig)
    rows = [('a', 'Pattern response', finite_cost['geometry_distribution_cost']),
            ('a', 'Internal excitation', finite_cost['genuine_conditional_orbital_excitation_cost']),
            ('a', 'Transfer', finite_cost['full_state_K']),
            ('a', 'Net', finite_cost['net_energy_difference_full_minus_P']),
            ('b', 'Static response', actual_cost['static_debt']),
            ('b', 'Transfer', actual_cost['movement']), ('b', 'Net', actual_cost['total'])]
    with (HERE/'source_data.csv').open('w', newline='') as handle:
        writer=csv.writer(handle); writer.writerow(['panel','quantity','energy'])
        writer.writerows(rows)
    (HERE/'interval_source_data.json').write_text(json.dumps({key:macro[key] for key in
        ('lambda_transfer','mu','static_endpoint_candidate_window','mobile_endpoint_candidate_window')},indent=2)+'\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--audit-dir', type=Path, required=True)
    main(parser.parse_args())
