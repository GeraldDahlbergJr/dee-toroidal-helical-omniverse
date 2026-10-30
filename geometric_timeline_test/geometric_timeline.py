#!/usr/bin/env python3
"""Geometric instrumentation for a prescribed ODE, not a DEE field solver.

The index observes rates; it never edits a trajectory or drives transitions.
ClosedBoundary is an explicit solver-adapter contract, not a physical wall.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np
import test_transition as model

OUT = Path(__file__).resolve().parent


class ExternalInputRejected(RuntimeError):
    pass


class ClosedBoundary:
    """No external source/incoming-data channel is allowed in this baseline.

    A real PDE adapter must additionally implement and check its characteristic
    boundary data, source terms and moving-boundary flux. This ODE has none.
    Rejecting an input aborts the attempted step; it is never silently clipped.
    """
    def __init__(self):
        self.events = []

    def rhs(self, t, y, external=None):
        if external is not None:
            v = np.asarray(external, dtype=float)
            if v.shape != (3,) or not np.all(np.isfinite(v)) or np.any(v != 0):
                self.events.append({'time': float(t), 'event': 'external_input_rejected',
                                    'input': repr(external), 'step_accepted': False})
                raise ExternalInputRejected('Closed baseline: attempted external input')
        return model.rhs(t, y)


def rk_step(boundary, t, y, h):
    a = boundary.rhs(t, y)
    b = boundary.rhs(t+h/2, y+h*a/2)
    c = boundary.rhs(t+h/2, y+h*b/2)
    d = boundary.rhs(t+h, y+h*c)
    return y+h*(a+2*b+2*c+d)/6


def integrate(h=.0125):
    boundary = ClosedBoundary()
    n = round(model.T/h)
    t = np.linspace(0, model.T, n+1)
    y = np.empty((n+1, 3)); y[0] = [.2, .35, 0.]
    defects = np.zeros(n+1)
    for i in range(n):
        dt = t[i+1]-t[i]
        full = rk_step(boundary, t[i], y[i], dt)
        half = rk_step(boundary, t[i], y[i], dt/2)
        halves = rk_step(boundary, t[i]+dt/2, half, dt/2)
        defects[i+1] = np.linalg.norm(full-halves)
        y[i+1] = full
    return t, y, defects


def index(t, y, rho=None):
    """Offline index using measured rates, never model.gate(t).

    Rates come from neighboring samples, so events have timestep uncertainty.
    Input angles must be lifted, trusted histories. Only the existing
    translating toroidal-tube geometry is covered here.
    """
    t, y = np.asarray(t), np.asarray(y)
    if t.ndim != 1 or len(t) < 3 or y.shape != (len(t), 3):
        raise ValueError('Expected increasing time and N by 3 states')
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError('Time must be finite and strictly increasing')
    rho = np.full(len(t), .55) if rho is None else np.asarray(rho)
    if rho.shape != t.shape:
        raise ValueError('One radius per sample is required')
    # Reject first, so NaNs are never passed to differentiation/classification.
    for i in range(len(t)):
        if not np.all(np.isfinite(y[i])) or not np.isfinite(rho[i]):
            return [], [{'time': float(t[i]), 'row': i, 'event': 'nonfinite_state'}]
    rates = np.gradient(y, t, axis=0, edge_order=2)
    radial_rate = np.gradient(rho, t, edge_order=2)
    q = (1+rho/model.R*np.cos(y[:, 0]))*rates[:, 0]/model.WP
    axial = rates[:, 2]/model.V
    dq = np.gradient(q, t, edge_order=2)
    eps = 1e-4  # Declared classification tolerance; not a DEE acceptance gate.
    phase = 0
    labels = ['helical', 'transition_to_toroidal', 'toroidal',
              'restore_axial', 'axial_restored']
    rows, events = [], []
    start_phi = y[0, 1]
    toroidal_phi = None
    for i, time in enumerate(t):
        turns = (y[i, 1]-start_phi)/(2*math.pi)
        event = None
        if not model.LO <= rho[i] <= model.HI:
            event = 'outside_tube'
        elif abs(radial_rate[i]) > 1e-8:
            event = 'radial_flux_detected'
        elif rates[i, 1] < -eps:
            event = 'clockwise_reversal'
        elif rates[i, 1] <= eps:
            event = 'circulation_stall'
        elif abs(q[i]-axial[i]) > 2e-3:
            event = 'winding_axial_mismatch'
        elif q[i] < -2e-3 or q[i] > 1+2e-3:
            event = 'rate_outside_baseline_range'
        if event:
            events.append({'time': float(time), 'row': i, 'event': event,
                           'last_stage': labels[phase], 'turns': float(turns)})
            break  # No anomaly is erased or permitted to advance the stage.
        candidate = (0 if phase < 3 else 4) if q[i] >= 1-eps else (
            2 if q[i] <= eps else (1 if dq[i] < 0 else 3))
        if candidate != phase:
            if candidate != phase+1:
                events.append({'time': float(time), 'row': i,
                               'event': 'stage_sequence_violation',
                               'from': labels[phase], 'to': labels[candidate]})
                break
            if candidate == 1 and turns < 1-1e-6:
                events.append({'time': float(time), 'row': i,
                               'event': 'incomplete_helical_pass'})
                break
            if candidate == 3 and (y[i, 1]-toroidal_phi)/(2*math.pi) < 3-.1:
                events.append({'time': float(time), 'row': i,
                               'event': 'incomplete_toroidal_hold'})
                break
            phase = candidate
            if phase == 2:
                toroidal_phi = y[i, 1]
            events.append({'time': float(time), 'row': i, 'event': 'stage_entered',
                           'stage': labels[phase], 'turns': float(turns)})
        rows.append(dict(sample=i, time=float(time), stage_index=phase,
                         stage=labels[phase], phi_lifted=float(y[i, 1]),
                         ccw_turns=float(turns), theta=float(y[i, 0]),
                         axial_position=float(y[i, 2]), rho=float(rho[i]),
                         jacobian=float(rho[i]*(model.R+rho[i]*np.cos(y[i, 0]))),
                         azimuthal_rate=float(rates[i, 1]), winding_fraction=float(q[i]),
                         axial_fraction=float(axial[i]), radial_rate=float(radial_rate[i])))
    return rows, events


def main():
    t, y, defects = integrate()
    rows, events = index(t, y)
    checks = []
    def check(name, ok, **evidence):
        checks.append(dict(name=name, passed=bool(ok), **evidence))
    check('baseline completes five observed stages',
          len(rows) == len(t) and {r['stage_index'] for r in rows} == set(range(5)),
          stages=list(dict.fromkeys(r['stage'] for r in rows)))
    check('counterclockwise throughout', min(r['azimuthal_rate'] for r in rows) > 0)
    check('step refinement defect', max(defects) < 1e-8, max_defect=float(max(defects)))
    error = np.max(np.linalg.norm(model.embed(y)-model.embed(model.exact(t)), axis=1))
    check('independent trajectory reference', error < 1e-7, max_position_error=float(error))
    check('toroidal axial pause is not a circulation stall',
          any(r['stage'] == 'toroidal' for r in rows) and
          all(e['event'] == 'stage_entered' for e in events))
    boundary = ClosedBoundary()
    caught = False
    try:
        boundary.rhs(1., y[80], external=[1., 0., 0.])
    except ExternalInputRejected:
        caught = True
    check('external forcing rejected during helix', caught and len(boundary.events) == 1,
          events=boundary.events)
    # Controlled faults are known injections, not diagnosed physical causes.
    faults = {}
    stalled = y.copy(); stalled[:, 1] = y[0, 1]
    _, faults['circulation_stall'] = index(t, stalled)
    reversed_y = y.copy(); reversed_y[:, 1] *= -1
    _, faults['clockwise_reversal'] = index(t, reversed_y)
    radial = np.full(len(t), .55); radial[80:] = .80
    _, faults['radial_flux_detected'] = index(t, y, radial)
    invalid = y.copy(); invalid[80, 0] = np.nan
    _, faults['nonfinite_state'] = index(t, invalid)
    mismatch = y.copy(); mismatch[:, 2] = 0
    _, faults['winding_axial_mismatch'] = index(t, mismatch)
    for expected, observed in faults.items():
        check('detect injected '+expected, any(e['event'] == expected for e in observed),
              first_event=observed[0] if observed else None)
    with (OUT/'geometric_timeline.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[*rows[0], 'rk4_step_defect'])
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, 'rk4_step_defect': float(defects[row['sample']])})
    summary = dict(scope='Instrumented prescribed ODE; no DEE PDE or Run 43 data',
                   thresholds=dict(classification_fraction=1e-4,
                                   winding_axial_mismatch=2e-3, radial_rate=1e-8,
                                   minimum_helical_turns=1, toroidal_turns=3,
                                   toroidal_turn_detection_allowance=.1),
                   event_time_uncertainty='At least one timestep plus classification tolerance; offline derivatives use neighboring samples',
                   checks=checks, events=events, total_rows=len(rows),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'timeline_results.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(passed=sum(c['passed'] for c in checks), total=len(checks),
                          failures=[c for c in checks if not c['passed']], events=events), indent=2))
    return int(not all(c['passed'] for c in checks))


if __name__ == '__main__':
    raise SystemExit(main())
