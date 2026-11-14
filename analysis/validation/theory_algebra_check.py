"""Independent algebra checks for the EIF and exact nuisance drift.

Run from any directory with Python >= 3.9 (standard library only):
    python theory_algebra_check.py
The output directory receives theory_algebra_results.json.
No manuscript or scientific output is modified.

References: sections/inference.tex (main:tceeif, main:gapeif);
supplement-fragments/appendices.tex (eq:exactnuisancedrift).
The derivative test perturbs probability laws and compares centered finite
changes in the target with inner products against the stated EIF. It does
not obtain target derivatives by differentiating the EIF implementation.
"""

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
import random

X_PROB = [0.4, 0.6]
SOURCE = [0.3, 0.7]
ALLOCATION = [0.25, 0.8]
OUTCOMES = [-2.0, 1.0, 4.0, 7.0]
# Group order: trial treated, internal control, external control.
OUTCOME_PROB = [
    [[0.1, 0.2, 0.4, 0.3], [0.3, 0.4, 0.2, 0.1], [0.4, 0.1, 0.2, 0.3]],
    [[0.4, 0.3, 0.2, 0.1], [0.15, 0.25, 0.45, 0.15], [0.3, 0.4, 0.2, 0.1]],
]


def share(p, e, vi, v0):
    return p * (1.0 - e) / (p * (1.0 - e) + (1.0 - p) * vi / v0)


def moments(outcome_prob):
    means = [
        [sum(pr * y for pr, y in zip(probs, OUTCOMES)) for probs in groups]
        for groups in outcome_prob
    ]
    variances = [
        [
            sum(pr * (y - means[x][g]) ** 2 for pr, y in zip(probs, OUTCOMES))
            for g, probs in enumerate(groups)
        ]
        for x, groups in enumerate(outcome_prob)
    ]
    return means, variances


def numerator(x, g, y, p, means, variances):
    """H_eta evaluated at one record, with e fixed at its known value."""
    e = ALLOCATION[x]
    m1, mi, m0 = means[x]
    vi, v0 = variances[x][1:]
    p = p[x]
    lam = share(p, e, vi, v0)
    s, i1, ii, i0 = int(g < 2), int(g == 0), int(g == 1), int(g == 2)
    residual = y - means[x][g]
    contrast = m1 - lam * mi - (1.0 - lam) * m0
    core = (
        s * contrast
        + i1 / e * residual
        - ii * lam / (1.0 - e) * residual
        - i0 * p * (1.0 - lam) / (1.0 - p) * residual
    )
    correction = (
        (mi - m0)
        * lam
        * (1.0 - lam)
        * (
            -(s - p) / (1.0 - p)
            + ii / (1.0 - e) * (residual**2 / vi - 1.0)
            - i0 * p / (1.0 - p) * (residual**2 / v0 - 1.0)
        )
    )
    return core + correction


def law(f, p, outcome_prob):
    q = sum(fx * px for fx, px in zip(f, p))
    means, variances = moments(outcome_prob)
    lam = [
        share(p[x], ALLOCATION[x], variances[x][1], variances[x][2]) for x in range(2)
    ]
    contrast = [
        means[x][0] - lam[x] * means[x][1] - (1.0 - lam[x]) * means[x][2]
        for x in range(2)
    ]
    tce = sum(f[x] * p[x] * contrast[x] for x in range(2)) / q
    trial = sum(f[x] * p[x] * (means[x][0] - means[x][1]) for x in range(2)) / q
    records = []
    for x in range(2):
        e = ALLOCATION[x]
        group_prob = [p[x] * e, p[x] * (1.0 - e), 1.0 - p[x]]
        for g in range(3):
            for j, y in enumerate(OUTCOMES):
                s = int(g < 2)
                pr = f[x] * group_prob[g] * outcome_prob[x][g][j]
                phi = (numerator(x, g, y, p, means, variances) - s * tce) / q
                rct = (
                    s * (means[x][0] - means[x][1] - trial)
                    + int(g == 0) / e * (y - means[x][0])
                    - int(g == 1) / (1.0 - e) * (y - means[x][1])
                ) / q
                records.append(
                    dict(
                        x=x,
                        group=g,
                        y=y,
                        probability=pr,
                        tce_eif=phi,
                        trial_eif=rct,
                        gap_eif=phi - rct,
                    )
                )
    return dict(
        q=q,
        tce=tce,
        trial=trial,
        gap=tce - trial,
        means=means,
        variances=variances,
        shares=lam,
        records=records,
    )


def exponential_tilt(probabilities, direction, amount):
    weights = [pr * math.exp(amount * h) for pr, h in zip(probabilities, direction)]
    total = sum(weights)
    return [w / total for w in weights]


def perturbed_law(component, amount):
    f, p, outcome_prob = list(X_PROB), list(SOURCE), deepcopy(OUTCOME_PROB)
    if component == -2:
        f = exponential_tilt(f, [0.5, -0.2], amount)
    elif component == -1:
        p = [
            1.0 / (1.0 + math.exp(-math.log(px / (1.0 - px)) - amount * h))
            for px, h in zip(p, [0.4, -0.3])
        ]
    else:
        x, g = divmod(component, 3)
        outcome_prob[x][g] = exponential_tilt(
            outcome_prob[x][g], [0.4, -0.6, 0.2, 0.8], amount
        )
    return law(f, p, outcome_prob)


def exact_drift_formula(base, fitted_p, fitted_means, fitted_variances):
    total = 0.0
    for x in range(2):
        p, pb, e = SOURCE[x], fitted_p[x], ALLOCATION[x]
        mi, m0 = base["means"][x][1:]
        mib, m0b = fitted_means[x][1:]
        vi, v0 = base["variances"][x][1:]
        vib, v0b = fitted_variances[x][1:]
        lam, lb = base["shares"][x], share(pb, e, vib, v0b)
        dp, dmi, dm0 = pb - p, mib - mi, m0b - m0
        dvi, dv0 = vib - vi, v0b - v0
        wb = (1.0 - p) * pb / (1.0 - pb)
        d = (
            -p * (mi - m0) * (lb - lam)
            + (1.0 - lb) * dp * dm0 / (1.0 - pb)
            + (mib - m0b)
            * lb
            * (1.0 - lb)
            * (
                dp / (1.0 - pb)
                - p * dvi / vib
                + wb * dv0 / v0b
                + p * dmi**2 / vib
                - wb * dm0**2 / v0b
            )
        )
        total += X_PROB[x] * d
    return total


def run():
    base = law(X_PROB, SOURCE, OUTCOME_PROB)
    step = 1e-6
    derivative_checks = []
    for component in range(-2, 6):
        plus, minus = perturbed_law(component, step), perturbed_law(component, -step)
        label = (
            "covariate exponential tilt"
            if component == -2
            else (
                "source logistic tilt"
                if component == -1
                else "outcome tilt: stratum {}, group {}".format(*divmod(component, 3))
            )
        )
        for target in ["tce", "trial", "gap"]:
            derivative = (plus[target] - minus[target]) / (2.0 * step)
            represented = sum(
                record[target + "_eif"]
                * (rp["probability"] - rm["probability"])
                / (2.0 * step)
                for record, rp, rm in zip(
                    base["records"], plus["records"], minus["records"]
                )
            )
            derivative_checks.append(
                dict(
                    path=label,
                    target=target,
                    finite_difference=derivative,
                    eif_inner_product=represented,
                    absolute_error=abs(derivative - represented),
                )
            )
    rng = random.Random(20261007)
    drift_errors = []
    for _ in range(500):
        fitted_p = [rng.uniform(0.1, 0.9) for _ in range(2)]
        fitted_means = [
            [m + rng.uniform(-1.5, 1.5) for m in ms] for ms in base["means"]
        ]
        fitted_variances = [
            [v * math.exp(rng.uniform(-0.8, 0.8)) for v in vs]
            for vs in base["variances"]
        ]
        direct = (
            sum(
                r["probability"]
                * numerator(
                    r["x"], r["group"], r["y"], fitted_p, fitted_means, fitted_variances
                )
                for r in base["records"]
            )
            - base["q"] * base["tce"]
        )
        displayed = exact_drift_formula(base, fitted_p, fitted_means, fitted_variances)
        drift_errors.append(abs(direct - displayed))
    max_derivative_error = max(row["absolute_error"] for row in derivative_checks)
    max_drift_error = max(drift_errors)
    eif_means = {
        target: sum(r["probability"] * r[target + "_eif"] for r in base["records"])
        for target in ["tce", "trial", "gap"]
    }
    assert max_derivative_error < 1e-8, max_derivative_error
    assert max_drift_error < 1e-10, max_drift_error
    assert max(abs(value) for value in eif_means.values()) < 1e-12, eif_means
    result = dict(
        status="passed",
        purpose="Independent finite-law EIF and exact nuisance-drift identity checks",
        limits="Algebraic verification on a finite law; not a proof of general regularity or asymptotics.",
        inputs=dict(
            covariate_probabilities=X_PROB,
            source_probabilities=SOURCE,
            known_allocation=ALLOCATION,
            outcomes=OUTCOMES,
            outcome_probabilities=OUTCOME_PROB,
        ),
        baseline={k: v for k, v in base.items() if k != "records"},
        finite_difference_step=step,
        factor_paths=8,
        target_path_checks=len(derivative_checks),
        eif_means=eif_means,
        max_derivative_absolute_error=max_derivative_error,
        exact_drift_random_seed=20261007,
        exact_drift_tuples=len(drift_errors),
        max_exact_drift_absolute_error=max_drift_error,
        derivative_checks=derivative_checks,
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "results/theory-check",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    destination = args.output_dir / "theory_algebra_results.json"
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "status",
                    "factor_paths",
                    "target_path_checks",
                    "eif_means",
                    "max_derivative_absolute_error",
                    "exact_drift_tuples",
                    "max_exact_drift_absolute_error",
                ]
            },
            indent=2,
        )
    )
    print("Results saved to " + str(destination.resolve()))


if __name__ == "__main__":
    run()
