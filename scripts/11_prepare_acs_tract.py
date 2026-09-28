#!/usr/bin/env python3
"""11_prepare_acs_tract.py --------------------------------------------------

Build a census-tract dataset from the ACS 5-year table-based summary files.

Why this dataset exists.  The BRFSS and quiron case studies both turned out to
contain almost no domain shift: the between-environment variance in the
covariate means is under 1% of the within-environment variance in each
(`scripts/10_environment_shift.R`).  The correction term
`K' Sigma_0^{-1} (X_0 - mu_0)` can only do work when the target domain's
covariate distribution departs from the training pool's, so in that regime
generative invariance is *expected* to coincide with pooled OLS, and neither
dataset can demonstrate an advantage.

Measuring several candidates showed the reason is structural rather than a
matter of which variables are used.  An individual inside a US state sits at a
ratio of ~0.007 whether the covariates are lifestyle behaviours (BRFSS 0.0074)
or human capital (ACS person-level 0.0073).  Individual heterogeneity simply
swamps any geographic mean.  Refining the geography helps a little (ACS at PUMA
level, 0.0298); changing the *unit* from a person to a small area helps by an
order of magnitude (UCI Communities and Crime, communities inside states,
0.4418).

This script applies that finding to clean, current, fully public data: the unit
is a census tract, the environment is the state, and every variable is a
documented ACS estimate rather than a pre-normalised extract.  Compared with
Communities and Crime it has ~85,000 units instead of 1,994 and a median of
about 1,300 tracts per state instead of 26, which is what the per-environment
plug-in covariances need.

It also removes a weakness the reviewers may press on in both existing case
studies: `X_ei ~ N(mu_e, Sigma_e)` applied to eleven binary dummies is a
stretch, and tract-level shares, medians and indices are genuinely continuous.

Inputs, from
  https://www2.census.gov/programs-surveys/acs/summary_file/2022/table-based-SF/data/5YRData/
    acsdt5y2022-b19013.dat   median household income   (response)
    acsdt5y2022-b15003.dat   educational attainment
    acsdt5y2022-b23025.dat   employment status
    acsdt5y2022-b25077.dat   median house value
    acsdt5y2022-b25064.dat   median gross rent
    acsdt5y2022-b25003.dat   tenure
    acsdt5y2022-b25010.dat   average household size
    acsdt5y2022-b19083.dat   Gini index
    acsdt5y2022-b08303.dat   travel time to work
    acsdt5y2022-b01002.dat   median age            (confounder)
    acsdt5y2022-b02001.dat   race                  (confounder)
    acsdt5y2022-b05002.dat   nativity              (confounder)

Usage:
    python3 scripts/11_prepare_acs_tract.py --in=DIR --out=PATH
"""

import argparse
import csv
import os
import sys

# Tract-level rows in the summary files carry this geography prefix; the eleven
# digits after "US" are state(2) + county(3) + tract(6).
TRACT_PREFIX = "1400000US"

FIPS_STATE = {
    "01": "Alabama", "02": "Alaska", "04": "Arizona", "05": "Arkansas",
    "06": "California", "08": "Colorado", "09": "Connecticut",
    "10": "Delaware", "11": "District of Columbia", "12": "Florida",
    "13": "Georgia", "15": "Hawaii", "16": "Idaho", "17": "Illinois",
    "18": "Indiana", "19": "Iowa", "20": "Kansas", "21": "Kentucky",
    "22": "Louisiana", "23": "Maine", "24": "Maryland",
    "25": "Massachusetts", "26": "Michigan", "27": "Minnesota",
    "28": "Mississippi", "29": "Missouri", "30": "Montana",
    "31": "Nebraska", "32": "Nevada", "33": "New Hampshire",
    "34": "New Jersey", "35": "New Mexico", "36": "New York",
    "37": "North Carolina", "38": "North Dakota", "39": "Ohio",
    "40": "Oklahoma", "41": "Oregon", "42": "Pennsylvania",
    "44": "Rhode Island", "45": "South Carolina", "46": "South Dakota",
    "47": "Tennessee", "48": "Texas", "49": "Utah", "50": "Vermont",
    "51": "Virginia", "53": "Washington", "54": "West Virginia",
    "55": "Wisconsin", "56": "Wyoming", "72": "Puerto Rico",
}


def read_table(path, cells):
    """Read one summary file, returning {geo_id: {cell: float or None}}.

    Suppressed and unavailable estimates are written as negative sentinels or
    as blanks; both become None so that they are dropped rather than silently
    treated as data.
    """
    out = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="|")
        header = next(reader)
        try:
            idx = {c: header.index(c) for c in cells}
        except ValueError as err:
            raise SystemExit("%s: %s (have %s)"
                             % (os.path.basename(path), err, header[:6]))
        for row in reader:
            if not row or not row[0].startswith(TRACT_PREFIX):
                continue
            vals = {}
            for c, i in idx.items():
                raw = row[i].strip() if i < len(row) else ""
                if raw in ("", ".", "null"):
                    vals[c] = None
                else:
                    try:
                        v = float(raw)
                    except ValueError:
                        vals[c] = None
                        continue
                    # ACS marks suppression with large negative values.
                    vals[c] = None if v <= -100000 else v
            out[row[0]] = vals
    return out


def ratio(num, den):
    if num is None or den is None or den <= 0:
        return None
    return num / den


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_dir", required=True,
                    help="directory holding the .dat files")
    ap.add_argument("--out", required=True, help="output CSV")
    args = ap.parse_args(argv)

    d = args.in_dir

    def load(tbl, cells):
        path = os.path.join(d, tbl + ".dat")
        if not os.path.exists(path):
            raise SystemExit("Missing %s" % path)
        sys.stderr.write("reading %s\n" % os.path.basename(path))
        return read_table(path, cells)

    inc = load("b19013", ["B19013_E001"])
    edu = load("b15003", ["B15003_E001", "B15003_E022", "B15003_E023",
                          "B15003_E024", "B15003_E025"])
    emp = load("b23025", ["B23025_E001", "B23025_E002", "B23025_E003",
                          "B23025_E005"])
    val = load("b25077", ["B25077_E001"])
    rnt = load("b25064", ["B25064_E001"])
    ten = load("b25003", ["B25003_E001", "B25003_E002"])
    hhs = load("b25010", ["B25010_E001"])
    gin = load("b19083", ["B19083_E001"])
    cmt = load("b08303", ["B08303_E001"] +
               ["B08303_E%03d" % i for i in range(8, 14)])
    age = load("b01002", ["B01002_E001"])
    rac = load("b02001", ["B02001_E001", "B02001_E002", "B02001_E003"])
    nat = load("b05002", ["B05002_E001", "B05002_E013"])

    cols = [
        "state", "geo_id", "median_income",
        # covariates: continuous tract characteristics
        "pct_bachelors_plus", "unemployment_rate", "labor_force_rate",
        "median_house_value", "median_rent", "pct_owner_occupied",
        "avg_household_size", "gini", "pct_commute_30plus",
        # confounders: present in the file, excluded from the design on purpose
        "median_age", "pct_white", "pct_black", "pct_foreign_born",
    ]

    n_written = 0
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for geo in inc:
            fips = geo[len(TRACT_PREFIX):len(TRACT_PREFIX) + 2]
            state = FIPS_STATE.get(fips)
            if state is None:
                continue
            e = edu.get(geo, {})
            p = emp.get(geo, {})
            t = ten.get(geo, {})
            c = cmt.get(geo, {})
            r = rac.get(geo, {})
            nb = nat.get(geo, {})

            ba = None
            if e.get("B15003_E001"):
                parts = [e.get("B15003_E0%d" % k) for k in (22, 23, 24, 25)]
                if all(v is not None for v in parts):
                    ba = sum(parts) / e["B15003_E001"]

            c30 = None
            if c.get("B08303_E001"):
                parts = [c.get("B08303_E%03d" % i) for i in range(8, 14)]
                if all(v is not None for v in parts):
                    c30 = sum(parts) / c["B08303_E001"]

            row = [
                state, geo, inc[geo].get("B19013_E001"),
                ba,
                ratio(p.get("B23025_E005"), p.get("B23025_E003")),
                ratio(p.get("B23025_E002"), p.get("B23025_E001")),
                val.get(geo, {}).get("B25077_E001"),
                rnt.get(geo, {}).get("B25064_E001"),
                ratio(t.get("B25003_E002"), t.get("B25003_E001")),
                hhs.get(geo, {}).get("B25010_E001"),
                gin.get(geo, {}).get("B19083_E001"),
                c30,
                age.get(geo, {}).get("B01002_E001"),
                ratio(r.get("B02001_E002"), r.get("B02001_E001")),
                ratio(r.get("B02001_E003"), r.get("B02001_E001")),
                ratio(nb.get("B05002_E013"), nb.get("B05002_E001")),
            ]
            w.writerow(["" if v is None else v for v in row])
            n_written += 1

    sys.stderr.write("wrote %d tracts to %s\n" % (n_written, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
