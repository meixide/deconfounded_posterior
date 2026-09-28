#!/usr/bin/env python3
"""06_prepare_brfss.py -----------------------------------------------------

Extract the case-study columns from a BRFSS annual LLCP file and write a CSV.

The BRFSS public-use files are distributed as SAS XPORT (version 5) transport
files.  This reads them with the standard library alone: the cluster's R
installation has neither `foreign` nor `haven`, and adding a dependency to the
replication package for one file-format step is a poor trade.  XPORT v5 is a
short, fully documented fixed-record format (SAS Technical Support document
TS-140), so a reader is about a hundred lines.

Only the columns named in `--columns` are decoded, which is what keeps this
fast: the 2023 file is 1.2 GB and 345 variables wide, and the case study needs
roughly twenty of them.

Usage:
    python3 scripts/06_prepare_brfss.py --xpt=PATH --out=PATH [--columns=A,B]

Download the input from
    https://www.cdc.gov/brfss/annual_data/annual_2023.html
    -> files/LLCP2023XPT.zip
"""

import argparse
import csv
import os
import struct
import sys

# ---------------------------------------------------------------------------
# XPORT v5 reader
# ---------------------------------------------------------------------------

_HEADER = b"HEADER RECORD*******"
_NAMESTR_LEN = 140

# A numeric missing value is written as a single missing code in the first byte
# followed by zeros: "." for ordinary missing, "_" and "A".."Z" for the special
# missing values.  No finite IBM float shares that representation, because a
# legitimate value with those exponent bits would need a non-zero mantissa.
_MISSING_CODES = frozenset(b"._ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# value = mantissa * 16 ** (exponent - 64) / 16 ** 14, precomputed per exponent.
_POW16 = tuple(16.0 ** (e - 78) for e in range(128))


def ibm_to_double(raw):
    """Convert an 8-byte big-endian IBM 360 hexadecimal float to a double.

    The format is one sign bit, seven exponent bits in excess-64 to base 16,
    and a 56-bit fraction interpreted as 14 hexadecimal digits after the point.
    """
    i = struct.unpack(">Q", raw)[0]
    if i == 0:
        return 0.0
    mantissa = i & 0x00FFFFFFFFFFFFFF
    value = mantissa * _POW16[(i >> 56) & 0x7F]
    return -value if i >> 63 else value


def _parse_namestrs(blob, n_vars):
    """Decode the NAMESTR block into a list of variable descriptors."""
    variables = []
    for k in range(n_vars):
        rec = blob[k * _NAMESTR_LEN:(k + 1) * _NAMESTR_LEN]
        ntype, _hash, length, number = struct.unpack(">hhhh", rec[:8])
        variables.append({
            "name": rec[8:16].decode("ascii", "replace").strip(),
            "label": rec[16:56].decode("ascii", "replace").strip(),
            "numeric": ntype == 1,
            "length": length,
            "number": number,
        })
    return variables


def read_xport_header(handle):
    """Read the library and member headers, returning the variable list.

    Leaves `handle` positioned at the first observation record.
    """
    first = handle.read(80)
    if not first.startswith(_HEADER + b"LIBRARY HEADER RECORD"):
        raise ValueError("Not a SAS XPORT v5 file (bad library header).")
    handle.read(160)                                    # two real header records

    if not handle.read(80).startswith(_HEADER + b"MEMBER"):
        raise ValueError("Expected a MEMBER header record.")
    if not handle.read(80).startswith(_HEADER + b"DSCRPTR"):
        raise ValueError("Expected a DSCRPTR header record.")
    handle.read(160)                                    # member name, labels

    namestr_header = handle.read(80)
    if not namestr_header.startswith(_HEADER + b"NAMESTR"):
        raise ValueError("Expected a NAMESTR header record.")
    n_vars = int(namestr_header[54:58])

    # NAMESTR records are packed contiguously and the block is then padded with
    # blanks to a multiple of 80 bytes.
    size = n_vars * _NAMESTR_LEN
    padded = size + (-size % 80)
    variables = _parse_namestrs(handle.read(padded), n_vars)

    if not handle.read(80).startswith(_HEADER + b"OBS"):
        raise ValueError("Expected an OBS header record.")
    return variables


def column_offsets(variables):
    """Byte offset of each variable within an observation record."""
    offsets, position = {}, 0
    for var in variables:
        offsets[var["name"]] = position
        position += var["length"]
    return offsets, position


def resolve_path(path):
    """Return `path`, tolerating the trailing space CDC ships in the archives."""
    if os.path.exists(path):
        return path
    if os.path.exists(path + " "):
        return path + " "
    raise IOError("No such file: %s" % path)


def read_columns(path, wanted, chunk_records=4096):
    """Stream a XPORT file, yielding one list of values per observation.

    `wanted` is a list of variable names; the yielded values follow that order.
    Character columns come back as stripped strings, numeric ones as floats or
    None for missing.
    """
    with open(resolve_path(path), "rb") as handle:
        variables = read_xport_header(handle)
        by_name = {v["name"]: v for v in variables}
        missing = [w for w in wanted if w not in by_name]
        if missing:
            raise KeyError("Not in file: %s\nAvailable: %s"
                           % (", ".join(missing),
                              ", ".join(sorted(by_name))))

        offsets, record_len = column_offsets(variables)
        # (offset, length, numeric) per requested column, in output order.
        plan = [(offsets[w], by_name[w]["length"], by_name[w]["numeric"])
                for w in wanted]

        block = chunk_records * record_len
        tail = b""
        while True:
            chunk = handle.read(block)
            if not chunk:
                # Whatever is left over is the blank padding that rounds the
                # observation block up to a multiple of 80 bytes.
                break
            buf = tail + chunk
            n_whole = len(buf) // record_len
            tail = buf[n_whole * record_len:]
            for r in range(n_whole):
                base = r * record_len
                rec = buf[base:base + record_len]
                row = []
                for off, length, numeric in plan:
                    raw = rec[off:off + length]
                    if numeric:
                        if raw[0] in _MISSING_CODES and not raw[1:].strip(b"\x00"):
                            row.append(None)
                        else:
                            row.append(ibm_to_double(raw.ljust(8, b"\x00")))
                    else:
                        row.append(raw.decode("latin-1").strip())
                yield row


# ---------------------------------------------------------------------------
# Case-study columns
# ---------------------------------------------------------------------------

# Variable names drifted between years (`_DRNKWK2` in 2023 became `_DRNKWK3` in
# 2024, `_RFDRHV8` became `_RFDRHV9`, and so on), so the set is resolved against
# what the file actually contains rather than hard-coded per year.
CANDIDATES = [
    # identifiers and outcome
    "_STATE", "_BMI5", "_LLCPWT", "_PSU", "_STSTR",
    # smoking and other tobacco
    "_SMOKER3", "SMOKE100", "SMOKDAY2", "USENOW3",
    "_CURECI1", "_CURECI2", "_CURECI3", "ECIGNOW1", "ECIGNOW2", "ECIGNOW3",
    # alcohol
    "_DRNKWK1", "_DRNKWK2", "_DRNKWK3", "DRNKANY5", "DRNKANY6",
    "_RFBING5", "_RFBING6", "_RFDRHV7", "_RFDRHV8", "_RFDRHV9", "DROCDY4_",
    # physical activity
    "_PACAT1", "_PACAT2", "_PACAT3", "PA3MIN_", "PA1MIN_", "_PASTRNG",
    "_TOTINDA", "_PAINDX3", "_PA150R4", "EXERANY2",
    # sleep, if the year has it
    "SLEPTIM1",
    # confounders: present in the file and excluded from the design on purpose
    "SEXVAR", "_SEX", "_AGE80", "_AGEG5YR", "INCOME3", "_INCOMG1",
    "EDUCA", "_EDUCAG", "DIABETE4", "_RACE", "GENHLTH", "MARITAL", "EMPLOY1",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xpt", required=True, help="path to LLCP20XX.XPT")
    parser.add_argument("--out", help="path to the output CSV")
    parser.add_argument("--columns", default="",
                        help="comma-separated column names; default is the "
                             "case-study set present in the file")
    parser.add_argument("--list", action="store_true",
                        help="print the variables in the file and exit")
    args = parser.parse_args(argv)

    try:
        path = resolve_path(args.xpt)
    except IOError as err:
        parser.error(str(err))

    with open(path, "rb") as handle:
        variables = read_xport_header(handle)
        header_bytes = handle.tell()
    available = [v["name"] for v in variables]
    # The observation block is a whole number of fixed-length records plus at
    # most 79 bytes of blank padding, so the row count is known before reading
    # any data.  Checking it against what the loop produces catches a silently
    # truncated download, which is otherwise very hard to notice.
    _, record_len = column_offsets(variables)
    expected_rows = (os.path.getsize(path) - header_bytes) // record_len

    if args.list:
        for name in available:
            sys.stdout.write(name + "\n")
        return 0
    if not args.out:
        parser.error("--out is required unless --list is given")

    if args.columns:
        wanted = [c.strip() for c in args.columns.split(",") if c.strip()]
        absent = [c for c in wanted if c not in available]
        if absent:
            parser.error("Not in file: %s" % ", ".join(absent))
    else:
        present = set(available)
        wanted = [c for c in CANDIDATES if c in present]

    sys.stderr.write("Reading %d of %d variables from %s\n"
                     % (len(wanted), len(available), os.path.basename(path)))
    sys.stderr.write("  %s\n" % ", ".join(wanted))

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    n = 0
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(wanted)
        for row in read_columns(path, wanted):
            writer.writerow(["" if v is None else v for v in row])
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("  %d rows\n" % n)

    sys.stderr.write("Wrote %d rows to %s\n" % (n, args.out))
    if n != expected_rows:
        sys.stderr.write("ERROR: expected %d records from the file geometry.\n"
                         % expected_rows)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
