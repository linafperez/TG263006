#!/usr/bin/env python3

import argparse
import csv


COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def ranking_key(row):
    return (
        float(row[10]),    # menor e-value
        -float(row[11]),  # mayor bitscore
        -float(row[2]),   # mayor identidad
        -float(row[3]),   # mayor longitud
    )


def main():
    parser = argparse.ArgumentParser(
        description="Keep the best DIAMOND hit per query."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    best_hits = {}

    with open(args.input, newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue

            if len(row) != len(COLUMNS):
                raise ValueError(
                    f"Expected {len(COLUMNS)} columns at line {line_number}, "
                    f"found {len(row)}"
                )

            qseqid = row[0]
            key = ranking_key(row)

            if qseqid not in best_hits or key < best_hits[qseqid][0]:
                best_hits[qseqid] = (key, row)

    with open(args.output, "w", newline="") as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(COLUMNS)

        for qseqid in sorted(best_hits):
            writer.writerow(best_hits[qseqid][1])


if __name__ == "__main__":
    main()
