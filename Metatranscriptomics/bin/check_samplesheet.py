#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
import sys


REQUIRED_COLUMNS = [
    "sample_id",
    "sra_id",
    "condition",
    "fastq_1",
    "fastq_2",
]

ALLOWED_CONDITIONS = {
    "PD",
    "HC",
}


def main():
    parser = argparse.ArgumentParser(
        description="Validate metatranscriptomics samplesheet."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input samplesheet CSV."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Validated output samplesheet CSV."
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {input_path}")

    seen_samples = set()
    rows = []

    with input_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError("Samplesheet is empty or has no header.")

        missing_columns = [
            column for column in REQUIRED_COLUMNS
            if column not in reader.fieldnames
        ]

        if missing_columns:
            raise ValueError(
                "Missing required columns: "
                + ", ".join(missing_columns)
            )

        for line_number, row in enumerate(reader, start=2):
            sample_id = row["sample_id"].strip()
            sra_id = row["sra_id"].strip()
            condition = row["condition"].strip()
            fastq_1 = row["fastq_1"].strip()
            fastq_2 = row["fastq_2"].strip()

            if not sample_id:
                raise ValueError(f"Line {line_number}: sample_id is empty")

            if sample_id in seen_samples:
                raise ValueError(f"Line {line_number}: duplicate sample_id: {sample_id}")

            seen_samples.add(sample_id)

            if not condition:
                raise ValueError(f"Line {line_number}: condition is empty for {sample_id}")

            if condition not in ALLOWED_CONDITIONS:
                raise ValueError(
                    f"Line {line_number}: invalid condition for "
                    f"{sample_id}: {condition}. "
                    f"Expected one of: "
                    f"{', '.join(sorted(ALLOWED_CONDITIONS))}"
                )

            if not fastq_1:
                raise ValueError(f"Line {line_number}: fastq_1 is empty for {sample_id}")

            if not fastq_2:
                raise ValueError(f"Line {line_number}: fastq_2 is empty for {sample_id}")

            if not Path(fastq_1).exists():
                raise FileNotFoundError(
                    f"Line {line_number}: fastq_1 not found for {sample_id}: {fastq_1}"
                )

            if not Path(fastq_2).exists():
                raise FileNotFoundError(
                    f"Line {line_number}: fastq_2 not found for {sample_id}: {fastq_2}"
                )

            rows.append({
                "sample_id": sample_id,
                "sra_id": sra_id,
                "condition": condition,
                "fastq_1": fastq_1,
                "fastq_2": fastq_2,
            })

    if not rows:
        raise ValueError("Samplesheet has no samples.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Validated samplesheet: {input_path}")
    print(f"Samples: {len(rows)}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
