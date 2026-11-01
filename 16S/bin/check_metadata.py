#!/usr/bin/env python3

import argparse
import csv
import shutil
import sys
from pathlib import Path


VALID_ID_HEADERS_CASE_INSENSITIVE = {
    "id",
    "sampleid",
    "sample id",
    "sample-id",
}

VALID_ID_HEADERS_CASE_SENSITIVE = {
    "#SampleID",
    "#Sample ID",
    "sample_name",
}


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def main():

    parser = argparse.ArgumentParser(
        description="Validate QIIME2 sample metadata against pipeline samplesheet."
    )

    parser.add_argument("--metadata", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)

    args = parser.parse_args()

    metadata = Path(args.metadata)
    samplesheet = Path(args.samplesheet)
    output = Path(args.output)
    report = Path(args.report)

    # ------------------------------------------------------------
    # READ METADATA
    # ------------------------------------------------------------

    if not metadata.exists():
        fail(f"Metadata file does not exist: {metadata}")

    with metadata.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if not rows:
        fail("Metadata file is empty.")

    header = rows[0]

    if len(header) < 1:
        fail("Metadata header is empty.")

    id_header = header[0].strip()

    valid_header = (
        id_header.lower() in VALID_ID_HEADERS_CASE_INSENSITIVE
        or id_header in VALID_ID_HEADERS_CASE_SENSITIVE
    )

    if not valid_header:
        fail(
            f"Invalid QIIME2 sample ID header: '{id_header}'. "
            "Use for example '#SampleID' or 'sample-id'."
        )

    # Column names
    cleaned_columns = [x.strip() for x in header]

    if any(not x for x in cleaned_columns):
        fail("Metadata contains an empty column name.")

    lowered = [x.lower() for x in cleaned_columns]

    if len(lowered) != len(set(lowered)):
        fail("Metadata contains duplicate column names.")

    # ------------------------------------------------------------
    # METADATA SAMPLE IDS
    # ------------------------------------------------------------

    metadata_ids = []

    for row_number, row in enumerate(rows[1:], start=2):

        if not row:
            continue

        first = row[0].strip() if row else ""

        # QIIME2 directives and comments
        if first.startswith("#"):
            continue

        if len(row) != len(header):
            fail(
                f"Metadata row {row_number} has {len(row)} columns; "
                f"expected {len(header)}."
            )

        if not first:
            fail(f"Metadata row {row_number} has an empty sample ID.")

        metadata_ids.append(first)

    if not metadata_ids:
        fail("Metadata contains no sample IDs.")

    if len(metadata_ids) != len(set(metadata_ids)):
        duplicates = sorted({
            x for x in metadata_ids
            if metadata_ids.count(x) > 1
        })
        fail(
            "Metadata contains duplicate sample IDs: "
            + ", ".join(duplicates)
        )

    # ------------------------------------------------------------
    # READ PIPELINE SAMPLESHEET
    # ------------------------------------------------------------

    if not samplesheet.exists():
        fail(f"Samplesheet does not exist: {samplesheet}")

    with samplesheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames or "sample" not in reader.fieldnames:
            fail("Samplesheet must contain a 'sample' column.")

        sample_ids = []

        for row in reader:
            sample = (row.get("sample") or "").strip()

            if sample:
                sample_ids.append(sample)

    sample_ids = sorted(set(sample_ids))

    if not sample_ids:
        fail("No sample IDs were found in the samplesheet.")

    # ------------------------------------------------------------
    # CROSS-CHECK
    # ------------------------------------------------------------

    metadata_set = set(metadata_ids)
    sample_set = set(sample_ids)

    missing = sorted(sample_set - metadata_set)
    extra = sorted(metadata_set - sample_set)

    if missing:
        fail(
            "Samples present in samplesheet but missing from metadata: "
            + ", ".join(missing)
        )

    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------

    shutil.copyfile(metadata, output)

    with report.open("w", encoding="utf-8") as handle:

        handle.write("QIIME2 METADATA VALIDATION\n")
        handle.write("=" * 60 + "\n")
        handle.write(f"Metadata file: {metadata.name}\n")
        handle.write(f"Samplesheet file: {samplesheet.name}\n")
        handle.write(f"ID column: {id_header}\n")
        handle.write(f"Metadata columns: {len(header)}\n")
        handle.write(f"Metadata sample IDs: {len(metadata_ids)}\n")
        handle.write(f"Pipeline sample IDs: {len(sample_ids)}\n")
        handle.write(f"Missing pipeline samples: {len(missing)}\n")
        handle.write(f"Extra metadata samples: {len(extra)}\n")

        if extra:
            handle.write(
                "Extra metadata IDs: " + ", ".join(extra) + "\n"
            )

        handle.write("STATUS: PASS\n")

    print("Metadata validation: PASS")
    print(f"Pipeline samples: {len(sample_ids)}")
    print(f"Metadata samples: {len(metadata_ids)}")

    if extra:
        print(
            "WARNING: metadata contains additional sample IDs: "
            + ", ".join(extra)
        )


if __name__ == "__main__":
    main()
