#!/usr/bin/env python3

import argparse


def read_fasta(path):
    records = {}
    order = []

    current_id = None
    current_lines = []

    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current_id is not None:
                    if current_id in records:
                        raise ValueError(
                            f"Duplicate FASTA identifier: {current_id}"
                        )

                    records[current_id] = current_lines
                    order.append(current_id)

                current_id = line[1:].split()[0]
                current_lines = [line]
            else:
                if current_id is not None:
                    current_lines.append(line)

    if current_id is not None:
        if current_id in records:
            raise ValueError(
                f"Duplicate FASTA identifier: {current_id}"
            )

        records[current_id] = current_lines
        order.append(current_id)

    return records, order


def write_selected_fasta(records, selected_ids, output_path):
    with open(output_path, "w") as out:
        for record_id in selected_ids:
            if record_id not in records:
                raise ValueError(
                    f"Selected identifier missing from FASTA: {record_id}"
                )

            out.write(
                "\n".join(records[record_id]) + "\n"
            )


def get_complete_ids_from_pep(pep_path):
    complete_ids = []

    with open(pep_path, "r") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            record_id = header.split()[0]

            if "type:complete" in header.lower():
                complete_ids.append(record_id)

    if len(complete_ids) != len(set(complete_ids)):
        raise ValueError(
            "Duplicate complete ORF identifiers found in peptide FASTA"
        )

    return complete_ids


def parse_attributes(attribute_text):
    attributes = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue

        attributes[key.strip()] = value.strip().strip('"')

    return attributes


def split_attribute_ids(value):
    if not value:
        return []

    return [
        item.strip().strip('"')
        for item in value.split(",")
        if item.strip()
    ]


def filter_gff3(gff3_path, selected_ids, output_path):
    selected_ids = set(selected_ids)

    entries = []
    records = []
    id_to_parents = {}

    with open(gff3_path, "r") as inp:
        for line_number, line in enumerate(inp, start=1):
            raw_line = line.rstrip("\n")

            if not raw_line or raw_line.startswith("#"):
                entries.append(("raw", raw_line, None))
                continue

            parts = raw_line.split("\t")

            if len(parts) != 9:
                raise ValueError(
                    f"Invalid GFF3 row at line {line_number}: "
                    f"expected 9 columns, found {len(parts)}"
                )

            attributes = parse_attributes(parts[8])

            record_id = attributes.get("ID")
            parents = split_attribute_ids(
                attributes.get("Parent")
            )

            record = {
                "line": raw_line,
                "id": record_id,
                "parents": parents,
                "feature": parts[2],
            }

            records.append(record)
            entries.append(("record", raw_line, record))

            if record_id:
                if record_id in id_to_parents:
                    existing = id_to_parents[record_id]

                    if existing != parents:
                        raise ValueError(
                            f"Conflicting GFF3 definitions for ID: "
                            f"{record_id}"
                        )
                else:
                    id_to_parents[record_id] = parents

    #
    # Ancestors:
    # Preserve gene/parent records belonging to each selected ORF.
    #
    ancestor_ids = set()
    frontier = list(selected_ids)

    while frontier:
        current = frontier.pop()

        for parent in id_to_parents.get(current, []):
            if (
                parent not in selected_ids
                and parent not in ancestor_ids
            ):
                ancestor_ids.add(parent)
                frontier.append(parent)

    #
    # Descendants:
    # Preserve CDS/exon/etc. records whose Parent belongs exactly
    # to a selected ORF. Recursively preserve lower-level children.
    #
    descendant_ids = set(selected_ids)

    changed = True

    while changed:
        changed = False

        for record in records:
            record_id = record["id"]

            if not record_id:
                continue

            if record_id in descendant_ids:
                continue

            if any(
                parent in descendant_ids
                for parent in record["parents"]
            ):
                descendant_ids.add(record_id)
                changed = True

    kept_records = 0
    matched_selected_ids = set()

    with open(output_path, "w") as out:
        for entry_type, raw_line, record in entries:
            if entry_type == "raw":
                out.write(raw_line + "\n")
                continue

            record_id = record["id"]
            parents = record["parents"]

            keep = (
                record_id in selected_ids
                or record_id in ancestor_ids
                or record_id in descendant_ids
                or any(
                    parent in descendant_ids
                    for parent in parents
                )
            )

            if not keep:
                continue

            out.write(raw_line + "\n")
            kept_records += 1

            if record_id in selected_ids:
                matched_selected_ids.add(record_id)

    missing_ids = sorted(
        selected_ids - matched_selected_ids
    )

    if missing_ids:
        preview = ", ".join(missing_ids[:10])

        raise ValueError(
            f"{len(missing_ids)} complete ORFs were not found as "
            f"exact GFF3 IDs. Examples: {preview}"
        )

    return kept_records


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Retain complete TransDecoder ORFs and their exact "
            "GFF3 feature hierarchy."
        )
    )

    parser.add_argument("--pep", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--gff3", required=True)

    parser.add_argument("--output-pep", required=True)
    parser.add_argument("--output-cds", required=True)
    parser.add_argument("--output-gff3", required=True)

    args = parser.parse_args()

    pep_records, pep_order = read_fasta(args.pep)
    cds_records, _ = read_fasta(args.cds)

    complete_ids = get_complete_ids_from_pep(
        args.pep
    )

    if not complete_ids:
        raise ValueError(
            "No complete TransDecoder isoforms found."
        )

    complete_id_set = set(complete_ids)

    complete_ids = [
        record_id
        for record_id in pep_order
        if record_id in complete_id_set
    ]

    missing_cds = [
        record_id
        for record_id in complete_ids
        if record_id not in cds_records
    ]

    if missing_cds:
        preview = ", ".join(missing_cds[:10])

        raise ValueError(
            f"{len(missing_cds)} complete peptide IDs are missing "
            f"from CDS FASTA. Examples: {preview}"
        )

    write_selected_fasta(
        pep_records,
        complete_ids,
        args.output_pep
    )

    write_selected_fasta(
        cds_records,
        complete_ids,
        args.output_cds
    )

    kept_gff_records = filter_gff3(
        args.gff3,
        complete_ids,
        args.output_gff3
    )

    print(f"Input proteins: {len(pep_records)}")
    print(f"Complete proteins kept: {len(complete_ids)}")
    print(f"GFF3 records kept: {kept_gff_records}")
    print(f"Output PEP: {args.output_pep}")
    print(f"Output CDS: {args.output_cds}")
    print(f"Output GFF3: {args.output_gff3}")


if __name__ == "__main__":
    main()
