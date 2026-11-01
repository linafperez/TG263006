#!/usr/bin/env python3

import argparse


def parse_attributes(attr_string):
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')

    return attrs


def read_diamond_hits(path):
    hits = {}

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            parts = line.split("\t")

            if parts[0] == "qseqid":
                continue

            if len(parts) < 12:
                raise ValueError(
                    f"DIAMOND row has fewer than 12 columns: {line}"
                )

            qseqid = parts[0]
            sseqid = parts[1]

            if qseqid in hits:
                raise ValueError(
                    f"Multiple best hits found for qseqid: {qseqid}"
                )

            hits[qseqid] = {
                "sseqid": sseqid,
                "gene_id": sseqid,
                "mag_id": sseqid.split("|", 1)[0],
                "pident": parts[2],
                "evalue": parts[10],
                "bitscore": parts[11],
            }

    return hits


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create transcript-to-MAG-gene mappings from "
            "TransDecoder GFF3 and DIAMOND best hits."
        )
    )

    parser.add_argument(
        "--diamond-best-hits",
        required=True,
    )

    parser.add_argument(
        "--transdecoder-gff3",
        required=True,
    )

    parser.add_argument(
        "--output-gtf",
        required=True,
    )

    parser.add_argument(
        "--output-map",
        required=True,
    )

    args = parser.parse_args()

    hits = read_diamond_hits(args.diamond_best_hits)

    if not hits:
        raise ValueError("No DIAMOND hits found.")

    mappings = {}
    gtf_records = []

    with open(args.transdecoder_gff3) as gff:
        for line in gff:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) != 9:
                continue

            (
                seqid,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attr_string,
            ) = parts

            # The TransDecoder mRNA ID is exactly the protein/ORF ID
            # used as qseqid by DIAMOND.
            if feature != "mRNA":
                continue

            attrs = parse_attributes(attr_string)
            qseqid = attrs.get("ID")

            if qseqid not in hits:
                continue

            hit = hits[qseqid]

            transcript_id = seqid
            gene_id = hit["gene_id"]
            mag_id = hit["mag_id"]

            if transcript_id in mappings:
                previous = mappings[transcript_id]

                if previous["gene_id"] != gene_id:
                    raise ValueError(
                        f"Transcript {transcript_id} maps to multiple MAG genes: "
                        f"{previous['gene_id']} and {gene_id}"
                    )

                continue

            mappings[transcript_id] = {
                "qseqid": qseqid,
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "mag_id": mag_id,
                "sseqid": hit["sseqid"],
                "pident": hit["pident"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
            }

            gtf_attrs = (
                f'gene_id "{gene_id}"; '
                f'transcript_id "{transcript_id}"; '
                f'orf_id "{qseqid}"; '
                f'mag_id "{mag_id}"; '
                f'diamond_target "{hit["sseqid"]}"; '
                f'pident "{hit["pident"]}"; '
                f'evalue "{hit["evalue"]}"; '
                f'bitscore "{hit["bitscore"]}";'
            )

            gtf_records.append(
                "\t".join([
                    seqid,
                    "create_transcript_mag_gtf",
                    "transcript",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    gtf_attrs,
                ])
            )

    matched_qseqids = {
        mapping["qseqid"]
        for mapping in mappings.values()
    }

    missing_qseqids = sorted(set(hits) - matched_qseqids)

    if missing_qseqids:
        preview = ", ".join(missing_qseqids[:10])

        raise ValueError(
            f"{len(missing_qseqids)} DIAMOND qseqids could not be matched "
            f"exactly to TransDecoder mRNA IDs. Examples: {preview}"
        )

    with open(args.output_gtf, "w") as out:
        for record in gtf_records:
            out.write(record + "\n")

    with open(args.output_map, "w") as out:
        out.write(
            "qseqid\ttranscript_id\tgene_id\tmag_id\t"
            "sseqid\tpident\tevalue\tbitscore\n"
        )

        for transcript_id in sorted(mappings):
            row = mappings[transcript_id]

            out.write(
                "\t".join([
                    row["qseqid"],
                    row["transcript_id"],
                    row["gene_id"],
                    row["mag_id"],
                    row["sseqid"],
                    row["pident"],
                    row["evalue"],
                    row["bitscore"],
                ]) + "\n"
            )

    print(f"Loaded DIAMOND best hits: {len(hits)}")
    print(f"Mapped transcripts: {len(mappings)}")
    print(f"Wrote GTF records: {len(gtf_records)}")
    print(f"Output GTF: {args.output_gtf}")
    print(f"Output map: {args.output_map}")


if __name__ == "__main__":
    main()
