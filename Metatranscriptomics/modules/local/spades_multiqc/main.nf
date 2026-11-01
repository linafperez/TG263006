process SPADES_MULTIQC {

    tag "${assembly_id}"
    label 'process_low'

    input:
    path transcripts
    val assembly_id

    output:
    path "${assembly_id}_spades_summary_mqc.yaml", emit: yaml

    script:
    """
    awk '
        /^>/ {
            if (seen) print len
            seen = 1
            len = 0
            next
        }
        {
            gsub(/[[:space:]]/, "")
            len += length(\$0)
        }
        END {
            if (seen) print len
        }
    ' ${transcripts} > transcript_lengths.txt

    transcript_count=\$(wc -l < transcript_lengths.txt | tr -d ' ')
    total_bases=\$(awk '{s += \$1} END {print s + 0}' transcript_lengths.txt)

    longest_transcript=\$(sort -nr transcript_lengths.txt | head -n 1)
    longest_transcript=\${longest_transcript:-0}

    half=\$(( (total_bases + 1) / 2 ))

    n50=\$(sort -nr transcript_lengths.txt | awk -v half="\$half" '
        {
            cumulative += \$1
            if (cumulative >= half) {
                print \$1
                exit
            }
        }
    ')
    n50=\${n50:-0}

    cat > ${assembly_id}_spades_summary_mqc.yaml << EOF_MQC
id: "spades_assembly_summary"
section_name: "rnaSPAdes - Assembly"
description: "Assembly metrics calculated directly from the transcript FASTA produced by rnaSPAdes."
plot_type: "table"
pconfig:
  id: "spades_assembly_summary_table"
  title: "rnaSPAdes Assembly Summary"
  col1_header: "Assembly"
headers:
  transcripts:
    title: "Transcripts"
    description: "Number of assembled transcript sequences"
  total_bases:
    title: "Total bases"
    description: "Total number of assembled transcript bases"
  longest_transcript:
    title: "Longest transcript"
    description: "Length of the longest assembled transcript in bases"
  n50:
    title: "N50"
    description: "Transcript assembly N50 in bases"
data:
  "${assembly_id}":
    transcripts: \${transcript_count}
    total_bases: \${total_bases}
    longest_transcript: \${longest_transcript}
    n50: \${n50}
EOF_MQC
    """

    stub:
    """
    cat > ${assembly_id}_spades_summary_mqc.yaml << EOF_MQC
id: "spades_assembly_summary"
section_name: "rnaSPAdes - Assembly"
description: "Assembly metrics calculated directly from the transcript FASTA produced by rnaSPAdes."
plot_type: "table"
pconfig:
  id: "spades_assembly_summary_table"
  title: "rnaSPAdes Assembly Summary"
  col1_header: "Assembly"
headers:
  transcripts:
    title: "Transcripts"
  total_bases:
    title: "Total bases"
  longest_transcript:
    title: "Longest transcript"
  n50:
    title: "N50"
data:
  "${assembly_id}":
    transcripts: 1
    total_bases: 100
    longest_transcript: 100
    n50: 100
EOF_MQC
    """
}
