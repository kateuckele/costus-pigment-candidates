#!/usr/bin/env python3ls
"""
blast_xml_to_table.py
---------------------
Convert BLAST XML output into a tab-delimited table.

Usage:
    ./blast_xml_to_table.py input.xml output_prefix

Example:
    ./blast_xml_to_table.py /results/blast_tables/0G799C6R013-Alignment.xml results/my_blast

Output:
    results/my_blast.tsv
"""

import sys
from Bio.Blast import NCBIXML

def main():
    if len(sys.argv) != 3:
        print("Usage: blast_xml_to_table.py input.xml output_prefix", file=sys.stderr)
        sys.exit(1)

    xml_file = sys.argv[1]
    out_prefix = sys.argv[2]
    out_file = f"{out_prefix}.tsv"

    with open(xml_file) as xml, open(out_file, "w") as out:
        blast_records = NCBIXML.parse(xml)
        out.write("query_id\tsubject_id\tpercent_identity\talignment_length\tevalue\tbit_score\n")

        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    pid = 100 * hsp.identities / hsp.align_length
                    out.write(
                        f"{record.query}\t{alignment.hit_def}\t{pid:.2f}\t"
                        f"{hsp.align_length}\t{hsp.expect}\t{hsp.bits}\n"
                    )

    print(f"[DONE] Wrote: {out_file}")

if __name__ == "__main__":
    main()
