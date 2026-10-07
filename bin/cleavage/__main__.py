"""Command line: python3 -m cleavage {samplesheet, call, finalize} ...

  samplesheet  check the input sheet and write samples.json
  call         measure chunk N of a sample
  finalize     merge a sample's chunks and write its outputs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .call import call_chunk
from .finalize import finalize
from .samplesheet import validate_samplesheet, write_samples_json
from .settings import load_settings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python3 -m cleavage")
    commands = parser.add_subparsers(dest="command", required=True)

    sheet = commands.add_parser("samplesheet", help="check the input samplesheet")
    sheet.add_argument("input_csv")
    sheet.add_argument("output_json")
    sheet.add_argument("--analysis", required=True, choices=["digenome", "ndigenome"])

    call = commands.add_parser("call", help="measure one chunk of one sample")
    call.add_argument("--settings", required=True)
    call.add_argument("--bam", required=True)
    call.add_argument("--sample", required=True)
    call.add_argument("--chunk", required=True, type=int, help="0-based chunk index")
    call.add_argument("--chunks", required=True, type=int, help="number of chunks")
    call.add_argument("--out-prefix", required=True)
    call.add_argument("--control-bam")
    call.add_argument("--control-sample", default="")
    call.add_argument("--vcf")
    call.add_argument("--blacklist")

    final = commands.add_parser("finalize", help="merge chunks and write outputs")
    final.add_argument("--settings", required=True)
    final.add_argument("--sample", required=True)
    final.add_argument("--control-sample", default="")
    final.add_argument("--out-prefix", required=True)
    final.add_argument("summaries", nargs="+", help="chunk summary JSON files")
    return parser


def run(args: argparse.Namespace) -> None:
    if args.command == "samplesheet":
        samples = validate_samplesheet(Path(args.input_csv), args.analysis)
        write_samples_json(samples, Path(args.output_json))
        print(f"Validated {len(samples)} sample(s) for {args.analysis} analysis.")
    elif args.command == "call":
        summary = call_chunk(
            load_settings(args.settings),
            bam_path=args.bam,
            sample=args.sample,
            index=args.chunk,
            count=args.chunks,
            out_prefix=args.out_prefix,
            control_path=args.control_bam,
            control_sample=args.control_sample,
            vcf_path=args.vcf,
            blacklist_path=args.blacklist,
        )
        print(f"{args.sample} chunk {args.chunk}/{args.chunks}: {summary['records']} record(s).")
    else:
        qc = finalize(
            load_settings(args.settings),
            args.summaries,
            args.sample,
            args.control_sample,
            args.out_prefix,
        )
        print(f"{args.sample}: {qc['rows']} candidate(s), {qc['tiers']['high_confidence']} high-confidence.")


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
