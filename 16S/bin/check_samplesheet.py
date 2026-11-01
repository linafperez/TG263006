#!/usr/bin/env python3

"""
check_samplesheet.py

Validate a samplesheet for the Parkinson-16S pipeline.

This script performs structural and semantic validation of the user-provided
samplesheet before the pipeline starts.

Main responsibilities
---------------------
- Validate the samplesheet structure.
- Validate required columns.
- Validate sample identifiers.
- Validate FASTQ file paths.
- Validate paired-end consistency.
- Report all detected errors and warnings.
- Write a validated samplesheet if no errors are found.

Author:
Santiago Castañeda

Pipeline:
parkinson-16s
"""

import argparse
import re
import sys

from pathlib import Path

import pandas as pd


# =============================================================================
# Constants
# =============================================================================

REQUIRED_COLUMNS = {
    "sample",
    "fastq_1",
    "fastq_2",
}

OPTIONAL_COLUMNS = {
    "run",
    "control",
    "quant_reading"
}

VALID_FASTQ_EXTENSIONS = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq"
)

SAMPLE_NAME_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9_]*$"
)

MAX_SAMPLE_NAME_LENGTH = 100
MIN_FASTQ_SIZE = 10 * 1024

class ValidationError:
    """
    Represents a validation error detected in the samplesheet.
    """

    def __init__(
        self,
        sample,
        field,
        message,
        suggestion=None,
    ):

        self.sample = sample
        self.field = field
        self.message = message
        self.suggestion = suggestion

    def __str__(self):

        text = []

        text.append(f"ERROR")

        if self.sample:
            text.append(f"Sample: {self.sample}")

        if self.field:
            text.append(f"Field : {self.field}")

        text.append(self.message)

        if self.suggestion:
            text.append(f"Suggestion: {self.suggestion}")

        return "\n".join(text)


class ValidationWarning:
    """
    Represents a non-fatal validation warning.
    """

    def __init__(self, sample, field, message):

        self.sample = sample
        self.field = field
        self.message = message

    def __str__(self):

        text = []

        text.append("WARNING")

        if self.sample:
            text.append(f"Sample: {self.sample}")

        if self.field:
            text.append(f"Field : {self.field}")

        text.append(self.message)

        return "\n".join(text)


class SamplesheetValidator:
    """
    Validate a Parkinson-16S samplesheet.
    """

    def __init__(self, samplesheet, base_dir):

        self.samplesheet = Path(samplesheet)
        self.base_dir = Path(base_dir).resolve()

        self.df = None

        self.errors = []
        self.warnings = []

        self.fastq_paths = []
        self.fastq_names = {}
        self.validated_df = None

 # ==========================================================
 # Private helper methods
 # ==========================================================

    def _add_error(
            self,
            sample,
            field,
            message,
            suggestion=None,
        ):
            """
            Store a validation error.
            """

            self.errors.append(

                ValidationError(

                    sample=sample,
                    field=field,
                    message=message,
                    suggestion=suggestion,

                )

            )

    def _add_warning(
        self,
        sample,
        field,
        message,
    ):
        """
        Store a validation warning.
        """

        self.warnings.append(

            ValidationWarning(

                sample=sample,
                field=field,
                message=message,

            )

        )

    def _can_validate(self):
        """
        Return True if validation can continue.
        """

        return (
            self.df is not None
            and len(self.errors) == 0
        )

    def _resolve_path(self, value):
        """
        Resolve relative FASTQ paths against base_dir.
        """

        path = Path(str(value))

        if not path.is_absolute():
            path = self.base_dir / path

        return path.resolve()

    def _normalize_fastq_name(self, path: Path) -> str:
        """
        Remove read identifiers from a FASTQ filename.
        """

        name = path.name

        # Remove FASTQ extension
        name = re.sub(
            r"\.(fastq|fq)(\.gz)?$",
            "",
            name,
            flags=re.IGNORECASE
        )

        # Remove common read suffixes
        name = re.sub(
            r"([._-]?)(R?)(1|2)([._-]?001)?$",
            "",
            name,
            flags=re.IGNORECASE
        )

        return name

 # =============================================================================
 # Read samplesheet
 # =============================================================================

    def read_samplesheet(self):
        """
        Read the samplesheet into a pandas DataFrame.
        """

        try:
            self.df = pd.read_csv(self.samplesheet)

        except FileNotFoundError:

            self._add_error(
                None,
                "samplesheet",
                f"Samplesheet not found:\n{self.samplesheet}",
                "Verify the path provided with --input."
            )

        except pd.errors.EmptyDataError:

            self._add_error(
                None,
                "samplesheet",
                "The samplesheet is empty.",
                "Provide a non-empty CSV file."
            )

        except pd.errors.ParserError:

            self._add_error(
                None,
                "samplesheet",
                "The samplesheet could not be parsed as a CSV file.",
                "Verify that the file is a valid comma-separated table."
            )

        except Exception as e:

            self._add_error(
                None,
                "samplesheet",
                str(e)
            )

 # =============================================================================
 # Validate structure
 # =============================================================================

    def validate_structure(self):
        """
        Validate the overall structure of the samplesheet.
        """

        if self.df is None:
            return
    
        if self.df.empty:

            self._add_error(
                None,
                "samplesheet",
                "The samplesheet contains no rows.",
                "Add at least one sample."
            )

 # =============================================================================
 # Validate columns
 # =============================================================================

    def validate_columns(self):
        """
        Validate the samplesheet header.
        """

        if self.df is None:
            return

        columns = list(self.df.columns)

        duplicated = self.df.columns[self.df.columns.duplicated()]

        for column in duplicated:

            self._add_error(
                None,
                column,
                f"Column '{column}' appears more than once.",
                "Remove duplicated columns."
            )

        missing = REQUIRED_COLUMNS - set(columns)

        for column in sorted(missing):

            self._add_error(
                None,
                column,
                f"Required column '{column}' is missing.",
                f"Add the '{column}' column."
            )

        unknown = set(columns) - REQUIRED_COLUMNS - OPTIONAL_COLUMNS

        for column in sorted(unknown):

            self._add_warning(
                None,
                column,
                f"Column '{column}' is not recognized and will be ignored."
            )

 # =============================================================================
 # Validate sample names
 # =============================================================================

    def validate_sample_names(self):
        """
        Validate sample identifiers.
        """

        if not self._can_validate():
            return

        if "sample" not in self.df.columns:
            return

        for row, value in enumerate(self.df["sample"], start=2):

            if pd.isna(value):

                self._add_error(
                    None,
                    "sample",
                    f"Row {row}: sample name is missing.",
                    "Provide a unique sample identifier."
                )

                continue

            sample = str(value)

            if sample.strip() == "":

                self._add_error(
                    None,
                    "sample",
                    f"Row {row}: sample name is empty.",
                    "Provide a unique sample identifier."
                )

                continue

            if sample != sample.strip():

                self._add_warning(
                    sample,
                    "sample",
                    "Leading or trailing whitespace detected."
                )

                sample = sample.strip()

            if len(sample) > MAX_SAMPLE_NAME_LENGTH:

                self._add_error(
                    sample,
                    "sample",
                    f"Sample name exceeds {MAX_SAMPLE_NAME_LENGTH} characters.",
                    "Use a shorter sample identifier."
                )

            if not SAMPLE_NAME_PATTERN.fullmatch(sample):

                self._add_error(
                    sample,
                    "sample",
                    "Sample names must start with a letter or number and contain only letters, numbers and underscores.",
                    "Example: Sample01, PD_001, Control_12, 016287d9"
                )

 #============================================================================
 # Validate fastq paths
 # =============================================================================

    def validate_fastq_paths(self):
        """
        Validate FASTQ paths and store resolved paths.
        """

        if not self._can_validate():
            return

        self.fastq_paths = []

        for _, row in self.df.iterrows():

            sample = str(row["sample"]).strip()

            sample_files = {}

            for field in ("fastq_1", "fastq_2"):

                value = row[field]

                if pd.isna(value) or str(value).strip() == "":

                    self._add_error(
                        sample,
                        field,
                        f"Missing path for '{field}'.",
                        f"Provide the path to the {field} FASTQ."
                    )

                    continue

                path = self._resolve_path(value)

                sample_files[field] = path

                if not path.exists():

                    self._add_error(
                        sample,
                        field,
                        f"FASTQ file not found:\n{path}",
                        "Verify that the path is correct."
                    )

                    continue

                if not path.is_file():

                    self._add_error(
                        sample,
                        field,
                        f"Path is not a file:\n{path}",
                        "Provide the path to a FASTQ file."
                    )

            self.fastq_paths.append(
                (sample, sample_files)
            )

 # =============================================================================
 # Validate fastq extensions
 # =============================================================================

    def validate_fastq_extensions(self):
        """
        Validate FASTQ file extensions.
        """

        if not self._can_validate():
            return
    
        for sample, files in self.fastq_paths:

            for field, path in files.items():

                if path.name.lower().endswith(VALID_FASTQ_EXTENSIONS):
                    continue

                self._add_error(
                    sample,
                    field,
                    f"Invalid FASTQ extension: '{path.name}'.",
                    "Accepted extensions are .fastq.gz, .fq.gz, .fastq and .fq."
                )

 # =============================================================================
 # Validate fastq pairs
 # =============================================================================

    def validate_fastq_pairs(self):
        """
        Validate that R1 and R2 FASTQ files belong to the declared sample.

        The sample identifier must be present in both FASTQ filenames.
        This supports conventional filenames as well as datasets exported
        from QIIME2 where R1 and R2 may contain different internal numeric
        identifiers.
        """

        if not self._can_validate():
            return

        for sample, files in self.fastq_paths:

            if "fastq_1" not in files or "fastq_2" not in files:
                continue

            r1_name = files["fastq_1"].name
            r2_name = files["fastq_2"].name

            r1_pair_id = self._normalize_fastq_name(
                files["fastq_1"]
            )

            r2_pair_id = self._normalize_fastq_name(
                files["fastq_2"]
            )

            if r1_pair_id == r2_pair_id:
                continue

            self._add_error(
                sample,
                "fastq_pair",
                (
                    "FASTQ files do not appear to form a matching "
                    f"paired-end set for sample '{sample}'.\n"
                    f"R1: {r1_name}\n"
                    f"R2: {r2_name}\n"
                    f"Normalized R1: {r1_pair_id}\n"
                    f"Normalized R2: {r2_pair_id}"
                ),
                (
                    "Verify that fastq_1 and fastq_2 originate from "
                    "the same sequencing run or technical replicate."
                )
            )
           

 # =============================================================================
 # Validate duplicate FASTQs
 # =============================================================================

    def validate_duplicate_fastqs(self):
        """
        Validate that each FASTQ file is used only once.
        """

        if not self._can_validate():
            return

        seen = {}

        for sample, files in self.fastq_paths:

            for field, path in files.items():

                path = str(path)

                if path not in seen:

                    seen[path] = (sample, field)
                    continue

                previous_sample, previous_field = seen[path]

                self._add_error(

                    sample,

                    field,

                    (
                        f"This FASTQ file is already assigned to sample "
                        f"'{previous_sample}' ({previous_field})."
                    ),
    
                    "Each FASTQ file should belong to only one sample."

                )

 # =============================================================================
 # Validate file sizes
 # =============================================================================

    def validate_file_sizes(self):
        """
        Warn if FASTQ files are unusually small.
        """

        if not self._can_validate():
            return
    
        for sample, files in self.fastq_paths:

            for field, path in files.items():

                size = path.stat().st_size

                if size >= MIN_FASTQ_SIZE:
                    continue
    
                self._add_warning(

                    sample,

                    field,

                    (
                        f"FASTQ file is very small "
                        f"({size} bytes)."
                    )

                )

 # =============================================================================
 # Print report
 # =============================================================================

    def print_report(self):
        """
        Print validation summary.
        """

        print("\n" + "=" * 70)
        print(" Parkinson-16S samplesheet validation")
        print("=" * 70)

        if self.warnings:

            print(f"\nWarnings ({len(self.warnings)})")

            print("-" * 70)

            for warning in self.warnings:

                print(warning)
                print()

        if self.errors:

            print(f"\nErrors ({len(self.errors)})")

            print("-" * 70)

            for error in self.errors:

                print(error)
                print()
    
            print("=" * 70)
            print("Validation FAILED")
            print("=" * 70)

            sys.exit(1)

        print("\nNo validation errors detected.")

        print("=" * 70)

        print("Validation PASSED")

        print("=" * 70)

 # =============================================================================
 # Write output
 # ============================================================================= 

    def write_output(self, output_file):
         """
         Write the validated samplesheet.
         """

         if self.df is None:
             return

         self.df.to_csv(
             output_file,
             index=False
         )

 # =============================================================================
 # Run
 # ============================================================================= 

    def run(self, output_file):
        """
        Execute the complete validation workflow.
        """

        self.read_samplesheet()

        self.validate_structure()
        self.validate_columns()

        self.validate_sample_names()

        self.validate_fastq_paths()
        self.validate_fastq_extensions()
        self.validate_fastq_pairs()
        self.validate_duplicate_fastqs()
        self.validate_file_sizes()

        self.print_report()

        if not self.errors:
            self.write_output(output_file)


def main():

    parser = argparse.ArgumentParser(
        description="Validate a Parkinson-16S samplesheet."
    )

    parser.add_argument(
        "samplesheet",
        help="Input samplesheet."
    )

    parser.add_argument(
        "--base-dir",
        default=".",
        help="Base directory used to resolve relative FASTQ paths."
    )

    parser.add_argument(
        "--output",
        default="samplesheet.valid.csv",
        help="Validated samplesheet."
    )

    args = parser.parse_args()

    validator = SamplesheetValidator(
        args.samplesheet,
        args.base_dir
    )

    validator.run(args.output)


if __name__ == "__main__":
    main()
