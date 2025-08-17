"""
Simple section extractor for PyFrag input files.

This module provides basic regex-based extraction of sections from input files.
"""

import re
from pathlib import Path
from typing import Dict, Optional


def extract_sections(input_file: str) -> Dict[str, str]:
    """Extract sections from a PyFrag input file using regex patterns.

    Args:
        input_file: Path to the input file

    Returns:
        Dictionary mapping section names (lowercase) to their content as strings
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    with open(input_path, "r") as f:
        content = f.read()

    sections = {}

    # Extract standard sections (case-insensitive)
    section_patterns = [
        ("jobsub", r"(?i)^JOBSUB\s*(?:#.*)?$.*?^JOBSUB\s+END\s*(?:#.*)?$"),
        ("ams", r"(?i)^AMS\s*(?:#.*)?$.*?^AMS\s+END\s*(?:#.*)?$"),
        ("orca", r"(?i)^ORCA\s*(?:#.*)?$.*?^ORCA\s+END\s*(?:#.*)?$"),
        ("gaussian", r"(?i)^GAUSSIAN\s*(?:#.*)?$.*?^GAUSSIAN\s+END\s*(?:#.*)?$"),
        ("turbomole", r"(?i)^TURBOMOLE\s*(?:#.*)?$.*?^TURBOMOLE\s+END\s*(?:#.*)?$"),
        ("adf", r"(?i)^ADF\s*(?:#.*)?$.*?^ADF\s+END\s*(?:#.*)?$"),
        ("pyfrag", r"(?i)^PyFrag\s*(?:#.*)?$.*?^PyFrag\s+END\s*(?:#.*)?$"),
        ("fragment1_extra", r"(?i)^fragment1\s+EXTRA\s*(?:#.*)?$.*?^fragment1\s+EXTRA\s+END\s*(?:#.*)?$"),
        ("fragment2_extra", r"(?i)^fragment2\s+EXTRA\s*(?:#.*)?$.*?^fragment2\s+EXTRA\s+END\s*(?:#.*)?$"),
        ("fragment1_open_extra", r"(?i)^fragment1\s+OPEN\s+EXTRA\s*(?:#.*)?$.*?^fragment1\s+OPEN\s+EXTRA\s+END\s*(?:#.*)?$"),
        ("fragment2_open_extra", r"(?i)^fragment2\s+OPEN\s+EXTRA\s*(?:#.*)?$.*?^fragment2\s+OPEN\s+EXTRA\s+END\s*(?:#.*)?$"),
        ("complex_extra", r"(?i)^complex\s+EXTRA\s*(?:#.*)?$.*?^complex\s+EXTRA\s+END\s*(?:#.*)?$"),
    ]

    for section_name, pattern in section_patterns:
        match = re.search(pattern, content, re.MULTILINE | re.DOTALL)
        if match:
            # Extract content between header and footer lines
            section_content = match.group(0)
            lines = section_content.split("\n")
            lines = [line for line in lines if line.strip()]  # Remove empty lines

            # Remove first line (header) and last line (footer)
            if len(lines) >= 2:
                content_lines = lines[1:-1]
            else:
                content_lines = []

            sections[section_name] = "\n".join(content_lines)

    return sections


def get_section_content(input_file: str, section_name: str) -> Optional[str]:
    """Get the content of a specific section from the input file.

    Args:
        input_file: Path to the input file
        section_name: Name of the section to extract (case-insensitive)

    Returns:
        Section content as string, or None if section not found
    """
    sections = extract_sections(input_file)
    return sections.get(section_name.lower())
