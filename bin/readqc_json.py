#!/usr/bin/env python
from datetime import datetime
from pathlib import Path
import os
import glob
import json
import re
import argparse

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--output", "-o")
parser.add_argument("--run_name", help="run name")


args = parser.parse_args()

# JSON keys we want to remove since they are too large and unnecessary
UNWANTED_KEYS = [
    "content_curves",
    "kmer_count",
    "quality_curves"
]


def dict_cleaner(data, unwanted_keys=UNWANTED_KEYS):
    if not isinstance(data, dict):
        return data if not isinstance(data, list) else list(map(dict_cleaner, data))
    return {a: dict_cleaner(b) for a, b in data.items() if a not in unwanted_keys}


def parse_json(lines, return_adress=None):
    """
    return the JSON as dict
    if return_adress is defined, returns the value at the key adress
    """
    data = json.loads(" ".join(lines))

    data = dict_cleaner(data)

    if return_adress and data:
        for key in return_adress:
            data = data[key]
        return data
    return data


def parse_flat(lines, sep="\t", no_header=False):
    if no_header:
        header = [f"col{i}" for i in range(len(lines[0].strip().split(sep)))]
    else:
        header = lines.pop(0).strip().split(sep)
    data = []
    for line in lines:
        this_data = {}
        elements = line.strip().split(sep)
        for idx, h in enumerate(header):
            if idx < len(elements):
                entry = elements[idx]
                # value is an integer
                if re.match(r"^[0-9]+$", entry):
                    entry = int(entry)
                # value is a float
                elif re.match(r"^[0-9]+\.[0-9]+$", entry):
                    entry = float(entry)
                # value is a file path (messes up md5 fingerprinting)
                elif re.match(r"^\/.*\/.*$", entry):
                    entry = entry.split("/")[-1]
                this_data[h] = entry
        data.append(this_data)

    return data


def parse_yaml(lines):
    data = {}
    key = ""

    for line in lines:

        line = line.replace("\"", "")
        if re.match(r"^\s+.*", line):
            line = line.replace(":", "")
            elements = line.strip().split(" ")
            tool = elements.pop(0)
            version = " ".join(elements)
            data[key][tool] = version
        else:
            key = line.strip()
            data[key] = {}

    return data


def parse_interop_csv(lines):

    bucket = {}

    run = lines[1]
    header = lines[2].split(",")

    keys = { "read_1": 3, "read_4": 6, "total": 8 }

    for item,line in keys.items():

        elements = lines[line].split(",")
        this_bucket = {}
        for index, h in enumerate(header):
            this_bucket[h] = elements[index]    

        bucket[item] = this_bucket

    return bucket


def main(run_name, output):
    # Mapping each JSON section to (json_key, file regex, parsing_function, kwargs)
    parser_mapper = {
        "fastp": ("fastp", ".fastp.json", parse_json, None),
        "interop": ("interop", "interop.csv", parse_interop_csv, None),
        "versions": ("versions", "versions.yml", parse_yaml, None),
        "kraken2": ("kraken", "report.txt", parse_flat, None),
    }

    files = [os.path.abspath(f) for f in glob.glob("*")]

    date = datetime.today().strftime('%Y-%m-%d')

    matrix = {
        "date": date,
        "kraken2": {},
        "run_date": datetime.now().strftime('%Y-%m-%d'),
        "run_name": run_name
    }

    # Iterating over Path objects
    for file_path in map(Path, files):

        with open(file_path, "r") as f:
            lines = [line.rstrip() for line in f]

        # keep track of matched keys to skip
        matched_keys = set()

        for k, (json_key, suffix, func, kwargs) in parser_mapper.items():
            if k in matched_keys:
                continue
            if file_path.name.endswith(suffix):
                
                kwargs = kwargs or {}
                matrix[json_key] = func(lines, **kwargs)
                matched_keys.add(k)
                break

    with open(output, "w") as fo:
        json.dump(matrix, fo, indent=4, sort_keys=True)


if __name__ == '__main__':
    main( args.run_name, args.output)
