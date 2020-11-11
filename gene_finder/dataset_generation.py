#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Copyright 2020 EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
"""


# standard library imports
import argparse
import pathlib
import sys

# third party imports
import pandas as pd

from Bio import SeqIO

# project imports


data_directory = pathlib.Path("data")


def fasta_sequences_to_dataframe(genomic_sequences_path, encoded_sequences_path):
    """
    """
    examples = []
    with open(genomic_sequences_path) as genomic_sequences, open(encoded_sequences_path) as encoded_sequences:
        for genomic_sequence, encoded_sequence in zip(SeqIO.FastaIO.SimpleFastaParser(genomic_sequences), SeqIO.FastaIO.SimpleFastaParser(encoded_sequences)):
            assert genomic_sequence[0] == encoded_sequence[0], f"{genomic_sequence=}, {encoded_sequence=}"
            examples.append({"description": genomic_sequence[0], "sequence": genomic_sequence[1], "encoding": encoded_sequence[1]})

    examples_dataframe = pd.DataFrame(examples)

    return examples_dataframe


def fasta_sequences_to_pickled_dataframe(genomic_sequences_path, encoded_sequences_path, dataframe_pickle_path):
    """
    """
    examples_dataframe = fasta_sequences_to_dataframe(genomic_sequences_path, encoded_sequences_path)

    examples_dataframe.to_pickle(dataframe_pickle_path)


def load_data():
    """
    """
    dataframe_pickle_path = data_directory / "sequences_dataframe.pickle"
    print("loading data...")
    data = pd.read_pickle(dataframe_pickle_path)
    print("data loaded")

    return data


def main():
    """
    main function
    """
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--fasta_sequences_to_pickled_dataframe", action="store_true")

    args = argument_parser.parse_args()

    if args.fasta_sequences_to_pickled_dataframe:
        genomic_sequences_path = data_directory / "run_1_genomic_seqs.fa"
        encoded_sequences_path = data_directory / "run_1_encoded_seqs.fa"
        dataframe_pickle_path = data_directory / "sequences_dataframe.pickle"

        fasta_sequences_to_pickled_dataframe(genomic_sequences_path, encoded_sequences_path, dataframe_pickle_path)
    else:
        print("nothing to do")


if __name__ == "__main__":
    main()
