# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
"""Config parser classes"""

import os
import configparser


class PairedEndReadsConfiguration:
    def __init__(self, args):
        config = configparser.ConfigParser()
        config.read(args.configuration)

        self.output_sample_name = config.get('general', 'output_sample_name')
        self.insert_size = int(config.get('general', 'insert_size'))
        self.insert_size_std = float(config.get('general', 'insert_size_std'))
        self.short_read_length = int(config.get('general', 'short_read_length'))
        self.error_rate = float(config.get('general', 'error_rate'))

        self.fasta_files = []
        self.fasta_files_dict = {}

        for section in [s for s in config.sections() if s != 'general']:
            alias = os.path.basename('.'.join(section.split('.')[:-1]))
            self.fasta_files.append(alias)
            self.fasta_files_dict[alias] = {'path': os.path.join(os.path.dirname(os.path.abspath(args.configuration)), section),
                                            'alias': alias,
                                            'coverage': int(config.get(section, 'coverage'))}


class SingleReadsConfiguration:
    def __init__(self, args):
        config = configparser.ConfigParser()
        config.read(args.configuration)

        self.output_file = config.get('general', 'output_file')
        self.short_read_length = int(config.get('general', 'short_read_length'))
        self.error_rate = float(config.get('general', 'error_rate'))

        self.fasta_files = []
        self.fasta_files_dict = {}

        for section in [s for s in config.sections() if s != 'general']:
            alias = os.path.basename('.'.join(section.split('.')[:-1]))
            self.fasta_files.append(alias)
            self.fasta_files_dict[alias] = {'path': os.path.join(os.path.dirname(os.path.abspath(args.configuration)), section),
                                            'alias': alias,
                                            'coverage': int(config.get(section, 'coverage'))}
