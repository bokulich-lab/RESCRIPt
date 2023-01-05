# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


def _validate_record_len(cells, current_line_number, exp_len):
    if len(cells) != exp_len:
        raise ValidationError(
            "Expected data record to be TSV with {0} "
            "fields. Detected {1} fields at line {2}:\n\n{3!r}".format(
                exp_len, len(cells), current_line_number, cells))


def _validate_is_numeric(inputvalue, valuedescription, line_number):
    try:
        float(inputvalue)
    except ValueError:
        raise ValidationError(
            "{0}must contain only numeric values. A non-numeric value "
            "({1!r}) was detected at line {2}.".format(
                valuedescription, inputvalue, line_number))


def _validate_file_not_empty(has_data):
    if not has_data:
        raise ValidationError(
            "There must be at least one tab-delimited data record present in "
            "the file (in addition to any header lines).")


def _validate_silva_taxonomy_format(inputvalue, columnnumber):
    if ';' not in inputvalue:
        raise ValidationError(
            "Column {0} (taxonomy) is not in SILVA taxonomy format. "
            "SILVA taxonomy consists of a semicolon-delimited "
            "string with a terminal semicolon.".format(columnnumber))


class SILVATaxonomyFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        with self.open() as fh:
            # note: there is no header line
            # validate body
            has_data = False
            for line_number, line in enumerate(fh, start=1):
                # we want to strip each cell, not the original line
                # otherwise empty cells are dropped, causing a TypeError
                cells = [c.strip() for c in line.split('\t')]
                # validate contents
                # must be 5 columns:
                # taxonomy, taxid, terminal rank, BLANK, silva version
                _validate_record_len(cells, line_number, 5)
                # first column must be semicolon-delimited taxonomy
                _validate_silva_taxonomy_format(cells[0], 1)
                # second column must be numeric (taxids)
                _validate_is_numeric(
                    cells[1], 'Column 2 (taxids) ', line_number)
                has_data = True
                if n_records is not None and (line_number - 1) >= n_records:
                    break

            _validate_file_not_empty(has_data)

    def _validate_(self, level):
        record_count_map = {'min': 5, 'max': None}
        self._validate(record_count_map[level])


SILVATaxonomyDirectoryFormat = model.SingleFileDirectoryFormat(
    'SILVATaxonomyDirectoryFormat', 'silva_taxonomy.tsv',
    SILVATaxonomyFormat)


class SILVATaxidMapFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        HEADER = ['primaryAccession', 'start', 'stop', 'path',
                  'organism_name', 'taxid']
        HEADER_VERSION132 = ['primaryAccession', 'start', 'stop', 'path',
                             'organismName', 'taxid']
        with self.open() as fh:
            # validate header
            # for now we will not validate any information in the header.
            line = [i.strip() for i in fh.readline().split('\t')]
            if line != HEADER and line != HEADER_VERSION132:
                raise ValidationError(
                    "Header line does not match SILVA format. Must consist of "
                    "the following values: " + ', '.join(HEADER) +
                    ".\n\nFound instead: " + ', '.join(line))
            # validate body
            has_data = False
            for line_number, line in enumerate(fh, start=1):
                # we want to strip each cell, not the original line
                # otherwise empty cells are dropped, causing a TypeError
                cells = [c.strip() for c in line.split('\t')]
                # validate contents
                # we will not validate the row index (primaryAccession), since
                # we assume that can be flexible. We do not validate start/stop
                # because we don't really use that anywhere.
                # must be 6 columns
                _validate_record_len(cells, line_number, 6)
                # fourth column must be semicolon-delimited taxonomy
                _validate_silva_taxonomy_format(cells[3], 4)
                # sixth column must be numeric (taxids)
                _validate_is_numeric(
                    cells[5], 'Column 6 (taxids) ', line_number)
                has_data = True
                if n_records is not None and (line_number - 1) >= n_records:
                    break

            _validate_file_not_empty(has_data)

    def _validate_(self, level):
        record_count_map = {'min': 5, 'max': None}
        self._validate(record_count_map[level])


SILVATaxidMapDirectoryFormat = model.SingleFileDirectoryFormat(
    'SILVATaxidMapDirectoryFormat', 'silva_taxmap.tsv',
    SILVATaxidMapFormat)
