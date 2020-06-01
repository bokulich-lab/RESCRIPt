# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="rescript",
    version='2020.6.0.dev0',
    packages=find_packages(),
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description="Reference sequence annotation and curation",
    license='BSD-3-Clause',
    url="https://github.com/nbokulich/RESCRIPt",
    entry_points={'qiime2.plugins': ['rescript=rescript.plugin_setup:plugin']},
    package_data={
        'rescript.tests': ['data/*'],
        'rescript.types.tests': ['data/*'],
        'rescript': ['citations.bib'],
    },
    zip_safe=False,
)
