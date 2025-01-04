# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'


importlib.import_module('rescript.types')
