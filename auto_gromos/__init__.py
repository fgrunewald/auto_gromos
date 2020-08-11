# Copyright 2020 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pbr.version

__version__ = pbr.version.VersionInfo('auto_gromos').release_string()

# Find the data directory once.
try:
    import pkg_resources
except ImportError:
    import os
    DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')
    del os
else:
    DATA_PATH = pkg_resources.resource_filename('auto_gromos', 'data')
    del pkg_resources

del pbr

from .assign_functional_groups import AssignFunctionalGroups
from .select_dihedrals import SelectDihedrals
from .gen_pairs import Gen14Pairs
from .remove_duplicates import RemoveDuplicates
