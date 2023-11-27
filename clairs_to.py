# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
from importlib import import_module

DATA_PREP_SCRIPTS_FOLDER = "src"
DEEP_LEARNING_FOLDER = "clairs"

deep_learning_folder = [
    "predict",
    "call_variants",
]

REPO_NAME = "clairs_to"
data_preprocess_folder = [
    "extract_candidates_calling",
    "create_tensor_pileup_calling",
    "sort_vcf",
    "postprocess_vcf",
    'concat_files',
    'add_back_missing_variants_in_genotyping',
    "nonsomatic_tagging",
    'realign_reads',
    'realign_variants',
    "haplotype_filtering",
    'cal_af_distribution',
    "compare_vcf",
]


def directory_for(submodule_name):
    if submodule_name in deep_learning_folder:
        return DEEP_LEARNING_FOLDER
    if submodule_name in data_preprocess_folder:
        return DATA_PREP_SCRIPTS_FOLDER
    return ""


def print_help_messages():
    from textwrap import dedent
    print(dedent("""\
        {0} submodule invocator:
            Usage: python clairs_to.py [submodule] [options of the submodule]
        Available data preparation submodules:\n{1}
        Available clairs_to submodules:\n{2}
        """.format(
            REPO_NAME,
            "\n".join("          - %s" % submodule_name for submodule_name in data_preprocess_folder),
            "\n".join("          - %s" % submodule_name for submodule_name in deep_learning_folder),
        )
    ))


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print_help_messages()
        sys.exit(0)

    submodule_name = sys.argv[1]
    if (
        submodule_name not in deep_learning_folder and
        submodule_name not in data_preprocess_folder
    ):
        sys.exit("[ERROR] Submodule %s not found." % (submodule_name))

    directory = directory_for(submodule_name)
    submodule = import_module("%s.%s" % (directory, submodule_name))

    sys.argv = sys.argv[1:]
    sys.argv[0] += (".py")

    # Note: need to make sure every submodule contains main() method
    submodule.main()

    sys.exit(0)


if __name__ == "__main__":
    main()