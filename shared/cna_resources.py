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
"""Locate the reference resources that the Verdict module needs.

Verdict (ClairS-TO's Python port of ASCAT) needs, for the assembly the input BAM was aligned
to, a set of common germline SNP sites (one loci file and one allele file per contig), the GC
content around those sites and, optionally, their replication timing. A resource directory
has this layout::

    <cna_resource_dir>/
        loci_files/<loci prefix>chr1.txt ... <loci prefix>chrX.txt
        allele_files/<alleles prefix>chr1.txt ... <alleles prefix>chrX.txt
        GC_<name>.txt
        RT_<name>.txt            (optional; LogR is corrected for GC content only when absent)

The file name prefixes are discovered from the directory, so a resource set built for another
assembly (for example T2T-CHM13) can carry names that say so. Exactly one candidate is allowed
for each item; when several match, the directory holds more than one resource set and the
caller is told instead of one being picked silently. The GRCh38 set shipped with ClairS-TO
(``G1000_loci_hg38_chr*.txt``, ``G1000_alleles_hg38_chr*.txt``, ``GC_G1000_hg38.txt``,
``RT_G1000_hg38.txt``) resolves to exactly the paths ClairS-TO used before this module existed.
"""
import os

# Verdict covers the autosomes and chrX only; chrY and unplaced contigs are never used.
VERDICT_CONTIGS = ["chr" + str(a) for a in list(range(1, 23)) + ["X"]]

_FIRST_CONTIG = 'chr1'


def normalize_contig(contig):
    """Drop a leading 'chr' so contig names can be compared across files.

    Loci files (and therefore LogR/BAF files) carry the BAM's contig names, whereas the GC
    content and replication timing files distributed with ASCAT spell contigs either bare
    ('1', 'X') or with the prefix ('chr1'). Matching is done on the normalised name only; the
    contig names written to any output are never changed.
    """
    contig = str(contig)
    return contig[3:] if contig.startswith('chr') else contig


def _single_candidate(directory, predicate, what):
    """Return (path, None) for the unique file in `directory` satisfying `predicate`, else (None, error)."""
    if not os.path.isdir(directory):
        return None, "directory {} not found".format(directory)
    matches = sorted(f for f in os.listdir(directory) if predicate(f))
    if len(matches) == 0:
        return None, "no {} found in {}".format(what, directory)
    if len(matches) > 1:
        return None, "{} candidates for the {} in {} ({}); keep one resource set per directory".format(
            len(matches), what, directory, ', '.join(matches[:4]) + (', ...' if len(matches) > 4 else ''))
    return os.path.join(directory, matches[0]), None


def _per_contig_prefix(resource_dir, subdir, what):
    """Derive '<dir>/<subdir>/<prefix>' from the unique file named '<prefix>chr1.txt'."""
    suffix = _FIRST_CONTIG + '.txt'
    path, error = _single_candidate(
        os.path.join(resource_dir, subdir),
        lambda f: f.endswith(suffix) and len(f) > len(suffix),
        "{} file ending in '{}'".format(what, suffix))
    if error is not None:
        return None, error
    return path[:-len(suffix)], None


def resolve_cna_resources(resource_dir):
    """Resolve the Verdict resource files under `resource_dir`.

    Returns a dict with keys:
      loci_prefix, alleles_prefix -- paths up to but excluding '<contig>.txt', so callers
                                     append the contig themselves (GNU parallel and
                                     get_logr_and_baf.py both work with prefixes)
      gc_file                     -- complete path
      rt_file                     -- complete path, or None when the directory holds no RT_*.txt
      errors                      -- list of human-readable problems; empty when usable
    Items that could not be resolved are None.
    """
    errors = []

    loci_prefix, error = _per_contig_prefix(resource_dir, 'loci_files', 'loci')
    if error is not None:
        errors.append(error)

    alleles_prefix, error = _per_contig_prefix(resource_dir, 'allele_files', 'allele')
    if error is not None:
        errors.append(error)

    gc_file, error = _single_candidate(
        resource_dir, lambda f: f.startswith('GC_') and f.endswith('.txt'), "GC content file 'GC_*.txt'")
    if error is not None:
        errors.append(error)

    # Replication timing is optional: only several candidates is an error.
    rt_file, error = _single_candidate(
        resource_dir, lambda f: f.startswith('RT_') and f.endswith('.txt'), "replication timing file 'RT_*.txt'")
    if error is not None and 'candidates' in error:
        errors.append(error)

    return {
        'loci_prefix': loci_prefix,
        'alleles_prefix': alleles_prefix,
        'gc_file': gc_file,
        'rt_file': rt_file,
        'errors': errors,
    }


def loci_file(resources, contig):
    return '{}{}.txt'.format(resources['loci_prefix'], contig)


def alleles_file(resources, contig):
    return '{}{}.txt'.format(resources['alleles_prefix'], contig)


def _resource_candidates(resources, contigs):
    """Every file Verdict needs for `contigs`, in one list."""
    candidates = [loci_file(resources, ctg) for ctg in contigs]
    candidates += [alleles_file(resources, ctg) for ctg in contigs]
    candidates.append(resources['gc_file'])
    if resources['rt_file'] is not None:
        candidates.append(resources['rt_file'])
    return candidates


def missing_resource_files(resources, contigs):
    """Return the resource files that Verdict would need for `contigs` but that do not exist."""
    if resources['errors']:
        return []
    return [fn for fn in _resource_candidates(resources, contigs) if not os.path.isfile(fn)]


def unreadable_resource_files(resources, contigs):
    """Return the resource files that exist but that this process cannot open for reading.

    Existence is not enough. The reference files ship mode 0640, so a container started under a
    different uid -- which is what `docker run -u $(id -u):$(id -g)` gives -- can stat them but
    not read them. Undetected that costs an hour: every
    other check here passes, Verdict stays enabled, and the run fails much later inside
    alleleCounter, which reports a file it cannot open as one that "does not appear to exist".
    """
    if resources['errors']:
        return []
    return [fn for fn in _resource_candidates(resources, contigs)
            if os.path.isfile(fn) and not os.access(fn, os.R_OK)]


def last_locus_position(loci_fn):
    """Return the highest position in a loci file, or None if no position could be read.

    None covers every reason alike -- unreadable, empty, or a malformed second column -- so it is
    not evidence about the file's contents. Callers that need to tell "could not read it" from
    "read it and found nothing wrong" check `unreadable_resource_files` first.

    Loci files are sorted by position, so reading the tail is enough.
    """
    try:
        with open(loci_fn, 'rb') as fp:
            fp.seek(0, os.SEEK_END)
            file_size = fp.tell()
            read_size = min(4096, file_size)
            fp.seek(file_size - read_size)
            lines = [line for line in fp.read(read_size).decode().split('\n') if line.strip()]
        return int(lines[-1].split('\t')[1])
    except Exception:
        return None


def loci_beyond_reference(resources, contigs, contig_length_by_name):
    """Return [(contig, last locus, contig length)] for contigs whose loci run past the reference.

    Contig names alone cannot tell two assemblies apart: T2T-CHM13 uses the same chr1..chrX
    names as GRCh38, so a GRCh38 resource set passes a name check and is then applied to
    coordinates it does not describe. A locus beyond the end of its contig proves such a
    mismatch (the converse does not hold, so this only ever reports true mismatches).
    """
    out_of_range = []
    for ctg in contigs:
        contig_length = contig_length_by_name.get(ctg)
        if not contig_length:
            continue
        last_position = last_locus_position(loci_file(resources, ctg))
        if last_position is not None and last_position > contig_length:
            out_of_range.append((ctg, last_position, contig_length))
    return out_of_range
