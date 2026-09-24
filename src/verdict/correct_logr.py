import os
import sys
from argparse import ArgumentParser

import numpy as np
from scipy.interpolate import BSpline
from sklearn.linear_model import LinearRegression

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from shared.cna_resources import normalize_contig

# Window names that the GC content file must provide, in ASCAT's order.
GC_WINDOWS = ['25bp', '50bp', '100bp', '200bp', '500bp', '1kb', '2kb', '5kb', '10kb', '20kb', '50kb', '100kb',
              '200kb', '500kb', '1Mb']


def create_bspline_basis(x, df, degree=3):
    """Create B-spline basis for given data"""
    n_knots = df - degree + 1
    knots = np.linspace(np.min(x), np.max(x), n_knots)
    knots = np.concatenate(([knots[0]] * degree, knots, [knots[-1]] * degree))
    spline = BSpline(knots, np.eye(len(knots) - degree - 1), degree)
    # BSpline evaluates the whole vector at once, giving bitwise the same matrix as a
    # per-element loop but far faster, which matters at the ~3M loci of a WGS resource set.
    return spline(x)


def read_covariate_file(file_name, description):
    """Read an ASCAT GC content or replication timing file.

    Column 1 is a row identifier, column 2 the contig (with or without 'chr'), column 3 the
    position and the remaining columns the covariate values. Returns the value column names
    and a dict {(normalised contig, position): tab-joined values}.
    """
    values = dict()
    with open(file_name, 'r') as fp:
        header = fp.readline().rstrip('\n').split('\t')
        if len(header) < 4:
            sys.exit("[ERROR] The {} file {} has {} columns, expected a row id, Chr, Position and at least one "
                     "value column.".format(description, file_name, len(header)))
        for line in fp:
            info = line.rstrip('\n').split('\t')
            if len(info) < 4:
                continue
            pos = info[2]
            if 'e' in pos or '.' in pos:
                # the ASCAT hg38 files spell one position in scientific notation ("8e+06")
                pos = str(int(float(pos)))
            key = (normalize_contig(info[1]), pos)
            values[key] = '\t'.join(info[3:])
    return header[3:], values


def check_full_coverage(values, lookup_keys, file_name, description):
    """Exit with an actionable message if a covariate file misses loci that Verdict needs.

    Every later stage lines up LogR, BAF and germline genotypes by row, so loci cannot simply
    be dropped here; an incomplete file has to be reported instead.
    """
    missing = [key for key, lookup_key in lookup_keys.items() if lookup_key not in values]
    if not missing:
        return
    sys.exit(
        "[ERROR] The {} file {} does not cover all loci used by Verdict: {} of {} loci are missing "
        "(e.g. {}:{}). The GC content file (and the replication timing file, if given) must hold one row "
        "for every locus in the loci and allele files of the same CNA resource set. Please check "
        "docs/verdict.md for the expected file formats.".format(
            description, file_name, len(missing), len(lookup_keys), missing[0][0], missing[0][1]))


def abs_correlation_with_logr(covariates, logr, autosome_mask):
    """|corr| of every covariate column with LogR, computed on autosomes only as ASCAT does.

    NaN correlations (constant columns) are treated as 0 so they are never selected.
    """
    if autosome_mask.sum() >= 2:
        covariates = covariates[autosome_mask]
        logr = logr[autosome_mask]
    corr = np.corrcoef(covariates, logr, rowvar=False)[-1, :-1]
    return np.nan_to_num(np.abs(corr), nan=0.0)


def correctLogR(tumor_logr_file, gc_content_file, replication_timing_file, tumor_logr_correction_output_file, sample_name):
    has_rt = replication_timing_file is not None

    tumor_logr_dict = dict()
    with open(tumor_logr_file, 'r') as fp:
        for idx, tumor_logr_line in enumerate(fp):
            if idx == 0:
                continue
            tumor_logr_info = tumor_logr_line.strip().split('\t')
            chr = tumor_logr_info[0]
            pos = tumor_logr_info[1]
            logr = tumor_logr_info[2]
            key = (str(chr), str(pos))
            tumor_logr_dict[key] = logr

    gc_windows, gc_content_dict = read_covariate_file(gc_content_file, 'GC content')
    for window in ['1kb', '100kb', '1Mb']:
        if window not in gc_windows:
            sys.exit("[ERROR] The GC content file {} has no '{}' column. Expected the ASCAT window columns {}."
                     .format(gc_content_file, window, ' '.join(GC_WINDOWS)))
    if gc_windows != GC_WINDOWS:
        print("[WARNING] The GC content file {} has window columns {} instead of the ASCAT layout {}; "
              "windows are selected by name.".format(gc_content_file, ' '.join(gc_windows), ' '.join(GC_WINDOWS)))

    overlap_keys = list(tumor_logr_dict.keys())
    # The LogR file keeps the contig names of the loci files, the GC and replication timing
    # files may spell them either way, so match on the normalised name.
    lookup_keys = {key: (normalize_contig(key[0]), key[1]) for key in overlap_keys}
    check_full_coverage(gc_content_dict, lookup_keys, gc_content_file, 'GC content')

    tumor_logr_dict_values = np.array([tumor_logr_dict[key] for key in overlap_keys]).astype(float)
    gc_content_dict_values = np.array(
        [gc_content_dict[lookup_keys[key]].split('\t') for key in overlap_keys]).astype(float)
    autosome_mask = np.array([lookup_keys[key][0] not in ('X', 'Y') for key in overlap_keys])

    corr_gc = abs_correlation_with_logr(gc_content_dict_values, tumor_logr_dict_values, autosome_mask)

    # As in ASCAT: the insert-size window is searched among 25bp..1kb, the amplicon-size window
    # among 5kb..100kb, or up to 500kb when no replication timing is available.
    index_1kb = gc_windows.index('1kb')
    index_max = gc_windows.index('100kb') if has_rt else gc_windows.index('1Mb') - 1
    maxGCcol_insert = int(np.argmax(corr_gc[:(index_1kb + 1)]))
    maxGCcol_amplic = int(np.argmax(corr_gc[(index_1kb + 2):(index_max + 1)])) + (index_1kb + 2)
    print("[INFO] GC correction: insert-size window {} (|corr| {:.3f}), amplicon-size window {} (|corr| {:.3f})"
          .format(gc_windows[maxGCcol_insert], corr_gc[maxGCcol_insert],
                  gc_windows[maxGCcol_amplic], corr_gc[maxGCcol_amplic]))

    bases = [create_bspline_basis(gc_content_dict_values[:, maxGCcol_insert], df=5),
             create_bspline_basis(gc_content_dict_values[:, maxGCcol_amplic], df=5)]

    if has_rt:
        rt_names, replication_timing_dict = read_covariate_file(replication_timing_file, 'replication timing')
        check_full_coverage(replication_timing_dict, lookup_keys, replication_timing_file, 'replication timing')
        replication_timing_dict_values = np.array(
            [replication_timing_dict[lookup_keys[key]].split('\t') for key in overlap_keys]).astype(float)
        corr_rep = abs_correlation_with_logr(replication_timing_dict_values, tumor_logr_dict_values, autosome_mask)
        maxreplic = int(np.argmax(corr_rep))
        print("[INFO] Replication timing correction: dataset {} (|corr| {:.3f})".format(
            rt_names[maxreplic], corr_rep[maxreplic]))
        bases.append(create_bspline_basis(replication_timing_dict_values[:, maxreplic], df=5))
    else:
        print("[INFO] No replication timing file given, proceeding with GC correction only.")

    X = np.hstack(bases)
    y = tumor_logr_dict_values.reshape(-1, 1)
    model = LinearRegression(fit_intercept=True).fit(X, y)
    residuals = y - model.predict(X)
    tumor_logr_dict_values_after = residuals.flatten()

    corr_after = abs_correlation_with_logr(gc_content_dict_values, tumor_logr_dict_values_after, autosome_mask)
    print("[INFO] |corr(LogR, GC)| before/after correction: {} {:.3f}/{:.3f}, {} {:.3f}/{:.3f}".format(
        gc_windows[maxGCcol_insert], corr_gc[maxGCcol_insert], corr_after[maxGCcol_insert],
        gc_windows[maxGCcol_amplic], corr_gc[maxGCcol_amplic], corr_after[maxGCcol_amplic]))

    output_header = 'Chromosome' + '\t' + 'Position' + '\t' + sample_name + '\n'
    with open(tumor_logr_correction_output_file, 'w') as tumor_logr_correction_output:
        tumor_logr_correction_output.write(output_header)
        for idx, key in enumerate(overlap_keys):
            tumor_logr_correction_string = '\t'.join(map(str, key)) + '\t' + str(tumor_logr_dict_values_after[idx]) + '\n'
            tumor_logr_correction_output.write(tumor_logr_correction_string)


def main():
    parser = ArgumentParser(description="Correct Tumor Sample LogR")

    parser.add_argument('--tumor_logr_file', type=str,
                        default=None,
                        help="Path of tumor sample LogR")

    parser.add_argument('--gc_content_file', type=str,
                        default=None,
                        help="Path of 1kG GC content file")

    parser.add_argument('--replication_timing_file', type=str,
                        default=None,
                        help="Path of 1kG replication timing file. Optional; LogR is corrected for GC content only when omitted")

    parser.add_argument('--tumor_logr_correction_output_file', type=str,
                        default=None,
                        help="Output path of tumor sample corrected LogR")

    parser.add_argument('--sample_name', type=str,
                        default="SAMPLE",
                        help="Tumor sample name")


    global args
    args = parser.parse_args()

    correctLogR(args.tumor_logr_file, args.gc_content_file, args.replication_timing_file, args.tumor_logr_correction_output_file, args.sample_name)


if __name__ == "__main__":
    main()
