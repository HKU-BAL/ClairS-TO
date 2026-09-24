About the Verdict module
---
A module named Verdict in ClairS-TO is applicable to the neural network called variants, and tags them as either a germline, somatic, or subclonal somatic variant. Verdict’s idea and algorithm are similar to and improved from SGZ [A computational approach to distinguish somatic vs. germline origin of genomic alterations from deep sequencing of cancer specimens without a matched normal, Sun et al., 2018]. Verdict has three steps: 1) find copy number segments, 2) estimate tumor purity, tumor ploidy, and the copy number profile, and 3) binomial tests.

SGZ suggested using ASCAT [Allele-specific copy number analysis of tumors, Van Loo et al., 2010] to estimate tumor purity, ploidy, and the copy number profile of each variant. We rewrote ASCAT in Python from R so that it could be integrated into Verdict and run reasonably fast. ASCAT uses the LogR (log ratios, representing log-transformed copy numbers derived from sequencing depth) and BAF (B allele frequencies, describing the allelic imbalance of variants) of germline heterozygous variants as input. LogR is calculated by the normalized read coverage of the tumor sample. BAF is calculated by dividing the signal intensities of minor alleles by those of major and minor alleles. ASCAT uses the LogR and BAF distributions to segment the genome into multiple regions with constant copy number states, identifying breakpoints based on LogR and BAF value changes. Next, ASCAT estimates tumor purity and ploidy by evaluating the goodness of fit for a grid of possible values for tumor purity and ploidy. Using the fact that true copy numbers are nonnegative whole numbers, ASCAT seeks values for tumor purity and ploidy such that the copy number estimates are as close as possible to nonnegative whole numbers for germline heterozygous variants. Finally, the allele frequency of each variant is used as input to two binomial tests that decide whether a variant is more likely to be germline, somatic, or subclonal somatic. The two binomial tests are calculated as follows using the tumor purity, tumor ploidy, and copy number as inputs. The p-value of the somatic hypothesis: p<sub>somatic</sub> = Binomial (n * f, n, AF<sub>somatic</sub>), where n is the read depth, f is the allele frequency, AF<sub>somatic</sub> is the expected allele frequency of the variant being somatic, calculated as p * V / (p * C + 2 * (1-p)), where p is the tumor purity, C is the copy number, and V is the variant allele count in the tumor. The p-value of the germline hypothesis: p<sub>germline</sub> = Binomial (n * f, n, AF<sub>germline</sub>), where AF<sub>germline</sub> is the expected allele frequency of the variant being a germline, calculated as (p * V + (1-p))/ (p * C + 2 * (1-p)). Verdict tags a variant as somatic if p<sub>somatic</sub> is greater than 0.001 and p<sub>germline</sub> is lower than 0.001. A variant is tagged as subclonal somatic if the p<sub>somatic</sub> and p<sub>germline</sub> are lower than 0.001, tumor purity is greater than 0.2, and f < AF<sub>somatic</sub> / gamma, where gamma is a tunable parameter with a default of 1.5.


## Using Verdict with another reference (e.g. T2T-CHM13)

Verdict itself is not tied to an assembly, but its resources are. It needs, for the reference the input
BAM was aligned to, a set of common germline SNP sites (the ASCAT "G1000" loci and allele files), the GC
content around those sites and, optionally, their replication timing. ClairS-TO ships such a set for
GRCh38 (`reference_files.tar.gz`, see the README). To run Verdict against another assembly, build the
equivalent set and point `--cna_resource_dir` at it. The assembly must name its contigs `chr1`..`chr22`
and `chrX`, as GRCh38 and T2T-CHM13 do; Verdict uses those contigs only.

Contig names alone cannot tell two assemblies apart: T2T-CHM13 uses the same `chr1`..`chrX` names as
GRCh38, so the GRCh38 resources would otherwise be accepted and applied to coordinates they do not
describe, leaving germline variants untagged. ClairS-TO therefore checks the loci against the contig
lengths of `--ref_fn` and disables Verdict, with a warning, when they cannot belong to that reference.

### Expected layout

```
<cna_resource_dir>/
  loci_files/<prefix>chr1.txt ... <prefix>chr22.txt, <prefix>chrX.txt
  allele_files/<prefix>chr1.txt ... <prefix>chr22.txt, <prefix>chrX.txt
  GC_<name>.txt
  RT_<name>.txt        (optional)
```

The prefixes are discovered from the file that ends in `chr1.txt` in each sub-directory, so the files can
carry any name that says which assembly they belong to (for example `G1000_loci_CHM13_chr1.txt`). Keep one
resource set per directory: if several files could be meant, ClairS-TO reports the ambiguity rather than
guessing. `chrY` and unplaced contigs are not used by Verdict.

When no `RT_*.txt` is present, LogR is corrected for GC content only, exactly as ASCAT does when no
replication timing file is given (the amplicon-size GC window is then searched up to 500kb instead of 100kb).

### File formats

All files are tab separated. Contig names in the GC content and replication timing files may be written
with or without the `chr` prefix; matching ignores the prefix.

* **loci** — no header, one row per site, 1-based positions, sorted by position, contig names as in the BAM:

      chr1	261839
      chr1	262034

* **alleles** — with header, the same sites in the same order, alleles encoded as `1=A, 2=C, 3=G, 4=T`,
  where `a0` is the reference allele of the assembly and `a1` the alternative allele:

      position	a0	a1
      261839	1	2

* **GC content** — with header. Column 1 is a row identifier, column 2 the contig, column 3 the position,
  followed by the 16 ASCAT window sizes in this order:

      	Chr	Position	25bp	50bp	100bp	200bp	500bp	1kb	2kb	5kb	10kb	20kb	50kb	100kb	200kb	500kb	1Mb
      1_261839	1	261839	0.6	0.529412	...

* **replication timing** (optional) — same first three columns, then one column per cell line:

      	Chr	Position	Bg02es	Bj	Gm06990	...
      1_261839	1	261839	54.427647	59.525208	58.798737	...

  The GC content file (and the replication timing file, if given) must cover **every** site present in
  the loci and allele files. ClairS-TO reports the first missing site if they do not.

### T2T-CHM13

ASCAT distributes CHM13 loci, allele and GC content files (`G1000_loci_WGS_CHM13.zip`,
`G1000_alleles_WGS_CHM13.zip`, `GC_G1000_WGS_CHM13.zip` on [Zenodo record 14008443](https://zenodo.org/records/14008443)),
built from the same 1000 Genomes SNPs as the GRCh38 set. Two things to know when laying them out:
the loci files are already `chr`-prefixed, and the allele files inside the zip are named `*_hg19_chr*.txt`
although they are CHM13 files, so rename them (e.g. to `G1000_alleles_CHM13_chr*.txt`). No replication
timing file is distributed for CHM13; Verdict then runs with GC-only correction.

Replication timing can be approximated for CHM13 by lifting the hg38 file over, but on a whole-genome
PacBio HiFi sample it changed almost nothing: identical tumour purity, ploidy 2.8738 against 2.8732,
92% of copy number segments identical, and a different tag on 106 of 4,385,912 variants (0.002%).
GC-only correction is therefore the recommended configuration for CHM13.

Everything else in a CHM13 run also has to match the assembly. The default non-somatic tagging databases
(gnomAD, dbSNP, 1000G PoN, CoLoRSdb) are GRCh38 coordinates, so pass `--disable_nonsomatic_tagging` or
CHM13-based files via `--panel_of_normals`. Do not pass the GRCh38 indel BED. Under Docker/Apptainer add
`--conda_prefix /opt/micromamba/envs/clairs-to` so that the models are found. A minimal example:

```bash
run_clairs_to -T tumor.chm13.bam -R chm13v2.0.fa -o output -t 24 -p hifi_revio \
    --cna_resource_dir /path/to/verdict_CHM13 --disable_nonsomatic_tagging
```
