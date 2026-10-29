#!/usr/bin/awk -f
# Reads `biscuit pileup` VCF (any number of samples) on stdin and prints one
# line per site whose alternative-allele evidence reaches min_alt reads,
# pooled over the samples:
#
#   chrom  pos (1-based)  ref  alt bases  depth  alt reads  samples with alt
#
# Only BISCUIT's definite base calls (A, C, G, T in the SP field) count. Its
# Y and R calls are reads on the strand where bisulfite conversion makes C/T
# (or G/A) indistinguishable, so they say nothing about a SNP and are left
# out of both depth and alt reads. Portable awk (busybox, BWK, gawk): no
# gensub, no arrays of arrays.
#
#   biscuit pileup ... | awk -v min_alt=2 -f pileup_sites.awk
BEGIN { FS = OFS = "\t"; if (min_alt == "") min_alt = 2 }
/^#/ { next }
$5 == "." { next }
{
    ref = toupper($4)
    if (ref !~ /^[ACGT]$/) next
    nf = split($9, fmt, ":")
    sp = 0
    for (i = 1; i <= nf; i++) if (fmt[i] == "SP") sp = i
    if (!sp) next
    depth = 0; alt = 0; nsamples = 0; seen = ""
    for (c = 10; c <= NF; c++) {
        n = split($c, v, ":")
        if (sp > n) continue
        s = v[sp]; here = 0
        while (match(s, /^[A-Z][0-9]+/)) {
            b = substr(s, 1, 1); k = substr(s, 2, RLENGTH - 1) + 0
            s = substr(s, RLENGTH + 1)
            if (b == "A" || b == "C" || b == "G" || b == "T") {
                depth += k
                if (b != ref) { alt += k; here += k; if (index(seen, b) == 0) seen = seen b }
            }
        }
        if (here > 0) nsamples++
    }
    if (alt >= min_alt) print $1, $2, ref, seen, depth, alt, nsamples
}
