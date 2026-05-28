• The source MGF has MERGED_STATS lines but no MERGED_SCANS, and this GitHub mascot-rs version requires both when merged-scan metadata is present. Since those fields are not part of the taxonomy enrichment and mascot-rs does not write them back in its
  canonical output anyway, I’m going to preprocess only the in-memory parse input to skip MERGED_STATS lines for this dataset, then still parse and write via mascot-rs.


• The error means mascot-rs was parsing an MSLEVEL=1 ion block and enforcing its MS1 invariant: for level 1, the precursor m/z must equal the first m/z after mascot-rs sorts the peaks. In your file, the first block is:

  MSLEVEL=1
  PEPMASS=593.2760620117188
  593.2760620117188 ...
  592.2678833007812 ...

  After sorting, the minimum peak is 592.267883..., so mascot-rs rejects it as an MS1 precursor mismatch.

  For this tool, those MS1 SPECTYPE=CORRELATED MS blocks are not the target. The file is features_ms2, and the useful spectra are the following MSLEVEL=2 fragmentation blocks. I fixed the loader to:

  - keep only MSLEVEL=2 ion blocks before handing data to mascot-rs
  - drop MERGED_STATS and MERGED_SCANS during parsing because some PF1600 files have incomplete/mismatched merged-scan metadata that mascot-rs validates strictly, and those fields are not needed for taxonomy enrichment