context("Generate genomic regions")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

expected_regions <- GenomicRanges::GRanges(
  seqnames = c(
    "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3", "chr3", "chr4",
    "chr4", "chr5", "chr5", "chr6", "chr6", "chr7", "chr7", "chr8", "chr8",
    "chr9", "chr9", "chr10", "chr10", "chr11", "chr11", "chr12", "chr12",
    "chr13", "chr13", "chr14", "chr14", "chr15", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrX", "chrY", "chrM"
  ),
  ranges = IRanges::IRanges(
    start = c(
      1L, 100000001L, 200000001L, 1L, 100000001L, 200000001L, 1L, 100000001L,
      1L, 100000001L, 1L, 100000001L, 1L, 100000001L, 1L, 100000001L, 1L,
      100000001L, 1L, 100000001L, 1L, 100000001L, 1L, 100000001L, 1L,
      100000001L, 1L, 100000001L, 1L, 100000001L, 1L, 100000001L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 100000001L, 1L, 1L
    ),
    end = c(
      100000000L, 200000000L, 248956422L, 100000000L, 200000000L,
      242193529L, 100000000L, 198295559L, 100000000L, 190214555L, 100000000L,
      181538259L, 100000000L, 170805979L, 100000000L, 159345973L, 100000000L,
      145138636L, 100000000L, 138394717L, 100000000L, 133797422L, 100000000L,
      135086622L, 100000000L, 133275309L, 100000000L, 114364328L, 100000000L,
      107043718L, 100000000L, 101991189L, 90338345L, 83257441L, 80373285L,
      58617616L, 64444167L, 46709983L, 50818468L, 100000000L, 156040895L,
      57227415L, 16569L
    )
  ),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

test_that(
  desc = "Generate an object with all genomic regions given parameters",
  code = {

    output_regions <- generateGenomicRegions(ref_genome, 1E8)
    GenomeInfoDb::seqinfo(output_regions) <- GenomeInfoDb::seqinfo(expected_regions)
    expect_equal(output_regions, expected_regions)

  }
)
