context("Genomic density")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr3", "chr3", "chr6", "chr9", "chr18"),
  ranges = IRanges::IRanges(
    start = c(40843844, 40843950, 117347610, 53157926, 61874684),
    width = 1
  ),
  strand = c("-", "-", "-", "-", "+"),
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

expected_regions <- GenomicRanges::GRanges(
  seqnames = c("chr3", "chr6", "chr9", "chr18"),
  ranges = IRanges::IRanges(
    start = c(40000001L, 110000001L, 50000001L, 60000001L),
    end = c(50000000L, 120000000L, 60000000L, 70000000L)
  ),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(ref_genome),
  count = c(2L, 1L, 1L, 1L),
  log.count = c(0.477, 0.301, 0.301, 0.301),
  norm.log.count = c(1, 0.631, 0.631, 0.631)
)

test_that(
  desc = "Determine the density of observations across the reference genome",
  code = {

    output <- genomicDensity(
      input_gr, res = 1E7, cutoff = 1L, drop.alt.chr = FALSE
    )

    output$log.count <- round(output$log.count, digits = 3)
    output$norm.log.count <- round(output$norm.log.count, digits = 3)
    expect_equal(output, expected_regions)

  }
)
