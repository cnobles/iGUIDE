context("Calculate genomic coverage")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

input_sites <- GenomicRanges::GRanges(
  seqnames = c("chr3", "chr4", "chr6", "chr9", "chr18"),
  ranges = IRanges::IRanges(
    start = c(40843844, 33583665, 117347610, 53157926, 61874684),
    width = 1
  ),
  strand = c("-", "-", "-", "-", "+"),
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

expected_coverage <- data.frame(
  seqnames = factor(
    c("chr3", "chr4", "chr6", "chr9", "chr18"),
    levels = GenomeInfoDb::seqlevels(ref_genome)
  ),
  start = c(40843844L, 33583665L, 117347610L, 53157926L, 61874684L),
  end = c(50843843L, 43583664L, 127347609L, 63157925L, 71874683L),
  width = c(10000000L, 10000000L, 10000000L, 10000000L, 10000000L),
  readCountsPos = c(0L, 0L, 0L, 0L, 1L),
  readCountsNeg = c(1L, 1L, 1L, 1L, 0L)
)

test_that(
  desc = "Calculate the coverage or frequency within genomic regions",
  code = {

    output_coverage <- calcCoverage(input_sites, 10E6)
    expect_equal(output_coverage, expected_coverage)

  }
)
