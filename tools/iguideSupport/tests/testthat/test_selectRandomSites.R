context("Select random sites")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

expected_sites <- GenomicRanges::GRanges(
  seqnames = c("chr3", "chr4", "chr6", "chr9", "chr18"),
  ranges = IRanges::IRanges(
    start = c(79310131L, 85229418L, 129367534L, 21612295L, 73949913L),
    width = 1
  ),
  strand = c("+", "+", "+", "-", "-"),
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

test_that(
  desc = "Select a number of random sites across the reference genome",
  code = {

    select_sites <- selectRandomSites(5, ref_genome, rnd.seed = 1)
    expect_equal(select_sites, expected_sites)

  }
)
