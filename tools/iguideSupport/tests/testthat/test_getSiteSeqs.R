context("Get Site specific sequences")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

sites <- GenomicRanges::GRanges(
  seqnames = c("chr4", "chr15", "chr14", "chr7", "chr7"),
  ranges = IRanges::IRanges(
    start = c(105272459, 44711569, 22547664, 142792020, 142801367),
    width = 1
  ),
  strand = c("+", "-", "+", "+", "+"),
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

target_seq <- c(
  paste0(
    "AATTTGTCTTATCTCCCTGTACCATTTTGTTGCTATTTTCATTAATAACAGGTAGGATGGTTTTATGGTA",
    "ATATATATGTCACTGATCTGGATCAACTAGGCCACCAACACAAATCTGAATACTGAGAGGA"
  ),
  paste0(
    "CAGAGGGTGCAGAGCGGGAGAGGAAGGACCAGAGCGGGAGGGTAGGAGAGACTCACGCTGGATAGCCTCC",
    "AGGCCAGAAAGAGAGAGTAGCGCGAGCACAGCTAAGGCCACGGAGCGAGACATCTCGGCCC"
  ),
  paste0(
    "GTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTAT",
    "ATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGG"
  ),
  paste0(
    "GAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAG",
    "TTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGA"
  ),
  paste0(
    "GAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAG",
    "TTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGA"
  )
)

test_that(
  desc = "Query sequences around a specific genomic site",
  code = {

    query <- getSiteSeqs(
      sites,
      upstream.flank = 100,
      downstream.flank = 30,
      ref.genome = ref_genome
    )

    expect_equal(as.character(query), target_seq)

  }
)
