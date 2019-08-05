context("Calculate edit cut site")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

seqs <- Biostrings::DNAStringSet(c(
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
))

sites <- GenomicRanges::GRanges(
  seqnames = c("chr15", "chr14", "chr7", "chr7"),
  ranges = IRanges::IRanges(
    start = c(44711569, 22547664, 142792020, 142801367),
    width = 1
  ),
  strand = c("-", "+", "+", "+"),
  seqs = seqs,
  siteID = c(2L, 3L, 4L, 5L),
  target.match = c(
    "B2M:(sense)", "TRAC:(sense)", "TRBC:(sense)", "TRBC:(sense)"
  ),
  target.mismatch = 0,
  aligned.sequence = c(
    "GAGTAGCGCGAGCACAGCTAAGG", "TGTGCTAGACATGAGGTCTATGG",
    "GGAGAATGACGAGTGGACCCAGG", "GGAGAATGACGAGTGGACCCAGG"
  ),
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

matched_seqs <- data.frame(
  "names" = as.character(c(2, 3, 4, 5)),
  "target" = c("B2M:(sense)", "TRAC:(sense)", "TRBC:(sense)", "TRBC:(sense)"),
  "target.mismatch" = 0,
  "start" = 85L,
  "end" = 107L,
  "width" = 20.0,
  "nt_width" = 131L,
  "aln.seq" = c(
    "GAGTAGCGCGAGCACAGCTAAGG", "TGTGCTAGACATGAGGTCTATGG",
    "GGAGAATGACGAGTGGACCCAGG", "GGAGAATGACGAGTGGACCCAGG"
  ),
  stringsAsFactors = FALSE
)

nuc_profile <- list(
  "PAM" = "NGG",
  "PAM_Loc" = "3p",
  "PAM_Tol" = 1,
  "Cut_Offset" = -4,
  "Insert_size" = FALSE
)

expected_output <- c(
  "chr15:-:44711569", "chr14:+:22547664", "chr7:+:142792020", "chr7:+:142801367"
)

test_that(
  desc = "Calculate the location of a target within a sequence reference",
  code = {

    output <- calcCutSite(
      sites = sites,
      matched.seqs = matched_seqs,
      nuc.profile = nuc_profile,
      upstream.flank = 100L,
      downstream.flank = 30L
    )

    expect_equal(output, expected_output)

  }
)

