context("Compare target sequences")

ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

seqs <- Biostrings::DNAStringSet(c(
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
))

sites <- GenomicRanges::GRanges(
  seqnames = c("chr4", "chr15", "chr14", "chr7", "chr7"),
  ranges = IRanges::IRanges(
    start = c(105272459, 44711569, 22547664, 142792020, 142801367),
    width = 1
  ),
  strand = c("+", "-", "+", "+", "+"),
  seqs = seqs,
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

tar_seqs <- c(
  "B2M" = "GAGTAGCGCGAGCACAGCTA",
  "TRAC" = "TGTGCTAGACATGAGGTCTA",
  "TRBC" = "GGAGAATGACGAGTGGACCC"
)

nuc_profile <- list(
  "PAM" = "NGG",
  "PAM_Loc" = "3p",
  "PAM_Tol" = 1,
  "Cut_Offset" = -4,
  "Insert_size" = FALSE
)

expected_output <- sites

expected_output$target.match <- c(
  "No_valid_match",
  "B2M:(sense)", "TRAC:(sense)", "TRBC:(sense)", "TRBC:(sense)"
)

expected_output$target.mismatch <- c(NA, 0, 0, 0, 0)

expected_output$aligned.sequence <- c(
  NA, "GAGTAGCGCGAGCACAGCTAAGG", "TGTGCTAGACATGAGGTCTATGG",
  "GGAGAATGACGAGTGGACCCAGG", "GGAGAATGACGAGTGGACCCAGG"
)

expected_output$edit.site <- c(
  NA, "chr15:-:44711569", "chr14:+:22547664", "chr7:+:142792020", "chr7:+:142801367"
)


test_that(
  desc = "Compare site-specific flanking sequences for targets",
  code = {

    output <- compareTargetSeqs(
      gr.with.sequences = sites,
      seq.col = "seqs",
      target.seqs = tar_seqs,
      nuc.profile = nuc_profile,
      upstream.flank = 100,
      downstream.flank = 30
    )

    expect_equal(output, expected_output)

  }
)

