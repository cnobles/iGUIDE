.. _configinfo:

.. contents::
   :depth: 4

Config - Nuclease Profiles
==========================

An additional component to the first part of the config file, is the Nuclease
Profiles. The user can specify which nuclease they are using and include
and profile to help identify edit sites. Nuclease can range from Cas9 to Cpf1
or TALEN based nickases. 

**Note:** For TALEN and dual flanking nickases / nucleases, each side will need
to be input as a different target. Specify in `Target_Sequences` the sequence
and `On_Target_Sites` the actual editing site. Make sure you include two 
distinct identifiers for the sequences on-target sites, then specify the 
target treatment as `{target_seq1};{target_seq2}.

Any name can be given in the `Nuclease` section, but that name needs to match
the profile name as well. So if you want to call it "Cas9v2", then just make 
sure you have a profile named "Cas9v2".

Below is some ascii art that indicates the differences between nucleases. 
Additionally, below the art are example profiles for input into the iGUIDE 
software.::

  Editing strategies by designer nucleases
  
  Cas9 :
                   ><   PAM
  ATGCATGCATGCATGCATGCA TGG (sense strand)

   TGCATGCATGCATGCATGCA NGG # gRNA
   |||||||||||||||||||| |||
  TACGTACGTACGTACGTACGT ACC (anti-sense strand)
                   ><       # Dominant cutpoint
  
  
  
  Cpf1 : Also known as Cas12a
                          ><         # Dominant cutpoint
  GTTTG ATGCATGCATGCATGCATGCATGCATGC (sense strand)
    PAM
   TTTV ATGCATGCATGCATGCATGCA        # gRNA, nuclease activity leave overhang
   |||| |||||||||||||||||||||
  CTAAC TACGTACGTACGTACGTACGTACGTACG (anti-sense strand)
                              ><     # Dominant cutpoint
  
  
  
  TALEN : Protin-DNA binding domain fused with FokI nickase
  
  ATATATATATATATATATAT GCATGCATGCATGCAT GCGCGCGCGCGCGCGCGCGC (sense strand)
  \\\\\\\\\\\\\\\\\\\\
                      |------->
                               <-------|
                                        \\\\\\\\\\\\\\\\\\\\
  TATATATATATATATATATA CGTACGTACGTACGTA CGCGCGCGCGCGCGCGCGCG (anti-sense strand)
  
  # Proteins bind flanking the cleavage site and cut in the "insert" sequence.
  
  
  
  CasCLOVER : Clo051 or another nickases with CRISPR-based binding domains
  
  ATCCT ATGCATGCATGCATGCATGC TTAACCGGTTAACCGG TACGTACGTACGTACGTACG CGGTC
    ||| ||||||||||||||||||||                              (sense strand)
    PAM    Target Sequence  \------->
                                     <-------\   Target Sequence   PAM
  (anti-sense strand)                         |||||||||||||||||||| |||
  TAGGA TACGTACGTACGTACGTACG AATTGGCCAATTGGCC ATGCATGCATGCATGCATGC GCCAG


Below are the example profiles:::
  
  Nuclease_Profiles :
  Cas9 :
    PAM : "NGG"
    PAM_Loc : "3p"
    PAM_Tol : 1
    Cut_Offset : -4
    Insert_size : FALSE

  Cpf1 :
    PAM : "TTTV"
    PAM_Loc : "5p"
    PAM_Tol : 1
    Cut_Offset : 26     #(Any where between 23 and 28)
    Insert_size : FALSE

  TALEN :
    PAM : FALSE
    PAM_Loc : FALSE
    PAM_Tol : 0
    Cut_Offset : Mid_insert
    Insert_size : "15:21"

  CasCLOVER :
    PAM : "NGG"
    PAM_Loc : "3p"
    PAM_Tol : 1
    Cut_Offset : Mid_insert
    Insert_size : "10:30"


Profile parameters
------------------

``PAM``
  protospacer adjacent motif - should be specified here and can contain 
  ambiguous nucleotides. 
  
``PAM_Loc`` 
  indicates the location of the PAM with respect to the pattern, either '5p', 
  '3p' or FALSE.
  
``PAM_Tol`` 
  indicates the tolerance for mismatches in the PAM sequence (ignorned if PAM 
  is FALSE). 
  
``Cut_Offset`` 
  indicates the offset from the 5' nucleotide of the PAM sequence where the 
  nuclease creates a double strand break, unless PAM is FALSE, then the 5' 
  position of the target sequence (also accepts "mid_insert" to specify middle 
  of region between paired alignments).
  
``Insert_size`` 
  is used if target sequences are expected to flank each other for editing, 
  such as with TALENs, and indicates the expected size of the insert. To input 
  a range, delimit the min and max by a colon, ie. 15:21. All names of 
  nucleases used to treat specimens need to have a profile. Additional profiles
  should be added under the 'Nuclease_Profiles' parameter.

