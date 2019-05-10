########## Uncertain gene name -> Uniprot naming ############
#######  Manually currated cases (these proteins are TFs, their concurrents are not) #######
###  These records are OK ###
# BACH1 {:uniprot_ac=>"O14867", :uniprot_id=>"BACH1_HUMAN", :gene_names=>["BACH1"]}
# CUX1  {:uniprot_ac=>"P39880", :uniprot_id=>"CUX1_HUMAN", :gene_names=>["CUX1", "CUTL1"]}
# DDIT3 {:uniprot_ac=>"P35638", :uniprot_id=>"DDIT3_HUMAN", :gene_names=>["DDIT3", "CHOP", "CHOP10", "GADD153"]}
# ERG {:uniprot_ac=>"P11308", :uniprot_id=>"ERG_HUMAN", :gene_names=>["ERG"]}
# HBP1  {:uniprot_ac=>"O60381", :uniprot_id=>"HBP1_HUMAN", :gene_names=>["HBP1"]}
# HES1  {:uniprot_ac=>"Q14469", :uniprot_id=>"HES1_HUMAN", :gene_names=>["HES1", "BHLHB39", "HL", "HRY"]}
# MAFA  {:uniprot_ac=>"Q8NHW3", :uniprot_id=>"MAFA_HUMAN", :gene_names=>["MAFA"]}
# MGA {:uniprot_ac=>"Q8IWI9", :uniprot_id=>"MGAP_HUMAN", :gene_names=>["MGA", "KIAA0518", "MAD5"]}
# PDX1  {:uniprot_ac=>"P52945", :uniprot_id=>"PDX1_HUMAN", :gene_names=>["PDX1", "IPF1", "STF1"]}
# RAX {:uniprot_ac=>"Q9Y2V3", :uniprot_id=>"RX_HUMAN", :gene_names=>["RAX", "RX"]}
# SP1 {:uniprot_ac=>"P08047", :uniprot_id=>"SP1_HUMAN", :gene_names=>["SP1", "TSFP1"]}
# STAG1 {:uniprot_ac=>"Q8WVM7", :uniprot_id=>"STAG1_HUMAN", :gene_names=>["STAG1", "SA1", "SCC3"]}
###################
###  These records are NOT OK ###
# BACH1 {:uniprot_ac=>"Q9BX63", :uniprot_id=>"FANCJ_HUMAN", :gene_names=>["BRIP1", "BACH1", "FANCJ"]}
# CUX1  {:uniprot_ac=>"Q13948", :uniprot_id=>"CASP_HUMAN", :gene_names=>["CUX1", "CUTL1"]}
# DDIT3 {:uniprot_ac=>"P0DPQ6", :uniprot_id=>"DT3UO_HUMAN", :gene_names=>["DDIT3"]}
# ERG {:uniprot_ac=>"Q12809", :uniprot_id=>"KCNH2_HUMAN", :gene_names=>["KCNH2", "ERG", "ERG1", "HERG"]}
# HBP1  {:uniprot_ac=>"Q8IV16", :uniprot_id=>"HDBP1_HUMAN", :gene_names=>["GPIHBP1", "HBP1"]}
# HES1  {:uniprot_ac=>"A0A0B4J2D5", :uniprot_id=>"GAL3B_HUMAN", :gene_names=>["GATD3B", "HES1", "KNPI"]}
# MAFA  {:uniprot_ac=>"Q96E93", :uniprot_id=>"KLRG1_HUMAN", :gene_names=>["KLRG1", "CLEC15A", "MAFA", "MAFAL"]}
# MGA {:uniprot_ac=>"O43451", :uniprot_id=>"MGA_HUMAN", :gene_names=>["MGAM", "MGA", "MGAML"]}
# PDX1  {:uniprot_ac=>"O00330", :uniprot_id=>"ODPX_HUMAN", :gene_names=>["PDHX", "PDX1"]}
# RAX {:uniprot_ac=>"O75569", :uniprot_id=>"PRKRA_HUMAN", :gene_names=>["PRKRA", "PACT", "RAX", "HSD-14", "HSD14"]}
# SP1 {:uniprot_ac=>"Q8N907", :uniprot_id=>"DAND5_HUMAN", :gene_names=>["DAND5", "CER2", "CKTSF1B3", "GREM3", "SP1"]}
# STAG1 {:uniprot_ac=>"Q969W9", :uniprot_id=>"PMEPA_HUMAN", :gene_names=>["PMEPA1", "STAG1", "TMEPAI"]}

MANUALLY_CURRATED_UNIPROT_AC = {
  'BACH1' => "O14867",
  'CUX1'  => "P39880",
  'DDIT3' => "P35638",
  'ERG'   => "P11308",
  'HBP1'  => "O60381",
  'HES1'  => "Q14469",
  'MAFA'  => "Q8NHW3",
  'MGA'   => "Q8IWI9",
  'PDX1'  => "P52945",
  'RAX'   => "Q9Y2V3",
  'SP1'   => "P08047",
  'STAG1' => "Q8WVM7",
  # TFs from SELEX
  'CART1'   => 'Q15699', # it's ALX1_HUMAN (Q15699), not TRAF4_HUMAN (Q9BUZ4)
  'HINFP1'  => 'Q9BQA5', # HINFP_HUMAN
  'MSX3'    => 'P70354', # MSX3_MOUSE
  'RHOX11'  => 'Q810N8', # Q810N8_MOUSE, Reproductive homeobox 11
  'TCFAP2A' => 'P34056', # AP2A_MOUSE
  'ZFP652'  => 'Q5DU09', # ZN652_MOUSE
  'ZFP740'  => 'Q6NZQ6', # ZN740_MOUSE
}

#####################################
###  Undecidable cases. Here we took TFs whose first name matches query ###
###  These records are CONTROVERSIAL BUT TAKEN ###
# HMX3  {:uniprot_ac=>"A6NHT5", :uniprot_id=>"HMX3_HUMAN", :gene_names=>["HMX3", "NKX-5.1", "NKX5-1"]}
# KLF14 {:uniprot_ac=>"Q8TD94", :uniprot_id=>"KLF14_HUMAN", :gene_names=>["KLF14", "BTEB5"]}
# MED1  {:uniprot_ac=>"Q15648", :uniprot_id=>"MED1_HUMAN", :gene_names=>["MED1", "ARC205", "CRSP1", "CRSP200", "DRIP205", "DRIP230", "PBP", "PPARBP", "PPARGBP", "RB18A", "TRAP220", "TRIP2"]}
# MTF1  {:uniprot_ac=>"Q14872", :uniprot_id=>"MTF1_HUMAN", :gene_names=>["MTF1"]}
# NRF1  {:uniprot_ac=>"Q16656", :uniprot_id=>"NRF1_HUMAN", :gene_names=>["NRF1"]}
# TCF3  {:uniprot_ac=>"P15923", :uniprot_id=>"TFE2_HUMAN", :gene_names=>["TCF3", "BHLHB21", "E2A", "ITF1"]}
# TCF4  {:uniprot_ac=>"P15884", :uniprot_id=>"ITF2_HUMAN", :gene_names=>["TCF4", "BHLHB19", "ITF2", "SEF2"]}
# ZFP36 {:uniprot_ac=>"P26651", :uniprot_id=>"TTP_HUMAN", :gene_names=>["ZFP36", "G0S24", "NUP475", "RNF162A", "TIS11A", "TTP"]}
###################
###  These records are CONTROVERSIAL AND DISCARDED ###
# HMX3  {:uniprot_ac=>"Q8IZL8", :uniprot_id=>"PELP1_HUMAN", :gene_names=>["PELP1", "HMX3", "MNAR"]}
# KLF14 {:uniprot_ac=>"Q3SY56", :uniprot_id=>"SP6_HUMAN", :gene_names=>["SP6", "KLF14"]}
# MED1  {:uniprot_ac=>"O95243", :uniprot_id=>"MBD4_HUMAN", :gene_names=>["MBD4", "MED1"]}
# MTF1  {:uniprot_ac=>"Q01538", :uniprot_id=>"MYT1_HUMAN", :gene_names=>["MYT1", "KIAA0835", "KIAA1050", "MTF1", "MYTI", "PLPB1"]}
# NRF1  {:uniprot_ac=>"Q14494", :uniprot_id=>"NF2L1_HUMAN", :gene_names=>["NFE2L1", "HBZ17", "NRF1", "TCF11"]}
# TCF3  {:uniprot_ac=>"Q9HCS4", :uniprot_id=>"TF7L1_HUMAN", :gene_names=>["TCF7L1", "TCF3"]}
# TCF4  {:uniprot_ac=>"Q9NQB0", :uniprot_id=>"TF7L2_HUMAN", :gene_names=>["TCF7L2", "TCF4"]}
# ZFP36 {:uniprot_ac=>"P16415", :uniprot_id=>"ZN823_HUMAN", :gene_names=>["ZNF823", "ZFP36"]}

SEMI_CURRATED_UNIPROT_AC = {
  "HMX3"  => "A6NHT5",
  "KLF14" => "Q8TD94",
  "MED1"  => "Q15648",
  "MTF1"  => "Q14872",
  "NRF1"  => "Q16656",
  "TCF3"  => "P15923",
  "TCF4"  => "P15884",
  "ZFP36" => "P26651",
}
############# Unknown TF ###############
### TFs which can't be found in uniprot by name ###
CURRATED_MISNAMED_UNIPROT_AC = {
  'TAL1_SCL' => 'P17542', # TAL1_HUMAN (SCL is an alias for TAL1)
  'NCOR' => 'O75376', # NCOR1_HUMAN (there're several TFs in N-CoR complex, but more often than others N-CoR is denoted as N-CoR)
  'RXR' => 'P19793', # RXRA_HUMAN (it's hard to understand which subunit of RXR was used: alpha, beta or gamma; RXRA is taken arbitrarily)
}

#####################################
CURRATED_UNIPROT_AC = [
  MANUALLY_CURRATED_UNIPROT_AC,
  SEMI_CURRATED_UNIPROT_AC,
  CURRATED_MISNAMED_UNIPROT_AC
].reduce(&:merge)
