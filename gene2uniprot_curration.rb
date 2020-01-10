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
  # TFs from CISBP
  'OSR1'    => 'Q8TAX0', # OSR1_HUMAN
  'ZNF225'  => 'Q9UK10', # ZN225_HUMAN
  # TFs from mouse uniprobe
  "Cart1"   => "Q8C8B0", # It's ALX1_MOUSE, not TRAF4_MOUSE
  "Elf1"    =>"Q60775",  # It's ELF1_MOUSE, not EFNA2_MOUSE
  "Elf2"    => "Q9JHC9", # It's ELF2_MOUSE, not EFNB2_MOUSE
  "Erg"     => "P81270", # It's ERG_MOUSE, not KCNH2_MOUSE
  "Fli1"    => "P26323", # It's FLI1_MOUSE, not FLII_MOUSE
  "Hbp1"    => "Q8R316", # It's HBP1_MOUSE, not HDBP1_MOUSE
  "Mrg2"    => "P97368", # It's MEIS3_MOUSE, not CITE4_MOUSE
  "Osr1"    => "Q9WVG7", # It's OSR1_MOUSE, not OXSR1_MOUSE
  "Rax"     => "O35602", # It's RX_MOUSE, not PRKRA_MOUSE
}

#####################################
###  Undecidable cases. Here we took TFs whose first name matches query (or there're other reasons to prefer one TF to another) ###
###  These records are CONTROVERSIAL BUT TAKEN ###
# HMX3  {:uniprot_ac=>"A6NHT5", :uniprot_id=>"HMX3_HUMAN", :gene_names=>["HMX3", "NKX-5.1", "NKX5-1"]}
# KLF14 {:uniprot_ac=>"Q8TD94", :uniprot_id=>"KLF14_HUMAN", :gene_names=>["KLF14", "BTEB5"]}
# MED1  {:uniprot_ac=>"Q15648", :uniprot_id=>"MED1_HUMAN", :gene_names=>["MED1", "ARC205", "CRSP1", "CRSP200", "DRIP205", "DRIP230", "PBP", "PPARBP", "PPARGBP", "RB18A", "TRAP220", "TRIP2"]}
# MTF1  {:uniprot_ac=>"Q14872", :uniprot_id=>"MTF1_HUMAN", :gene_names=>["MTF1"]}
# NRF1  {:uniprot_ac=>"Q16656", :uniprot_id=>"NRF1_HUMAN", :gene_names=>["NRF1"]}
# TCF3  {:uniprot_ac=>"P15923", :uniprot_id=>"TFE2_HUMAN", :gene_names=>["TCF3", "BHLHB21", "E2A", "ITF1"]}
# TCF4  {:uniprot_ac=>"P15884", :uniprot_id=>"ITF2_HUMAN", :gene_names=>["TCF4", "BHLHB19", "ITF2", "SEF2"]}
# ZFP36 {:uniprot_ac=>"P26651", :uniprot_id=>"TTP_HUMAN", :gene_names=>["ZFP36", "G0S24", "NUP475", "RNF162A", "TIS11A", "TTP"]}
# SOX12 {:uniprot_ac=>"O15370", :uniprot_id=>"SOX12_HUMAN", :gene_names=>["SOX12", "SOX22"]}
# Cutl1 {:uniprot_ac=>"P53564", :uniprot_id=>"CUX1_MOUSE", :gene_names=>["Cux1", "Cutl1", "Cux", "Kiaa4047"]}
# Irx6  {:uniprot_ac=>"Q9ER75", :uniprot_id=>"IRX6_MOUSE", :gene_names=>["Irx6", "Irxb3"]}
# Mrg1  {:uniprot_ac=>"P97367", :uniprot_id=>"MEIS2_MOUSE", :gene_names=>["Meis2", "Mrg1", "Stra10"]}
# Sox21 {:uniprot_ac=>"Q811W0", :uniprot_id=>"SOX21_MOUSE", :gene_names=>["Sox21"]}
# Tcf1  {:uniprot_ac=>"Q00417", :uniprot_id=>"TCF7_MOUSE", :gene_names=>["Tcf7", "Tcf-1", "Tcf1"]}
# Tcf3  {:uniprot_ac=>"P15806", :uniprot_id=>"TFE2_MOUSE", :gene_names=>["Tcf3", "Alf2", "Me2", "Tcfe2a"]}
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
# SOX12 {:uniprot_ac=>"O60248", :uniprot_id=>"SOX15_HUMAN", :gene_names=>["SOX15", "SOX12", "SOX20", "SOX26", "SOX27"]}
# Cutl1 {:uniprot_ac=>"P70403", :uniprot_id=>"CASP_MOUSE", :gene_names=>["Cux1", "Cutl1"]}
# Irx6  {:uniprot_ac=>"P81066", :uniprot_id=>"IRX2_MOUSE", :gene_names=>["Irx2", "Irx6", "Irxa2"]}
# Mrg1  {:uniprot_ac=>"O35740", :uniprot_id=>"CITE2_MOUSE", :gene_names=>["Cited2", "Mrg1", "Msg2"]}
# Sox21 {:uniprot_ac=>"Q04888", :uniprot_id=>"SOX10_MOUSE", :gene_names=>["Sox10", "Sox-10", "Sox21"]}
# Tcf1  {:uniprot_ac=>"P22361", :uniprot_id=>"HNF1A_MOUSE", :gene_names=>["Hnf1a", "Hnf-1", "Hnf-1a", "Tcf1"]}
# Tcf3  {:uniprot_ac=>"Q9Z1J1", :uniprot_id=>"TF7L1_MOUSE", :gene_names=>["Tcf7l1", "Tcf3"]}

SEMI_CURRATED_UNIPROT_AC = {
  "HMX3"  => "A6NHT5",
  "KLF14" => "Q8TD94",
  "MED1"  => "Q15648",
  "MTF1"  => "Q14872",
  "NRF1"  => "Q16656",
  "TCF3"  => "P15923",
  "TCF4"  => "P15884",
  "ZFP36" => "P26651",
  "SOX12" => "O15370",
  "Cutl1" =>"P53564",
  "Irx6"  =>"Q9ER75",
  "Mrg1"  =>"P97367",
  "Sox21" =>"Q811W0",
  "Tcf1"  =>"Q00417",
  "Tcf3"  =>"P15806",
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
