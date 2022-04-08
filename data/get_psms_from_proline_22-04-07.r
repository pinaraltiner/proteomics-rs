library("RPostgres")
#setting working directory
set_wd <- setwd("D:/dev/rust/proteomics-rs/data/benchmarks/")

project_id <- 481

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = paste0('msi_db_project_',project_id), 
  host = 'pippin', # i.e. 'ec2-54-83-201-96.compute-1.amazonaws.com'
  port = 5432, # or any other port specified by your DBA
  user = 'proline_db_user',
  password = 'proline'
)

### EXECUTE ONLY ONCE PER PROJECT ###
#dbSendQuery(con, "CREATE INDEX peptide_readable_ptm_string_pep_idx ON peptide_readable_ptm_string ( peptide_id ASC )")
#dbSendQuery(con, "CREATE INDEX peptide_readable_ptm_string_rs_idx ON peptide_readable_ptm_string ( result_set_id ASC )")


### QUERY TO GET THE FULL DATASET ###
# df <- dbGetQuery(con, statement = paste(
#   "SELECT count(*)",
#   "FROM peptide_match WHERE (result_set_id =",
#   paste(
#     sep = " OR result_set_id =",
#       470,
#       466,
#       472,
#       468,
#       480,
#       478,
#       484,
#       482,
#       490,
#       486,
#       492,
#       488,
#       496,
#       494,
#       502,
#       500,
#       504,
#       498,
#       528,
#       526
#     ),
#   ") AND score >= 13"
#   )
# )


peptide_df <- dbGetQuery(con, statement = paste(
  "SELECT id as peptide_id, sequence, calculated_mass ",
  "FROM peptide WHERE sequence IN(",
  paste(
    sep = ", ",
    "'ADENYYK'",
    "'AGGKPSQSPSQEAAGEAVLGAK'",
    "'AVGMPSPVSPK'",
    "'DKSPSSLLEDAK'",
    "'ESKSSPRPTAEK'",
    "'FGESDTENQNNK'",
    "'GGFFSSFMK'",
    "'GHLSEGLVTK'",
    "'GQEYLILEK'",
    "'IGEGTYGVVYK'",
    "'IKSYSFPK'",
    "'KTSPLNFK'",
    "'LMTGDTYTAHAGAK'",
    "'LQTVHSIPLTINK'",
    "'NFSAAKSLLNK'",
    "'NSFNNPAYYVLEGVPHQLLPPEPPSPAR'",
    "'RLSSFVTK'",
    "'SDGGVKPQSNNK'",
    "'SEGSPVLPHEPAK'",
    "'SESPPPLSDPK'",
    "'SNSTSSMSSGLPEQDR'",
    "'SQSTSEQEK'",
    "'SSSFREMDGQPER'",
    "'SSSFREMENQPHK'",
    "'SSSPTQYGLTK'",
    "'STLVLHDLLK'",
    "'TSPLNFK'",
    "'TVSTSSQPEENVDR'",
    "'VGSLTPPSSPK'",
    "'VKEEGYELPYNPATDDYAVPPPR'",
    "'VPASPLPGLER'",
    "'VQTTPPPAVQGQK'",
    "'WTAPESLAYNK'",
    "'YIEDEDYYK'",
    "'YMEDSTYYKASK'",
    "'YSLTVAVK'"
  ),
  ")"
)
)

### QUERY TO GET PEPTIDES OF INTEREST ###
psm_df <- dbGetQuery(con, statement = paste(
  "SELECT id, charge, peptide_id",
  "FROM peptide_match WHERE (peptide_id =",
  paste(peptide_df$peptide_id, collapse=" OR peptide_id="),
  ")"
)
)

ptm_df <- dbGetQuery(con, statement = paste(
  "SELECT DISTINCT peptide_id, readable_ptm_string AS modifications ",
  "FROM peptide_readable_ptm_string WHERE (peptide_id =",
  paste(peptide_df$peptide_id, collapse=" OR peptide_id="),
  ")"
)
)

spectrum_df <- dbGetQuery(
  con,
  statement = paste(
    "SELECT peptide_match.id AS peptide_match_id, spectrum.title AS spectrum_title, peaklist.raw_file_identifier AS peaklist_name",
    "FROM peptide_match, ms_query, spectrum, peaklist WHERE (peptide_match.peptide_id =",
    paste(peptide_df$peptide_id, collapse=" OR peptide_match.peptide_id="),
    ") AND peptide_match.ms_query_id = ms_query.id AND ms_query.spectrum_id = spectrum.id",
    "AND spectrum.peaklist_id = peaklist.id"
  )
)


#fix the merging system to get rid of multiplication of same psm match id for the different peptides
psm_peptide <- merge(x=psm_df, y=peptide_df, by="peptide_id", x.all=FALSE, y.all=FALSE)
psm_peptide_ptm <- merge(x=psm_peptide, y=ptm_df, by="peptide_id", x.all=FALSE, y.all=FALSE)
psm_peptide_ptm_spectrum <- merge(x=psm_peptide_ptm, y=spectrum_df, by.x ="id", by.y = "peptide_match_id",x.all=FALSE, y.all=FALSE)

write.table(psm_peptide_ptm_spectrum, file = "./eyers_phsopho_peptides_mascot.tsv", row.names = FALSE, sep="\t")

dbDisconnect(con)
rm(con)