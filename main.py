import ufe_calculator as u
FASTA_HEADER = "Isolate name | Isolate ID | Type | Passage details/history | Lineage | Clade | Collection date | Submitter | Sample ID by sample provider | Sample ID by submitting lab | Last modified | Originating lab | Submitting lab | Gene name | Protein Accession no. | Protein INSDC"
clades_to_use = ["6B.1A.1", 
                "6B.1A.5b",
                "6B.1A.3",
                "6B.1A.6",
                "6B.1A.5a",
                "6B.1A.2",
                "6B.1A.7",
                "6B.1A.5a.1",
                "6B.1A.5a.2a",
                "6B.1A.5a.2a.1"]

logger = u.Logger()
number_of_first_sites_to_drop=17
reindex_fun=lambda x: x-number_of_first_sites_to_drop+1
filter_fun=lambda x: reindex_fun(x)>0
data_filter = lambda df: df["Clade"] not in clades_to_use

columns_raw = [u.ColumnAdder("Year", lambda df_row: df_row['File'][-10:-6], logger, calculate_rowwise=True),
               u.WeightAdder("Weight", "", {}, logger)]
columns_clade_based = [u.ColumnAdder("Year", lambda df_row: df_row['File'][-10:-6], logger, calculate_rowwise=True),
                       u.WeightAdder("Weight", "Clade", {clade:1 for clade in clades_to_use}, logger)]
columns_year_based = [u.ColumnAdder("Year", lambda df_row: df_row['File'][-10:-6], logger, calculate_rowwise=True),
                      u.WeightAdder("Weight", "Year", {str(year):1 for year in range(2016,2024)}, logger)]

columns=[columns_raw, columns_clade_based, columns_year_based]
titles=["no_weighting", "clade_averaging", "year_averaging"]

for i, title in enumerate(titles):
    calculator = u.Pipeline(job_title=title, job_saves_folder="./out", logger=logger, 
                            fasta_colnames=FASTA_HEADER.split(" | "), input_folders=['./H1N1_2016_2023'], input_files=[], fasta_header_splitter="|",
                            records_filters=[data_filter], new_columns=columns[i], columns_to_keep=["File", "Year", "Clade", "Isolate name", "Seq", "Weight"], columns_to_drop=[],
                            weight_colname="Weight", sites_reindex_fun=reindex_fun, sites_filter_fun=filter_fun, polymorphic_sites_cut=.05,
                            ufe_epinet_cut=.5, ufe0_epinet_cut=.5, ufe_plot_limits=(-1.5,1.5),
                            sites_annotation={}, protein_len=566-number_of_first_sites_to_drop+1,
                            title_to_load="", folder_to_load_from="", stage_to_load=u.Stage.NONE, run_untill_stage=u.Stage.DONE)
    calculator.run()

