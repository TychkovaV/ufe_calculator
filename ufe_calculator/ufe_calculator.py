from .data_preprocessing import AaData
from .logger import Logger
from .method import MethodApplier
from matplotlib.figure import Figure
from enum import IntEnum
from typing import Callable, Any
from pandas import DataFrame, Series
from pathlib import Path



class Stage(IntEnum):
    NONE = 0
    LOAD_FROM_FASTA = 1
    EDIT_RAW_DATA = 2
    CALCULATE_FREQUNCES = 3
    CALCULATE_UFE = 4
    DONE = 5

    def __str__(self)->str:
        return self.name
    def name_lowercase(self)->str:
        return self.name.lower()
    
class ColumnAdder:
    logger:Logger
    colname:str
    reference_colname:str|None
    calculate_rowwise:bool
    formula:Callable[[DataFrame|Series],Any]
    def __init__(self, colname:str, formula:Callable[[DataFrame|Series],Any], 
                 logger:Logger,
                 reference_colname:str|None=None, calculate_rowwise=False) -> None:
        self.logger=logger
        self.colname=colname
        self.reference_colname=reference_colname
        self.calculate_rowwise=calculate_rowwise
        self.formula=formula
        if reference_colname==None and calculate_rowwise==False:
            self.logger.error("No reference column set for adding new column. Provide source column name or calculate values row wise")
        if colname=="":
            self.logger.error("Column name can't be empty string")
    def apply(self, aa_data:AaData) -> None:
        if self.calculate_rowwise:
            aa_data.add_column_based_on_row(self.colname, self.formula)
        else:
            aa_data.add_column_based_on_other(self.colname, self.reference_colname, self.formula)
            
class WeightAdder(ColumnAdder):
    def __init__(self, colname:str, categories_colname:str, categories_weights:dict[str,float], logger:Logger) -> None:
        self.logger=logger
        self.colname=colname
        self.calculate_rowwise=False
        self.reference_colname=categories_colname
        self.categories_weights=categories_weights
        self.formula=lambda _: self.logger.error("WeightAdder should never have its formula called directly. Use WeightAdder.apply() instead")
    def apply(self, aa_data: AaData) -> None:
        if self.categories_weights=={}:
            aa_data.add_column_based_on_row(self.colname, lambda _: 1)
        else:
            weights = {val : (self.categories_weights[val] if val in self.categories_weights.keys() else 0)/count 
                    for val, count in aa_data.value_counts()[self.reference_colname].items()}
            aa_data.add_column_based_on_other(self.colname, self.reference_colname, lambda x: weights[x])


class Pipeline :
    #### data ####
    aa:AaData
    ufe:MethodApplier
    logger:Logger
    stage:Stage
    run_untill:Stage
    #### parameters ####
    # common fields
    job_title:str
    job_saves_folder:str
    # parsing parameters
    fasta_colnames:list[str]
    input_folders:list[str]
    input_files:list[str]
    fasta_header_splitter:str
    # raw data editing parameters
    records_filters:list[Callable[[DataFrame],bool]]
    new_columns:list[ColumnAdder]
    columns_to_keep:list[str]
    columns_to_drop:list[str]
    # mutation frequence calculating parameters
    weight_colname:str
    sites_reindex_fun:Callable[[int],int]
    sites_filter_fun:Callable[[int],bool]
    # filtering, binarizing and UFE calculation parameters
    polymorphic_sites_cut:float
    # result plotting parameters
    sites_annotation:dict[str,tuple[int]]
    protein_len:int
    ufe_epinet_cut:float
    ufe0_epinet_cut:float
    ufe_plot_limits:tuple[float]


    def load(self, job_title, folder, stage:Stage) -> None:
        if stage == Stage.NONE:
            self.logger.error("Trying to load empty data (stage = NONE), aborted")
            return
        title = folder+"/"+job_title+" "+stage.name_lowercase()
        self.aa.load(title)
        self.ufe.load(title)

    def get_job_and_folder_prefix(self) -> str:
        return self.job_saves_folder+"/"+self.job_title

    def save(self) -> None:
        if self.stage == Stage.NONE:
            self.logger.error("Trying to save empty data (stage = NONE), aborted")
            return
        Path(self.job_saves_folder).mkdir(parents=True, exist_ok=True)
        title = self.get_job_and_folder_prefix()+" "+self.stage.name_lowercase()
        self.aa.save(title)
        self.ufe.save(title)

    def __init__(self, job_title:str="", job_saves_folder:str=".", logger:Logger=Logger(),
                fasta_colnames:list[str]=[], input_folders:list[str]=[], input_files:list[str]=[], fasta_header_splitter:str="|",
                records_filters:list[Callable[[DataFrame],bool]]=[], new_columns:list[ColumnAdder]=[], 
                columns_to_keep:list[str]=[], columns_to_drop:list[str]=[], 
                weight_colname:str="", sites_reindex_fun:Callable[[int],int]=lambda i: i, sites_filter_fun:Callable[[int],bool]=lambda i : True,
                polymorphic_sites_cut:float=.05,
                ufe_epinet_cut:float=.5, ufe0_epinet_cut:float=.5, ufe_plot_limits:tuple[float]=(-1.5,1.5),
                sites_annotation:dict[str,tuple[int]]={}, protein_len=-1,
                title_to_load:str="", folder_to_load_from:str=".", stage_to_load:Stage=Stage.NONE, run_untill_stage:Stage=Stage.DONE) -> None:
        self.job_title=job_title
        self.job_saves_folder=job_saves_folder
        self.logger=logger
        self.aa=AaData(self.logger)
        self.ufe=MethodApplier(self.logger)
        self.stage=stage_to_load
        self.run_untill=run_untill_stage
        if stage_to_load > Stage.NONE:
            self.load(title_to_load, folder_to_load_from, stage_to_load)
            self.stage=stage_to_load
        if self.stage <= Stage.NONE:
            self.fasta_colnames=fasta_colnames
            self.input_folders=input_folders
            self.input_files=input_files
            self.fasta_header_splitter=fasta_header_splitter
        if self.stage <= Stage.LOAD_FROM_FASTA:
            self.records_filters=records_filters
            self.new_columns=new_columns
            self.columns_to_keep=columns_to_keep
            self.columns_to_drop=columns_to_drop
        if self.stage <= Stage.EDIT_RAW_DATA:
            self.weight_colname=weight_colname
        if self.stage <= Stage.CALCULATE_FREQUNCES:
            self.polymorphic_sites_cut=polymorphic_sites_cut
        if self.stage <= Stage.CALCULATE_UFE:
            ...
        if self.stage <= Stage.DONE:
            self.sites_reindex_fun=sites_reindex_fun
            self.sites_filter_fun=sites_filter_fun
            self.sites_annotation=sites_annotation
            self.protein_len=protein_len
            self.ufe_epinet_cut=ufe_epinet_cut
            self.ufe0_epinet_cut=ufe0_epinet_cut
            self.ufe_plot_limits=ufe_plot_limits

        
    def run(self):
        if self.stage <= Stage.NONE and self.stage < self.run_untill:
            for f in self.input_folders:
                self.aa.load_from_folder(self.fasta_colnames, f, self.fasta_header_splitter)
            if self.input_files!=[]:
                self.aa.load_from_files(self.fasta_colnames, self.input_files, self.fasta_header_splitter)
            self.aa.drop_odd_lengths()
            self.stage=Stage.LOAD_FROM_FASTA
            self.save()

        if self.stage <= Stage.LOAD_FROM_FASTA and self.stage < self.run_untill:
            for filter in self.records_filters:
                self.aa.drop_rows_cond(filter)
            for column_adder in self.new_columns:
                column_adder.apply(self.aa)
            self.aa.drop_columns(to_keep=self.columns_to_keep, to_drop=self.columns_to_drop)
            self.aa.save_to_excel(self.get_job_and_folder_prefix())
            self.stage=Stage.EDIT_RAW_DATA
            self.save()

        if self.stage <= Stage.EDIT_RAW_DATA and self.stage < self.run_untill:
            self.logger.info(f"Filtering out rows with {self.weight_colname} <= 0")
            self.aa.drop_rows_cond(lambda df_row: df_row[self.weight_colname] <= 0)
            self.aa.calculate_aa_frequences(weight_col=self.weight_colname)
            self.aa.logoplot(self.get_job_and_folder_prefix()+" logo full.png", reindex_fun=self.sites_reindex_fun)

            self.aa.drop_non_polymorphic(cut=self.polymorphic_sites_cut)
            self.aa.calculate_mutation_ratio()

            self.aa.plot_column_hist("Mutation ratio", self.get_job_and_folder_prefix()+" mutation rate hist.png")
            self.aa.plot_column_kde("Mutation ratio", self.get_job_and_folder_prefix()+" mutation rate kde.png")
            self.aa.logoplot(self.get_job_and_folder_prefix()+" logo polymorphic.png", reindex_fun=self.sites_reindex_fun)

            self.stage=Stage.CALCULATE_FREQUNCES
            self.save()

        if self.stage <= Stage.CALCULATE_FREQUNCES and self.stage < self.run_untill:           
            data = self.aa.binarize()
            self.ufe.set(data)
            self.ufe.calculate_ufe()
            self.ufe.calculate_ufe0()

            self.stage=Stage.CALCULATE_UFE
            self.save()

        if self.stage <= Stage.CALCULATE_UFE and self.stage < self.run_untill:
            self.ufe.plot_circular_graph(self.get_job_and_folder_prefix()+" circular graph.png",
                                         sites_annotation=self.sites_annotation, 
                                         sites_count=self.protein_len,
                                         ufe0_cut=self.ufe0_epinet_cut,
                                         label_filter_fun=self.sites_filter_fun,
                                         label_rename_fun=self.sites_reindex_fun)
            self.ufe.save_to_excel(self.get_job_and_folder_prefix())
            self.ufe.plot_epistatic_net(self.get_job_and_folder_prefix()+" epistatic net.png", 
                                        ufe_cut=self.ufe_epinet_cut, ufe0_cut=self.ufe0_epinet_cut, plot_lims=self.ufe_plot_limits,
                                        label_filter_fun=self.sites_filter_fun, 
                                        label_rename_fun=self.sites_reindex_fun,
                                        bottom_left_filler=lambda x: (self.aa.logoplot("", reindex_fun=self.sites_reindex_fun, figure=x), 
                                                                     Figure.suptitle(x, 'logo plot for proccessed sites')))