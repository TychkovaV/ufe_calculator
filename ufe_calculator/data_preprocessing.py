import pandas as pd
from .logger import Logger
import os.path
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import logomaker as lm

class GisaidFastaParser:
    def __init__(self, files:list[str], fasta_colnames:list[str], header_splitter:str, logger:Logger) -> None:
        self.files:list[str] = files
        self.headers:list[str] = fasta_colnames
        self.logger:Logger = logger
        self.fasta_spitter = header_splitter

    class _myseq:
        id = ""
        seq = ""
        def __init__(self,i,s) -> None:
            self.id = i
            self.seq = s           

    def parse(self) -> pd.DataFrame:
        def parse_file(file_handler):
            f = file_handler
            res = GisaidFastaParser._myseq("","")
            for line in f:
                if line[0] == ">":
                    if res.id != "":
                        yield GisaidFastaParser._myseq(res.id,res.seq)
                        res = GisaidFastaParser._myseq("","")
                    res.id = line[1:-1]
                else:
                    res.seq += line[:-1]
            yield GisaidFastaParser._myseq(res.id,res.seq)

        res = pd.DataFrame(dict(), columns = ["File"]+self.headers+["Seq"])
        self.logger.info("Parsing fasta files")
        for f in tqdm(self.files):
            with open(f) as fasta_file:
                for s in parse_file(fasta_file):
                    res.loc[len(res.index)] = [f,*s.id.split(self.fasta_spitter),str(s.seq)]
        return res


def get_fasta_files_in_folder(folder_path:str) -> list[str]:
    res = []
    extension=".fasta"
    for filename in os.listdir(folder_path):
        if filename.lower().endswith(extension):
            res.append(folder_path+"/"+filename)
    return res


class AaData:
    def __init__(self, logger) -> None:
        self.logger:Logger = logger
        self.data:None|pd.DataFrame = None
        self.aa_frequences:None|pd.DataFrame = None
        self.data_parser:None|GisaidFastaParser = None


    def load_from_folder(self, fasta_colnames:list[str], folder:str, header_splitter:str="|") -> None:
        if folder=="":
            self.logger.warning("Folder to load sequence data wasn't provided. Aborting")
            return  

        self.load_from_files(fasta_colnames, get_fasta_files_in_folder(folder), header_splitter)
    
    def load_from_files(self, fasta_colnames:list[str], files:list[str], header_splitter:str="|") -> None:
        if files==[]:
            self.logger.warning("No .fasta file were found to load. Aborting")
            return

        self.data_parser = GisaidFastaParser(files, fasta_colnames, header_splitter, self.logger)
        
        new_data = self.data_parser.parse()
        if self.data==None:
            self.data = new_data
        else:
            if not set(new_data.columns)<=set(self.data.columns):
                self.logger.error(f"Tried to use AAData.add() whith new data not having mandatory columns: \n\tneeded \n\t{self.data.columns}\n\tbut provided\n\t{new_data.columns}")
                return
            self.data = pd.concat([self.data, new_data[self.data.columns]])
            
    def drop_columns(self, to_drop:list[str]=[], to_keep:list[str]=[]) -> None:
        to_remove = list(set(self.data.columns)-set(to_keep)) if to_keep!=[] else to_drop
        if "Seq" in to_remove:
            self.logger.warning("Removing 'Seq' column from the data in AAData.drop_column(). Some methods will be unavailable from now on")
        self.data = self.data.drop(columns=to_remove)
        self.logger.info(f"Dropped columns \n\t{''.join(i+', ' for i in to_remove)[:-2]}\n\tkept\n\t{''.join(i+', ' for i in self.data.columns)[:-2]}")
    def drop_rows_cond(self, condition_fun) -> None:
        fun_res = self.data.apply(condition_fun, axis=1)
        self.data = self.data[fun_res != True]
        self.logger.info(f"Droped {sum(fun_res)} out of {len(fun_res)} rows")

    def value_counts(self)->dict[str,pd.Series]:
        return {name:self.data[name].value_counts() for name in self.data.columns}

    def add_column_based_on_row(self, title:str, value_calculator) -> None:
        if title in self.data.columns:
            self.logger.error(f"Column '{title}' already exists, AAData.add_column() failed")
            return
        self.data[title] = self.data.apply(value_calculator, axis=1)
        self.logger.info(f"Added '{title}' column")
    def add_column_based_on_other(self, title:str, reference_title:str, value_calculator) -> None:
        if title in self.data.columns:
            self.logger.error(f"Column '{title}' already exists, AAData.add_column() failed")
            return
        if reference_title not in self.data.columns:
            self.logger.error(f"Reference column '{reference_title}' doesn't exist, AAData.add_column() failed")
            return
        self.data[title] = self.data[reference_title].apply(value_calculator)
        self.logger.info(f"Added '{title}' column")
    def edit_column(self, column:str, series_fun):
        if column not in self.data.columns:
            self.logger.error(f"Column {column} does not exist, AAData.edit_column() failed")
            return
        self.data[column] = series_fun(self.data[column])

    def calculate_aa_frequences(self, weight_col:str="") -> None:
        protein_len:int = len(self.data['Seq'].iloc[0])
        alphabet:list[str] = list(set([i for a in self.data['Seq'] for i in a]))
        data_len:int = len(self.data.index)
        freq:list[dict[str,float]] = [{a:0. for a in alphabet} for i in range(protein_len)]

        if weight_col not in self.data.columns:
            self.logger.error(f"Column {weight_col} does not exist, AAData.calculate_aa_frequences() proceeds without weights")
            weight_col = ""
        weight_sum = data_len if weight_col=="" else sum(self.data[weight_col])
        delta_freq:list[float] = [1/weight_sum]*data_len if weight_col=="" else self.data[weight_col].apply(lambda x: x/weight_sum).values.flatten().tolist()
        row_idx = 0
        self.logger.info("Calculating sequence frequences")
        for s in tqdm(self.data['Seq']):
            for i,l in zip(range(protein_len),s):
                freq[i][l]+=delta_freq[row_idx]
            row_idx += 1
        freq_dict = {k: [dic[k] for dic in freq] for k in freq[0]}
        self.aa_frequences = pd.DataFrame(freq_dict)
        self.logger.info(f"Successefully calculated AA frequences")
    
    def drop_odd_lengths(self) -> None:
        mfl = (self.data['Seq'].str.len().value_counts().nlargest(1).index.tolist()[0]) # most frequent length
        filtered_data = self.data[self.data['Seq'].str.len() == mfl]
        self.logger.info(f"Most frequent length is {mfl}, filtering all other values, {round(len(filtered_data.index)*10000/len(self.data.index))/100}% sequences were kept")
        self.data = filtered_data
    def drop_non_polymorphic(self, cut:float=.05) -> None:
        if not isinstance(self.aa_frequences, pd.DataFrame):
            self.logger.error(f"AA frequences should be calculated before AAData.drop_non_polymorphic() usage. Aborting")
            return
        polymorphic_indexes = self.aa_frequences.max(axis = "columns") < 1-cut
        self.aa_frequences = self.aa_frequences[polymorphic_indexes]
        self.logger.info(f"Kept {sum(polymorphic_indexes)} sites out of {len(self.data['Seq'].iloc[0])}")
        self.data["Seq"] = self.data["Seq"].apply(lambda x: "".join(x[i] for i in polymorphic_indexes[polymorphic_indexes==True].index))
    def get_consensus(self) -> str:
        if not isinstance(self.aa_frequences, pd.DataFrame):
            self.logger.error("Can't get consensus sequence without AA frequences calculated")
            return ''
        consensus = ''
        freqs_ranked_tmp = self.aa_frequences.rank(axis = 1, ascending = False, method = 'max')
        for i in freqs_ranked_tmp.index:
            tmp = freqs_ranked_tmp.loc[i,:].values.flatten().tolist()
            consensus += freqs_ranked_tmp.columns[tmp.index(1.)]
        self.logger.info(f"Consensus sequence is {consensus}")
        return consensus
    def calculate_mutation_ratio(self) -> None:
        consensus = self.get_consensus()
        delta = 1/len(consensus)
        self.data["Mutation ratio"] = self.data["Seq"].apply(lambda str: sum(delta for (a, b) in zip(str, consensus) if a != b))

    def binarize(self) -> pd.DataFrame:
        consensus = self.get_consensus()
        res = pd.DataFrame(dict(), columns = ['Isolate name']+
                                             self.aa_frequences.index.tolist()+
                                             ['Weight']) 
        def seq_to_bin(s):
            res = []
            for i,l in zip(range(len(s)),s):
                if l == consensus[i]:
                    res.append(0)
                else:
                    res.append(1)
            return res
        self.logger.info("Binarizing sequence data")
        for i in tqdm(self.data.index):
            s = self.data.loc[i]
            bin_seq = seq_to_bin(s['Seq'])
            new_row = [s['Isolate name']]+bin_seq+[s['Weight']]
            res.loc[len(res.index)] = new_row
        return res
        

    def save_to_excel(self, title) -> None:
        if not isinstance(self.data, pd.DataFrame): 
            self.logger.warning("No data, unable to write to excel")    
        else:
            self.data.to_excel(title+" samples.xlsx")
            self.logger.info(f"AA data table saved to {title+' samples.xlsx'}")
       
    def plot_column_hist(self, colname:str, title:str, bins:int=20, group_by:str='') -> None:
        if colname not in self.data.columns:
            self.logger.error(f"Column {colname} does not exist, can't plot")
            return
        if group_by!='' and (group_by not in self.data.columns):
            self.logger.error(f"Trying to plot {colname} with grouping according to {group_by} but such column doesn't exist. Plotting without group_by")
            group_by=''

        if group_by=='':
            fig = plt.figure(figsize=(8, 24))
            ax = plt.subplot()
            self.data[colname].plot.hist(bins = bins, 
                                        ax = ax)
        else:
            fig, axs = plt.subplots(len(self.data[group_by].unique())+1,
                                    figsize=(8, 24), tight_layout=True,
                                    sharex = True, sharey = True)
            i = 0
            for val in self.data[group_by].unique():
                df = self.data[self.data[group_by] == val]
                df[colname].plot.hist(bins = bins, 
                                     ax = axs[i], 
                                     weights = np.ones_like(df[colname]) / len(df[colname]))
                axs[i].set_title(val)
                i += 1
            self.data[colname].plot.hist(bins = 20, 
                                        ax = axs[i], 
                                        weights = np.ones_like(self.data[colname]) / len(self.data[colname]))
            axs[i].set_title("Total")
        fig.savefig(title)

    def plot_column_kde(self, colname:str, title:str, group_by:str='') -> None:
        if colname not in self.data.columns:
            self.logger.error(f"Column {colname} does not exist, can't plot kde")
            return
        if group_by!='' and (group_by not in self.data.columns):
            self.logger.error(f"Trying to plot kde of {colname} with grouping according to {group_by} but such column doesn't exist. Plotting without group_by")
            group_by=''
        
        fig = plt.figure(figsize=(24, 12))
        ax = plt.subplot()
        if group_by=='':
            self.data[colname].plot.kde(ax = ax, 
                                       bw_method = 0.35)
        else:
            self.data.groupby(group_by).plot.kde(y=colname, 
                                                bw_method = 0.35,
                                                ax=ax, 
                                                legend = True,
                                                )
            ax.legend(self.data[group_by].unique())
        ax.set_xlim((0,ax.get_xlim()[1]))
        fig.savefig(title)
    
    def logoplot(self, title:str, width:int=17, reindex_fun=lambda x: x, figure=None) -> None:
        if not isinstance(self.aa_frequences, pd.DataFrame):
            self.logger.error("Can't draw logplot without AA frequences calculated")
            return
        count = len(self.aa_frequences.index)//width + (0 if len(self.aa_frequences.index)%width==0 else 1)
        if figure==None:
            fig, axs = plt.subplots(count, figsize=(np.ceil(12*width/17), np.ceil(3*count*width/17)), tight_layout=True)
        else:
            fig = figure
            axs = fig.subplots(count)
        for i in range(count):
            lm.Logo(self.aa_frequences.iloc[i*width:(i+1)*width].reset_index().iloc[: , 1:], 
                    color_scheme='dmslogo_funcgroup', 
                    vpad=.05,
                    ax=axs[i])

            sites_count=min(width, len(self.aa_frequences.index)-i*width)
            axs[i].set_xticks([i for i in range(sites_count)], minor=False)
            axs[i].set_xticklabels([reindex_fun(self.aa_frequences.index[j]) for j in range(i*width, i*width+sites_count)], 
                                   fontdict=None, 
                                   minor=False)
        if title!="":
            fig.savefig(title)

    def clear(self) -> None:
        self.logger.info("AAData cleared")
        self.data = None
        self.aa_frequences = None
        
    def save(self, data_title:str) -> None:
        if isinstance(self.data, pd.DataFrame):
            self.data.to_pickle(f'{data_title} data.pkl')

        if isinstance(self.aa_frequences, pd.DataFrame):
            self.aa_frequences.to_pickle(f'{data_title} freqs.pkl')

        self.logger.info(f"AAData saved as {data_title}")

    def load(self, data_title:str) -> None:
        if os.path.isfile(f'{data_title} data.pkl'):
            self.clear()
            self.data = pd.read_pickle(f'{data_title} data.pkl')
            self.logger.info(f"Loaded AAData for '{data_title}'")
        else:
            self.logger.warning(f"File '{data_title} data.pkl' doesn't exist, failed to load AAData")

        if os.path.isfile(f'{data_title} freqs.pkl'):
            self.aa_frequences = pd.read_pickle(f'{data_title} freqs.pkl')
            self.logger.info(f"Loaded frequences for AAData for '{data_title}'")
        else:
            self.logger.info(f"No frequences saved for '{data_title}'")