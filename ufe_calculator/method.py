import pandas as pd
import os.path
from .logger import Logger
import numpy as np
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from pycirclize import Circos



class MethodApplier:
    def __init__(self, logger:Logger) -> None:
        self.logger:Logger = logger
        self.data:None|pd.DataFrame = None
        self.ufe:None|pd.DataFrame = None
        self.ufe0:None|pd.DataFrame = None

    def set(self, binarized_data:pd.DataFrame) -> None:
        self.data = binarized_data.copy()

    def clear(self) -> None:
        self.logger.info("UFE data cleared")
        self.data = None
        self.ufe = None
        self.ufe0 = None
        
    def save(self, data_title:str) -> None:
        if isinstance(self.data, pd.DataFrame):
            self.data.to_pickle(f'{data_title} bin_data.pkl')
        if isinstance(self.ufe, pd.DataFrame):
            self.ufe.to_pickle(f'{data_title} ufe_data.pkl')
        if isinstance(self.ufe0, pd.DataFrame):
            self.ufe0.to_pickle(f'{data_title} ufe0_data.pkl')

        self.logger.info(f"UFE data saved as {data_title}")

    def load(self, data_title:str) -> None:
        if os.path.isfile(f'{data_title} bin_data.pkl'):
            self.clear()
            self.data = pd.read_pickle(f'{data_title} bin_data.pkl')
            self.logger.info(f"Loaded binarized sequence data for '{data_title}'")
        else:
            self.logger.warning(f"File '{data_title} bin_data.pkl' doesn't exist, failed to load binarized sequence data")
        
        if os.path.isfile(f'{data_title} ufe_data.pkl'):
            self.ufe = pd.read_pickle(f'{data_title} ufe_data.pkl')
            self.logger.info(f"Loaded ufe data for for '{data_title}'")
        else:
            self.logger.info(f"No ufe data saved for '{data_title}'")

        if os.path.isfile(f'{data_title} ufe0_data.pkl'):
            self.ufe0 = pd.read_pickle(f'{data_title} ufe0_data.pkl')
            self.logger.info(f"Loaded ufe0 data for for '{data_title}'")
        else:
            self.logger.info(f"No ufe0 data saved for '{data_title}'")

    def ufe_ij(self, i:int, j:int, min_freq_to_calc:float) -> float:
        f11 = self.data["Weight"][(self.data[i] == 1) & (self.data[j] == 1)].sum()
        f10 = self.data["Weight"][(self.data[i] == 1) & (self.data[j] == 0)].sum()
        f01 = self.data["Weight"][(self.data[i] == 0) & (self.data[j] == 1)].sum()
        f00 = self.data["Weight"][(self.data[i] == 0) & (self.data[j] == 0)].sum()
        tmp =  (1-(np.log(f11/f00)/np.log(f01*f10/f00/f00))) if (min(f11,f00,f10,f01)/sum([f11,f00,f10,f01]) > min_freq_to_calc) else np.nan
        return tmp
    
    def ufe_ij0(self, i:int, j:int, min_freq_to_calc:float) -> float:
        tmp = np.nan
        for k in self.data.columns[1:-1]:
            if k != i and k != j:
                f11 = self.data["Weight"][(self.data[i] == 1) & (self.data[j] == 1) & (self.data[k] == 0)].sum()
                f10 = self.data["Weight"][(self.data[i] == 1) & (self.data[j] == 0) & (self.data[k] == 0)].sum()
                f01 = self.data["Weight"][(self.data[i] == 0) & (self.data[j] == 1) & (self.data[k] == 0)].sum()
                f00 = self.data["Weight"][(self.data[i] == 0) & (self.data[j] == 0) & (self.data[k] == 0)].sum()
                if sum([f11,f00,f10,f01]) == 0:
                    #self.logger.error("ASDASDASD internal error if UFE0")
                    continue
                elif min(f11,f00,f10,f01)/sum([f11,f00,f10,f01]) < min_freq_to_calc:
                    continue
                elif np.isnan(tmp) or (abs(1-(np.log(f11/f00)/np.log(f01*f10/f00/f00))) < abs(tmp)):
                    tmp = (1-(np.log(f11/f00)/np.log(f01*f10/f00/f00)))
        return tmp

    def calculate_ufe(self, min_freq_to_calc:float=.035) -> None:
        if not isinstance(self.data, pd.DataFrame):
            self.logger.error("No data set, unable to calculate UFE. Stopping")
        sites = self.data.columns[1:-1]
        self.ufe = pd.DataFrame(index = sites, columns = sites)
        self.logger.info("Calculating UFE")
        for i in tqdm(sites):
            for j in sites:
                if i<j:
                    tmp = self.ufe_ij(i, j, min_freq_to_calc)
                    self.ufe.loc[i,j] = tmp
                    self.ufe.loc[j,i] = tmp
    
    def calculate_ufe0(self, min_freq_to_calc:float=.035) -> None:
        if not isinstance(self.data, pd.DataFrame):
            self.logger.error("No data set, unable to calculate UFE0. Stopping")
        if not isinstance(self.ufe, pd.DataFrame):
            self.logger.error("No ufe set, unable to calculate UFE0. Stopping")
        sites = self.data.columns[1:-1]
        self.ufe0 = pd.DataFrame(index = sites, columns = sites)
        self.logger.info("Calculating UFE0")
        for i in tqdm(sites):
            for j in sites:
                if i<j:
                    tmp = self.ufe_ij0(i, j, min_freq_to_calc)
                    self.ufe0.loc[i,j] = tmp
                    self.ufe0.loc[j,i] = tmp
        self.ufe0[pd.isna(self.ufe)]=np.nan
    
    def save_to_excel(self, title) -> None:
        if not isinstance(self.ufe, pd.DataFrame): 
            self.logger.warning("No UFE calculated, unable to write to excel")    
        else:
            self.ufe.to_excel(title+"_ufe.xlsx")
            self.logger.info(f"UFE table saved to {title+'_ufe.xlsx'}")
        if not isinstance(self.ufe0, pd.DataFrame): 
            self.logger.warning("No UFE0 calculated, unable to write to excel")    
        else:
            self.ufe0.to_excel(title+"_ufe0.xlsx")
            self.logger.info(f"UFE0 table saved to {title+'_ufe0.xlsx'}")

    def plot_circular_graph(self, title:str,
                           sites_annotation:dict[str,tuple[int]]={},
                           sites_count:int=-1,
                           ufe0_cut:float=.5,
                           label_filter_fun=lambda x: True, 
                           label_rename_fun=lambda x: x) -> None:
        if sites_count==-1 and sites_annotation=={}:
            self.logger.error("Provide either sites count or sites annotation to plot circular plot")
            return
        if sites_annotation!={}:
            ...
            # TODO: implement sectors split according to functional parts of protein
            return
        else:
            circos = Circos(sectors={"protein": (1,sites_count+1)})
            sector = circos.sectors[0]
            track = sector.add_track(r_lim=(100, 100))
            minor_xticks = np.arange(1, sites_count+1, 1)
            track.xticks(minor_xticks, outer=False, tick_length=1.2)
            major_xticks = np.arange(10, sites_count+1, 10)
            track.xticks(major_xticks, outer=False, show_bottom_line=True)
            for x in major_xticks:
                track.text(str(x), x=x, r=105, size=10, adjust_rotation=False)

            labels_to_use = [label_filter_fun(c) for c in self.ufe0.columns]
            sites = self.ufe0.columns[labels_to_use]

            connections = set()
            
            

            for i in sites:
                for j in sites:
                    if abs(self.ufe0[i][j])>ufe0_cut:
                        connections.add(label_rename_fun(i))
                        connections.add(label_rename_fun(j))
                        circos.link_line(("protein", label_rename_fun(i)), ("protein", label_rename_fun(j)), color="green" if self.ufe0[i][j]>0 else "red")
            for x in connections:
                track.text(str(x), x=x, r=96, size=5, adjust_rotation=False)

            circos.plotfig().savefig(title)





    def plot_epistatic_net(self, title:str, 
                           ufe_cut:float=.5, ufe0_cut:float=.5, 
                           plot_lims = None, 
                           label_filter_fun=lambda x: True, 
                           label_rename_fun=lambda x: x, 
                           bottom_left_filler=None,
                           annotate_standouts=False) -> None:

        def plot_ufe(ufe:pd.DataFrame, cut:float, ax:Axes):
            labels_to_use = [label_filter_fun(c) for c in ufe.columns]
            table = ufe.loc[labels_to_use, labels_to_use]
            
            def pack_index(df:pd.DataFrame, label_rename_fun):
                return {i:j for i,j in zip(range(len(df.columns)), 
                                           [label_rename_fun(label) for label in df.index.values.flatten().tolist()])}
            adjacency_matrix = table.to_numpy()
            mylabels = pack_index(table, label_rename_fun)

            rows, cols = np.where(np.logical_not(pd.isna(adjacency_matrix))) and np.where(abs(adjacency_matrix) > cut)
            edges = zip(rows.tolist(), cols.tolist())
            gr = nx.Graph()
            all_rows = np.unique(rows)
            mylabels = {l:mylabels[l] for l in all_rows}
            for n in all_rows:
                gr.add_node(n)
            for e in edges:
                gr.add_edge(*e,
                            weight=adjacency_matrix[e[0]][e[1]], 
                            color = "r" if adjacency_matrix[e[0]][e[1]] < 0 else "g", 
                            label = round(adjacency_matrix[e[0]][e[1]],2))
            colors = nx.get_edge_attributes(gr,'color').values()
            weights = [abs(i) for i in nx.get_edge_attributes(gr,'weight').values()]
            
            pos = nx.nx_agraph.graphviz_layout(gr)    
                
            nx.draw(gr, pos,
                    ax = ax,
                    node_size=220, labels=mylabels, with_labels=True, font_size = 8,
                    width=[min(1.5*(i-.2),2) for i in weights],
                    edge_color=colors)
            nx.draw_networkx_edge_labels(gr, pos, ax=ax,
                                        edge_labels=nx.get_edge_attributes(gr, "label"),
                                        font_color='black', font_size=7
                                        )

        fig = plt.figure(figsize=(20, 20))    
        subfigs = fig.subfigures(2,2)#, figsize=(20, 20)), tight_layout=True)
        #fig, axs = plt.subplots(2,2, figsize=(20, 20), tight_layout=True)
        axs = [[(subfigs[i][j].subplots() if not (i==1 and j ==0) else None) for j in range(2)] for i in range(2)]
        plot_ufe(self.ufe, ufe_cut, axs[1][1])
        subfigs[1][1].suptitle('UFE ij epistatic network')
        #axs[1][1].set_title("UFE ij epistatic network")
        plot_ufe(self.ufe0, ufe0_cut, axs[0][0])
        subfigs[0][0].suptitle('UFE ij0 epistatic network')
        #axs[0][0].set_title("UFE ij0 epistatic network")

        xs = self.ufe.values.flatten()
        ys = self.ufe0.values.flatten()
        axs[0][1].scatter(xs, 
                       ys,
                       marker=".")
        dim = len(self.ufe)
        if annotate_standouts:
            for i in range(dim):
                for j in range(dim):
                    if i<j and \
                    not np.isnan(xs[i*dim+j]) and \
                    not np.isnan(ys[i*dim+j]) and \
                    ys[i*dim+j] > ufe0_cut: 
                        axs[0][1].annotate(f"{label_rename_fun(self.ufe.columns[i])}-{label_rename_fun(self.ufe.columns[j])}",
                                        (xs[i*dim+j]*1.03-.05, ys[i*dim+j]*1.02))
        axs[0][1].set_aspect('equal', adjustable='box')
        axs[0][1].set_xlabel("UFE ij")
        axs[0][1].set_ylabel("UFE ij0")
        if plot_lims != None:
            axs[0][1].set_ylim(plot_lims)
            axs[0][1].set_xlim(plot_lims)
        else:
            axs[0][1].set_ylim(axs[0][1].get_xlim())
            axs[0][1].set_xlim(axs[0][1].get_xlim())
        axs[0][1].grid(axis='both', which = 'both')
        subfigs[0][1].suptitle('UFE metrics for each site pair')        
        #axs[0][1].set_title("UFE metrics for each site pair")
        axs[0][1].hlines([ufe0_cut, -ufe0_cut], *axs[0][1].get_xlim(), colors='black', linestyles='dashed', linewidth=.7).set_dashes((0, (7,10)))
        axs[0][1].vlines([ufe_cut, -ufe_cut], *axs[0][1].get_ylim(), colors='black', linestyles='dashed', linewidth=.7).set_dashes((0, (7,10)))

        if bottom_left_filler!=None:
            bottom_left_filler(subfigs[1][0])       
            
        fig.savefig(title)

