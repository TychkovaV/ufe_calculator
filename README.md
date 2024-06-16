## UFE calculator by Varvara Tychkova

## About

Originating from Pedruzzi and Rouzine article (https://doi.org/10.1371/journal.ppat.1009669), this python module provides convenient high-level interface for UFE metrics calculation using any reads (.fasta) provided.

Pipeline class takes all method parameters on init and calculate metrics on run(), saving all results including detected epistatic network.
Stages class allows user to start the pipeline from saved stage, saving time on repeating time-consuming stages (e.g. reads data parsing while proccessing same data with different parameters).

See main.py for example of module use
