class Params:
    def __init__(self):
        print("ver 1.01")
        self.folder = "/osc-fs_home/ramzan/ALS/"
        self.file_name = "als_meta.h5ad"
        self.file_path = self.folder + self.file_name
        self.deseq = True
        self.cell_type = "manual_anno_L0"
        self.gene_sets = [
        "Reactome_2022",
        "MSigDB_Hallmark_2020",
        "KEGG_2021_Human",
        "Jensen_DISEASES",
        "Human_Phenotype_Ontology",
        "GWAS_Catalog_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "GO_Biological_Process_2023",
        "Elsevier_Pathway_Collection",
        "BioPlanet_2019",
        "WikiPathway_2023_Human"
        ]