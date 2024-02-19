import pandas as pd
from params import Params
import seaborn as sns
import textwrap
import numpy as np
from utils.genes_ncbi_homo_sapiens_proteincoding import GENEID2NT as GeneID2nt_hs
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('macosx')

params = Params()

if params.deseq:
    res = pd.read_csv("/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/STAGE_2/deseq_results.tsv",
                      sep="\t", index_col=0)
    res = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
else:
    df = pd.read_csv('/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/R_scripts/diff_expr/output/all_cells'
                     '/glmQLFTest_result.tsv', sep="\t")
    res = df[(df['logCPM'] > 5) & (df['FDR'] < 0.05) & (abs(df['logFC']) > 0.5)]

res = res.index.tolist()
print(len(res))
# run one time to initialize
obo_fname = download_go_basic_obo()
fin_gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")
# run one time to initialize
mapper = {}

for key in GeneID2nt_hs:
    mapper[GeneID2nt_hs[key].Symbol] = GeneID2nt_hs[key].GeneID

inv_map = {v: k for k, v in mapper.items()}
# run one time to initialize

# Read NCBI's gene2go. Store annotations in a list of namedtuples
objanno = Gene2GoReader(fin_gene2go)
# Get namespace2association where:
#    namespace is:
#        BP: biological_process
#        MF: molecular_function
#        CC: cellular_component
#    assocation is a dict:
#        key: NCBI GeneID
#        value: A set of GO IDs associated with that gene
ns2assoc = objanno.get_ns2assc()
# run one time to initialize
goeaobj = GOEnrichmentStudyNS(
    GeneID2nt_hs.keys(),  # List of mouse protein-coding genes
    ns2assoc,  # geneid/GO associations
    obodag,  # Ontologies
    propagate_counts=False,
    alpha=0.05,  # default significance cut-off
    methods=['fdr_bh'])  # defult multipletest correction method
# run one time to initialize
GO_items = []

temp = goeaobj.ns2objgoea['BP'].assoc
for item in temp:
    GO_items += temp[item]

temp = goeaobj.ns2objgoea['CC'].assoc
for item in temp:
    GO_items += temp[item]

temp = goeaobj.ns2objgoea['MF'].assoc
for item in temp:
    GO_items += temp[item]


# pass list of gene symbols
def go_it(test_genes):
    print(f'input genes: {len(test_genes)}')

    mapped_genes = []
    for gene in test_genes:
        try:
            mapped_genes.append(mapper[gene])
        except:
            pass
    print(f'mapped genes: {len(mapped_genes)}')

    goea_results_all = goeaobj.run_study(mapped_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,
                                          x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),
                                          list(map(lambda y: inv_map[y], x.study_items)), ], goea_results_sig)),
                      columns=['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',
                               'n_study', 'n_go', 'study_genes'])

    GO = GO[GO.n_genes > 1]
    return GO


df = go_it(res)
df = df.sort_values(by='p_corr', ascending=True).head(100)
df['per'] = df.n_genes / df.n_go

fig, (ax, cax) = plt.subplots(figsize=(12, 12), ncols=2, gridspec_kw={'width_ratios': [20, 1], 'wspace': 0.05})
cmap = mpl.cm.bwr_r
norm = mpl.colors.Normalize(vmin=df.p_corr.min(), vmax=df.p_corr.max())
mapper = cm.ScalarMappable(norm=norm, cmap=cm.bwr_r)
sns.barplot(data=df, x='per', y='term', palette=mapper.to_rgba(df.p_corr.values), ax=ax)
ax.set_yticklabels([textwrap.fill(e, 30) for e in df['term']])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
cb.set_label('p_corr')
plt.subplots_adjust(left=0.4)
plt.tight_layout()
plt.savefig("GOATOOLS.png", dpi=300)