

### Import libraries
import sys, getopt


### Load arguments
def main(argv):
    ### DataSettings options
    inputType = '' # [-i, required, no default] Nucleotide (NC), Nucleotide (C), Protein, Distance Matrix (.meg), Tree (newick) 
    MBS = '?' # [-M, default ?] Missing Base Symbol 
    IBS = '.' # [-I, default .] Identical Base Symbol 
    GS = '-' # [-G, default -] Gap Symbol
    LS = 'All Sites' # [-L, default All Sites] Labelled sites
    LtI = '' # [-Z] Labels to include 
    ### ProcessType options
    ppXX = 'Infer' # [-X, default Infer] infer tree
    ppYY = 'NJ' # [-Y, default NJ] Type of analysis
    ### AnalysisSettings    
    Analysis = 'Phylogeny Reconstruction' # [-A, default Phylogeny Reconstruction] Type of analysis
    Scope = 'All Selected Taxa' # [-W, default All Selected Taxa]
    StatisticalMethod = 'Neighbor-joining' # [-P, default NJ] default Type of phylogeny 
    TestOfPhy = 'Bootstrap method' # [-B, default boostrap] Statistical test used
    BootstrapReplicates = 100 # [-R, default 100] Number of bs replicates 
    SubstType = 'Nucleotide' # [-s, default Nucleotide] Type of substitutions to include 
    Model = '' # [-m, required no default] Evolutionary model to use on the distance matrix
    TiTvRatio = 'Not Applicable' # [-t, default NA] Transition/Transversion Ratio for specific models 
    SubstToInclude = 'All' # [-a, default all] default = All
    RatesAmongSites = 'Uniform Rates' # [-r, default uniform] Substituion rates among sites 
    GammaP = 'Not Applicable' # [-g, default NA] Gamma parameter if rates sampled from a gamma distribution 
    PatternAmongLineages = 'Same (Homogeneous)' # [-P, default Same] Pattern among lineages 
    GapsNAdata = 'Pairwise deletion' # [-D, default Pairwise]
    SiteCovCutoff = 'Not Applicable' # [-C, default NA] Coverage cutoff
    CodonPos = '1st, 2nd, 3rd, Non-Coding' # [-c, default 123NC]
    Threads = '1' # [-T, default 1] change if bootstrap = true 
    GeneticCodeTable = 'Not Applicable' # [-g, default NA]
    GeneticCode = 'Not Applicable' # [-j, default NA]
    TimeLimit = 'False' # [w, default False]
    MaximumExecTime = '-1' # [m, default -1]
    arg_list = sys.argv[1:]
    short_options = 'i:M:I:G:L:Z:X:Y:A:W:P:B:R:s:m:t:a:r:g:P:D:C:c:T:g:j:w:m:'
    long_options = ['inputType=', 'MBS=', 'IBS=', 'GS=', 'LS=', 'LtI=', 'ppInfer=', 'ppXX=', 'Analysis=', 'Scope=', 'StatisticalMethod=', 'TestOfPhy=', 'BootstrapReplicates=', 'SubstType=', 'RatesAmongSites=', 'GammaP=', 'PatternAmongLineages=', 'GapsNAdata=', 'SiteCovCutoff=', 'CodonPos=', 'Threads=', 'GeneticCodeTable=', 'GeneticCode=', 'TimeLimit=', 'MaximumExecTime=']
    usage = ''

    try:
        options, args = getopt.getopt(arg_list,short_options,long_options)
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)

    if len(arg_list)<1 or not options:
        print(usage)
        sys.exit(2)

    for arg, val in options:
        if arg in ('-i', '--infile'):
            inputType = val
        elif arg in ('-M', '--MBS'):
            MBS = val
        elif arg in ('-I', '--IBS'):
            IBS = val
        elif arg in ('-G', '--GS'):
            GS = val
        






