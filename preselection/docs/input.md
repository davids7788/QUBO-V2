# Configuration file preselection
The preselection script takes YAML configuration file. 
It consists of the sections [doublet](#doublet), [triplet](#triplet), [binning](#binning), [qubo](#qubo) and [scaling](#scaling).

# doublet
Criterion for a preselection on doublet level. Creating a Doublet object, if a set of two hits fulfill the criterion.
* `dx/x0:` (x2 - x1) / x0 with x1 and x2 coordinates of hits and x0 as linearly interpolated hit on reference layer (float)
* `eps:`. maximum deviation from the d/dx criteria (float)

Criteria for a preselection on triplet level. Creating a Doublet object, if a set of two doublets fulfill both criteria.
# triplet
* `angle diff x:` maximum allowed angle in x (float)
* `angle diff y:` maximum allowed angle in y (float)
  
# binning
To reduce computational costs, separating hits into bins. 
* `num bins x:` number of bins in x-direction (int)

# qubo
QUBO parameters. Separating parameters for quadratic terms to match and conflict.
* `b_ij conflict:` linear term  (conflict) of triplets not creating a track, (float) or name of function from self written class in src/functions
* `b_ij match:` linear term (interaction) of triplets creating a track candidate, (float) or name of function from self written class in src/functions
* `a_i:` quadratic term (quality) of triplet term of 

# scaling:
Rescaling QUBO parameters
* `z_scores:` True if using z-scores on quadratic term (quality), else False
* `quality:` Interval [min (float), max (float)] if rescaling, else False
* `interaction:` Interval [min (float), max (float)] if rescaling, else False