# Configuration file for `make_triplets.py`
The `make_triplets.py` script takes a YAML configuration file as first input argument. 
It consists of the sections [doublet](#doublet), [triplet](#triplet), [binning](#binning), [qubo](#qubo) and [scaling](#scaling).

# doublet
Criterion for a preselection on doublet level. Creating a Doublet object, if a set of two hits fulfill the criterion.
* `dx/x0:` $\frac{(x_2 - x_1)}{x_0}$ with $x_1$ (float) and $x_2$ (float) as x-coordinates of hits ond 
  (consecutive layers) and $x_0$ as linearly interpolated hit on reference layer (float)
* `eps:` maximum deviation from the d/dx criteria (float)
* `dy/x0:` $\frac{(x_y - x_y)}{x_0}$ with $y_1$ (float) and $y_2$ (float) as y-coordinates of hits ond 
  (consecutive layers) and $x_0$ as linearly interpolated hit on reference layer (float)

Criteria for a preselection on triplet level. Creating a Doublet object, if a set of two doublets fulfill both criteria.
# triplet
* `max scattering:` maximum allowed value for scattering, calculated via $\sqrt{angle_{xz}^2 + angle_{yz}^2}$ (float)
  
# binning
To reduce computational costs, separating hits into bins. 
* `num bins x:` number of bins in x (int)
* `num bins y:` number of bins in y (int)

# qubo
QUBO parameters. Separating parameters for quadratic terms to match and conflict. The QuboCoefficients class administers the functions. 
To add a function, ensure that it is added to the function dictionary with a name as a key and the function as the value.
* `b_ij conflict:` quadratic tern  (conflict) of triplets not creating a track, (float) or name of function
* `b_ij match:` quadratic term (interaction) of triplets creating a track candidate, (float) or name of function
* `a_i:` linear term (quality) of triplet term, (float) or name of function

# scaling:
Rescaling QUBO parameters
* `z_scores:` True if using z-scores on linear term (quality), else False
* `quality:` interval [min (float), max (float)] if rescaling, else null
* `interaction:` interval [min (float), max (float)] if rescaling, else null