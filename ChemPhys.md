---
name: ChemPhys
topic: Chemometrics and Computational Physics
maintainer: Katharine Mullen
email: katharine.mullen@stat.ucla.edu
version: 2021-12-27
source: https://github.com/cran-task-views/ChemPhys/
---

Chemometrics and computational physics are concerned with the analysis
of data arising in chemistry and physics experiments, as well as the
simulation of physico-chemico systems. Many of the functions in base R
are useful for these ends.

The second edition of *Chemometrics with R: Multivariate Data Analysis
in the Natural and Life Sciences* by Ron Wehrens, ISBN
978-3-662-62027-4, Springer, 2020, provides an introduction to
multivariate statistics in the life sciences, as well as coverage of
several specific topics from the area of chemometrics. The associated
package `r github("rwehrens/ChemometricsWithR")` facilitates
reproduction of the examples in the book.

The book *Modern Statistical Methods for Astronomy With R Applications*
by Eric D. Feigelson and G. Jogesh Babu, ISBN-13: 9780521767279,
Cambridge, 2012, provides an introduction to statistics for astronomers
and an overview of the foremost methods being used in astrostatistical
analysis, illustrated by examples in R.

The book by Kurt Varmuza and Peter Filzmoser, *Introduction to
Multivariate Statistical Analysis in Chemometrics,* ISBN
978-1-420-05947-2, CRC Press, 2009, is associated with the package
`r pkg("chemometrics")`.

A special issue of R News with a focus on [R in
Chemistry](http://CRAN.R-project.org/doc/Rnews/Rnews_2006-3.pdf) was
published in August 2006. A special volume of Journal of Statistical
Software (JSS) dedicated to [oscopy and Chemometrics in
R](http://www.jstatsoft.org/v18/) was published in January 2007.

Please e-mail the maintainer, submit an issue or pull request in the GitHub
repository linked above, if we have omitted something of importance, or
if a new package or function should be mentioned here.

### Linear Regression Models

-   Linear models can be fitted (via OLS) with `lm()` (from stats). A
    least squares solution for `x` in `Ax = b` can also be computed as
    `qr.coef(qr(A), b)`.
-   The package `r pkg("nnls", priority = "core")` provides a
    means of constraining `x` to non-negative or non-positive values;
    the package `r pkg("bvls")` allows other bounds on `x`
    to be applied.
-   Functions for isotonic regression are available in the package
    `r pkg("Iso", priority = "core")`, and are useful to
    determine the unimodal vector that is closest to a given vector `x`
    under least squares criteria.
-   Heteroskedastic linear models can be fit using the `gls()` function
    of the `r pkg("nlme")` package.

### Nonlinear Regression Models

-   The `nls()` function (from stats) as well as the package
    `r pkg("minpack.lm")` allow the solution of nonlinear
    least squares problems.
-   Correlated and/or unequal variances can be modeled using the
    `gnls()` function of the `r pkg("nlme")` package and by
    `r pkg("nlreg")`.

### Curve Resolution

-   The `r pkg("PTAk", priority = "core")` package provides
    functions for Principal Tensor Analysis on k modes. The package
    includes also some other multiway methods: PCAn (Tucker-n) and
    PARAFAC/CANDECOMP.
-   Multivariate curve resolution alternating least squares (MCR-ALS) is
    implemented in the package
    `r pkg("ALS", priority = "core")`.
-   The `r bioc("alsace")` package provides MCR-ALS support
    for Liquid chromatography with PhotoDiode Array Detection (LC-DAD)
    data with many injections, with features for peak alignment and
    identification.
-   The package `r pkg("drc")` provides functions for the
    analysis of one or multiple non-linear curves with focus on models
    for concentration-response, dose-response and time-response data.

### Partial Least Squares

-   The package `r pkg("pls", priority = "core")` implements
    Partial Least Squares Regression (PLSR) and Principal Component
    Regression (PCR).
-   The package `r pkg("lspls")` implements the least
    squares-partial least squares (LS-PLS) method.
-   Sparse PLS is implemented in the package `r pkg("spls")`
    package.
-   The `r bioc("gpls")` package implements generalized
    partial least squares, based on the Iteratively ReWeighted Least
    Squares (IRWLS) method of Brian Marx.
-   The package `r pkg("enpls")` implements ensemble partial
    least squares, a framework for measuring feature importance, outlier
    detection, and ensemble modeling based on (sparse) partial least
    squares regressions.

### Principal Component Analysis

-   Principal component analysis (PCA) is in the package stats as
    functions `princomp()`. Some graphical PCA representations can be
    found in the `r pkg("psy")` package.
-   The `r pkg("homals")` package provides nonlinear PCA
    and, by defining sets, nonlinear canonical correlation analysis
    (models of the Gifi-family).
-   A desired number of robust principal components can be computed with
    the `r pkg("pcaPP")` package. The package
    `r pkg("elasticnet")` is applicable to sparse PCA. The
    package `r pkg("fpca")` can be applied to restricted MLE
    for functional PCA.
-   The `r pkg("subselect")` provides a collection of
    functions which assess the quality of variable subsets as surrogates
    for a full data set.
-   See the `r view("Multivariate")` task view for further
    packages dealing with PCA and other projection methods.

### Factor Analysis

-   Factor analysis (FA) is in the package stats as functions
    `factanal()`; see `r view("Psychometrics")` task view
    for details on extensions.

### Compositional Data Analysis

-   The package `r pkg("compositions")` provides functions
    for the consistent analysis of compositional data (e.g. portions of
    substances) and positive numbers (e.g. concentrations). See also the
    book, *Analyzing Compositional Data with R* by K. Gerald von den
    Boogaart und Raimon Tolosana-Delgado, ISBN: 978-3-642-36808-0,
    Springer, 2013.

### Independent Component Analysis

-   Independent component analysis (ICA) can be computed using
    `r pkg("fastICA")`.

### Clustering

-   The `r view("Cluster")` task view provides a list of
    packages that can be used for clustering problems.

### Variable Selection

-   Stepwise variable selection for linear models, using AIC, is
    available in function `step()`; package `r pkg("leaps")`
    implements leaps-and-bounds variable selection, by default using
    Mallow's Cp. `r pkg("stepPlr")` provides stepwise
    variable selection for penalized logistic regression.
-   Package `r pkg("varSelRF")` provides variable selection
    methods for random forests. Cross-validation-based variable
    selection using Wilcoxon rank sum tests is available in package
    `r pkg("WilcoxCV")`, focused on binary classification in
    microarrays. Package `r pkg("clustvarsel")` implements
    variable selection for model-based clustering.
-   The `r pkg("BioMark")` package implements two
    meta-methods for variable selection: stability selection (applying a
    primary selection method like a t-test, VIP value or PLSDA
    regression coefficient) to different subsets of the data, and higher
    criticism, which provides a data-driven choice of significance
    cutoffs in statistical testing.

### Self-Organizing Maps

-   The `r pkg("kohonen", priority = "core")` package
    implements self-organizing maps as well as some extensions for
    supervised pattern recognition and data fusion. The
    `r pkg("som")` package provides functions for
    self-organizing maps.

### Differential Equations

-   See the `r view("DifferentialEquations")` task view
    packages dealing with differential equations.

### Metrology

-   The `r pkg("units")` package attaches unit metadata to
    vectors, matrices and arrays, providing automatic propagation,
    conversion, derivation and simplification of units.
-   The `r pkg("errors")` attaches uncertainty metadata to
    vectors, matrices and arrays, providing automatic propagation and
    reporting.
-   The `r pkg("constants")` package provides values of the
    fundamental physical constants based on values reported by the
    Committee on Data for Science and Technology (CODATA), an
    interdisciplinary committee of the International Council for
    Science.
-   `r pkg("NISTunits")` also provides values of the
    fundamental physical constants. The values it contains are based on
    the values reported by the National Institute of Standards and
    Technology, (NIST).
-   The `r pkg("measurements")` contains tools to make
    working with physical measurements easier, such as functions to
    convert between metric and imperial units, or to calculate a
    dimension's unknown value from other dimensions' measurements.
-   The `r pkg("metRology")` package provides support for
    metrology applications, including measurement uncertainty estimation
    and inter-laboratory metrology comparison studies.
-   The `r pkg("ATmet")` package provides functions for
    smart sampling and sensitivity analysis for metrology applications,
    including computationally expensive problems.

### Calibration

-   The `r pkg("investr")` package facilitates
    calibration/inverse estimation with linear and nonlinear regression
    models.
-   The `r pkg("chemCal", priority = "core")` package
    provides functions for plotting linear calibration functions and
    estimating standard errors for measurements.
-   The `r pkg("nlreg")` package is useful for nonlinear
    calibration models.
-   The package `r pkg("represent")` calculates the
    'representativity' of two multidimensional datasets, which
    involves comparison of the similarity of principal component
    analysis loading patterns, variance-covariance matrix structures,
    and data set centroid locations.

### Cellular Automata

-   The `r pkg("simecol")` package includes functions for
    cellular automata modeling.

### Thermodynamics

-   The `r pkg("CHNOSZ")` package provides functions for
    calculating the standard Gibbs energies and other thermodynamic
    properties, and chemical affinities of reactions between species
    contained in a thermodynamic database.

### Interfaces to External Libraries

-   The package `r pkg("rcdk")` allows the user to access
    functionality in the [Chemistry Development Kit
    (CDK),](http://sourceforge.net/projects/cdk/) a Java framework for
    cheminformatics. This allows the user to load molecules, evaluate
    fingerprints (via the package `r pkg("fingerprint")`),
    calculate molecular descriptors and so on. In addition, the CDK API
    allows the user to view structures in 2D. The
    `r pkg("rcdklibs")` package provides the CDK libraries
    for use in R.
-   `r bioc("ChemmineR")` is a cheminformatics toolkit for
    analyzing small molecules in R. Its add-on packages include
    `r bioc("fmcsR")` for mismatch tolerant maximum common
    substructure matching, `r bioc("eiR")` for accelerated
    structure similarity searching; `r bioc("bioassayR")`
    for analyzing bioactivity data, and
    `r bioc("ChemmineOB")` for accessing
    [OpenBabel](http://openbabel.org/wiki/Main_Page) functionalities
    from R.
-   The `r pkg("webchem")` package allows users to retrieve
    chemical information from various sources on the web and to interact
    with various APIs. Sources include: [Chemical Identifier
    Resolver](http://cactus.nci.nih.gov/chemical/structure) ,
    [ChemSpider](http://www.chemspider.com/) ,
    [PubChem](https://pubchem.ncbi.nlm.nih.gov/) , [Chemical Translation
    Service](http://cts.fiehnlab.ucdavis.edu/) , [PAN Pesticide
    Database](http://www.pesticideinfo.org/) , [Alan Wood's Compendium
    of Pesticide Common Names](http://www.alanwood.net/pesticides/) ,
    [PHYSPROP
    Database](http://www.srcinc.com/what-we-do/environmental/scientific-databases.html)
    , [ETOX](http://webetox.uba.de/webETOX/index.do) ,
    [PPDB](http://sitem.herts.ac.uk/aeru/iupac/search.htm) , and
    [ChemIDplus](http://chem.sis.nlm.nih.gov/chemidplus/) .

### Spectroscopy

-   Bryan Hanson has compiled a broad range of [Free and Open Source
    Software (FOSS) for
    Spectroscopy](https://bryanhanson.github.io/FOSS4Spectroscopy/) ,
    much of which is in the form of R packages.
-   The `r pkg("spectralAnalysis")` package allows users to
    pre-process, visualize and analyze spectroscopy data. Non-negative
    matrix factorization analysis is included.
-   The `r pkg("ChemoSpec")` package collects user-friendly
    functions for plotting spectra (NMR, IR, etc) and carrying top-down
    exploratory data analysis, such as HCA, PCA and model-based
    clustering.
-   The `r github("Chathurga/HyperChemoBridge")`
    interconverts `r pkg("ChemoSpec")` (and hyperSpec)
    objects
-   The `r pkg("speaq")` package implements the hierarchical
    Cluster-based Peak Alignment (CluPA) and may be used for aligning
    NMR spectra.
-   The package `r pkg("TIMP")` provides a problem solving
    environment for fitting separable nonlinear models in physics and
    chemistry applications, and has been extensively applied to
    time-resolved spectroscopy data.
-   The package `r pkg("ChemoSpec2D")` allows exploratory
    chemometrics of 2D spectroscopic data sets such as COSY (correlated
    spectroscopy) and HSQC (heteronuclear single quantum coherence) 2D
    NMR (nuclear magnetic resonance) spectra.
-   The `r pkg("spectrino")` package provides tools for
    spectra viewing and organization.

### Mass Spectrometry

-   The `r bioc("MSnbase")` defines infrastructure for mass
    spectrometry-based proteomics data handling, plotting, processing
    and quantification.
-   The `r pkg("MALDIquant")` provides tools for
    quantitative analysis of MALDI-TOF mass spectrometry data, with
    support for baseline correction, peak detection and plotting of mass
    spectra.
-   The `r pkg("OrgMassSpecR")` package is for
    organic/biological mass spectrometry, with a focus on graphical
    display, quantification using stable isotope dilution, and protein
    hydrogen/deuterium exchange experiments.
-   The Bioconductor packages `r bioc("MassSpecWavelet")`,
    `r bioc("PROcess")`, and `r bioc("xcms")`
    are designed for the analysis of mass spectrometry data.
-   The [apLCMS](http://web1.sph.emory.edu/apLCMS/) package is designed
    for the processing of LC/MS based metabolomics data.
-   The [xMSanalyzer](http://sourceforge.net/projects/xmsanalyzer/)
    package allows merging [apLCMS](http://web1.sph.emory.edu/apLCMS/)
    sample processing results from multiple sets of parameter settings,
    among other features.
-   The [MSPrep](http://sourceforge.net/projects/msprep/) package is for
    post-processing of metabolomic data, including summarization of
    replicates, filtering, imputation, and normalization.
-   The `r bioc("metaMS")` package is an MS-based
    metabolomics data processing and compound annotation pipeline.

### Functional Magnetic Resonance Imaging

-   Functions for I/O, visualization and analysis of functional Magnetic
    Resonance Imaging (fMRI) datasets stored in the ANALYZE or NIFTI
    format are available in the package
    `r pkg("AnalyzeFMRI")`. The package
    `r pkg("fmri")` contains functions to analyze fMRI data
    using adaptive smoothing procedures.

### Fluorescence Lifetime Imaging Microscopy

-   Functions for visualization and analysis of Fluorescence Lifetime
    Imaging Microscopy (FLIM) datasets are available in the package
    `r pkg("TIMP")`.

### Fluorescence Excitation-Emission Matrix (EEM)

-   The `r pkg("EEM")` reads raw EEM data and prepares it
    for further analysis.

### Carbon Dating

-   The package `r pkg("Bchron")` creates chronologies based
    on radiocarbon and non-radiocarbon dated depths.

### Astronomy and astrophysics

-   The `r pkg("astrodatR")` package collects 19 datasets
    from contemporary astronomy research, many of which are described in
    the aforementioned textbook 'Modern Statistical Methods for
    Astronomy with R Applications'.
-   The `r pkg("astrolibR")` package presents an R interface
    to low-level utilities and codes from the [Interactive Data Language
    (IDL) Astronomy Users Library](http://idlastro.gsfc.nasa.gov) .
-   The `r pkg("CRAC")` collects R functions for
    cosmological research, with its main functions being similar to the
    python library, cosmolopy.
-   The `r pkg("RobPer")` package calculates periodograms
    based on (robustly) fitting periodic functions to light curves.
-   The package `r pkg("snapshot")` contains functions for
    reading and writing N-body snapshots from the GADGET code for
    cosmological N-body/SPH simulations.
-   The package `r pkg("UPMASK")` performs unsupervised
    photometric membership assignment in stellar clusters using, e.g.,
    photometry and spatial positions.
-   The `r pkg("solaR")` package provides functions to
    determine the movement of the sun from the earth and to determine
    incident solar radiation.
-   The `r pkg("FITSio")` package provides utilities to read
    and write files in the FITS (Flexible Image Transport System)
    format, a standard format in astronomy.
-   The `r pkg("stellaR")` package manages and displays
    stellar tracks and isochrones from the Pisa low-mass database.
-   The `r pkg("astroFns")` provides miscellaneous astronomy
    functions, utilities, and data.
-   The `r pkg("cosmoFns")` contains standard expressions
    for distances, times, luminosities, and other quantities useful in
    observational cosmology, including molecular line observations.
-   The `r pkg("celestial")` package includes a number of
    common astronomy conversion routines, particularly the HMS and
    degrees schemes.
-   The `r pkg("SCEPtER")` package is used to estimate
    stellar mass and radius given observational data of effective
    temperature, \[Fe/H\], and astroseismic parameters.
-   The `r pkg("lira")` package performs Bayesian linear
    regression and forecasting in Astronomy, accounting for all kinds of
    errors and correlations in the data.
-   The `r pkg("SPADAR")` package provides functions to
    create all-sky grid plots of widely used astronomical coordinate
    systems (equatorial, ecliptic, galactic) and scatter plots of data
    on any of these systems, including on-the-fly system conversion.
-   The `r pkg("SCEPtERbinary")` allows for estimating the
    stellar age for double-lined detached binary systems, adopted from
    the effective temperature, the metallicity \[Fe/H\], the mass, and
    the radius of the two stars.
-   The [Astrostatistics and Astroinformatics
    Portal](https://asaip.psu.edu) is an R-centric collection of
    information regarding statistical analysis in astronomy.
-   Hans Werner Borchers has a page on [Astronomy modules and links for
    R, Python, and
    Julia](https://github.com/hwborchers/zaRastro/blob/master/README.md)
    .

### Optics and Scattering Approximations

-   The `r pkg("planar")` package provides code to simulate
    reflection and transmission at a multilayer planar interface.
-   The `r pkg("dielectric")` package defines some physical
    constants and dielectric functions commonly used in optics and
    plasmonics.

### Energy Modeling

-   The `r pkg("solaR")` package provides functions to
    simulate and model systems involved in the capture and use of solar
    energy, including photovoltaics.

### Water and Soil Chemistry

-   The `r pkg("AquaEnv")` package is a toolbox for aquatic
    chemical modelling focused on (ocean) acidification and CO2
    air-water exchange.
-   See the `r view("Environmetrics")` task view for further
    related packages related to water and soil chemistry.

### Titration Curves

-   The `r pkg("titrationCurves")` package provides
    functions to plot acid/base, complexation, redox, and precipitation
    titration curves.

### Electrochemistry

-   The `r pkg("eChem")` package provides functions to
    simulate voltammetry, chronoamperometry and chronocoulometry
    experiments, which may be useful in courses in analytical chemistry.

### Health Physics

-   The package `r pkg("radsafer")` provides functions for
    radiation safety; the package `r pkg("RadData")`
    provides nuclear decay data for dosimetric calculations from the
    International Commission on Radiological Protection.



### Links
-   [R News: R in Chemistry](http://CRAN.R-project.org/doc/Rnews/Rnews_2006-3.pdf)
-   [Journal of Statistical Software: Spectroscopy and Chemometrics in R](http://www.jstatsoft.org/v18/)
-   [apLCMS](http://web1.sph.emory.edu/apLCMS/)
-   [Interactive Data Language (IDL) Astronomy Users Library](http://idlastro.gsfc.nasa.gov)
-   [Astrostatistics and Astroinformatics Portal Software Forum](https://asaip.psu.edu/forums/software-forum)
-   [Chemistry Development Kit (CDK)](http://sourceforge.net/projects/cdk/)
-   [MSPrep](http://sourceforge.net/projects/msprep/)
-   [OpenBabel](http://openbabel.org/wiki/Main_Page)
-   [PubChem](http://pubchem.ncbi.nlm.nih.gov/)
-   [xMSanalyzer](http://sourceforge.net/projects/xmsanalyzer/)
-   [Chemical Identifier Resolver](http://cactus.nci.nih.gov/chemical/structure)
-   [ChemSpider](http://www.chemspider.com/)
-   [Free and Open Source Software (FOSS) for Spectroscopy](https://bryanhanson.github.io/FOSS4Spectroscopy/)
-   [Chemical Translation Service](http://cts.fiehnlab.ucdavis.edu/)
-   [PAN Pesticide Database](http://www.pesticideinfo.org/)
-   [Alan Wood's Compendium of Pesticide Common Names](http://www.alanwood.net/pesticides/)
-   [PHYSPROP Database](http://www.srcinc.com/what-we-do/environmental/scientific-databases.html)
-   [ETOX](http://webetox.uba.de/webETOX/index.do)
-   [PPDB](http://sitem.herts.ac.uk/aeru/iupac/search.htm)
-   [ChemIDplus](http://chem.sis.nlm.nih.gov/chemidplus/)
-   [Astronomy modules and links for R, Python, and Julia](https://github.com/hwborchers/zaRastro/blob/master/README.md)
