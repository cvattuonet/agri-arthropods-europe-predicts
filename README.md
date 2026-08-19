# Arthropod Biodiversity and Agricultural Land Use in Europe

# Abstract

Agriculture is often described as one of the main drivers of arthropod decline, 
but the contribution of the individual practices behind it has not been disentangled 
at large scales. Analyses at that scale represent agriculture through broad 
land-use classes and generally reduce management intensity to a single proxy. 
This is partly because the global spatial data describing what is actually applied
in the field are coarse or absent, and partly because the biodiversity records 
available are compiled from many heterogeneous sources. In Europe, some products 
describing agricultural inputs and the structure of the agricultural landscape are 
available at a finer resolution than their global equivalents. This study asked whether
a global biodiversity database, combined with independent large-scale agronomic datasets,
can be used to separate the effects of individual agricultural practices and of the surrounding 
landscape composition on arthropods in Europe, and whether those effects differ between habitat types.

Records of Insecta and Arachnida for the European Union and the United Kingdom were 
extracted from the 2023 release of PREDICTS and curated reference by reference against the 
primary publications. The dataset comprised 235,661 records from 60 studies 
covering 2,456 sites sampled between 2000 and 2011. Each site was assigned a local 
habitat category crossing land-use classes with use intensity, and was characterized
within a radius of 1 kilometer by land cover, pesticide and nitrogen application, crop diversity, 
field size and climate, drawn from public spatial products. Total abundance and taxa 
richness responses were then modelled with mixed-effects models.

Habitats under little or no management held higher taxa richness than intensively managed 
ones, with abundance following the same direction but with weaker support. The least disturbed
habitat was not the richest: minimally used pasture held 96% more individuals and 21% more taxa 
than natural vegetation, and lightly used semi-natural vegetation 72% and 36% more. A clear 
gradient across use intensity appeared only in pastures. Pesticide application in the 
surrounding landscape reduced abundance by 12.9% and richness by 7.5%, with no evidence that 
the effect differed between habitats. Effect of crop diversity and the prevalence of 
small fields depends on local habitat, with greater crop diversity associated 
with higher abundance and richness in intensively used cropland, but lower in natural vegetation.
Contrary to what was expected, natural vegetation cover in the landscape did not buffer 
these pressures. Marginal R² was low throughout (0.034 to 0.051), a pattern common to studies 
using this database that nonetheless calls for testing alternative modelling approaches.

Nitrogen application could not be separated from pesticides, and no spatial products were 
available for practices such as tillage, rotation or organic management. Separating individual 
practices at this scale therefore remains limited by the agronomic data, and by the biases of 
the biodiversity database. The analysis nonetheless provides, to our knowledge, the first 
estimate of the effect of pesticide application on the abundance and richness of insects 
and arachnids across habitat types at European scale, together with a curated dataset that 
can be used in further work.


# Overview
This repository contains the full R pipeline used to build a site-level
dataset of arthropod biodiversity across European PREDICTS sites, extract
agronomic and landscape predictors around each site, and fit mixed-effects
models (Gaussian, Poisson, negative binomial, Tweedie, zero-inflated NB, and
LASSO-penalized) relating biodiversity outcomes to those predictors.


# Pipeline 
Run the scripts in order; each one reads the outputs of the previous steps.

1. `1_Biodiversity_Dataset_Creation.R` - build the initial site-level
   biodiversity dataset from the PREDICTS database.
2. `2_Biodiversity_Dataset_Curation.R` - clean and curate site coordinates and
   metadata.
3. `3_Agronomic_Databases_Preparation.R` - prepare agronomic raster layers
   (pesticides, fertilizers, crop diversity) from raw source databases.
4. `4_Fixed_effects_calculation.R` - extract landscape/agronomic predictors
   around each site (land cover, climate, field size, pesticides,
   fertilizers, crop diversity, water, etc.).
5. `5_Dataset_final_preparation.R` - merge biodiversity and predictor data
   into the final abundance/richness model-ready tables.
6. `6_Null_models.R` - fit null (intercept-only) models for comparison.
7. `7_Models.R` - fit the main mixed-effects models (Gaussian, Poisson,
   negative binomial, Tweedie, zero-inflated NB) and run model selection
   (`MuMIn::dredge`) for abundance and richness.
8. `8_Lasso_Models.R` - fit LASSO-penalized mixed models (`glmmLasso`) as a
   complementary variable-selection approach.
9. `9_Results_Extraction.R` - extract model diagnostics, averaged effects,
   and LASSO results into a formatted Excel workbook.
10. `10_Figures.R` - reproduces the manuscript/thesis figures.


# Data availability
The raw source databases used by this pipeline are too large to host on
GitHub (several hundred GB in total). Each subfolder under `Databases/`
contains a `README.md` stub with the dataset name and source/download link.
Download each dataset and place its files in the corresponding
folder before running the scripts.

One exception is `Databases/Smith2006_coordinates/`, which is committed
directly rather than left as a download task. These are site coordinates
shared privately by one of the primary authors of the underlying study (see
that folder's `README.md`).

`Intermediate_dataset/` is regenerated automatically by running scripts 1-5.
Three key site-level summary tables are committed directly so results can be explored without
re-running the full pipeline.

`Models_results/` (fitted models and `dredge()` model-selection objects) is
small enough to be included directly in the repository although it can be obtained 
by running scripts 7-8. 
