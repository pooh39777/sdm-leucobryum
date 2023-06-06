# README
## Code version

### SDM_Leucobryum_aduncum_20210822.R
- R version

### SDM_Leucobryum_aduncum_20220506.Rmd
- Move to R markdown

### SDM_Leucobryum_aduncum_20220621.Rmd
- Future prediction test (wrong method)

### SDM_Leucobryum_aduncum_20220808.Rmd
- Clean code
- Remove wrong future

### SDM_Leucobryum_aduncum_20220815.Rmd
- Using SDMtune package

### SDM_Leucobryum_aduncum_20221002.Rmd
- Spatial cross validation with blockCV (Valavi et al., 2019)

### SDM_Leucobryum_aduncum_20221003.Rmd
- Spatial cross validation with ENMeval (Muscarella et al., 2014) by Block method partition
- Change in SEA domain follow to Ge et al. (10° S–23° N, 95° E–140° E)
- Area Calculation from https://caucasus-spiders.info/r-spatial/raster-basics-3/
- Selected variable model (about 28 min.)
- Tune hyperparameters with optimizeModel (about 7 min.)
- Plot species occurrence with mapview()
- Shortening bioclim data path with "./"
- Fix bug (upgrade moving from raster to terra)
- Change map color
- Change ggsave map size
- Started using GitHub with this version on 2023-06-06
