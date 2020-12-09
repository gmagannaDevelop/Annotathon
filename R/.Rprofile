# Renv :
source("renv/activate.R")

# Chacking our BioConductor installation is valid:
print("Checking that we have a valid BioConductor instance")
BiocManager::valid()


# Parse config
config <- RcppTOML::parseToml("config.toml")

# If willing to load  packages from the config.toml file
# uncomment the following : 
# lapply(   , require, character.only = TRUE)

# Enable reticulate
#Sys.setenv(RETICULATE_PYTHON=config$paths$python)

# Enhance reproducibility
set.seed(config$random$seed)
#py_set_seed(config$random$seed)
#data.lock <- parseToml(config$data$lockfile)

