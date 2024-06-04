# invoke the build scripts in the respective folders
if(!file.exists("scripts/build-curated-datasets.R")) rlang::abort(
  "this script should be invoked from the project folder",
)


cat("Building curated GBInd datasets\n")
source("scripts/GBInd/create-full-datasets.R")
source("scripts/GBInd/densify-datasets.R")
source("scripts/GBInd/test-dependencies.R")


cat("Building curated TypLinkInd datasets\n")
source("scripts/TypLinkInd/compile-external-input.R")
source("scripts/TypLinkInd/create-full-datasets.R")
source("scripts/TypLinkInd/densify-datasets.R")
source("scripts/TypLinkInd/test-dependencies.R")
