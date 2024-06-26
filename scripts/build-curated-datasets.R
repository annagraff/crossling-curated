# invoke the build scripts in the respective folders
if(!file.exists("scripts/build-curated-datasets.R")) rlang::abort(
  "this script should be invoked from the project folder",
)


cat("Building curated GBI datasets\n")
source("scripts/GBI/build-datasets.R")
source("scripts/GBI/test-dependencies.R")
source("scripts/GBI/densify-datasets.R")

cat("Building curated TLI datasets\n")
source("scripts/TLI/compile-external-input.R")
source("scripts/TLI/build-datasets.R")
source("scripts/TLI/test-dependencies.R")
source("scripts/TLI/densify-datasets.R")
