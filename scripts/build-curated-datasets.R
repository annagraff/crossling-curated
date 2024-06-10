# invoke the build scripts in the respective folders
if(!file.exists("scripts/build-curated-datasets.R")) rlang::abort(
  "this script should be invoked from the project folder",
)


cli::cli_h1("Building GBInd")
source("scripts/GBInd/build-datasets.R")