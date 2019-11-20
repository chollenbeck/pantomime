library(devtools)
library(usethis)

#create("pantomime")

use_data_raw()

##### DESCRIPTION Setup #####

# Remove default DESC
unlink("DESCRIPTION")

# Create a new description object
my_desc <- desc::description$new("!new")

# Set your package name
my_desc$set("Package", "pantomime")
#Set your name

my_desc$set("Authors@R", "person('Chris', 'Hollenbeck', email =
'chollenbeck07@gmail.com', role = c('cre', 'aut'))")

# Remove some author fields
my_desc$del("Maintainer")

# Set the version
my_desc$set_version("0.0.0.9000")

# The title of your package
my_desc$set(Title = "Tools for tidy population genetics")

# The description of your package
my_desc$set(Description = "An R package containing tools for population genetics analysis in R.")

# Save everyting
my_desc$write(file = "DESCRIPTION")

##### Misc Setup #####

use_roxygen_md()
use_pipe()
use_testthat()
use_git()
