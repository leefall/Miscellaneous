# Install from CRAN
install.packages("ukbtools")

#Install latest development version
#devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE)

library(ukbtools)

#my_ukb_data<-ukb_df("",path="/mnt/towel/UKBiobank")


# To load the example data
path_to_example_data <- system.file("extdata", package = "ukbtools")
path_to_example_data
df <- ukb_df("ukbxxxx", path = path_to_example_data)

# To create a field code to name key
df_field <- ukb_df_field("ukbxxxx", path = path_to_example_data)



ukb_context(my_ukb_data, nonmiss.var = "my_variable_of_interest")



