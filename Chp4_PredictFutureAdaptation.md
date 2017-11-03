# Predict genetic offset with future climates

Since I used some pond-specific variables, it's going to be impossible to get widespread current and future data for these variables. 

Instead I will download all the bioclim data and find which variables are most correlated with Temp and Season within each gradient. 
I'll use these to predict the genotypes across the region under current climate, and to predict genotypes needed under future climates. 

i.e. in the env file, I'll change the BioClim variables to Temp and Season for the function to work. 

BioClim data resolution: 30second resolution 0.86km at the equator. 
