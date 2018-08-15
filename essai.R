# script using function , that is the "main" 

# note : matrix : genetData coordinate + where there are 

# donnees environnementale : environmentalData

# BUT : comparer les deux méthodes ( matrice absorbante et chaine de markov)

# matrice absorbante

# créer le modèle démographique : fonction ReactNorm mais que signifie p et choisir "share

source("initialisation.R")
source("Laurianne.R") 

ReactNorm(environmentalData,835,shapes="enveloppe")


# matrice de migration 

migrationMatrix(environmentalData,prior$dispersion$distribution,prior$dispersion$p)
