from cobra.io import load_model

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion, flux_variability_analysis)

model = load_model("iML1515")
# Default objective is biomass growth
# BIOMASS_Ec_iML1515_WT_75p37M <- Reaction name

print(f"Default solution value: {model.optimize().objective_value}")

# Test1 : Environment without o2
with model:
    medium = model.medium
    medium["EX_o2_e"] = 0.0
    model.medium = medium
    print(f"No o2 in the medium solution: {model.optimize().objective_value}")


# Test2 : CYTK2 knockout
with model:
    model.reactions