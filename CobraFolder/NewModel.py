from cobra import Model, Reaction, Metabolite

model = Model("NewModel")

# --- Helper Functions ---
def create_reaction(id_num, name, genes, metabolites_dict, low_bound=0, up_bound=1000, subsystem=""):
    reaction = Reaction(
        id=id_num,
        name=name,
        subsystem=subsystem, # Added ="" so that it doesnt create an error if i forget to add it
        lower_bound=low_bound,
        upper_bound=up_bound
    )
    reaction.gene_reaction_rule = genes

    reaction.add_metabolites(metabolites_dict)

    return reaction

def create_metabolites(id_name, formula, name, compartment=""):
    metabolite = Metabolite(
        id=id_name,
        formula=formula,
        name=name,
        compartment=compartment
    )

    return metabolite

# --- Metabolite Definitions ---
glc_EX = create_metabolites("glc_EX", "C6H12O6", "Glucose", "e")
glc_p = create_metabolites("glc_p", "C6H13O9P", "Alpha-D-Glucose-6-Phosphate", "c")
frc_p = create_metabolites("frc_p", "C6H13O9P", "Fructose-6-Phosphate", "c")
atp = create_metabolites("atp", "C10H16N5O13P3", "Adenosine-5-triphosphate", "c")
adp = create_metabolites("adp", "C10H15N5O10P2", "Adenosine-5-diphosphate", "c")
frc_dp = create_metabolites("frc_dp", "C6H14O12P2", "D-Fructose-1,6-bisphosphate", "c")
h2o = create_metabolites("h2o", "H2O", "Water", "c")
pi = create_metabolites("pi", "H3PO4", "Phosphate", "c")
glyc_p = create_metabolites("glyc_p", "C3H7O6P", "Glycerone phosphate", "c")
glycal_p = create_metabolites("glycal_p", "C3H7O6P", "D-Glyceraldehyde 3-phosphate", "c")
nad = create_metabolites("nad", "C21H28N7O14P2", "NAD+", "c")
nadh = create_metabolites("nadh", "C21H29N7O14P2", "NADH", "c")
h = create_metabolites("h", "H", "Proton", "c")
glyclr_p = create_metabolites("glyclr", "C3H8O10P2", "3-Phospho-D-glyceroyl phosphate", "c")
glyclr3_p = create_metabolites("glyclr3_p", "C3H7O7P", "3-Phospho-D-glycerate", "c")
glyclr2_p = create_metabolites("glyclr2_p", "C3H7O7P", "2-Phospho-D-glycerate", "c")
pep = create_metabolites("pep", "C3H5O6P", "Phosphoenolpyruvate", "c")
pyr = create_metabolites("pyr", "C3H4O3", "Pyruvate", "c")
amp = create_metabolites("amp", "C10H14N5O7P", "Adenosine-5-monophosphate", "c")

# --- Reaction Definitions ---
R_27199 = create_reaction("R02738", "protein-N(pi)-phosphohistidine:D-glucose 6-phosphotransferase", "(b1101 or b1621)", {glc_EX: -1.0, glc_p: 1.0})
R_5319 = create_reaction("R13199", "alpha-D-glucose-6-phosphate aldose-ketose-isomerase", "(b4025)", {glc_p: -1.0, frc_p: 1.0}, low_bound=-1000)
R_27111 = create_reaction("R00756", "ATP:D-fructose-6-phosphate 1-phosphotransferase", "(b1723 or b3916)", {atp: -1.0, frc_p: -1.0, adp: 1.0, frc_dp: 1.0})
R_31311 = create_reaction("R00762", "D-fructose-1,6-bisphosphate 1-phosphohydrolase", "(b2930 or b3925 or b4323)", {frc_dp: -1.0, h2o: -1.0, frc_p: 1.0, pi: 1.0})
R_41213 = create_reaction("R01068", "D-fructose-1,6-bisphosphate D-glyceraldehyde-3-phosphate-lyase", "(b2097 or b2925)", {frc_dp: -1.0, glyc_p: 1.0, glycal_p: 1.0}, low_bound=-1000)
R_5311 = create_reaction("R01015", "D-glyceraldehyde-3-phosphate aldose-ketose-isomerase", "(b3919)", {glycal_p: -1.0, glyc_p: 1.0}, low_bound=-1000)
R_12112 = create_reaction("R01061", "D-glyceraldehyde-3-phosphate:NAD+ oxidoreductase (phosphorylating)", "(b1779)", {glycal_p: -1.0, pi: -1.0, nad: -1.0, glyclr_p: 1.0, nadh: 1.0, h: 1.0}, low_bound=-1000)
R_2723 = create_reaction("R01512", "ATP:3-phospho-D-glycerate 1-phosphotransferase", "(b2926)", {glyclr_p: -1.0, adp: -1.0, glyclr3_p: 1.0, atp: 1.0}, low_bound=-1000)
R_54211 = create_reaction("R01518", "D-phosphoglycerate 2,3-phosphomutase", "(b3612 or b0755 or b4395)", {glyclr2_p: -1.0, glyclr3_p: 1.0}, low_bound=-1000)
R_42111 = create_reaction("R00658", "2-phospho-D-glycerate hydro-lyase (phosphoenolpyruvate-forming)", "(b2779)", {glyclr2_p: -1.0, pep: 1.0, h2o: 1.0}, low_bound=-1000)
R_27140 = create_reaction("R00200", "ATP:pyruvate 2-O-phosphotransferase", "(b1676 or b1854)", {pep: -1.0, adp: -1.0, h: -1.0, atp: 1.0, pyr: 1.0})
R_2792 = create_reaction("R00199", "ATP:pyruvate,water phosphotransferase", "(b1702)", {h2o: -1.0, atp: -1.0, pyr: -1.0, amp: 1.0, pi: 1.0, pep: 1.0})

# Needed Reactions Because The Cell Doesn't Get Enough
atp_reaction = create_reaction("ATP_M", "ADP Generation for balance", "", {atp: -1.0, adp: 1.0, pi: 1.0})
# R5319 -> Uses 1 ATP, produces 1 ADP
# R2723 -> Uses 1 ADP, prodcues 1 ATP
# R27140 -> Uses 1 ADP, produces 1 ATP
# 1 ADP is needed for the whole cycle to continue

nadh_oxidation = create_reaction("NADH_OX", "NADH Oxidation", "", {nadh: -1.0, nad: 1.0})
# R12112 -> Uses 1 NAD, produces 1 NADH
# No NAD source to get NAD

model.add_reactions([
    R_27199, 
    R_5319, 
    R_27111, 
    R_31311, 
    R_41213, 
    R_5311, 
    R_12112,
    R_2723, 
    R_54211, 
    R_42111, 
    R_27140, 
    R_2792,
    atp_reaction, 
    nadh_oxidation
])

model.add_boundary(model.metabolites.get_by_id("glc_EX"), type="exchange", lb=-10.0)
model.add_boundary(model.metabolites.get_by_id("pyr"), type="demand", lb=0.0)
# There is too much proton
model.add_boundary(model.metabolites.get_by_id("pi"), type="demand", lb=0.0)

model.objective = "R00200"

