import pandas as pd

from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion,
    double_gene_deletion, double_reaction_deletion,
    flux_variability_analysis
)

pd.set_option("display.float_format", "{:.10f}".format)

model = Model("ArginineModel")

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

glt = create_metabolites("glt","C5H9NO4","L-Glutamate","c")
asp = create_metabolites("asp","C4H7NO4","L-Aspartate","c") # Comes from Alanine, glutamate, aspartate metabolism
oxa = create_metabolites("oxa","C4H4O5","Oxaloacetate","c")
oxa_glt = create_metabolites("oxa_glt","C5H6O5","2-Oxoglutarate","c") # Comes from TCA Cycle
pyr = create_metabolites("pyr", "C3H4O3", "Pyruvate", "c")
ala = create_metabolites("ala","C3H7NO2","L-Alanine","c")
ac_coa = create_metabolites("ac_coa","C23H38N7O17P3S","Acetyl-CoA","c")
ac_glt = create_metabolites("ac_glt","C7H11NO5","N-Acetyl-L-glutamate","c")
coa = create_metabolites("coa","C21H36N7O16P3S","Coenzyme A")
h = create_metabolites("h", "H", "Proton", "c")
atp = create_metabolites("atp", "C10H16N5O13P3", "Adenosine-5-triphosphate", "c")
adp = create_metabolites("adp", "C10H15N5O10P2", "Adenosine-5-diphosphate", "c")
ac_glt_p = create_metabolites("ac_glt_p","C7H12NO8P","N-Acetyl-L-glutamate 5-phosphate","c")
nadph = create_metabolites("nadph","C21H30N7O17P3","NADPH","c")
nadp = create_metabolites("nadp","C21H29N7O17P3","NADP+","c")
pi = create_metabolites("pi", "H3PO4", "Phosphate", "c")
ac_glt_sa = create_metabolites("ac_glt_sa","C7H11NO4","N-Acetyl-L-glutamate 5-semialdehyde","c")
ac_ort = create_metabolites("ac_ort","C7H14N2O3","N-Acetylornithine","c")
ort = create_metabolites("ort","C5H12N2O2","L-Ornithine","c")
h2o = create_metabolites("h2o", "H2O", "Water", "c")
acet = create_metabolites("acet","C2H4O2","Acetate","c")
cp = create_metabolites("cp","CH4NO5P","Carbamoyl phosphate","c") # Can come from pyrimidine metabolism and also be consumed by it
cit = create_metabolites("cit","C6H13N3O3","L-Citrulline","c")
amp = create_metabolites("amp", "C10H14N5O7P", "Adenosine-5-monophosphate", "c")
dpi = create_metabolites("dpi","H4P2O7","Diphosphate","c")
arg_suc = create_metabolites("arg_suc","C10H18N4O6","L-Argininosuccinate","c")
arg = create_metabolites("arg","C6H14N4O2","L-Arginine","c")
fum = create_metabolites("fum","C4H4O4","Fumarate","c") # Can be consumed by TCA cycle
gln = create_metabolites("gln","C5H10N2O3","L-Glutamine","c")
hco3 = create_metabolites("hco3","HCO3","HCO3-","c")

# --- Side Metabolites
nh3 = create_metabolites("nh3","NH3","Ammonia","c") # Can come from Nitrogen Metabolism
ac_cit = create_metabolites("ac_cit","C8H15N3O4","N-Acetyl-L-citrulline","c")

# --- Main Road 
R_2611 = create_reaction("R00355","L-aspartate:2-oxoglutarate aminotransferase","(b0928)",{asp:-1.0,oxa_glt:-1.0,oxa:1.0,glt:1.0},subsystem="Arginine biosynthesis")
R_2612 = create_reaction("R00258","L-alanine:2-oxoglutarate aminotransferase","(b2290)",{ala:-1.0,oxa_glt:-1.0,pyr:1.0,glt:1.0},subsystem="Arginine biosynthesis")
R_2311 = create_reaction("R00259","acetyl-CoA:L-glutamate N-acetyltransferase","(b2818)",{ac_coa:-1.0,glt:-1.0,ac_glt:1.0,coa:1.0,h:1.0},subsystem="Arginine biosynthesis")
R_2728 = create_reaction("R02649","ATP:N-acetyl-L-glutamate 5-phosphotransferase","(b3959)",{atp:-1.0,ac_glt:-1.0,adp:1.0,ac_glt_p:1.0},-1000,subsystem="Arginine biosynthesis")
R_12138 = create_reaction("R03443","N-acetyl-L-glutamate-5-semialdehyde:NADP+ 5-oxidoreductase (phosphrylating)","(b3958)",{nadph:-1.0,h:-1.0,ac_glt_p:-1.0,pi:1.0,nadp:1.0,ac_glt_sa:1.0},-1000,subsystem="Arginine biosynthesis")
R_26111 = create_reaction("R02283","N2-acetyl-L-ornithine:2-oxoglutarate aminotransferase","(b3359)",{glt:-1.0,ac_glt_sa:-1.0,oxa_glt:1.0,ac_ort:1.0},-1000,subsystem="")
R_35116 = create_reaction("R00669","N2-acetyl-L-ornithine amidohydrolase","(b3957)",{ac_ort:-1.0,h2o:-1.0,ort:1.0,acet:1.0},subsystem="")
R_2133 = create_reaction("R01398","carbamoyl-phosphate:L-ornithine carbamoyltransferase","(b0273 or b4254)",{ort:-1.0,cp:-1.0,cit:1.0,pi:1.0,h:1.0},-1000,subsystem="")
R_6345 = create_reaction("R01954","L-citrulline:L-aspartate ligase (AMP-forming)","(b3172)",{asp:-1.0,cit:-1.0,atp:-1.0,amp:1.0,dpi:1.0,arg_suc:1.0,h:1.0},subsystem="")
R_4321 = create_reaction("R01086","2-(Nomega-L-arginino)succinate arginine-lyase (fumarate-forming)","(b3960)",{arg_suc:-1.0,arg:1.0,fum:1.0},-1000,subsystem="")
R_6355 = create_reaction("R00575","HCO3-:L-glutamine amido-ligase (ADP-forming, carbamate-phosphorylating)","(b0032 or b0033)",{atp:-2.0,gln:-1.0,hco3:-1.0,h2o:-1.0,adp:2.0,pi:1.0,glt:1.0,cp:1.0,h:2.0},subsystem="") # Probably uses R_2722 and gets the hco3 from pyr-hco3

# Total missing:
# asp -2 | oxa_glt -1 | ala -1 | ac_coa -1
# atp -4 | nadph -1 | h2o -2 | gln -1 | hco3 -1

# Total overflow:
# glt +1 | oxa +1 | pyr +1 | coa +1
# h +4 | adp +3 | nadp +1 | pi +3
# acet +1 | amp +1 | dpi +1 | fum +1 

# --- Side Reactions (Not on MetaCyc)
R_1414 = create_reaction("R00248","L-glutamate:NADP+ oxidoreductase (deaminating)","(b1761)",{glt:-1.0,nadp:-1.0,h2o:-1.0,oxa_glt:1.0,nadph:1.0,h:1.0,nh3:1.0},-1000,subsystem="")
R_3512 = create_reaction("R00256","L-glutamine amidohydrolase","(b0485 or b1524)",{gln:-1.0,glt:1.0,nh3:1.0,h:1.0,h2o:-1.0},subsystem="")
R_2722 = create_reaction("R00150","ATP:carbamate phosphotransferase","(b0323 or b0521 or b2874)",{nh3:-1.0,hco3:-1.0,atp:-1.0,cp:1.0,adp:1.0,h2o:1.0,h:1.0},-1000,subsystem="")
R_6312 = create_reaction("R00253","L-glutamate:ammonia ligase (ADP-forming)","(b3870)",{nh3:-1.0,glt:-1.0,atp:-1.0,gln:1.0,adp:1.0,pi:1.0,h:1.0},subsystem="") #gln is produced by gln-acet
R_35116_b = create_reaction("R09107","N-acetyl-L-citrulline amidohydrolase","(b3957)",{h2o:-1.0,ac_cit:-1.0,cit:1.0,acet:1.0},subsystem="") # Doesnt work bcs of ac_cit

# Total missing:
# (s)asp -2 | ~ala -1 | ~ac_coa -1 | ~atp -6 
# (s)h2o -4 | gln -1 | ~hco3 -2 | ac_cit -1

# Total overflow:
# ~oxa +1 | ~pyr +1 | ~coa +1 | (d)h +8
# ~adp +5 | ~pi +4 | acet +2 | ~cp +1
# cit +1 | ~amp +1 | ~dpi +1 | fum +1 

# --- Balance Reactions
atp_reaction = create_reaction("ATP-M","ATP from ADP Reaction for balance","(b0000)",{adp:-1.0,pi:-1.0,atp:1.0}) # atp -2 | adp +1 | pi 0
amp_reaction = create_reaction("AMP-M","ATP from AMP Reaction for balance","(b0000)",{amp:-1.0,dpi:-1.0,atp:1.0}) # atp -1 | amp 0 | dpi 0

# --- Non-Realistic Reactions
adp_reaction = create_reaction("ADP-M","Non-accurate, here for atp balance","(b0000)",{atp:1.0,adp:-1.0}) # atp 0 | adp 0
ac_coa_reaction = create_reaction("CoA-AceCoA", "AceCoA Reaction for balance","(b0000)",{coa:-1.0,ac_coa:1.0}) # ac_coa 0 | coa 0
oxa_ala_reaction = create_reaction("Oxa-Ala","Balance Reaction for Oxa and Ala","(b0000)",{oxa:-1.0,ala:1.0}) # oxa 0 | ala 0
hco3_cp_reaction = create_reaction("Cp-HCO3","Balance for cp","(b0000)",{cp:-1.0,hco3:1.0}) # cp 0 | hco3 -1
hco3_pyr_reaction = create_reaction("Pyr-HCO3","Balance for hco3 and pyr","(b0000)",{hco3:1.0,pyr:-1.0}) # hco3 0 | pyr 0
acet_gln_reaction = create_reaction("Acet-Gln","Balance for gln and acet","(b0000)",{acet:-1.0,gln:1.0}) # acet 0/(c) +1 | gln 0

model.add_reactions([
    R_2611,
    R_2612,
    R_2311,
    R_2728,
    R_12138,
    R_26111,
    R_35116,
    R_2133,
    R_6345,
    R_4321,
    R_6355,
    R_1414,
    R_3512,
    R_2722,
    R_6312,
    R_35116_b,
    atp_reaction,
    amp_reaction,
    adp_reaction,
    ac_coa_reaction,
    oxa_ala_reaction,
    hco3_cp_reaction,
    hco3_pyr_reaction,
    acet_gln_reaction
])

model.objective = "R01086"

# --- Exchanges and Sinks
model.add_boundary(h2o, "sink")
model.add_boundary(asp, "sink", lb=-30.0) # If lower bound not set, solution becomes 166.667 . If set, solution becomes 15.000 . So 2:1

model.add_boundary(h, "demand")
model.add_boundary(arg, "demand")
model.add_boundary(fum, "demand") # Consumption by TCA

with model:
    print("---Acet-Gln Closed---")
    model.add_boundary(acet, "demand")
    model.add_boundary(glt, "sink")

    model.reactions.get_by_id("Acet-Gln").knock_out()

    solution = model.optimize()
    print(solution)
    print("\n")
    print(solution.fluxes)
    print("\n")


# with model:
#     print("---Hco3-Pyr Closed---")
#     model.add_boundary(pyr, "demand")
#     model.add_boundary(hco3, "sink")

#     model.reactions.get_by_id("Pyr-HCO3").knock_out()

#     solution = model.optimize()
#     print(solution)
#     print("\n")
#     print(solution.fluxes)
#     print("\n")

# with model:
#     print("---Hco3-Cp Closed---")
#     model.add_boundary(cp, "demand")
#     model.add_boundary(hco3, "sink")

#     model.reactions.get_by_id("Cp-HCO3").knock_out()

#     solution = model.optimize()
#     print(solution)
#     print("\n")
#     print(solution.fluxes)
#     print("\n")

# with model:
#    print("---Oxa-Ala Closed---")
#    model.add_boundary(oxa, "demand")
#    model.add_boundary(ala, "sink")

#    model.reactions.get_by_id("Oxa-Ala").knock_out()

#    solution = model.optimize()
#    print(solution)
#    print("\n")
#    print(solution.fluxes)
#    print("\n")

    

print("---Normal Model---")
solution = model.optimize()
print(solution)
print("\n")
print(solution.fluxes)
print("\n")