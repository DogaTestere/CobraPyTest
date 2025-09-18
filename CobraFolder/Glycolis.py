import pandas as pd

from cobra import Model, Reaction, Metabolite

pd.set_option("display.float_format", "{:.10f}".format)

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
glyc_p = create_metabolites("dhap", "C3H7O6P", "Glycerone phosphate", "c")
gap = create_metabolites("gap", "C3H7O6P", "D-Glyceraldehyde 3-phosphate", "c")
nad = create_metabolites("nad", "C21H28N7O14P2", "NAD+", "c")
nadh = create_metabolites("nadh", "C21H29N7O14P2", "NADH", "c")
h = create_metabolites("h", "H", "Proton", "c")
gly_p = create_metabolites("gly_p", "C3H8O10P2", "3-Phospho-D-glyceroyl phosphate", "c")
g3p = create_metabolites("g3p", "C3H7O7P", "3-Phospho-D-glycerate", "c")
g2p = create_metabolites("g2p", "C3H7O7P", "2-Phospho-D-glycerate", "c")
pep = create_metabolites("pep", "C3H5O6P", "Phosphoenolpyruvate", "c")
pyr = create_metabolites("pyr", "C3H4O3", "Pyruvate", "c")
amp = create_metabolites("amp", "C10H14N5O7P", "Adenosine-5-monophosphate", "c")

# --- Extra Pathways
glc_a = create_metabolites("glc_d","C6H12O6","Alpha-D-Glucose","c")
oxa = create_metabolites("oxa","C4H4O5","Oxaloacetate","c")
co2 = create_metabolites("co2","CO2","Carbondioxide","c")
glc_ap = create_metabolites("glc_ap","C6H13O9P","Alpha-D-Glucose 1-phosphate","c")
glc_bp = create_metabolites("glc_bp","C6H11O9P","Beta-D-glucose 6-phosphate","c")
glc_b = create_metabolites("glc_b","C6H12O6","beta-D-Glucose","c")
salc_p = create_metabolites("salc_p","C13H19O10P","Salicin 6-phosphate","c")
salc_oh = create_metabolites("salc_oh","C7H8O2","Salicyl alcohol","c")
arbt_p = create_metabolites("arbt_p","C12H17O10P","Arbutin 6-phosphate","c")
hqui = create_metabolites("hqui","C6H6O2","Hydroquinone","c")
salc = create_metabolites("salc","C13H18O7","Salicin","e")
arbt = create_metabolites("arbt","C12H16O7","Arbutin","e")
ace_coa = create_metabolites("ace_coa","C23H38N7O17P3S","Acetyl-CoA")
coa = create_metabolites("coa","C21H36N7O16P3S","Coenzyme A")
thipp = create_metabolites("thipp","C12H19N4O7P2S","Thiamin diphosphate","c")
h_thipp = create_metabolites("h_thipp","C14H23N4O8P2S","2-(alpha-Hydroxyethyl)thiamine diphosphate","c")

# --- Reaction Definitions ---
R_27199 = create_reaction("R02738", "protein-N(pi)-phosphohistidine:D-glucose 6-phosphotransferase", "(b1101 or b1621)", {glc_EX: -1.0, glc_p: 1.0})
R_5319 = create_reaction("R13199", "alpha-D-glucose-6-phosphate aldose-ketose-isomerase", "(b4025)", {glc_p: -1.0, frc_p: 1.0}, -1000,subsystem="Glycolysis (Embden-Meyerhof pathway)")
R_27111 = create_reaction("R00756", "ATP:D-fructose-6-phosphate 1-phosphotransferase", "(b1723 or b3916)", {atp: -1.0, frc_p: -1.0, adp: 1.0, frc_dp: 1.0},subsystem="Glycolysis (Embden-Meyerhof pathway)")
R_31311 = create_reaction("R00762", "D-fructose-1,6-bisphosphate 1-phosphohydrolase", "(b2930 or b3925 or b4323)", {frc_dp: -1.0, h2o: -1.0, frc_p: 1.0, pi: 1.0})
R_41213 = create_reaction("R01068", "D-fructose-1,6-bisphosphate D-glyceraldehyde-3-phosphate-lyase", "(b2097 or b2925)", {frc_dp: -1.0, glyc_p: 1.0, gap: 1.0}, -1000,subsystem="Glycolysis (Embden-Meyerhof pathway)")
R_5311 = create_reaction("R01015", "D-glyceraldehyde-3-phosphate aldose-ketose-isomerase", "(b3919)", {gap: -1.0, glyc_p: 1.0}, -1000,subsystem="Glycolysis, core module involving three-carbon compounds")
R_12112 = create_reaction("R01061", "D-glyceraldehyde-3-phosphate:NAD+ oxidoreductase (phosphorylating)", "(b1779)", {gap: -1.0, pi: -1.0, nad: -1.0, gly_p: 1.0, nadh: 1.0, h: 1.0}, low_bound=-1000, subsystem="Glycolysis, core module involving three-carbon compounds")
R_2723 = create_reaction("R01512", "ATP:3-phospho-D-glycerate 1-phosphotransferase", "(b2926)", {gly_p: -1.0, adp: -1.0, g3p: 1.0, atp: 1.0}, -1000,subsystem="Glycolysis, core module involving three-carbon compounds")
R_54211 = create_reaction("R01518", "D-phosphoglycerate 2,3-phosphomutase", "(b3612 or b0755 or b4395)", {g2p: -1.0, g3p: 1.0}, -1000,subsystem="Glycolysis, core module involving three-carbon compounds")
R_42111 = create_reaction("R00658", "2-phospho-D-glycerate hydro-lyase (phosphoenolpyruvate-forming)", "(b2779)", {g2p: -1.0, pep: 1.0, h2o: 1.0}, -1000,subsystem="Glycolysis, core module involving three-carbon compounds")
R_27140 = create_reaction("R00200", "ATP:pyruvate 2-O-phosphotransferase", "(b1676 or b1854)", {pep: -1.0, adp: -1.0, h: -1.0, atp: 1.0, pyr: 1.0},subsystem="Glycolysis, core module involving three-carbon compounds")
R_2792 = create_reaction("R00199", "ATP:pyruvate,water phosphotransferase", "(b1702)", {h2o: -1.0, atp: -1.0, pyr: -1.0, amp: 1.0, pi: 1.0, pep: 1.0}) #Consider non functional

# --- From Extra Pathways
R_2712 = create_reaction("R01786","ATP:alpha-D-glucose 6-phosphotransferase","(b2388)",{atp:-1.0,glc_a:-1.0,glc_p:1.0,adp:1.0,h:1.0,},subsystem="Glycolysis (Embden-Meyerhof pathway)")
R_5422 = create_reaction("R00959","alpha-D-glucose 1-phosphate 1,6-phosphomutase","(b0688)",{glc_p:1.0,glc_ap:-1.0},-1000,subsystem="Glycolysis / Gluconeogenesis")
R_41149 = create_reaction("R00341", "ATP:oxaloacetate carboxy-lyase (transphosphorylating;phosphoenolpyruvate-forming)","(b3403)",{oxa:-1.0,atp:-1.0,pep:1.0,co2:1.0,adp:1.0},-1000,subsystem="Glycolysis / Gluconeogenesis")
R_51315 = create_reaction("R02739","alpha-D-glucose 6-phosphate ketol-isomerase","(b1780)",{glc_bp:-1.0,glc_p:1.0},-1000,subsystem="Glycolysis / Gluconeogenesis")
R_5133 = create_reaction("R01602","D-glucose 1-epimerase","(b0756 or b3879)",{glc_b:-1.0,glc_a:1.0},-1000,subsystem="Glycolysis / Gluconeogenesis")
R_2712_b = create_reaction("R01600","ATP:beta-D-glucose 6-phosphotransferase","(b2388)",{atp:-1.0,adp:1.0,glc_b:-1.0,glc_bp:1.0},subsystem="Glycolysis / Gluconeogenesis")
R_31310 = create_reaction("R00947","alpha-D-glucose-1-phosphate phosphohydrolase","(b1002)",{h2o:-1.0,glc_ap:-1.0,glc_a:1.0,pi:1.0},subsystem="Glycolysis / Gluconeogenesis")
R_32186 = create_reaction("R05134","salicin 6-phosphate glucohydrolase","(b1734 or b2716 or b2901 or b3721)",{h2o:-1.0,glc_bp:1.0,salc_p:-1.0,salc_oh:1.0},subsystem="Glycolysis / Gluconeogenesis")
R_32186_a = create_reaction("R05133","arbutin 6-phosphate glucohydrolase","(b1734 or b2716 or b2901 or b3721)",{h2o:-1.0,glc_bp:1.0,arbt_p:-1.0,hqui:1.0},subsystem="Glycolysis / Gluconeogenesis")

R_1271 = create_reaction("R01196","pyruvate:ferredoxin 2-oxidoreductase (CoA-acetylating)","(b1378)",{pyr:-1.0,ace_coa:1.0,coa:-1.0,co2:1.0,h:-2.0},-1000,subsystem="Glycolysis / Gluconeogenesis")
R_1241 = create_reaction("R00014","pyruvate:thiamin diphosphate acetaldehydetransferase (decarboxylating)","(b0114)",{pyr:-1.0,h:-1.0,co2:1.0,thipp:-1.0,h_thipp:1.0})
R_1241_e = create_reaction("R03270","2-(alpha-Hydroxyethyl)thiamine diphosphate + Enzyme N6-(lipoyl)lysine <=> [Dihydrolipoyllysine-residue acetyltransferase] S-acetyldihydrolipoyllysine + Thiamin diphosphate","(b0114)",{h_thipp:-1.0,thipp:1.0},subsystem="Glycolysis / Gluconeogenesis")
R_23112 = create_reaction("R02569","acetyl-CoA:enzyme N6-(dihydrolipoyl)lysine S-acetyltransferase","(b0115)",{coa:-1.0,ace_coa:1.0})

# --- From External Pathways
R_271165 = create_reaction("R08572", "Pentose Way Simulation","(b0514 or b3124 or b1850)",{glc_p:-1.0,gap:1.0,g2p:1.0,atp:-1.0,adp:1.0},subsystem="Pentose Phosphate Pathway Simulation")
R_2331 = create_reaction("R00351", "TCA Cycle Simulation(Not Accurate)","(b0720)",{oxa:1.0,ace_coa:-1.0},-1000,subsystem="TCA Cycle Simulation")

# --- Needed Reactions
atp_reaction = create_reaction("ATP_M", "ATP/ADP Generation for balance", "(b0000)", {atp: 1.0, adp: -1.0, pi: -1.0})
adp_reaction = create_reaction("ADP_M", "ADP generation for balance","(b0000)",{atp:-1.0,adp:1.0,pi:1.0})
nadh_oxidation = create_reaction("R07618", "enzyme N6-(dihydrolipoyl)lysine:NAD+ oxidoreductase", "(b0116)", {nadh: -1.0, nad: 1.0, h:1.0},-1000)

# --- Reaction Adding
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
    # Extra Pathways
    R_2712,
    R_5422,
    R_41149,
    R_51315,
    R_2712_b,
    R_31310,
    R_32186,
    R_32186_a,
    R_1271,
    R_1241,
    R_1241_e,
    R_23112,
    # External Pathways
    R_271165,
    R_2331,
    # Enegry Circulation
    #atp_reaction,
    #adp_reaction,
    nadh_oxidation
])

model.add_boundary(glc_EX, type="exchange", lb=-10.0)
model.add_boundary(salc, type="exchange", lb=-10.0)
model.add_boundary(arbt, type="exchange", lb=-10.0)
model.add_boundary(pyr, type="demand", lb=0.0) # Added incase pyr amount stops the other ways

model.add_boundary(pi, type="demand", lb=0.0)

# Other Metabolic Pathways
# model.add_boundary(model.metabolites.get_by_id("glc_ap"), type="sink") # Comes from starch and sucrose metabolism
# Model uses this directly instead of importing molecules if open

model.objective = "R00200"

with model:
    model.add_boundary(h,type="demand", lb=0.0)
    model.add_reactions([adp_reaction]) # Becomes zero if not added
    print("---PPP way shut down---")
    model.reactions.R08572.knock_out()
    solution = model.optimize()
    print(solution)
    print("\n")
    print(solution.fluxes)
    print("\n")

with model:
    print("---TCA way shut down---")
    model.reactions.R00351.knock_out()
    solution = model.optimize()
    print(solution)
    print("\n")
    print(solution.fluxes)
    print("\n")

print("---Normal Model---")
solution = model.optimize()
print(solution)
print("\n")
print(solution.fluxes)
print("\n")