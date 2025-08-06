from cobra import Model, Reaction, Metabolite

model = Model("test_model")

glc_Intake = Reaction("EX_glc__D_e") # These are for manual adding of the reactions
glc_Intake.bounds = (-18.5, 1000)

o2_Intake = Reaction("EX_o2_e")
o2_Intake.bounds = (0, 1000)
# Metabolites and genes aren't added into this

print(f'{len(model.reactions)} reactions initially')
print(f'{len(model.metabolites)} metabolites initially')
print(f'{len(model.genes)} genes initially')

model.add_reactions([glc_Intake, o2_Intake])

print(f'{len(model.reactions)} reactions now')
print(f'{len(model.metabolites)} metabolites now')
print(f'{len(model.genes)} genes now')