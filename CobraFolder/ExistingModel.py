from cobra.io import load_model
model = load_model("textbook")

o2_Intake = model.reactions.get_by_id("EX_o2_e") 
old_bounds = model.reactions.get_by_id("EX_o2_e").bounds
print(f"Old bounds: {old_bounds}")

o2_Intake.bounds = (0, 1000)
print(f"New Bounds: {o2_Intake.bounds}")