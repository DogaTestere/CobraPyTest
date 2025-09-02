
import cobra
import pandas as pd
from cobra import io
from cobra.manipulation import knock_out_model_genes
import copy

# Modeli yükle
model = cobra.io.read_sbml_model("iML1515.xml")

print(f"Model adı: {model.name}")
print(f"Toplam reaksiyon: {len(model.reactions)}")
print(f"Toplam metabolit: {len(model.metabolites)}")
print(f"Toplam gen: {len(model.genes)}")

# Büyüme ortamı
model.medium = model.medium  
model.medium["EX_glc__D_e"] = 10.0  # glukoz uptake (mmol/gDW/h)
model.medium["EX_o2_e"] = 20.0      # oksijen uptake

# FBA çözümü
solution = model.optimize()
print(f"Maksimum büyüme hızı: {solution.objective_value:.4f} 1/h")

solution.fluxes.sort_values(ascending=False).head(10)

model.medium["EX_o2_e"] = 0.0  # anaerobik koşul
anaerobic_solution = model.optimize()
print(f"Anaerobik büyüme hızı: {anaerobic_solution.objective_value:.4f}")
print("Anaerobik koşulda en yüksek flux değerleri:")
print(anaerobic_solution.fluxes.sort_values(ascending=False).head(10))  

# Aerobik koşul
model.medium["EX_o2_e"] = 20.0
aerobic_sol = model.optimize()

# Anaerobik koşul
model.medium["EX_o2_e"] = 0.0
anaerobic_sol = model.optimize()

comparison = pd.DataFrame({
    "Koşul": ["Aerobik", "Anaerobik"],
    "Büyüme hızı (1/h)": [aerobic_sol.objective_value, anaerobic_sol.objective_value]
})

print(comparison)

# Karşılaştırma tablosu
flux_comparison = pd.DataFrame({
    "Reaksiyon": [rxn.id for rxn in model.reactions],
    "Aerobik Flux": aerobic_sol.fluxes.values,
    "Anaerobik Flux": anaerobic_sol.fluxes.values
})

# Fark sütunu ekleme (Aerobik - Anaerobik)
flux_comparison["ΔFlux"] = flux_comparison["Aerobik Flux"] - flux_comparison["Anaerobik Flux"]

# Reaksiyon ID'sini index yap
flux_comparison = flux_comparison.set_index("Reaksiyon")

# En çok değişen reaksiyonlar (mutlak değere göre sıralama)
flux_diff_sorted = flux_comparison.reindex(flux_comparison["ΔFlux"].abs().sort_values(ascending=False).index)
print(flux_diff_sorted.head(10))

#Konkout işlemi
solution_before = model.optimize()
print("Knockout öncesi büyüme hızı:", round(solution_before.objective_value, 4), "1/h")
model_copy = copy.deepcopy(model)

gene_to_knockout = ['b1779']  # Buraya istediğin genin ID'sini yazabiliriz
knocked_out_rxns = knock_out_model_genes(model_copy, gene_to_knockout)

# Knockout sonrası büyüme hızı
solution_after = model_copy.optimize()
print(f"{gene_to_knockout} gen knockout sonrası büyüme hızı:", round(solution_after.objective_value, 4), "1/h")
print("Devre dışı kalan reaksiyon sayısı:", len(knocked_out_rxns))

# devre dışı kalan reaksiyonlar
knocked_out_reactions = [rxn.id for rxn in knocked_out_rxns]
print("Devre dışı kalan reaksiyonlar:")
print(knocked_out_reactions)    

#değişen flux değerleri
flux_change = solution_after.fluxes - solution_before.fluxes
print("Knockout sonrası değişen flux değerleri:")
print(flux_change.sort_values(ascending=False).head(10))

flux_change.to_csv("flux_change_after_knockout.csv")

#Metabolit üretimi için FBA çözümü
target_rxn = "EX_succ_e"  # Dışarıya verilen sukkinat reaksiyonu
if target_rxn not in model.reactions:
    raise ValueError(f"{target_rxn} modelde bulunamadı!")
model.objective = target_rxn
solution = model.optimize()

print(f"Hedef metabolit: {target_rxn}")
print(f"Maksimum üretim miktarı: {solution.objective_value:.4f} mmol/gDW/h")
