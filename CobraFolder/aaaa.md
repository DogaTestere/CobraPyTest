Genlerde ismi 'eco' olan = Escherichia coli K-12 MG1655 <- Buna bak

glyc exchange yerine

glc-p ----> gap + g2p yapsak
glukoz-p -> Pentoz fosfat yolu var zaten.

oxa üretimi de benzer şekilde yapılabilir.
pirüvat ----> asetilCoA ----> citrat ----> oxa
pirüvat ----> ThPP ----> acetyl-dihidro ----> asetilCoA ----> citrat ----> oxa

R_41119 -> b3403 | ATP + Oxaloacetate <=> ADP + Phosphoenolpyruvate + CO2 | R00341 

R_1271 -> b1378 | 2 Reduced ferredoxin + Acetyl-CoA + CO2 + 2 H+ <=> 2 Oxidized ferredoxin + Pyruvate + CoA | R01196

R_1241 -> b0114 | Pyruvate + Thiamin diphosphate <=> 2-(alpha-Hydroxyethyl)thiamine diphosphate + CO2 | R00014

R_1241_e ->  b0114 | 2-(alpha-Hydroxyethyl)thiamine diphosphate + Enzyme N6-(lipoyl)lysine <=> [Dihydrolipoyllysine-residue acetyltransferase] S-acetyldihydrolipoyllysine + Thiamin diphosphate | R03270

Yeni yolda moleküller eksik özellikle coa, ama oluşturacak yer yok. Sadece aceCoa tüketilerek oluşturuluyor. O da, pruvat harcanarak oluşuyor. Pürivat amaç olduğu sürecede de yeni eklenen yollar çalışmıyor.

---
## Yan yollar

alpha-D-glucose-6-p:
- [x] 5.4.2.2
- [ ] 5.1.3.15
- [x] 2.7.1.2
- [x] (Ters) 5.3.1.9

gliserolaldehit-3-p:
- [x] pentoz-fosfat yolu
- [x] (Ters) 4.1.2.13

glycreate-3-p:
- [ ] calvin döngüsünden ( Yok bulamıyom)
- [x] (Ters) 2.7.2.3

glycreate-2-p:
- [x] pentoz-fosfat yolundan
- [x] (Ters) 4.2.1.11

Phosphoenolpruvate:
- [x] 4.1.1.49

# checklist
(M) atp -1 atp +1 atp +1 **atp -1** atp -1 atp -1 atp -1
(P) adp +1 adp -1 adp -1 **adp +1** adp +1 adp +1 adp +1
(D) h2o -1 h2o +1 
(P) pi +1 pi -1 **pi +1**
(D) nad -1 **nad +1**
(D) nadh +1 **nadh -1**
(P) h +1 h -1 h +1 h +1
( ) co2 +1

glc -1 <-Unlimited
(P) glc-p +1 glc-p -1 glc_p +1 glp-p +1 glc-p +1
(P) frc-p +1 frc-p -1 frc-p +1
(M) frc-dp +1 frc-dp -1 frc-dp -1
(P) dhap +1 dhap +1
(D) gap +1 gap -1 gap -1 gap +1
(D) dpg +1 dpg -1
(P) g3p +1 g3p +1
(M) g2p -1 g2p -1 g2p +1
(P) pep +1 pep -1 pep +1
(P) pyr +1 pyr +1
kdpg <-Unlimited from pentose phosphate
(M) glc_d -1.0
glyc <- Unlimited from pentose phosphate
oxa <- Unlimited from citrate cycle
glc-ap -1 <- from Startch and sucrose metabolism
glc-bp -1 
