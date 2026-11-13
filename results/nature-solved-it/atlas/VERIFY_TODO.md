# What the checker must open

Written 2026-10-19 by the agent that built the atlas. A second agent that has not
seen README.md does this list. Expect some rows to fail; report them, do not fix
them silently.

## What to check in every row of atlas.csv

1. Open the citation in column 4 and confirm the number in column 3 is in that
   paper, under that condition. A number that is in the body and not the abstract
   counts only if you opened the body.
2. Open the phylogeny paper in column 5 and confirm the count of separate origins.
3. Read the UniProt accession in column 6 today and confirm the protein and the
   species match the row.
4. Open the citation in column 7 and confirm it shows the molecule acting in a
   human cell, on a human receptor, or in a human assay.
5. Open the agency page in column 8 and confirm the approval and the agency.
6. Recompute the five funnel counts and the negative-row count from atlas.csv with
   scripts/build_tolerance_atlas.py, and compare them with counts.json and the
   numbers in README.md.
7. Check that every row marked 'none found' really has at least three queries in
   queries.csv. The atlas names the rows where it has fewer.

## Identifiers, row by row

- **D01 Snakebite envenoming [toxin]** (Didelphis aurita): https://doi.org/10.1371/journal.pntd.0005829.; PMID 28759578; https://doi.org/10.3390/ijms26157510.; PMID 40806638; UniProt Q8MIS3
- **D02 Scorpion sting envenoming [toxin]** (Onychomys torridus (southern grasshopper mouse)): https://doi.org/10.1126/science.1236451.; PMID 24159039
- **D03 Cardiac glycoside poisoning (digoxin and plant or toad cardenolides) [toxin]** (Danaus plexippus): https://doi.org/10.1007/BF02028512.; PMID 24227116; https://doi.org/10.1073/pnas.1202111109.; PMID 22826239; UniProt Q2XN70
- **D04 Tetrodotoxin poisoning [toxin]** (none found with a number): https://doi.org/10.1007/s10886-005-1345-x,; PMID 15856788; https://doi.org/10.1038/s41437-020-0300-5,; PMID 32111979; https://doi.org/10.1007/s10886-005-1345-x.
- **D05 Paralytic shellfish poisoning (saxitoxin) [toxin]** (Mya arenaria): https://doi.org/10.1038/nature03415.; PMID 15815630; https://doi.org/10.1016/j.aquatox.2018.06.017.; PMID 30075308
- **D06 Ciguatera poisoning (ciguatoxin) [toxin]** (none found): no identifier to open
- **D07 Nicotine poisoning [toxin]** (none found with a tolerance number): https://doi.org/10.1073/pnas.1314848111.; PMID 24379363
- **D08 Ethanol poisoning and alcohol-related liver injury [toxin]** (Drosophila melanogaster): https://doi.org/10.1111/j.1558-5646.1991.tb04378.x.; PMID 28564174; UniProt P00334
- **D09 Cadmium poisoning [toxin]** (none found with a cadmium number): https://doi.org/10.1002/tox.20074.; PMID 15712333
- **D10 Arsenic poisoning [toxin]** (none found): https://doi.org/10.1016/j.mrgentox.2013.08.008.; PMID 24004878
- **D11 Methylmercury poisoning [toxin]** (none found): https://doi.org/10.1186/s13028-021-00607-w; PMID 35086557; https://doi.org/10.1016/j.scitotenv.2019.136325; PMID 31927288
- **D12 Lead poisoning [toxin]** (none found): https://doi.org/10.1002/etc.1762; PMID 22278886; https://doi.org/10.1897/08-359.1; PMID 19379016
- **D13 Cyanide poisoning [toxin]** (none found): https://doi.org/10.1038/s41598-018-29148-5; PMID 30018390; https://doi.org/10.1093/molbev/msag197; PMID 42572112; UniProt A0A345BJE4; PMID 28629163
- **D14 Hydrogen sulfide poisoning [toxin]** (Poecilia mexicana): https://doi.org/10.1073/pnas.2004223117; PMID 32586956; https://doi.org/10.1111/j.1558-5646.2011.01298.x; PMID 21790570; UniProt A0A023IMR0
- **D15 Carbon monoxide poisoning [toxin]** (Mirounga angustirostris): https://doi.org/10.1242/jeb.100677; PMID 24829326
- **D16 Organophosphate poisoning [toxin]** (Culex quinquefasciatus): https://doi.org/10.1093/jmedent/42.5.815; PMID 16363165; https://doi.org/10.1111/j.1365-2583.2004.00452.x; PMID 14728661; https://doi.org/10.1098/rspb.2002.2122; PMID 12396499; UniProt Q867X2
- **D17 Hyperammonaemia (urea cycle disorders and hepatic encephalopathy) [metabolite that accumulates]** (Periophthalmodon schlosseri): https://doi.org/10.1242/jeb.00788; PMID 14747411; https://doi.org/10.1086/374281; PMID 12794674; UniProt A0A191B6Q3; https://doi.org/10.1242/jeb.203.15.2297; PMID 10887068
- **D18 Uraemia in kidney failure [metabolite that accumulates]** (Squalus acanthias): https://doi.org/10.1016/s0166-445x(00)00180-6; PMID 11489311; UniProt Q9I9F0
- **D19 Type 2 diabetes (chronic high blood glucose) [metabolite that accumulates]** (Calypte anna): https://doi.org/10.1016/s1095-6433(98)10039-9; PMID 9787825; UniProt A0A7H1CRJ6; https://doi.org/10.1242/jeb.229989; PMID 32895327
- **D20 Unconjugated hyperbilirubinaemia (jaundice) [metabolite that accumulates]** (none found): no identifier to open
- **D21 Hyperuricaemia and gout [metabolite that accumulates]** (none found): https://doi.org/10.1073/pnas.1320393111; PMID 24550457
- **D22 HSAN1 neuropathy (1-deoxysphingolipid accumulation) [metabolite that accumulates]** (none found): https://doi.org/10.1016/j.ejphar.2008.02.011; PMID 18343365
- **D23 Iron overload (haemochromatosis) [metabolite that accumulates]** (Desmodus rotundus): https://doi.org/10.1016/0300-9629(82)90178-5; PMID 6128125; https://doi.org/10.1126/sciadv.abm6494; PMID 35333583; https://doi.org/10.1007/s00360-018-1155-4; PMID 29594459
- **D24 Copper overload (Wilson disease) [metabolite that accumulates]** (none found): no identifier to open
- **D25 Fatty liver disease (hepatic steatosis) [metabolite that accumulates]** (Anser anser): https://doi.org/10.1016/j.psj.2025.106052; PMID 41205434; https://doi.org/10.1016/j.psj.2025.105515; PMID 40651204
- **D26 Lactic acidosis [metabolite that accumulates]** (Chrysemys picta bellii): https://doi.org/10.1086/316754; PMID 10893168; https://doi.org/10.1152/ajpregu.00519.2006; PMID 17008457
- **D27 Hypercapnia (carbon dioxide retention and respiratory acidosis) [metabolite that accumulates]** (Heterocephalus glaber): https://doi.org/10.1152/ajpregu.00251.2019; PMID 31967860; https://doi.org/10.1007/s00359-009-0420-0; PMID 19219614; https://doi.org/10.1016/j.ejphar.2016.09.006; PMID 27609608
- **D28 Protein aggregation neurodegeneration [metabolite that accumulates]** (Heterocephalus glaber): https://doi.org/10.1093/gerona/glaa069; PMID 32515471
- **D29 Hyperoxaluria and oxalate kidney stones [metabolite that accumulates]** (Neotoma albigula): https://doi.org/10.1128/AEM.00216-16; PMID 26896138; https://doi.org/10.3109/21681805.2014.996251; PMID 25645091
- **D30 Hypercholesterolaemia and atherosclerosis [metabolite that accumulates]** (Columba livia): https://doi.org/10.1016/j.yexmp.2004.08.010; PMID 15596055; https://doi.org/10.1093/jhered/92.5.439; PMID 11773255
- **D31 Ischaemic stroke [loss of oxygen]** (Heterocephalus glaber): https://doi.org/10.1126/science.aab3896; PMID 28428423; https://doi.org/10.1007/s00239-018-9870-8; PMID 30327830
- **D31 Ischaemic stroke [loss of oxygen]** (Chrysemys picta bellii): https://doi.org/10.1242/jeb.200.17.2295; PMID 9320212
- **D32 Myocardial infarction (ischaemic heart injury) [loss of oxygen]** (Urocitellus parryii): https://doi.org/10.1097/ALN.0000000000001113; PMID 27187119; https://doi.org/10.1007/s00239-018-9870-8; PMID 30327830; https://doi.org/10.1186/s12872-023-03203-0; PMID 36991356
- **D33 Ischaemic acute kidney injury [loss of oxygen]** (none found): no identifier to open
- **D34 Neonatal hypoxic-ischaemic encephalopathy [loss of oxygen]** (Heterocephalus glaber): https://doi.org/10.1371/journal.pone.0031568; PMID 22363676; https://doi.org/10.1007/s00239-018-9870-8; PMID 30327830; https://doi.org/10.1016/j.neulet.2011.11.042; PMID 22155615
- **D35 Brain injury after cardiac arrest [loss of oxygen]** (Urocitellus parryii): https://doi.org/10.1111/jnc.13996; PMID 28222226
- **D36 Chronic mountain sickness (high-altitude hypoxia) [loss of oxygen]** (Anser indicus): https://doi.org/10.1371/journal.pone.0094015; PMID 24710001; https://doi.org/10.1093/molbev/msq205; PMID 20685719
- **D37 Hypoxic pulmonary hypertension [loss of oxygen]** (Lama glama): https://doi.org/10.1113/JP274340; PMID 29369354
- **D38 Obstructive sleep apnoea (intermittent hypoxia) [loss of oxygen]** (Mirounga angustirostris): https://doi.org/10.1242/jeb.018887; PMID 18840667; https://doi.org/10.1152/japplphysiol.00846.2019; PMID 32324472
- **D39 Ischaemia-reperfusion injury in organ transplantation [loss of oxygen]** (Ictidomys tridecemlineatus): https://doi.org/10.1097/TP.0000000000000937; PMID 26457601
- **D40 Haemorrhagic shock [loss of oxygen]** (none found): no identifier to open
- **D41 Frostbite [freezing]** (Tenebrio molitor): https://doi.org/10.1046/j.1432-1327.2000.01734.x; PMID 11029589; https://doi.org/10.1073/pnas.94.8.3817; PMID 9108061; UniProt O16119; https://doi.org/10.1016/s0014-5793(02)03355-0; PMID 12372611; https://doi.org/10.1371/journal.pone.0116562; PMID 25714402; https://doi.org/10.1021/acsomega.5c11565; PMID 41939373
- **D42 Freezing injury in cell and organ cryopreservation [freezing]** (Rana sylvatica (also written Lithobates sylvaticus)): https://doi.org/10.1242/jeb.181.1.245; PMID 8409827
- **D42 Freezing injury in cell and organ cryopreservation [freezing]** (Flavobacterium frigoris): https://doi.org/10.1093/humrep/dev170; PMID 26202918; https://doi.org/10.1073/pnas.94.8.3817; PMID 9108061
- **D43 Accidental hypothermia [freezing]** (Urocitellus parryii): https://doi.org/10.1086/679591; PMID 25590595; https://doi.org/10.1016/j.cell.2018.03.010; PMID 29576452
- **D44 Dry eye disease [drying]** (none found): no identifier to open
- **D45 Severe dehydration [drying]** (Artemia franciscana): https://doi.org/10.1016/j.cryobiol.2005.04.007; PMID 15963489; UniProt O44112; https://doi.org/10.1074/jbc.C900001200; PMID 19228698
- **D46 Acute radiation syndrome [radiation]** (Ramazzottius varieornatus): https://doi.org/10.1038/ncomms12808; PMID 27649274; https://doi.org/10.1016/j.dnarep.2010.08.003; PMID 20817622; UniProt P0DOW4; https://doi.org/10.3390/biology10100970; PMID 34681069; https://doi.org/10.1016/j.mcn.2023.103826; PMID 36858083
- **D47 Ultraviolet skin damage (sunburn) [radiation]** (Deinococcus radiodurans): https://doi.org/10.1016/j.lssr.2024.01.006; PMID 38670653; UniProt Q9RXY7; https://doi.org/10.18632/oncotarget.13893; PMID 27974687; PMID 34681069
- **D48 Disuse osteoporosis (bone loss during immobility) [long immobility]** (Ursus arctos horribilis): https://doi.org/10.1016/j.bone.2007.10.010; PMID 18037367; https://doi.org/10.1016/j.bone.2011.08.017; PMID 21893223
- **D49 Disuse muscle atrophy [long immobility]** (Ursus americanus): https://doi.org/10.1086/515429; PMID 9678502; https://doi.org/10.1371/journal.pone.0263085; PMID 35077510
- **D50 Venous thromboembolism during immobility [long immobility]** (none found): UniProt P50454; https://doi.org/10.1126/science.abo5044; PMID 37053338
- **D51 Limb loss (no regrowth after amputation) [injury that scars instead of regrowing]** (Ambystoma mexicanum): https://doi.org/10.1002/reg2.8.; PMID 27499857
- **D52 Heart scarring after myocardial infarction [injury that scars instead of regrowing]** (Danio rerio): https://doi.org/10.1126/science.1077857.; PMID 12481136; UniProt Q8AYG3
- **D53 Spinal cord injury [injury that scars instead of regrowing]** (Danio rerio): https://doi.org/10.1126/science.aaf2679.; PMID 27811277; UniProt Q5RI33
- **D54 Skin scarring after full-thickness wounds [injury that scars instead of regrowing]** (Acomys spp. (two wild-caught Kenyan spiny mouse species; the PubMed full-text tool strips italic species names, so the two binomials were not read today)): https://doi.org/10.1038/nature11499.; PMID 23018966
- **D55 Liver fibrosis and cirrhosis [injury that scars instead of regrowing]** (none found): no identifier to open
- **D56 Sensorineural hearing loss (hair cell loss) [injury that scars instead of regrowing]** (Gallus gallus (domestic chicken; the abstract says 'chick', MeSH term Chickens)): https://doi.org/10.1016/0378-5955(87)90135-3.; PMID 3680064
- **D57 Retinal degeneration (photoreceptor loss) [injury that scars instead of regrowing]** (Danio rerio): https://doi.org/10.1002/1097-4695(20000905)44:3<289::aid-neu1>3.0.co;2-h.; PMID 10942883
- **D58 Tooth loss [injury that scars instead of regrowing]** (Alligator mississippiensis): https://doi.org/10.1002/aja.1001870407.; PMID 2353674; PMID 23671090; https://doi.org/10.1073/pnas.1213202110)
- **D59 Pulmonary fibrosis [injury that scars instead of regrowing]** (Acomys cahirinus): https://doi.org/10.1186/s41232-025-00391-4.; PMID 40883787
- **D60 Cancer in a large or long-lived body (Peto's paradox) [cancer in a large or long-lived body]** (Loxodonta africana (TP53 copy count); the cancer mortality is from zoo necropsy data on elephants, n = 644, African and Asian not split in the abstract): https://doi.org/10.1001/jama.2015.13134.; PMID 26447779; https://doi.org/10.7554/eLife.11994.; PMID 27642012

## Counts to recompute

```json
{
  "in_scope": 60,
  "named_species_with_a_number": 42,
  "molecule_named": 29,
  "shown_on_a_human_target": 7,
  "approved_drug": 0
}
```

Negative rows (diseases with no species and no number): 18.


## Already spot-checked by the agent that built the atlas

Three identifiers from the blocks collected in parallel were opened again with the PubMed
tool and matched their rows: PMID 42572112 (cyanogenic glucoside biosynthesis in
Heliconiini, which carries the independent-origins claim in the cyanide row), PMID
16363165 (chlorpyrifos resistance in Culex quinquefasciatus) and PMID 27649274 (the
tardigrade genome paper with the Dsup result in human cells). Everything else in the list
above is unchecked.
