# L01 stage2 view: MODEL 1 poses (receptor + ligand) aligned onto 7POU/7POL.
# Both experimental proteins (7POU, 7POL) shown; all pose receptors + ligands shown. No capture.
# run:  @/path/to/casp17-local/L01/make_figures.pml

reinitialize
bg_color white
set ray_shadows, 0
set cartoon_transparency, 0.8
set stick_radius, 0.25

fetch 7pou, 7pou, async=0
fetch 7pol, 7pol, async=0

# both experimental proteins (7POL aligned onto 7POU frame)
super 7pol and polymer, 7pou and polymer
hide everything
show cartoon, 7pou
color grey70, 7pou
show cartoon, 7pol
color lightblue, 7pol

# experimental drugs (green) + catalytic zinc (purple)
select drugs, (7pou or 7pol) and resn 6JP+7X9
show sticks, drugs
set stick_radius, 0.3, drugs
color green, drugs
util.cnc drugs
select zn, 7pou and elem Zn
show spheres, zn
color purple, zn
set sphere_scale, 0.5, zn

# poses: receptor cartoon (wheat, transparent) + ligand sticks (colored)
load /path/to/casp17-local/L01/consensus/L010016/05_final/model_4.cif, m_016
super m_016 and polymer, 7pou and polymer
show cartoon, m_016
set cartoon_transparency, 0.7, m_016
color wheat, m_016
select lig_016, m_016 and hetatm and not solvent
show sticks, lig_016
color cyan, lig_016
util.cnc lig_016

load /path/to/casp17-local/L01/consensus/L010078/05_final/model_1.cif, m_078
super m_078 and polymer, 7pou and polymer
show cartoon, m_078
set cartoon_transparency, 0.7, m_078
color wheat, m_078
select lig_078, m_078 and hetatm and not solvent
show sticks, lig_078
color yellow, lig_078
util.cnc lig_078

load /path/to/casp17-local/L01/consensus/L010223/05_final/model_1.cif, m_223
super m_223 and polymer, 7pou and polymer
show cartoon, m_223
set cartoon_transparency, 0.7, m_223
color wheat, m_223
select lig_223, m_223 and hetatm and not solvent
show sticks, lig_223
color orange, lig_223
util.cnc lig_223

load /path/to/casp17-local/L01/consensus/L010695/05_final/model_1.cif, m_695
super m_695 and polymer, 7pou and polymer
show cartoon, m_695
set cartoon_transparency, 0.7, m_695
color wheat, m_695
select lig_695, m_695 and hetatm and not solvent
show sticks, lig_695
color magenta, lig_695
util.cnc lig_695

load /path/to/casp17-local/L01/consensus/L010319/05_final/model_1.cif, m_319
super m_319 and polymer, 7pou and polymer
show cartoon, m_319
set cartoon_transparency, 0.7, m_319
color wheat, m_319
select lig_319, m_319 and hetatm and not solvent
show sticks, lig_319
color red, lig_319
util.cnc lig_319

# BFT1 co-folding templates (aligned onto 7POU frame), transparent reference
load /path/to/casp17-local/templates/L01/3p24.cif, t_3p24
super t_3p24 and polymer, 7pou and polymer
hide everything, t_3p24
show cartoon, t_3p24
set cartoon_transparency, 0.8, t_3p24
color palegreen, t_3p24

load /path/to/casp17-local/L01/inputs/A1+L010016/af3_msa/template_A_01_8wen.cif, t_8wen
super t_8wen and polymer, 7pou and polymer
hide everything, t_8wen
show cartoon, t_8wen
set cartoon_transparency, 0.8, t_8wen
color lightpink, t_8wen

load /path/to/casp17-local/L01/inputs/A1+L010016/af3_msa/template_A_02_4on1.cif, t_4on1
super t_4on1 and polymer, 7pou and polymer
hide everything, t_4on1
show cartoon, t_4on1
set cartoon_transparency, 0.8, t_4on1
color paleyellow, t_4on1

load /path/to/casp17-local/L01/inputs/A1+L010016/af3_msa/template_A_03_8weo.cif, t_8weo
super t_8weo and polymer, 7pou and polymer
hide everything, t_8weo
show cartoon, t_8weo
set cartoon_transparency, 0.8, t_8weo
color lightorange, t_8weo

# template het/ligands (only 3p24 has any: azide + glycerol/PEG; 8wen/4on1/8weo are protein-only)
select tmpl_het, (t_3p24 or t_8wen or t_4on1 or t_8weo) and not polymer and not solvent and not elem Zn
show sticks, tmpl_het
set stick_radius, 0.2, tmpl_het
color grey40, tmpl_het
util.cnc tmpl_het
select tmpl_zn, (t_3p24 or t_8wen or t_4on1 or t_8weo) and elem Zn
show spheres, tmpl_zn
color purple, tmpl_zn
set sphere_scale, 0.4, tmpl_zn

deselect
orient drugs
zoom drugs, 12
