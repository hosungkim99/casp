# L01 (BFT1) binding pocket candidates on L010016 receptor
# A=red(consensus 188-268) B=blue(p2rank1) C=green(p2rank4) E=orange(new 194-275) D=purple(p2rank2)
bg_color white
load /path/to/casp17-local/L01/outputs/L010016/05_final/score/P1_L1_prot.pdb, rec
hide everything
show cartoon, rec
color grey80, rec
set cartoon_transparency, 0.6, rec
select pA, rec and resi 188+190+192+230+232+233+254+257+259+260+263+268
show sticks, pA
show spheres, p{k} and name CA
set sphere_scale, 0.5, pA and name CA
color red, pA
pseudoatom cenA, pos=[11.3,2.2,11.9], label=A
select pB, rec and resi 105+106+109+110+142+151+153+154+273+306+308
show sticks, pB
show spheres, p{k} and name CA
set sphere_scale, 0.5, pB and name CA
color blue, pB
pseudoatom cenB, pos=[5.0,-3.1,-15.3], label=B
select pC, rec and resi 61+62+67+68+69+80+82+210+213+341+342
show sticks, pC
show spheres, p{k} and name CA
set sphere_scale, 0.5, pC and name CA
color green, pC
pseudoatom cenC, pos=[-10.6,-5.0,-1.8], label=C
select pE, rec and resi 194+195+196+197+198+199+200+272+273+274+275
show sticks, pE
show spheres, p{k} and name CA
set sphere_scale, 0.5, pE and name CA
color orange, pE
pseudoatom cenE, pos=[12.8,1.4,-3.7], label=E
select pD, rec and resi 8+33+34+36+81+83+84+86+88+90+119+120+123+141+143+146
show sticks, pD
show spheres, p{k} and name CA
set sphere_scale, 0.5, pD and name CA
color purple, pD
pseudoatom cenD, pos=[-10.6,-1.9,-16.5], label=D
set label_size, 24
set label_color, black
deselect
zoom rec
