# PyMOL 一键出图：load egfr_pocket.pml（在该文件所在目录打开 PyMOL 后 @egfr_pocket.pml）
load 03_1M17_protein.pdb, crystal
load 04_1M17_ligand_AQ4.pdb, ligand
load 02_AF_superposed_on_1M17.pdb, afmodel
load 07_pocket_geo_AF_full.pdb, pocket_af
load 06_pocket_geo_crystal.pdb, pocket_cry

hide everything
# 晶体：灰色卡通；AF：按 B 因子(pLDDT)彩虹着色
show cartoon, crystal;  color grey40, crystal
show cartoon, afmodel;  spectrum b, blue_white_red, afmodel, 50, 100
# 配体：绿色球棍
show sticks, ligand;    show spheres, ligand;  color green, ligand
# 检测到的口袋：半透明球壳
show spheres, pocket_af;   set sphere_scale, 1.0, pocket_af;  color orange, pocket_af;  set transparency, 0.55, pocket_af
show spheres, pocket_cry;  set sphere_scale, 1.0, pocket_cry; color marine,  pocket_cry; set transparency, 0.55, pocket_cry
zoom ligand, 12
ray 1600, 1200; png egfr_pocket.png
