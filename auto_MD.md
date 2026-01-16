echo 6 | gmx pdb2gmx -f 7alv-8etr.pdb -o receptor_processed.gro -water spce

gmx editconf -f receptor_processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

touch ions.mdp
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 200
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -quiet

touch minim.mdp
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 200
gmx mdrun -v -deffnm em -nb gpu

touch nvt.mdp
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro -maxwarn 200
gmx mdrun -deffnm nvt -v -nb gpu

touch npt.mdp
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -r nvt.gro -maxwarn 200
gmx mdrun -deffnm npt -v -nb gpu

touch md.mdp
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr -maxwarn 200
gmx mdrun -deffnm md_0_1 -v -nb gpu