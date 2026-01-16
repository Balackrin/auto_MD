import argparse
import subprocess
import os
import sys
from .analysis import GMXAnalysis

def run_command(cmd, input_text=None, verbose=True):
    """Execute command line command"""
    if verbose:
        print(f"Executing command: {' '.join(cmd)}")
    
    if input_text:
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=input_text)
    else:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()
    
    if verbose:
        if stdout:
            print(stdout)
        if stderr:
            print(stderr)
    
    if process.returncode != 0:
        print(f"Command execution failed: {' '.join(cmd)}")
        print(f"Error message: {stderr}")
        sys.exit(1)
    
    return stdout

def get_mdp_path(filename, custom_mdp_dir=None, simulation_type='protein'):
    """Get the path of mdp file"""
    if custom_mdp_dir:
        return os.path.join(custom_mdp_dir, filename)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if simulation_type == 'protein':
        return os.path.join(current_dir, 'protein_mdp', filename)
    else:  # protein_ligand
        return os.path.join(current_dir, 'protein_ligand_mdp', filename)

def modify_mdp_file(file_path, steps=None, dt=None):
    """Modify mdp file steps and step size"""
    if not os.path.exists(file_path):
        print(f"Error: MDP file not found at {file_path}")
        sys.exit(1)
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Modify steps if provided
    if steps is not None:
        import re
        content = re.sub(r'nsteps\s*=\s*\d+', f'nsteps      = {steps}', content)
    
    # Modify dt if provided
    if dt is not None:
        import re
        content = re.sub(r'defdt\s*=\s*\d+\.?\d*', f'defdt        = {dt}', content)
    
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"Modified mdp file: {file_path}")
    if steps is not None:
        print(f"  Updated nsteps to {steps}")
    if dt is not None:
        print(f"  Updated defdt to {dt}")

def simulate_protein(pdb_file, custom_mdp_dir=None, steps=None, dt=None):
    """Molecular dynamics simulation of single protein"""
    print("=== Starting Protein Molecular Dynamics Simulation ===")
    
    # Modify md.mdp file if steps or dt is provided
    if steps is not None or dt is not None:
        print("\nModifying md.mdp file...")
        mdp_path = get_mdp_path('md.mdp', custom_mdp_dir, simulation_type='protein')
        modify_mdp_file(mdp_path, steps, dt)
    
    # Build protein topology
    print("\n1. Building protein topology...")
    run_command(['gmx', 'pdb2gmx', '-f', pdb_file, '-o', 'receptor_processed.gro', '-water', 'spce', '-ignh'], 
                input_text='6\n')
    
    # Add box
    print("\n2. Adding box...")
    run_command(['gmx', 'editconf', '-f', 'receptor_processed.gro', '-o', 'newbox.gro', '-c', '-d', '1.0', '-bt', 'cubic'])
    
    # Add solvent
    print("\n3. Adding solvent...")
    run_command(['gmx', 'solvate', '-cp', 'newbox.gro', '-cs', 'spc216.gro', '-o', 'solv.gro', '-p', 'topol.top'])
    
    # Neutralize charge
    print("\n4. Neutralizing charge...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('ions.mdp', custom_mdp_dir, simulation_type='protein'), '-c', 'solv.gro', '-p', 'topol.top', '-o', 'ions.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'genion', '-s', 'ions.tpr', '-o', 'solv_ions.gro', '-p', 'topol.top', '-pname', 'NA', '-nname', 'CL', '-neutral'], 
                input_text='13\n')
    
    # Energy minimization
    print("\n5. Performing energy minimization...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('minim.mdp', custom_mdp_dir, simulation_type='protein'), '-c', 'solv_ions.gro', '-p', 'topol.top', '-o', 'em.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-v', '-deffnm', 'em', '-nb', 'gpu'])
    
    # NVT equilibration
    print("\n6. Performing NVT equilibration...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('nvt.mdp', custom_mdp_dir, simulation_type='protein'), '-c', 'em.gro', '-p', 'topol.top', '-o', 'nvt.tpr', '-r', 'em.gro', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-deffnm', 'nvt', '-v', '-nb', 'gpu'])
    
    # NPT equilibration
    print("\n7. Performing NPT equilibration...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('npt.mdp', custom_mdp_dir, simulation_type='protein'), '-c', 'nvt.gro', '-t', 'nvt.cpt', '-p', 'topol.top', '-o', 'npt.tpr', '-r', 'nvt.gro', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-deffnm', 'npt', '-v', '-nb', 'gpu'])
    
    # Run MD simulation
    print("\n8. Running MD simulation...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('md.mdp', custom_mdp_dir, simulation_type='protein'), '-c', 'npt.gro', '-t', 'npt.cpt', '-p', 'topol.top', '-o', 'md_0_1.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-deffnm', 'md_0_1', '-v', '-nb', 'gpu'])
    
    # Process trajectory
    print("\n9. Processing trajectory...")
    run_command(['gmx', 'trjconv', '-s', 'md_0_1.tpr', '-f', 'md_0_1.xtc', '-o', 'md_0_1_noPBC.xtc', '-pbc', 'mol', '-ur', 'compact'], 
                input_text='0\n')
    
    # Copy files to match analysis script expectations
    import shutil
    if os.path.exists('md_0_1.tpr') and not os.path.exists('md.tpr'):
        shutil.copy('md_0_1.tpr', 'md.tpr')
    if os.path.exists('md_0_1_noPBC.xtc') and not os.path.exists('md_noPBC.xtc'):
        shutil.copy('md_0_1_noPBC.xtc', 'md_noPBC.xtc')
    
    # Perform analysis
    print("\n10. Performing trajectory analysis...")
    gmx_analysis = GMXAnalysis(work_dir='.', output_dir='analysis_results')
    gmx_analysis.analyze_rmsd()
    gmx_analysis.analyze_rmsf()
    gmx_analysis.analyze_rg()
    gmx_analysis.analyze_sasa()
    gmx_analysis.analyze_hbonds()
    
    print("\n=== Protein Molecular Dynamics Simulation Completed ===")

def simulate_protein_ligand(protein_file, ligand_file, custom_mdp_dir=None, steps=None, dt=None):
    """Molecular dynamics simulation of protein-ligand complex"""
    print("=== Starting Protein-Ligand Complex Molecular Dynamics Simulation ===")
    
    # Modify md.mdp file if steps or dt is provided
    if steps is not None or dt is not None:
        print("\nModifying md.mdp file...")
        mdp_path = get_mdp_path('md.mdp', custom_mdp_dir, simulation_type='protein_ligand')
        modify_mdp_file(mdp_path, steps, dt)
    
    # 1. Generate ligand topology using acpype
    print("\n1. Generating ligand topology with acpype...")
    run_command(['acpype', '-i', ligand_file, '-n', '0', '-a', 'gaff'])
    
    # 2. Copy ligand files
    print("\n2. Copying ligand files...")
    import shutil
    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
    acpype_dir = f"{ligand_name}.acpype"
    
    # Copy ligand_GMX.gro and ligand_GMX.itp
    shutil.copy(os.path.join(acpype_dir, f"{ligand_name}_GMX.gro"), f"{ligand_name}_GMX.gro")
    shutil.copy(os.path.join(acpype_dir, f"{ligand_name}_GMX.itp"), f"{ligand_name}_GMX.itp")
    
    # 3. Build protein topology
    print("\n3. Building protein topology...")
    run_command(['gmx', 'pdb2gmx', '-f', protein_file, '-o', 'receptor_processed.gro', '-water', 'spce', '-ignh'], 
                input_text='6\n')
    
    # 4. Combine protein and ligand into complex.gro
    print("\n4. Combining protein and ligand into complex.gro...")
    
    # Copy protein coordinates
    shutil.copy('receptor_processed.gro', 'complex.gro')
    
    # Read ligand coordinates
    with open(f"{ligand_name}_GMX.gro", 'r') as f:
        ligand_lines = f.readlines()
    
    # Read complex coordinates
    with open('complex.gro', 'r') as f:
        complex_lines = f.readlines()
    
    # Update atom count in complex.gro
    complex_atom_count = int(complex_lines[1].strip())
    ligand_atom_count = int(ligand_lines[1].strip())
    total_atoms = complex_atom_count + ligand_atom_count
    complex_lines[1] = f"{total_atoms:5d}\n"
    
    # Add ligand coordinates to complex.gro
    ligand_coords = ligand_lines[2:-1]  # Third line to second last line
    complex_lines = complex_lines[:-1] + ligand_coords + complex_lines[-1:]
    
    # Write updated complex.gro
    with open('complex.gro', 'w') as f:
        f.writelines(complex_lines)
    
    # 5. Update topol.top with ligand information
    print("\n5. Updating topol.top with ligand information...")
    
    # Read topol.top
    with open('topol.top', 'r') as f:
        topol_lines = f.readlines()
    
    # Add ligand topology include
    insert_index = None
    for i, line in enumerate(topol_lines):
        if '#include "amber99sb-ildn.ff/forcefield.itp"' in line:
            insert_index = i + 1
            break
    
    if insert_index is not None:
        ligand_include = [
            '; Include ligand topology\n',
            f'#include "{ligand_name}_GMX.itp"\n',
            '\n',
            '; Icariin position restraints\n',
            '#ifdef POSRES\n',
            '#include "posre_lig.itp"\n',
            '#endif\n'
        ]
        topol_lines = topol_lines[:insert_index] + ligand_include + topol_lines[insert_index:]
    
    # Add ligand molecule count
    topol_lines.append(f'\n{ligand_name.upper():<20}1\n')
    
    # Write updated topol.top
    with open('topol.top', 'w') as f:
        f.writelines(topol_lines)
    
    # 6. Add simulation box
    print("\n6. Adding simulation box...")
    run_command(['gmx', 'editconf', '-f', 'complex.gro', '-o', 'newbox.gro', '-c', '-d', '1.0', '-bt', 'cubic'])
    
    # 7. Add solvent
    print("\n7. Adding solvent...")
    run_command(['gmx', 'solvate', '-cp', 'newbox.gro', '-cs', 'spc216.gro', '-o', 'solv.gro', '-p', 'topol.top'])
    
    # 8. Neutralize charge
    print("\n8. Neutralizing charge...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('ions.mdp', custom_mdp_dir, simulation_type='protein_ligand'), '-c', 'solv.gro', '-p', 'topol.top', '-o', 'ions.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'genion', '-s', 'ions.tpr', '-o', 'solv_ions.gro', '-p', 'topol.top', '-pname', 'NA', '-nname', 'CL', '-neutral', '-rmin', '0.5'], 
                input_text='15\n')
    
    # 9. Energy minimization
    print("\n9. Performing energy minimization...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('em.mdp', custom_mdp_dir, simulation_type='protein_ligand'), '-c', 'solv_ions.gro', '-p', 'topol.top', '-o', 'em.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-v', '-deffnm', 'em', '-nb', 'gpu'])
    
    # 10. Generate position restraints for ligand
    print("\n10. Generating position restraints for ligand...")
    run_command(['gmx', 'genrestr', '-f', f"{ligand_name}_GMX.gro", '-o', 'posre_lig.itp', '-fc', '1000', '1000', '1000'], 
                input_text='2\n')
    
    # 11. Create index groups
    print("\n11. Creating index groups...")
    run_command(['gmx', 'make_ndx', '-f', 'em.gro', '-o', 'index.ndx'], 
                input_text='1|13\nq\n')
    
    # 12. NVT equilibration
    print("\n12. Performing NVT equilibration...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('nvt.mdp', custom_mdp_dir, simulation_type='protein_ligand'), '-c', 'em.gro', '-r', 'em.gro', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'nvt.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-v', '-deffnm', 'nvt', '-nb', 'gpu'])
    
    # 13. NPT equilibration
    print("\n13. Performing NPT equilibration...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('npt.mdp', custom_mdp_dir, simulation_type='protein_ligand'), '-c', 'nvt.gro', '-r', 'nvt.gro', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'npt.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-deffnm', 'npt', '-v', '-nb', 'gpu'])
    
    # 14. Run MD simulation
    print("\n14. Running MD simulation...")
    run_command(['gmx', 'grompp', '-f', get_mdp_path('md.mdp', custom_mdp_dir, simulation_type='protein_ligand'), '-c', 'npt.gro', '-t', 'npt.cpt', '-p', 'topol.top', '-n', 'index.ndx', '-o', 'md_0_1.tpr', '-maxwarn', '200'])
    run_command(['gmx', 'mdrun', '-deffnm', 'md_0_1', '-v', '-nb', 'gpu'])
    
    # 15. Process trajectory for analysis
    print("\n15. Processing trajectory for analysis...")
    run_command(['gmx', 'trjconv', '-s', 'md_0_1.tpr', '-f', 'md_0_1.xtc', '-o', 'md_0_1_noPBC.xtc', '-pbc', 'mol', '-ur', 'compact'],
                input_text='0\n')
    
    # Copy files to match analysis script expectations
    import shutil
    if os.path.exists('md_0_1.tpr') and not os.path.exists('md.tpr'):
        shutil.copy('md_0_1.tpr', 'md.tpr')
    if os.path.exists('md_0_1_noPBC.xtc') and not os.path.exists('md_noPBC.xtc'):
        shutil.copy('md_0_1_noPBC.xtc', 'md_noPBC.xtc')
    
    # Perform analysis
    print("\n16. Performing trajectory analysis...")
    gmx_analysis = GMXAnalysis(work_dir='.', output_dir='analysis_results')
    gmx_analysis.analyze_rmsd()
    gmx_analysis.analyze_rmsf()
    gmx_analysis.analyze_rg()
    gmx_analysis.analyze_sasa()
    gmx_analysis.analyze_hbonds()
    
    print("\n=== Protein-Ligand Complex Molecular Dynamics Simulation Completed ===")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='autoMD - Automated Molecular Dynamics Simulation Tool')
    parser.add_argument('--pdb', help='Protein PDB file path (for single protein simulation)')
    parser.add_argument('--protein', help='Protein PDB file path (for protein-ligand complex simulation)')
    parser.add_argument('--ligand', help='Ligand PDB file path (for protein-ligand complex simulation)')
    parser.add_argument('--mdp_files', help='Custom mdp files directory path')
    parser.add_argument('--steps', type=int, help='Number of steps for MD simulation')
    parser.add_argument('--dt', type=float, help='Step size (ps) for MD simulation')
    
    args = parser.parse_args()
    
    # Check parameters
    if args.pdb and (args.protein or args.ligand):
        print("Error: --pdb parameter cannot be used with --protein or --ligand")
        sys.exit(1)
    
    if args.protein and not args.ligand:
        print("Error: --ligand parameter must be provided when using --protein")
        sys.exit(1)
    
    if args.ligand and not args.protein:
        print("Error: --protein parameter must be provided when using --ligand")
        sys.exit(1)
    
    if not any([args.pdb, args.protein, args.ligand]):
        parser.print_help()
        sys.exit(1)
    
    # Execute simulation
    if args.pdb:
        simulate_protein(args.pdb, args.mdp_files, args.steps, args.dt)
    elif args.protein and args.ligand:
        simulate_protein_ligand(args.protein, args.ligand, args.mdp_files, args.steps, args.dt)

if __name__ == '__main__':
    main()
