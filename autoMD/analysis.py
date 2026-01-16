import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5

class GMXAnalysis:
    def __init__(self, work_dir="gmx_output", output_dir="gmx_analysis"):
        self.work_dir = work_dir
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def run_cmd(self, cmd):
        """Run command"""
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=self.work_dir)
        return result

    def analyze_rmsd(self):
        """RMSD analysis"""
        print("Calculating RMSD...")

        # Calculate backbone RMSD
        self.run_cmd("echo '4\n4' | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns")

        # Read data
        data = np.loadtxt(f"{self.work_dir}/rmsd.xvg", comments=['@', '#'])
        time = data[:, 0]
        rmsd = data[:, 1] * 10  # nm to Å

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        ax.plot(time, rmsd, linewidth=1.5, color='#2E86AB')
        ax.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
        ax.set_ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
        ax.set_title('Backbone RMSD', fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/rmsd.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/rmsd.pdf', bbox_inches='tight')
        plt.close()

        print("✓ RMSD completed")

    def analyze_rmsf(self):
        """RMSF analysis"""
        print("Calculating RMSF...")

        # Calculate C-alpha RMSF
        self.run_cmd("echo '3' | gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg -res")

        # Read data
        data = np.loadtxt(f"{self.work_dir}/rmsf.xvg", comments=['@', '#'])
        residues = data[:, 0]
        rmsf = data[:, 1] * 10  # nm to Å

        # Plot
        fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
        ax.plot(residues, rmsf, linewidth=2, color='#A23B72')
        ax.fill_between(residues, rmsf, alpha=0.3, color='#A23B72')
        ax.set_xlabel('Residue Number', fontsize=14, fontweight='bold')
        ax.set_ylabel('RMSF (Å)', fontsize=14, fontweight='bold')
        ax.set_title('C-alpha RMSF', fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/rmsf.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/rmsf.pdf', bbox_inches='tight')
        plt.close()

        print("✓ RMSF completed")

    def analyze_rg(self):
        """Radius of gyration"""
        print("Calculating Rg...")

        # Calculate Rg
        self.run_cmd("echo '1' | gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg")

        # Read data
        data = np.loadtxt(f"{self.work_dir}/gyrate.xvg", comments=['@', '#'])
        time = data[:, 0]
        rg = data[:, 1] * 10  # nm to Å

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        ax.plot(time, rg, linewidth=1.5, color='#F18F01')
        ax.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Radius of Gyration (Å)', fontsize=14, fontweight='bold')
        ax.set_title('Radius of Gyration', fontsize=16, fontweight='bold')
        ax.axhline(y=np.mean(rg), color='red', linestyle='--', linewidth=1.5,
                   label=f'Mean: {np.mean(rg):.2f} Å', alpha=0.7)
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/rg.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/rg.pdf', bbox_inches='tight')
        plt.close()

        print("✓ Rg completed")

    def analyze_sasa(self):
        """Solvent Accessible Surface Area"""
        print("Calculating SASA...")

        # Calculate SASA
        self.run_cmd("echo '1' | gmx sasa -s md.tpr -f md_noPBC.xtc -o sasa.xvg -tu ns")

        # Read data
        data = np.loadtxt(f"{self.work_dir}/sasa.xvg", comments=['@', '#'])
        time = data[:, 0]
        sasa = data[:, 1] * 100  # nm^2 to Å^2

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        ax.plot(time, sasa, linewidth=1.5, color='#06A77D')
        ax.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
        ax.set_ylabel('SASA (Å²)', fontsize=14, fontweight='bold')
        ax.set_title('Solvent Accessible Surface Area', fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/sasa.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/sasa.pdf', bbox_inches='tight')
        plt.close()

        print("✓ SASA completed")
    
    def analyze_hbonds(self):
        """Hydrogen bond analysis"""
        print("Analyzing hydrogen bonds...")

        # Protein internal hydrogen bonds
        self.run_cmd("echo '1\n1' | gmx hbond -s md.tpr -f md_noPBC.xtc -num hbond.xvg -tu ns")

        # Read data
        data = np.loadtxt(f"{self.work_dir}/hbond.xvg", comments=['@', '#'])
        time = data[:, 0]
        hbonds = data[:, 1]

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
        ax.plot(time, hbonds, linewidth=1.5, color='#D62839')
        ax.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Number of H-bonds', fontsize=14, fontweight='bold')
        ax.set_title('Hydrogen Bonds', fontsize=16, fontweight='bold')
        ax.axhline(y=np.mean(hbonds), color='blue', linestyle='--', linewidth=1.5,
                   label=f'Mean: {np.mean(hbonds):.1f}', alpha=0.7)
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/hbonds.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{self.output_dir}/hbonds.pdf', bbox_inches='tight')
        plt.close()

        print("✓ Hydrogen bond analysis completed")