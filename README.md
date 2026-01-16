# autoMD

基于GROMACS的全自动分子动力学模拟Python三方包。

## 功能

- **单独蛋白模拟**：通过单个命令实现蛋白的全自动分子动力学模拟及分析
- **蛋白-配体复合物模拟**：支持蛋白-配体复合物的分子动力学模拟及分析
- **自动化流程**：包含拓扑构建、盒子添加、溶剂化、离子添加、能量最小化、NVT平衡、NPT平衡和MD模拟等完整流程
- **GPU加速**：使用GPU加速分子动力学模拟计算

## 安装

```bash
pip install -e .
```

## 依赖

- Python 3.7+
- GROMACS (已安装并配置环境变量)

## 使用方法

### 单独蛋白模拟

```bash
run_gromacs --pdb protein.pdb
```

### 蛋白-配体复合物模拟

```bash
run_gromacs --protein protein.pdb --ligand ligand.pdb
```

## 模拟流程

### 单独蛋白模拟流程

1. **拓扑构建**：使用gmx pdb2gmx构建蛋白拓扑
2. **盒子添加**：使用gmx editconf添加立方盒子
3. **溶剂化**：使用gmx solvate添加水溶剂
4. **离子添加**：使用gmx genion添加离子中和电荷
5. **能量最小化**：使用gmx mdrun进行能量最小化
6. **NVT平衡**：在恒温下进行NVT系综平衡
7. **NPT平衡**：在恒温恒压下进行NPT系综平衡
8. **MD模拟**：执行生产级分子动力学模拟
9. **轨迹处理**：处理轨迹，去除PBC影响

## 配置文件

所有的mdp配置文件存放在`autoMD/protein_mdp`目录下：

- `ions.mdp`：离子添加配置
- `minim.mdp`：能量最小化配置
- `nvt.mdp`：NVT平衡配置
- `npt.mdp`：NPT平衡配置
- `md.mdp`：MD模拟配置

## 注意事项

1. 确保GROMACS已正确安装并配置环境变量
2. 确保PDB文件格式正确
3. 建议在Linux系统上运行，获得更好的性能
4. 首次运行可能需要下载GROMACS的力场文件

## 开发

### 项目结构

```
autoMD/
├── __init__.py
├── main.py              # 主程序入口
└── protein_mdp/         # MDP配置文件目录
    ├── ions.mdp
    ├── minim.mdp
    ├── nvt.mdp
    ├── npt.mdp
    └── md.mdp
scripts/
├── run_gromacs          # Unix命令行入口
└── run_gromacs.bat      # Windows命令行入口
setup.py                 # 安装配置
README.md                # 项目说明
```

## 许可证

MIT
