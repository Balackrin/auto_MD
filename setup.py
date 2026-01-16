from setuptools import setup, find_packages
import os

# 读取README文件
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='autoMD',
    version='0.1.0',
    author='autoMD Team',
    author_email='autoMD@example.com',
    description='基于GROMACS的全自动分子动力学模拟包',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/autoMD/autoMD',
    packages=find_packages(),
    package_data={
        'autoMD': ['protein_mdp/*'],
    },
    install_requires=[
        # 依赖项（如果有）
    ],
    entry_points={
        'console_scripts': [
            'run_gromacs=autoMD.main:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.7',
)
