import sys
import os

# Add the current directory to the path
sys.path.append('.')

from autoMD.main import modify_mdp_file

# Test modifying the mdp file
modify_mdp_file('test.mdp', steps=50000, dt=0.001)

# Print the modified file
print("\nModified mdp file content:")
with open('test.mdp', 'r') as f:
    print(f.read())