#!/usr/bin/env python3
import os,shutil
from SARTS import setupStack

sartPath = setupStack.__file__
sartPath = os.path.dirname(sartPath)
paramPath = os.path.join(sartPath, 'docs', 'params_template.yaml')

if not os.path.isfile('./params.yaml'):
    print('Copying from ' + paramPath)
    shutil.copy(paramPath, './params.yaml')
else:
    print('params.yaml already exists.')

# Define the path to your file
file_path = os.path.join(sartPath, 'docs', 'sarts.txt')

# Open the file and print its contents
with open(file_path, 'r') as file:
    contents = file.read()
    print(contents)

print('\n')
print('Next, edit the params.yaml file. Then you can download data with downloadData.py')
