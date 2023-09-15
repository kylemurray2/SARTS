#!/usr/bin/env python3
import os,shutil
from SARTS import setupStack

sartPath = setupStack.__file__
sartPath = os.path.dirname(sartPath)

# If you want to ensure the output always ends with a '/'
paramPath = os.path.join(sartPath, 'docs', 'params_template.yaml')

print('Copying from ' + paramPath)

shutil.copy(paramPath, './params.yaml')

print('Next, edit the params.yaml file. Then you can download data with downloadData.py')
