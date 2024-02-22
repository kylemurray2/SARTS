#!/usr/bin/env python3
import os,shutil,argparse
from SARTS import setupStack

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--alos', action='store_true', dest='alos',help='Use this flag if processing ALOS-1 stripmap')

    return parser.parse_args()


def main(inps):
    sartPath = setupStack.__file__
    sartPath = os.path.dirname(sartPath)
    if inps.alos:
        print('Using params.yaml for ALOS-1 strip map')
        paramPath = os.path.join(sartPath, 'docs', 'params_alos_template.yaml')
    else:
        print('Using params.yaml for Sentinel-1 tops processing')
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

if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = argparse.Namespace()
    inps = cmdLineParser()
    main(inps)