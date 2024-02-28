#!/usr/bin/env python3
# David Bekaert

import argparse, zipfile, sys, shutil, tarfile,zipfile,os,glob
from stripmapStack import uncompressFile


def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser(description='Prepare ALOS raw processing (unzip/untar files, '
                                     'organize in date folders, generate script to unpack into isce formats).')
    parser.add_argument('-i', '--input', dest='inputDir', type=str, required=True,
            help='directory with the raw data')
    parser.add_argument('-rmfile', '--rmfile', dest='rmfile',action='store_true', default=False,
            help='Optional: remove zip/tar/compressed files after unpacking into date structure '
                 '(default is to keep in archive fo  lder)')
    parser.add_argument('-o', '--output', dest='outputDir', type=str, required=False,
            help='output directory where data needs to be unpacked into isce format (for script generation).')
    parser.add_argument('-t', '--text_cmd', dest='text_cmd', type=str, default='',
            help='text command to be added to the beginning of each line of the run files. ')
    parser.add_argument('--dual2single','--fbd2fbs', dest='fbd2fbs', action='store_true',
            help='resample the FBD acquisitions to FBS. Recommended for "interferogram" workflow without ionosphere.')
    return parser


def cmdLineParse(iargs=None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args = iargs)

    # parsing required inputs
    inps.inputDir = os.path.abspath(inps.inputDir)
    # parsing optional inputs
    if inps.outputDir:
        inps.outputDir = os.path.abspath(inps.outputDir)
    return inps


def get_Date(ALOSfolder):

    # will search for different version of workreport to be compatible with ASf, WInSAR etc
    workreport_files = ('*workreport','summary.txt')
    for workreport_file in workreport_files:
        workreports = glob.glob(os.path.join(ALOSfolder,workreport_file))

        # if nothing is found return a failure
        if len(workreports) > 0:
            for workreport in workreports:
                template_dict = {}
                with open(workreport) as openfile:
                    for line in openfile:
                        c = line.split("=")
                        template_dict[c[0].strip()] = c[1].strip()
                acquisitionDate = (str(template_dict['Img_SceneCenterDateTime'][1:9]))
                if acquisitionDate:
                    successflag = True
                    return successflag, acquisitionDate

    # if it reached here it could not find the acqusiitionDate
    successflag = False                                                                                                                                                    
    acquisitionDate = 'FAIL'
    return successflag, acquisitionDate


def get_ALOS_ALP_name(infile):
    """Get the ALPSRP075780620 name from compress file in various format."""
    outname = None
    fbase = os.path.basename(infile)
    if fbase.startswith("ALP"):
        outname = fbase.split("-")[0]
    else:
        fext = os.path.splitext(infile)[1]
        if fext in ['.tar', '.gz']:
            with tarfile.open(infile, 'r') as tar:
                file_list = tar.getnames()
        elif fext in ['.zip']:
            with zipfile.ZipFile(infile, 'r') as z:
                file_list = z.namelist()
        else:
            raise ValueError('unrecognized file extension: {}'.format(fext))
        led_file = [i for i in file_list if 'LED' in i][0]
        led_file = os.path.basename(led_file)
        outname = [i for i in led_file.split("-") if 'ALP' in i][0]
    return outname




def check_and_delete_zip(filename):
    flag = False
    try:
        # Attempt to open the zip file
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            print(f"{filename} is a valid zip file and can be opened.")
    except zipfile.BadZipFile:
        flag =True
        # If the file is not a valid zip file, delete it
        os.remove(filename)
        print(f"{filename} is not a valid zip file and has been deleted.")
    return flag

def main(iargs=None):
    '''
    The main driver.
    '''
    # inps = argparse.Namespace()
    # inps.inputDir = 'SLCS'
    # inps.outputDir = 'SLCS'
    # inps.rmfile = False
    # inps.fbd2fbs = True
    # inps.text_cmd = ''
    
    inps = cmdLineParse(iargs)

    zipFiles = glob.glob('SLCS/*zip')
    nbad = 0
    for z in zipFiles:
        flag = check_and_delete_zip(z)
        if flag:
            nbad+=1
    
    if nbad>0:
        print('Rerun downloadData.py')
        sys.exit(1)

    # filename of the runfile
    run_unPack = 'run_unPackALOS'   


    ALOS_extension = os.path.join(inps.inputDir, '*.zip')
    # loop over zip/tar files
    ALOS_filesfolders = sorted(glob.glob(ALOS_extension))
    for ALOS_infilefolder in ALOS_filesfolders:
        ## the path to the folder/zip
        workdir = os.path.dirname(ALOS_infilefolder)

        ## get the output name folder without any extensions
        ALOS_outfolder = ALOS_infilefolder.split('/')[-1].split('.')[0]
        # add the path back in
        ALOS_outfolder = os.path.join(inps.outputDir, ALOS_outfolder)

        # loop over two cases (either file or folder): 
        ### this is a file, try to unzip/untar it
        if os.path.isfile(ALOS_infilefolder):
            # unzip the file in the outfolder
            successflag_unzip = uncompressFile.uncompressfile(ALOS_infilefolder, ALOS_outfolder)

            # put failed files in a seperate directory
            if not successflag_unzip:
                os.makedirs(os.path.join(workdir,'FAILED_FILES'), exist_ok=True)
                os.rename(ALOS_infilefolder,os.path.join(workdir,'FAILED_FILES','.'))
            else:
                # check if file needs to be removed or put in archive folder
                if inps.rmfile:
                    os.remove(ALOS_infilefolder)
                    print('Deleting: ' + ALOS_infilefolder)
                else:
                    os.makedirs(os.path.join(workdir,'ARCHIVED_FILES'), exist_ok=True)
                    cmd  = 'mv ' + ALOS_infilefolder + ' ' + os.path.join(workdir,'ARCHIVED_FILES','.')
                    os.system(cmd)



    # loop over the different ALOS folders and organize in date folders
    ALOS_folders = glob.glob(os.path.join(inps.outputDir, 'ALP*'))                        
    for ALOS_folder in ALOS_folders:
        # get the date
        successflag, imgDate = get_Date(ALOS_folder)       
        
        workdir = os.path.dirname(ALOS_folder)
        if successflag:
            # move the file into the date folder
            SLC_dir = os.path.join(workdir,imgDate,'')
            os.makedirs(SLC_dir, exist_ok=True)

            # check if the folder already exist in that case overwrite it
            ALOS_folder_out = os.path.join(SLC_dir,os.path.basename(ALOS_folder))
            if os.path.isdir(ALOS_folder_out):
                shutil.rmtree(ALOS_folder_out)
            # move the ALOS acqusition folder in the date folder
            cmd  = 'mv ' + ALOS_folder + ' ' + SLC_dir + '.' 
            os.system(cmd)

            print ('Succes: ' + imgDate)
        else:
            print('Failed: ' + ALOS_folder)
        

    # now generate the unpacking script for all the date dirs
    dateDirs = sorted(glob.glob(os.path.join(inps.outputDir,'2*')))
    if inps.outputDir is not None:
        for dateDir in dateDirs:
            AlosFiles = glob.glob(os.path.join(dateDir, 'ALP*'))
            if len(AlosFiles)>0:
                acquisitionDate = os.path.basename(dateDir)
                slcDir = os.path.join(inps.outputDir, acquisitionDate)
                os.makedirs(slcDir, exist_ok=True)
                cmd = 'unpackFrame_ALOS_raw.py -i ' + os.path.abspath(dateDir) + ' -o ' + slcDir      
                IMG_files = glob.glob(os.path.join(AlosFiles[0],'IMG*'))
                if inps.fbd2fbs:
                    #recommended for regular interferometry to use all FBS bandwidth
                    if len(IMG_files) == 2:
                        cmd += ' -f fbd2fbs '
                else:
                    #used for ionosphere workflow for simplicity
                    if len(IMG_files) == 1:
                        cmd = cmd + ' -f  fbs2fbd ' 
                if len(AlosFiles) > 1:
                    cmd = cmd + ' -m' 
                print (cmd)
                # Execute the command
                # unpackFrame_ALOS_raw.py needs to be in system path
                os.system(inps.text_cmd + cmd)
    return


if __name__ == '__main__':

    main()
