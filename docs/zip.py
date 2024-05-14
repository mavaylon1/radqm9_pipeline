import zipfile
import os
import sys

def zipfolder(foldername, target_dir):            
    zipobj = zipfile.ZipFile(foldername + '.zip', 'w', zipfile.ZIP_DEFLATED)
    rootlen = len(target_dir) + 1
    for base, dirs, files in os.walk(target_dir):
        for file in files:
            fn = os.path.join(base, file)
            zipobj.write(fn, fn[rootlen:])

zipfolder('/pscratch/sd/m/mavaylon/chem_final_data/SP', '/pscratch/sd/m/mavaylon/chem_final_data/SP') #insert your variables here
sys.exit()