import os
import shutil
import subprocess
import linecache


def mysplit(s):
    """
    It reads in a string+number and separate them.
    
    Args: 
        s:  is the input string+number.
    
    return head:string, tail:number
    """
    head = s.rstrip('-0123456789.0')
    tail = s[len(head):]
    #print(head, tail)
    return head, tail

def readable_fatbands(path_eigfat2plot, eigfat_file, dat_file):
    """ 
    It makes the output fatbands readable and write them into a .dat file.
    
    Args:
        path_eigfat2plot:   is the path to the eigfat2plot binary.
        eigfat_file:        is the outfile produced by eig2plot binary.
        dat_file:           is the readable file in dat format.
        
    return:
        None    
    """
    myoutput = open(dat_file, 'w')
    output, errors = subprocess.Popen([path_eigfat2plot, eigfat_file], stdout=myoutput, stderr=subprocess.PIPE, universal_newlines=True).communicate()
    myoutput.close()
    
    return None

def from_datfile(path_datfile, number):
    """
    It reads the datfile containing fatbands and produce the amount of energy and fatband for specific a band.
    
    Args:
        path_datfile: the path to the datfile containing the amount of energy and fatbands.
        number:       it expecifies the number of the line from which the fatband for specific band starts.
    
    return:
        energy, fatband
    """
    
    energy = linecache.getline(path_datfile, number).rstrip().split()[1]
    linecache.clearcache()
    fatband = linecache.getline(path_datfile, number).rstrip().split()[2]
    linecache.clearcache()
    return energy, fatband

PATH = './'

# ***************************************************************************************************
# Start producing the datfiles.
for fat_dir in next(os.walk(PATH))[1]:
    print(fat_dir)
    os.chdir(fat_dir+'/fatbands/')
    shutil.copy('../../../strain0.0/fatbands/MoS2.mpr', './')
    shutil.copy('MoS2.bands.WFSX','MoS2.WFSX')
    subprocess.run(['/home/mb1988/opt/siesta-4.1-b3/Util/COOP/fat', 'MoS2'])
    #print("The fatbands are calculating ...")
    #wait = input("PRESS ENTER TO CONTINUE.")
    #print("Make them readable and store them in a .dat foramt.")
    readable_fatbands('/home/mb1988/opt/siesta-4.1-b3/Util/Bands/eigfat2plot', 'MoS2.fatbands_Mo_4d.EIGFAT', 'Mo.4d.dat')
    readable_fatbands('/home/mb1988/opt/siesta-4.1-b3/Util/Bands/eigfat2plot', 'MoS2.fatbands_Mo_4dx2-y2.EIGFAT', 'Mo.4dx2-y2.dat')
    readable_fatbands('/home/mb1988/opt/siesta-4.1-b3/Util/Bands/eigfat2plot', 'MoS2.fatbands_Mo_4dxy.EIGFAT','Mo.4dxy.dat')
    readable_fatbands('/home/mb1988/opt/siesta-4.1-b3/Util/Bands/eigfat2plot', 'MoS2.fatbands_S_3p.EIGFAT', 'S.3p.dat')
    readable_fatbands('/home/mb1988/opt/siesta-4.1-b4/Util/Bands/gnubands', 'MoS2.bands', 'bs.dat')
    os.chdir('../../')

# ***************************************************************************************************    
# Start reading in the desired properties and writing them into the appropriate output for the case of 4dx2-y2.
fhand = open('v_2stop_fatbands_4dx2-y2.dat', 'w')

# First read the energy and its corresponding fatbands for the case of strain0.0
print(linecache.getline('../strain0.0/fatbands/Mo.4dx2-y2.dat', 96))
energy_lower, fatbands_lower = from_datfile('../strain0.0/fatbands/Mo.4dx2-y2.dat', 96)
print(linecache.getline('../strain0.0/fatbands/Mo.4dx2-y2.dat', 188))
energy_upper, fatbands_upper = from_datfile('../strain0.0/fatbands/Mo.4dx2-y2.dat', 188)
fhand.write('   ' + str(0.0) + '        ' + energy_lower + '        ' + fatbands_lower + '\n')
fhand.write('   ' + str(0.0) + '        ' + energy_upper + '        ' + fatbands_upper + '\n')

for fat_dir in next(os.walk(PATH))[1]:
    print(fat_dir)
    #print(mysplit(fat_dir)[0], mysplit(fat_dir)[1])
    strain = float(mysplit(fat_dir)[1])
    print(strain)
    #os.chdir(fat_dir)
    print(linecache.getline(fat_dir+'/fatbands/Mo.4dx2-y2.dat', 96))
    energy_lower, fatbands_lower = from_datfile(fat_dir+'/fatbands/Mo.4dx2-y2.dat', 96)
    print(linecache.getline(fat_dir+'/fatbands/Mo.4dx2-y2.dat', 188))
    energy_upper, fatbands_upper = from_datfile(fat_dir+'/fatbands/Mo.4dx2-y2.dat', 188)
    fhand.write('   ' + str(strain) + '        ' + energy_lower + '        ' + fatbands_lower + '\n')
    fhand.write('   ' + str(strain) + '        ' + energy_upper + '        ' + fatbands_upper + '\n')
    #os.chdir('../')
  
fhand.close()

# ***************************************************************************************************
# Start reading in the desired properties and writing them into the appropriate output for the case of 4dxy.
fhand = open('v_2stop_fatbands_4dxy.dat', 'w')

# First read the energy and its corresponding fatbands for the case of strain0.0
print(linecache.getline('../strain0.0/fatbands/Mo.4dxy.dat', 96))
energy_lower, fatbands_lower = from_datfile('../strain0.0/fatbands/Mo.4dxy.dat', 96)
print(linecache.getline('../strain0.0/fatbands/Mo.4dxy.dat', 188))
energy_upper, fatbands_upper = from_datfile('../strain0.0/fatbands/Mo.4dxy.dat', 188)
fhand.write('   ' + str(0.0) + '        ' + energy_lower + '        ' + fatbands_lower + '\n')
fhand.write('   ' + str(0.0) + '        ' + energy_upper + '        ' + fatbands_upper + '\n')

for fat_dir in next(os.walk(PATH))[1]:
    print(fat_dir)
    #print(mysplit(fat_dir)[0], mysplit(fat_dir)[1])
    strain = float(mysplit(fat_dir)[1])
    print(strain)
    #os.chdir(fat_dir)
    print(linecache.getline(fat_dir+'/fatbands/Mo.4dxy.dat', 96))
    energy_lower, fatbands_lower = from_datfile(fat_dir+'/fatbands/Mo.4dxy.dat', 96)
    print(linecache.getline(fat_dir+'/fatbands/Mo.4dxy.dat', 188))
    energy_upper, fatbands_upper = from_datfile(fat_dir+'/fatbands/Mo.4dxy.dat', 188)
    fhand.write('   ' + str(strain) + '        ' + energy_lower + '        ' + fatbands_lower + '\n')
    fhand.write('   ' + str(strain) + '        ' + energy_upper + '        ' + fatbands_upper + '\n')
    #os.hdir('../')
    
fhand.close()
