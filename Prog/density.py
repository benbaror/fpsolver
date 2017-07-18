import os
def main():
    path = '/home/ben/mp/Code/FP_solver/Runs/'
    bin_path = '/home/ben/mp/Code/FP_solver/Prog/'
    Dirs = os.listdir(path)
    for Dir in Dirs:
        os.chdir(path+Dir)
        files = os.listdir('./')
        if ('density.out' not in files and 'FP_last.out' in files):
            f = open('./FP.inp')
            line = f.readline()
            while (line and line[1:5] != 'NPDE'):
                line = f.readline()
            f.close()
            nM = int(line[-3])
            if nM == 1:
                os.system(bin_path + 'density < ' + bin_path + 'density_sm.inp')
            if nM == 4:
                os.system(bin_path + 'density < ' + bin_path + 'density.inp')

if __name__ == '__main__':
    main()
