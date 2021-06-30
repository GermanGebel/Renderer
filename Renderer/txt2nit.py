txt_file = 'out.txt'
nit_file = 'out.nit'

r, g, b = [], [], []

def get_color_array(file, color_array):
    line = file.readline()
    while (line != '' and line != '\n'):
        color_array.append([float(i) for i in line.split()])
        line = file.readline()
    
with open(txt_file, 'r') as file:
    line = file.readline()
    while line != '':
        if 'RED' in line:
            get_color_array(file, r)
            print('RED')
        if 'GREEN' in line:
            get_color_array(file, g)
            print('GREEN')
        if 'BLUE' in line:
            get_color_array(file, b)
            print('BLUE')
        line = file.readline()
    file.close()


pp = PostProcessor(PPDataUnits.LUMINANCE, [], r, g, b) 
pp.SaveToHDR(nit_file, overwrite = OverwriteMode.OVERWRITE)
