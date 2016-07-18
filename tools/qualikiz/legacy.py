import os
def convert_current_to_2_3_0(inputdir):
    print ("go!")
    # Magic to fix ordering of the p-files
    os.rename(os.path.join(inputdir, 'p14.bin'), 'temp')
    for i in range(15, 21):
        os.rename(os.path.join(inputdir, 'p' + str(i) + '.bin'),
                  os.path.join(inputdir, 'p' + str(i - 1) + '.bin'))
    os.rename('temp', os.path.join(inputdir, 'p20.bin'))


    os.rename(os.path.join(inputdir, 'p43.bin'), 'temp')
    os.rename(os.path.join(inputdir, 'p44.bin'), 'temp2')
    for i in reversed(range(36, 43)):
        os.rename(os.path.join(inputdir, 'p' + str(i) + '.bin'),
                  os.path.join(inputdir, 'p' + str(i + 2) + '.bin'))
    os.rename('temp', os.path.join(inputdir, 'p36.bin'))
    os.rename('temp2', os.path.join(inputdir, 'p37.bin'))
