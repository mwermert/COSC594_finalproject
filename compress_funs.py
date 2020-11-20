import os

def nt2int(nt):
    if nt == 'A':
        return 0
    elif nt == 'T':
        return 1
    elif nt == 'C':
        return 2
    elif nt == 'G':
        return 3
    else:
        return 0



def compress(uncompressed, base):
    base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"
    compseq = 0
    if type(uncompressed) == str:
        for i in range(len(uncompressed)):
            val = nt2int(uncompressed[i]) * pow(4, i)  # multiplying by power-4 converts to base10
            compseq += val
        uncompressed = compseq
    compreturn = str()
    while uncompressed >= base:
        rem = uncompressed%base
        uncompressed = int(uncompressed/base)
        compreturn = base_array_64[rem] + compreturn
    compreturn = base_array_64[uncompressed] + compreturn
    return compreturn



def compress_guides(seq_line):
    ### Assumes that each seq line is a list of [seq, score, location, pam, strand] (agnostic of gene)
    OT_tmp = os.getcwd() + "/temp.txt" ###open file for writing compressed guides to
    with open(OT_tmp, "w+") as f:
        for item in seq_line:
            print(item)
            tmp = str(compress(str(item[0]),64) + "," + compress(str(item[1]),64) + str(item[2]) + compress(str(item[3]),64) + "," + compress(str(item[4]),64))
            f.write(tmp + "\n")
