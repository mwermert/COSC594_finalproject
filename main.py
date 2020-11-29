import gzip
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from Bio import SeqIO
import os, sys, platform, time

### GUI Declarations/Setup
class MyMainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MyMainWindow, self).__init__()
        uic.loadUi(os.getcwd() + '/MyMainWindow.ui', self)

        ###Initialize any necessary variables here
        self.crRNA_dict = {} #Stores crRNA info, key = sequence
        self.switcher = [1,1,1,1,1,1,1]  #Keeps track of where we are in the sorting clicking for each column
        self.casper_info = os.getcwd() + "/bin/CASPERinfo" #Supplementary file for Off-Target and Index Builder algorithms
        self.path1 = "" #Path to crRNA file
        self.path2 = "" #Path to ref file
        self.path3 = "" #Path to query file
        self.query_dict = {}
        self.mash_dict = {}

        ### crRNA table settings
        self.crRNA_table.setColumnCount(7)  # hardcoded because there will always be nine columns
        self.crRNA_table.setShowGrid(False)
        self.crRNA_table.setHorizontalHeaderLabels("Gene;Sequence;Avg. Off-Target;On-Target;Location;Strand;PAM".split(";"))
        self.crRNA_table.horizontalHeader().setSectionsClickable(True)
        self.crRNA_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.crRNA_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.crRNA_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.crRNA_table.horizontalHeader().sectionClicked.connect(self.table_sorting)

        ###Make sure window is centered upon start-up
        self.mwfg = self.frameGeometry()  ##Create frame geometry
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Create center point
        self.mwfg.moveCenter(self.cp)  ##Assign center of frame geometry to center point
        self.move(self.mwfg.topLeft())  ##Move window to center

        ###Connect push buttons to appropriate functions
        self.browse_1.clicked.connect(self.get_file1)
        self.browse_2.clicked.connect(self.get_file2)
        self.browse_3.clicked.connect(self.get_file3)
        self.run_OT.clicked.connect(self.run_all)
        self.show()

    """ Start defining functions here """
    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.crRNA_table.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.crRNA_table.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)

    def get_file1(self):
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose crRNA Sequence File (.csv only)")
        if (myFile[0] != "" and '.csv' in myFile[0]):
            self.path_1.setText(myFile[0].split("/")[-1])
            self.path1 = myFile[0]
            self.fill_table(myFile[0])
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .csv file.\n\nTry again!",QtWidgets.QMessageBox.Ok)




    def get_file2(self):
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose Reference File (.fasta only)")
        if (myFile[0] != "" and '.fasta' in myFile[0]):
            self.path_2.setText(myFile[0].split("/")[-1])
            self.path2 = myFile[0]
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .fasta file.\n\nTry again!",QtWidgets.QMessageBox.Ok)


    def get_file3(self):
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose Query File (.fasta only)")
        if (myFile[0] != "" and '.fasta' in myFile[0]):
            self.path_3.setText(myFile[0].split("/")[-1])
            self.path3 = myFile[0]
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .fasta file.\n\nTry again!",QtWidgets.QMessageBox.Ok)

    def fill_table(self, input_path):
        with open(input_path, "r") as f: #Table won't populate with file open
            readlines = f.readlines()[1:] #Assumes first line contains column names

        ###Assumes formata is gene, sequence, on-target score, location, PAM, strand
        self.crRNA_table.setRowCount(len(readlines))
        self.crRNA_table.setColumnCount(7)
        for index, item in enumerate(readlines):
            item = item.strip().split(",")
            self.crRNA_dict[str(item[1])] = [item[0], "--", item[2], item[3], item[4], item[5]] ###Add sequence to dictionary
#                gene = item[0]
#                seq = item[1]
#                on = item[2]
#                loc = item[3]
#                pam = item[4]
#                strand = item[5]
            ###Initilize table items
            gene = QtWidgets.QTableWidgetItem(str(item[0]))
            seq = QtWidgets.QTableWidgetItem(str(item[1]))
            avg_off = QtWidgets.QTableWidgetItem("--") ###Set to -- initially before off-target analysis done
            score = QtWidgets.QTableWidgetItem(str(item[2]))
            loc = QtWidgets.QTableWidgetItem(str(item[3]))
            pam = QtWidgets.QTableWidgetItem(str(item[4]))
            strand = QtWidgets.QTableWidgetItem(str(item[5]))

            ###Set table items
            self.crRNA_table.setItem(index, 0, gene)
            self.crRNA_table.setItem(index, 1, seq)
            self.crRNA_table.setItem(index, 2, avg_off)
            self.crRNA_table.setItem(index, 3, score)
            self.crRNA_table.setItem(index, 4, loc)
            self.crRNA_table.setItem(index, 5, pam)
            self.crRNA_table.setItem(index, 6, strand)

        self.crRNA_table.resizeColumnsToContents()

    def run_all(self):
        """
        This is the master command that runs when GUI button is clicked for running the entire analysis.
        It calls the functions that run Mash, Index Builder, Off-Target Algorithm, and requisite parsing
        functions.
        """
        self.run_mash()
        self.filter_query()
        self.build_index(self.filtered_out_path)
        self.Compress_gRNA(self.path1)
        self.OTF()

    def run_mash(self):
        """
        This function rans Mash and saves the output to self.mash_out_path
        """

        if (self.path2 != "" and self.path3 != ""):
            ###Set up command###
            # <exe path> dist <options> <ref seq> <query seq>
            # Sketch size by default 10000, kmer = 16
            ####################

            ###Initialize variables
            ref_file = self.path2
            query_file = self.path3
            self.mash_out_path = os.getcwd() + "/mash_out.txt"
            path_to_exe = "/Users/ddooley/bioinformatics_packages/Mash/mash"
            cmd = ""
            flag_str = "-s 10000 -k 16 -i"
            cmd += path_to_exe + " dist " + " " + flag_str + " " + ref_file + " " + query_file + " > " + self.mash_out_path
            os.system(cmd) ###Run the command
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select reference and query .fasta files.",QtWidgets.QMessageBox.Ok)

    def build_index(self, input_fasta):
    ###"path to exe" "endo" "PAM" "org_code" "FALSE" "output dir" "CASPERinfo" "query file" "org_name "guide length" "seed length" " "
        if self.mash_out_path: ###Makes sure that Mash has been run first
            path_to_exe = os.getcwd() + "/bin/index_builder_mac"
            org_code = "temp"
            endo = "spCas9"
            pam = "NGG"
            output_dir = os.getcwd()
            self.index = output_dir + "/temp_spCas9.cspr"
            org_name = "temp_index"
            guide_len = "20"
            seed_len = "16"
            index_cmd = '"' + path_to_exe + '" "' + endo + '" "' + pam + '" "' + org_code + '" "' + "FALSE" + '" "' + output_dir + '" "' + self.casper_info + '" "' + input_fasta + '" "' + org_name + '" "' + guide_len + '" "'  + seed_len + '" " "'
            print(index_cmd)
            os.system(index_cmd)

        else:
            QtWidgets.QMessageBox.question(self, "Error!","Run Mash first!",QtWidgets.QMessageBox.Ok)
    
 
    def Compress_gRNA(self, input_file):
        """
        This function compresses the gRNA's for off-target input
        """

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
        
        with open(input_file, "r") as file:
            self.guides = os.getcwd() + "/temp.txt"
            with open(self.guides, "w") as f:
                for index in file:
                    item = index.strip().split(",")
                    if item[1] == "Sequence":
                        continue
                    list = [item[1], item[2], item[3], item[4], item[5]]
                    tmp = str(compress(str(list[0]),64) + "," + compress(str(list[1]),64) + str(list[2]) + compress(str(list[3]),64) + "," + compress(str(list[4]),64))
                    f.write(tmp + '\n')
            f.close()
            file.close()
        
    
    def filter_query(self):
        """
        This function filters out query sequences with p-values under 0.001
        """
        with open(self.mash_out_path) as f:
            for line in f:
                line = " ".join(line.split())
                hold = []
                for x in line.split(' '):
                    hold.append(x)
                if float(hold[3]) > 0.001:
                    self.mash_dict[str(hold[1])] = ["--", hold[2], hold[3]] ## key = fasta id = [name, distance, p-val]
        f.close()

        q = SeqIO.parse(open(self.path3), 'fasta')
        self.filtered_out_path = os.getcwd() + "/filtered_query.fasta"  ## Holds sequences for off-target algorithm
        with open(self.filtered_out_path, 'w') as outfile:
            for fasta in q:
                id, description = fasta.id, fasta.description
                if id in self.mash_dict.keys():
                    line = description.replace(",","")
                    hold = []
                    name = ""
                    for x in line.split(' '):
                        hold.append(x)
                    for i in range(1,len(hold)-2):
                        name += hold[i] + " "
                    self.mash_dict[id][0] = name
                    SeqIO.write(fasta, outfile, 'fasta')
        outfile.close()
        
            
    def OTF(self):
        """
        This function rans the Off-Target Algorithm and saves it to OT_results.txt
        """
        ### "path to exe" "compressed guides" True "CSPR path" "output path" "CASPERinfo path" [max mismatches = (3-5)] [tolerance value = 0.05] [average output] [detailed output]
        path_to_exe = os.getcwd() + "/executables/OT_mac"
        path_to_guides = self.guides
        CSPR_path = self.index
        self.ot_out = os.getcwd() + "/OT_results.txt"
        index_cmd = '"' + path_to_exe + '" "' + path_to_guides + '" ' + "True" + ' "' + CSPR_path + '" "' + self.ot_out + '" "' + self.casper_info + '" ' + "3" + ' ' + "0.05" + ' ' + "True" + ' ' + "False"
        print(index_cmd)
        os.system(index_cmd)

        
###Off-target + mash output parsing function


if __name__ == '__main__':
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    #app = Qt.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    app.setApplicationName("crRNA Viewer")
    sys.exit(app.exec_())
