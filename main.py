import gzip
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from Bio import SeqIO
import os, sys, platform, time

import matplotlib
matplotlib.use('QT5Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

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
        self.mash_all = {} #Holds all mash output 
        self.mash_dict = {} #Holds filtered mash output data
        self.output = [] #List to hold all output data

        ### crRNA table settings
        self.crRNA_table.setColumnCount(7)  # hardcoded because there will always be nine columns
        self.crRNA_table.setShowGrid(False)
        self.crRNA_table.setHorizontalHeaderLabels("Gene;Sequence;Avg. Off-Target;On-Target;Location;PAM;Strand".split(";"))
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
        #self.build_index(self.filtered_out_path)
        self.format_gRNA(self.path1)
        self.OTF()
        self.output_parse()

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
            path_to_exe = os.getcwd() + "/bin/index_builder_linux"
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
            os.remove(input_fasta) ###Clean up intermediate fasta file

        else:
            QtWidgets.QMessageBox.question(self, "Error!","Run Mash first!",QtWidgets.QMessageBox.Ok)


    def format_gRNA(self, input_file):
        """
        This function formats the gRNA's for off-target input
        """
        ###Assumes that gRNA csv file comes in format of gene, sequence, score, location, pam, strand
        ###Needs to be in format of loc, seq, pam, score, strand w/ semicolon sep.
        with open(input_file, "r") as file:
            self.guides = os.getcwd() + "/temp.txt"
            with open(self.guides, "w") as f:
                for index in file:
                    item = index.strip().split(",")
                    if item[1] == "Sequence":
                        continue
                    list = [item[3], item[1], item[4], item[2], item[5]]
                    ###Compressed sequences must be in format: loc;seq;pam;score;strand
                    tmp = str(list[0] + ";" + list[1] + ";" + list[2] + ";" + list[3] + ";" + list[4])
                    f.write(tmp + '\n')


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
                self.mash_all[str(hold[1])] = ["--", hold[2], hold[3]] ## Holds entire mash output
                if (float(hold[3]) > 0.001 or float(hold[3]) == 0):
                    self.mash_dict[str(hold[1])] = ["--", hold[2], hold[3]] ## key = fasta id = [name, distance, p-val], holds filtered mash output
        f.close()
        os.remove(self.mash_out_path) ###Delete intermediate Mash file

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

    def OTF(self):
        """
        This function runs the Off-Target Algorithm and saves it to OT_results.txt
        """
        ### "path to exe" "compressed guides" True "CSPR path" "output directory" "CASPERinfo path" [max mismatches = (3-5)] [tolerance value = 0.05] [average output] [detailed output]
        path_to_exe = os.getcwd() + "/bin/OT_linux"
        path_to_guides = self.guides
        self.index =  os.getcwd() + "/bin/skin_spCas9.cspr" ###Index file containing unique sequences and organism info
        self.db =  os.getcwd() + "/bin/skin_spCas9_repeats.db"
        CSPR_path = self.index
        self.ot_out = os.getcwd() + "/OT_results.txt"
        index_cmd = '"' + path_to_exe + '" "' + path_to_guides + '" "' + CSPR_path + '" "' + self.db + '" "' + self.ot_out + '" "' + self.casper_info + '" ' + "3" + " " + "0.05" + ' ' + "True False"
        print(index_cmd)
        os.system(index_cmd)
        os.remove(self.guides) ###Removes intermediate gRNA file


    def output_parse(self):
        """
        This function finds average off-target scores and populates the table, and creates a list with output data for plotting
        """

        off_target = open(self.ot_out, 'r')


        for x in off_target:
            if x[0] != '0':
                if x[0] == "D":
                    continue
                split_1 = x.replace("\n","")
                split_1 = split_1.replace(" ","")
                split_1 = split_1.split(":") # split_1[0] holds sequence, split_1[1] holds average OT score
                self.crRNA_dict[str(split_1[0])][1] = split_1[1]
            else:
                split_2 = x.replace("\n","")
                split_2 = split_2.split(",") # split_2[0] holds off-target score, split_2[1] holds index in index_file
                ### output list = sequence, off-target score, off-target organism, distance, gene, location, PAM, strand, on-target score
                self.output.append([str(split_1[0]), split_2[0], split_2[1], "distance", self.crRNA_dict[split_1[0]][0], self.crRNA_dict[split_1[0]][3], self.crRNA_dict[split_1[0]][4], self.crRNA_dict[split_1[0]][5], self.crRNA_dict[split_1[0]][2]])
            
        ## Error generated if off-target output file shows all scores as 0
        if not self.output:
            QtWidgets.QMessageBox.question(self, "Error!","Check off-target output file",QtWidgets.QMessageBox.Ok)

        with open(self.index, 'rb') as f:
            index_file = gzip.GzipFile(fileobj=f)
            index_dict = {}

            for x in index_file:
                x = x.decode("utf-8")
                if x[0] == '>':
                    line = x.replace(">","")
                    line = line.replace("(","")
                    line = line.replace(")","")
                    line = line.replace("\n","")
                    line = line.replace(",","")
                    hold = line.split(" ")
                    num = hold[len(hold)-1]
                    id = hold[0]
                    name = ""
                    i = 1
                    while (hold[i] != "complete"):
                        name += hold[i] + " "
                        i += 1
                    num = num.replace("\r","")
                    index_dict[num] = [id, name]
        f.close()

        for x in self.output:
            name = index_dict[x[2]][1]
            x[3] = self.mash_all[index_dict[x[2]][0]][1] # Distance
            x[2] = name # Off-target organism name

        ### Populates table with average off-target scores
        for index, (x,y) in enumerate(self.crRNA_dict.items()):
            avg_off = QtWidgets.QTableWidgetItem(str(y[1]))
            self.crRNA_table.setItem(index, 2, avg_off)
        self.crRNA_table.resizeColumnsToContents()
        
    
       #with open('./test_files/example_plottable.csv', 'r') as f:
    
       #     count = 0
       #     for x in f:
       #         if count > 0:
       #             line = x.replace('\n', '')
       #             arr = line.split(',')
       #             self.output.append(arr)
       #         count += 1

        markers = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4']
        
        organismList = {}

        for row in  self.output:
            if row[2] != '':
                if not row[2] in organismList:
                    organismList[row[2]] = []
                organismList[row[2]].append((row[1], row[3]))
            

        fig, axs = plt.subplots()

        count = 0
        for key in organismList:
            x_vals = []
            y_vals = []
            for i in range(0, len(organismList[key])):
                x_vals.append(float(organismList[key][i][0]))
                y_vals.append(float(organismList[key][i][1]))

            scatter = axs.scatter(x_vals, y_vals, s = 30, label = key, marker=markers[count])

            count += 1

        # produce a legend with the unique colors from the scatter
        fig.set_size_inches(3, 3, forward=True)
        legend1 = axs.legend(loc="upper right", title="Off Target Organism")
        legend1 = axs.legend( prop={'size': 4})
        axs.set_ylabel('Organism Distance')
        axs.set_xlabel('Off-Target Score')
        axs.set_title('gRNA selection')
        axs.add_artist(legend1)
        plt.tight_layout()
        
        self.plotWidget = FigureCanvas(fig)
        lay = QtWidgets.QVBoxLayout(self.total_crRNA)
        lay.setContentsMargins(0,0,0,0)
        lay.addWidget(self.plotWidget)
        

        ## self.output  ---list of lists containing data for plotting
        ## Format for each index in the list: (sequence, off-target score, off-target organism, distance, gene, location, PAM, strand, on-target score


###Plot data


if __name__ == '__main__':
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    #app = Qt.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    app.setApplicationName("crRNA Viewer")
    sys.exit(app.exec_())
