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
        self.dbpath = ""
        self.crRNA_dict = {}
        self.switcher = [1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column
        #self.info_path = info_path

        ### crRNA table settings
        self.crRNA_table.setColumnCount(7)  # hardcoded because there will always be nine columns
        self.crRNA_table.setShowGrid(False)
        self.crRNA_table.setHorizontalHeaderLabels("Gene;Sequence;Avg. Off-Target;On-Target;Location;Strand;PAM".split(";"))
        self.crRNA_table.horizontalHeader().setSectionsClickable(True)
        self.crRNA_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.crRNA_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.crRNA_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
#        self.back_button.clicked.connect(self.goBack)
        self.crRNA_table.horizontalHeader().sectionClicked.connect(self.table_sorting)
#        self.actionSave.triggered.connect(self.save_data)
#        self.actionOpen.triggered.connect(self.open_data)
#        self.actionOff_Target_Analysis.triggered.connect(self.Off_Target_Analysis)
#        self.actionCoTargeting.triggered.connect(self.open_cotarget)







        ###Make sure window is centered upon start-up
        self.mwfg = self.frameGeometry()  ##Create frame geometry
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Create center point
        self.mwfg.moveCenter(self.cp)  ##Assign center of frame geometry to center point
        self.move(self.mwfg.topLeft())  ##Move window to center

        ###Connect push buttons to appropriate functions
        self.browse_1.clicked.connect(self.get_file1)
        self.browse_2.clicked.connect(self.get_file2)
        self.browse_3.clicked.connect(self.get_file3)
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
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .csv file.\n\nTry again!",QtWidgets.QMessageBox.Ok)

        self.fill_table(myFile[0])


    def get_file2(self):
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose crRNA Sequence File (.csv only)")
        if (myFile[0] != "" and '.csv' in myFile[0]):
            self.path_2.setText(myFile[0].split("/")[-1])
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .csv file.\n\nTry again!",QtWidgets.QMessageBox.Ok)


    def get_file3(self):
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose crRNA Sequence File (.csv only)")
        if (myFile[0] != "" and '.csv' in myFile[0]):
            self.path_3.setText(myFile[0].split("/")[-1])
        else:
            QtWidgets.QMessageBox.question(self, "Error!","You must select a .csv file.\n\nTry again!",QtWidgets.QMessageBox.Ok)

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
#                loc.setData(QtCore.Qt.EditRole, num)
#                score.setData(QtCore.Qt.EditRole, num1)
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



if __name__ == '__main__':
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    #app = Qt.QApplication(sys.argv)
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    app.setApplicationName("crRNA Viewer")
    sys.exit(app.exec_())

#
#    def getData(self):
#        mypath = os.getcwd()
#        found = False;
#        self.dbpath = mypath
#        onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
#        orgsandendos = {}
#        shortName = {}
#        for file in onlyfiles:
#            if file.find('.cspr') != -1:
#                found = True
#                newname = file[0:-4]
#                s = newname.split('_')
#                hold = gzip.open(file, 'r')
#                buf = (hold.readline())
#                buf = str(buf)
#                buf = buf.strip("'b")
#                buf = buf[:len(buf) - 4]
#                species = buf[8:]
#                endo = str(s[1][:len(s[1]) - 1])
#                if species not in shortName:
#                    shortName[species] = s[0]
#                if species in orgsandendos:
#                    orgsandendos[species].append(endo)
#                else:
#                    orgsandendos[species] = [endo]
##                    if self.orgChoice.findText(species) == -1:
#        # auto fill the kegg search bar with the first choice in orgChoice
##        self.Search_Input.setText(self.orgChoice.currentText())
##        if found == False:
##            return False
#        self.data = orgsandendos
#        self.shortHand = shortName
##        self.endoChoice.clear()
##        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
##        self.orgChoice.currentIndexChanged.connect(self.changeEndos)
#        GlobalSettings.mainWindow.mwfg.moveCenter(GlobalSettings.mainWindow.cp)  ##Center window
#        GlobalSettings.mainWindow.move(GlobalSettings.mainWindow.mwfg.topLeft())  ##Center window
#        GlobalSettings.mainWindow.show()
#
#
#    def select_ref(self):
#        filed = QtWidgets.QFileDialog()
#        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose reference virus genome")
#        if (myFile[0] != ""):
#            self.refEdit.setText(myFile[0])
#
#    def select_anno(self):
#        filed = QtWidgets.QFileDialog()
#        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose reference annotation")
#        if (myFile[0] != ""):
#            self.annoEdit.setText(myFile[0])
#

