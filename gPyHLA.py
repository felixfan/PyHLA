#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore
import pandas as pd
import sys
import os

 ##############################3
class PyHLAWin(QtGui.QWidget):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.width = 700
        self.height =600
        self.setFixedSize(self.width, self.height)
        self.move(200, 200)
        self.setWindowTitle('PyHLA v1.1.1')
        #######################################################
        self.gfile = ''
        self.covfile = ''
        self.outfile = ''
        #######################################################
        self.radio1 = QtGui.QRadioButton("&Summary")
        self.radio2 = QtGui.QRadioButton("&Association")
        self.radio3 = QtGui.QRadioButton("&Alignment")
        self.radio4 = QtGui.QRadioButton("&Zygosity")
        self.radio5 = QtGui.QRadioButton("&Interaction")

        self.gLab = QtGui.QLabel('Genotype file')
        self.gEdit = QtGui.QLineEdit()
        self.gButton = QtGui.QPushButton("Open")
        self.gvButton = QtGui.QPushButton("View")

        self.vLab = QtGui.QLabel('Covariates file')
        self.vEdit = QtGui.QLineEdit()
        self.vButton = QtGui.QPushButton("Open")
        self.vvButton = QtGui.QPushButton("View")

        self.oLab = QtGui.QLabel('Output file')
        self.oEdit = QtGui.QLineEdit()
        self.oButton = QtGui.QPushButton("Open")

        self.traitLab = QtGui.QLabel('Trait')
        self.traitCombo = QtGui.QComboBox()
        self.traitCombo.addItem("disease trait/case-control study")
        self.traitCombo.addItem("quantitative trait")

        self.levelLab = QtGui.QLabel('Level')
        self.levelCombo = QtGui.QComboBox()
        self.levelCombo.addItem("allele")
        self.levelCombo.addItem("residue")

        self.testLab = QtGui.QLabel('Test')
        self.testCombo = QtGui.QComboBox()
        self.testCombo.addItem("Fisher's exact test")
        self.testCombo.addItem("Pearson chi-squared test")
        self.testCombo.addItem("Logistic regression")

        self.digitLab = QtGui.QLabel('Digit')
        self.digitCombo = QtGui.QComboBox()
        self.digitCombo.addItems(["4","2","6"])

        self.modelLab = QtGui.QLabel('Model')
        self.modelCombo = QtGui.QComboBox()
        self.modelCombo.addItems(["allelic", "dom", "rec"])
        # self.modelCombo.addItems(["allelic", "dominant", "recessive"])

        self.adjLab = QtGui.QLabel('adjustment')
        self.adjCombo = QtGui.QComboBox()
        self.adjCombo.addItems(["FDR", "FDR_BY","Bonferroni", "Holm"])

        self.frqLab = QtGui.QLabel('allele frequency')
        self.frqEdit = QtGui.QLineEdit()

        self.permLab = QtGui.QLabel('Number of permutation')
        self.permEdit = QtGui.QLineEdit()

        self.covLab = QtGui.QLabel('Covariates name')
        self.covEdit = QtGui.QLineEdit()

        self.consensusLabel = QtGui.QLabel('consensus')
        self.consensusCB = QtGui.QCheckBox('--consensus')

        self.runButton = QtGui.QPushButton("Run")
        #######################################################
        ### layout
        self.grid = QtGui.QGridLayout()
        self.grid.setSpacing(10)

        self.grid.addWidget(self.radio1,0,0) # first row
        self.grid.addWidget(self.radio2,0,1)
        self.grid.addWidget(self.radio3,0,2)
        self.grid.addWidget(self.radio4,0,3)
        self.grid.addWidget(self.radio5,0,4)

        self.grid.addWidget(self.gLab, 1, 0)  # second row
        self.grid.addWidget(self.gEdit, 1, 1, 1, 2)
        self.grid.addWidget(self.gButton, 1, 3)
        self.grid.addWidget(self.gvButton, 1, 4)

        self.grid.addWidget(self.vLab, 2, 0)  # third row
        self.grid.addWidget(self.vEdit, 2, 1,1,2)
        self.grid.addWidget(self.vButton, 2, 3)
        self.grid.addWidget(self.vvButton, 2, 4)

        self.grid.addWidget(self.oLab, 3, 0) # fourth row
        self.grid.addWidget(self.oEdit, 3, 1, 1, 3)
        self.grid.addWidget(self.oButton, 3, 4)

        self.grid.addWidget(self.traitLab, 4, 0) # 5
        self.grid.addWidget(self.traitCombo, 4, 1, 1, 4)

        self.grid.addWidget(self.levelLab, 5, 0) # 6
        self.grid.addWidget(self.levelCombo, 5, 1, 1, 4)


        self.grid.addWidget(self.testLab, 6, 0) # 7
        self.grid.addWidget(self.testCombo, 6, 1, 1, 4)

        self.grid.addWidget(self.digitLab, 7, 0) # 8
        self.grid.addWidget(self.digitCombo, 7, 1, 1, 4)

        self.grid.addWidget(self.modelLab, 8, 0) # 9
        self.grid.addWidget(self.modelCombo, 8, 1, 1, 4)

        self.grid.addWidget(self.adjLab, 9, 0) # 10
        self.grid.addWidget(self.adjCombo, 9, 1, 1, 4)

        self.grid.addWidget(self.frqLab, 10, 0) # 11
        self.grid.addWidget(self.frqEdit, 10, 1, 1, 4)

        self.grid.addWidget(self.permLab, 11, 0) # 12
        self.grid.addWidget(self.permEdit, 11, 1, 1, 4)

        self.grid.addWidget(self.covLab, 12, 0) # 13
        self.grid.addWidget(self.covEdit, 12, 1, 1, 4)

        self.grid.addWidget(self.consensusLabel, 13, 0) # 14
        self.grid.addWidget(self.consensusCB, 13, 1,1,4)

        self.grid.addWidget(self.runButton, 14, 2) # 15
        ########################################################
        ### init (same as radio1_clicked)
        self.radio1.setChecked(True)
        self.frqEdit.setText('0.05')
        # self.oEdit.setText('output.txt')

        self.vLab.setEnabled(False)
        self.vEdit.setEnabled(False)
        self.vButton.setEnabled(False)
        self.vvButton.setEnabled(False)
        self.traitLab.setEnabled(True)
        self.traitCombo.setEnabled(True)
        self.traitCombo.clear()
        self.traitCombo.addItem("disease trait/case-control study")
        self.traitCombo.addItem("quantitative trait")
        self.levelLab.setEnabled(False)
        self.levelCombo.setEnabled(False)
        self.testLab.setEnabled(False)
        self.testCombo.setEnabled(False)
        self.digitLab.setEnabled(True)
        self.digitCombo.setEnabled(True)
        self.modelLab.setEnabled(False)
        self.modelCombo.setEnabled(False)
        self.adjLab.setEnabled(False)
        self.adjCombo.setEnabled(False)
        self.frqLab.setEnabled(False)
        self.frqEdit.setEnabled(False)
        self.permLab.setEnabled(False)
        self.permEdit.setEnabled(False)
        self.covLab.setEnabled(False)
        self.covEdit.setEnabled(False)
        self.consensusLabel.setEnabled(False)
        self.consensusCB.setEnabled(False)
        ### Events and signals
        self.connect(self.radio1, QtCore.SIGNAL("toggled(bool)"),self.radio1_clicked)
        self.connect(self.radio2, QtCore.SIGNAL("toggled(bool)"),self.radio2_clicked)
        self.connect(self.radio3, QtCore.SIGNAL("toggled(bool)"),self.radio3_clicked)
        self.connect(self.radio4, QtCore.SIGNAL("toggled(bool)"),self.radio4_clicked)
        self.connect(self.radio5, QtCore.SIGNAL("toggled(bool)"),self.radio4_clicked) # radio5 checked is the same (parameters used) as radio4 checked

        self.gButton.clicked.connect(self.gButtonClicked)
        self.gvButton.clicked.connect(self.gvButtonClicked)

        self.vButton.clicked.connect(self.vButtonClicked)
        self.vvButton.clicked.connect(self.vvButtonClicked)

        self.oButton.clicked.connect(self.oButtonClicked)

        self.connect(self.traitCombo, QtCore.SIGNAL('activated(QString)'), self.traitCombo_chosen)

        self.connect(self.levelCombo, QtCore.SIGNAL('activated(QString)'), self.levelCombo_chosen)

        self.connect(self.testCombo, QtCore.SIGNAL('activated(QString)'), self.testCombo_chosen)

        self.runButton.clicked.connect(self.runButtonClicked)
        #######################################################
        self.setLayout(self.grid)
        #######################################################
    def radio1_clicked(self, enabled):
        if enabled:
            self.vLab.setEnabled(False)
            self.vEdit.setEnabled(False)
            self.vButton.setEnabled(False)
            self.vvButton.setEnabled(False)
            self.traitLab.setEnabled(True)
            self.traitCombo.setEnabled(True)
            self.traitCombo.clear()
            self.traitCombo.addItem("disease trait/case-control study")
            self.traitCombo.addItem("quantitative trait")
            self.levelLab.setEnabled(False)
            self.levelCombo.setEnabled(False)
            self.testLab.setEnabled(False)
            self.testCombo.setEnabled(False)
            self.digitLab.setEnabled(True)
            self.digitCombo.setEnabled(True)
            self.modelLab.setEnabled(False)
            self.modelCombo.setEnabled(False)
            self.adjLab.setEnabled(False)
            self.adjCombo.setEnabled(False)
            self.frqLab.setEnabled(False)
            self.frqEdit.setEnabled(False)
            self.permLab.setEnabled(False)
            self.permEdit.setEnabled(False)
            self.covLab.setEnabled(False)
            self.covEdit.setEnabled(False)
            self.consensusLabel.setEnabled(False)
            self.consensusCB.setEnabled(False)
    def radio2_clicked(self, enabled):
        if enabled:
            self.vLab.setEnabled(False)
            self.vEdit.setEnabled(False)
            self.vButton.setEnabled(False)
            self.vvButton.setEnabled(False)
            self.traitLab.setEnabled(True)
            self.traitCombo.setEnabled(True)
            self.traitCombo.clear()
            self.traitCombo.addItem("disease trait/case-control study")
            self.traitCombo.addItem("quantitative trait")
            self.levelLab.setEnabled(True)
            self.levelCombo.setEnabled(True)
            self.testLab.setEnabled(True)
            self.testCombo.setEnabled(True)
            self.testCombo.clear()
            self.testCombo.addItem("Fisher's exact test")
            self.testCombo.addItem("Pearson chi-squared test")
            self.testCombo.addItem("Logistic regression")
            self.digitLab.setEnabled(True)
            self.digitCombo.setEnabled(True)
            self.modelLab.setEnabled(True)
            self.modelCombo.setEnabled(True)
            self.adjLab.setEnabled(True)
            self.adjCombo.setEnabled(True)
            self.frqLab.setEnabled(True)
            self.frqEdit.setEnabled(True)
            self.permLab.setEnabled(True)
            self.permEdit.setEnabled(True)
            self.covLab.setEnabled(False)
            self.covEdit.setEnabled(False)
            if self.levelCombo.currentText() == 'residue':
                self.consensusLabel.setEnabled(True)
                self.consensusCB.setEnabled(True)
            else:
                self.consensusLabel.setEnabled(False)
                self.consensusCB.setEnabled(False)
    def radio3_clicked(self, enabled):
        if enabled:
            self.vLab.setEnabled(False)
            self.vEdit.setEnabled(False)
            self.vButton.setEnabled(False)
            self.vvButton.setEnabled(False)
            self.traitLab.setEnabled(False)
            self.traitCombo.setEnabled(False)
            self.traitCombo.clear()
            self.traitCombo.addItem("disease trait/case-control study")
            self.traitCombo.addItem("quantitative trait")
            self.levelLab.setEnabled(False)
            self.levelCombo.setEnabled(False)
            self.testLab.setEnabled(False)
            self.testCombo.setEnabled(False)
            self.digitLab.setEnabled(False)
            self.digitCombo.setEnabled(False)
            self.modelLab.setEnabled(False)
            self.modelCombo.setEnabled(False)
            self.adjLab.setEnabled(False)
            self.adjCombo.setEnabled(False)
            self.frqLab.setEnabled(False)
            self.frqEdit.setEnabled(False)
            self.permLab.setEnabled(False)
            self.permEdit.setEnabled(False)
            self.covLab.setEnabled(False)
            self.covEdit.setEnabled(False)
            self.consensusLabel.setEnabled(True)
            self.consensusCB.setEnabled(True)
    def radio4_clicked(self, enabled):
        if enabled:
            self.vLab.setEnabled(False)
            self.vEdit.setEnabled(False)
            self.vButton.setEnabled(False)
            self.vvButton.setEnabled(False)
            self.traitLab.setEnabled(True)
            self.traitCombo.setEnabled(True)
            self.traitCombo.clear()
            self.traitCombo.addItem("disease trait/case-control study")
            self.levelLab.setEnabled(True)
            self.levelCombo.setEnabled(True)
            self.testLab.setEnabled(True)
            self.testCombo.setEnabled(True)
            self.testCombo.clear()
            self.testCombo.addItem("Fisher's exact test")
            self.testCombo.addItem("Pearson chi-squared test")
            self.digitLab.setEnabled(True)
            self.digitCombo.setEnabled(True)
            self.modelLab.setEnabled(False)
            self.modelCombo.setEnabled(False)
            self.adjLab.setEnabled(False)
            self.adjCombo.setEnabled(False)
            self.frqLab.setEnabled(True)
            self.frqEdit.setEnabled(True)
            self.permLab.setEnabled(False)
            self.permEdit.setEnabled(False)
            self.covLab.setEnabled(False)
            self.covEdit.setEnabled(False)
            if self.levelCombo.currentText() == 'residue':
                self.consensusLabel.setEnabled(True)
                self.consensusCB.setEnabled(True)
            else:
                self.consensusLabel.setEnabled(False)
                self.consensusCB.setEnabled(False)
    def gButtonClicked(self):
        self.gfile = QtGui.QFileDialog.getOpenFileName(self,'Open File', '.')
        self.gEdit.setText(self.gfile)
        if self.gfile:
            self.gvButton.setEnabled(True)
    def gvButtonClicked(self):
        self.gdialog = QtGui.QDialog()
        self.gdialog.resize(self.width, self.height)
        self.gdialog.setWindowTitle("view of genotype data")
        df  = pd.read_csv(str(self.gfile), delim_whitespace= True, index_col = None, header = None)
        gdatatable = QtGui.QTableWidget(parent=self.gdialog)
        gdatatable.setColumnCount(len(df.columns))
        gdatatable.setRowCount(len(df.index))
        labels = list(df.iloc[0])
        labels = labels[2:]
        lheader = ['IID', 'PHT']
        lindex = 0
        for label in labels:
            temp = label.split('*')
            tlabel = 'HLA-' + temp[0]
            if lindex == 0:
                tlabel += '-1'
                lindex = 1
            else:
                tlabel += '-2'
                lindex = 0
            lheader.append(tlabel)
        gdatatable.setHorizontalHeaderLabels(lheader)
        for i in range(len(df.index)):
            for j in range(len(df.columns)):
                gdatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iat[i, j])))
        gdatatable.resizeColumnsToContents()
        ghbox = QtGui.QHBoxLayout()
        ghbox.addWidget(gdatatable)
        self.gdialog.setLayout(ghbox)
        self.gdialog.show()
    def vButtonClicked(self):
        self.covfile = QtGui.QFileDialog.getOpenFileName(self,'Open File', '.')
        self.vEdit.setText(self.covfile)
        if self.covfile:
            self.vvButton.setEnabled(True)
    def vvButtonClicked(self):
        self.vdialog = QtGui.QDialog()
        self.vdialog.resize(self.width, self.height)
        self.vdialog.setWindowTitle("view of covariates data")
        df  = pd.read_csv(str(self.covfile), delim_whitespace= True, index_col = None, header = 0)
        vdatatable = QtGui.QTableWidget(parent=self.vdialog)
        vdatatable.setColumnCount(len(df.columns))
        vdatatable.setRowCount(len(df.index))
        vdatatable.setHorizontalHeaderLabels(list(df.columns.values))
        for i in range(len(df.index)):
            for j in range(len(df.columns)):
                vdatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iat[i, j])))
        vdatatable.resizeColumnsToContents()
        vhbox = QtGui.QHBoxLayout()
        vhbox.addWidget(vdatatable)
        self.vdialog.setLayout(vhbox)
        self.vdialog.show()
    def oButtonClicked(self):
        self.outfile = QtGui.QFileDialog.getSaveFileName(self,'Open File', '.')
        self.oEdit.setText(self.outfile)
    def traitCombo_chosen(self, text):
        if self.radio2.isChecked(): # for assoc only
            if text == 'quantitative trait':
                self.vLab.setEnabled(True)
                self.vEdit.setEnabled(True)
                self.vButton.setEnabled(True)
                self.vvButton.setEnabled(True)
                self.levelLab.setEnabled(False)
                self.levelCombo.setEnabled(True)
                self.levelCombo.clear()
                self.levelCombo.addItems(["allele"])
                self.levelCombo.setEnabled(False)
                self.testLab.setEnabled(True)
                self.testCombo.setEnabled(True)
                self.testCombo.clear()
                self.testCombo.addItem("Linear regression")
                self.modelLab.setEnabled(True)
                self.modelCombo.setEnabled(True)
                self.modelCombo.clear()
                self.modelCombo.addItems(["additive", "dominant", "recessive"])
                self.covLab.setEnabled(True)
                self.covEdit.setEnabled(True)
            else:
                self.vLab.setEnabled(False)
                self.vEdit.setEnabled(False)
                self.vButton.setEnabled(False)
                self.vvButton.setEnabled(False)
                self.levelLab.setEnabled(True)
                self.levelCombo.setEnabled(True)
                self.levelCombo.clear()
                self.levelCombo.addItems(["allele", "residue"])
                self.testLab.setEnabled(True)
                self.testCombo.setEnabled(True)
                self.testCombo.clear()
                self.testCombo.addItem("Fisher's exact test")
                self.testCombo.addItem("Pearson chi-squared test")
                if str(self.levelCombo.currentText()) == 'allele':
                    self.testCombo.addItem("Logistic regression")
                self.modelLab.setEnabled(True)
                self.modelCombo.setEnabled(True)
                self.modelCombo.clear()
                # self.modelCombo.addItems(["allelic", "dominant", "recessive"])
                self.modelCombo.addItems(["allelic", "dom", "rec"])
                self.covLab.setEnabled(False)
                self.covEdit.setEnabled(False)
    def levelCombo_chosen(self, text):
        if self.radio2.isChecked(): # for assoc only
            if text == 'allele':
                self.testCombo.clear()
                self.testCombo.addItem("Fisher's exact test")
                self.testCombo.addItem("Pearson chi-squared test")
                self.testCombo.addItem("Logistic regression")
                self.digitLab.setEnabled(True)
                self.digitCombo.setEnabled(True)
                self.modelLab.setEnabled(True)
                self.modelCombo.setEnabled(True)
                self.modelCombo.clear()
                # self.modelCombo.addItems(["allelic", "dominant", "recessive"])
                self.modelCombo.addItems(["allelic", "dom", "rec"])
                self.adjLab.setEnabled(True)
                self.adjCombo.setEnabled(True)
                self.frqLab.setEnabled(True)
                self.frqEdit.setEnabled(True)
                self.permLab.setEnabled(True)
                self.permEdit.setEnabled(True)
                self.consensusLabel.setEnabled(False)
                self.consensusCB.setEnabled(False)
                self.covLab.setEnabled(False)
                self.covEdit.setEnabled(False)
            else: # residue
                self.testCombo.clear()
                self.testCombo.addItem("Fisher's exact test")
                self.testCombo.addItem("Pearson chi-squared test")
                self.digitLab.setEnabled(False)
                self.digitCombo.setEnabled(False)
                self.modelLab.setEnabled(False)
                self.modelCombo.setEnabled(False)
                self.adjLab.setEnabled(False)
                self.adjCombo.setEnabled(False)
                self.frqLab.setEnabled(False)
                self.frqEdit.setEnabled(False)
                self.permLab.setEnabled(False)
                self.permEdit.setEnabled(False)
                self.covLab.setEnabled(False)
                self.covEdit.setEnabled(False)
                self.consensusLabel.setEnabled(True)
                self.consensusCB.setEnabled(True)
        elif self.radio4.isChecked() or self.radio5.isChecked(): # zyg & inter
            if text == 'allele':
                self.digitLab.setEnabled(True)
                self.digitCombo.setEnabled(True)
                self.frqLab.setEnabled(True)
                self.frqEdit.setEnabled(True)
                self.consensusLabel.setEnabled(False)
                self.consensusCB.setEnabled(False)
            else: # residue
                self.digitLab.setEnabled(False)
                self.digitCombo.setEnabled(False)
                self.frqLab.setEnabled(False)
                self.frqEdit.setEnabled(False)
                self.consensusLabel.setEnabled(True)
                self.consensusCB.setEnabled(True)
    def testCombo_chosen(self, text):
        if self.radio2.isChecked(): # for assoc only
            if text == 'Logistic regression' or text == 'Linear regression':
                self.modelLab.setEnabled(True)
                self.modelCombo.setEnabled(True)
                self.modelCombo.clear()
                # self.modelCombo.addItems(["additive", "dominant", "recessive"])
                self.modelCombo.addItems(["additive", "dom", "rec"])
                self.vLab.setEnabled(True)
                self.vEdit.setEnabled(True)
                self.vButton.setEnabled(True)
                self.covLab.setEnabled(True)
                self.covEdit.setEnabled(True)
            else:
                if str(self.levelCombo.currentText()) == 'allele':
                    self.modelLab.setEnabled(True)
                    self.modelCombo.setEnabled(True)
                    self.modelCombo.clear()
                    # self.modelCombo.addItems(["allelic", "dominant", "recessive"])
                    self.modelCombo.addItems(["allelic", "dom", "rec"])
                else:
                    self.modelLab.setEnabled(False)
                    self.modelCombo.setEnabled(False)
                self.vLab.setEnabled(False)
                self.vEdit.setEnabled(False)
                self.vButton.setEnabled(False)
                self.vvButton.setEnabled(False)
                self.covLab.setEnabled(False)
                self.covEdit.setEnabled(False)
    def runButtonClicked(self):
        if self.radio1.isChecked():
            comm = 'python PyHLA.py --summary --input '
            comm += str(self.gfile)
            digi = str(self.digitCombo.currentText())
            comm += ' --digit '
            comm += digi
            comm += ' --out '
            comm += str(self.outfile)
            os.system(comm)
        elif self.radio2.isChecked(): # for assoc only
            if self.levelCombo.currentText() == 'allele':
                comm = 'python PyHLA.py --assoc --input '
                comm += str(self.gfile)
                digi = str(self.digitCombo.currentText())
                comm += ' --digit '
                comm += digi
                testStr = str(self.testCombo.currentText())
                testM = ''
                if testStr == 'Pearson chi-squared test':
                    testM = 'chisq'
                elif testStr == "Fisher's exact test":
                    testM = 'fisher'
                elif testStr == 'Logistic regression':
                    testM = 'logistic'
                elif testStr == 'Linear regression':
                    testM = 'linear'
                comm += ' --test '
                comm += testM
                model = str(self.modelCombo.currentText())
                comm += ' --model '
                comm += model
                freq = str(self.frqEdit.text())
                comm += ' --freq '
                comm += freq
                adjust = str(self.adjCombo.currentText())
                comm += ' --adjust '
                comm += adjust
                if testM == 'logistic' or testM == 'linear':
                    if self.covfile:
                        comm += ' --covar '
                        comm += str(self.covfile)
                        if str(self.covEdit.text()):
                            comm += ' --covarname '
                            comm += str(self.covEdit.text())
                if str(self.permEdit.text()) and int(str(self.permEdit.text())) > 0:
                    comm += ' --perm '
                    comm += str(self.permEdit.text())
                comm += ' --out '
                comm += str(self.outfile)
                os.system(comm)
            else: # residue
                comm = 'python PyHLA.py --assoc-AA --input '
                comm += str(self.gfile)
                testStr = str(self.testCombo.currentText())
                testM = ''
                if testStr == 'Pearson chi-squared test':
                    testM = 'chisq'
                elif testStr == "Fisher's exact test":
                    testM = 'fisher'
                comm += ' --test '
                comm += testM
                if self.consensusCB.isChecked():
                    comm += ' --consensus'
                comm += ' --out '
                comm += str(self.outfile)
                os.system(comm)
        elif self.radio3.isChecked():
            comm = 'python PyHLA.py --align --input '
            comm += str(self.gfile)
            if self.consensusCB.isChecked():
                comm += ' --consensus'
            comm += ' --out '
            comm += str(self.outfile)
            os.system(comm)
        elif self.radio4.isChecked():
            comm = 'python PyHLA.py --zygosity --input '
            comm += str(self.gfile)
            comm += ' --level '
            comm += str(self.levelCombo.currentText())
            testStr = str(self.testCombo.currentText())
            testM = ''
            if testStr == 'Pearson chi-squared test':
                testM = 'chisq'
            elif testStr == "Fisher's exact test":
                testM = 'fisher'
            comm += ' --test '
            comm += testM
            if self.levelCombo.currentText() == 'allele':
                digi = str(self.digitCombo.currentText())
                comm += ' --digit '
                comm += digi
                freq = str(self.frqEdit.text())
                comm += ' --freq '
                comm += freq
            else:
                if self.consensusCB.isChecked():
                    comm += ' --consensus'
            comm += ' --out '
            comm += str(self.outfile)
            os.system(comm)
        elif self.radio5.isChecked():
            comm = 'python PyHLA.py --interaction --input '
            comm += str(self.gfile)
            comm += ' --level '
            comm += str(self.levelCombo.currentText())
            testStr = str(self.testCombo.currentText())
            testM = ''
            if testStr == 'Pearson chi-squared test':
                testM = 'chisq'
            elif testStr == "Fisher's exact test":
                testM = 'fisher'
            comm += ' --test '
            comm += testM
            if self.levelCombo.currentText() == 'allele':
                digi = str(self.digitCombo.currentText())
                comm += ' --digit '
                comm += digi
                freq = str(self.frqEdit.text())
                comm += ' --freq '
                comm += freq
            else:
                if self.consensusCB.isChecked():
                    comm += ' --consensus'
            comm += ' --out '
            comm += str(self.outfile)
            os.system(comm)
        ### show results
        self.odialog = QtGui.QDialog()
        self.odialog.resize(self.width, self.height)
        self.odialog.setWindowTitle("view of output")
        if self.radio4.isChecked() or self.radio5.isChecked() or (self.radio2.isChecked() and testM != 'raw') or (self.radio2.isChecked() and testM != 'score'):
            df  = pd.read_csv(str(self.outfile), delim_whitespace= True, index_col = None, header = 0)
            odatatable = QtGui.QTableWidget(parent=self.odialog)
            odatatable.setColumnCount(len(df.columns))
            odatatable.setRowCount(len(df.index))
            odatatable.setHorizontalHeaderLabels(list(df.columns.values))
            for i in range(len(df.index)):
                for j in range(len(df.columns)):
                    odatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iat[i, j])))
            odatatable.resizeColumnsToContents()
            ohbox = QtGui.QHBoxLayout()
            ohbox.addWidget(odatatable)
            self.odialog.setLayout(ohbox)
            self.odialog.show()
        else:
            odatatable = QtGui.QTextEdit(parent=self.odialog)
            odatatable.setText('')
            fff = open(str(self.outfile))
            for line in fff:
                odatatable.append(line)
            fff.close()
            odatatable.setAlignment(QtCore.Qt.AlignLeft)
            ohbox = QtGui.QHBoxLayout()
            ohbox.addWidget(odatatable)
            self.odialog.setLayout(ohbox)
            self.odialog.show()
#####################################################

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    gui = PyHLAWin()
    gui.show()
    sys.exit(app.exec_())

