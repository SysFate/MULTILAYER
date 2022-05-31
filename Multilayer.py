__author__ = 'Julien Moehlin'

import os
from tkinter import *
import tkinter as tk
import tkinter.messagebox
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename 
from tkinter import ttk
import numpy as np
import pandas as pd
import time
import math
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
import scipy.stats as stats
import seaborn as sb
import threading
import multiprocessing
from sklearn.cluster import AgglomerativeClustering
import networkx as nx
import community as louvain
import webbrowser
import copy
import PIL
from PIL import ImageTk
from PIL import Image
from PIL import ImageDraw
import sys
import re

class GeneSample(object):
	"""
	Gene with information from the whole sample.
	GeneSample(geneID, totalCount, maxCount, minCount, dic_coordinate_count)
		geneID : 					Gene ID
		totalCount : 				The sum of counts for all coordinates
		maxCount : 					Highest counts
		minCount : 					Lower counts
			##numberUp :				Number of gexel where this gene is up regulated
			##numberDown : 				Number of gexel where this gene is down regulated
		dic_coordinate_count :		Dictionnary[coordinate ID] = count. This dic has as key the name of coordinate and as value the count for this coordinate.
		dic_pattern :				Dictionnary[index pattern] = [coordinates]. This dic has as key a number (pattern index) and as value the list of coordinates.
	"""
	def __init__(self, geneID):
		self.geneID = geneID
		self.totalCount = 0
		self.maxCount = 0
		self.minCount = 0
		self.dic_coordinate_count = {}
		self.dic_pattern = {}

class Gexel(object):
	"""
	Gexel(gexelID, totalGeneCount, genesValueGexel)
		gexelID :			Gexel name
		totalGeneCount : 	Total number of counts for all genes
		genesValueGexel :	Dictionary[gene] = count 
	"""
	def __init__(self, gexelID):
		self.gexelID = gexelID
		self.totalGeneCount = 0
		self.genesValueGexel = {}

class Run(object):
	"""
	Run has all parameters for the analysis provide by the user
	"""
	def __init__(self):
		self.sizeX = 0
		self.sizeY = 0
		self.totalGexel = 0
		self.mode = None
		self.raw = ['']
		self.norm = ['']
		self.diff = ['']
		self.transpose = False
		self.round = False
		self.method = None
		self.indexGetParameters = None
		self.dicSample = {}
		self.numberCPU = None
		self.indexRun = None
		self.minPattern = None
		self.saveMatrix = None
		self.upDiffThreshold = 1
		self.downDiffThreshold = -1
		self.calculatedX = 0
		self.calculatedY = 0

class Sample(object):
	"""
	Sample 
	"""
	def __init__(self, name):
		self.name = name
		self.raw = None
		self.norm = None
		self.diff = None
		self.numberGexelCreate = None
		self.maxCountRaw = None
		self.minCountRaw = None
		self.maxCountNorm = None
		self.minCountNorm = None
		self.dicGexelRaw = None
		self.dicSampleRaw = None
		self.dicGexelNorm = None
		self.dicSampleNorm = None
		self.listGexelDiff = None
		self.dicSampleDiff = None
		self.infoPattern = []
		self.infoUp = []
		self.infoDown = []
		self.simMatrix = None
		self.matrixLouvain = None
		self.currentLouvain = None

class StartGUI(tk.Tk):
	"""
	main window
	"""
	def __init__(self, runObject):
		tk.Tk.__init__(self)
		self.minsize(580, 300)
		self.title('MULTILAYER - Molecular Tissu Digitalization Analyser')
		self.container = tk.Frame(self)
		self.container.pack(expand=1)
		self.dic_frames = {}
		for frame_name in [StartPageGUI, InputGUIExplore, InputGUIBatch, ParamGUI]:
			page_name = frame_name.__name__
			frame = frame_name(self.container, self, runObject)
			self.dic_frames[frame_name.__name__] = frame
			frame.grid(row=0, column=0, sticky="nsew")
		self.show_frame("StartPageGUI", runObject)

	def show_frame(self, page_name, runObject):
		"""
		Function for switch panel
		"""
		frame = self.dic_frames[page_name]
		frame.tkraise()

	def roundVal(self, x):
		"""
		Round values
		"""
		new = x.split('x')
		return(f'{round(float(new[0]))}x{round(float(new[1]))}')

	def create_matrix(self, sampleData, runObject):
		"""
		Prepare raw, norm, diff. genes and patterns
		"""
		sample = Sample(os.path.basename(sampleData))
		dirNamePath = os.path.join(os.path.dirname(os.path.abspath(runObject.raw[0])), 'Similarity_Matrix_Sample')
		if not os.path.exists(dirNamePath):
			os.mkdir(dirNamePath)
		runObject.dicSample[sample.name] = sample
		if runObject.transpose == True:
			sample.raw = pd.read_csv(sampleData, sep='\t', index_col=[0]).astype(float).transpose()
		else:
			sample.raw = pd.read_csv(sampleData, sep='\t', index_col=[0]).astype(float)
		if runObject.round == True:
			sample.raw.columns = sample.raw.columns.to_series().apply(self.roundVal)
		dic_gexel_raw = {}
		listGexelCountRaw = []
		sample.raw = sample.raw.loc[(sample.raw.sum(axis=1) != 0), (sample.raw.sum(axis=0) != 0)]
		sample.raw.index = sample.raw.index.str.upper()
		if runObject.mode == 'explore':
			dic_gexel_raw = dict(sample.raw.apply(lambda x : self.createGexel(x, listGexelCountRaw, runObject), axis=0))
			sample.maxCountRaw =  max(listGexelCountRaw)
			sample.minCountRaw =  min(listGexelCountRaw)
			dic_geneSample_raw = dict(sample.raw.apply(self.createGeneSample, axis=1))
			runObject.dicSample[sample.name].dicGexelRaw = dic_gexel_raw
			runObject.dicSample[sample.name].dicSampleRaw = dic_geneSample_raw
		listGexelCountNorm = []
		if runObject.mode == 'cluster':
			self.normalization(sampleData, sample, runObject, dirNamePath)
		else:
			if runObject.norm == [''] and runObject.diff == ['']:
				self.normalization(sampleData, sample, runObject, dirNamePath)
			elif runObject.norm != ['']:
				print(f'### - Loading Normalisation file - {sample.name}')
				if runObject.transpose == True:
					sample.norm = pd.read_csv(runObject.norm[0], sep='\t', index_col=[0]).astype(float).transpose()
					sample.norm.index = sample.norm.index.str.upper()
					if runObject.diff == ['']:
						sample.diff = pd.read_csv(runObject.norm[0], sep='\t', index_col=[0]).astype(float).transpose()
						sample.diff += 1
				else:
					sample.norm = pd.read_csv(runObject.norm[0], sep='\t', index_col=[0]).astype(float)
					sample.norm.index = sample.norm.index.str.upper()
					if runObject.diff == ['']:
						sample.diff = pd.read_csv(runObject.norm[0], sep='\t', index_col=[0]).astype(float)
						sample.diff += 1
		if runObject.mode == 'explore' and runObject.norm != ['']:
			dic_geneSampleNorm = dict(runObject.dicSample[sample.name].norm.apply(self.createGeneSample, axis=1))
			runObject.dicSample[sample.name].dicSampleNorm = dic_geneSampleNorm
			dic_gexelNorm = dict(runObject.dicSample[sample.name].norm.apply(lambda x : self.createGexel(x, listGexelCountNorm, runObject), axis=0))
			runObject.dicSample[sample.name].dicGexelNorm = dic_gexelNorm
			sample.maxCountNorm =  max(listGexelCountNorm)
			sample.minCountNorm =  min(listGexelCountNorm)
		#print(f'Diff gene {sampleData}')
		self.diffgene(sample, runObject, dirNamePath)
		runObject.dicSample[sample.name].listGexelDiff	= list(runObject.dicSample[sample.name].diff.columns.values)
		dic_geneSampleDiff = dict(runObject.dicSample[sample.name].diff.apply(self.createGeneSample, axis=1))
		runObject.dicSample[sample.name].dicSampleDiff = dic_geneSampleDiff
		#print(f'Patterns detection {sampleData}')
		self.agglomerative(dic_geneSampleDiff, sample, runObject,  runObject.minPattern)
		print(f'### - Similarity Matrix {sample.name}')
		self.generateMatrixSimilarity(dic_geneSampleDiff, runObject, dirNamePath, sample.name)

	def clusterMatrix(self, runObject):
		startTime = time.time()
		to_process = []
		listThread = []
		index = 0
		indexStr = ''
		if index <= len(runObject.raw):
			index = len(runObject.raw)
			indexStr = 'raw'
		elif index <= len(runObject.norm):
			index = len(runObject.norm)
			indexStr = 'norm'
		elif index <= len(runObject.diff):
			index = len(runObject.diff)
			indexStr = 'diff'
		while len(to_process) < index:
			if threading.active_count() <= runObject.numberCPU:
				file = runObject.raw[len(to_process)]
				newThread = threading.Thread(target=self.create_matrix, args=[file, runObject])
				newThread.start()
				listThread.append(newThread)
				to_process.append(file)
			else:
				time.sleep(0.01)
			for thread in listThread:
				thread.join()
		stopTime = time.time()
		print(f'Time {stopTime-startTime} seconds')
		runObject.indexRun = 1
		self.destroy()

	def destroyFrame(self, run_Object):
		runObject = run_Object
		self.destroy()

	def createGexel(self, x, listGexelCount, run):
		"""
		Create Gexel object (column from matrix)
		"""
		sumCount = sum(x)
		if sumCount != 0:
			gexel = Gexel(x.name)
			if int(x.name.split('x')[0]) > run.calculatedX:
				run.calculatedX = int(x.name.split('x')[0])
			if int(x.name.split('x')[1]) > run.calculatedY:
				run.calculatedY = int(x.name.split('x')[1])
			gexel.totalGeneCount = sumCount
			gexel.genesValueGexel = x.to_dict()
			listGexelCount.append(sumCount)
			return gexel

	def createGeneSample(self, x):
		"""
		Create GeneSample object (row from matrix)
		"""
		geneSample = GeneSample(x.name.upper())
		geneSample.maxCount = max(x)
		geneSample.minCount = min(x)
		geneSample.dic_coordinate_count = x.to_dict()
		return geneSample

	def agglomerative(self, dicGeneSample, sample, runObject, indexMinPattern):
		"""
		Agglomerative clustering for define contigous pattern
		"""
		print(f'### - Patterns detection {sample.name}')
		sample.infoPattern = []
		for geneSample in dicGeneSample.keys():
			tempListCoordinate = []
			indexFind = 0
			for coordinate in dicGeneSample[geneSample].dic_coordinate_count.keys():
				if dicGeneSample[geneSample].dic_coordinate_count[coordinate] >= runObject.upDiffThreshold:
					x = int(coordinate.split('x')[0])
					y = int(coordinate.split('x')[1])
					tempListCoordinate.append([x, y])
			tempListCoordinate = np.asarray(tempListCoordinate)
			if len(tempListCoordinate) >= runObject.minPattern:
				clustering = AgglomerativeClustering(n_clusters=None, affinity='euclidean',
					linkage='single', distance_threshold=1.5)
				labels = clustering.fit_predict(tempListCoordinate)			
				(unique, counts) = np.unique(labels, return_counts=True)
				frequencies = np.asarray((unique, counts)).T
				for i in frequencies:
					index = i[0]
					count = i[1]
					if count >= int(indexMinPattern):
						indexFind += 1
						listTemp = []
						for o in tempListCoordinate[labels==index]:
							listTemp.append(f'{o[0]}x{o[1]}')
						dicGeneSample[geneSample].dic_pattern[indexFind] = listTemp
						sample.infoPattern.append((len(listTemp), dicGeneSample[geneSample].geneID))
	
	def generateMatrixSimilarity(self, dicGeneSample, runObject, pathDirSave, name):
		"""
		Compare all patterns and generate a similarity matrix
		"""
		startTime_sim = time.time()
		list_similarity = []
		list_column = []
		for geneQuery in dicGeneSample.keys():
			for numberPatternQuery in dicGeneSample[geneQuery].dic_pattern.keys():
				set_gene_query = set(dicGeneSample[geneQuery].dic_pattern[int(numberPatternQuery)])
				list_temp = []
				list_column.append(f'{geneQuery} | {numberPatternQuery}')
				index_find = 0
				for geneComp in dicGeneSample.keys():
					for numberPatternComp in dicGeneSample[geneComp].dic_pattern.keys():
						if geneQuery == geneComp:
							if numberPatternQuery == numberPatternComp:
								list_temp.append(100)
							else:
								list_temp.append(np.NaN)
						else:
							if runObject.method == 'Percent':
								inter = len(set_gene_query.intersection(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)])))
								similarityValue = float(inter/len(set_gene_query))*100
								if round(similarityValue, 2) >= 0:
									list_temp.append(round(similarityValue, 2))
									index_find += 1
								else:
									list_temp.append(np.NaN)
							elif runObject.method == 'Tanimoto':
								inter = len(set_gene_query.intersection(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)])))
								union = len(set_gene_query.union(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)])))
								similarityValue = float(inter/union)*100
								if round(similarityValue, 2) >= 0:
									list_temp.append(round(similarityValue, 2))
									index_find += 1
								else:
									list_temp.append(np.NaN)
							elif runObject.method == 'Dice':
								inter = len(set_gene_query.intersection(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)])))
								union = len(set_gene_query.union(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)])))
								similarityValue = float((2*inter)/(len(set_gene_query) + len(set(dicGeneSample[geneComp].dic_pattern[int(numberPatternComp)]))))*100
								if round(similarityValue, 2) >= 0:
									list_temp.append(round(similarityValue, 2))
									index_find += 1
								else:
									list_temp.append(np.NaN)
				list_similarity.append(list_temp)			
		dfSim = pd.DataFrame(list_similarity, columns=list_column, index=list_column)
		dfSim = dfSim.transpose()
		runObject.dicSample[name].simMatrix = dfSim.copy()
		dfSim.to_csv(os.path.join(pathDirSave, f'{name}_size{runObject.minPattern}_{runObject.method}_similarity_matrix.tsv'), sep='\t')
		stopTime_sime = time.time()
		print(f'Time to make similarity  {stopTime_sime-startTime_sim} seconds')
		StartGUI.matrixForLouvain(StartGUI, runObject.dicSample[name])

	def matrixForLouvain(self, sample):
		"""
		Create file for Louvain analysis
		"""
		dico_tot ={}
		sample.simMatrix.apply(lambda x: StartGUI.dico_df(StartGUI, x, dico_tot))
		dico_fin = StartGUI.dico_final(StartGUI, dico_tot)
		df = pd.DataFrame(dico_fin, columns=['TF', 'TG', 'Weight'])
		sample.matrixLouvain = df

	def dico_df(self, x, dic):
		dic[x.name] = x[x.notna()].to_dict()

	def dico_final(self, dico):
		dico_final = []
		for i in dico.keys():
			for j in dico[i].keys():
				dico_final.append([i, j, dico[i][j]])
		return dico_final

	def normalization(self, data, sample, runObject, dirNamePath):
		"""
		Quantile normalization
		"""
		print(f'### - Normalisation - {sample.name}')
		if runObject.transpose == True:
			df = pd.read_csv(data, sep='\t', index_col=[0]).astype(float).transpose()
		else:
			df = pd.read_csv(data, sep='\t', index_col=[0]).astype(float)
		if runObject.round == True:
			df.columns = df.columns.to_series().apply(self.roundVal)
		df = df.loc[(df.sum(axis=1) != 0), (df.sum(axis=0) != 0)]
		df += 1
		rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
		dfQN = df.rank(method='first').stack().astype(int).map(rank_mean).unstack()
		dfDiff = dfQN.copy()
		dfQN -= 1
		sample.norm = dfQN
		sample.diff = dfDiff
		sample.norm.index = sample.norm.index.str.upper()
		runObject.norm.append(sample.norm)
		if runObject.saveMatrix == True:
			print('SAVE NORM MATRIX')
			print(os.path.join(dirNamePath, f'{sample.name}_size{runObject.minPattern}_{runObject.method}_Norm.tsv'))
			sample.norm.to_csv(os.path.join(dirNamePath, f'{sample.name}_size{runObject.minPattern}_{runObject.method}_Norm.tsv'), sep='\t')

	def diffgene(self, sample, runObject, dirNamePath):
		"""
		Differentially expressed genes
		"""
		print(f'### - Differentially expression - {sample.name}')
		if runObject.diff == [''] or runObject.mode == 'cluster':
			sample.diff['mean'] = sample.diff.mean(axis=1)
			sample.diff = sample.diff.drop(['mean'], axis=1).apply(lambda x: np.log2(x/sample.diff['mean']))
		else:
			print(f'### - Loading Differentially expression file - {sample.name}')
			if runObject.transpose == True:
				sample.diff = pd.read_csv(runObject.diff[0], sep='\t', index_col=[0]).astype(float).transpose()
			else:
				sample.diff = pd.read_csv(runObject.diff[0], sep='\t', index_col=[0]).astype(float)
		df = sample.diff.copy()
		sample.diff.index = sample.diff.index.str.upper()
		df['positive'] = df[df >= runObject.upDiffThreshold].count(axis=1).values
		df['negative'] = df[df <= runObject.downDiffThreshold].count(axis=1).values
		runObject.dicSample[sample.name].infoUp = list(zip(df['positive'].values.tolist(), df.index.values.tolist()))
		runObject.dicSample[sample.name].infoDown = list(zip(df['negative'].values.tolist(), df.index.values.tolist()))
		df.drop(['positive'], axis=1)
		df.drop(['negative'], axis=1)
		if runObject.saveMatrix == True:
			print('SAVE DIFF MATRIX')
			sample.diff.to_csv(os.path.join(dirNamePath, f'{sample.name}_size{runObject.minPattern}_{runObject.method}_Diff.tsv'), sep='\t')
		runObject.diff.append(sample.diff)

class StartPageGUI(tk.Frame):
	"""
	Home page
	"""
	def __init__(self, parent, controller, runObject):
		tk.Frame.__init__(self, parent)
		self.controller = controller
		self.runObject = runObject
		self.canvasLogo = Canvas(self, width=580, height=95)
		self.canvasLogo.grid(row=1, columnspan=3, sticky=NW, pady=0)
		self.logo = PhotoImage(file=os.path.join('logo', 'multilayer2.png'))
		self.canvasLogo.create_image(0, 0, anchor=NW, image=self.logo)
		self.canvasLogo.bind("<Button-1>", lambda event : webbrowser.open_new("https://www.sysfate.org/"))
		self.canvasLogo.bind("<Enter>",  lambda event : self.canvasLogo.config(cursor="hand2"))
		self.canvasLogo.bind("<Leave>", lambda event : self.canvasLogo.config(cursor=""))
		self.labelExplore = tk.Label(self, text="Explore one dataset :", font= "Arial 14")#, font=controller.title_font)
		self.labelExplore.grid(row=2, column=1, sticky=NW, pady=30)
		self.buttonExplore = ttk.Button(self, width=20, text="Explore Mode", command=lambda: self._next('explore'))
		self.buttonExplore.grid(row=2, column=2, sticky=NW, pady=30)
		self.labelBatch = tk.Label(self, text="Explore multiple datasets :", font= "Arial 14")#, font=controller.title_font)
		self.labelBatch.grid(row=3, column=1, sticky=NW)
		self.buttonBatch = ttk.Button(self, width=20, text="Batch Mode", command=lambda: self._next('cluster'))
		self.buttonBatch.grid(row=3, column=2, sticky=NW)
		self.labelProvide = Label(self, text="This tool is provided by :", font= "Arial 8")
		self.labelProvide.grid(row=4, column=1, sticky=SE, pady=30)
		self.labelSysfate = Label(self, text="SysFate Team", fg="forestgreen", font= "Arial 8")
		self.labelSysfate.bind("<Button-1>", lambda event : webbrowser.open_new("https://www.sysfate.org/"))
		self.labelSysfate.bind("<Enter>",  lambda event : self.labelSysfate.configure(font = "Arial 8 underline", fg="red"))
		self.labelSysfate.bind("<Leave>", lambda event : self.labelSysfate.configure(font = "Arial 8", fg="forestgreen"))
		self.labelSysfate.grid(row=4, column=2, sticky=SW, pady=30)

	def _validate(self, P):
		if str.isdigit(P) or P == '':
			return True
		else:
			return False

	def _next(self, mode):
		self.runObject.mode = mode
		if mode == 'explore':
			self.controller.show_frame("InputGUIExplore", self.runObject)
			print('explore')
		elif mode == 'cluster':
			self.controller.show_frame("InputGUIBatch", self.runObject)
			print('cluster')

class InputGUIBatch(tk.Frame):
	"""
	Input Batch
	"""
	def __init__(self, parent, controller, runObject):
		tk.Frame.__init__(self, parent)
		self.controller = controller
		self.runObject = runObject
		self.onlyDigit = (self.register(self._validate))
		self.labelInput = tk.Label(self, text='Select dataset(s) :', font= "Arial 14")
		self.labelInput.pack(side=TOP)
		self.frameInputFile = tk.Frame(self)
		self.listRawFile = []
		self.frameMasterRaw = tk.Frame(self)
		self.frameRaw = tk.Frame(self.frameMasterRaw)
		self.varRaw = IntVar()
		self.namePathRaw = StringVar()
		self.rawBut = ttk.Button(self.frameRaw, text='Raw Matrix',  width = 20, command=lambda: self.browse_button(self.namePathRaw, 'raw'))
		self.rawBut.grid(row=1, column=1)
		self.frameRaw.grid(row=1, column=0, sticky="w")
		self.frameMasterRaw.pack(pady=5)
		self.frameInformation = tk.Frame(self)
		self.labelInformation = tk.Label(self.frameInformation, text='Matrix Informations :', font= "Arial 14")
		self.labelInformation.grid(row=0, columnspan=3)
		self.varTranspose = IntVar()
		self.checkTranspose = Checkbutton(self.frameInformation, width=2, variable=self.varTranspose)
		self.checkTranspose.deselect()
		self.checkTranspose.grid(row=1, column=0, sticky="w")
		self.labelTranspose = tk.Label(self.frameInformation, text=': Transpose matrix')
		self.labelTranspose.grid(row=1, column=1, sticky="w")
		self.varRound = IntVar()
		self.checkRound = Checkbutton(self.frameInformation, width=2, variable=self.varRound)
		self.checkRound.deselect()
		self.checkRound.grid(row=2, column=0, sticky="w")
		self.labelRound = tk.Label(self.frameInformation, text=': Round matrix')
		self.labelRound.grid(row=2, column=1, sticky="w")
		self.frameInformation.pack()
		self.frameInformationSize = tk.Frame(self)
		self.labelX = tk.Label(self.frameInformationSize, text='Size X of matrix :')
		self.labelX.grid(row=1, column=0, sticky="w")
		self.sizeXVar = IntVar()
		self.sizeX = Entry(self.frameInformationSize, textvar=self.sizeXVar, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.sizeXVar.set(32)
		self.sizeX.grid(row=1, column=1, sticky="w")
		self.labelY = tk.Label(self.frameInformationSize, text='Size Y of matrix:')
		self.labelY.grid(row=2, column=0, sticky="w")
		self.sizeYVar = IntVar()
		self.sizeY = Entry(self.frameInformationSize, textvar=self.sizeYVar, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.sizeYVar.set(32)
		self.sizeY.grid(row=2, column=1, sticky="w")
		self.frameInformationSize.pack()
		self.frameBackNext = Frame(self)
		self.buttonBack = ttk.Button(self.frameBackNext, width=25, text="Back", command=lambda: controller.show_frame('StartPageGUI', self.runObject))
		self.buttonBack.pack(side=LEFT, fill=X)
		self.buttonNext = ttk.Button(self.frameBackNext, width=25, text="Next", command=lambda: self._next())
		self.buttonNext.pack(side=RIGHT, fill=X)
		self.frameBackNext.pack(side=BOTTOM, fill=X)

	def _validate(self, P):
		if str.isdigit(P) or P == '':
			return True
		else:
			return False

	def browse_button(self, namePath, typeOfMatrix):
		listFile = []
		if self.runObject.mode == 'cluster':
			self.fileDir = askdirectory()
			namePath.set(self.fileDir)
			if self.fileDir != '':
				for file in os.listdir(self.fileDir):
					fileName, fileExtension = os.path.splitext(file)
					if fileExtension == '.tsv' or fileExtension == '.csv':
						listFile.append(os.path.join(self.fileDir, file))
		self.runObject.raw = listFile

	def _next(self):
		if self.runObject.raw[0] == '':
			tk.messagebox.showwarning('No raw file(s) selected !', 'You have to select at least one raw matrix.', icon='warning')
		else:
			try:
				if self.varTranspose.get() != 0:
					self.runObject.transpose = True
				else:
					self.runObject.transpose = False
				if self.varRound.get() != 0:
					self.runObject.round = True
				else:
					self.runObject.round = False
				self.runObject.sizeX = self.sizeXVar.get()
				self.runObject.sizeY = self.sizeYVar.get()
				self.runObject.totalGexel = self.runObject.sizeX * self.runObject.sizeY
				self.controller.show_frame('ParamGUI', self.runObject)
			except:
				tk.messagebox.showwarning('Error with Matrix size', 'You have to provide an integer as X and Y values.', icon='warning')

class InputGUIExplore(tk.Frame):
	"""
	Input page
	"""
	def __init__(self, parent, controller, runObject):
		tk.Frame.__init__(self, parent)
		self.controller = controller
		self.runObject = runObject
		self.onlyDigit = (self.register(self._validate))
		self.labelInput = tk.Label(self, text='Select dataset(s) :', font= "Arial 14")
		self.labelInput.pack(side=TOP)
		self.frameInputFile = tk.Frame(self)
		self.listRawFile = []
		self.frameMasterRaw = tk.Frame(self)
		self.frameRaw = tk.Frame(self.frameMasterRaw)
		self.varRaw = IntVar()
		self.namePathRaw = StringVar()
		self.rawBut = ttk.Button(self.frameRaw, text='Raw Matrix',  width = 20, command=lambda: self.browse_button(self.namePathRaw, 'raw'))
		self.rawBut.grid(row=1, column=1)
		self.namePathNorm = StringVar()
		self.normBut = ttk.Button(self.frameRaw, text='Normalized', width = 20, command=lambda: self.browse_button(self.namePathNorm, 'norm'))
		self.normBut.grid(row=1, column=2)
		self.namePathDiff = StringVar()
		self.diffBut = ttk.Button(self.frameRaw, text='Differential Expressed', width = 20, command=lambda: self.browse_button(self.namePathDiff, 'diff'))
		self.diffBut.grid(row=1, column=3)
		self.resetBut = ttk.Button(self.frameRaw, text='Reset Selection', width = 20, command=lambda: self.resetSelection())
		self.resetBut.grid(row=2, column=2)
		self.frameRaw.grid(row=1, column=0, sticky="w")
		self.frameMasterRaw.pack(pady=5)
		self.frameInformation = tk.Frame(self)
		self.labelInformation = tk.Label(self.frameInformation, text='Matrix Informations :', font= "Arial 14")
		self.labelInformation.grid(row=0, columnspan=3)
		self.varTranspose = IntVar()
		self.checkTranspose = Checkbutton(self.frameInformation, width=2, variable=self.varTranspose)
		self.checkTranspose.deselect()
		self.checkTranspose.grid(row=1, column=0, sticky="w")
		self.labelTranspose = tk.Label(self.frameInformation, text=': Transpose matrix')
		self.labelTranspose.grid(row=1, column=1, sticky="w")
		self.varRound = IntVar()
		self.checkRound = Checkbutton(self.frameInformation, width=2, variable=self.varRound)
		self.checkRound.deselect()
		self.checkRound.grid(row=2, column=0, sticky="w")
		self.labelRound = tk.Label(self.frameInformation, text=': Round matrix')
		self.labelRound.grid(row=2, column=1, sticky="w")
		self.frameInformation.pack()
		self.frameInformationSize = tk.Frame(self)
		self.labelX = tk.Label(self.frameInformationSize, text='Size X of matrix :')
		self.labelX.grid(row=1, column=0, sticky="w")
		self.sizeXVar = IntVar()
		self.sizeX = Entry(self.frameInformationSize, textvar=self.sizeXVar, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.sizeXVar.set(32)
		self.sizeX.grid(row=1, column=1, sticky="w")
		self.labelY = tk.Label(self.frameInformationSize, text='Size Y of matrix:')
		self.labelY.grid(row=2, column=0, sticky="w")
		self.sizeYVar = IntVar()
		self.sizeY = Entry(self.frameInformationSize, textvar=self.sizeYVar, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.sizeYVar.set(32)
		self.sizeY.grid(row=2, column=1, sticky="w")
		self.labelSave = tk.Label(self.frameInformationSize, text="Save matrix (norm, diff) :")
		self.labelSave.grid(row=3, column=0, sticky="w")
		self.varSave = IntVar()
		self.saveMatrix = Checkbutton(self.frameInformationSize, variable=self.varSave)
		self.saveMatrix.deselect()
		self.saveMatrix.grid(row=3, column=1, sticky="w")
		self.frameInformationSize.pack()
		self.frameBackNext = Frame(self)
		self.buttonBack = ttk.Button(self.frameBackNext, width=25, text="Back", command=lambda: controller.show_frame('StartPageGUI', self.runObject))
		self.buttonBack.pack(side=LEFT, fill=X)
		self.buttonNext = ttk.Button(self.frameBackNext, width=25, text="Next", command=lambda: self._next())
		self.buttonNext.pack(side=RIGHT, fill=X)
		self.frameBackNext.pack(side=BOTTOM, fill=X)

	def _validate(self, P):
		if str.isdigit(P) or P == '':
			return True
		else:
			return False

	def browse_button(self, namePath, typeOfMatrix):
		listFile = []
		if self.runObject.mode == 'cluster':
			self.fileDir = askdirectory()
			namePath.set(self.fileDir)
			if self.fileDir != '':
				for file in os.listdir(self.fileDir):
					fileName, fileExtension = os.path.splitext(file)
					if fileExtension == '.tsv' or fileExtension == '.csv':
						listFile.append(os.path.join(self.fileDir, file))
		else:
			self.fileDir = askopenfilename(multiple=False, filetypes=[('Tabulation-separated values', '*.tsv')])#,
			listFile.append(self.fileDir)
		if typeOfMatrix == 'raw':
			self.runObject.raw = listFile
		elif typeOfMatrix == 'norm':
			self.runObject.norm = listFile
		elif typeOfMatrix == 'diff':		
			self.runObject.diff = listFile
		print('\n')
		print('RAW')
		print(self.runObject.raw)
		print('NORM')
		print(self.runObject.norm)
		print('DIFF')
		print(self.runObject.diff)

	def resetSelection(self):
		self.runObject.raw = ['']
		self.runObject.norm = ['']
		self.runObject.diff = ['']
		print('All selected files are reset !')

	def _next(self):
		if self.runObject.raw[0] == '':
			tk.messagebox.showwarning('No raw file(s) selected !', 'You have to select at least one raw matrix.', icon='warning')
		else:
			try:
				if self.varTranspose.get() != 0:
					self.runObject.transpose = True
				else:
					self.runObject.transpose = False
				if self.varRound.get() != 0:
					self.runObject.round = True
				else:
					self.runObject.round = False
				self.runObject.sizeX = self.sizeXVar.get()
				self.runObject.sizeY = self.sizeYVar.get()
				self.runObject.totalGexel = self.runObject.sizeX * self.runObject.sizeY
				if self.varSave.get() != 0:
					self.runObject.saveMatrix = True
				else:
					self.runObject.saveMatrix = False
				self.controller.show_frame('ParamGUI', self.runObject)
			except:
				tk.messagebox.showwarning('Error with Matrix size', 'You have to provide an integer as X and Y values.', icon='warning')

class ParamGUI(tk.Frame):
	"""
	Parameters page
	"""
	def __init__(self, parent, controller, runObject):
		tk.Frame.__init__(self, parent)
		self.controller = controller
		self.runObject = runObject
		self.onlyDigit = (self.register(self._validate))
		listThresholdValuesUp = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
		listThresholdValuesDown = [-0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4]

		self.indexRadio = IntVar()
		self.indexRadio.set(1)

		self.frameDiff = tk.Frame(self)
		self.labelDiff = tk.Label(self.frameDiff, text='Arguments for differential gene expression :', font= "Arial 14")
		self.labelDiff.grid(row=0, columnspan=3)
		self.labelUp = tk.Label(self.frameDiff, text='Threshold up regulated : ', font= "Arial 10")
		self.labelUp.grid(row=1, column=1, sticky='e')
		self.comboUp = ttk.Combobox(self.frameDiff,
				values=listThresholdValuesUp, state='readonly', height=3, width=3)
		self.comboUp.grid(row=1, column=2 , sticky='w')
		self.comboUp.current(1)
		self.labelDown = tk.Label(self.frameDiff, text='Threshold down regulated : ', font= "Arial 10")
		self.labelDown.grid(row=2, column=1, sticky='e')
		self.comboDown = ttk.Combobox(self.frameDiff,
				values=listThresholdValuesDown, state='readonly', height=3, width=3)
		self.comboDown.grid(row=2, column=2 , sticky='w')
		self.comboDown.current(1)
		self.frameDiff.pack(pady=0.5)

		self.framePattern = tk.Frame(self)
		self.sizeMinGexelVar = IntVar()
		self.labelPattern = tk.Label(self.framePattern, text='Patterns detection :', font= "Arial 14")
		self.labelPattern.grid(row=0, columnspan=3)
		self.labelMinGexel = tk.Label(self.framePattern, text='Minimum n° of contiguous gexels :')
		self.labelMinGexel.grid(row=1, column=0, sticky="w")
		self.minGexel = Entry(self.framePattern, textvar=self.sizeMinGexelVar, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.sizeMinGexelVar.set(10)
		self.minGexel.grid(row=1, column=1, sticky="w")
		self.framePattern.pack(pady=0.5)

		self.frameMethods = tk.Frame(self)
		self.varTanimoto = IntVar()
		self.varDice = IntVar()
		self.labelMethods = tk.Label(self.frameMethods, text='Similarity Methods :', font= "Arial 14")
		self.labelMethods.grid(row=0, columnspan=3)
		self.radioTanimoto = tk.Radiobutton(self.frameMethods, variable=self.varTanimoto.get(), value='radioButton2',
			command = lambda: self.indexRadio.set(1))
		self.radioTanimoto.select()
		self.radioTanimoto.grid(row=1, column=0, sticky="e")

		self.labelTanimoto = tk.Label(self.frameMethods, text=' Tanimoto', font= "Arial 10")
		self.labelTanimoto.grid(row=1, column=1, sticky="w")

		self.radioDice = tk.Radiobutton(self.frameMethods, variable=self.varDice.get(), value='radioButton3',
			command = lambda: self.indexRadio.set(2))
		self.radioDice.deselect()
		self.radioDice.grid(row=2, column=0, sticky="e")
		self.labelDice = tk.Label(self.frameMethods, text=' Dice', font= "Arial 10")
		self.labelDice.grid(row=2, column=1, sticky="w")
		self.frameMethods.pack(pady=0.5)

		self.frameCPU = tk.Frame(self)
		self.numberCPU = IntVar()
		self.labelCPU = tk.Label(self.frameCPU, text='Multiprocessing :', font= "Arial 14")
		self.labelCPU.grid(row=0, columnspan=3)
		
		self.labelNumberCPU = tk.Label(self.frameCPU,
			text=f'Number of thread : ')
		self.labelNumberCPU.grid(row=1, column=0, sticky="e")
		self.numberCPU.set(1)
		self.cpuEntry = Entry(self.frameCPU, textvar=self.numberCPU, width=3, validate='all', validatecommand=(self.onlyDigit, '%P'))
		self.cpuEntry.grid(row=1, column=1, sticky="e")
		self.labelMaxCPU = tk.Label(self.frameCPU, text=f'max ({multiprocessing.cpu_count()}).')
		self.labelMaxCPU.grid(row=1, column=2)
		self.frameCPU.pack(pady=0.5)

		self.frameBackNext = Frame(self)
		self.buttonBack = ttk.Button(self.frameBackNext, text="Back", width=25, command=lambda: self._back())
		self.buttonBack.pack(side=LEFT, fill=X)
		self.buttonNext = ttk.Button(self.frameBackNext, text="Next", width=25, command=lambda: self._next())
		self.buttonNext.pack(side=RIGHT, fill=X)
		self.frameBackNext.pack(side=BOTTOM, fill=X)

	def _validate(self, P):
		if str.isdigit(P) or P == '':
			return True
		else:
			return False

	def _back(self):
		if self.runObject.mode == 'explore':
			self.controller.show_frame('InputGUIExplore', self.runObject)
		elif self.runObject.mode == 'cluster':
			self.controller.show_frame('InputGUIBatch', self.runObject)

	def _next(self):
		try:
			if self.numberCPU.get() > multiprocessing.cpu_count():
				tk.messagebox.showwarning('Error with Threads', f'The number of threads is higher than those available ({multiprocessing.cpu_count()}).', icon='warning')
			else:
				if self.indexRadio.get() == 1:
					self.runObject.method = 'Tanimoto'
				elif self.indexRadio.get() == 2:
					self.runObject.method = 'Dice'
				self.runObject.numberCPU = self.numberCPU.get()
				self.runObject.minPattern = self.sizeMinGexelVar.get()
				self.runObject.indexGetParameters = 1
				self.runObject.upDiffThreshold = float(self.comboUp.get())
				self.runObject.downDiffThreshold = float(self.comboDown.get())
				print('==============================')
				print('Selected Parameters')
				print('==============================')
				print(f'sizeX : {self.runObject.sizeX}')
				print(f'sizeY : {self.runObject.sizeY}')		
				print(f'mode : {self.runObject.mode}')
				print(f'raw : {self.runObject.raw}')
				print(f'norm : {self.runObject.norm}')
				print(f'diff : {self.runObject.diff}')
				print(f'transpose : {self.runObject.transpose}')
				print(f'round : {self.runObject.round}')
				print(f'method : {self.runObject.method}')
				print(f'indexGetParameters : {self.runObject.indexGetParameters}')
				print(f'size for pattern : {self.runObject.minPattern}')
				print(f'dicSample : {self.runObject.dicSample}')
				print(f'up : {self.runObject.upDiffThreshold}')
				print(f'down : {self.runObject.downDiffThreshold}')
				print(f'Save : {self.runObject.saveMatrix}')
				print('==============================\n')
				self.controller.clusterMatrix(self.runObject)
		except:
			tk.messagebox.showwarning('Error with Values', 'You have to provide a value for each field.', icon='warning')

class AutocompleteEntry(Entry):
	"""
	Autocomplete list with genes names
	AutocompleteEntry(list_input_gene_name, frame_Entry, frame_ListBox)
	"""
	def __init__(self, list_input_gene_name, frame_Entry, frame_ListBox):
		Entry.__init__(self, frame_Entry)
		self.list_input_gene_name = list_input_gene_name		
		self.var = self["textvariable"]
		self.frame_ListBox = frame_ListBox
		if self.var == '':
			self.var = self["textvariable"] = StringVar()
		self.var.trace('w', self.changed)
		self.lb_up = False

	def setEntry(self, stringToAdd):
		self.var.set(stringToAdd)
		if self.lb_up:
			self.lb.destroy()
			self.lb_up = False

	def erase(self):
		self.var.set('')

	def changed(self, name, index, mode):  
		if self.var.get() == '':
			self.lb.destroy()
			self.lb_up = False
		else:
			words = self.comparison()
			if words:			
				if not self.lb_up:
					self.lb = Listbox(self.frame_ListBox, height=5)
					self.lb.bind("<Double-Button-1>", self.selection)
					self.lb.bind("<Right>", self.selection)
					self.lb.bind("<Down>", self.down)
					self.lb.bind("<Return>", self.selection)
					self.lb.pack(side=TOP)
					self.lb_up = True
				self.lb.delete(0, END)
				for w in words:
					self.lb.insert(END,w)
			else:
				if self.lb_up:
					self.lb.destroy()
					self.lb_up = False
		
	def selection(self, event):
		if self.lb_up:
			self.var.set(self.lb.get(ACTIVE))
			self.lb.destroy()
			self.lb_up = False
			self.icursor(END)

	def goToList(self, event):
		if self.lb_up:
			self.lb.curselection()[0]

	def selectAll(self, event):
		self.var.select_range(0, len(self.var))
		self.var.icursor(5)

	def up(self, event):
		if self.lb_up:
			if self.lb.curselection() == ():
				index = '0'
			else:
				index = self.lb.curselection()[0]
			if index != '0':				
				self.lb.selection_clear(first=index)
				index = str(int(index)-1)				
				self.lb.selection_set(first=index)
				self.lb.activate(index) 

	def down(self, event):
		if self.lb_up:
			if self.lb.curselection() == ():
				index = '0'
			else:
				index = self.lb.curselection()[0]
			if index != END:						
				self.lb.selection_clear(first=index)
				index = str(int(index)+1)		
				self.lb.selection_set(first=index)
				self.lb.activate(index) 

	def comparison(self):
		pattern = re.compile(f'{self.var.get().upper()}.*')
		return [w for w in self.list_input_gene_name if re.match(pattern, w)]

class GoGUI(object):
	"""
	Gene Ontology GUI
	"""
	def __init__(self, master, listGeneQuery, listGOBank, title):
		self.master = master
		self.dic_Com_ListGene = {}
		sb.set(font_scale=0.75)
		self.listGeneQuery = listGeneQuery
		self.listGOBank = listGOBank
		self.titleGO = title
		if 'All - Communities' not in self.titleGO:
			self.master.minsize(width=600, height=500)
		else:
			self.master.minsize(width=600, height=100)
		
		if 'All - Communities' not in self.titleGO:
			self.frameGOWindow = tk.Frame(self.master)
			self.listBoxMultiple  = ttk.Treeview(self.frameGOWindow, columns=[0, 'p-values', 'Common genes'],
				displaycolumns=[0, 'p-values', 'Common genes'], show="headings", height=15)
			self.listBoxMultiple.bind('<Double-Button>', self.clipboard)
			self.listBoxMultiple.heading(0, text='GO terms')
			self.listBoxMultiple.heading('p-values', text='P-Values')
			self.listBoxMultiple.heading('Common genes', text='Genes')
			self.vsb = ttk.Scrollbar(self.frameGOWindow, orient="vertical",
				command=self.listBoxMultiple.yview)
			self.hsb = ttk.Scrollbar(self.frameGOWindow, orient="horizontal",
				command=self.listBoxMultiple.xview)
			self.listBoxMultiple.configure(yscrollcommand=self.vsb.set, xscrollcommand=self.hsb.set)
			self.listBoxMultiple.grid(column=0, row=0, sticky='nsew', in_=self.frameGOWindow)
			self.vsb.grid(column=1, row=0, sticky='ns', in_=self.frameGOWindow)
			self.hsb.grid(column=0, row=1, sticky='ew', in_=self.frameGOWindow)
			self.frameGOWindow.grid_columnconfigure(0, weight=1)
			self.frameGOWindow.grid_rowconfigure(0, weight=1)
			self.frameGOWindow.pack(fill=X, side=TOP)
		self.frameGOText = tk.Frame(self.master)
		if 'All - Communities' not in self.titleGO:
			self.textOverview = Label(self.frameGOText, text='Select a GO terms database and press \'Run\'.\nGenes can be copied to clipboard by double-click on the table.')
		else:
			self.textOverview = Label(self.frameGOText, text='Select a GO terms database and press \'Run\'.')
		self.textOverview.pack()
		self.frameGOText.pack()
		self.frameGOButton = tk.Frame(self.master)
		self.frameGOButtonTop = tk.Frame(self.frameGOButton)
		self.comboGOBank = ttk.Combobox(self.frameGOButtonTop,
				values=self.listGOBank, state='readonly', height=5)
		self.comboGOBank.pack(side=LEFT)
		self.comboGOBank.current(0)
		self.runButton = ttk.Button(self.frameGOButtonTop, text="Run", command= lambda : self.analysisGO())
		self.runButton.pack(side=LEFT)
		if 'All - Communities' not in self.titleGO:
			self.saveGOButton = ttk.Button(self.frameGOButtonTop, text="Save", command= lambda : self.saveGO())
			self.saveGOButton.pack(side=RIGHT)
			self.saveGOButton.config(state = 'disabled')
		self.frameGOButtonTop.pack(side=TOP)
		self.frameGOButtonBot = tk.Frame(self.frameGOButton)
		if 'All - Communities' not in self.titleGO:
			self.barplotButton = ttk.Button(self.frameGOButtonBot, text="Barplot",
				command= lambda : Barplot(tk.Toplevel(self.master), self.table, self.titleGO, self.comboGOBank.get()))
			self.barplotButton.pack(side=LEFT)
			self.barplotButton.config(state = "disabled")

			self.heatmapButton = ttk.Button(self.frameGOButtonBot, text="Heatmap",
				command= lambda : HeatmapGO(tk.Toplevel(self.master), self.listGeneQuery, self.table, self.titleGO, self.comboGOBank.get()))
			self.heatmapButton.pack(side=LEFT)
			self.heatmapButton.config(state = "disabled")
		self.quitGOButton = ttk.Button(self.frameGOButtonBot, text="Quit", command= lambda : self.master.destroy())
		self.quitGOButton.pack(side=RIGHT)
		self.frameGOButtonBot.pack(side=BOTTOM)
		self.frameGOButton.pack(fill=X, pady=10, side=BOTTOM)

	def clipboard(self, event):
		self.master.clipboard_clear()
		self.master.clipboard_append(self.listBoxMultiple.item(self.listBoxMultiple.focus(), 'values')[2])

	def fillGOListBox(self):
		for i in range(len(self.table)):
			if self.table['p-values'][i]!= 1:
				if len(self.table['Common genes'][i]) != 0:
					self.listBoxMultiple.insert('', 'end', values=[self.table.index[i],
						self.table['p-values'][i], self.table['Common genes'][i]])

	def fisher(self, input, disease, nb_totgen):
		input = set(input)
		disease = set(disease)
		inputxdisease = len(disease & input)
		inputxtot = len(input) - inputxdisease
		diseasewoinput = len(disease) - inputxdisease
		totwoinput = nb_totgen - len(disease) - inputxtot
		odd, pvalues = stats.fisher_exact([[inputxdisease, diseasewoinput], [inputxtot, totwoinput]])
		return pvalues

	def correction_Benjamini(self, pvalues):
		print('Ca tourne')
		pvalues = np.array(pvalues)
		n = float(pvalues.shape[0])
		print(n)
		new_pvalues = np.empty(int(n))
		values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
		values.sort()
		values.reverse()
		new_values = []
		for i, vals in enumerate(values):
			rank = n - i
			pvalue, index = vals
			new_values.append((n / rank) * pvalue)
		for i in range(0, int(n) - 1):
			if new_values[i] < new_values[i + 1]:
				new_values[i + 1] = new_values[i]
		for i, vals in enumerate(values):
			pvalue, index = vals
			new_pvalues[index] = new_values[i]
		return new_pvalues

	def inter(self, liste1, liste2):
		liste1 = set(liste1)
		liste2 = set(liste2)
		return list(liste2 & liste1)

	def analysisGO(self):
		"""
		Perform gene ontology enrichment analysis
		"""
		if 'All - Communities' in self.titleGO:
			list_GO = []
			dic = {}
			for i in self.listGeneQuery.keys():
				if self.listGeneQuery[i] not in self.dic_Com_ListGene.keys():
					self.dic_Com_ListGene[self.listGeneQuery[i]] = [i.split(' ')[0]]
				else:
					self.dic_Com_ListGene[self.listGeneQuery[i]].append(i.split(' ')[0])
			for i in self.dic_Com_ListGene.keys():
				self.dfBank = pd.read_csv(os.path.join('GO_DB', f'{self.comboGOBank.get()}.tsv'), sep='\t', index_col=[0])
				geneGOList = self.dfBank.iloc[:, 0]
				geneGOList = geneGOList.apply(lambda x: x.replace('[\'', '').replace('\']', '').split('\', \''))
				self.dfBank['Gene Signature'] = geneGOList
				tot_genes = geneGOList.tolist()
				tot_genes = [item for sublist in tot_genes for item in sublist]
				tot_genes = set(tot_genes)
				nb_totgen = len(tot_genes)
				similarity = geneGOList.apply(lambda x: self.fisher(self.dic_Com_ListGene[i], x, nb_totgen))
				print(similarity)
				similarity = similarity.sort_values()
				similarity.name = 'p-values'
				sub_similarity = similarity.index.values
				dise_comm = geneGOList.loc[sub_similarity]
				dise_comm = dise_comm.apply(lambda x: self.inter(x, self.dic_Com_ListGene[i]))
				dise_comm = dise_comm.dropna()
				dise_comm.name = 'Common genes'
				self.table = pd.concat([similarity, dise_comm], axis=1)
				self.table = self.table.dropna()
				print(f'### TOP5 GO TERMS - Communities {i}')
				print(self.table.index[:5].tolist())
				dic[i] = self.table
				list_GO.extend(self.table.index[:5].tolist())
			print(list(set(list_GO)))
			print(list_GO)
			l_list_GO =  list(set(list_GO))
			l = []
			minV = 1
			maxV = 0
			for i in dic.keys():
				ll = []
				for u in list(set(list_GO)):
					ll.append(dic[i].loc[str(u), 'p-values'])
				if minV > min(ll):
					minV = min(ll)
				if maxV < max(ll):
					maxV = max(ll) 
				l.append(ll)
			ar = np.array(l)
			dfT = pd.DataFrame(ar, index=dic.keys(), columns=list(set(list_GO)))
			HeatmapGOAll(tk.Toplevel(self.master), self.listGeneQuery, dfT, f'{self.titleGO} - TOP5', self.comboGOBank.get(), minV, maxV, True)
		else:
			for item in self.listBoxMultiple.get_children():
				self.listBoxMultiple.delete(item)
			self.dfBank = pd.read_csv(os.path.join('GO_DB', f'{self.comboGOBank.get()}.tsv'), sep='\t', index_col=[0])
			geneGOList = self.dfBank.iloc[:, 0]
			geneGOList = geneGOList.apply(lambda x: x.replace('[\'', '').replace('\']', '').split('\', \''))
			self.dfBank['Gene Signature'] = geneGOList
			tot_genes = geneGOList.tolist()
			tot_genes = [item for sublist in tot_genes for item in sublist]
			tot_genes = set(tot_genes)
			nb_totgen = len(tot_genes)
			self.listGeneQuery = [x.split(' ')[0].upper() for x in self.listGeneQuery]
			similarity = geneGOList.apply(lambda x: self.fisher(self.listGeneQuery, x, nb_totgen))
			similarity = similarity.sort_values()
			similarity.name = 'p-values'
			sub_similarity = similarity.index.values
			dise_comm = geneGOList.loc[sub_similarity]
			dise_comm = dise_comm.apply(lambda x: self.inter(x, self.listGeneQuery))
			dise_comm = dise_comm.dropna()
			dise_comm.name = 'Common genes'
			self.table = pd.concat([similarity, dise_comm], axis=1)
			self.table = self.table.dropna()
			print(self.table)
			self.fillGOListBox()
			try:
				self.saveGOButton.config(state = 'normal')
				self.barplotButton.config(state = 'normal')
				self.heatmapButton.config(state = 'normal')
				self.correctedscoreButton(state = 'normal')
			except:
				pass

	def saveGO(self):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			self.table.to_csv(filename, sep='\t')

class HeatmapGOAll(object):
	"""
	Heatmap for communities
	"""
	def __init__(self, master, listGeneQuery, table, titleGO, comboGOBank, minV, maxV, annotation):
		self.master = master
		self.frameTOP = tk.Frame(self.master)
		self.listGeneQuery = listGeneQuery
		self.fig = Figure()
		ax = self.fig.add_subplot(111)
		if annotation == True:
			plot = sb.heatmap(table.transpose(), norm=LogNorm(), linewidths=0.3, cmap=sb.light_palette(color='green', as_cmap = True, reverse = True),
				linecolor='black', ax=ax, xticklabels=True, yticklabels=True, vmin=minV, vmax=0.10,
				annot=True, annot_kws={'size': 7.5})
		else:
			plot = sb.heatmap(table.transpose(), norm=LogNorm(), linewidths=0.3, cmap=sb.light_palette(color='green', as_cmap = True, reverse = True),
				linecolor='black', ax=ax, xticklabels=True, yticklabels=True, vmin=minV, vmax=0.10,
				annot=False, annot_kws={'size': 7.5})
		plot.xaxis.set_ticks_position('top')
		plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
		plot.set_yticklabels(plot.get_yticklabels(), rotation=0)
		plot.set(xlabel='Gene Ontology Terms', ylabel='TOP GO terms per community')
		plot.set_title(f'Heatmap GO Terms - {titleGO} - {comboGOBank}')
		canvas = FigureCanvasTkAgg(self.fig, master=self.frameTOP)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(row=0, column=0)
		self.frameTOP.pack()
		self.framBOT = tk.Frame(self.master)
		ttk.Button(self.framBOT,text="Save Heatmap",command = lambda : self.saveImage()).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Save Matrix",command = lambda : self.saveMatrix(table.transpose())).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Quit",command = lambda : self.master.destroy()).pack(side=LEFT)
		self.framBOT.pack(pady=10, side=BOTTOM)

	def saveImage(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			self.fig.savefig(filename, format='png', bbox_inches='tight')

	def saveMatrix(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.iloc[::-1].to_csv(filename, sep='\t')

class HeatmapGO(object):
	"""
	Heatmap for gene ontology enrichment analysis
	"""
	def __init__(self, master, listGeneQuery, table, titleGO, comboGOBank):
		self.master = master
		self.frameTOP = tk.Frame(self.master)
		self.listGeneQuery = listGeneQuery
		self.master.title(f'Top40 genes - Heatmap GO Terms - {titleGO} - {comboGOBank}')
		self.table = table.iloc[:15]
		matric_hit = self.table.dropna().iloc[:, 1]
		self.dico_hit = matric_hit.to_dict()
		matric_hit = pd.DataFrame([self.listGeneQuery for i in range(len(self.table))],
			index=self.table.index.tolist(), columns=self.listGeneQuery)
		matric_hit = matric_hit.apply(lambda x: self.intersec_serie(x), axis=1)
		self.fig = Figure()
		ax = self.fig.add_subplot(111)
		matric_hit = matric_hit.transpose()
		matric_hit['sum'] = matric_hit.apply(lambda x : sum(x), axis=1)
		matric_hit = matric_hit.sort_values(by ='sum', ascending=True)
		df2 = matric_hit.drop(['sum'], axis = 1)
		cmap = LinearSegmentedColormap.from_list('Custom', ['#D3D3D3', '#319819'], 2)
		df3 = None
		print(f'count : {df2.shape[0]}')
		if df2.shape[0] >= 40:
			df3 = df2[-40:]
		else:
			df3 = df2
		plot = sb.heatmap(df3, cmap=cmap, linewidths=0.3,
			linecolor='black', ax=ax, xticklabels=True, yticklabels=True)
		plot.set_ylim(0 , len(df3))
		plot.xaxis.set_ticks_position('top')
		plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
		colorbar = ax.collections[0].colorbar
		colorbar.set_ticks([0.25, 0.75])
		colorbar.set_ticklabels(['Absence', 'Presence'])
		plot.set_yticklabels(plot.get_yticklabels(), rotation=35)
		plot.set(xlabel='Gene Ontology Terms', ylabel='TOP40 - Genes')
		plot.set_title(f'Heatmap GO Terms - {titleGO} - {comboGOBank}')
		canvas = FigureCanvasTkAgg(self.fig, master=self.frameTOP)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(row=0, column=0)
		self.frameTOP.pack()
		self.framBOT = tk.Frame(self.master)
		ttk.Button(self.framBOT,text="Save Heatmap",command = lambda : self.saveImage()).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Save Matrix",command = lambda : self.saveMatrix(df2)).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Quit",command = lambda : self.quit()).pack(side=LEFT)
		self.framBOT.pack(pady=10, side=BOTTOM)

	def saveImage(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			self.fig.savefig(filename, format='png', bbox_inches='tight')

	def saveMatrix(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.iloc[::-1].to_csv(filename, sep='\t')

	def quit(self):
		self.master.destroy()

	def inter(self, liste1, liste2):
		liste1 = set(liste1)
		liste2 = set(liste2)
		return list(liste2 & liste1)

	def intersec_serie(self, serie):
		name = serie.name
		serie = serie.apply(lambda x: self.intersect(x, name))
		return serie

	def intersect(self, x, name):
		if x in self.dico_hit[name]:
			return 1
		else:
			return 0

class Barplot(object):
	"""
	Barplot for gene ontology enrichment analysis
	"""
	def __init__(self, master, df, titleGO, comboGOBank):
		self.master = master
		self.master.title(f'Barplot GO Terms - {titleGO} - {comboGOBank}')
		df2 = df.iloc[:15].drop(['Common genes'], axis = 1)
		self.frameInputFile = tk.Frame(self.master)
		self.fig = Figure()
		ax = self.fig.add_subplot(111)
		plot = sb.barplot(x=1/df2['p-values'] , y=[i for i in range(1, len(df2.index)+1)],
			data=df2, orient='h', ax=ax)
		plot.set(xscale="log")
		plot.set(xlabel='1/P-values (log10)', ylabel='TOP 15 Gene Ontology Terms')
		plot.set_title(f'Barplot GO Terms - {titleGO} - {comboGOBank}')
		ind = 0
		_x = ax.patches[len(df2)-1].get_width()
		for p in ax.patches:
			_y = p.get_y()+ 0.6
			value = int(p.get_width())
			ax.text(_x, _y, df2.index[ind], ha="left")
			ind += 1
		canvas = FigureCanvasTkAgg(self.fig, master=self.frameInputFile)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(row=0, column=0)
		self.frameInputFile.pack()
		self.framBOT = tk.Frame(self.master)
		ttk.Button(self.framBOT,text="Save Barplot",command = lambda : self.saveImage()).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Save Matrix",command = lambda : self.saveMatrix(df.drop(['Common genes'], axis = 1))).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Quit",command = lambda : self.quit()).pack(side=LEFT)
		self.framBOT.pack(pady=10)
	
	def saveImage(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			self.fig.savefig(filename, format='png', bbox_inches='tight')

	def saveMatrix(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.to_csv(filename, sep='\t')

	def quit(self):
		self.master.destroy()

class mainGUI(object):
	"""
	Main GUI panel
	"""
	def __init__(self, master, runObject):
		self.master = master
		self.master.protocol('WM_DELETE_WINDOW', self.close_windows)
		self.runObject = runObject
		self.textStringVar = StringVar()
		self.list_GO_Bank = self.loadGOBank()
		self.frame = tk.Frame(self.master)
		self.frameText = tk.Frame(self.frame,  height=10, width=10)
		self.master.title('MULTILAYER')
		if self.runObject.mode == 'explore':
			self.master.minsize(260, 405)
			self.textOverview = Label(self.frameText, text='Informations :')
			self.textOverview.pack(side=TOP)
			self.textStatPanel = Text(self.frameText, height=9, width=30)
			self.nameSample = list(self.runObject.dicSample.keys())[0]
			self.textString = f'Sample ID:\n\t{self.nameSample}\n\
				Matrix Size:\n\t{self.runObject.sizeX}x{self.runObject.sizeY} ({self.runObject.totalGexel} gexels)\n\
				Gexels Containing Signal:\n\t{len(self.runObject.dicSample[self.nameSample].raw.columns.values)} ({round((len(self.runObject.dicSample[self.nameSample].diff.columns.values)/self.runObject.totalGexel)*100, 2)}%)\n\
				Minimun n° of Contiguous Gexels\
				per Pattern:\n\t{self.runObject.minPattern}\n\
				Similarity Method:\n\t{self.runObject.method}\n\
				Differential Genes\nexpression:\nThreshold up :\t{self.runObject.upDiffThreshold}\nThreshold down :\t{self.runObject.downDiffThreshold}'
			self.textStatPanel.insert(END, f'{self.textString}\n')
			self.textStatPanel.config(state = DISABLED)
			self.textStatPanel.pack(side = LEFT, fill = X)
			self.scrollbarText = Scrollbar(self.frameText, orient=VERTICAL)
			self.scrollbarText.pack(fill = Y, side = RIGHT)
			self.scrollbarText.config(command = self.textStatPanel.yview)
			self.textStatPanel.config(yscrollcommand = self.scrollbarText.set)
			self.textStatPanel.pack(fill=X)
		elif self.runObject.mode == 'cluster':
			self.master.minsize(260, 450)
			self.displaySampleSelected = StringVar()
			self.textOverview = Label(self.frameText, text='Select a sample :')
			self.textOverview.pack(side=TOP)
			self.listCluster = Listbox(self.frameText, height=9, width=43)
			self.fillListbox(self.listCluster, self.runObject.dicSample.keys())
			self.listCluster.bind('<Double-Button>', self.selection_listbox)
			self.listCluster.bind('<Return>', self.selection_listbox)
			self.listCluster.bind('<Button>', self.selection_listbox)
			self.listCluster.pack(fill=X)
		self.frameText.pack(fill = X)
		if self.runObject.mode == 'cluster':
			self.frameDisplaySample = tk.Frame(self.frame)
			self.labelDisplaySampleText = Label(self.frameDisplaySample, text='Sample : ')
			self.labelDisplaySample = Label(self.frameDisplaySample, textvariable=self.displaySampleSelected)
			self.labelDisplaySampleText.pack(side=LEFT)
			self.labelDisplaySample.pack(side=LEFT, fill = X)
			self.frameDisplaySample.pack(fill=X)
		self.frameDisplayCoor = tk.Frame(self.frame)
		self.countDisplayStringCoor = StringVar()
		self.labelDisplayCoorText = Label(self.frameDisplayCoor, text='Coordinate : ')
		self.labelDisplayCoor = Label(self.frameDisplayCoor, textvariable=self.countDisplayStringCoor)
		self.labelDisplayCoorText.pack(side=LEFT)
		self.labelDisplayCoor.pack(side=LEFT, fill = X)
		self.frameDisplayCoor.pack(fill=X)
		self.frameCountDisplay = tk.Frame(self.frame)
		self.countDisplayString = StringVar()
		self.labelDisplayCountText = Label(self.frameCountDisplay, text='Number of count(s) : ')
		self.labelDisplayCount = Label(self.frameCountDisplay, textvariable=self.countDisplayString)
		self.labelDisplayCountText.pack(side=LEFT)
		self.labelDisplayCount.pack(side=LEFT, fill = X)
		self.frameCountDisplay.pack(fill=X)
		self.frameCountDisplayDiff = tk.Frame(self.frame)
		self.countDisplayStringDiff = StringVar()
		self.labelDisplayCountTextDiff = Label(self.frameCountDisplayDiff, text='Level of expression : ')
		self.labelDisplayCountDiff = Label(self.frameCountDisplayDiff, textvariable=self.countDisplayStringDiff)
		self.labelDisplayCountTextDiff.pack(side=LEFT)
		self.labelDisplayCountDiff.pack(side=LEFT, fill = X)
		self.frameCountDisplayDiff.pack(fill=X)
		if self.runObject.mode == 'cluster':
			self.countMatrixButton = ttk.Button(self.frame, text = 'Raw Counts', width = 25, 
				command = lambda : self.matrixGUIwrapper(self.runObject.dicSample[self.displaySampleSelected.get()].raw,
				f'{self.displaySampleSelected.get()} - Counts Matrix', 'rainbow', self.countDisplayString, self.displaySampleSelected.get(), []))
			self.normalizedCountMatrixButton = ttk.Button(self.frame, text = 'Normalized Counts', width = 25, 
				command = lambda : self.matrixGUIwrapper(self.runObject.dicSample[self.displaySampleSelected.get()].norm,
				f'{self.displaySampleSelected.get()} - Normalized Counts Matrix', 'rainbow', self.countDisplayString, self.displaySampleSelected.get(), []))
			self.diffGeneButton = ttk.Button(self.frame, text = 'Gene Expression Matrix', width = 25, 
				command = lambda : self.matrixWindow(f'{self.displaySampleSelected.get()} - Gene Expression Levels', self.runObject.dicSample[self.displaySampleSelected.get()].listGexelDiff,
				self.runObject.dicSample[self.displaySampleSelected.get()].dicSampleDiff, 4, -4, 'diffGeneGradient', self.countDisplayStringDiff, self.displaySampleSelected.get(), []))
			self.patternButton = ttk.Button(self.frame, text = 'Genes Co-Expression Patterns', width = 25, 
				command = lambda : self.matrixWindow(f'{self.displaySampleSelected.get()} - Genes Patterns Matrix', self.runObject.dicSample[self.displaySampleSelected.get()].listGexelDiff,
				self.runObject.dicSample[self.displaySampleSelected.get()].dicSampleDiff, 0, 0, 'patternGradient', self.countDisplayStringDiff, self.displaySampleSelected.get(), []))
			self.communitiesButton = ttk.Button(self.frame, text = 'Gexel Communities',
				command = lambda : self.wrapperCommunities(self.displaySampleSelected.get()))
			self.classButton = ttk.Button(self.frame, text = 'Comparison datasets',
				command = lambda : self.classAnalysis())
		else:
			self.countMatrixButton = ttk.Button(self.frame, text = 'Raw Counts', width = 25, 
				command = lambda : self.matrixWindow('Counts Matrix', self.runObject.dicSample[self.nameSample].dicGexelRaw,
				self.runObject.dicSample[self.nameSample].dicSampleRaw, self.runObject.dicSample[self.nameSample].maxCountRaw,
				self.runObject.dicSample[self.nameSample].minCountRaw, 'rainbow', self.countDisplayString, self.nameSample, []))
			self.normalizedCountMatrixButton = ttk.Button(self.frame, text = 'Normalized Counts', width = 25, 
				command = lambda : self.matrixWindow('Normalized Counts Matrix', self.runObject.dicSample[self.nameSample].dicGexelNorm,
				self.runObject.dicSample[self.nameSample].dicSampleNorm, self.runObject.dicSample[self.nameSample].maxCountNorm,
				self.runObject.dicSample[self.nameSample].minCountNorm, 'rainbow', self.countDisplayString, self.nameSample, []))
			if self.runObject.dicSample[self.nameSample].dicGexelNorm == None:
				self.normalizedCountMatrixButton.config(state = "disabled")
			self.diffGeneButton = ttk.Button(self.frame, text = 'Gene Expression Levels', width = 25, 
				command = lambda : self.matrixWindow('Gene Expressed Matrix', self.runObject.dicSample[self.nameSample].listGexelDiff,
				self.runObject.dicSample[self.nameSample].dicSampleDiff, 4, -4, 'diffGeneGradient', self.countDisplayStringDiff, self.nameSample, []))
			self.patternButton = ttk.Button(self.frame, text = 'Genes Co-Expression Patterns', width = 25, 
				command = lambda : self.wrapperPatterns())
			self.communitiesButton = ttk.Button(self.frame, text = 'Gexel Communities',
				command = lambda : self.wrapperCommunities(self.nameSample))
		self.countMatrixButton.pack(fill = X)
		self.normalizedCountMatrixButton.pack(fill = X)
		self.diffGeneButton.pack(fill = X)
		self.patternButton.pack(fill = X)
		self.communitiesButton.pack(fill=X)
		if self.runObject.mode == 'cluster':
			self.classButton.pack(fill=X)
			self.countMatrixButton.config(state = 'disabled')
			self.normalizedCountMatrixButton.config(state = 'disabled')
			self.diffGeneButton.config(state = 'disabled')
			self.patternButton.config(state = 'disabled')	
			self.communitiesButton.config(state = 'disabled')
		self.frameQuit = Frame(self.frame)
		self.buttonHome = ttk.Button(self.frameQuit, text = 'Home', command = self.goHome)
		self.buttonHome.pack(side = TOP, fill=X)
		self.buttonQuit = ttk.Button(self.frameQuit, text = 'Quit', command = self.close_windows)
		self.buttonQuit.pack(side = BOTTOM, fill=X)
		self.frameQuit.pack(fill = X)
		self.frame.pack(fill=X, side=LEFT)

	def goHome(self):
		"""
		Back to home
		"""
		closeGUI = tk.messagebox.askquestion('Go Home', 'All loaded data will be lost.', icon='warning')
		if closeGUI == 'yes':
			self.master.destroy()
			delattr(self.runObject, 'sizeX')
			delattr(self.runObject, 'sizeY')
			delattr(self.runObject, 'totalGexel')
			delattr(self.runObject, 'mode')
			delattr(self.runObject, 'raw')
			delattr(self.runObject, 'norm')
			delattr(self.runObject, 'diff')
			delattr(self.runObject, 'transpose')
			delattr(self.runObject, 'round')
			delattr(self.runObject, 'method')
			delattr(self.runObject, 'indexGetParameters')
			delattr(self.runObject, 'dicSample')
			delattr(self.runObject, 'numberCPU')
			delattr(self.runObject, 'indexRun')
			delattr(self.runObject, 'minPattern')
			delattr(self.runObject, 'saveMatrix')
			delattr(self.runObject, 'upDiffThreshold')
			delattr(self.runObject, 'downDiffThreshold')
			del self.runObject
			main()		
		
	def fillListbox(self, listBox, listData):
		listBox.delete(0, END)
		for i in listData:
			listBox.insert(END, i)

	def selection_listbox(self, event):
		"""
		Select sample in batch mode
		"""
		self.sampleSelected = self.listCluster.get(ACTIVE)
		self.displaySampleSelected.set(self.sampleSelected)
		print(f'You selected {self.sampleSelected}')
		self.countMatrixButton.config(state = 'normal')
		self.normalizedCountMatrixButton.config(state = 'normal')
		self.diffGeneButton.config(state = 'normal')
		self.patternButton.config(state = 'normal')
		self.communitiesButton.config(state = 'normal')

	def communitiesArgs(self, master):
		"""
		Communities arguments panel
		"""
		self.frameArg = tk.Frame(master)
		self.frameWeight = tk.Frame(self.frameArg)
		self.labelArg = tk.Label(self.frameWeight, text='Arguments for communities')
		self.labelArg.pack(side=TOP)
		self.varWeight = IntVar()
		self.labelWeight = tk.Label(self.frameWeight, text='Weight : ')
		self.labelWeight.pack(side=LEFT)
		self.checkWeight = Checkbutton(self.frameWeight, variable=self.varWeight)
		self.checkWeight.select()
		self.checkWeight.pack(side=LEFT)
		self.frameWeight.pack(side=TOP)
		self.thresholdCommunities = IntVar()
		self.frameThresoldCommunities = tk.Frame(self.frameArg)
		self.frameRandom = tk.Frame(self.frameArg)
		self.varRandom = IntVar()
		self.labelRandom = tk.Label(self.frameRandom, text='Multiple iterations (15 times) : ')
		self.labelRandom.pack(side=LEFT)
		self.checkRandom = Checkbutton(self.frameRandom, variable=self.varRandom)
		self.checkRandom.select()
		self.checkRandom.pack(side=LEFT)
		self.frameRandom.pack(side=TOP)
		self.thresoldCommunities = Label(self.frameThresoldCommunities, text='Threshold <= :')
		self.thresoldCommunities.pack(side=LEFT)
		self.entryThresholdCommunities = Entry(self.frameThresoldCommunities, textvar=self.thresholdCommunities, width=3)
		self.entryThresholdCommunities.pack(side=LEFT)
		self.thresholdCommunities.set(0)
		self.percentThresoldCommunities = Label(self.frameThresoldCommunities, text='%.')
		self.percentThresoldCommunities.pack(side=LEFT)
		self.frameThresoldCommunities.pack(side=BOTTOM)
		self.frameArg.pack(side=TOP)
		self.runcommunitiesButton = ttk.Button(master, text='Run Communities',
			command = lambda : self.louvain(self.runObject.dicSample[self.nameSample].matrixLouvain, self.varWeight.get(), self.varRandom.get(), master))
		self.runcommunitiesButton.pack()

	def louvain(self, dfLouvain, weightIndex, randomIndex, frame):
		"""
		Perform communities detection with Louvain algorithm
		"""
		dic_valCount_df = {}
		dic_valCount_index = {}
		indexRandom = 1
		if randomIndex == 1:
			indexRandom = 15
		for i in range(indexRandom):
			df = dfLouvain.copy()
			df['Community'] = df['TF']
			weightString = ''
			indexNames = df[(df['TF'] == df['TG'])].index
			df.drop(indexNames, inplace=True)
			indexNames2 = df[(df['Weight'] <= self.thresholdCommunities.get())].index
			df.drop(indexNames2, inplace=True)
			if weightIndex == 1:
				nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr='Weight')
				partition = louvain.best_partition(nw, weight='Weight')#, random_state=100)
				weightString = 'Weight'
			elif weightIndex == 0:
				nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr=None)
				partition = louvain.best_partition(nw, weight='None')#, random_state=100)
				weightString = 'without Weight'
			df['Community'] = df['Community'].map(partition)
			tmp = df['Community'].value_counts().to_string()
			if tmp not in dic_valCount_df.keys():
				dic_valCount_df[tmp] = df
				dic_valCount_index[tmp] = 1
			else:
				dic_valCount_index[tmp] += 1
		valCount_index_sort = sorted(dic_valCount_index.items(), key=lambda x: x[1], reverse=True)
		if indexRandom != 1:
			self.logFileRandomCommunities(valCount_index_sort, indexRandom)
		self.runObject.dicSample[self.nameSample].currentLouvain = dic_valCount_df[valCount_index_sort[0][0]]
		frame.destroy()
		title = f'Communities - {weightString} - {self.thresholdCommunities.get()} %'
		self.matrixWindow(title, self.runObject.dicSample[self.nameSample].listGexelDiff,
				self.runObject.dicSample[self.nameSample].dicSampleDiff, 0, 0, 'communitiesGradient', self.countDisplayStringDiff, self.nameSample, self.runObject.dicSample[self.nameSample].currentLouvain)

	def logFileRandomCommunities(self, tupleSort, total):
		"""
		Results of Louvain communities detection
		"""
		print('##########\n#  LOG   #\n##########')
		for i in range(len(tupleSort)):
			if i == 0:
				print(f'{i+1} [*] {(tupleSort[i][1]/total)*100}%')
			else:
				print(f'{i+1} [ ] {(tupleSort[i][1]/total)*100}%')
		print('\n')

	def classAnalysis(self):
		self.classWindow = tk.Toplevel(self.master)
		wrapperClass(self.classWindow, self.runObject)
		self.classWindow.title('Arguments for class')
		self.classWindow.minsize(260, 150)

	def wrapperCommunities(self, nameOfSample):
		self.nameSample = nameOfSample
		self.communitiesWindow = tk.Toplevel(self.master)
		self.communitiesArgs(self.communitiesWindow)
		self.communitiesWindow.title('Arguments for communities')
		self.communitiesWindow.minsize(260, 150)

	def wrapperPatterns(self):
		self.patternWindow = tk.Toplevel(self.master)
		self.patternsArgs(self.patternWindow)
		self.patternWindow.title('Arguments for patterns')
		self.patternWindow.minsize(260, 60)

	def patternsArgs(self, master):
		self.frameArg = tk.Frame(master)
		self.runObject.minPattern
		self.thresholdMinPattern = IntVar()
		self.labelthresholdMinPattern = Label(self.frameArg, text='Min. Gexel(s) per pattern :')
		self.labelthresholdMinPattern.pack(side=LEFT)
		self.entryMinPattern = Entry(self.frameArg, textvar=self.thresholdMinPattern, width=3)
		self.entryMinPattern.pack(side=LEFT)
		self.thresholdMinPattern.set(int(self.runObject.minPattern))
		self.frameArg.pack(side=TOP)
		self.runcommunitiesButton = ttk.Button(master, text='Run Patterns',
			command = lambda : self.wrapperMinPattern(master))
		self.runcommunitiesButton.pack()

	def wrapperMinPattern(self, master):
		if self.thresholdMinPattern.get() != self.runObject.minPattern:
			self.runObject.minPattern = self.thresholdMinPattern.get()
			print(self.runObject.minPattern)
			StartGUI.agglomerative(StartGUI, self.runObject.dicSample[self.nameSample].dicSampleDiff,
				self.runObject.dicSample[self.nameSample], self.runObject,  self.runObject.minPattern)
			StartGUI.generateMatrixSimilarity(StartGUI, self.runObject.dicSample[self.nameSample].dicSampleDiff,
				self.runObject, 'matrix_prostate/Similarity_Matrix_Sample', self.runObject.dicSample[self.nameSample].name)
		master.destroy()
		self.matrixWindow('Patterns Matrix', self.runObject.dicSample[self.nameSample].listGexelDiff,
			self.runObject.dicSample[self.nameSample].dicSampleDiff, 0, 0, 'patternGradient', self.countDisplayStringDiff, self.nameSample, [])

	def matrixWindow(self, title, gexel, geneSample, maxCount, minCount, gradient, stringDisplay, nameSample, currentLouvain):
		self.windowGUI = tk.Toplevel(self.master)
		self.title = f'All - {title}'
		matrixGUI(self.title, self.windowGUI, self.runObject, self.countDisplayStringCoor, stringDisplay, nameSample, gexel, geneSample, maxCount, minCount, gradient, self.list_GO_Bank, currentLouvain)#, dicCoordianteGexel, minGeneCountValue, maxGeneCountValue, dic_geneSample, 'rainbowGradient', self.countDisplayString, self.countDisplayStringCoor, 'menu', self.list_pattern)
		self.windowGUI.title(self.title)		

	def close_windows(self):
		closeGUI = tk.messagebox.askquestion('Quit now ?', 'Do you want to quit now ? Are you sure ?', icon='warning')
		if closeGUI == 'yes':
			self.master.destroy()

	def matrixGUIwrapper(self, sampleMatrix, title, gradient, stringDisplay, nameSample, currentLouvain):
		self.windowGUI = tk.Toplevel(self.master)
		self.title = title
		listGexelCount = []
		dicGexel = dict(sampleMatrix.apply(lambda x : StartGUI.createGexel(StartGUI, x, listGexelCount, self.runObject), axis=0))
		maxCountGexel =  max(listGexelCount)
		minCountGexel =  min(listGexelCount)
		dicGeneSample = dict(sampleMatrix.apply(lambda x : StartGUI.createGeneSample(StartGUI, x), axis=1))
		matrixGUI(self.title, self.windowGUI, self.runObject, self.countDisplayStringCoor, stringDisplay, nameSample, dicGexel, dicGeneSample, maxCountGexel, minCountGexel, gradient, [], currentLouvain)
		self.windowGUI.title(self.title)

	def loadGOBank(self):
		"""
		Load GO bank file. The Directory is 'GO_DB'.
		"""
		list_GO_Bank_temp = []
		for file in os.listdir('GO_DB'):
			if not file.startswith('.'):
				list_GO_Bank_temp.append(os.path.splitext(file)[0])
		return(sorted(list_GO_Bank_temp))

class resultGOClass(object):
	"""
	datasets communities heatmap for batch mode
	"""
	def __init__(self, master, runObject, df, dic):
		self.master = master
		self.runObject = runObject
		self.df = df
		self.widthGUI  = self.master.winfo_screenwidth()
		self.heightGUI = self.master.winfo_screenheight()
		self.dic = dic
		self.indexMiniSizeGexelX = round((self.widthGUI/2)/(len(self.df.columns.values)+4))
		self.indexMiniSizeGexelY = round((3*self.heightGUI/4)/(len(self.df.index.values)+4))
		self.textTop = Label(self.master, text='Dataset(s) communities', font= "Arial 14")
		self.textTop.grid(row=0, column=1, sticky=S, pady=0)
		self.textLeft = Label(self.master, text='C\nl\na\ns\ns\ne\ns', font= "Arial 14")
		self.textLeft.grid(row=1, column=4, sticky=E, pady=0)		
		self.can = Canvas(self.master, width=((len(self.df.columns.values)+1)*(self.indexMiniSizeGexelX+1)), 
			height = ((len(self.df.index.values)+1)*(self.indexMiniSizeGexelY+1)))
		self.can.grid(row=1, columnspan=3, sticky=NW, pady=0)
		self.frameDisplay = tk.Frame(self.master)
		self.frameCommunities = tk.Frame(self.frameDisplay)
		self.labelDisplayCommunities = Label(self.frameCommunities, text='Community : ',  font= "Arial 12")
		self.labelDisplayCommunities.pack(side=LEFT)
		self.countDisplayCommunities = StringVar()
		self.labelcountDisplayCommunities = Label(self.frameCommunities, textvariable=self.countDisplayCommunities, font= "Arial 12")
		self.labelcountDisplayCommunities.pack(side=LEFT)
		self.frameCommunities.pack(side=TOP)
		self.frameClass = tk.Frame(self.frameDisplay)
		self.labelDisplayClass = Label(self.frameDisplay, text='Class : ',  font= "Arial 12")
		self.labelDisplayClass.pack(side=LEFT)
		self.countDisplayClass = StringVar()
		self.labelcountDisplayClass = Label(self.frameClass, textvariable=self.countDisplayClass, font= "Arial 12")
		self.labelcountDisplayClass.pack(side=LEFT)
		self.frameClass.pack(side=BOTTOM)
		self.frameDisplay.grid(row=2, column=1, sticky=NW, pady=0)
		self.frameButton = tk.Frame(self.master)
		self.buttonSave = ttk.Button(self.frameButton, text="Save", command= lambda : self.wrapperHeatmapClass(master))
		self.buttonSave.pack(side=LEFT)
		self.buttonSave = ttk.Button(self.frameButton, text="GO Terms", command= lambda : self.geneOntologyClass(master, self.dic, mainGUI.loadGOBank(mainGUI), 'All'))
		self.buttonSave.pack(side=LEFT)
		self.buttonQuit = ttk.Button(self.frameButton, text="Quit", command= lambda : master.destroy())
		self.buttonQuit.pack(side=LEFT)
		self.frameButton.grid(row=3, column=1, sticky=NW)
		self.drawingHeatmap()

	def geneOntologyClass(self, master, dicQuery, listGOBank, title):
		self.master = master
		self.topLvl = tk.Toplevel(self.master)
		self.topLvl.title(f'Gene Ontology Analysis - Class')
		self.guiGOClass(self.topLvl, dicQuery, listGOBank, title)

	def guiGOClass(self, master, dicQuery, listGOBank, title):
		self.dic_Com_ListGene = {}
		sb.set(font_scale=0.75)
		self.dicQuery = dicQuery
		self.listGOBank = listGOBank
		self.titleGO = title
		master.minsize(width=600, height=100)
		self.frameGOText = tk.Frame(master)
		self.textOverview = Label(self.frameGOText, text='Select a GO terms database and press \'Run\'.')
		self.textOverview.pack()
		self.frameListLabal = tk.Frame(self.frameGOText)
		self.textTopGO = Label(self.frameListLabal, text='Top GO terms :')
		self.textTopGO.pack(side=LEFT)
		self.listCombobox = ttk.Combobox(self.frameListLabal,
				values=[5, 10, 15], state='readonly', height=3, width=3)
		self.listCombobox.current(0)
		self.listCombobox.pack(side=RIGHT)
		self.frameListLabal.pack()
		self.frameGOText.pack()
		self.frameGOButton = tk.Frame(master)
		self.frameGOButtonTop = tk.Frame(self.frameGOButton)
		self.comboGOBank = ttk.Combobox(self.frameGOButtonTop,
				values=self.listGOBank, state='readonly', height=5)
		self.comboGOBank.pack(side=LEFT)
		self.comboGOBank.current(0)
		self.runButton = ttk.Button(self.frameGOButtonTop, text="Run", command= lambda : self.analysisGO(master, int(self.listCombobox.get())))
		self.runButton.pack(side=LEFT)
		self.frameGOButtonTop.pack(side=TOP)
		self.frameGOButtonBot = tk.Frame(self.frameGOButton)
		self.quitGOButton = ttk.Button(self.frameGOButtonBot, text="Quit", command= lambda : master.destroy())
		self.quitGOButton.pack(side=RIGHT)
		self.frameGOButtonBot.pack(side=BOTTOM)
		self.frameGOButton.pack(fill=X, pady=10, side=BOTTOM)

	def wrapperHeatmapClass(self, master):
		self.master = master
		self.heatmapClassFrame = tk.Toplevel(self.master)
		self.heatmapClass(self.heatmapClassFrame)

	def heatmapClass(self, master):
		self.master = master
		self.frameTOP = tk.Frame(self.master)
		self.master.title(f'Heatmap Class')
		self.fig = Figure()
		ax = self.fig.add_subplot(111)
		cmap = LinearSegmentedColormap.from_list('Custom', ['#D3D3D3', '#319819'], 2)
		plot = sb.heatmap(self.df, cmap=cmap, linewidths=0.3,
			linecolor='black', ax=ax, xticklabels=True, yticklabels=True)
		plot.set_ylim(0 , len(self.df))
		plot.xaxis.set_ticks_position('top')
		plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
		colorbar = ax.collections[0].colorbar
		colorbar.set_ticks([0.25, 0.75])
		colorbar.set_ticklabels(['Absence', 'Presence'])
		plot.set_yticklabels(plot.get_yticklabels(), rotation=0)
		plot.set(xlabel='Class', ylabel='Communities')
		plot.set_title(f'Heatmap - Class')
		canvas = FigureCanvasTkAgg(self.fig, master=self.frameTOP)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(row=0, column=0)
		self.frameTOP.pack()
		self.framBOT = tk.Frame(self.master)
		ttk.Button(self.framBOT,text="Save Heatmap",command = lambda : self.saveImage()).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Save Matrix",command = lambda : self.saveMatrix(self.df)).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Quit",command = lambda : self.master.destroy()).pack(side=LEFT)
		self.framBOT.pack(pady=10, side=BOTTOM)

	def saveImage(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			self.fig.savefig(filename, format='png', bbox_inches='tight')

	def saveMatrix(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.iloc[::-1].to_csv(filename, sep='\t')
			
	def analysisGO(self, master, indexTop):
		list_GO = []
		dic = {}
		for i in self.dicQuery.keys():
			self.dfBank = pd.read_csv(os.path.join('GO_DB', f'{self.comboGOBank.get()}.tsv'), sep='\t', index_col=[0])
			geneGOList = self.dfBank.iloc[:, 0]
			geneGOList = geneGOList.apply(lambda x: x.replace('[\'', '').replace('\']', '').split('\', \''))
			self.dfBank['Gene Signature'] = geneGOList
			tot_genes = geneGOList.tolist()
			tot_genes = [item for sublist in tot_genes for item in sublist]
			tot_genes = set(tot_genes)
			nb_totgen = len(tot_genes)
			similarity = geneGOList.apply(lambda x: GoGUI.fisher(GoGUI, self.dicQuery[i], x, nb_totgen))
			print(similarity)
			similarity = similarity.sort_values()
			similarity.name = 'p-values'
			sub_similarity = similarity.index.values
			dise_comm = geneGOList.loc[sub_similarity]
			dise_comm = dise_comm.apply(lambda x: GoGUI.inter(GoGUI, x, self.dicQuery[i]))
			dise_comm = dise_comm.dropna()
			dise_comm.name = 'Common genes'
			self.table = pd.concat([similarity, dise_comm], axis=1)
			self.table = self.table.dropna()
			print(f'### TOP{indexTop} GO TERMS - Communities {i}')
			print(self.table.index[:indexTop].tolist())
			dic[i] = self.table
			list_GO.extend(self.table.index[:indexTop].tolist())
		print(list(set(list_GO)))
		print(list_GO)
		l_list_GO =  list(set(list_GO))
		l = []
		minV = 1
		maxV = 0
		for i in dic.keys():
			ll = []
			for u in list(set(list_GO)):
				ll.append(dic[i].loc[str(u), 'p-values'])
			if minV > min(ll):
				minV = min(ll)
			if maxV < max(ll):
				maxV = max(ll) 
			l.append(ll)
		ar = np.array(l)
		dfT = pd.DataFrame(ar, index=dic.keys(), columns=list(set(list_GO)))
		HeatmapGOAll(tk.Toplevel(master), self.dicQuery, dfT, f'{self.titleGO} - TOP{indexTop}', self.comboGOBank.get(), minV, maxV, False)
	
	def fisher(self, input, disease, nb_totgen):
		input = set(input)
		disease = set(disease)
		inputxdisease = len(disease & input)
		inputxtot = len(input) - inputxdisease
		diseasewoinput = len(disease) - inputxdisease
		totwoinput = nb_totgen - len(disease) - inputxtot
		odd, pvalues = stats.fisher_exact([[inputxdisease, diseasewoinput], [inputxtot, totwoinput]])
		return pvalues

	def inter(self, liste1, liste2):
		liste1 = set(liste1)
		liste2 = set(liste2)
		return list(liste2 & liste1)

	def gexel_active_heatmap(self, event, communityDisplay, classDisplay):
		"""
		When mouse over a gexel
		"""
		self.countDisplayClass.set(classDisplay)
		self.countDisplayCommunities.set(communityDisplay)

	def drawingHeatmap(self):
		self.can.delete('all')
		x1 = 1
		x2 = self.indexMiniSizeGexelX
		y1 = 1
		y2 = self.indexMiniSizeGexelY
		print(self.df.index.values)
		print(self.df.columns.values)
		l = []
		for index, row in self.df.iterrows():
			x1 = 1
			x2 = self.indexMiniSizeGexelX
			for x in self.df.columns.values.tolist():
				currentValue = int(row[x])
				if currentValue == 1:
					can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexelX, y1+self.indexMiniSizeGexelY, x2+self.indexMiniSizeGexelX, y2+self.indexMiniSizeGexelY,
						fill = 'green', activefill = 'white', activeoutline = 'white')
					self.can.tag_bind(can_temp, "<Enter>", lambda event,
						communityDisplay = x,
						classDisplay = index 
						: self.gexel_active_heatmap(event, communityDisplay, classDisplay))
					self.can.tag_bind(can_temp, "<Leave>", lambda event,
						communityDisplay = '',
						classDisplay = '' 
						: self.gexel_active_heatmap(event, communityDisplay, classDisplay))
				else:
					can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexelX, y1+self.indexMiniSizeGexelY, x2+self.indexMiniSizeGexelX, y2+self.indexMiniSizeGexelY,
						fill = 'grey')
				x1 += self.indexMiniSizeGexelX
				x2 += self.indexMiniSizeGexelX
			y1 += self.indexMiniSizeGexelY
			y2 += self.indexMiniSizeGexelY
			
class wrapperClass(object):
	def __init__(self, master, runObject):
		self.master = master
		self.runObject = runObject
		self.frameArg = tk.Frame(self.master)
		self.frameTOP = tk.Frame(self.frameArg)
		self.frameListboxClass = tk.Frame(self.frameTOP)
		self.frameListboxClass_text = tk.Frame(self.frameListboxClass)
		self.textOverview = Label(self.frameListboxClass_text, text='Select datasets for comparison :')
		self.textOverview.pack(side=TOP)
		listClass = Listbox(self.frameListboxClass_text, height=5, width=43, selectmode='multiple')
		mainGUI.fillListbox(mainGUI, listClass, self.runObject.dicSample.keys())
		listClass.pack(fill=X, side=BOTTOM)
		self.frameListboxClass_text.pack(side=TOP)
		self.frameButtonListBox = tk.Frame(self.frameListboxClass)
		self.selectAllButton = ttk.Button(self.master, text='Select All',
			command = lambda : listClass.selection_set(0, END))
		self.selectAllButton.pack(fill=X)
		self.resetAllButton = ttk.Button(self.master, text='Reset All',
			command = lambda : listClass.selection_clear(0, END))
		self.resetAllButton.pack(fill=X)
		self.frameButtonListBox.pack(side=BOTTOM)
		self.frameListboxClass.pack(side=TOP)
		self.frameWeight = tk.Frame(self.frameTOP)
		self.varWeight = IntVar()
		self.labelWeight = tk.Label(self.frameWeight, text='Weight : ')
		self.labelWeight.pack(side=LEFT)
		self.checkWeight = Checkbutton(self.frameWeight, variable=self.varWeight)
		self.checkWeight.select()
		self.checkWeight.pack(side=LEFT)
		self.frameWeight.pack(side=BOTTOM)
		self.frameTOP.pack(side=TOP)
		self.thresholdCommunities = IntVar()
		self.frameThresoldCommunities = tk.Frame(self.frameArg)
		self.frameRandom = tk.Frame(self.frameArg)
		self.varRandom = IntVar()
		self.labelRandom = tk.Label(self.frameRandom, text='Multiple iterations (15 times) : ')
		self.labelRandom.pack(side=LEFT)
		self.checkRandom = Checkbutton(self.frameRandom, variable=self.varRandom)
		self.checkRandom.select()
		self.checkRandom.pack(side=LEFT)
		self.frameRandom.pack(side=TOP)
		self.thresoldCommunities = Label(self.frameThresoldCommunities, text='Threshold <= :')
		self.thresoldCommunities.pack(side=LEFT)
		self.entryThresholdCommunities = Entry(self.frameThresoldCommunities, textvar=self.thresholdCommunities, width=3)
		self.entryThresholdCommunities.pack(side=LEFT)
		self.thresholdCommunities.set(0)
		self.percentThresoldCommunities = Label(self.frameThresoldCommunities, text='%.')
		self.percentThresoldCommunities.pack(side=LEFT)
		self.frameThresoldCommunities.pack(side=BOTTOM)
		self.frameArg.pack(side=TOP)
		self.runcommunitiesButton = ttk.Button(self.master, text='Run Communities',
			command = lambda : self.wrapperGexelMinBatch([listClass.get(idx) for idx in listClass.curselection()]))
		self.runcommunitiesButton.pack()

	def mergeTG(self, x):
		tf_filter = x['TF'].split(' | ')[0]
		tg = f'{tf_filter}'
		return tg

	def mergeTF(self, x, name):
		community = x['Community']
		tf = f'{name}_{community}'
		return tf

	def wrapperGexelMinBatch(self, list_selected):
		if len(list_selected) == 0:
			tk.messagebox.showwarning('No sample selected !', 'You have to select a sample.', icon='warning')
		else:
			self.selectAllButton.destroy()
			self.resetAllButton.destroy()
			self.frameArg.destroy()
			self.runcommunitiesButton.destroy()
			self.gexelMinBatch(self.master, list_selected)
			self.master.title('Arguments for class')

	def gexelMinBatch(self, master, list_selected):
		self.framegexelMinBatch = tk.Frame(master)
		l = list_selected
		list_val = []
		for i in l:
			list_val.append(i)
		for row in range(len(l)):
			for column in range(2):
				if column == 0:
					cell = Label(self.framegexelMinBatch, text=f'{l[row]} Minimal  n° of gexels : ')
				else:
					list_val[row] = IntVar()
					cell = Entry(self.framegexelMinBatch, textvar=list_val[row], width=3)
					list_val[row].set(int(self.runObject.minPattern))
				cell.grid(row=row, column=column)
		self.framegexelMinBatch.pack(side=TOP)
		self.runcommunitiesButton = ttk.Button(master, text='Run',
			command = lambda : self.dic_generate(list_val, l))
		self.runcommunitiesButton.pack(side=BOTTOM, fill=X)

	def dic_generate(self, listL, l):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':	
			self.dicClass = {}
			for i in range(len(listL)):
				print(f'{l[i]} = {listL[i].get()}')
				self.dicClass[l[i]]=listL[i].get()
			print(self.dicClass)
			self.framegexelMinBatch.destroy()
			self.louvainClass(self.varWeight.get(), self.varRandom.get(), filename)
		else:
			pass

	def louvainClass(self, weightIndex, randomIndex, filename):
		list_ = []
		for name in self.dicClass.keys():
			if int(self.dicClass[name]) == self.runObject.minPattern:
				print(f'{name} - {len(self.runObject.dicSample[name].diff.columns.values)} - Gexel min {int(self.dicClass[name])}')
				dic_valCount_df = {}
				dic_valCount_index = {}
				indexRandom = 1
				if randomIndex == 1:
					indexRandom = 15
				for i in range(indexRandom):
					df = self.runObject.dicSample[name].matrixLouvain.copy()
					df['Community'] = df['TF']
					weightString = ''
					indexNames = df[(df['TF'] == df['TG'])].index
					df.drop(indexNames, inplace=True)
					indexNames2 = df[(df['Weight'] <= self.thresholdCommunities.get())].index
					df.drop(indexNames2, inplace=True)
					if weightIndex == 1:
						nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr='Weight')
						partition = louvain.best_partition(nw, weight='Weight')
						weightString = 'Weight'
					elif weightIndex == 0:
						nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr=None)
						partition = louvain.best_partition(nw, weight='None')
						weightString = 'without Weight'
					df['Community'] = df['Community'].map(partition)
					tmp = df['Community'].value_counts().to_string()
					if tmp not in dic_valCount_df.keys():
						dic_valCount_df[tmp] = df
						dic_valCount_index[tmp] = 1
					else:
						dic_valCount_index[tmp] += 1
				valCount_index_sort = sorted(dic_valCount_index.items(), key=lambda x: x[1], reverse=True)
				if indexRandom != 1:
					print(name)
					mainGUI.logFileRandomCommunities(mainGUI, valCount_index_sort, indexRandom)
				df_class = dic_valCount_df[valCount_index_sort[0][0]]
				df_class['TG'] = df_class.apply(self.mergeTG, axis=1)#"{}{}".format(x, 'str'))
				df_class['TF'] = df_class.apply(lambda x : self.mergeTF(x, name), axis=1)
				list_.append(df_class)
			else:
				print(f'Agglomerative - {int(self.dicClass[name])}')
				dic_geneSampleDiffCopy = copy.deepcopy(self.runObject.dicSample[name].dicSampleDiff)
				for geneSample in dic_geneSampleDiffCopy.keys():
					tempListCoordinateCopy = []
					indexFind = 0
					for coordinate in dic_geneSampleDiffCopy[geneSample].dic_coordinate_count.keys():
						if dic_geneSampleDiffCopy[geneSample].dic_coordinate_count[coordinate] >= 1:
							x = int(coordinate.split('x')[0])
							y = int(coordinate.split('x')[1])
							tempListCoordinateCopy.append([x, y])
					tempListCoordinateCopy = np.asarray(tempListCoordinateCopy)
					if len(tempListCoordinateCopy) >= int(self.dicClass[name]):
						clustering_temp = AgglomerativeClustering(n_clusters=None, affinity='euclidean',
							linkage='single', distance_threshold=1.5)
						labels = clustering_temp.fit_predict(tempListCoordinateCopy)
						(unique, counts) = np.unique(labels, return_counts=True)
						frequencies = np.asarray((unique, counts)).T
						for i in frequencies:
							index = i[0]
							count = i[1]
							if count >= int(self.dicClass[name]):
								indexFind += 1
								listTempCopy = []
								for o in tempListCoordinateCopy[labels==index]:
									listTempCopy.append(f'{o[0]}x{o[1]}')
								dic_geneSampleDiffCopy[geneSample].dic_pattern[indexFind] = listTempCopy
				print('Similarity Matrix')
				startTime_sim = time.time()
				list_similarity = []
				list_column = []
				for geneQuery in dic_geneSampleDiffCopy.keys():
					for numberPatternQuery in dic_geneSampleDiffCopy[geneQuery].dic_pattern.keys():
						set_gene_query = set(dic_geneSampleDiffCopy[geneQuery].dic_pattern[int(numberPatternQuery)])
						list_temp = []
						list_column.append(f'{geneQuery} | {numberPatternQuery}')
						index_find = 0
						for geneComp in dic_geneSampleDiffCopy.keys():
							for numberPatternComp in dic_geneSampleDiffCopy[geneComp].dic_pattern.keys():
								if geneQuery == geneComp:
									if numberPatternQuery == numberPatternComp:
										list_temp.append(100)
									else:
										list_temp.append(np.NaN)
								else:
									if self.runObject.method == 'Tanimoto':
										inter = len(set_gene_query.intersection(set(dic_geneSampleDiffCopy[geneComp].dic_pattern[int(numberPatternComp)])))
										union = len(set_gene_query.union(set(dic_geneSampleDiffCopy[geneComp].dic_pattern[int(numberPatternComp)])))
										similarityValue = float(inter/union)*100
										if round(similarityValue, 2) >= 0:
											list_temp.append(round(similarityValue, 2))
											index_find += 1
										else:
											list_temp.append(np.NaN)
									elif self.runObject.method == 'Dice':
										inter = len(set_gene_query.intersection(set(dic_geneSampleDiffCopy[geneComp].dic_pattern[int(numberPatternComp)])))
										union = len(set_gene_query.union(set(dic_geneSampleDiffCopy[geneComp].dic_pattern[int(numberPatternComp)])))
										similarityValue = float((2*inter)/(len(set_gene_query) + len(set(dic_geneSampleDiffCopy[geneComp].dic_pattern[int(numberPatternComp)]))))*100
										if round(similarityValue, 2) >= 0:
											list_temp.append(round(similarityValue, 2))
											index_find += 1
										else:
											list_temp.append(np.NaN)
						list_similarity.append(list_temp)
				dfSim = pd.DataFrame(list_similarity, columns=list_column, index=list_column)
				dfSim = dfSim.transpose()
				stopTime_sime = time.time()
				print(f'Time to make similarity  {stopTime_sime-startTime_sim} seconds')
				print('Generate matrix for louvain')
				dico_tot ={}
				dfSim.apply(lambda x: StartGUI.dico_df(StartGUI, x, dico_tot))
				dico_fin = StartGUI.dico_final(StartGUI, dico_tot)
				dfLouvain = pd.DataFrame(dico_fin, columns=['TF', 'TG', 'Weight'])
				print(dfLouvain)
				print(f'{name} - {len(self.runObject.dicSample[name].diff.columns.values)} - Gexel min {int(self.dicClass[name])}')
				dic_valCount_df = {}
				dic_valCount_index = {}
				indexRandom = 1
				if randomIndex == 1:
					indexRandom = 15
				for i in range(indexRandom):
					dfLouvain['Community'] = dfLouvain['TF']
					weightString = ''
					indexNames = dfLouvain[(dfLouvain['TF'] == dfLouvain['TG'])].index
					dfLouvain.drop(indexNames, inplace=True)
					indexNames2 = dfLouvain[(dfLouvain['Weight'] <= self.thresholdCommunities.get())].index
					dfLouvain.drop(indexNames2, inplace=True)
					if weightIndex == 1:
						nw = nx.from_pandas_edgelist(dfLouvain, source='TF', target='TG', edge_attr='Weight')
						partition = louvain.best_partition(nw, weight='Weight')
						weightString = 'Weight'
					elif weightIndex == 0:
						nw = nx.from_pandas_edgelist(dfLouvain, source='TF', target='TG', edge_attr=None)
						partition = louvain.best_partition(nw, weight='None')
						weightString = 'without Weight'
					dfLouvain['Community'] = dfLouvain['Community'].map(partition)
					tmp = dfLouvain['Community'].value_counts().to_string()
					if tmp not in dic_valCount_df.keys():
						dic_valCount_df[tmp] = dfLouvain
						dic_valCount_index[tmp] = 1
					else:
						dic_valCount_index[tmp] += 1
				valCount_index_sort = sorted(dic_valCount_index.items(), key=lambda x: x[1], reverse=True)
				if indexRandom != 1:
					print(name)
					mainGUI.logFileRandomCommunities(mainGUI,valCount_index_sort, indexRandom)
				df_class = dic_valCount_df[valCount_index_sort[0][0]]
				df_class['TG'] = df_class.apply(self.mergeTG, axis=1)#"{}{}".format(x, 'str'))
				df_class['TF'] = df_class.apply(lambda x : self.mergeTF(x, name), axis=1)
				maxCommunitiesSample = max(df_class['Community'])
				print(f'For {name} Number of communities {maxCommunitiesSample+1}')
				list_.append(df_class)
		print('########################\nAll Communities\n########################')
		dic_valCount_df_all = {}
		dic_valCount_index_all = {}
		indexRandom = 1
		if randomIndex == 1:
			indexRandom = 15
		for i in range(indexRandom):
			df = pd.concat(list_, axis=0)
			df['Community'] = df['TF']
			weightString = ''
			indexNames = df[(df['TF'] == df['TG'])].index
			df.drop(indexNames, inplace=True)
			indexNames2 = df[(df['Weight'] <= self.thresholdCommunities.get())].index
			df.drop(indexNames2, inplace=True)
			if weightIndex == 1:
				nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr='Weight')
				partition = louvain.best_partition(nw, weight='Weight')#, random_state=100)
				weightString = 'Weight'
			elif weightIndex == 0:
				nw = nx.from_pandas_edgelist(df, source='TF', target='TG', edge_attr=None)
				partition = louvain.best_partition(nw, weight='None')#, random_state=100)
				weightString = 'without Weight'
			df['Community'] = df['Community'].map(partition)
			tmp = df['Community'].value_counts().to_string()
			if tmp not in dic_valCount_df_all.keys():
				dic_valCount_df_all[tmp] = df
				dic_valCount_index_all[tmp] = 1
			else:
				dic_valCount_index_all[tmp] += 1
		valCount_index_sort = sorted(dic_valCount_index_all.items(), key=lambda x: x[1], reverse=True)
		if indexRandom != 1:
			mainGUI.logFileRandomCommunities(mainGUI, valCount_index_sort, indexRandom)
		print('Write')
		dic_valCount_df_all[valCount_index_sort[0][0]].to_csv(filename, sep='\t')
		df_heatmap = dic_valCount_df_all[valCount_index_sort[0][0]].copy()
		df_heatmap = df_heatmap.sort_values(by ='Community', ascending=False)
		print(df_heatmap.head())
		print(len(df_heatmap))
		dic_valCount_df_all[valCount_index_sort[0][0]] = dic_valCount_df_all[valCount_index_sort[0][0]].sort_values(by ='Weight', ascending=False)
		dic_valCount_df_all[valCount_index_sort[0][0]] = dic_valCount_df_all[valCount_index_sort[0][0]].drop_duplicates(subset =['TF', 'TG', 'Community'], keep = 'first')#.set_index('TF')['Community'].to_dict()
		print('Write filtered')
		dic_valCount_df_all[valCount_index_sort[0][0]].to_csv(f'{filename[:-3]}_filtered.tsv', sep='\t')
		maxCommunities = max(dic_valCount_df_all[valCount_index_sort[0][0]]['Community'])
		print(f'Number of communities {maxCommunities+1}')
		print('End')
		dicTest = {}
		for index, row in dic_valCount_df_all[valCount_index_sort[0][0]].iterrows():
			if row['Community'] not in dicTest.keys():
				dicTest[row['Community']] = [row['TG'].split(' ')[0]]
			else:
				dicTest[row['Community']].append(row['TG'].split(' ')[0])
		dicTest2 = {}
		print(df_heatmap.head())
		print(len(df_heatmap))
		for index, row in df_heatmap.iterrows():
			if row['TF'] not in dicTest2.keys():
				dicTest2[row['TF']] = [row['Community']]
			else:
				dicTest2[row['TF']].append(row['Community'])
		print(len(dicTest2.keys()))
		l = list(set(df_heatmap['Community']))
		ll = []
		print(dicTest2.keys())
		for i in dicTest2.keys():
			lll = []
			for element in l:
				if element in dicTest2[i]:
					lll.append(1)
				else:
					lll.append(0)
			ll.append(lll)
		ar = np.array(ll)
		self.dfTT = pd.DataFrame(ar, index=dicTest2.keys(), columns=l)
		print(self.dfTT)
		self.runcommunitiesButton.destroy()
		resultGOClass(self.master, self.runObject, self.dfTT.transpose(), dicTest)

	def wrapperHeatmapClass(self, master):
		self.heatmapClassFrame = tk.Toplevel(master)
		self.heatmapClass(self.heatmapClassFrame)

	def heatmapClass(self, master):
		self.master = master
		self.frameTOP = tk.Frame(self.master)
		self.master.title(f'Heatmap Class')
		self.fig = Figure()
		ax = self.fig.add_subplot(111)
		cmap = LinearSegmentedColormap.from_list('Custom', ['#D3D3D3', '#319819'], 2)
		plot = sb.heatmap(self.dfTT, cmap=cmap, linewidths=0.3,
			linecolor='black', ax=ax, xticklabels=True, yticklabels=True)
		plot.set_ylim(0 , len(self.dfTT))
		plot.xaxis.set_ticks_position('top')
		plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
		colorbar = ax.collections[0].colorbar
		colorbar.set_ticks([0.25, 0.75])
		colorbar.set_ticklabels(['Absence', 'Presence'])
		plot.set_yticklabels(plot.get_yticklabels(), rotation=0)
		plot.set(xlabel='Class', ylabel='Communities')
		plot.set_title(f'Heatmap - Class')
		canvas = FigureCanvasTkAgg(self.fig, master=self.frameTOP)
		plot_widget = canvas.get_tk_widget()
		plot_widget.grid(row=0, column=0)
		self.frameTOP.pack()
		self.framBOT = tk.Frame(self.master)
		ttk.Button(self.framBOT,text="Save Heatmap",command = lambda : self.saveImage()).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Save Matrix",command = lambda : self.saveMatrix(self.dfTT)).pack(side=LEFT)
		ttk.Button(self.framBOT,text="Quit",command = lambda : self.master.destroy()).pack(side=LEFT)
		self.framBOT.pack(pady=10, side=BOTTOM)

	def geneOntologyClass(self, dicQuery, listGOBank, title):
		self.master.title(f'Gene Ontology Analysis - Class')
		self.guiGOClass(self.master, dicQuery, listGOBank, title)

	def saveImage(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			self.fig.savefig(filename, format='png', bbox_inches='tight')

	def saveMatrix(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.iloc[::-1].to_csv(filename, sep='\t')

	def geneOntologyClass(self, dicQuery, listGOBank, title):
		self.master.title(f'Gene Ontology Analysis - Class')
		self.guiGOClass(self.master, dicQuery, listGOBank, title)

	def guiGOClass(self, master, dicQuery, listGOBank, title):
		self.dic_Com_ListGene = {}
		sb.set(font_scale=0.75)
		self.dicQuery = dicQuery
		self.listGOBank = listGOBank
		self.titleGO = title
		master.minsize(width=600, height=100)
		self.frameGOText = tk.Frame(master)
		self.textOverview = Label(self.frameGOText, text='Select a GO terms database and press \'Run\'.')
		self.textOverview.pack()
		self.frameListLabal = tk.Frame(self.frameGOText)
		self.textTopGO = Label(self.frameListLabal, text='Top GO terms :')
		self.textTopGO.pack(side=LEFT)
		self.listCombobox = ttk.Combobox(self.frameListLabal,
				values=[5, 10, 15], state='readonly', height=3, width=3)
		self.listCombobox.current(0)
		self.listCombobox.pack(side=RIGHT)
		self.frameListLabal.pack()
		self.frameGOText.pack()
		self.frameGOButton = tk.Frame(master)
		self.frameGOButtonTop = tk.Frame(self.frameGOButton)
		self.comboGOBank = ttk.Combobox(self.frameGOButtonTop,
				values=self.listGOBank, state='readonly', height=5)
		self.comboGOBank.pack(side=LEFT)
		self.comboGOBank.current(0)
		self.runButton = ttk.Button(self.frameGOButtonTop, text="Run", command= lambda : self.analysisGO(master, int(self.listCombobox.get())))
		self.runButton.pack(side=LEFT)
		self.frameGOButtonTop.pack(side=TOP)
		self.frameGOButtonBot = tk.Frame(self.frameGOButton)
		self.heatmapButton = ttk.Button(self.frameGOButtonBot, text="Heatmap",
			command= lambda : self.wrapperHeatmapClass(master))
		self.heatmapButton.pack(side=LEFT)
		self.quitGOButton = ttk.Button(self.frameGOButtonBot, text="Quit", command= lambda : master.destroy())
		self.quitGOButton.pack(side=RIGHT)
		self.frameGOButtonBot.pack(side=BOTTOM)
		self.frameGOButton.pack(fill=X, pady=10, side=BOTTOM)

	def analysisGO(self, master, indexTop):
		list_GO = []
		dic = {}
		for i in self.dicQuery.keys():
			self.dfBank = pd.read_csv(os.path.join('GO_DB', f'{self.comboGOBank.get()}.tsv'), sep='\t', index_col=[0])
			geneGOList = self.dfBank.iloc[:, 0]
			geneGOList = geneGOList.apply(lambda x: x.replace('[\'', '').replace('\']', '').split('\', \''))
			self.dfBank['Gene Signature'] = geneGOList
			tot_genes = geneGOList.tolist()
			tot_genes = [item for sublist in tot_genes for item in sublist]
			tot_genes = set(tot_genes)
			nb_totgen = len(tot_genes)
			similarity = geneGOList.apply(lambda x: GoGUI.fisher(GoGUI, self.dicQuery[i], x, nb_totgen))
			print(similarity)
			similarity = similarity.sort_values()
			similarity.name = 'p-values'
			sub_similarity = similarity.index.values
			dise_comm = geneGOList.loc[sub_similarity]
			dise_comm = dise_comm.apply(lambda x: GoGUI.inter(GoGUI, x, self.dicQuery[i]))
			dise_comm = dise_comm.dropna()
			dise_comm.name = 'Common genes'
			self.table = pd.concat([similarity, dise_comm], axis=1)
			self.table = self.table.dropna()
			print(f'### TOP{indexTop} GO TERMS - Communities {i}')
			print(self.table.index[:indexTop].tolist())
			dic[i] = self.table
			list_GO.extend(self.table.index[:indexTop].tolist())
		print(list(set(list_GO)))
		print(list_GO)
		l_list_GO =  list(set(list_GO))
		l = []
		minV = 1
		maxV = 0
		for i in dic.keys():
			ll = []
			for u in list(set(list_GO)):
				ll.append(dic[i].loc[str(u), 'p-values'])
			if minV > min(ll):
				minV = min(ll)
			if maxV < max(ll):
				maxV = max(ll) 
			l.append(ll)
		ar = np.array(l)
		dfT = pd.DataFrame(ar, index=dic.keys(), columns=list(set(list_GO)))
		HeatmapGOAll(tk.Toplevel(master), self.dicQuery, dfT, f'{self.titleGO} - TOP{indexTop}', self.comboGOBank.get(), minV, maxV, False)
	
	def fisher(self, input, disease, nb_totgen):
		input = set(input)
		disease = set(disease)
		inputxdisease = len(disease & input)
		inputxtot = len(input) - inputxdisease
		diseasewoinput = len(disease) - inputxdisease
		totwoinput = nb_totgen - len(disease) - inputxtot
		odd, pvalues = stats.fisher_exact([[inputxdisease, diseasewoinput], [inputxtot, totwoinput]])
		return pvalues

	def inter(self, liste1, liste2):
		liste1 = set(liste1)
		liste2 = set(liste2)
		return list(liste2 & liste1)


class matrixGUI(object):
	"""
	Panel allows user choose which matrix to display
	"""
	def __init__(self, title, master, runObject, countDisplayStringCoor, countDisplayString, nameSample, dicGexel, dicGeneSample, maxCountGexel, minCountGexel, gradient, listGOBank, currentLouvain):
		self.master = master
		self.runObject = runObject
		self.nameSample = nameSample
		self.dicGexel = dicGexel
		self.dicGeneSample = dicGeneSample
		self.maxCountGexel = round(maxCountGexel, 2)
		self.minCountGexel = round(minCountGexel, 2)
		self.listGeneCommunities = []
		if 'All - Communities' in title:
			self.title = title[6:]
		else:
			self.title = title
		self.indexLog = 0
		self.master.title(self.title)
		self.listGOBank = listGOBank
		self.frame = tk.Frame(self.master)
		self.widthGUI  = self.master.winfo_screenwidth()
		self.heightGUI = self.master.winfo_screenheight()
		self.minValueGradient = IntVar()
		self.maxValueGradient = IntVar()
		self.countDisplayStringCoor = countDisplayStringCoor
		self.countDisplayString = countDisplayString
		self.currentLouvain = currentLouvain
		if (3*self.heightGUI/4)/(self.runObject.sizeY+4) > (self.widthGUI/2)/(self.runObject.sizeX+4):
			self.indexMiniSizeGexel = (self.widthGUI/2)/(self.runObject.sizeX+4)
		else:
			self.indexMiniSizeGexel = (3*self.heightGUI/4)/(self.runObject.sizeY+4)
		self.thresholdMin = StringVar()
		self.thresholdMax = StringVar()
		self.thresholdDownDiff = DoubleVar()
		self.thresholdDownDiff.set(self.runObject.downDiffThreshold)
		self.thresholdUpDiff = DoubleVar()
		self.thresholdUpDiff.set(self.runObject.upDiffThreshold)
		self.thresholdSimilarity = StringVar()
		self.gradient = gradient
		if int(self.indexMiniSizeGexel*self.runObject.sizeX) > 510:
			self.master.minsize(int(self.indexMiniSizeGexel*self.runObject.sizeX)+270, int(self.indexMiniSizeGexel*self.runObject.sizeY)+100)
		else:
			self.master.minsize(780, 610)
		self.frameMatrix = tk.Frame(self.frame)
		self.can = Canvas(self.frameMatrix, width=((self.runObject.sizeX+1)*self.indexMiniSizeGexel)+1, 
			height = ((self.runObject.sizeY+1)*self.indexMiniSizeGexel)+1)
		w = Canvas.winfo_width(self.can)
		h = Canvas.winfo_height(self.can)		
		self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
			(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
		self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
		self.colorGradientFrame = tk.Frame(self.frameMatrix)
		self.colorScale = Canvas(self.colorGradientFrame, width=255*2, height=15)
		if self.gradient == 'patternGradient' or self.gradient == 'communitiesGradient':
			self.colorList = []
		else:
			self.colorList = self.colorGradient()
		if self.gradient == 'communitiesGradient':
			self.comMatrix('All')
		else:
			self.drawingMatrix(self.colorList)
		self.can.pack(padx = 10)
		self.colorScale.pack()
		self.colorGradientFrame.pack()
		self.index_hide = 1
		if self.gradient != 'communitiesGradient':
			self.colorGradientFrame2 = tk.Frame(self.frameMatrix)
			self.labelMinValue = Label(self.colorGradientFrame2, textvariable=self.minValueGradient)
			self.labelMinValue.pack(side=LEFT)
			if self.gradient == 'patternGradient':
				self.minValueGradient.set('')
				self.maxValueGradient.set('')
			self.tempFrame = tk.Frame(self.colorGradientFrame2)
			self.tempFrame.pack(padx=250)
			self.labelMaxValue = Label(self.colorGradientFrame2, textvariable=self.maxValueGradient)
			self.labelMaxValue.pack(side=RIGHT)
			self.colorGradientFrame2.pack()
		self.buttonFrame = tk.Frame(self.frameMatrix)
		if self.gradient == 'communitiesGradient':
			self.list_number_communities = []
			self.list_number_communities.append('All')
			self.list_number_communities.extend(self.currentLouvain['Community'].value_counts().index.tolist())
			self.menu()
		else:
			self.showButton = ttk.Button(self.buttonFrame, text = 'Menu', width = 12, command = lambda : self.hideMenu())
			self.showButton.pack(side=LEFT)
		self.quitButton = ttk.Button(self.buttonFrame, text = 'Quit', width = 12, command = self.quit)#self.close_windows)
		self.quitButton.pack(side=RIGHT)
		self.saveCanvasButton = ttk.Button(self.buttonFrame, text = 'Save Image', width = 12, command = lambda : self.saveCanvas())
		self.saveCanvasButton.pack(side=RIGHT)
		self.frameMatrix.pack(side=LEFT, padx=10)
		self.buttonFrame.pack()
		self.frame.pack()

	def saveCanvas(self):
		filename = asksaveasfilename(defaultextension=".png", filetypes=(('Portable Network Graphics', '*.png'),))
		if filename != '':
			if self.gradient != 'patternGradient' and self.gradient != 'communitiesGradient':
				self.colorGradient()
				self.draw.text((self.indexMiniSizeGexel, 4), str(self.minValueGradient.get()),(0,0,0))
				self.draw.text((550+(2*self.indexMiniSizeGexel), 4), str(self.maxValueGradient.get()),(0,0,0))
			self.saveCanvasImage.save(filename)

	def saveCommunities(self):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			self.currentLouvain.to_csv(filename, sep='\t')

	def menu(self):
		"""
		frame menu content
		"""
		self.frameMenu = tk.Frame(self.frame)

		if self.gradient == 'communitiesGradient':
			self.temp_title = ''
			self.displayNumberPerCommunities = StringVar()
			self.frameCommunities = tk.Frame(self.frameMenu, highlightthickness=1, highlightbackground='black',
				highlightcolor='black')
			self.frameCommuSelected = tk.Frame(self.frameCommunities)
			self.labelCommunities = Label(self.frameCommuSelected, text='Communities selected : ')
			self.labelCommunities.pack(side=LEFT)
			self.comboCommunities = ttk.Combobox(self.frameCommuSelected,
				values=self.list_number_communities, state='readonly', width=2, height=5)
			self.comboCommunities.pack(side=LEFT)
			self.comboCommunities.current(0)
			self.frameCommuSelected.pack(side=TOP)
			self.searchButtonCommunities = ttk.Button(self.frameCommunities, text = 'Search',
				command = lambda : self.comMatrix(self.comboCommunities.get()))
			self.searchButtonCommunities.pack(side=BOTTOM)
			self.frameNumberGenePerCommunities = tk.Frame(self.frameCommunities)
			self.labelGenePerCommunities = Label(self.frameNumberGenePerCommunities, text='Number of genes : ')
			self.labelDisplayGenePerCommunities = Label(self.frameNumberGenePerCommunities,
				textvariable=self.displayNumberPerCommunities)
			self.labelGenePerCommunities.pack(side=LEFT)
			self.labelDisplayGenePerCommunities.pack(side=LEFT, fill = X)
			self.frameNumberGenePerCommunities.pack(side=BOTTOM)
			self.listBoxCommunities = Listbox(self.frameCommunities, height=10, width=15)
			self.listBoxCommunities.bind('<Double-Button>', self.clipboard)
			self.listBoxCommunities.pack(side=BOTTOM, pady = 10)
			self.frameCommunities.pack(side=TOP)
			self.frameCommunities3 = tk.Frame(self.frameMenu)
			self.GOButton = ttk.Button(self.frameCommunities3, text = 'Gene Ontology', command = lambda : self.geneOntology(self.listGeneCommunities, self.listGOBank, self.master.title()))
			self.GOButton.pack(side=BOTTOM)
			self.saveComunitiesButton = ttk.Button(self.frameCommunities3, text = 'Save', command = lambda : self.saveCommunities())
			self.saveComunitiesButton.pack(side=BOTTOM)
			self.frameCommunities3.pack(side=BOTTOM)
		else:
			self.frameSearchTool = tk.Frame(self.frameMenu)
			self.frameSearchButton = tk.Frame(self.frameSearchTool)
			self.frameSearch = tk.Frame(self.frameSearchTool)
			self.searchGene = AutocompleteEntry(self.dicGeneSample.keys(), self.frameSearchButton, self.frameSearch)
			self.searchGene.pack(side = TOP, fill=BOTH)
			if self.gradient == 'patternGradient':
				indexPatternSelected = ''
				self.temp_title = ''
				self.searchButton = ttk.Button(self.frameSearchButton, text = 'Search', command = lambda : self.patternMatrix())
				self.searchButton.pack(side = TOP, fill=X)
				self.frameDisplaySelectedPattern = tk.Frame(self.frameSearchButton)
				self.patternDisplayString = StringVar()
				self.labelPattern = Label(self.frameDisplaySelectedPattern, text='Pattern selected : ')
				self.labelPatternDisplay = Label(self.frameDisplaySelectedPattern, textvariable=self.patternDisplayString)
				self.labelPattern.pack(side=LEFT)
				self.labelPatternDisplay.pack(side=LEFT, fill = X)
				self.labelPattern = Label(self.frameDisplaySelectedPattern, text='Pattern selected : ')
				self.frameDisplaySelectedPattern.pack(side = TOP, fill=X)
				self.geneQuery = StringVar()
				self.frameThresoldSimilarity = tk.Frame(self.frameSearchButton)
				self.labelThresholdSimilirity = Label(self.frameThresoldSimilarity, text='Similarity > :')
				self.labelThresholdSimilirity.pack(side=LEFT)
				self.listValue = [90, 80, 70, 60, 50, 40, 30, 20, 10, 0]
				self.comboThresholdSimilarity = ttk.Combobox(self.frameThresoldSimilarity,
					values=self.listValue, state='readonly', height=5, width=8)
				self.comboThresholdSimilarity.current(0)
				self.comboThresholdSimilarity.pack(side=LEFT)
				self.labelThresholdSimilirityPercent = Label(self.frameThresoldSimilarity, text=' %.')
				self.labelThresholdSimilirityPercent.pack(side=LEFT)
				self.frameThresoldSimilarity.pack(side = TOP, fill=X)
				self.similarityButton = ttk.Button(self.frameSearchButton, text = 'Similarity', command = lambda : self.similarity())
				self.similarityButton.pack(side = TOP, fill=X)
			else:
				self.searchButton = ttk.Button(self.frameSearchButton, text = 'Search', command = lambda : self.matrixForGene(self.colorList))
				self.searchButton.pack(side = TOP, fill=X)
			self.reset = ttk.Button(self.frameSearchButton, text = 'Reset', command = lambda : self.resetMatrix())
			self.reset.pack(side = TOP, fill=X)
			if self.gradient == 'patternGradient':
				self.listGeneGO = []
				self.frameInfoPattern = tk.Frame(self.frameSearchButton)
				self.infoText = Label(self.frameInfoPattern, text='Patterns Informations :\nNumber of gexel | Gene')
				self.infoText.pack(fill=X)
				self.listBoxPatterns = Listbox(self.frameInfoPattern, height=10, width=15)
				self.fillListboxSort(self.listBoxPatterns, self.runObject.dicSample[self.nameSample].infoPattern)
				self.listBoxPatterns.bind('<Double-Button>', lambda event, listBox = self.listBoxPatterns, index = 1 : self.selection_gene(event, listBox, index))
				self.listBoxPatterns.bind('<Return>', lambda event, listBox = self.listBoxPatterns, index = 1 : self.selection_gene(event, listBox, index))
				self.listBoxPatterns.pack(fill=X)
				self.GOButton = ttk.Button(self.frameInfoPattern, text = 'Gene Ontology', command = lambda : self.geneOntology(self.listGeneGO, self.listGOBank, self.temp_title))
				self.GOButton.pack(side=BOTTOM, fill = X, pady=10)
				self.frameInfoPattern.pack(side = BOTTOM, fill=X)
			else:
				self.frameThresold = tk.Frame(self.frameSearchButton)
				self.labelThresholdMinMax = Label(self.frameThresold, text='Threshold :')
				self.labelThresholdMinMax.pack(side = TOP, fill=X)
				self.frameThresoldmin = tk.Frame(self.frameThresold)
				self.labelThresholdMin = Label(self.frameThresoldmin, text='Min :')
				self.labelThresholdMin.pack(side = LEFT)
				self.minTresholdEntry = Entry(self.frameThresoldmin, textvar=self.thresholdMin, width=5)
				self.minTresholdEntry.pack(side = LEFT)
				self.frameThresoldmin.pack(side=LEFT)
				self.frameThresoldmax = tk.Frame(self.frameThresold)
				self.labelThresholdMax = Label(self.frameThresoldmax, text='Max :')	
				self.labelThresholdMax.pack(side = LEFT)
				self.maxTresholdEntry = Entry(self.frameThresoldmax, textvar=self.thresholdMax, width=5)
				self.maxTresholdEntry.pack(side = LEFT)
				self.frameThresoldmax.pack(side=RIGHT)
				if self.gradient == 'diffGeneGradient':
					self.frameThresoldDiff = tk.Frame(self.frameSearchButton)
					self.labelThresholdDiff = Label(self.frameThresoldDiff, text='Threshold expression :')
					self.labelThresholdDiff.pack(side = TOP, fill=X)
					self.frameThresoldminDiff = tk.Frame(self.frameThresoldDiff)
					self.labelThresholdDownDiff = Label(self.frameThresoldminDiff, text='Down :')
					self.labelThresholdDownDiff.pack(side = LEFT)
					self.minTresholdEntryDiff = Entry(self.frameThresoldminDiff, textvar=self.thresholdDownDiff, width=5)
					self.minTresholdEntryDiff.pack(side = LEFT)
					self.frameThresoldminDiff.pack(side=LEFT)
					self.frameThresoldmaxDiff = tk.Frame(self.frameThresoldDiff)
					self.labelThresholdUpDiff = Label(self.frameThresoldmaxDiff, text='Up :')	
					self.labelThresholdUpDiff.pack(side = LEFT)
					self.maxTresholdEntryDiff = Entry(self.frameThresoldmaxDiff, textvar=self.thresholdUpDiff, width=5)
					self.maxTresholdEntryDiff.pack(side = LEFT)
					self.frameThresoldmaxDiff.pack(side=RIGHT)
					self.panelDiff = tk.Frame(self.frameMenu)
					self.infoTextDiff = Label(self.panelDiff, text='Genes Informations :')
					self.infoTextDiff.pack(fill=X)
					self.textButtonSwap =StringVar()
					self.textButtonSwap.set('Show repressed genes')
					self.index_hidePanel = 1
					self.swapUpDown = ttk.Button(self.panelDiff, textvariable = self.textButtonSwap, command = lambda : self.hidePanelGene())
					self.swapUpDown.pack(fill=X)
					self.panelGeneUp()
					self.panelDiff.pack(side = BOTTOM, fill=X)
					self.frameThresoldDiff.pack(side = BOTTOM, fill=X)
				elif self.gradient == 'rainbow':
					self.frameLogValue = tk.Frame(self.frameSearchButton)
					self.varLog = IntVar()
					self.checkLog = Checkbutton(self.frameLogValue, variable=self.varLog)
					self.checkLog.deselect()
					self.varLog.trace('w', lambda name, index, mode : self.selectedLog())
					self.checkLog.pack(side=LEFT, padx=10)
					self.labelLog = tk.Label(self.frameLogValue, text=': Log values.')
					self.labelLog.pack(side=LEFT)
					self.frameLogValue.pack(side = BOTTOM, fill=X)
				self.frameThresold.pack(side = BOTTOM, fill=X)
			self.frameSearchButton.pack(side=BOTTOM)
			self.frameSearch.pack(side=BOTTOM, fill=Y)
			self.frameSearchTool.pack(side=TOP)
		self.frameMenu.pack(side = RIGHT, fill=X, anchor = CENTER, expand=1)

	def clipboard(self, event):
		self.master.clipboard_clear()
		self.master.clipboard_append(''.join([item for item in self.listGeneCommunities]))

	def selectedLog(self):
		if self.indexLog == 1:
			self.indexLog = 0
		else:
			self.indexLog = 1

	def panelGeneUp(self):
		self.panelUpDiff = tk.Frame(self.panelDiff)
		self.infoTextDiffUp = Label(self.panelUpDiff, text='Upregulated Genes')
		self.infoTextDiffUp.pack(side = TOP, fill=X)
		self.listBoxUp = Listbox(self.panelUpDiff, height=10, width=10)
		self.fillListboxSort(self.listBoxUp, self.runObject.dicSample[self.nameSample].infoUp)
		self.listBoxUp.bind('<Double-Button>', lambda event, listBox = self.listBoxUp, index = 0 : self.selection_gene(event, listBox, index))
		self.listBoxUp.bind('<Return>', lambda event, listBox = self.listBoxUp, index = 0 : self.selection_gene(event, listBox, index))
		self.listBoxUp.pack(fill=X)
		self.infoLabelUp = Label(self.panelUpDiff, text='Number of gexels | Gene')
		self.infoLabelUp.pack(side = TOP, fill=X)
		self.panelUpDiff.pack(fill=X)

	def panelGeneDown(self):
		self.panelDownDiff = tk.Frame(self.panelDiff)
		self.infoTextDiffDown = Label(self.panelDownDiff, text='Downregulated Genes')
		self.infoTextDiffDown.pack(side = TOP, fill=X)
		self.listBoxDown = Listbox(self.panelDownDiff, height=10, width=10)
		self.fillListboxSort(self.listBoxDown, self.runObject.dicSample[self.nameSample].infoDown)
		self.listBoxDown.bind('<Double-Button>', lambda event, listBox = self.listBoxDown, index = 0 : self.selection_gene(event, listBox, index))
		self.listBoxDown.bind('<Return>', lambda event, listBox = self.listBoxDown, index = 0 : self.selection_gene(event, listBox, index))
		self.listBoxDown.pack(fill=X)
		self.infoLabelDown = Label(self.panelDownDiff, text='Number of gexels | Gene')
		self.infoLabelDown.pack(side = TOP, fill=X)
		self.panelDownDiff.pack(fill=X)

	def hidePanelGene(self):
		if self.index_hidePanel == 0:
			self.textButtonSwap.set('Show repressed genes')
			self.panelDownDiff.pack_forget()
			self.panelGeneUp()
			self.index_hidePanel = 1
		else:
			self.textButtonSwap.set('Show induced genes')
			self.panelUpDiff.pack_forget()
			self.panelGeneDown()		
			self.index_hidePanel = 0


	def selection_gene(self, event, listBox, index):
		"""
		index = 0 --> Diff genes
		index = 1 --> Pattern
		"""
		self.geneSelected = listBox.get(ACTIVE).split(' | ')[1]
		self.searchGene.setEntry(self.geneSelected)
		if index == 0:
			self.matrixForGene(self.colorList)
		else:
			self.patternMatrix()

	def fillListboxSort(self, listBox, listData):
		listData.sort(key=lambda x: x[0], reverse = True)
		for i in listData:
			if i[0] >= self.runObject.minPattern:
				listBox.insert(END, f'{i[0]} | {i[1].upper()}')

	def hideMenu(self):
		"""
		hide the frame menu
		"""
		if self.index_hide == 0:
			self.frameMenu.pack_forget()
			self.index_hide = 1
		else:
			self.menu()
			self.index_hide = 0

	def colorGradient(self):
		"""
		generate colors gradient
		"""
		self.colorScale.delete('all')
		list_color_hex = []
		for x in range(0, 256):
			#Green to red gradient
			if self.gradient == 'diffGeneGradient':
				b = 0
				if x > 128:
					r = (x-128)*2
					g = 0
				else:
					r = 0
					if x != 128:
						g = 255-(x*2)
					else:
						g = 255-((x*2)-1)
			#Rainbow grandient
			else:
				if x < 64:
					r=0
					g = x*2
					b = 255
				elif x >= 64 and x < 128:
					r = 0
					g = x*2
					b = 255 - (((x-64)*4)+3)
				elif x >= 128 and x < 192:
					r = ((x-127)*4)-1
					g = 255 - (x-128)*2
					b = 0
				else:
					r=255
					g = 255 - (x-128)*2
					b = 0
			self.colorScale.create_rectangle(x*2, 0, x*2 + 2, 16, fill=self.rgb(r, g, b), outline=self.rgb(r, g, b))
			self.draw.rectangle([(x*2)+(3*self.indexMiniSizeGexel), 0, (x*2)+(3*self.indexMiniSizeGexel) + 2, 16], fill = self.rgb(r, g, b), outline = self.rgb(r, g, b))
			list_color_hex.append(self.rgb(r, g, b))
		return list_color_hex

	def colors(self):
		"""
		Colors for similarity patterns 
		"""
		self.colorScale.delete('all')
		list_color_hex = [self.rgb(0, 0, 255), self.rgb(75, 0, 230), self.rgb(244, 114, 208), self.rgb(0, 191, 255), self.rgb(50, 205, 50), \
			self.rgb(96, 169, 23), self.rgb(255, 255, 0), self.rgb(240, 163, 10), self.rgb(250, 104, 0), \
			self.rgb(255, 0, 0)]
		x1 = 0
		x2 = 50
		for color in list_color_hex:
			self.colorScale.create_rectangle(x1, 0, x2, 16, fill=color, outline=color)
			self.draw.rectangle([x1+(3*self.indexMiniSizeGexel), 0, x2+(3*self.indexMiniSizeGexel), 16], fill = color, outline = color)
			x1 += 50
			x2 += 50
		self.draw.text((self.indexMiniSizeGexel, 4), '< 0%',(0,0,0))
		self.draw.text((x2+(2*self.indexMiniSizeGexel), 4), '100%',(0,0,0))
		return list_color_hex

	def rgb(self, r, g, b):
		"""
		RGB to hexadecimal
		"""
		return "#%s%s%s" % tuple([hex(color)[2:].rjust(2, "0") for color in (r, g, b)])

	def quit(self):
		"""
		clode window
		"""
		self.master.destroy()

	def similarity(self):
		"""
		Co-Expression Pattern analysis matrix
		"""
		if self.patternDisplayString.get() == '':
			tk.messagebox.showwarning('No pattern selected !', 'You have to select a pattern before performing a similarity analysis.', icon='warning')
		else:
			self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
				(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
			self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
			dic = {}
			listSort = set()
			finalDic = {}
			self.listGeneGO = []
			list_header = []
			list_genes = []
			indexMaxSize = 0
			for i in self.listValue:
				if int(self.comboThresholdSimilarity.get()) <= i:
					list_header.append(i)
					listTemp = set()
					print(f'From {i+10} to {i}% :')
					nameGenePattern = f'{self.geneQuery.get()} | {self.patternDisplayString.get()}'
					listTemp = set(self.runObject.dicSample[self.nameSample].simMatrix[self.runObject.dicSample[self.nameSample].simMatrix[nameGenePattern] > i].index.tolist()).difference(set(listSort))
					listSort.update(listTemp)
					if len(listTemp) != 0:
						print(listTemp)
						if indexMaxSize < len(listTemp):
							indexMaxSize = len(listTemp)
						list_genes.append(list(listTemp))
					else:
						list_genes.append([])
						print('No genes')
					if len(listTemp) != 0:
						for gene in listTemp:
							listCoor = []
							geneQ = gene.split(' | ')[0]
							self.listGeneGO.append(geneQ)
							patternNumber = int(gene.split(' | ')[1])
							listCoor.extend(self.runObject.dicSample[self.nameSample].dicSampleDiff[geneQ].dic_pattern[patternNumber])
							for coor in listCoor:
								if coor not in finalDic.keys():
									if i == 0:
										finalDic[coor] = 0
									else:
										finalDic[coor] = int(i/10)
					print('###################')
			newL = []
			for i in list_genes:
				if not indexMaxSize == len(i):
					i.extend(['']*(indexMaxSize-len(i)))
					newL.append(i)
				else:
					newL.append(i)
			ar =  np.array(newL)
			df = pd.DataFrame(ar, index=list_header)
			self.can.delete('all')
			self.colorList = self.colors()
			self.maxValueGradient.set('100%')
			self.minValueGradient.set('< 0%')
			x1 = 1
			x2 = self.indexMiniSizeGexel
			y1 = 1
			y2 = self.indexMiniSizeGexel
			for y in range(1, self.runObject.sizeY+1):
				x1 = 1
				x2 = self.indexMiniSizeGexel
				for x in range(1, self.runObject.sizeX+1):
					coordinateTemp = f'{x}x{y}'
					if coordinateTemp in self.dicGexel:
						if coordinateTemp in finalDic.keys():
							can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = self.colorList[finalDic[coordinateTemp]], activefill = 'white', activeoutline = 'white')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = self.colorList[finalDic[coordinateTemp]], outline = 'black')
							self.can.tag_bind(can_temp, "<Enter>", lambda event,
								countDisplay = '',
								coordinate = coordinateTemp,
								countDisplayString = self.countDisplayString
								: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
							self.can.tag_bind(can_temp, "<Leave>", lambda event,
								countDisplay = '',
								coordinate = '', 
								countDisplayString = self.countDisplayString 
								: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
						else:
							can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = 'black')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = 'black', outline = 'black')
					else:
						self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'grey')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'grey', outline = 'black')
					x1 += self.indexMiniSizeGexel
					x2 += self.indexMiniSizeGexel
				y1 += self.indexMiniSizeGexel
				y2 += self.indexMiniSizeGexel
			self.scatterFrame = tk.Toplevel(self.master)
			self.displayScatter(self.scatterFrame, df)
			self.scatterFrame.title(f'Result of Co-Expression Pattern analysis for {self.geneQuery.get()} - {self.comboThresholdSimilarity.get()}%')
			self.scatterFrame.minsize(500, 380)
			self.scatterFrame.maxsize(500, 380)

	def displayScatter(self, master, df):
		"""
		Co-Expression Pattern analysis table
		"""
		listCol = [i for i in df.index]
		self.frameScatterWindow = tk.Frame(master)
		self.listBoxScatterResult = ttk.Treeview(self.frameScatterWindow, columns=listCol,
			displaycolumns=listCol, show="headings", height=15)
		for i in listCol:
			self.listBoxScatterResult.heading(i, text=i)
			self.listBoxScatterResult.column(i, width=100)
		self.vsb = ttk.Scrollbar(self.frameScatterWindow, orient="vertical",
			command=self.listBoxScatterResult.yview)
		self.hsb = ttk.Scrollbar(self.frameScatterWindow, orient="horizontal",
			command=self.listBoxScatterResult.xview)
		self.listBoxScatterResult.configure(yscrollcommand=self.vsb.set, xscrollcommand=self.hsb.set)
		self.listBoxScatterResult.grid(column=0, row=0, sticky='nsew', in_=self.frameScatterWindow)
		self.vsb.grid(column=1, row=0, sticky='ns', in_=self.frameScatterWindow)
		self.hsb.grid(column=0, row=1, sticky='ew', in_=self.frameScatterWindow)
		self.frameScatterWindow.grid_columnconfigure(0, weight=1)
		self.frameScatterWindow.grid_rowconfigure(0, weight=1)
		self.frameScatterWindow.pack(fill=X, side=TOP)
		self.frameButtonScatter = tk.Frame(master)
		self.buttonSaveScatterMatrix = ttk.Button(self.frameButtonScatter, text="Save", command= lambda : self.saveMatrixScatter(df.transpose()))
		self.buttonSaveScatterMatrix.pack(side=LEFT)
		self.buttonQuit = ttk.Button(self.frameButtonScatter, text="Quit", command= lambda : master.destroy())
		self.buttonQuit.pack(side=RIGHT)
		self.frameButtonScatter.pack()
		for i in df.columns:
			self.listBoxScatterResult.insert('', 'end', values=list(df[i]))

	def saveMatrixScatter(self, df):
		filename = asksaveasfilename(defaultextension=".tsv", filetypes=(('Tabulation-separated values', '*.tsv'),))
		if filename != '':
			df.to_csv(filename, sep='\t')
	
	def geneOntology(self, listQuery, listGOBank, title):
		if self.gradient == 'patternGradient' and len(self.listGeneGO) == 0:
			tk.messagebox.showwarning('No Co-Expression Patterns', 'Perform a Co-Expression Pattern analysis before asking for GO terms.', icon='warning')
		else:
			self.go_GUI = tk.Toplevel(self.master)
			self.go_GUI.title(f'Gene Ontology Analysis - {title}')
			GoGUI(self.go_GUI, listQuery, listGOBank, title)
			print(listQuery)

	def selectPattern(self, event, indexSelected, gene):
		"""
		When click on pattern, get the gene and index of pattern
		"""
		self.patternDisplayString.set(indexSelected)
		self.geneQuery.set(gene)

	def resetMatrix(self):
		"""
		Reset matrix
		"""
		if self.gradient == 'patternGradient':
			self.colorScale.delete('all')
			self.listGeneGO = []
			self.patternDisplayString.set('')
			self.maxValueGradient.set('')
			self.minValueGradient.set('')
		if self.searchGene.get() != '':
			self.searchGene.erase()
		if self.thresholdUpDiff.get() != 1:
			self.thresholdUpDiff.set(1)
		if self.thresholdDownDiff.get() != -1:
			self.thresholdDownDiff.set(-1)
		if self.indexLog == 1:
			self.checkLog.deselect()
		self.indexLog = 0
		self.thresholdMin.set('')
		self.thresholdMax.set('')
		self.master.title(self.title)
		self.drawingMatrix(self.colorList)

	def drawingMatrix(self, colorList):
		self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
			(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
		self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
		self.can.delete('all')
		if self.indexLog == 1:
			self.master.title(f'{self.title} - Log')
		else:
			self.master.title(self.title)
		x1 = 1
		x2 = self.indexMiniSizeGexel
		y1 = 1
		y2 = self.indexMiniSizeGexel
		if self.gradient != 'patternGradient':
			if self.indexLog == 1:
				if self.thresholdMax.get() == '':
					self.maxValueGradient.set(int(np.log2(self.maxCountGexel)))
				else:
					self.maxValueGradient.set(float(self.thresholdMax.get()))

				if self.thresholdMin.get() == '':
					self.minValueGradient.set(int(np.log2(self.minCountGexel)))
				else:
					self.minValueGradient.set(float(self.thresholdMin.get()))
				
			else:
				if self.thresholdMax.get() == '':
					self.maxValueGradient.set(int(self.maxCountGexel))
				else:
					self.maxValueGradient.set(float(self.thresholdMax.get()))
				

				if self.thresholdMin.get() == '':
					self.minValueGradient.set(int(self.minCountGexel))
				else:
					self.minValueGradient.set(float(self.thresholdMin.get()))
				
		for y in range(1, self.runObject.sizeY+1):
			x1 = 1
			x2 = self.indexMiniSizeGexel
			for x in range(1, self.runObject.sizeX+1):
				coordinateTemp = f'{x}x{y}'
				if self.gradient == 'diffGeneGradient' or self.gradient == 'patternGradient':
					if coordinateTemp in self.dicGexel:
						can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'black')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'black', outline = 'black')
					else:
						self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'grey')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'grey', outline = 'black')
				elif coordinateTemp in self.dicGexel.keys():
					try:
						if self.indexLog == 1:
							if np.log2(self.dicGexel[coordinateTemp].totalGeneCount) <= float(self.minValueGradient.get()) and self.thresholdMin.get() != '':
								indexColor = 'black'
							elif np.log2(self.dicGexel[coordinateTemp].totalGeneCount) >= float(self.maxValueGradient.get()):
								indexColor = 1
							elif (np.log2(self.maxCountGexel)-np.log2(self.minCountGexel)) == 0:
								indexColor = 1
							else:
								indexColor = ((np.log2(self.dicGexel[coordinateTemp].totalGeneCount)-np.log2(self.minCountGexel))/(np.log2(self.maxCountGexel)-np.log2(self.minCountGexel)))
						else:
							if self.dicGexel[coordinateTemp].totalGeneCount <= float(self.minValueGradient.get()) and self.thresholdMin.get() != '':
								indexColor = 'black'
							elif self.dicGexel[coordinateTemp].totalGeneCount >= float(self.maxValueGradient.get()):
								indexColor = 1
							elif (self.maxCountGexel-self.minCountGexel) == 0:
								indexColor = 1
							else:
								indexColor = (((self.dicGexel[coordinateTemp].totalGeneCount)-(self.minCountGexel))/((self.maxCountGexel)-(self.minCountGexel)))
						if indexColor == 0:
							color = colorList[0]
						elif indexColor == 'black':
							color = 'black'
						else:
							color = colorList[math.ceil(indexColor*len(colorList))-1]
					except ValueError:
						color = 'white'
					can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
						fill = color, activefill = 'white', activeoutline = 'white')
					self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
						fill = color, outline = 'black')
					if self.indexLog == 1:
						self.can.tag_bind(can_temp, "<Enter>", lambda event,
							countDisplay = np.log2(self.dicGexel[coordinateTemp].totalGeneCount),
							coordinate = coordinateTemp,
							countDisplayString = self.countDisplayString 
							: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
					else:
						self.can.tag_bind(can_temp, "<Enter>", lambda event,
							countDisplay = self.dicGexel[coordinateTemp].totalGeneCount,
							coordinate = coordinateTemp,
							countDisplayString = self.countDisplayString 
							: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
					self.can.tag_bind(can_temp, "<Leave>", lambda event,
						countDisplay = '',
						coordinate = '', 
						countDisplayString = self.countDisplayString 
						: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
				else:
					self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
						fill = 'grey')
					self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
						fill = 'grey', outline = 'black')
				x1 += self.indexMiniSizeGexel
				x2 += self.indexMiniSizeGexel
			y1 += self.indexMiniSizeGexel
			y2 += self.indexMiniSizeGexel

	def comMatrix(self, value):
		"""
		Communities matrix
		"""
		try:
			self.listBoxCommunities.delete(0, END)
			self.displayNumberPerCommunities.set('')
		except:
			pass
		df = self.currentLouvain
		self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
			(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
		self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
		if value == 'All':
			self.listGeneCommunities = df.drop_duplicates(subset =['TF','Community'], keep = 'first').set_index('TF')['Community'].to_dict()
		else:	
			self.listGeneCommunities = []
		df =  df.drop(['TG'], axis=1)
		df = df.sort_values(by ='Weight', ascending=False)
		dicCom = df.drop_duplicates(subset =['TF','Community'], keep = 'first').set_index('TF')['Community'].to_dict()
		dic_Coordiante_Communities = {}
		for gene in dicCom.keys():
			geneName = gene.split(' | ')[0]
			patternNumber = int(gene.split(' | ')[1])
			for coordinate in self.dicGeneSample[geneName].dic_pattern[patternNumber]:
				if value == 'All':
					if coordinate in dic_Coordiante_Communities.keys():
						pass
					else:	
						dic_Coordiante_Communities[coordinate] = dicCom[gene]
				elif dicCom[gene] == int(value):
					dic_Coordiante_Communities[coordinate] = dicCom[gene]
					if geneName not in self.listBoxCommunities.get(0, END):
						self.listBoxCommunities.insert(END, geneName)
						self.listGeneCommunities.append(geneName)
		if value != 'All':
			self.displayNumberPerCommunities.set(self.listBoxCommunities.size())
		self.temp_title = f'{value} - {self.title}'
		self.master.title(f'{value} - {self.title}')
		listColor = ['green', 'blue', 'red', 'gold', 'purple',
			'pink', 'orange', 'yellow', 'darkgreen', 'darkred', 
			'darkblue', 'limegreen', 'darkorange', 'deeppink', 'cyan']
		x1 = 1
		x2 = self.indexMiniSizeGexel
		y1 = 1
		y2 = self.indexMiniSizeGexel
		if len(dicCom.keys()) == 0:
			tk.messagebox.showwarning('No communities found !', 'We did not detected any community.', icon='warning')
		for y in range(1, self.runObject.sizeY+1):
			x1 = 1
			x2 = self.indexMiniSizeGexel
			for x in range(1, self.runObject.sizeX+1):
				coordinateTemp = f'{x}x{y}'
				if len(dicCom.keys()) != 0:
					if coordinateTemp in self.dicGeneSample[geneName].dic_coordinate_count.keys():
						if coordinateTemp in dic_Coordiante_Communities.keys():
							can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = listColor[((dic_Coordiante_Communities[coordinateTemp]+1)%(len(listColor)))-1], activefill = 'white', activeoutline = 'white')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = listColor[((dic_Coordiante_Communities[coordinateTemp]+1)%(len(listColor)))-1], outline = 'black')
						else:			
							can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = 'black')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = 'black', outline = 'black')
					else:
						can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'grey')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'grey', outline = 'black')
				else:
					if coordinateTemp in self.dicGexel:
						self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'black')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'black', outline = 'black')
					else:
						self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'grey')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'grey', outline = 'black')
				x1 += self.indexMiniSizeGexel
				x2 += self.indexMiniSizeGexel
			y1 += self.indexMiniSizeGexel
			y2 += self.indexMiniSizeGexel

	def patternMatrix(self):
		"""
		Display the selected gene's pattern
		"""
		if self.searchGene.get() in self.dicGeneSample.keys():
			self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
				(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
			self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
			self.colorScale.delete('all')
			self.maxValueGradient.set('')
			self.minValueGradient.set('')
			geneName = self.searchGene.get()
			indexGexelCreate = 0
			self.patternDisplayString.set('')
			self.can.delete('all')
			self.master.title(f'{self.title} - {geneName}')
			self.temp_title = f'{self.title} - {geneName}'
			listColor = ['green', 'blue', 'red', 'gold', 'purple',
				'pink', 'orange', 'yellow', 'darkgreen', 'darkred', 
				'darkblue', 'limegreen', 'darkorange', 'deeppink', 'cyan']
			x1 = 1
			x2 = self.indexMiniSizeGexel
			y1 = 1
			y2 = self.indexMiniSizeGexel
			for y in range(1, self.runObject.sizeY+1):
				x1 = 1
				x2 = self.indexMiniSizeGexel
				for x in range(1, self.runObject.sizeX+1):
					coordinateTemp = f'{x}x{y}'
					if coordinateTemp in self.dicGeneSample[geneName].dic_coordinate_count.keys():
						indexFind = 0
						for index in self.dicGeneSample[geneName].dic_pattern.keys():
							if coordinateTemp in self.dicGeneSample[geneName].dic_pattern[index]:
								can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
									fill = listColor[(index%(len(listColor)))-1], activefill = 'white', activeoutline = 'white')
								self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
									fill = listColor[(index%(len(listColor)))-1], outline = 'black')
								self.can.tag_bind(can_temp, "<Button-1>", lambda event,
									indexSelected = index,
									gene = geneName,
									: self.selectPattern(event, indexSelected, gene))
								self.can.tag_bind(can_temp, "<Enter>", lambda event,
									countDisplay = '',
									coordinate = coordinateTemp,
									countDisplayString = self.countDisplayString
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
								self.can.tag_bind(can_temp, "<Leave>", lambda event,
									countDisplay = '',
									coordinate = '', 
									countDisplayString = self.countDisplayString 
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
								indexFind += 1
								indexGexelCreate += 1
						if indexFind == 0:
							can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = 'black')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = 'black', outline = 'black')
					else:
						can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'grey')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'grey', outline = 'black')
					x1 += self.indexMiniSizeGexel
					x2 += self.indexMiniSizeGexel
				y1 += self.indexMiniSizeGexel
				y2 += self.indexMiniSizeGexel
			if indexGexelCreate == 0:
				tk.messagebox.showwarning('No pattern detected', 'No pattern is detected for the selected gene.', icon='warning')
		
	def matrixForGene(self, colorList):
		"""
		Display a selected gene
		"""
		indexGexelCreate = 0
		self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
			(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
		self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
		if self.searchGene.get() in self.dicGeneSample.keys():
			if self.gradient == 'patternGradient':
				self.colorScale.delete('all')
			if self.indexLog == 1 and self.gradient == 'rainbow':
				self.matrixLog(colorList)
			else:
				geneName = self.searchGene.get()
				self.thresholdDown = 0
				self.thresholdUp = 0
				self.can.delete('all')
				self.master.title(f'{self.title} - {geneName}')
				x1 = 1
				x2 = self.indexMiniSizeGexel
				y1 = 1
				y2 = self.indexMiniSizeGexel
				if self.gradient != 'diffGeneGradient':
					if self.thresholdMax.get() == '':
						self.maxValueGradient.set(int(self.dicGeneSample[geneName].maxCount))
					else:
						self.maxValueGradient.set(float(self.thresholdMax.get()))
					if self.thresholdMin.get() == '':
						self.minValueGradient.set(int(self.dicGeneSample[geneName].minCount))
					else:
						self.minValueGradient.set(float(self.thresholdMin.get()))
				else:
					if self.thresholdUpDiff.get() == '':
						self.thresholdUp = self.runObject.upDiffThreshold
					else:
						self.thresholdUp = float(self.thresholdUpDiff.get())

					if self.thresholdDownDiff.get() == '':
						self.thresholdDown = self.runObject.downDiffThreshold
					else:
						self.thresholdDown = float(self.thresholdDownDiff.get())

					if self.thresholdMax.get() == '':
						self.maxValueGradient.set(int(4))
					else:
						self.maxValueGradient.set(float(self.thresholdMax.get()))

					if self.thresholdMin.get() == '':
						self.minValueGradient.set(int(-4))
					else:
						self.minValueGradient.set(float(self.thresholdMin.get()))
				for y in range(1, self.runObject.sizeY+1):
					x1 = 1
					x2 = self.indexMiniSizeGexel
					for x in range(1, self.runObject.sizeX+1):
						coordinateTemp = f'{x}x{y}'
						if coordinateTemp in self.dicGeneSample[geneName].dic_coordinate_count.keys():
							if self.gradient == 'diffGeneGradient':
								if self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] >= self.maxValueGradient.get():
									indexColor = 255
								elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] <= self.minValueGradient.get():
									indexColor = 0
								elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] < self.maxValueGradient.get() and self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] >= float(self.thresholdUp):
									indexColor = round((float(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])*110/(self.maxValueGradient.get())))+145
								elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] > self.minValueGradient.get() and self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] <= float(self.thresholdDown):
									indexColor = 110-round((float(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])*110/(self.minValueGradient.get())))
								else:
									indexColor = 127
								can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
									fill = colorList[indexColor], activefill = 'white', activeoutline = 'white')
								self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
									fill = colorList[indexColor], outline = 'black')
								self.can.tag_bind(can_temp, "<Enter>", lambda event,
									countDisplay = self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp],
									coordinate = coordinateTemp,
									countDisplayString = self.countDisplayString 
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
								self.can.tag_bind(can_temp, "<Leave>", lambda event,
									countDisplay = '',
									coordinate = '', 
									countDisplayString = self.countDisplayString
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
								indexGexelCreate += 1
							elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] != 0:
								if self.thresholdMin.get() != '' or self.thresholdMax.get() != '':

									if self.thresholdMax.get() == '':
										indexMax = self.dicGeneSample[geneName].maxCount
									else:
										indexMax = self.thresholdMax.get()

									if self.thresholdMin.get() == '':
										indexMin = self.dicGeneSample[geneName].minCount
									else:
										indexMin = self.thresholdMin.get()

									if self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] < float(indexMin):
										color = 'black'
									elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] >= float(indexMax):
										color = colorList[256-1]
									else:
										indexColor = (((self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])-(float(indexMin)))/(float(indexMax)-float(indexMin)))
										if math.ceil(indexColor*len(colorList)) == 0:
											color = colorList[0]
										else:
											color = colorList[math.ceil(indexColor*len(colorList))-1]
								else:
									indexMin = self.dicGeneSample[geneName].minCount
									indexMax = self.dicGeneSample[geneName].maxCount
									indexColor = (((self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])-(float(indexMin)))/(float(indexMax)-float(indexMin)))
									if indexColor == 0:
										color = colorList[0]
									else:
										color = colorList[math.ceil(indexColor*len(colorList))-1]
								can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
									fill = color, activefill = 'white', activeoutline = 'white')
								self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
									fill = color, outline = 'black')
								self.can.tag_bind(can_temp, "<Enter>", lambda event,
									countDisplay = self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp],
									coordinate = coordinateTemp,
									countDisplayString = self.countDisplayString 
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
								self.can.tag_bind(can_temp, "<Leave>", lambda event,
									countDisplay = '',
									coordinate = '', 
									countDisplayString = self.countDisplayString 
									: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
							else:
								self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
									fill = 'black')
								self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
									fill = 'black', outline = 'black')
						else:
							self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
								fill = 'grey')
							self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
								fill = 'grey', outline = 'black')
						x1 += self.indexMiniSizeGexel
						x2 += self.indexMiniSizeGexel
					y1 += self.indexMiniSizeGexel
					y2 += self.indexMiniSizeGexel
		elif self.searchGene.get() =='':
			self.drawingMatrix(self.colorList)
		else:
			if indexGexelCreate == 0 and self.gradient == 'diffGeneGradient':
				tk.messagebox.showwarning('Selected gene is not differentially expressed', 'Selected gene is not differentially expressed.', icon='warning')
			else:
				tk.messagebox.showwarning('Name of gene is not found', 'Please provide a correct gene name.', icon='warning')
	
	def matrixLog(self, colorList):
		"""
		Display matrix with log values
		"""
		self.saveCanvasImage = PIL.Image.new('RGBA', (int(((self.runObject.sizeX+1)*(self.indexMiniSizeGexel))+16),
			(int((self.runObject.sizeY+1)*(self.indexMiniSizeGexel))+16)), (0,0,0,0))
		self.draw = PIL.ImageDraw.Draw(self.saveCanvasImage)
		geneName = self.searchGene.get()
		self.can.delete('all')
		self.master.title(f'{self.title} - {geneName} - Log2')
		x1 = 1
		x2 = self.indexMiniSizeGexel
		y1 = 1
		y2 = self.indexMiniSizeGexel
		if self.thresholdMax.get() == '':
			if self.dicGeneSample[geneName].maxCount == 0:
				self.maxValueGradient.set(0)
			else:
				self.maxValueGradient.set(round(np.log2(self.dicGeneSample[geneName].maxCount)))
		elif self.thresholdMax.get() == 0:
			self.maxValueGradient.set(0)
		else:
			self.maxValueGradient.set(float(self.thresholdMax.get()))

		if self.thresholdMin.get() == '':
			if self.dicGeneSample[geneName].minCount == 0:
				self.minValueGradient.set(0)
			else:
				self.minValueGradient.set(round(np.log2(self.dicGeneSample[geneName].minCount)))
		elif self.thresholdMin.get() == 0:
			self.minValueGradient.set(0)
		else:
			self.minValueGradient.set(float(self.thresholdMin.get()))
		for y in range(1, self.runObject.sizeY+1):
			x1 = 1
			x2 = self.indexMiniSizeGexel
			for x in range(1, self.runObject.sizeX+1):
				coordinateTemp = f'{x}x{y}'
				if coordinateTemp in self.dicGeneSample[geneName].dic_coordinate_count.keys():						
					if self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] != 0:
						if self.thresholdMin.get() != '' or self.thresholdMax.get() != '':
							if self.thresholdMax.get() == '':
								indexMax = np.log2(self.dicGeneSample[geneName].maxCount)
							else:
								indexMax = self.thresholdMax.get()
							if self.thresholdMin.get() == '':
								if self.dicGeneSample[geneName].minCount == 0:
									indexMin = self.dicGeneSample[geneName].minCount
								else:
									indexMin = np.log2(self.dicGeneSample[geneName].minCount)
							else:
								indexMin = self.thresholdMin.get()
							if np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp]) < float(indexMin):
								color = 'black'

							elif np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp]) >= float(indexMax):
								color = colorList[255]
							else:
								indexColor = ((np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])-float(indexMin))/(float(indexMax)-float(indexMin)))
								if indexColor == 0:
									color = colorList[0]
								else:
									color = colorList[math.ceil(indexColor*len(colorList))-1]
						else:
							indexMin = self.dicGeneSample[geneName].minCount
							indexMax = self.dicGeneSample[geneName].maxCount
							if self.dicGeneSample[geneName].minCount == 0: 
								indexMini = 1
							else: 
								indexMini = np.log2(self.dicGeneSample[geneName].minCount)
							if (np.log2(self.dicGeneSample[geneName].maxCount)-indexMini) == 0:
								indexColor = 1
							elif indexMin == 0:
								indexColor = np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])/np.log2(indexMax)
							elif self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp] == 0:
								indexColor = 'black'
							else:
								indexColor = (np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp])-np.log2(indexMin))/(np.log2(indexMax)-np.log2(indexMin))
							if indexColor == 'black':
								color = 'black'
							elif indexColor == 0 or indexColor < 0:
								color = colorList[0]
							else:
								color = colorList[math.ceil(indexColor*len(colorList))-1]
						can_temp = self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = color, activefill = 'white', activeoutline = 'white')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = color, outline = 'black')
						self.can.tag_bind(can_temp, "<Enter>", lambda event,
							countDisplay = np.log2(self.dicGeneSample[geneName].dic_coordinate_count[coordinateTemp]),
							coordinate = coordinateTemp,
							countDisplayString = self.countDisplayString 
							: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
						self.can.tag_bind(can_temp, "<Leave>", lambda event,
							countDisplay = '',
							coordinate = '', 
							countDisplayString = self.countDisplayString 
							: self.gexel_active(event, countDisplay, coordinate, countDisplayString))
					else:
						self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
							fill = 'black')
						self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
							fill = 'black', outline = 'black')
				else:
					self.can.create_rectangle(x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel,
						fill = 'grey')
					self.draw.rectangle([x1+self.indexMiniSizeGexel, y1+self.indexMiniSizeGexel, x2+self.indexMiniSizeGexel, y2+self.indexMiniSizeGexel],
						fill = 'grey', outline = 'black')
				x1 += self.indexMiniSizeGexel
				x2 += self.indexMiniSizeGexel
			y1 += self.indexMiniSizeGexel
			y2 += self.indexMiniSizeGexel

	def gexel_active(self, event, countDisplay, coordinate, countDisplayString):
		"""
		When mouse over a gexel
		"""
		self.countDisplayStringCoor.set(coordinate)
		try:
			countDisplayString.set(round(countDisplay, 2))
		except:
			countDisplayString.set(countDisplay)

def main():
	run = Run()
	runGUI = StartGUI(run)
	#runGUI.iconbitmap(os.path.join('logo', 'multilaye_ico.ico'))
	runGUI.mainloop()
	if run.indexRun == 1:
		root = tk.Tk()
		if run.mode == 'explore':
			#print(f'X {run.sizeX} X cal {run.calculatedX}')
			#print(f'Y {run.sizeY} Y cal {run.calculatedY}')
			if run.sizeX < run.calculatedX or run.sizeY < run.calculatedY:
				tempX = 0
				tempY = 0
				if run.sizeX < run.calculatedX:
					tempX = run.calculatedX
				else:
					tempX = run.sizeX
				if run.sizeY < run.calculatedY:
					tempY = run.calculatedY
				else:
					tempY = run.sizeY
				warningGUI = tk.messagebox.askquestion('Size warning', f'At least 1 coordinate is out of the given map size {run.sizeX}x{run.sizeY}.\nDo you want to change the current size {run.sizeX}x{run.sizeY} by {tempX}x{tempY} for the display ?', icon='warning')
				if warningGUI == 'yes':
					if run.sizeX < run.calculatedX:
						run.sizeX = run.calculatedX
					if run.sizeY < run.calculatedY:
						run.sizeY = run.calculatedY
					run.totalGexel = run.sizeX * run.sizeY
		app = mainGUI(root, run) 
		#root.iconbitmap(os.path.join('logo', 'multilaye_ico.ico'))
		root.mainloop()

if __name__ == "__main__":
	main()