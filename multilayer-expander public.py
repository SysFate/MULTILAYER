
import argparse
from xmlrpc.client import boolean
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from time import time
import warnings

e=None
m=None
d=None

def quincux_detection(dataframe):
	colnames=list(dataframe.columns)
	for colname in colnames:
		x, y = colname.split('x')
		if f"{int(x)+1}x{int(y)}" in dataframe.columns or f"{int(x)}x{int(y)+1}" in dataframe.columns:
			print('Quincux = False')
			return False
	print('Quincux = True')
	return True	

def first_gexel_parity(dataframe):
	dict_min_max = dict_min_max_y_x(dataframe)
	ymin = min(dict_min_max.keys())
	if ymin%2==0:
		return 0
	return 1

def dict_min_max_y_x(dataframe): # dictionnaire coordonnée_y -> [min_x, max_x] 
	colnames=list(dataframe.columns)
	colnames=[[int(col.split('x')[0]), int(col.split('x')[1])] for col in colnames]
	dict_min_max = {}
	for index in colnames:
		if index[1] not in dict_min_max:
			all_x = [col[0] for col in colnames if col[1]==index[1]]
			dict_min_max[index[1]]=[min(all_x), max(all_x)]
	return dict_min_max

def dict_min_max_x_y(dataframe): # dictionnaire coordonnée_x -> [min_y, max_y] 
	colnames=list(dataframe.columns)
	colnames=[[int(col.split('x')[0]), int(col.split('x')[1])] for col in colnames]
	dict_min_max = {}
	for index in colnames:
		if index[0] not in dict_min_max:
			all_y = [col[1] for col in colnames if col[0]==index[0]]
			dict_min_max[index[0]]=[min(all_y), max(all_y)]
	return dict_min_max

def fill_interstices(dataframe,parity):
	colnames=list(dataframe.columns)
	dict_min_max = dict_min_max_x_y(dataframe)
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if y!=dict_min_max[x][0]: #y n'est pas le minimum
			dataframe[f"{x}x{y-1}"]=calc_mean_straight(dataframe, x, y-1)
	smooth_border(dataframe,dict_min_max)
	if m:
		mod_true_straight(dataframe, colnames)

def smooth_border(dataframe,dict_min_max):
	colnames=list(dataframe.columns)
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if f"{x+2}x{dict_min_max[x][0]}" in colnames and (y==dict_min_max[x+2][0] or y==dict_min_max[x][0]) and f"{x+1}x{dict_min_max[x][0]}" not in colnames:
			dataframe[f"{x+1}x{dict_min_max[x][0]}"]=calc_mean_straight(dataframe, x+1, dict_min_max[x][0])
		if f"{x+2}x{dict_min_max[x][1]}" in colnames and (y==dict_min_max[x+2][1] or y==dict_min_max[x][1]) and f"{x+1}x{dict_min_max[x][1]}" not in colnames:
			dataframe[f"{x+1}x{dict_min_max[x][1]}"]=calc_mean_straight(dataframe, x+1, dict_min_max[x][1])


def calc_mean_straight(dataframe, x, y):
	global e
	columns = [f"{x-1}x{y}", f"{x+1}x{y}", f"{x}x{y+1}", f"{x}x{y+1}"]
	columns = [n for n in columns if n in dataframe.columns]
	if e:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/4
	else:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/len(columns)
	return counts

def calc_mean_diag(dataframe, x, y):
	global e
	columns = [f"{x-1}x{y-1}", f"{x+1}x{y-1}", f"{x-1}x{y+1}", f"{x+1}x{y+1}"]
	columns = [n for n in columns if n in dataframe.columns]
	if e:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/4
	else:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/len(columns)
	return counts

def mod_true_diag(dataframe, colnames):
	global e
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		columns = [f"{x}x{y}", f"{x-1}x{y-1}", f"{x+1}x{y-1}", f"{x-1}x{y+1}", f"{x+1}x{y+1}"]
		columns = [n for n in columns if n in dataframe.columns]
		if e:
			counts = round(dataframe[columns].sum(axis=1, numeric_only=True)/5, 2)
		else:
			counts = round(dataframe[columns].sum(axis=1, numeric_only=True)/len(columns), 2)
		dataframe[colname]=counts

def mod_true_straight(dataframe, colnames):
	global e
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		columns = [f"{x}x{y}", f"{x-1}x{y}", f"{x+1}x{y}", f"{x}x{y+1}", f"{x}x{y-1}"]
		columns = [n for n in columns if n in dataframe.columns]
		if e:
			counts = round(dataframe[columns].sum(axis=1, numeric_only=True)/5, 2)
		else:
			counts = round(dataframe[columns].sum(axis=1, numeric_only=True)/len(columns), 2)
		dataframe[colname]=counts
	
def delete_true(dataframe,colnames):
	dataframe.drop(labels = colnames, axis = 1, inplace = True, errors = "ignore")

def connect_solitary(dataframe): #TODO improve
	colnames=list(dataframe.columns)
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if f'{x+1}x{y+1}' not in colnames and f'{x-1}x{y+1}' and f'{x+1}x{y-1}' and f'{x-1}x{y-1}':
			dataframe[f"{x+1}x{y+1}"]=calc_mean_diag(dataframe, x+1, y+1)
	
def add_interstices(dataframe):
	global d
	global m
	colnames=list(dataframe.columns)
	new_colnames = []
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		x = x*2 - 1
		y = y*2 - 1
		new_colnames.append(f"{x}x{y}")
	dataframe.columns = new_colnames
	dict_min_max = dict_min_max_y_x(dataframe)
	maxY = max(dict_min_max.keys())
	for colname in new_colnames:
		x, y = [int(index) for index in colname.split('x')]
		if x!=dict_min_max[y][1] and y!=maxY:
			dataframe[f"{x+1}x{y+1}"]=calc_mean_diag(dataframe, x+1, y+1)
	if d:
		delete_true(dataframe, new_colnames)
	elif m:
		mod_true_diag(dataframe, new_colnames)
	

def save(dataframe):
	fileDir = asksaveasfilename(filetypes=[('Tabulation-separated values', '*.tsv')])
	dataframe.to_csv(fileDir, sep='\t', index=True)

def main():
	global e
	global m
	global d
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-n","--n_iter", metavar='n_iter', type=int, default=1, help="number of iterations")
	parser.add_argument("-t","--transpose", metavar='transpose', type=boolean, default=False, help="transpose matrix")
	args = parser.parse_args()
	fileDir = askopenfilename(multiple=False, filetypes=[('Tabulation-separated values','*.tsv')])
	begin = time()
	if args.transpose:
		dataframe = pd.read_csv(fileDir, sep='\t', header=0, index_col=0).transpose()
	else:
		dataframe = pd.read_csv(fileDir, sep='\t', header=0, index_col=0)
	print(f"Shape : {dataframe.shape}")
	n = args.n_iter
	e = True
	m = False
	d = False
	parity = first_gexel_parity(dataframe)
	is_quincux = quincux_detection(dataframe)
	warnings.filterwarnings("ignore")
	if is_quincux:
		while n>0:
			print(f"Fill_interstices : {time()-begin} sec")
			fill_interstices(dataframe, parity)
			print(f"Shape : {dataframe.shape}")
			n-=1
			if n==0:
				break
			print(f"Add_interstices : {time()-begin} sec")
			add_interstices(dataframe)
			print(f"Shape : {dataframe.shape}")
			n-=1
	else:
		while n>0:
			print(f"Add_interstices : {time()-begin} sec")
			add_interstices(dataframe)
			print(f"Shape : {dataframe.shape}")
			n-=1
			if n==0:
				break
			print(f"Fill_interstices : {time()-begin} sec")
			fill_interstices(dataframe,parity)
			print(f"Shape : {dataframe.shape}")
			n-=1
	print(f"Save : {time()-begin} sec")
	save(dataframe)
	print(f"Time : {time()-begin} sec")
	
if __name__ == "__main__":
	main()
