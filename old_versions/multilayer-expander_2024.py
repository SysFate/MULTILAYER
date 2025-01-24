import argparse
from xmlrpc.client import boolean
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from time import time
import warnings

e=None


def quincux_detection(dataframe):
	dict_min_max = min_max_detection(dataframe)
	y = min(dict_min_max.keys())
	x = dict_min_max[y][0]
	if f"{int(x)+1}x{int(y)}" not in dataframe.columns and f"{int(x)}x{int(y)+1}" not in dataframe.columns:
		print('Quincux = True')
		return True
	print('Quincux = False')
	return False

def first_gexel_parity(dataframe):
	dict_min_max = min_max_detection(dataframe)
	ymin = min(dict_min_max.keys())
	if ymin%2==0:
		return 0
	return 1


def min_max_detection(dataframe): # dictionnaire coordonnée_y -> [min_x, max_x] 
	colnames=list(dataframe.columns)
	colnames=[[int(col.split('x')[0]), int(col.split('x')[1])] for col in colnames]
	dict_min_max = {}
	for index in colnames:
		if index[1] not in dict_min_max:
			all_x = [col[0] for col in colnames if col[1]==index[1]]
			dict_min_max[index[1]]=[min(all_x), max(all_x)]
	return dict_min_max


def fill_interstices(dataframe,parity):
	colnames=list(dataframe.columns)
	max_x=max([int(index.split('x')[0]) for index in colnames])
	#max_y=max([int(index.split('x')[1]) for index in colnames])
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if x!=max_x:
			dataframe[f"{x+1}x{y}"]=calc_mean_straight(dataframe, x+1, y)
	smooth_border(dataframe,parity)

def smooth_border(dataframe,parity):
	colnames=list(dataframe.columns)
	max_x=max([int(index.split('x')[0]) for index in colnames])
	min_x=min([int(index.split('x')[0]) for index in colnames])
	dict_min_max = min_max_detection(dataframe)
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if x!=min_x and y%2==parity and dict_min_max[y][0]==x:
			dataframe[f"{x-1}x{y}"]=calc_mean_straight(dataframe,x-1,y)
		if x!=max_x and y%2!=parity and dict_min_max[y][1]==x:
			dataframe[f"{x+1}x{y}"]=calc_mean_straight(dataframe,x+1,y)


def calc_mean_straight(dataframe, x, y):
	global e
	columns = [f"{x-1}x{y}", f"{x+1}x{y}", f"{x}x{y+1}", f"{x}x{y+1}"]
	columns = [n for n in columns if n in dataframe.columns]
	if e==True:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/4
	else:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/len(columns)
	return counts


def calc_mean_diag(dataframe, x, y):
	global e
	columns = [f"{x-1}x{y-1}", f"{x+1}x{y-1}", f"{x-1}x{y+1}", f"{x+1}x{y+1}"]
	columns = [n for n in columns if n in dataframe.columns]
	if e==True:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/4
	else:
		counts = dataframe[columns].sum(axis=1, numeric_only=True)/len(columns)
	return counts
	

def add_interstices(dataframe):
	colnames=list(dataframe.columns)
	new_colnames = []
	for colname in colnames:
		x, y = [int(index) for index in colname.split('x')]
		if x != 1:
			x = x*2 - 1
		if y !=1:
			y = y*2 - 1
		new_colnames.append(f"{x}x{y}")
	max_x=max([int(index.split('x')[0]) for index in new_colnames])
	max_y=max([int(index.split('x')[1]) for index in new_colnames])
	dataframe.columns = new_colnames
	for colname in new_colnames:
		x, y = [int(index) for index in colname.split('x')]
		if x!=max_x and y!=max_y:
			dataframe[f"{x+1}x{y+1}"]=calc_mean_diag(dataframe, x+1, y+1)
	
def save(dataframe):
	fileDir = asksaveasfilename(filetypes=[('Tabulation-separated values', '*.tsv')])
	dataframe.to_csv(fileDir, sep='\t', index=True)


def main():
	global e
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-n","--n_iter", metavar='n_iter', type=int, default=1, help="number of iterations")
	parser.add_argument("-e","--empty_as_zero", metavar='empty_as_zero', type=boolean, default=True, help="Count empty pixels as zeros")
	args = parser.parse_args()
	fileDir = askopenfilename(multiple=False, filetypes=[('Tabulation-separated values','*.tsv')])
	begin = time()
	dataframe = pd.read_csv(fileDir, sep='\t', header=0, index_col=0)
	n = args.n_iter
	e = args.empty_as_zero
	parity = first_gexel_parity(dataframe)
	is_quincux = quincux_detection(dataframe)
	warnings.filterwarnings("ignore")
	if is_quincux:
		while n>0:
			fill_interstices(dataframe, parity)
			n-=1
			if n==0:
				break
			add_interstices(dataframe)
			n-=1
	else:
		while n>0:
			add_interstices(dataframe)
			n-=1
			if n==0:
				break
			fill_interstices(dataframe,parity)
			n-=1
	print(f"Time : {time()-begin} sec")
	save(dataframe)
	
if __name__ == "__main__":
	main()
