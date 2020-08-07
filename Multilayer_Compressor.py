__author__ = 'Julien Moehlin'

import argparse
import math
import pandas as pd
import time
import numpy as np

def parseArguments():
	parser = argparse.ArgumentParser(description='Multilayer Compressor provided by SysFate Team.')
	parser.add_argument('-i', help='Input matrix with : - bacordes (XxY) with header \'bc\'. - Genes with header \'gene\'. - Genes counts with header \'count\'', required=True)
	parser.add_argument('-o', help='Output matrix formatted for Multilayer tool.', required=True)
	parser.add_argument('-cx', help='Compressor factor.', required=True)
	parser.add_argument('-cy', help='Compressor factor for Y (if cy is not assign, cx will be apply to x and y).')
	args = parser.parse_args()
	return args

def compressor():
	args = 	parseArguments()
	inputFile = args.i
	outputFile = args.o
	compression_x = int(args.cx)
	if args.cy == None:
		compression_y = compression_x
	else:
		compression_y = args.cy
	startTime = time.time()
	df = pd.read_csv(inputFile, sep='\t', index_col=None)
	dic = {}
	l = []
	xMax = 0
	yMax = 0
	for index, row in df.iterrows():
		x, y = row['bc'].split('x')
		coordinate = f'{math.ceil((float(x))/compression_x)}x{math.ceil((float(y))/compression_y)}'
		if math.ceil((float(x))/compression_x) > xMax:
			xMax = math.ceil((float(x))/compression_x)
		if math.ceil((float(y))/compression_y) > yMax:
			yMax = math.ceil((float(y))/compression_y)
		if coordinate not in dic.keys():
			dicGexel = {}
			dicGexel[row['gene'].upper()] = float(row['count'])
			l.append(str(row['gene']).upper())
			dic[coordinate] = dicGexel
		else:
			if str(row['gene']).upper() in dic[coordinate].keys():
				dic[coordinate][str(row['gene']).upper()] = float(dic[coordinate][str(row['gene']).upper()])+int(row['count'])
				l.append(str(row['gene']).upper())
			else:
				dic[coordinate][str(row['gene']).upper()] = float(row['count'])
				l.append(str(row['gene']).upper())
	stopTime = time.time()
	l = set(l)
	ll = []
	for i in dic.keys():
		lll = []
		for o in l:
			if o in dic[i].keys():
				lll.append(dic[i][o])
			else:
				lll.append(0)
		ll.append(lll)
	print(f'Time {stopTime-startTime} secondes')
	ar = np.array(ll)
	dfT = pd.DataFrame(ar, index=dic.keys(), columns=l)
	print('######')
	print('Write compressed Matrix')
	dfT.transpose().to_csv(outputFile, sep='\t')
	print(f'x max : {xMax}')
	print(f'y max : {yMax}')

if __name__ == '__main__':
	compressor()