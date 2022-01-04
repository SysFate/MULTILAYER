__author__ = 'Julien Moehlin'

import argparse
import pandas as pd
import os

def parseArguments():
	parser = argparse.ArgumentParser(description='Enrichr converter provided by SysFate Team.')
	parser.add_argument('-i', help='Input : enrichr librairy.', required=True)
	#parser.add_argument('-o', help='Output : GO DB for Multilayer.', required=True)
	args = parser.parse_args()
	return args

def enrichr_converter():
	args = parseArguments()
	#outputFile = args.o
	outputFile = f'{args.i[:-4]}_converted.tsv'
	read_file = pd.read_csv(args.i, na_values='True')
	inputFile = f'{args.i[:-3]}tsv'
	read_file.to_csv(inputFile, index=None)
	df = pd.read_csv(inputFile, index_col=None, header=None)
	list_col = []
	list_data = []
	list_col.append('0')
	list_data.append('Gene Signature')
	for index, row in df.iterrows():
		filter_list = filter(lambda x: x != '', row.str.split('\t')[0][2:])
		final_list = list(filter_list)
		list_col.append(row.str.split('\t')[0][0])
		list_data.append(final_list)
	final_df = pd.DataFrame(list_data, index=list_col)
	final_df.to_csv(outputFile, sep='\t', header=None)
	os.remove(inputFile)
	print(f'{args.i} has been converted')

if __name__ == '__main__':
	enrichr_converter()