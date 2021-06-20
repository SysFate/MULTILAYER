'''
Created on 26 mai 2021

@author: francois

'''

import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script will convert a raw csv matrix from visium to tsv matrix for multilayer', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m","--matrix", metavar='matrix', type=str, required=True, help="Path to the raw csv matrix")
    parser.add_argument("-p","--position", metavar='position', type=str, required=True, help="Path to the spacial matrix")
    parser.add_argument("-o","--out", metavar='out', type=str, required=True, help="Path where the out file are saved")
    parser.add_argument("-g","--gene", metavar='gene', type=str, required=True, help="Path to genes features")
    parser.add_argument("--compressor", dest='compressor', default=False, action='store_true', help="If you want the matrix in multilayer compressor format")
    args = parser.parse_args()
    
    matrix = pandas.read_csv(args.matrix,index_col=0,header=0)
    posDict = {}
    with open(args.position) as posFile:
        for line in posFile.readlines():
            cutedLine = line.split(",")
            posDict[cutedLine[0]] = str(int(cutedLine[2])+1)+"x"+str(int(cutedLine[3])+1)
    geneDict = {}
    if args.gene[-3:] == ".gz":
        import gzip
        with gzip.open(args.gene,'rt') as geneFile:
            for line in geneFile:
                cutedLine = line.split("\t")
                geneDict[cutedLine[0]] = cutedLine[1]+'_'+cutedLine[0]
    else:
        with open(args.gene) as geneFile:
            for line in geneFile.readlines():
                cutedLine = line.split("\t")
                geneDict[cutedLine[0]] = cutedLine[1]+'_'+cutedLine[0]
    newMatrix = matrix.rename(columns = posDict, index = geneDict)
    newMatrix = newMatrix[~newMatrix.index.duplicated(keep='first')]
    newMatrix.to_csv(args.out,sep="\t")
    if args.compressor:
        compressorVersion = pandas.melt(newMatrix.assign(index=newMatrix.index), id_vars=['index'])
        compressorVersion = compressorVersion[compressorVersion.value != 0]
        compressorVersion.columns = ['gene', 'bc', 'count']
        compressorVersion = compressorVersion[['bc', 'gene', 'count']]
        compressorVersion.to_csv(args.out+".compressor",sep="\t",index=False)