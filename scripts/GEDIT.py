#!/usr/bin/python
import sys
import getSigGenesModal
import MatrixTools
import HandleInput
import os
import random
import numpy as np

def main():
  """
  usage default:
  python ThisScript.py -mix SamplesMat.tsv -ref RefMat.tsv

  user selected parameters:
  python ThisScript.py -mix SamplesMat.tsv -ref RefMat.tsv -numSigs SigsPerCT -method SigMethod -RS rowscaling
  """
  #where to write results 
  scratchSpace = os.path.join(os.getcwd(), "scratch","")
  try: 
    os.mkdir(scratchSpace) 
  except OSError as error: 
    print(error)  


  myArgs = HandleInput.checkInputs(sys.argv[1:])
  if myArgs[0] == False:
    print (myArgs[1:])
    return

  rawMix = myArgs[0]
  rawRef = myArgs[1]
  SigsPerCT = myArgs[2]
  #minSigs = myArgs[3]
  SigMethod = myArgs[4]
  RowScaling = myArgs[5]
  MixFName = myArgs[6].split("/")[-1]
  RefFName = myArgs[7].split("/")[-1]
  outFile = myArgs[8]

  numCTs = len(rawRef[0])-1
  TotalSigs = int(SigsPerCT*numCTs)
 
  stringParams = [str(m) for m in \
          [MixFName,RefFName,SigsPerCT,SigMethod,RowScaling]]
  
  scratchSpace = scratchSpace + "_".join(stringParams) + "_"

  SampleNames = rawMix[0]
  CTNames = rawRef[0]

  betRef = MatrixTools.remove0s(rawRef)
  
  normMix, normRef = MatrixTools.qNormMatrices(rawMix,betRef)
  sharedMix, sharedRef = MatrixTools.getSharedRows(normMix,betRef)

  if len(sharedMix) < 1 or len(sharedRef) < 1:
      print ("error: no gene names match between reference and mixture")
      return
  if len(sharedMix) < numCTs or len(sharedRef) < numCTs:
      print ("warning: only ", len(sharedMix) , " gene names match between reference and mixture")

  #write normalized matrices
  MatrixTools.writeMatrix([CTNames] + normRef, scratchSpace + "NormRef.tsv") 
  MatrixTools.writeMatrix([SampleNames] + normMix, scratchSpace + "NormMix.tsv")

  SigRef = getSigGenesModal.returnSigMatrix([CTNames] + sharedRef, \
  SigsPerCT, TotalSigs, SigMethod)
  
  SigMix, SigRef = MatrixTools.getSharedRows(sharedMix, SigRef)

  """
  write matrices with only sig genes. files are not used by this program,
  but potentially informative to the user
  """
  
  MatrixTools.writeMatrix([CTNames] + SigRef, scratchSpace + "SigRef.tsv") 
  ScaledRef, ScaledMix = MatrixTools.RescaleRows(SigRef[1:], SigMix[1:], RowScaling)
 
  ScaledRef = [CTNames] + ScaledRef
  ScaledMix = [SampleNames] + ScaledMix 

  refFile = scratchSpace + "ScaledRef.tsv"
  mixFile = scratchSpace + "ScaledMix.tsv"

  MatrixTools.writeMatrix(ScaledRef, refFile)
  MatrixTools.writeMatrix(ScaledMix, mixFile)
   
  print (mixFile)
  print(refFile)
  #predictions = Regression(scratchSpace, CTNames, SampleNames,refFile, mixFile,strDescr)
  #if predictions == False:
  #   return

  #for line in predictions:
  #  print "\t".join([str(el) for el in line])
  return


main()
