#load files
datafile = './data/Xdata.txt'
infofile = './data/XvarInfo.xlsx'
addfile = './data/posAdductsTable_RP.xlsx'
anno = 'annotClass'

#other datafile
datafile = './data/MESA_Data_LPOS.csv'
infofile = './data/MESA_VarInfo_LPOS.xlsx'
addfile = './data/posAdductsTable_RP.xlsx'
anno = 'RefMet_Sub_class'

datafile = './data/Airwave1_Data_LPOS.csv'
infofile = './data/Airwave1_VarInfo_LPOS.xlsx'
addfile = './data/posAdductsTable_RP.xlsx'
anno = 'RefMet_Sub_class'

datafile = './data/Rotterdam_Data_LPOS.csv'
infofile = './data/Rotterdam_VarInfo_LPOS.xlsx'
addfile = './data/posAdductsTable_RP.xlsx'
anno = 'RefMet_Sub_class'

#load variables
'''
rtThresh = 1/200
corrThresh = 0.75
mzThresh = 0.002235
'''
rtThresh = 1/30
corrThresh = 0.75
mzThresh = 0.002235

ratrt = 0.5
ratmz = 0.5
sdnum = 1.96
corrThreshadd = 0.7
mzCH = 28.0313
rtCH = 1.5
cluster = 3

#load packages

import bdnw
import numpy as np

#load information
dataMat = load.loaddata(datafile)
corrMat = np.corrcoef(dataMat.T)#1848 * 1848


#draw.dr(finIso)
#isotopes

#find all isotope pairs
finIso = bdnw.findIso(dataMat, corrMat, infofile, rtThresh, corrThresh, mzThresh, ratrt, ratmz)

#bdnw.drawhist(finIso, infofile, corrMat)
#visual the network
#bdnw.ig(finIso)

#get connected compounents list
conList = bdnw.getConList(finIso)
#output ./network without label.png
#bdnw.drawIsoNet(finIso)



#output ./edgeInfo.xls
edgeInfo = bdnw.edInfo(infofile, finIso, conList)

#get information for scatter plot
scat, resList, tarList = bdnw.infoOfScat(conList, infofile)


#output ./scatComp.png
#bdnw.drawScatOfPair(scat)
#lineInfo = bdnw.drawScatOfPair_line(scat)

lineInfo = bdnw.drawScatOfPair_intercept0(scat)
#get new prediction for scatter plot (according to standard deviation)
new = bdnw.findWrong(scat, lineInfo)
#bdnw.dishist(withOutlier, lineInfo)
#output ./noWrong.png
#bdnw.drawNoWrong(new)
#wrongP: MZRT list of wrongList
#wrongP = findwrong.wrongPoint(infofile, wrongList)
#correctIso: list of each isotope type, finIso[0]: list of +0 dictionary
#correctIso[0][0]: key:resource isotope MZRT, value: targetMZRT & pair-type & correlation
#correctIso[0][0]: SLPOS_1001.7229_10.8676:SLPOS_1002.7260_10.8665&0_1&0.9055846806697808
correctIso = bdnw.isoAll(new, resList, tarList, infofile, corrMat)

#sdList = [sd, lineA, lineB, mean]
#newList: corrected resource isotope type of each isotope pair(no outlier)
#accurate: accurate of prediction, outP: amount of outliers
#withOutlier = bdnw.findOut_sd(new, sdnum)

withOutlier = bdnw.findOut_remain(new, lineInfo)
#output ./Outlier.png
#bdnw.drawWithOut(withOutlier)

isoEdge = bdnw.isoEdgeInfo(infofile, withOutlier, resList, tarList, corrMat)
#isoList, isoMzrt: list of mz, MZRT of isotope in groups
isoNodes, pr = bdnw.nodeInfo(infofile, withOutlier, resList, tarList, conList)
#get list of mz and MZRT of each group
isoList, isoMzrt = bdnw.isoGroup(pr, infofile)


#adduct


#finAdd: [resMZRT, tarMZRT, adduct type, correlation]
finAdd, addNode, addEdge = bdnw.findAdd(infofile, addfile, dataMat, corrMat, isoList, isoMzrt, mzThresh, rtThresh, corrThreshadd)
#output ./adducts network.png
#bdnw.drawInLayer(finAdd)


#C2H4

#chdict: dictionary, key:resMZRT, value:list of targetMZRT
finCH, chEdge, chNode = bdnw.findCH(infofile, dataMat, corrMat, isoList, isoMzrt, mzThresh, mzCH, rtCH)
#output ./outfile/CH.png
#bdnw.drawCH(finCH, infofile)

allEdge = bdnw.edgeInfo_all(isoEdge, addEdge, chEdge)

#bdnw.drawAll(allEdge)

unprediction = bdnw.nodelabel(infofile, anno)


prediction, name, predict = bdnw.makePredict(allEdge, infofile, anno)
#bdnw.drawPred(allEdge, unprediction)
#bdnw.drawPred(allEdge, prediction)





