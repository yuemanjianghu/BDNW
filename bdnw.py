# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 12:02:25 2021

@author: 白藏
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import xlwt
import operator
from itertools import cycle
import math
import pandas as pd
import seaborn as sns
from scipy import stats
import igraph
from igraph import Graph
from igraph import *
import statsmodels.api as sm
import matplotlib.ticker as ticker


import sys
import json


#find all isotope pairs
#output:[resMZRT, tarMZRT, correlation]
def findIso(dataMat, corrMat, infofile, rtThresh, corrThresh, mzThresh, ratrt, ratmz):
    file = pd.read_excel(infofile)
    info = file[['MZRT_str', 'mzmed', 'rtmed', 'fimed']].values.tolist()
    #info = load.loadallinfo(infofile).values.tolist()
    #info = file[['MZRT_str', 'mzmed', 'rtmed', 'fimed']]
    info.sort(key = lambda i:i[1])
    temp = []
    
    for i in range(len(info) - 1):
        posList = []
        indexList = []
        #set threshold
        rtmin = float(float(info[i][2]) - ratrt * rtThresh)
        rtmax = float(float(info[i][2]) + (1 - ratrt) * rtThresh)
        mzmin = float(float(info[i][1]) + 1.003055 - ratmz * mzThresh)
        mzmax = float(float(info[i][1]) + 1.003055 + (1 - ratmz) * mzThresh)
        
        for j in range(i + 1, len(info)): 
            #mz threshold
            if info[j][1] >= mzmin and info[j][1] <= mzmax:
                posList.append(info[j])
                indexList.append(j)
            
            #rt threshold
            if posList !=[]:
                delList = []
                for k in range(len(posList)):
                    if posList[k][2] < rtmin or posList[k][2] > rtmax:
                        delList.append(k)
                posList = [posList[x] for x in range(len(posList)) if x not in delList]
                indexList = [indexList[x] for x in range(len(indexList)) if x not in delList]
            
            #corr threshold
            if posList != []:
                delList = []
                for m in range(len(posList)):
                    if corrMat[i][indexList[m]] <= corrThresh:
                        delList.append(m)
                posList = [posList[x] for x in range(len(posList)) if x not in delList]
                indexList = [indexList[x] for x in range(len(indexList)) if x not in delList]
        
        #output        
        for n in range(len(posList)):
            temp.append([info[i][0], posList[n][0], corrMat[i][indexList[n]]])
    finIso = pd.DataFrame(temp, columns=['resource', 'target', 'correlation'])

        # -*- coding: UTF-8 -*-
    return finIso

def drawhist(finIso, infofile, corrMat):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    mzDelta = []
    rtDelta = []
    corr = []
    
    for pair in finIso.iterrows():
        resInd = mzrt.index(pair[1]['resource'])
        tarInd = mzrt.index(pair[1]['target'])
        
        mzDelta.append(mz[tarInd] - mz[resInd])
        rtDelta.append(rt[tarInd] - rt[resInd])
        corr.append(corrMat[resInd][tarInd])
    mzplt = plt.subplot(1,3,1)
    rtplt = plt.subplot(1,3,2)
    coplt = plt.subplot(1,3,3)
    
    plt.sca(mzplt)
    mzbin = 201
    mzhist = plt.hist(mzDelta, bins = mzbin, facecolor = 'red', edgecolor='black')
    #tick_spacing = 0.0005
    #plt.tick_params(labelsize=12)
    #mzplt.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.title('mzdifferences', fontsize = 15)
    plt.xlabel('mzdifference', fontsize = 12)
    plt.ylabel('count', fontsize = 12)
    '''
    for i in range(mzbin):
        if int(mzhist[0][i]) != 0:
            plt.text(mzhist[1][i], mzhist[0][i], str(int(mzhist[0][i])))
    '''
    plt.sca(rtplt)
    rtbin = 51
    rthist = plt.hist(rtDelta, bins = rtbin, facecolor = 'blue', edgecolor='black')
    #tick_spacing = 0.0005
    #plt.tick_params(labelsize=12)
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.title('rtdifferences', fontsize = 15)
    plt.xlabel('rtdifference', fontsize = 12)
    plt.ylabel('count', fontsize = 12)
    
    plt.sca(coplt)
    cobin = 20
    cohist = plt.hist(corr, bins = cobin, facecolor = 'green', edgecolor='black')
    #tick_spacing = 0.0005
    #plt.tick_params(labelsize=12)
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.title('correlation', fontsize = 15)
    plt.xlabel('correlation', fontsize = 12)
    plt.ylabel('count', fontsize = 12)
    
    plt.savefig('./outfile/hist.png')
    plt.show() 



def ig(finIso):
    g = Graph.TupleList(finIso.itertuples(index=False), directed=False, edge_attrs="correlation")
    layout = g.layout('random')
    visual_style = {"vertex_label": g.vs["name"]}
    igraph.plot(g, './outfile/n1.png', vertex_size = 3)

#get connected compounents list
#output:(point1, point2,...)
def getConList(finIso):
    G = nx.from_pandas_edgelist(finIso, 'resource', 'target')
    G.remove_nodes_from(list(nx.isolates(G)))
    conList = list(nx.connected_components(G))
    return conList

def getConDict(conList):
    conDict = {}
    for i in range(len(conList)):
        for j in conList[i]:
            conDict[j] = i
    return conDict

def drawIsoNet(finIso):
    G = nx.from_pandas_edgelist(finIso, 'resource', 'target')
    G.remove_nodes_from(list(nx.isolates(G)))
    nx.draw(G, node_size = 5, width = 1, edge_color='r', with_labels = False)
    plt.savefig('./outfile/network without label.png')
    plt.show()

#output edge-info file
#['res', 'tar', 'component', 'mzDelta', 'rtDelta', 'correlation']
def edInfo(infofile, finIso, conList):
    conDict = getConDict(conList)
    #conDict: dictionary, key:pointMZRT, value: component number
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    edgeInfo = []
    #finIso: [resMZRT, tarMZRT, correlation]
    for pair in finIso.iterrows():
        comp = conDict[pair[1]['resource']]
        mzDelta = abs(mz[mzrt.index(pair[1]['resource'])] - mz[mzrt.index(pair[1]['target'])])
        rtDelta = abs(rt[mzrt.index(pair[1]['resource'])] - rt[mzrt.index(pair[1]['target'])])
        corr = pair[1]['correlation']
        edgeInfo.append([pair[1]['resource'], pair[1]['target'], comp, mzDelta, rtDelta, corr])

    file = xlwt.Workbook()
    sheet1 = file.add_sheet(u'sheet1', cell_overwrite_ok=True)
    m = 0
    for edge in edgeInfo:
        for n in range(len(edge)):
            sheet1.write(m, n, str(edge[n]))
        m += 1
    file.save('./outfile/edgeInfo.xls')
    edgeInfo = pd.DataFrame(edgeInfo, columns = ['resource', 'target', 'component', 'mzDelta', 'rtDelta', 'correlation'])
    return edgeInfo

def infoOfScat(conList, infofile):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    info = file[['MZRT_str', 'mzmed', 'rtmed', 'fimed']].values.tolist()
    #['MZRT_str', 'mzmed', 'rtmed', 'fimed']
    scat = []
    resList = []
    tarList = []
    forIso = {}
    for comp in conList:
        forSort = []
        for samp in comp:
            index = mzrt.index(samp)
            sampmz = [samp, float(info[index][1])]
            forSort.append(sampmz)
        forSort.sort(key = operator.itemgetter(1))
        for i in range (len(forSort)):
            forIso[forSort[i][0]] = i
    #forIso: dictionary, key: metaMZRT, value: predicted isotope type
    for comp in conList:
        forSort = []
        for samp in comp:
            index = list(mzrt).index(samp)
            sampmz = [samp, float(info[index][1])]
            forSort.append(sampmz)
        forSort.sort(key = operator.itemgetter(1))
        compSort = []
        for name in forSort:
            compSort.append(name[0])
        for i in range (len(compSort)- 1):
            reference = list(compSort)[i]
            target = list(compSort)[i + 1]
            resList.append(reference)
            tarList.append(target)
            refIndex = list(mzrt).index(reference)
            tarIndex = list(mzrt).index(target)
            scat.append([info[refIndex][1], info[tarIndex][3] / info[refIndex][3], forIso[reference]])
    scat = pd.DataFrame(scat, columns = ['mzmed', "fimed_t / fimed_r", "type"])
    return scat, resList, tarList

def drawScatOfPair(scat):
    sns.lmplot(x='mzmed', y="fimed_t / fimed_r", data = scat, scatter_kws={'s':5}, hue='type', ci = 95, fit_reg = False)
    plt.savefig('./outfile/scatComp.png')
    plt.show()

def drawScatOfPair_line(scat):
    fig,axes=plt.subplots(1,scat['type'].max()) 
    lineInfo = []
    names = locals()
    for i in range(scat['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
    for pair in scat.iterrows():
        for j in range(scat['type'].max() + 1):
            if pair[1]['type'] == j:
                names['x' + str(j)].append(pair[1]['mzmed'])
                names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
    for k in range(scat['type'].max() + 1):
        slope, intercept, r_value, p_value, std_err = stats.linregress(names['x' + str(k)], names['y' + str(k)])
        lineInfo.append([slope, intercept])
    lineInfo = pd.DataFrame(lineInfo, columns = ['slope', 'intercept'])
    sns.lmplot(x='mzmed', y="fimed_t / fimed_r", data = scat, scatter_kws={'s':5}, hue='type', 
               ci = 95, line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
    plt.savefig('./outfile/scatComp robust.png')
    plt.show()
    return lineInfo

def drawScatOfPair_intercept0(scat):
    lineInfo = []
    names = locals()
    for i in range(scat['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
    
    for pair in scat.iterrows():
        for j in range(scat['type'].max() + 1):
            if pair[1]['type'] == j:
                names['x' + str(j)].append(pair[1]['mzmed'])
                names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
    
    for k in range(scat['type'].max() + 1):
        y = names['y' + str(k)]
        x = names['x' + str(k)]
        model = sm.RLM(y, x, M=sm.robust.norms.HuberT()).fit()
        predicts = model.predict()
        plt.scatter(x, y, s = 5, label='M+' + str(k) + ' - M+' + str(k + 1))
        plt.plot(x, predicts, label='y = ' + str(format(model.params[0], '.6f')) + 'x')
        plt.legend(fontsize=15)
        plt.xlabel('mzmed_r', fontsize=15)
        plt.ylabel('fimed_t / fimed_r', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()
        lineInfo.append([model.params[0], 0])
    lineInfo = pd.DataFrame(lineInfo, columns = ['slope', 'intercept'])
    plt.savefig('./outfile/scatComp robust.png')
    return lineInfo
    
def distance(pointX, pointY, lineA, lineB):
    lineY = lineA * pointX + lineB
    lineX = (pointY - lineB) / lineA
    xDelta = lineX - pointX
    yDelta = lineY - pointY
    pro = lineX * lineY - pointX * pointY
    dis=(math.fabs(yDelta*pointX+xDelta*pointY+pro))/(math.pow(yDelta * yDelta + xDelta * xDelta, 0.5))
    return dis

def findWrong(scat, lineInfo):
    names = globals()        
    for i in range(scat['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
        names['i' + str(i)] = []
        
    for pair in scat.iterrows():
        for j in range(scat['type'].max() + 1):
            if pair[1]['type'] == j:
                names['x' + str(j)].append(pair[1]['mzmed'])
                names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
                names['i' + str(j)].append(pair[1]['type'])
    
    for k in range(scat['type'].max() + 1):
        if len(names['x' + str(k)]) > 1:
            names['line' + str(k)] = np.polyfit( names['x' + str(k)], 
                                             names['y' + str(k)], 1)
            names['lineA' + str(k)] = names['line' + str(k)][0]
            names['lineB' + str(k)] = names['line' + str(k)][1]

    wrongList = []
    realList = []
    #realList: corrected resource isotope type of each isotope pair
    for p in range(scat.shape[0]):
        pointX = scat[p:p + 1]['mzmed']
        pointY = scat[p:p + 1]['fimed_t / fimed_r']
        pointIso = scat[p:p + 1]['type']
        disList = []
        for m in range(scat['type'].max() + 1):
            dis = distance(pointX, pointY, float(lineInfo[m:m + 1]['slope']), float(lineInfo[m:m + 1]['intercept']))
            disList.append(dis)
        realiso = disList.index(min(disList))
        realList.append(realiso)
    new = []
    for q in range(scat.shape[0]):
        new.append([float(scat[q:q + 1]['mzmed']), float(scat[q:q + 1]['fimed_t / fimed_r']), realList[q]])
    new = pd.DataFrame(new, columns = ['mzmed', "fimed_t / fimed_r", "type"])
    return new
'''
def drawNoWrong(new):
    sns.lmplot(x='mzmed', y="fimed_t / fimed_r", data = new, scatter_kws={'s':5}, hue='type', ci = 95)
    plt.savefig('./outfile/noWrong.png')
'''

def drawNoWrong(new):
    names = locals()
    for i in range(new['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
    
    for pair in new.iterrows():
        for j in range(new['type'].max() + 1):
            if pair[1]['type'] == j:
                names['x' + str(j)].append(pair[1]['mzmed'])
                names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
    
    for k in range(new['type'].max() + 1):
        y = names['y' + str(k)]
        x = names['x' + str(k)]
        model = sm.RLM(y, x, M=sm.robust.norms.HuberT()).fit()
        predicts = model.predict()
        plt.scatter(x, y, s = 5, label='M+' + str(k) + ' - M+' + str(k + 1))
        plt.plot(x, predicts, label='y = ' + str(format(model.params[0], '.6f')) + 'x')
        plt.legend()
        plt.xlabel('mzmed_r')
        plt.ylabel('fimed_t / fimed_r')
        plt.show()
    plt.savefig('./outfile/no Wrong.png')
    
def dishist(new, lineInfo):
    disList = []
    for point in new.iterrows():
        if point[1]['type'] == 0:
            pointX = point[1]['mzmed']
            pointY = point[1]['fimed_t / fimed_r']
            dis = distance(pointX, pointY, float(lineInfo.iat[0, 0]), float(lineInfo.iat[0, 1]))
            disList.append(dis)
    disbin = 50
    dishist = plt.hist(disList, bins = disbin, facecolor = 'red', edgecolor='black')
    plt.title('distance')
    plt.xlabel('distance')
    plt.ylabel('count')
    plt.savefig('./outfile/dishist.png')
    
def isoAll(new, resList, tarList, infofile, corrMat):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    correctIso = []
    names = locals()
    for i in range (new['type'].max() + 1):
        names['is' + str(i)] = {}
    for j in range(len(resList)):
        for k in range(new['type'].max() + 1):
            if int(new[j:j + 1]['type']) == k:
                res = mzrt.index(resList[j])
                tar = mzrt.index(tarList[j])
                names['is' + str(k)][resList[j]] = tarList[j] + "&" + str(k) + "_" + str(k + 1) + "&" + str(corrMat[res][tar])
    for m in range (new['type'].max() + 1):
        correctIso.append(names['is' + str(m)])
    
    return correctIso

def getSD(new):
    sdList = []
    names = globals()        
    for i in range(new['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
        names['i' + str(i)] = []
        names['d' + str(i)] = []
        
    for pair in new.iterrows():
        for j in range(new['type'].max() + 1):
            if pair[1]['type'] == j:
                names['x' + str(j)].append(pair[1]['mzmed'])
                names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
                names['i' + str(j)].append(pair[1]['type'])            
    
    
    
    #get standard deviation list
    for k in range(new['type'].max() + 1):
        if len(names['x' + str(k)]) > 1:
            names['line' + str(k)] = np.polyfit( names['x' + str(k)], 
                                             names['y' + str(k)], 1)
            names['lineA' + str(k)] = names['line' + str(k)][0]
            names['lineB' + str(k)] = names['line' + str(k)][1]
            
        for m in range(len(names['x' + str(k)])):
            dis = distance(names['x' + str(k)][m], names['y' + str(k)][m],
                           names['lineA' + str(k)], names['lineB' + str(k)])
            names['d' + str(k)].append(dis)
        #sdList = [sd, lineA, lineB, mean]
        sdList.append([np.std(names['d' + str(k)]), names['lineA' + str(k)],
                       names['lineB' + str(k)], np.mean(names['d' + str(k)])])
    return sdList
   
def findOut_sd(new, sdnum):
    newList = []
    accurate = []
    sdList = getSD(new)
    names = locals() 
    
    for pair in new.iterrows():
        disList = []
        pointX = pair[1]['mzmed']
        pointY = pair[1]['fimed_t / fimed_r']
        pointI = pair[1]['type']
        for j in range(len(sdList)):
            sd = sdList[j][0]
            lineA = sdList[j][1]
            lineB = sdList[j][2]
            mean = sdList[j][3]
            dis = distance(pointX, pointY, lineA, lineB)
            if abs(dis - mean) <= (sdnum * sd):
                disList.append(j)
        if disList == []:
            newList.append([pair[1]['mzmed'], pair[1]['fimed_t / fimed_r'], -1])
        elif pointI in disList:
            newList.append([pair[1]['mzmed'], pair[1]['fimed_t / fimed_r'], int(pair[1]['type'])])
        else:
            newList.append([pair[1]['mzmed'], pair[1]['fimed_t / fimed_r'], int(disList[0])])
    withOutlier = pd.DataFrame(newList, columns = ['mzmed', "fimed_t / fimed_r", "type"])
    return withOutlier

def findOut_remain(new, lineInfo):
    remain = []
    names = locals()
    for i in range(new['type'].max() + 1):
        names['remain' + str(i)] = []
    for point in new.iterrows():
        for j in range(new['type'].max() + 1):
            if point[1]['type'] == j:
                names['remain' + str(j)].append(float(abs(lineInfo[j : j + 1]['slope'] * point[1]['mzmed'] + 
                                                    lineInfo[j : j + 1]['intercept'] - point[1]['fimed_t / fimed_r'])))
                remain.append(float(abs(lineInfo[j : j + 1]['slope'] * point[1]['mzmed'] + 
                                                    lineInfo[j : j + 1]['intercept'] - point[1]['fimed_t / fimed_r'])))
    for k in range (new['type'].max() + 1):
        names['lower' + str(k)] = np.quantile(names['remain' + str(k)],0.25,interpolation='lower')#下四分位数
        names['higher' + str(k)] = np.quantile(names['remain' + str(k)],0.75,interpolation='higher')#上四分位数
        names['int_r' + str(k)] = names['higher' + str(k)] - names['lower' + str(k)]
        names['top' + str(k)] = names['higher' + str(k)] + 1.5 * names['int_r' + str(k)]
        names['bottom' + str(k)] = names['lower' + str(k)] - 1.5 * names['int_r' + str(k)]
    newList = []
    for m in range (len(remain)):
        typem = int(new[m:m + 1]['type'])
        if remain[m] >= names['bottom' + str(typem)] and remain[m] <= names['top' + str(typem)]:
            newList.append([new.iat[m, 0], new.iat[m, 1], new.iat[m, 2]])
        else:
            newList.append([new.iat[m, 0], new.iat[m, 1], -1])
    withOutlier = pd.DataFrame(newList, columns = ['mzmed', "fimed_t / fimed_r", "type"])
    return withOutlier
     
    

def drawWithOut(withOutlier):
    outlierx = []
    outliery = []
    names = locals()
    for i in range(0, withOutlier['type'].max() + 1):
        names['x' + str(i)] = []
        names['y' + str(i)] = []
    
    for pair in withOutlier.iterrows():
        if pair[1]['type'] == -1:
            outlierx.append(pair[1]['mzmed'])
            outliery.append(pair[1]['fimed_t / fimed_r'])
        else:
            for j in range(0, withOutlier['type'].max() + 1):
                if pair[1]['type'] == j:
                    names['x' + str(j)].append(pair[1]['mzmed'])
                    names['y' + str(j)].append(pair[1]['fimed_t / fimed_r'])
    plt.scatter(outlierx, outliery, s = 20, label='outlier', c = 'black')
    for k in range(0, withOutlier['type'].max() + 1):
        y = names['y' + str(k)]
        x = names['x' + str(k)]
        model = sm.RLM(y, x, M=sm.robust.norms.HuberT()).fit()
        predicts = model.predict()
        plt.scatter(x, y, s = 5, label='M+' + str(k) + ' - M+' + str(k + 1))
        plt.plot(x, predicts, label='y = ' + str(format(model.params[0], '.6f')) + 'x')
        plt.legend(fontsize=15)
        plt.xlabel('mzmed_r', fontsize=15)
        plt.ylabel('fimed_t / fimed_r', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()
    plt.savefig('./outfile/outlier.png')
    

def nodeInfo(infofile, withOutlier, resList, tarList, conList):
    conDict = getConDict(conList)
    file = pd.read_excel(infofile)
    info = file[['MZRT_str', 'mzmed', 'rtmed', 'fimed']].values.tolist()
    #info = ['MZRT_str', 'mzmed', 'rtmed', 'fimed']
    nodes = []
    #nodes = ['MZRT_str', 'mzmed', 'rtmed', 'fimed', 'component', 'predictIso']
    pr = []
    
    for i in info:
        if i[0] in resList:
            ind = int(resList.index(i[0]))
            temp = int(withOutlier[ind:ind + 1]['type'])
            if temp != -1:
                pos = 'M+' + str(temp)
                pr.append(temp)
            elif temp == -1:
                pos = 'outlier-res'
                pr.append(0)
        elif i[0] in tarList:
            ind = tarList.index(i[0])
            if temp != -1:
                pos = 'M+' + str(temp + 1)
                pr.append(temp + 1)
            elif temp == -1:
                pos = 'outlier-tar'
                pr.append(0)
        else:
            pos = 'isolated'
            pr.append(0)
            
        if i[0] in conDict.keys():
            comp = conDict[i[0]]
        else:
            comp = 'isolated'
        nodes.append([i[0], i[1], i[2], i[3], comp, pos])
    
    isoNodes = pd.DataFrame(nodes, columns = ['MZRT_str', 'mzmed', 'rtmed', 'fimed', 'component', 'predictIso'])
    return isoNodes, pr

def isoEdgeInfo(infofile, new, resList, tarList, corrMat):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    isoEdge = []
    for i in range(len(resList)):
        restype = int(new.iat[i , 2])  
        corr = corrMat[mzrt.index(resList[i])][mzrt.index(tarList[i])]
        isoEdge.append([resList[i], tarList[i], 'isotope', '+' + str(restype) + ' - +' + str(restype + 1), corr])
    isoEdge = pd.DataFrame(isoEdge, columns = ['resource', 'target', 'type', 'detail', 'correlation'])
    return isoEdge
    
      
def isoGroup(pr, infofile):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    isoList = []
    isoMzrt = []
    
    names = locals()
    for i in range (max(pr) + 1):
        names['isomz' + str(i)] = []
        names['isomzrt' + str(i)] = []
    
    for p in range(len(pr)):
        for k in range (max(pr) + 1):
            if pr[p] == k:
                names['isomz' + str(k)].append(mz[p])
                names['isomzrt' + str(k)].append(mzrt[p])
    
    for j in range (max(pr) + 1):
        isoList.append(names['isomz' + str(j)])
        isoMzrt.append(names['isomzrt' + str(j)])
    
    return isoList, isoMzrt


def findAdd(infofile, addfile, dataMat, corrMat, isoList, isoMzrt, mzThresh, rtThresh, corrThreshadd):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    addfile = pd.read_excel(addfile)
    addinfo = addfile[['ID', 'likelihoodGroup', 'chargeFactor', 'massToSubtract', 'massToAdd', 'oligomerIndex']].values.tolist()

    finAdd = []
    addNode = []
    addEdge = []
    char = []
    mas = []
    count = 0
    for ad in addinfo:
        #if ad[1] == 1:
        if count <= 8:
            char.append(ad[2])
            mas.append(ad[4] - ad[3])
        count += 1
    
    names = locals()
    for i in range(len(isoList)):
        names['isomz' + str(i)] = isoList[i]
        names['isomzrt' + str(i)] = isoMzrt[i]
        names['addList' + str(i)] = []
    
        for k in range(len(names['isomz' + str(i)]) - 1):
            for j in range (k + 1, len(names['isomz' + str(i)])):
                refInd = mzrt.index(names['isomzrt' + str(i)][k])
                tarInd = mzrt.index(names['isomzrt' + str(i)][j])
                if abs(rt[refInd] - rt[tarInd]) <= rtThresh and corrMat[refInd][tarInd] >= corrThreshadd:
                    refmass = []
                    tarmass = []
                    for p in range(len(char)):
                        refmass.append(names['isomz' + str(i)][k] * char[p] + mas[p])
                        tarmass.append(names['isomz' + str(i)][j] * char[p] + mas[p])
                        
                    same = []
                    for ref in range (len(refmass)):
                        mzmin = refmass[ref] - mzThresh
                        mzmax = refmass[ref] + mzThresh
                        for tar in range(len(tarmass)):
                            if tarmass[tar] >= mzmin and tarmass[tar] <= mzmax:
                                same.append((ref, tar))

                    if len(same) != 0:
                        for m in same:
                            refID = addinfo[m[0]][0]
                            tarID = addinfo[m[1]][0]
                            names['addList' + str(i)].append([names['isomzrt' + str(i)][k], names['isomzrt' + str(i)][j], 
                                                              str(refID) + ' - ' + str(tarID), abs(mz[refInd] - mz[tarInd]), 
                                                              abs(rt[refInd] - rt[tarInd]), corrMat[refInd][tarInd]])
                            addNode.append([names['isomzrt' + str(i)][k], mz[refInd], rt[refInd], 'adduct', str(refID)])
                            addNode.append([names['isomzrt' + str(i)][j], mz[tarInd], rt[tarInd], str(tarID)])
                            addEdge.append([names['isomzrt' + str(i)][k], names['isomzrt' + str(i)][j], 'adduct', 
                                           str(refID) + ' - ' + str(tarID), corrMat[refInd][tarInd]])
        names['addList' + str(i)] = pd.DataFrame(names['addList' + str(i)], columns = ['resource', 'target', 'adduct type', 'mzDelta', 'rtDelta', 'correlation'])
        finAdd.append(names['addList' + str(i)])
    addNode = pd.DataFrame(addNode, columns = ['MZRT_str', 'mzmed', 'rtmed', 'type', 'detail'])
    addEdge = pd.DataFrame(addEdge, columns = ['resource', 'target', 'type', 'detail', 'correlation'])
    return finAdd, addNode, addEdge
    

def drawInLayer(finAdd):
    names = locals()
    for i in range (len(finAdd)):
        names['pic' + str(i)] = finAdd[i]
        #plt.subplot(1, len(finAdd), i + 1)
        g = Graph.TupleList(names['pic' + str(i)].itertuples(index=False), directed=False, edge_attrs="correlation")
        layout = g.layout('random')
        visual_style = {"vertex_label": g.vs["name"]}
        igraph.plot(g, './outfile/add' + str(i) + '.png', vertex_size = 3)

    
    
def findCH(infofile, dataMat, corrMat, isoList, isoMzrt, mzThresh, mzCH, rtCH):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    finCH = []
    chNode = []
    chEdge = []
    
    names = locals()
    for i in range(len(isoList)):
        names['isomz' + str(i)] = isoList[i]
        names['isomzrt' + str(i)] = isoMzrt[i]
        names['chList' + str(i)] = []
        names['isomz' + str(i)].sort()
    
        for k in range(len(names['isomz' + str(i)]) - 1):
            indexList = []
            resInd = mzrt.index(names['isomzrt' + str(i)][k])
            rtmin = float(float(rt[resInd]))
            rtmax = float(float(rt[resInd]) + rtCH)
            mzmin = float(float(names['isomz' + str(i)][k]) + mzCH - mzThresh)
            mzmax = float(float(names['isomz' + str(i)][k]) + mzCH + mzThresh)
            
            for j in range (k + 1, len(names['isomz' + str(i)])):
                tarInd = mzrt.index(names['isomzrt' + str(i)][j])
                #mz threshold
                if names['isomz' + str(i)][j] >= mzmin and names['isomz' + str(i)][j] <= mzmax:
                    if rt[tarInd] >= rtmin and rt[tarInd] <= rtmax:
                        names['chList' + str(i)].append([mzrt[resInd], mzrt[tarInd], corrMat[resInd][tarInd]])
                        chEdge.append([mzrt[resInd], mzrt[tarInd], 'C2H4', 'C2H4 - C2H4', corrMat[resInd][tarInd]])
                        chNode.append([mzrt[resInd], mz[resInd], rt[resInd], 'C2H4', 'C2H4'])
                        chNode.append([mzrt[tarInd], mz[tarInd], rt[tarInd], 'C2H4', 'C2H4'])
        names['chList' + str(i)] = pd.DataFrame(names['chList' + str(i)], columns = ['resource', 'target', 'correlation'])
        finCH.append(names['chList' + str(i)])
    chEdge = pd.DataFrame(chEdge, columns = ['resource', 'target', 'type', 'detail', 'correlation'])
    chNode = pd.DataFrame(chEdge, columns = ['MZRT_str', 'memed', 'rtmed', 'type', 'detail'])
    return finCH, chEdge, chNode

def drawCH(finCH, infofile):
    file = pd.read_excel(infofile)
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    names = locals()
    '''
    for i in range (len(finCH)):
        names['pic' + str(i)] = finCH[i]
        names['x' + str(i)] = []
        names['y' + str(i)] = []
        
        if len(names['pic' + str(i)]) != 0:
            for pair in names['pic' + str(i)].iterrows():
                y = [mz[mzrt.index(pair[1]['resource'])], mz[mzrt.index(pair[1]['target'])]]
                x = [rt[mzrt.index(pair[1]['resource'])], rt[mzrt.index(pair[1]['target'])]]
                #plt.subplot(1, len(finCH), i + 1)
                plt.scatter(x, y, s = 5)
                plt.plot(x, y)
                plt.xlabel('rtmed', fontsize=15)
                plt.ylabel('mzmed', fontsize=15)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.show()
                
                plt.savefig('./outfile/ch' + str(i) + '.png')
    '''
    i = 3
    names['pic' + str(i)] = finCH[i]
    names['x' + str(i)] = []
    names['y' + str(i)] = []
        
    if len(names['pic' + str(i)]) != 0:
        for pair in names['pic' + str(i)].iterrows():
            y = [mz[mzrt.index(pair[1]['resource'])], mz[mzrt.index(pair[1]['target'])]]
            x = [rt[mzrt.index(pair[1]['resource'])], rt[mzrt.index(pair[1]['target'])]]
            #plt.subplot(1, len(finCH), i + 1)
            plt.scatter(x, y, s = 5, c = 'blue')
            plt.plot(x, y, c = 'darkseagreen')
            plt.title('M+3 layer', fontsize = 15)
            plt.xlabel('rtmed', fontsize=15)
            plt.ylabel('mzmed', fontsize=15)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.show()

            


def edgeInfo_all(isoEdge, addEdge, chEdge):
    frames = [isoEdge, addEdge, chEdge]
    allEdge = pd.concat(frames)
    return allEdge

def drawAll(allEdge):
    
    allEdge['type']=pd.Categorical(allEdge['type'])
    allEdge['type'].cat.codes
    '''
    G = nx.Graph()
    for pair in allEdge.iterrows():
        G.add_node(pair[1]['resource'], node_size = 1)
        G.add_node(pair[1]['target'], node_size = 1)
        if pair[1]['type'] == 'isotope':
            G.add_edge(pair[1]['resource'],pair[1]['target'],color='r')
        elif pair[1]['type'] == 'adduct':
            G.add_edge(pair[1]['resource'],pair[1]['target'],color='g')
        else:
            G.add_edge(pair[1]['resource'],pair[1]['target'],color='b')
    '''
    G = nx.from_pandas_edgelist(allEdge, 'resource', 'target')
    nx.draw(G, node_size = 8, width = 2, with_labels = False, 
            edge_color = allEdge['type'].cat.codes, edge_cmap=plt.cm.rainbow)
    '''
    pos = nx.random_layout(G)
    nx.draw(G, pos)
    plt.savefig('./outfile/alledge.png')
    plt.show()
    g = Graph.TupleList(allEdge.itertuples(index=False), directed=False, edge_attrs="type")
    layout = g.layout('random')
    visual_style = {"vertex_label": g.vs["name"]}
    igraph.plot(g, './outfile/111.png', vertex_size = 3)
    '''
def nodelabel(infofile, anno):
    file = pd.read_excel(infofile)
    annot = file[anno].values.tolist()
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    
    unprediction = []
    for n in range(len(mzrt)):
        unprediction.append([mzrt[n], mz[n], rt[n], annot[n]])
    unprediction = pd.DataFrame(unprediction, columns=(['MZRT_str', 'mzmed', 'rtmed', 'prediction']))
    return unprediction

def drawPred(allEdge, prediction):
    name = []
    G = nx.from_pandas_edgelist(allEdge, 'resource', 'target')
    conList = list(nx.connected_components(G))
    for i in conList:
        for j in i:
            name.append(j)
    edge_labels = {}
    ed = 0
    for edge in G.edges:
        edge_labels[edge] = allEdge.iat[ed, 2]
        ed += 1
    
    allEdge['type'] = pd.Categorical(allEdge['type'])
    #allEdge['type'].cat.codes
    color = []
    node_labels = {}
    nd = -1
    edgeiso = []
    edgeadd = []
    edgech = []
    
    for edge in allEdge.iterrows():
        if edge[1]['type'] == 'isotope':
            edgeiso.append((edge[1]['resource'], edge[1]['target']))
        elif edge[1]['type'] == 'adduct':
            edgeadd.append((edge[1]['resource'], edge[1]['target']))
        else:
            edgech.append((edge[1]['resource'], edge[1]['target']))
            
    for edg in edgeiso:
        G[edg[0]][edg[1]]['color'] = 'salmon'
    for edg in edgeadd:
        G[edg[0]][edg[1]]['color'] = 'darkorange'
    for edg in edgech:
        G[edg[0]][edg[1]]['color'] = 'limegreen'
    edgecolors = nx.get_edge_attributes(G,'color').values()
    
    for node in prediction.iterrows():
        if node[1]['MZRT_str'] in name:
            nd += 1
            if node[1]['prediction'] != 'isolated':
                color.append(node[1]['prediction'])
                node_labels[node[1]['MZRT_str']] = prediction.iat[nd, 3]
    color = pd.DataFrame(color, columns=(['prediction']))
    color['prediction'] = pd.Categorical(color['prediction'])
    color['prediction'].cat.codes
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_size = 12, width = 2, with_labels = False, node_color = color['prediction'].cat.codes, 
            edge_color = edgecolors, node_cmap = plt.cm.rainbow)
    #nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)
    nx.draw_networkx_labels(G, pos, node_labels, font_size=12)
    plt.savefig('./outfile/alledge.png')
    plt.show()
    
def makePredict(allEdge, infofile, anno):
    file = pd.read_excel(infofile)
    annot = file[anno].values.tolist()
    mzrt = file['MZRT_str'].values.tolist()
    mz = file['mzmed'].values.tolist()
    rt = file['rtmed'].values.tolist()
    name = []
    predict = []
    G = nx.from_pandas_edgelist(allEdge, 'resource', 'target')
    conList = list(nx.connected_components(G))
    nodeType = {}
    for ind in range(len(allEdge)):
        if allEdge.iat[ind, 0] not in nodeType.keys():
            nodeType[allEdge.iat[ind, 0]] = [allEdge.iat[ind, 3].split(' - ')[0]]
        else:
            nodeType[allEdge.iat[ind, 0]].append(allEdge.iat[ind, 3].split(' - ')[0])
        if allEdge.iat[ind, 1] not in nodeType.keys():
            nodeType[allEdge.iat[ind, 1]] = [allEdge.iat[ind, 3].split(' - ')[1]]
        else:
            nodeType[allEdge.iat[ind, 1]].append(allEdge.iat[ind, 3].split(' - ')[1])
    for i in range(len(conList)):
        temp = []
        for j in conList[i]:
            ind = mzrt.index(j)
            name.append(j)
            if pd.isnull(annot[ind]):
                temp.append(0)
            else:
                temp.append(annot[ind])
        for k in conList[i]:
            if any(temp) == False:
                predict.append('nan')
            elif len(set(temp)) == 2:
                for m in temp:
                    if m != 0:
                        predict.append(m)
            elif len(set(temp)) > 2:
                pre = 'possible:'
                for m in temp:
                    if m != 0 and pre.split(':')[1] != m:
                        pre += str(m) + ', '
                predict.append(pre)
                
    prediction = []
    for n in range(len(mzrt)):
        if pd.isnull(annot[n]):
            if mzrt[n] in name:
                if nodeType[mzrt[n]][0].startswith("+") and len(nodeType[mzrt[n]]) == 1:
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], predict[name.index(mzrt[n])], nodeType[mzrt[n]][0], 'nan'])
                elif nodeType[mzrt[n]][0].startswith("+") and nodeType[mzrt[n]][1] != 'C2H4':
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], predict[name.index(mzrt[n])], nodeType[mzrt[n]][0], nodeType[mzrt[n]][1]])
                elif nodeType[mzrt[n]][0] != 'C2H4':
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], predict[name.index(mzrt[n])], 'isolated', nodeType[mzrt[n]][0]])
                else:
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], predict[name.index(mzrt[n])], 'isolated', 'nan'])
            else:
                prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], annot[n], 'isolated', 'nan'])
        else:
            if mzrt[n] in name:
                if nodeType[mzrt[n]][0].startswith("+") and len(nodeType[mzrt[n]]) == 1:
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], annot[n], nodeType[mzrt[n]][0], 'nan'])
                elif nodeType[mzrt[n]][0].startswith("+") and nodeType[mzrt[n]][1] != 'C2H4':
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], annot[n], nodeType[mzrt[n]][0], nodeType[mzrt[n]][1]])
                elif nodeType[mzrt[n]][0] != 'C2H4':
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], annot[n], 'isolated', nodeType[mzrt[n]][0]])
                else:
                    prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], predict[name.index(mzrt[n])], 'isolated', 'nan'])
            else:
                prediction.append(['N' + str(n), mzrt[n], mz[n], rt[n], annot[n], annot[n], 'isolated', 'nan'])

    file = xlwt.Workbook()
    sheet1 = file.add_sheet(u'MESA', cell_overwrite_ok=True)
    nam = [['node', 'MZRT_str', 'mzmed', 'rtmed', 'annotation', 'prediction', 'isotope', 'adduct']] + prediction
    m = 0
    for pre in nam:
        for n in range(len(pre)):
            sheet1.write(m, n, str(pre[n]))
        m += 1
    file.save('./outfile/prediction_MESA.xls')
    prediction = pd.DataFrame(prediction, columns=(['node', 'MZRT_str', 'mzmed', 'rtmed', 'annotation', 'prediction', 'isotope', 'adduct']))
    return prediction, name, predict
