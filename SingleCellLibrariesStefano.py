#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:45:46 2019

@author: stefano.magni
"""

import numpy as np
import math
import statsmodels.stats.weightstats as wstats
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.stats import pearsonr

from genes_1 import *

mpl.style.use('classic')

font05 = {'family': 'serif',
'color':  'darkred',
'weight': 'normal',
'size': 10,
}
            
font01 = {'family': 'serif',
'color':  'orange',
'weight': 'normal',
'size': 10,
}
            
font001 = {'family': 'serif',
'color':  'green',
'weight': 'normal',
'size': 10,
}

font0001 = {'family': 'serif',
'color':  'darkgreen',
'weight': 'normal',
'size': 10,
}

font2 = {'family': 'serif',
'color':  'blue',
'weight': 'normal',
'size': 15,
}

font3 = {'family': 'serif',
'color':  'black',
'weight': 'normal',
'size': 20,
}

font3s = {'family': 'serif',
'color':  'black',
'weight': 'normal',
'size': 10,
}

font3ss = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 6,
    }

def ComputeFoldChange(GenesExpressionA,GenesExpressionB):
            
            MeanA = sum(GenesExpressionA)/len(GenesExpressionA)
            MeanB = sum(GenesExpressionB)/len(GenesExpressionB)
            
            StdDevA = np.std(GenesExpressionA)
            StdDevB = np.std(GenesExpressionB)
            
            StandardErrorOfTheMeanA = StdDevA / math.sqrt(len(GenesExpressionA))
            StandardErrorOfTheMeanB = StdDevB / math.sqrt(len(GenesExpressionB))
            
            RatioMutOnCtrl = MeanA/MeanB
            MyBase = 2
            if RatioMutOnCtrl != 0:
                FoldChangeCurrentGene = math.log(RatioMutOnCtrl, MyBase) 
            else:
                FoldChangeCurrentGene = 0
            
            Temp1 = StandardErrorOfTheMeanA ** 2 / (MeanA ** 2 * math.log(2) ** 2)
            Temp2 = StandardErrorOfTheMeanB ** 2 / (MeanB ** 2 * math.log(2) ** 2)
            
            STDDEVofFoldChangeCurrentGene = math.sqrt(Temp1 + Temp2)
            
            return FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene
        
        
def ApplyZtestOnMeansOfDistributions(DataPopulation1, DataPopulation2):
    ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
    # of 2 populations are statistically significantly different ######
    DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
    DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
    MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
    TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

    return pvalue   


def autolabel_forOneComparison(ax,rects, values, ListOfFoldChangesThisDay, GeneStd, LabelUpBy, LabelDownBy):
    # attach pvalue text labels, for one comparison (one day)
    i = 0
    for rect in rects:
        height = ListOfFoldChangesThisDay
        Displacement = 0
        if ListOfFoldChangesThisDay > 0:
            Displacement = GeneStd + LabelUpBy
        elif ListOfFoldChangesThisDay < 0:
            Displacement = - GeneStd - LabelDownBy
        ax.text(rect.get_x()+rect.get_width()/2., height + Displacement, values,
                ha='center', va='bottom')
        i=i+1
       

def autolabel_forMultipleComparisons(ax, rects, values, GeneFoldChanges, GeneStd, LabelUpBy, LabelDownBy):
    # attach pvalue text labels, for multiple comparisons (days)
    i = 0
    for rect in rects:
        if GeneFoldChanges[i] > 0:
            height = GeneFoldChanges[i]
            Displacement = GeneStd[i] + LabelUpBy
        elif GeneFoldChanges[i] <= 0:
            height = GeneFoldChanges[i]
            Displacement = - GeneStd[i] - LabelDownBy
        ax.text(rect.get_x()+rect.get_width()/2., height + Displacement, values[i], #'%d'%int(height),
                ha='center', va='bottom')
        i=i+1         


def PlotFoldChangesManyGenesOneComparison(CurrentNamesList, 
                                              ListOfFoldChanges, 
                                              ListOfSTDEVforFoldChanges, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              CurrentDaysList="None",
                                              MyTitle = "", 
                                              BarsColors = "Rainbow",
                                              NamesOnEachRectangle = 'No'):
    
    N = Ncols # Number of days
    ind = np.arange(N)  # the x locations for the groups
    width = 1/(len(CurrentNamesList)+2)# 0.15        # the width of the bars
    
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)

    ColorsList = []
    if BarsColors == "Rainbow":
        for i in range(10):
            ColorsList.append('r') 
            ColorsList.append('y') 
            ColorsList.append('b') 
            ColorsList.append('c') 
            ColorsList.append('g') 
            ColorsList.append('m')
    elif BarsColors == "Bicolor":
        if  Ncols == 1:
            for FoldChange in ListOfFoldChanges:
                # print(ListOfFoldChanges)
                # print(FoldChange)
                if FoldChange >= 0:
                    ColorsList.append('Red') 
                elif FoldChange < 0:
                    ColorsList.append('Gray')                
                else:
                    ColorsList.append('Black')
                    print("PROBLEM HERE!!!")
        elif  Ncols == 2:
            ColorsMatrix = []
            for List in ListOfFoldChanges:
                ColorsList = []
                for FoldChange in List:
                    if FoldChange >= 0:
                        ColorsList.append('Red') 
                    elif FoldChange < 0:
                        ColorsList.append('Gray')                
                    else:
                        ColorsList.append('Black')
                        print("PROBLEM HERE!!!")
                ColorsMatrix.append(ColorsList)
            print(ColorsMatrix)
            print(ColorsList)
        
    ListOfRects = []
    i = 0
    if Ncols == 1:
        for IndexGenes in range(len(CurrentNamesList)):
            # print(CurrentNamesList)
            # print(IndexGenes)
            Gene0FoldChanges = ListOfFoldChanges[IndexGenes]
            Gene0Std = ListOfSTDEVforFoldChanges[IndexGenes]
            # print(len(ColorsList))
            # print(ColorsList)
            ListOfRects.append(ax.bar(ind+i*width, Gene0FoldChanges, width, color=ColorsList[i], yerr=Gene0Std))  # color=ColorsList[i],
            i = i+1
            
        ListOfSteps = []
        for step in range(len(CurrentNamesList)):
            ListOfSteps.append(ind+(step+0.5)*width)
        ax.set_xticks(ListOfSteps)
        ax.set_xticklabels(CurrentNamesList, rotation=90)
        
        ListOfAsterisks = []
        for current_pvalue in ListOfPvalues:
            if current_pvalue < 0.0001 / NBonferroni:
                Asterisks = '****'
            elif current_pvalue < 0.001 / NBonferroni:
                Asterisks = '***'
            elif current_pvalue < 0.01 / NBonferroni:
                Asterisks = '**'
            elif current_pvalue < 0.05 / NBonferroni:
                Asterisks = '*'
            else:
                Asterisks = ' '
            ListOfAsterisks.append(Asterisks)
                
        i = 0
        for i in range(len(CurrentNamesList)):
            print(i)
            autolabel_forOneComparison(ax,ListOfRects[i],ListOfAsterisks[i], ListOfFoldChanges[i], ListOfSTDEVforFoldChanges[i], LabelUpBy[i], LabelDownBy[i])
            i= i + 1
    
    elif Ncols == 2:
        for IndexGenes in range(len(CurrentNamesList)):
            Gene0FoldChanges = (ListOfFoldChanges[0][IndexGenes], ListOfFoldChanges[1][IndexGenes])
            Gene0Std = (ListOfSTDEVforFoldChanges[0][IndexGenes], ListOfSTDEVforFoldChanges[1][IndexGenes])
            if BarsColors == "Rainbow":
                MyColors = ColorsList[i]
            elif BarsColors == "Bicolor":
                MyColors = [ColorsMatrix[0][i],ColorsMatrix[1][i]]
            ListOfRects.append(ax.bar(ind+i*width, Gene0FoldChanges, width, color=MyColors, yerr=Gene0Std))
            i = i+1

        if NamesOnEachRectangle == 'No':
            ListOfSteps = []
            i = 0
            for step in range(N):
                ListOfSteps.append(ind[i]+len(CurrentNamesList)/2*width)
                i = i+1
            ax.set_xticks(ListOfSteps)
            # ax.set_xticklabels(CurrentDaysList, rotation=90)
            ax.set_xticklabels(CurrentDaysList, rotation=0, fontsize=20)
            
        elif NamesOnEachRectangle == 'Yes':
            ListOfSteps = []
            for step in range(len(CurrentNamesList)):
                ListOfSteps.append(ind+(step+0.5)*width)
            print(ListOfSteps)
            ax.set_xticks(np.array(ListOfSteps).flatten())
            
            MatrixOfNames = []
            for Day in CurrentDaysList:
                CurrentNamesList_MOD = []
                for Name in CurrentNamesList:
                    Name_MOD = Name + ' ' + Day
                    CurrentNamesList_MOD.append(Name_MOD)
                MatrixOfNames.append(CurrentNamesList_MOD)
            ax.set_xticklabels(np.array(MatrixOfNames).flatten(), rotation=45)
        
        MatrixOfAsterisks = []
        
        for i in range(len(CurrentNamesList)):
            ListOfAsterisks = []
            for j in range(Ncols):
                current_pvalue = ListOfPvalues[j][i]
                if current_pvalue < 0.0001 / NBonferroni:
                    Asterisks = '****'
                elif current_pvalue < 0.001 / NBonferroni:
                    Asterisks = '***'
                elif current_pvalue < 0.01 / NBonferroni:
                    Asterisks = '**'
                elif current_pvalue < 0.05 / NBonferroni:
                    Asterisks = '*'
                else:
                    Asterisks = ' '
                ListOfAsterisks.append(Asterisks)
            MatrixOfAsterisks.append(ListOfAsterisks)
        
        for i in range(len(CurrentNamesList)):
            MatrixOfFoldChanges = ListOfFoldChanges
            MatrixOfSEMforFoldChanges = ListOfSTDEVforFoldChanges
            CurrentGeneFoldChange = (MatrixOfFoldChanges[0][i], MatrixOfFoldChanges[1][i])
            CurrentGeneStd = (MatrixOfSEMforFoldChanges[0][i], MatrixOfSEMforFoldChanges[1][i])
            autolabel_forMultipleComparisons(ax, ListOfRects[i],
                                             (MatrixOfAsterisks[i][0],MatrixOfAsterisks[i][1]), 
                                             CurrentGeneFoldChange, CurrentGeneStd, LabelUpBy[i], LabelDownBy[i])
        # ax.legend(CurrentNamesList, loc='best')
                
    ax.set_ylabel(My_ylabel)
    ax.set_title(MyTitle)

    ax.set_axis_bgcolor('white')
    ax.axis('on')
    ax.axhline(linewidth=1, color='k')

    plt.show()

   
def ComputeCumulativeGeneExpressionFromListAndCompareObject(C, C_norm, N_cells, CurrentGenesList):
    
    CumulativeGeneExpressionA = np.zeros(N_cells)
    CumulativeGeneExpressionB = np.zeros(N_cells)
    
    for i in range(len(CurrentGenesList)): 
        GeneName = CurrentGenesList[i]
        if GeneName in C.c_names:
            CumulativeGeneExpressionA = CumulativeGeneExpressionA + C_norm[C.c_names.index(GeneName),0:N_cells]
            CumulativeGeneExpressionB = CumulativeGeneExpressionB + C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
        elif GeneName not in C.c_names:
            print('Gene ' + GeneName + ' was not found in compare object C, when computing Cumulative Gene Expression (on list with ' + str(len(CurrentGenesList)) + ' genes)')
            
    return CumulativeGeneExpressionA, CumulativeGeneExpressionB
        

def PlotHistogramsCumulatives(Nlines, 
                              Ncols, 
                              ListOfCumulatives35, 
                              ListOfCumulatives70, 
                              ListOfListsNames, 
                              ListOfPvalues, 
                              NBonferroni, 
                              DistribLabels,
                              ylabel,
                              ylimVector=[], 
                              FontListLabel=font2,
                              Colors=['Red','Blue'],
                              BinWidth = 0.05):

    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols, facecolor="white")
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide ticks and labels on X
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True) # ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_visible(True)
        
    if Ncols == 1:
        for i in range(Nlines):            
            MaxX = max(max(ListOfCumulatives35[i]),max(ListOfCumulatives70[i]))
            CumulativeGeneExpression35_NORM = ListOfCumulatives35[i] / MaxX
            CumulativeGeneExpression70_NORM = ListOfCumulatives70[i] / MaxX
        
            binwidth = BinWidth # 0.05
            axes[i].hist(CumulativeGeneExpression35_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Colors[0], normed=0, label=DistribLabels[0], alpha=0.5)        
            axes[i].hist(CumulativeGeneExpression70_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Colors[1], normed=0, label=DistribLabels[1], alpha=0.5)
                
            axes[i].set_ylabel("N Cells")
            axes[Nlines-1].set_xlabel(ylabel)
            # axes[i].set_ylim(0,220)
            # axes[2].set_ylim(0,450)
            if ylimVector: # this returns True if the vecor is not empty
                axes[i].set_ylim(0,ylimVector[i])
            axes[i].text(1.2, 0, ListOfListsNames[i], fontdict=FontListLabel) 
    
            if i == Nlines-1:
                plt.legend()
                                
            BONF = 'Cumulative Gene Expressions,\n Bonferroni correction applied for ' + str(round(NBonferroni)) + ' z-tests repetitions'
            plt.suptitle(BONF)
            
            if ListOfPvalues[i] < 0.0001 / NBonferroni:
                PValueLabel = 'p-val < 0.0001'
                font = font0001
            elif ListOfPvalues[i] < 0.001 / NBonferroni:
                PValueLabel = 'p-val < 0.001'
                font = font001
            elif ListOfPvalues[i] < 0.01 / NBonferroni:
                PValueLabel = 'p-val < 0.01'
                font = font01
            elif ListOfPvalues[i] < 0.05 / NBonferroni:
                PValueLabel = 'p-val < 0.05'
                font = font05
            else:
                PValueLabel = ' '
                font = font05
                
            axes[i].text(0.65, 15, PValueLabel, fontdict=font)
            
    elif Ncols > 1:
        for i in range(Nlines):   
            for j in range(Ncols):
                MaxX = max(max(ListOfCumulatives35[j][i]),max(ListOfCumulatives70[j][i]))
                CumulativeGeneExpression35_NORM = ListOfCumulatives35[j][i] / MaxX
                CumulativeGeneExpression70_NORM = ListOfCumulatives70[j][i] / MaxX
            
                binwidth = BinWidth # 0.05
                axes[i,j].hist(CumulativeGeneExpression35_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Colors[0], normed=0, label=DistribLabels[0], alpha=0.5)        
                axes[i,j].hist(CumulativeGeneExpression70_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color=Colors[1], normed=0, label=DistribLabels[1], alpha=0.5)
                    
                axes[i,j].set_ylabel("N Cells")
                axes[Nlines-1,j].set_xlabel(ylabel)
                # axes[i].set_ylim(0,220)
                # axes[2].set_ylim(0,450)
                if ylimVector: # this returns True if the vecor is not empty
                    axes[i,j].set_ylim(0,ylimVector[i])
                axes[i,j].text(1.2, 0, ListOfListsNames[i], fontdict=FontListLabel) 
        
                if i == Nlines-1:
                    plt.legend()
                                    
                BONF = str(ylabel) + ',\n Bonferroni correction applied for ' + str(round(NBonferroni)) + ' z-tests repetitions'
                plt.suptitle(BONF)
                
                if ListOfPvalues[j][i] < 0.0001 / NBonferroni:
                    PValueLabel = 'p-val < 0.0001'
                    font = font0001
                elif ListOfPvalues[j][i] < 0.001 / NBonferroni:
                    PValueLabel = 'p-val < 0.001'
                    font = font001
                elif ListOfPvalues[j][i] < 0.01 / NBonferroni:
                    PValueLabel = 'p-val < 0.01'
                    font = font01
                elif ListOfPvalues[j][i] < 0.05 / NBonferroni:
                    PValueLabel = 'p-val < 0.05'
                    font = font05
                else:
                    PValueLabel = ' '
                    font = font05
                    
                axes[i,j].text(0.65, 15, PValueLabel, fontdict=font)

    plt.show()
    
    
def ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                        C_norm,
                                                        List_Gene_Names_for_Correlation,
                                                        N_cells,
                                                        ThrasoldForPval,
                                                        RemoveNotSignificantCorr,
                                                        WhiteDiagonal,
                                                        SampleIndex):

    CurrentList = []
    N = len(List_Gene_Names_for_Correlation)
    pvalPearsonCorr_Matrix = 100 * np.ones([N, N])
    MyCorrelationMatrix_Full = 100 * np.ones([N, N])
    MyCorrelationMatrix_Masked = 100 * np.ones([N, N])
    
    i = 0
    for GeneName1 in List_Gene_Names_for_Correlation:
        j = 0
        CurrentList.append(GeneName1)
        for GeneName2 in List_Gene_Names_for_Correlation:
            
            Gene_1_Expression = C_norm[C.c_names.index(GeneName1), 0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            Gene_2_Expression = C_norm[C.c_names.index(GeneName2), 0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            
            PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
            pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
            
#            if RemoveNotSignificantCorr == False:
#                MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
#                              
#                if WhiteDiagonal == True:
#                    MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
#                    if i == j:
#                        MyCorrelationMatrix_Masked[i,j] = 0
                        
#            elif RemoveNotSignificantCorr == True:
            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff

            if pvalPearsonCorr <= ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
            elif pvalPearsonCorr > ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = 0
                          
            if WhiteDiagonal == True:
#                MyCorrelationMatrix_Masked[i,j] = MyCorrelationMatrix_Full[i,j]
                if i == j:
                    MyCorrelationMatrix_Masked[i,j] = 0 
                    MyCorrelationMatrix_Full[i,j] = 0 
            j = j + 1
        i = i + 1
    
    return MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked


def ComputeGeneGeneCorrelationOnSelectedCellsOfCompareObject(C,
                                                             C_norm,
                                                             List_Gene_Names_for_Correlation,
                                                             N_cells,
                                                             ThrasoldForPval,
                                                             RemoveNotSignificantCorr,
                                                             WhiteDiagonal,
                                                             IndexesOfSelectedCells):

    CurrentList = []
    N = len(List_Gene_Names_for_Correlation)
    pvalPearsonCorr_Matrix = 100 * np.ones([N, N])
    MyCorrelationMatrix_Full = 100 * np.ones([N, N])
    MyCorrelationMatrix_Masked = 100 * np.ones([N, N])
    
    i = 0
    for GeneName1 in List_Gene_Names_for_Correlation:
        j = 0
        CurrentList.append(GeneName1)
        for GeneName2 in List_Gene_Names_for_Correlation:
            
            Gene_1_Expression = C_norm[C.c_names.index(GeneName1),IndexesOfSelectedCells] # 0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            Gene_2_Expression = C_norm[C.c_names.index(GeneName2),IndexesOfSelectedCells] # 0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            
            PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
            pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
            
#            if RemoveNotSignificantCorr == False:
#                MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
#                              
#                if WhiteDiagonal == True:
#                    MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
#                    if i == j:
#                        MyCorrelationMatrix_Masked[i,j] = 0
                        
#            elif RemoveNotSignificantCorr == True:
            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff

            if pvalPearsonCorr <= ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
            elif pvalPearsonCorr > ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = 0
                          
            if WhiteDiagonal == True:
#                MyCorrelationMatrix_Masked[i,j] = MyCorrelationMatrix_Full[i,j]
                if i == j:
                    MyCorrelationMatrix_Masked[i,j] = 0 
                    MyCorrelationMatrix_Full[i,j] = 0 
            j = j + 1
        i = i + 1
    
    return MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked


def ComputeGeneGeneCorrelation_Between2DifferentGenesLists_On2SamplesOfCompareObject(C,
                                                                                     C_norm,
                                                                                     List_1_Gene_Names_for_Correlation,
                                                                                     List_2_Gene_Names_for_Correlation,
                                                                                     N_cells,
                                                                                     ThrasoldForPval,
                                                                                     RemoveNotSignificantCorr,
                                                                                     WhiteDiagonal,
                                                                                     SampleIndex):

    CurrentList = []
    N_1 = len(List_1_Gene_Names_for_Correlation)
    N_2 = len(List_2_Gene_Names_for_Correlation)
    pvalPearsonCorr_Matrix = 100 * np.ones([N_1, N_2])
    MyCorrelationMatrix_Full = 100 * np.ones([N_1, N_2])
    MyCorrelationMatrix_Masked = 100 * np.ones([N_1, N_2])
    
    i = 0
    for GeneName1 in List_1_Gene_Names_for_Correlation:
        j = 0
        CurrentList.append(GeneName1)
        for GeneName2 in List_2_Gene_Names_for_Correlation:
            
            Gene_1_Expression = C_norm[C.c_names.index(GeneName1),0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            Gene_2_Expression = C_norm[C.c_names.index(GeneName2),0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            
            PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
            pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff

            if pvalPearsonCorr <= ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
            elif pvalPearsonCorr > ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = 0
                          
            if WhiteDiagonal == True:
                if i == j:
                    MyCorrelationMatrix_Masked[i,j] = 0 
                    MyCorrelationMatrix_Full[i,j] = 0 
            j = j + 1
        i = i + 1
    
    return MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked


def ClusterDataWithGaussianMixtureModel(samples, Ncomponents, MyTitle, SaveAs='DoNotSave'):
    #fit the gaussian mixture model to the data
    gmix = mixture.GaussianMixture(n_components=Ncomponents, covariance_type='full')
    gmix.fit(samples)
    #print(gmix.means_)
    #print(' ')
    #print(gmix.covariances_)
    
    # plot the data coloured differently for each cluster
    colors = ['r' if i==0 else ('g' if i==1 else ('b' if i==2 else 
                                                  ('y' if i==3 else 
                                                   ('c' if i==4 else 
                                                    ('m' if i==5 else 
                                                     ('k' if i==6 else
                                                      ('tab:orange' if i==7 else
                                                       ('tab:purple' if i==8 else
                                                        ('tab:brown' if i==9 else
                                                         ('tab:pink' if i==10 else
                                                          ('tab:gray'))))))))))) 
                                                        for i in gmix.predict(samples)]
    figGMM = plt.figure()
    plt.title(MyTitle)
    ax = plt.gca()
    ax.scatter(samples[:,0], samples[:,1], c=colors, alpha=0.8)
    
    # Make ellipses around each claster
    ListOfAngles = []
    ListOfVs = []
    for Covar in gmix.covariances_:
        v, w = linalg.eigh(Covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ListOfAngles.append(angle)
        ListOfVs.append(v)
    ColorsList = ['r','g','b','y','c','m','k','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray']
    i=0
    for Center in gmix.means_:
        v=ListOfVs[i]
        ellipse = Ellipse(xy=Center, width=v[0], height=v[1], angle=180+ListOfAngles[i], edgecolor='black')
        ellipse.set_alpha(0.3)
        ellipse.set_facecolor(ColorsList[i])
        ax.add_artist(ellipse)
        i=i+1
        
    plt.show()
    if SaveAs != 'DoNotSave':
        plt.savefig(SaveAs)
    return gmix, figGMM


def Intersect4listsOfGenePrintNumbers(L_Aall, L_Ball, L_Call, L_Dall):
    
    L_all = []
    for GeneName in L_Aall:
        L_all.append(GeneName)
    for GeneName in L_Ball:
        L_all.append(GeneName)
    for GeneName in L_Call:
        L_all.append(GeneName)
    for GeneName in L_Dall:
        L_all.append(GeneName)

    L_all = list(set(L_all))

    L_A = []
    L_B = []
    L_C = []
    L_D = []
    
    L_AB = []
    L_AC = [] 
    L_AD = []
    L_BC = []
    L_BD = []
    L_CD = []

    L_ABC = []
    L_ABD = [] 
    L_ACD = []
    L_BCD = []
            
    L_ABCD = []

    for GeneName in L_all:
        if (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_ABCD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_BCD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Dall):
            L_ABD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball) and (GeneName in L_Call):
            L_ABC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Call) and (GeneName in L_Dall):
            L_ACD.append(GeneName)
        elif (GeneName in L_Call) and (GeneName in L_Dall):
            L_CD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Dall):
            L_BD.append(GeneName)
        elif (GeneName in L_Ball) and (GeneName in L_Call):
            L_BC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Dall):
            L_AD.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Call):
            L_AC.append(GeneName)
        elif (GeneName in L_Aall) and (GeneName in L_Ball):
            L_AB.append(GeneName)
        elif (GeneName in L_Aall):
            L_A.append(GeneName)
        elif (GeneName in L_Ball):
            L_B.append(GeneName)
        elif (GeneName in L_Call):
            L_C.append(GeneName)
        elif (GeneName in L_Dall):
            L_D.append(GeneName)
            
    print(" ")
    print(len(L_Aall))
    print(len(L_Ball))
    print(len(L_Call))
    print(len(L_Dall))
    print(" ")                
    print(len(L_all))
    print(" ")
    print(len(L_Aall))
    print(len(L_Ball))
    print(len(L_Call))
    print(len(L_Dall))
    print(" ")
    print(len(L_A))
    print(len(L_B))
    print(len(L_C))
    print(len(L_D))
    print(" ")
    print(len(L_AB))
    print(len(L_AC)) 
    print(len(L_AD))
    print(len(L_BC))
    print(len(L_BD))
    print(len(L_CD))
    print(" ")
    print(len(L_ABC))
    print(len(L_ABD)) 
    print(len(L_BCD))
    print(len(L_ACD))
    print(" ")
    print(len(L_ABCD))
    print(" ")
    print(L_AB)
    
    
def PCA_SelectingNumberOfPCsForGivenPercExplVarianceThrashold(M,ThrasholdForVarianceExplained, MaxComponents):
    print("Here the summary of performing PCA, underlying the fractions of explained variance for each component, and the total explained variance with the maximum number of components allowed by user")
    print(np.shape(M)[0])
    y_pca_temp, y_tsne_temp, PercVarExplainPCs = get_DR(M,g_ind=-1, d=MaxComponents, off=[])
    print(PercVarExplainPCs)
    print(" ")
    print(sum(PercVarExplainPCs))
    print(" ")
    TotalExplVar = 0
    IndexPCs = 0
    for PercExplVar1comp in PercVarExplainPCs:
        TotalExplVar = TotalExplVar + PercExplVar1comp
        print(TotalExplVar)
        if TotalExplVar >= ThrasholdForVarianceExplained:
            break
        else:
            IndexPCs = IndexPCs + 1
    if IndexPCs == MaxComponents:
        print("Reached the maximum components number without reaching the desired percentage of explained variance")
    return y_pca_temp, IndexPCs
