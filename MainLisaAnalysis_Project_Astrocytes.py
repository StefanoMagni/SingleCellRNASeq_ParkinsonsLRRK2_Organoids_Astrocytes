#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Mon Feb 18 2019
# @authors: stefano.magni
    
from scipy.stats.stats import pearsonr
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.ticker as plticker
import seaborn as sns
import copy
import math
import statsmodels.stats.weightstats as wstats
import matplotlib as mpl
import matplotlib.colors as colors
import statistics
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.patches as mpatches
    
from Script_data import *
from Script_information import *
from genes_1 import *
import Script_data as sd

from SingleCellLibrariesStefano import *

# Tags to decide what to run
# Pre Processing
GetDataFromFiles = 0                                                               # Pre-processing
ComputeKneePlots = 0                                                                # Pre-processing
PlotKneePlots = 0                                                                   # Pre-processing
CutCellNumber = 0                                                                   # Pre-processing 
Normalize = 0                                                                       # Pre-processing
ComparePreliminaries = 0                                                            # Pre-processing

ProduceListsOfDifferentiallyExpressedGenes = 0                                      # DEGs Mu VS Wt(Differentially Expressed Genes)
IdentifyGenesActuallypresentInData = 0                                              # Count DEGs in lists and isolate absent genes

DoComparisonDR = 0                                                                  # DR day by day
DoComparisonDR_ColoredByScores = 0                                                  # DR day by day

##### NEW METHOD TO SEPARATE ASTROCYTES BASED ON CLUSTERING #####
Compare_All_Wt_Mu_35_70_Preliminaries = 0
GenesToUseInALLComparison = "SelectedGenes" # "AllGenes" # Selceted genes are those relevant for cell type identification
DoComparisonDR_ALL = 0
UsePCApreviousTotSNE_ALL = 0
UsePCApreviousTotSNE_ALL_AlsoColoredByScores = 0
DetermineOptimalNumberOfClustersBySilhouetteAnalysis = 0
ClusteringOnFullDEMorPCs = "PCs" # "AllDEM" = run clustering on the DEM after cuts but before DR "PCs" = do it after DR by PCA
DoKmeansClusteringOnHDdata_PlotItAgainstPCApreviousTotSNE_ALL = 0 
ClustersIdentificationByVisualizationOfMarkersGeneExpression = 0
UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_NFIA = 0
UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_TGFBI = 0
UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_NR2F1 = 0

##### OLD METHOD TO SEPARATE ASTROCYTES BASED ON CUMULATIVE SCORES #####
Plot3d_Astro_Neurons_Stemness_SCORES = 0                                            # 3D SCORES
Separate_Cells_By_Cell_Type_From_3D_Plot3d = 0                                      # USE MEAN 3D SCORES TO SEPARATE CELLS
Threshold_Celltype = "Median" # Specify to use "Mean" or "Median" or "Max" as Threshold for determining cell type. Mutually exclusives.

SeparateAstrocytesUsingCumulativesOrClustering = "Clustering" # "Cumulatives"

HistoAndFoldchangesGeneListsGeneralCellTypesProcesses_ASTRO_ONLY = 0                # LISTS OF GENES ASTRO ONLY
HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES_ASTRO_ONLY=0  # CRC GENES ONLY ASTRO

HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES = 0           # CRC GENES ALL CELLS
HistoAndFoldChangesIndividualGenes_AstroList = 0
HistoAndFoldChangesIndividualGenes_AstroList_ASTRO_ONLY = 0
HistoAndFoldChangesGeneListsGeneralCellTypesProcesses = 0                           # LISTS OF GENES ALL CELLS
HistoAndFoldChangesGeneListsGeneralCellTypesProcesses_PaperSelection_Astro = 0      # FEWER LISTS OF GENES ALL CELLS  
HistoAndFoldChangesGeneListsGeneralCellTypesProcesses_PaperSelection_OtherLists = 0 # 
DoGeneGeneCorrelation_Astro_Genes = 0                                               # Gene-Gene correlation among astrocytes specific genes, all cells
DoGeneGeneCorrelation_Astro_Genes_ASTRO_ONLY = 0                                               # Gene-Gene correlation among astrocytes specific genes, all cells

ComputeCorrelationTo_NR2F1_CRCgenes = 0                                             # Correlation each gene (in genome) vs NR2F1
ComputeIntersectionGenesLists_CORRELATION_NR2F1 = 0                                 # Intersect results

# Senescence Analysis
Use_Ohashi2018_UPDATED = True                                                       # Select updated genes list
HistoAndFoldChangesGeneLists_SenescenceGenes = 0                                    # Senescence Cumulatives
HistoAndFoldChangesIndividualGenes_Senescence = 0                                   # Senescence Individual Genes Fold Changes
HistoAndFoldChangesIndividualGenes_Senescence_ASTRO_ONLY = 0                        # Senescence Individual Genes Fold Changes
DoGeneGeneCorrelation_SenescenceGenes = 0                                           # Senescence Correlation Matrices
DoGeneGeneCorrelation_SenescenceGenes_ASTRO_ONLY = 0                                           # Senescence Correlation Matrices
DoGeneGeneCorrelation_Ohashi_Genes = 0                                              # Senescence Correlation Matrices only Ohashi List
DoGeneGeneCorrelation_Ohashi_Genes_ASTRO_ONLY = 0                                              # Senescence Correlation Matrices only Ohashi List
PlotGeneExpressionCOL1A1vsCOL14A1toSeeCorrelation = 0                               # Plot an example Gene-Gene Correlation
PlotGeneExpressionCOL1A1vsCOL3A1toSeeCorrelation = 0                                # Plot an example Gene-Gene Correlation

ExpressionNFIA = 0
ExpressionTGFBI = 0
Expression_NFIA_TGFBI = 0
Expression_NFIA_TGFBI_with_ZOOM = 'Yes'

ExpressionERN1 = 0

PlotSomeExpressionsForMasterEqProj = 1


######### START ##########

ListOfDays = [35,70]
    
N_cells = 500

L_CoreRegulatoryGenesLasse = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLasse/20180427-CRC-CoreTFs-unique-AllCellTypes_5_TOP_GENES.txt')[0]

L_stem = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0] 
L_CellCycle = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
L_ProApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Pro-Apoptotic.txt')[0]
L_AntiApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Anti-Apoptotic.txt')[0]
L_Caspases = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Caspases.txt')[0]
L_Astro = sd.get_from_txt('./ListsOfGenes/ListsGenesLisa/Lastro2-1_GFAP_S100b_Added.txt')[0]
L_Neurons_Lisa = sd.get_from_txt('./ListsOfGenes/ListsGenesLisa/EnrichedNeuronList.txt')[0]
L_DA_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisaNeuroSubtypes/L_DA_Neurons.txt')[0]

# L_Mito = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Mito.txt')[0]


ListsOfListsOfGenesCellTypesProcesses = [L_stem, L_Neurons_Lisa, L_Astro, L_CellCycle, L_DA_Neuro, L_ProApoptotic, L_AntiApoptotic, L_Caspases]

ListOfListsNames_CellTypesProcesses = ['List of\nGenes for\nStemness',
                        'List of\nGenes for\nNeurons',
                        'List of\nGenes for\nAstrocytes',
                        'List of\nGenes for\nCell Cycle',
                        'List of\nGenes for\nDopaminergic\nNeurons',
                        'List of\nGenes\nPro Apoptotic',
                        'List of\nGenes\nAnti Apoptotic',
                        'List of\nGenes\nCaspases']

ListOfListsNames_CellTypesProcesses_Short = ['Stemness',
                        'Neurons',
                        'Astrocytes',
                        'CellCycle',
                        'DANeurons',
                        'ProApoptotic',
                        'AntiApoptotic',
                        'Caspases']

if Use_Ohashi2018_UPDATED == False:
    L_Senescence_Ohashi_2018 = sd.get_from_txt('./ListsOfGenes/ListsOfGenesSenescenceSilvia/List1_SASP_Genes_Ohashi2018StemCellReports.txt')[0]
elif Use_Ohashi2018_UPDATED == True:
    L_Senescence_Ohashi_2018 = sd.get_from_txt('./ListsOfGenes/ListsOfGenesSenescenceSilvia/List1_SASP_Genes_Ohashi2018_UPDATED.txt')[0]
else:
    print("PROBLEM WITH LIST OHASHI!!!")
L_Senescence_Wileyetal_2017 = sd.get_from_txt('./ListsOfGenes/ListsOfGenesSenescenceSilvia/List2_SASP_genes_Wileyetal2017AgingCells.txt')[0]

ListOfListsOfGenesSenescence = [L_Senescence_Ohashi_2018, L_Senescence_Wileyetal_2017]
ListOfNamesOfListsOfGenesSenescence = ["Senescence genes Ohashi et al. 2018", 
                                       "Senescence genes Wiley et al. 2017"]

My_ylabel_Fold_Changes = r'Fold Change G2019S / WT'
# My_ylabel_Fold_Changes = r'Fold change ($log_2(\langle LRRK2G2019S \rangle / \langle LRRK2WT \rangle)$)'
            
if GetDataFromFiles == 1:
    print("GetDataFromFiles")
    # A = Mutant, B = Wild Tipe
    A_35_genes, A_35_expr, A_35_expr_np = sd.get_from_txt('./DataLisa/org4Mu35_S2_DGE.txt') 
    B_35_genes, B_35_expr, B_35_expr_np = sd.get_from_txt('./DataLisa/org3wt35_S1_DGE.txt')
    A_70_genes, A_70_expr, A_70_expr_np = sd.get_from_txt('./DataLisa/org2Mu70_S1_DGE.txt')
    B_70_genes, B_70_expr, B_70_expr_np = sd.get_from_txt('./DataLisa/org1wt70II_S1_DGE.txt')
    
    genesAlist = [A_35_genes,A_70_genes]
    genesBlist = [B_35_genes,B_70_genes]
    print("Done.")
    
    
if ComputeKneePlots == 1:
    mu_a35,mu_b35 = sd.get_knee_plot(A_35_expr_np)
    mu_a70,mu_b70 = sd.get_knee_plot(A_70_expr_np)
    
    wt_a35,wt_b35 = sd.get_knee_plot(B_35_expr_np)
    wt_a70,wt_b70 = sd.get_knee_plot(B_70_expr_np)
    print("Done.")

    
if PlotKneePlots == 1:
    plt.clf()    
    plt.xlim([0,3000]) 
    plt.xlabel('Cell Number') 
    plt.ylabel('Cumulative fraction of transcripts')
    plt.plot(mu_a35, label='Mu Day 35')
    plt.plot(mu_a70, label='Mu Day 70')  
    plt.plot(wt_a35, label='Wt Day 35')  
    plt.plot(wt_a70, label='Wt Day 70')
    plt.legend(loc='center right')
    plt.title('Knee Plot Cumulative')
    plt.savefig('./KneePlotCumulative.pdf')
    
    plt.clf() 
    plt.xlim([-10,3000]) 
    plt.xlabel('Cell Number') 
    plt.ylabel('Number of transcripts')  
    plt.plot(mu_b35, label='Mu Day 35') 
    plt.plot(mu_b70, label='Mu Day 70') 
    plt.plot(wt_b35, label='Wt Day 35') 
    plt.plot(wt_b70, label='Wt Day 70') 
    plt.legend(loc='center right')
    plt.title('Knee Plot Absolute')
    plt.savefig('./KneePlotAbsolute.pdf')    
    print("Done.")
    
    
if CutCellNumber == 1:
    print("CutCellNumber")
    A_35 = sd.get_top_cells_matrix(A_35_expr_np,N_cells)[0]
    B_35 = sd.get_top_cells_matrix(B_35_expr_np,N_cells)[0]
    A_70 = sd.get_top_cells_matrix(A_70_expr_np,N_cells)[0]
    B_70 = sd.get_top_cells_matrix(B_70_expr_np,N_cells)[0]
    print("Done.")

    
if Normalize == 1:
    print("Normalize")
    A_35_norm, A_35_log, A_35_sca = sd.get_normalization(A_35)
    B_35_norm, B_35_log, B_35_sca = sd.get_normalization(B_35)
    A_70_norm, A_70_log, A_70_sca = sd.get_normalization(A_70)
    B_70_norm, B_70_log, B_70_sca = sd.get_normalization(B_70)
    
    normAlist = [A_35_norm,A_70_norm]
    normBlist = [B_35_norm,B_70_norm]

    scaAlist = [A_35_sca,A_70_sca]
    scaBlist = [B_35_sca,B_70_sca]
    print("Done.")

    
if ComparePreliminaries == 1:
    print("ComparePreliminaries")
    C_35 = Compare((A_35_norm,A_35_genes),(B_35_norm,B_35_genes)) 
    C_70 = Compare((A_70_norm,A_70_genes),(B_70_norm,B_70_genes)) 
    C_35_merged = C_35.merge()[0] # aux35
    C_70_merged = C_70.merge()[0]
    C_35_norm, C_35_log, C_35_sca = sd.get_normalization(C_35_merged)
    C_70_norm, C_70_log, C_70_sca = sd.get_normalization(C_70_merged)
    normClist = [C_35_norm, C_70_norm]
    scaClist = [C_35_sca, C_70_sca]
    Clist = [C_35, C_70]
    print("Done.")


if ProduceListsOfDifferentiallyExpressedGenes == True:
    print("ProduceListsOfDifferentiallyExpressedGenes")
    
    ListOfDictionariesOfDEGsWithinLists = []
    j = 0
    DEGs_AllGenome = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
        
        #Saves list of genes that are different between MU and CT
        aux = C.avg() 
        a = np.array(aux) #some transformations
        a = np.transpose(a) #some transformations
        N_Bonferroni_Correction = len(C_norm)
        Corrected_Pvalue = 0.05/N_Bonferroni_Correction
        L_min, L_min_stats = C.avg_list(Corrected_Pvalue, aux[3], a)
        sd.save_array(L_min_stats, './Results/DiffExprGenes_Day_' + str(DAY) + '_005_Bonf.csv')
    
        DictionaryOfGenesListsToCompareWithDEGs = {"Stemness" : L_stem,
                                                   "Neurons" : L_Neurons_Lisa,
                                                   "Astrocytes" : L_Astro,
                                                   "CRC" : L_CoreRegulatoryGenesLasse,
                                                   "Senescence_Ohashi_2018" : L_Senescence_Ohashi_2018, 
                                                   "Senescence_Wileyetal_2017" : L_Senescence_Wileyetal_2017}

        DictionaryOfDEGsWithinLists = {}
        for KeyThisList in DictionaryOfGenesListsToCompareWithDEGs.keys():
            GenesList = DictionaryOfGenesListsToCompareWithDEGs[KeyThisList]
            ListOfDEGsWithinCurrentList = []
            for GeneName in L_min:
                if GeneName in GenesList:
                    ListOfDEGsWithinCurrentList.append(GeneName)
            DictionaryOfDEGsWithinLists.update({KeyThisList:ListOfDEGsWithinCurrentList})
            sd.save_from_list(ListOfDEGsWithinCurrentList, './Results/DiffExprGenes_Day_' + str(DAY) + '_' + str(KeyThisList) + '_Genes.txt')

        ListOfDictionariesOfDEGsWithinLists.append(DictionaryOfDEGsWithinLists)
        DEGs_AllGenome.append(L_min)
        j = j + 1
    print("Done.")
    

if IdentifyGenesActuallypresentInData == 1:
    print("IdentifyGenesActuallypresentInData")
    
    MatrixOfListsOfLisaGenesActuallyPresentInData = []
    j = 0
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
        
        ListOfListsOfLisaGenesActuallyPresentInData = []
        ListGenesC = C.comm_genes() 
        i = 0
        LisOfListsOfGenes_NoAbsentGenes = []
        for KeyThisList in DictionaryOfGenesListsToCompareWithDEGs.keys():
            CurrentList = DictionaryOfGenesListsToCompareWithDEGs[KeyThisList]
            ListGenesPresentInDataset = []
            for GeneInMyList in CurrentList:
                if GeneInMyList in ListGenesC:
                    ListGenesPresentInDataset.append(GeneInMyList)

            print(' ')
            print('List of Genes Considered: ' + str(KeyThisList))
            print(len(CurrentList))
            print(len(ListGenesPresentInDataset))
            print("The following genes are absent from the common data:")
            print([item for item in CurrentList if item not in ListGenesC])
            ListOfListsOfLisaGenesActuallyPresentInData.append(ListGenesPresentInDataset)
            print("The following genes are DEGs within the current list of genes:")
            ListOfDictionariesOfDEGs_WithinLists = ListOfDictionariesOfDEGsWithinLists[j]
            ListOfDEGs_WithinLists = ListOfDictionariesOfDEGs_WithinLists[KeyThisList]
            print(ListOfDEGs_WithinLists)
            print(len(ListOfDEGs_WithinLists))
            print(len(ListOfDEGs_WithinLists) / len(ListOfListsOfLisaGenesActuallyPresentInData[i]))
            i = i+1
                        
        print(" ")
        print(len(ListGenesC))
        print(len(DEGs_AllGenome[j]))
        print(len(DEGs_AllGenome[j])/len(ListGenesC))
        print(" ")
        print(" ")

        MatrixOfListsOfLisaGenesActuallyPresentInData.append(ListOfListsOfLisaGenesActuallyPresentInData)

        j = j + 1
    print("Done.")
    

if DoComparisonDR == 1:
    print("DoComparisonDR")

    mpl.style.use('default')

    j = 0
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        # C = Clist[j]
        C_norm = normClist[j]
        C_sca = scaClist[j]

        Myt1 = 'COMPARISON,sca - PCA - Day ' + str(DAY)
        Myt2 = 'COMPARISON,sca - tSNE - Day ' + str(DAY)
        MySaveAs1 = './Results/COMPARISON_PCA_sca'+ '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
        MySaveAs2 = './Results/COMPARISON_tSNE_sca'+ '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
        a,b = get_DR_bis(C_sca[:,0:N_cells],C_sca[:,N_cells:N_cells+N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

        Myt1 = 'COMPARISON,norm - PCA - Day ' + str(DAY)
        Myt2 = 'COMPARISON,norm - tSNE - Day ' + str(DAY)
        MySaveAs1 = './Results/COMPARISON_PCA_norm'+ '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
        MySaveAs2 = './Results/COMPARISON_tSNE_norm'+ '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
        a,b = get_DR_bis(C_norm[:,0:N_cells],C_norm[:,N_cells:N_cells+N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

        j = j + 1
        
#    Myt1 = 'COMPARISON,norm - PCA 3D'
#    Myt2 = 'COMPARISON,norm - tSNE 3D'
#    MySaveAs1 = './Results/COMPARISON_PCA_3D_norm'+ '_' + str(N_cells) + '.pdf'
#    MySaveAs2 = './Results/COMPARISON_tSNE_3D_norm'+ '_' + str(N_cells) + '.pdf'
#    a,b = get_DR_bis(C_35_70_norm[:,0:N_cells],C_35_70_norm[:,N_cells:N_cells+N_cells],com=3,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
    print("Done.")
    

if DoComparisonDR_ColoredByScores == 1:
    print("DoComparisonDR_ColoredByScores")
    
    j = 0
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = Clist[j]
        C_norm = normClist[j]
        C_sca = scaClist[j]
        
        k = 0
        for ListOfListsOfGenes in [ListsOfListsOfGenesCellTypesProcesses, [L_CoreRegulatoryGenesLasse]]:
        
            ListOfListsOfGenesForCumulatives = ListOfListsOfGenes
            ListOfCumulativesA = []
            ListOfCumulativesB = []
            
            if k == 0:
                ListNames_Short = ListOfListsNames_CellTypesProcesses_Short
            elif k == 1:
                ListNames_Short = L_CoreRegulatoryGenesLasse
        
            for CurrentList in ListOfListsOfGenesForCumulatives:
        
                CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                        C, C_norm, N_cells, CurrentList)
                
                ListOfCumulativesA.append(CumulativeGeneExpressionA)
                ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            i = 0
            for i in range(len(ListOfListsOfGenesForCumulatives)):
                
                Myt1 = 'COMPARISON,norm - PCA - Day '+ str(DAY)+'\npoints colored by expression of ' + str(ListNames_Short[i])
                Myt2 = 'COMPARISON,norm - tSNE - Day '+ str(DAY)+'\npoints colored by expression of ' + str(ListNames_Short[i])
                MySaveAs1 = './Results/COMPARISON_PCA_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
                MySaveAs2 = './Results/COMPARISON_tSNE_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
                ColorMap = 'rainbow'
                CumulativeGeneExpressionOneList = np.concatenate((ListOfCumulativesA[i], ListOfCumulativesB[i]), axis=0)
                a,b = get_DR_Colored_By_User_Provided_Vector(C_norm,CumulativeGeneExpressionOneList,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
    
                Myt1 = 'COMPARISON,sca - PCA - Day '+ str(DAY)+'\npoints colored by expression of ' + str(ListNames_Short[i])
                Myt2 = 'COMPARISON,sca - tSNE - Day '+ str(DAY)+'\npoints colored by expression of ' + str(ListNames_Short[i])
                MySaveAs1 = './Results/COMPARISON_PCA_sca_' + str(ListNames_Short[i]) + '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
                MySaveAs2 = './Results/COMPARISON_tSNE_sca_' + str(ListNames_Short[i]) + '_' + str(N_cells) + 'Day_' + str(DAY) + '.pdf'
                ColorMap = 'rainbow'
                CumulativeGeneExpressionOneList = np.concatenate((ListOfCumulativesA[i], ListOfCumulativesB[i]), axis=0)
                a,b = get_DR_Colored_By_User_Provided_Vector(C_sca,CumulativeGeneExpressionOneList,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
                
            k = k + 1
    
        ListOfCumulativesA = []
        ListOfCumulativesB = []
        
        j = j + 1
    print("Done.") 


if Compare_All_Wt_Mu_35_70_Preliminaries == 1:
    print("Compare_All_Wt_Mu_35_70_Preliminaries")
    
    C_ALL = Compare((C_35_norm,C_35.comm_genes()),(C_70_norm,C_70.comm_genes()))
    C_ALL_merged = C_ALL.merge()[0]
    C_ALL_norm, C_ALL_log, C_ALL_sca = sd.get_normalization(C_ALL_merged)
    
    if GenesToUseInALLComparison == "AllGenes":
        print("Ok - AllGenes")
        LabelWhichGenesDR = "_All_Genes_"
    elif GenesToUseInALLComparison == "SelectedGenes":
        print("Ok - SelectedGenes")
        LabelWhichGenesDR = "_Selected_Genes_"
        L_CellTypesGenes = L_stem + L_Neurons_Lisa + L_Astro
        C_ALL_norm = sd.get_submatrix(L_CellTypesGenes, C_ALL_norm, C_ALL.comm_genes())
        C_ALL_sca = sd.get_submatrix(L_CellTypesGenes, C_ALL_sca, C_ALL.comm_genes())
        C_ALL_log = sd.get_submatrix(L_CellTypesGenes, C_ALL_log, C_ALL.comm_genes())
    else:
        PRINT("ERROR!!! - IN DEFINING WHICH GENES SHOULD BE CONSIDERED IN DR")

    print("Done.")


if DoComparisonDR_ALL == 1:
    print("DoComparisonDR_ALL" + LabelWhichGenesDR)

    mpl.style.use('default')

    C_norm = C_ALL_norm
    C_sca = C_ALL_sca

    Myt1 = 'COMPARISON,sca - PCA - ALL conditions and days ' + LabelWhichGenesDR
    Myt2 = 'COMPARISON,sca - tSNE - ALL conditions and days ' + LabelWhichGenesDR
    MySaveAs1 = './Results/COMPARISON_PCA_sca'+ '_' + str(N_cells) + '_ALL_conditions_days' + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_tSNE_sca'+ '_' + str(N_cells) + '_ALL_conditions_days' + LabelWhichGenesDR + '.pdf'
    a,b = get_DR_bis(C_sca[:,0:N_cells],C_sca[:,N_cells:2*N_cells],C_sca[:,2*N_cells:3*N_cells],C_sca[:,3*N_cells:4*N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

    Myt1 = 'COMPARISON,norm - PCA - ALL conditions and days ' + LabelWhichGenesDR
    Myt2 = 'COMPARISON,norm - tSNE - ALL conditions and days ' + LabelWhichGenesDR
    MySaveAs1 = './Results/COMPARISON_PCA_norm'+ '_' + str(N_cells) + '_ALL_conditions_days' + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_tSNE_norm'+ '_' + str(N_cells) + '_ALL_conditions_days' + LabelWhichGenesDR + '.pdf'
    a,b = get_DR_bis(C_norm[:,0:N_cells],C_norm[:,N_cells:2*N_cells],C_norm[:,2*N_cells:3*N_cells],C_norm[:,3*N_cells:4*N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
    

if UsePCApreviousTotSNE_ALL == 1:

    if GenesToUseInALLComparison == "AllGenes":
        ThrasholdForVarianceExplained = 0.86 # This is the percentage (1 = 100%) of total explained variance of selected PCs
        MaxNumberOfPCs = 3
    elif GenesToUseInALLComparison == "SelectedGenes":
        ThrasholdForVarianceExplained = 0.90 # 0.86 # This is the percentage (1 = 100%) of total explained variance of selected PCs
        MaxNumberOfPCs = 12

    ##### Normalization: NORM #####
    M = C_ALL_norm

    ##### TAKE THIS OUT WHEN NOT NEEDED AS VERY TIME-CONSUMING
    y_pca_temp, IndexPCs = PCA_SelectingNumberOfPCsForGivenPercExplVarianceThrashold(M,ThrasholdForVarianceExplained, MaxNumberOfPCs) # 12) # 3)  

    NumberOfPCAbeforetSNE = IndexPCs + 1
    Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR
    Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '.pdf'

    a,b = get_DR_bis(np.transpose(y_pca_temp[0:N_cells,:]),np.transpose(y_pca_temp[N_cells:2*N_cells,:]),np.transpose(y_pca_temp[2*N_cells:3*N_cells,:]),np.transpose(y_pca_temp[3*N_cells:4*N_cells,:]),t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

    # cmap = matplotlib.cm.get_cmap(ColorMap)
    PatchMu35 = mpatches.Patch(color='red', label='G2019S d35')
    PatchCt35 = mpatches.Patch(color='blue', label='WT d35')
    PatchMu70 = mpatches.Patch(color='orange', label='G2019S d70')
    PatchCt70 = mpatches.Patch(color='cyan', label='WT d70')
    ListofPatches = [PatchMu35, PatchCt35, PatchMu70, PatchCt70]
    plt.legend(handles=ListofPatches, loc='best')
    
    NumberOfPCAbeforetSNE = IndexPCs + 1
    Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA 3D, norm - ALL conditions and days ' + LabelWhichGenesDR
    Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE 3D, norm - ALL conditions and days ' + LabelWhichGenesDR
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_3D_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_3D_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '.pdf'

    a,b = get_DR_bis(np.transpose(y_pca_temp[0:N_cells,:]),np.transpose(y_pca_temp[N_cells:2*N_cells,:]),np.transpose(y_pca_temp[2*N_cells:3*N_cells,:]),np.transpose(y_pca_temp[3*N_cells:4*N_cells,:]),com=3,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)
    
#    ##### Normalization: SCA #####
#    M = C_ALL_sca
#    Copy-paste here the same code as for NORM above and write sca instead of norm, if you wish to get the analysis for sca

    if UsePCApreviousTotSNE_ALL_AlsoColoredByScores == 1:
        k = 0
        Up_To_Which_lists = 3
        ListOfCumulatives_C_ALL = []
        for ListOfListsOfGenes in [ListsOfListsOfGenesCellTypesProcesses[0:Up_To_Which_lists]]: # 0:3 means [L_stem, L_Neurons_Lisa, L_Astro]:
        
            ListOfListsOfGenesForCumulatives = ListOfListsOfGenes
            ListOfCumulativesALL = []
            
            ListNames_Short = ListOfListsNames_CellTypesProcesses_Short[0:Up_To_Which_lists]

        
            for CurrentList in ListOfListsOfGenesForCumulatives:
        
                CumulativeGeneExpressionA_35, CumulativeGeneExpressionB_35 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                        C_ALL, C_ALL_norm[:,0:2*N_cells], N_cells, CurrentList)
                CumulativeGeneExpressionA_70, CumulativeGeneExpressionB_70 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                        C_ALL, C_ALL_norm[:,2*N_cells:4*N_cells], N_cells, CurrentList)
                
                ListOfCumulatives_C_ALL.append(np.concatenate((CumulativeGeneExpressionA_35,CumulativeGeneExpressionB_35,CumulativeGeneExpressionA_70,CumulativeGeneExpressionB_70), axis=0))
            
            i = 0
            for i in range(len(ListOfListsOfGenesForCumulatives)):
                
                Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(ListNames_Short[i])
                Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(ListNames_Short[i])
                MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + '_' + str(ListNames_Short[i]) + LabelWhichGenesDR + '.pdf'
                MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + '_' + str(ListNames_Short[i]) + LabelWhichGenesDR + '.pdf'

                ColorMap = 'rainbow'
                a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),ListOfCumulatives_C_ALL[i],ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
                
                Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA 3D, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(ListNames_Short[i])
                Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE 3D, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(ListNames_Short[i])
                MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_3D_norm'+ '_' + str(N_cells) + '_' + str(ListNames_Short[i]) + LabelWhichGenesDR + '.pdf'
                MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_3D_norm'+ '_' + str(N_cells) + '_' + str(ListNames_Short[i]) + LabelWhichGenesDR + '.pdf'

                ColorMap = 'rainbow'
                a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),ListOfCumulatives_C_ALL[i],ColorMap,3,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
                
            k = k + 1


if UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_NFIA == 1:
    GeneName = 'NFIA'
    #ExpressionOfNFIA = C_ALL_norm[C_ALL.c_names.index("NFIA"),:]
    GenesExpressionA_35 = C_35_norm[C_35.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_35 = C_35_norm[C_35.c_names.index(GeneName),N_cells:N_cells+N_cells]
    GenesExpressionA_70 = C_70_norm[C_70.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_70 = C_70_norm[C_70.c_names.index(GeneName),N_cells:N_cells+N_cells]
    Temp1 = np.concatenate((GenesExpressionA_35, GenesExpressionB_35), axis=None)
    Temp2 = np.concatenate((GenesExpressionA_70, GenesExpressionB_70), axis=None)

    ExpressionOfNFIA = np.concatenate((Temp1, Temp2), axis=None)
    
    Myt1 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of NFIA'
    Myt2 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of NFIA'
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + '_' + 'NFIA' + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + '_' + 'NFIA' + LabelWhichGenesDR + '.pdf'

    ColorMap = 'rainbow'
    a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),ExpressionOfNFIA,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)


if UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_TGFBI == 1:
    GeneName = 'TGFBI'
    #ExpressionOfNFIA = C_ALL_norm[C_ALL.c_names.index("NFIA"),:]
    GenesExpressionA_35 = C_35_norm[C_35.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_35 = C_35_norm[C_35.c_names.index(GeneName),N_cells:N_cells+N_cells]
    GenesExpressionA_70 = C_70_norm[C_70.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_70 = C_70_norm[C_70.c_names.index(GeneName),N_cells:N_cells+N_cells]
    Temp1 = np.concatenate((GenesExpressionA_35, GenesExpressionB_35), axis=None)
    Temp2 = np.concatenate((GenesExpressionA_70, GenesExpressionB_70), axis=None)

    ExpressionOfNFIA = np.concatenate((Temp1, Temp2), axis=None)
    
    Myt1 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(GeneName)
    Myt2 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(GeneName)
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + '_' + str(GeneName) + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + '_' + str(GeneName) + LabelWhichGenesDR + '.pdf'

    ColorMap = 'rainbow'
    a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),ExpressionOfNFIA,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
      

if UsePCApreviousTotSNE_ALL_AlsoColoredBy_ExpressionOf_NR2F1 == 1:
    GeneName = 'SOX2'
    #ExpressionOfNFIA = C_ALL_norm[C_ALL.c_names.index("NFIA"),:]
    GenesExpressionA_35 = C_35_norm[C_35.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_35 = C_35_norm[C_35.c_names.index(GeneName),N_cells:N_cells+N_cells]
    GenesExpressionA_70 = C_70_norm[C_70.c_names.index(GeneName),0:N_cells]
    GenesExpressionB_70 = C_70_norm[C_70.c_names.index(GeneName),N_cells:N_cells+N_cells]
    Temp1 = np.concatenate((GenesExpressionA_35, GenesExpressionB_35), axis=None)
    Temp2 = np.concatenate((GenesExpressionA_70, GenesExpressionB_70), axis=None)

    ExpressionOfNFIA = np.concatenate((Temp1, Temp2), axis=None)
    
    Myt1 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(GeneName)
    Myt2 = ' ' # 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by expression of ' + str(GeneName)
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + '_' + str(GeneName) + LabelWhichGenesDR + '.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + '_' + str(GeneName) + LabelWhichGenesDR + '.pdf'

    ColorMap = 'rainbow'
    a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),ExpressionOfNFIA,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)
        

if DetermineOptimalNumberOfClustersBySilhouetteAnalysis == True:
    
    if ClusteringOnFullDEMorPCs == "AllDEM":
        DataMatrixForDR = C_ALL_norm
    elif ClusteringOnFullDEMorPCs == "PCs": # "AllDEM"
        DataMatrixForDR = np.transpose(y_pca_temp)

    SilhouetteScores = []
    CellLabelsInClusters = []
    ClustersNumber = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    for N_Clusters in ClustersNumber: 
        ### RUN KMEANS with each clusters number ###
        kmeans = KMeans(n_clusters=N_Clusters, random_state=0).fit(np.transpose(DataMatrixForDR))
        CellLabelsInClusters.append(kmeans)
        MySilhouette = silhouette_score(np.transpose(DataMatrixForDR), kmeans.labels_, metric='euclidean')  
        SilhouetteScores.append(MySilhouette)
        print("\n" + "Kmeans using " + str(N_Clusters) + " clusters gives a silhouette score of: " + str(MySilhouette))
        print("the numbers of cells per cluster were:")
        for j in range(N_Clusters):
            print(list(kmeans.labels_).count(j))
    
    fig, ax = plt.subplots()
    ax.plot(ClustersNumber, SilhouetteScores,'-o')
    ax.set(xlabel='Clusters Number', ylabel='Average Silhouette Score',
    title='Identification of Optimal Cluster Number by Silhouette Analysis - ' + LabelWhichGenesDR)
    fig.savefig('SilhouetteAnalysisOnNumberOfClustersForKMeans' + LabelWhichGenesDR + '.pdf')
    plt.show()


if DoKmeansClusteringOnHDdata_PlotItAgainstPCApreviousTotSNE_ALL == True:
    
    if ClusteringOnFullDEMorPCs == "AllDEM":
        DataMatrixForDR = C_ALL_norm
    elif ClusteringOnFullDEMorPCs == "PCs": # "AllDEM"
        DataMatrixForDR = np.transpose(y_pca_temp)
        
    if GenesToUseInALLComparison == "AllGenes":
        N_Clusters = 4
    elif GenesToUseInALLComparison == "SelectedGenes":
        N_Clusters = 9 # 6
        
    ### RUN KMEANS with selected Clusters Number ###
    kmeans = KMeans(n_clusters=N_Clusters, random_state=0).fit(np.transpose(DataMatrixForDR))
    print("\n" + "Kmeans using " + str(N_Clusters))

    Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by Kmeans Clustering'
    Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by Kmeans Clustering'
    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '_Kmeans_Clustering_' + str(N_Clusters) + '_clusters.pdf'
    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '_Kmeans_Clustering_' + str(N_Clusters) + '_clusters.pdf'

    ColorMap = 'rainbow'
    a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),kmeans.labels_,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2,InsertColorBar=False)

    cmap = matplotlib.cm.get_cmap(ColorMap)
    ListofPatches = []
    for i in range(N_Clusters):
        print(list(kmeans.labels_).count(i))
        CurrentClustersPatch = mpatches.Patch(color=cmap(i/(N_Clusters-1)), label='C' + str(i) + ', ' + str(list(kmeans.labels_).count(i)) + ' cells')
        ListofPatches.append(CurrentClustersPatch)
    plt.legend(handles=ListofPatches, loc='best')
    print(' ')
    print(list(kmeans.labels_[0:500]).count(2))
    print(list(kmeans.labels_[500:1000]).count(2))
    print(list(kmeans.labels_[1000:1500]).count(2))
    print(list(kmeans.labels_[1500:2000]).count(2))


#    Myt1 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before PCA 3D, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by Kmeans Clustering'
#    Myt2 = 'COMPARISON,' + str(NumberOfPCAbeforetSNE) + ' PCA before tSNE 3D, norm - ALL conditions and days ' + LabelWhichGenesDR +'\npoints colored by Kmeans Clustering'
#    MySaveAs1 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_PCA_3D_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '_Kmeans_Clustering_' + str(N_Clusters) + '_clusters.pdf'
#    MySaveAs2 = './Results/COMPARISON_' + str(NumberOfPCAbeforetSNE) + '_PCA_before_tSNE_3D_norm'+ '_' + str(N_cells) + LabelWhichGenesDR + '_Kmeans_Clustering_' + str(N_Clusters) + '_clusters.pdf'
#
#    ColorMap = 'rainbow'
#    a,b = get_DR_Colored_By_User_Provided_Vector(np.transpose(y_pca_temp),kmeans.labels_,ColorMap,3,[],Myt1,Myt2,MySaveAs1,MySaveAs2,InsertColorBar=False)

    
#if ClustersIdentificationByVisualizationOfMarkersGeneExpression == True:
#    for 
    
                
if Plot3d_Astro_Neurons_Stemness_SCORES == 1:
    print("Plot3d_Astro_Neurons_Stemness_SCORES")

    CellTypesDictionary = {"Stem_Cells" : 0, "Neurons" : 1, "Astrocytes" : 2}
                           
    MatrixOfMU_Selected_Scores = []
    MatrixOfWT_Selected_Scores = []
    for DAY in ListOfDays:
        DAY_index = 0 if DAY is 35 else 1
       
        CurrentDataMU_Scores = []
        CurrentDataWT_Scores = []
        i = 0
        for List in ListsOfListsOfGenesCellTypesProcesses[0:3]: # Should be L_stem, L_Neurons_Lisa, L_Astro
            # THE NORMALIZATION TO EACH COUPLE (35/70) WAS PERFORMED PREVIOUSLY BUT I WOULD RATHER TAKE IT OUAT AS MEANINGLESS HERE
            # MaxXXX = max(max(MatrixOfCumulatives_A[DAY_index][i]),max(MatrixOfCumulatives_B[DAY_index][i]))
            CurrentDataMU_Scores.append(MatrixOfCumulatives_A[DAY_index][i]) # /MaxXXX)   # /len(List))
            CurrentDataWT_Scores.append(MatrixOfCumulatives_B[DAY_index][i]) # /MaxXXX)   # /len(List))
            i = i+ 1
        MatrixOfMU_Selected_Scores.append(CurrentDataMU_Scores)
        MatrixOfWT_Selected_Scores.append(CurrentDataWT_Scores)

    xdata_MU35 = MatrixOfMU_Selected_Scores[0][CellTypesDictionary["Stem_Cells"]]
    ydata_MU35 = MatrixOfMU_Selected_Scores[0][CellTypesDictionary["Neurons"]]
    zdata_MU35 = MatrixOfMU_Selected_Scores[0][CellTypesDictionary["Astrocytes"]]

    xdata_WT35 = MatrixOfWT_Selected_Scores[0][CellTypesDictionary["Stem_Cells"]]
    ydata_WT35 = MatrixOfWT_Selected_Scores[0][CellTypesDictionary["Neurons"]]
    zdata_WT35 = MatrixOfWT_Selected_Scores[0][CellTypesDictionary["Astrocytes"]]
    
    xdata_MU70 = MatrixOfMU_Selected_Scores[1][CellTypesDictionary["Stem_Cells"]]
    ydata_MU70 = MatrixOfMU_Selected_Scores[1][CellTypesDictionary["Neurons"]]
    zdata_MU70 = MatrixOfMU_Selected_Scores[1][CellTypesDictionary["Astrocytes"]]
    
    xdata_WT70 = MatrixOfWT_Selected_Scores[1][CellTypesDictionary["Stem_Cells"]]
    ydata_WT70 = MatrixOfWT_Selected_Scores[1][CellTypesDictionary["Neurons"]]
    zdata_WT70 = MatrixOfWT_Selected_Scores[1][CellTypesDictionary["Astrocytes"]]
    
    MaxX = max(max(xdata_MU35), max(xdata_WT35), max(xdata_MU70), max(xdata_WT70))
    MaxY = max(max(ydata_MU35), max(ydata_WT35), max(ydata_MU70), max(ydata_WT70))
    MaxZ = max(max(zdata_MU35), max(zdata_WT35), max(zdata_MU70), max(zdata_WT70))
    
    fig1 = plt.figure(figsize=plt.figaspect(0.5), facecolor="white")
    ax0 = fig1.add_subplot(2, 2, 1, projection='3d')

    ax0.set_xlabel('Stem Cells Cumulative Score (norm)')
    ax0.set_ylabel('Neurons Cumulative Score (norm)');
    ax0.set_zlabel('Astrocytes Cumulative Score (norm)')

    MyColors = ['red', 'blue', 'orange', 'cyan']

    ax0.scatter3D(xdata_MU35/MaxX, ydata_MU35/MaxY, zdata_MU35/MaxZ, c=MyColors[0]);
    ax0.scatter3D(xdata_WT35/MaxX, ydata_WT35/MaxY, zdata_WT35/MaxZ, c=MyColors[1]);
    ax0.scatter3D(xdata_MU70/MaxX, ydata_MU70/MaxY, zdata_MU70/MaxZ, c=MyColors[2]);
    ax0.scatter3D(xdata_WT70/MaxX, ydata_WT70/MaxY, zdata_WT70/MaxZ, c=MyColors[3]);

    binwidth = 0.0125 # 0.005  
    
    ax1= fig1.add_subplot(2, 2, 2)
    ax1.set_xlabel('Stem Cells Cumulative Score (norm)')
    ax1.set_ylabel('Number of Cells')
    ax1.hist([xdata_MU35/MaxX,xdata_WT35/MaxX,xdata_MU70/MaxX,xdata_WT70/MaxX],
             bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', 
             normed=0, stacked=True, color = MyColors) 
    
    ax2= fig1.add_subplot(2, 2, 3)
    ax2.set_xlabel('Neurons Cumulative Score (norm)')
    ax2.set_ylabel('Number of Cells')
    ax2.hist([ydata_MU35/MaxY,ydata_WT35/MaxY,ydata_MU70/MaxY,ydata_WT70/MaxY],
             bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', 
             normed=0, stacked=True, color = MyColors)
    
    ax3= fig1.add_subplot(2, 2, 4)
    ax3.set_xlabel('Astrocytes Cumulative Score (norm)')
    ax3.set_ylabel('Number of Cells')
    ax3.hist([zdata_MU35/MaxZ,zdata_WT35/MaxZ,zdata_MU70/MaxZ,zdata_WT70/MaxZ],
             bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', 
             normed=0, stacked=True, color = MyColors, label = ("LRRK2G2019S Day 35","LRRK2WT Day 35","LRRK2G2019S Day 70","LRRK2WT Day 70"))    

    plt.legend()
    plt.show()

    fig2 = plt.figure(figsize=plt.figaspect(0.5), facecolor="white")
    ax0 = fig2.add_subplot(2, 2, 1, projection='3d')

    ax0.set_xlabel('Stem Cells Cumulative Score (norm)')
    ax0.set_ylabel('Neurons Cumulative Score (norm)');
    ax0.set_zlabel('Astrocytes Cumulative Score (norm)')

    ax0.scatter3D(xdata_MU35/MaxX, ydata_MU35/MaxY, zdata_MU35/MaxZ, c=MyColors[0]);
    ax0.scatter3D(xdata_WT35/MaxX, ydata_WT35/MaxY, zdata_WT35/MaxZ, c=MyColors[1]);
    ax0.scatter3D(xdata_MU70/MaxX, ydata_MU70/MaxY, zdata_MU70/MaxZ, c=MyColors[2]);
    ax0.scatter3D(xdata_WT70/MaxX, ydata_WT70/MaxY, zdata_WT70/MaxZ, c=MyColors[3]);

    plt.legend()
    
    ax1= fig2.add_subplot(2, 2, 2)
    ax1.set_xlabel('Stem Cells Score Norm')
    ax1.set_ylabel('Neurons Cells Score Norm')
    ax1.hist2d(np.column_stack([xdata_MU35/MaxX,xdata_WT35/MaxX,xdata_MU70/MaxX,xdata_WT70/MaxX]).ravel(), 
              np.column_stack([ydata_MU35/MaxY,ydata_WT35/MaxY,ydata_MU70/MaxY,ydata_WT70/MaxY]).ravel(), (40, 40),
                   norm=colors.LogNorm(), cmap=plt.cm.jet)     

    ax2= fig2.add_subplot(2, 2, 3)
    ax2.set_xlabel('Stem Cells Score Norm')
    ax2.set_ylabel('Astrocytes Score Norm')
    ax2.hist2d(np.column_stack([xdata_MU35/MaxX,xdata_WT35/MaxX,xdata_MU70/MaxX,xdata_WT70/MaxX]).ravel(), 
              np.column_stack([zdata_MU35/MaxZ,zdata_WT35/MaxZ,zdata_MU70/MaxZ,zdata_WT70/MaxZ]).ravel(), (40, 40),
                   norm=colors.LogNorm(), cmap=plt.cm.jet)
    
    ax3= fig2.add_subplot(2, 2, 4)
    ax3.set_xlabel('Neurons Cells Score Norm')
    ax3.set_ylabel('Astrocytes Score Norm')
    ax3.hist2d(np.column_stack([ydata_MU35/MaxY,ydata_WT35/MaxY,ydata_MU70/MaxY,ydata_WT70/MaxY]).ravel(),
                np.column_stack([zdata_MU35/MaxZ,zdata_WT35/MaxZ,zdata_MU70/MaxZ,zdata_WT70/MaxZ]).ravel(), (40, 40),
                   norm=colors.LogNorm(), cmap=plt.cm.jet)
        
    plt.show()
    print("Done.")
    

if Separate_Cells_By_Cell_Type_From_3D_Plot3d == 1:
    print("Separate_Cells_By_Cell_Type_From_3D_Plot3d")
    
    Lists_Stemness_Scores = [xdata_MU35/MaxX, xdata_WT35/MaxX, xdata_MU70/MaxX, xdata_WT70/MaxX]
    Lists_Neurons_Scores = [ydata_MU35/MaxY, ydata_WT35/MaxY, ydata_MU70/MaxY, ydata_WT70/MaxY]
    Lists_Astro_Scores = [zdata_MU35/MaxZ, zdata_WT35/MaxZ, zdata_MU70/MaxZ, zdata_WT70/MaxZ]

    Num = (sum(xdata_MU35)+sum(xdata_WT35)+sum(xdata_MU70)+sum(xdata_WT70)) / MaxX
    Denom = len(xdata_MU35)+len(xdata_WT35)+len(xdata_MU70)+len(xdata_WT70)
    Mean_Stem_Score = Num / Denom
    
    Num = (sum(ydata_MU35)+sum(ydata_WT35)+sum(ydata_MU70)+sum(ydata_WT70)) / MaxY
    Denom = len(ydata_MU35)+len(ydata_WT35)+len(ydata_MU70)+len(ydata_WT70)
    Mean_Neuro_Score = Num / Denom
    
    Num = (sum(zdata_MU35/MaxZ)+sum(zdata_WT35/MaxZ)+sum(zdata_MU70/MaxZ)+sum(zdata_WT70/MaxZ))
    Denom = len(zdata_MU35/MaxZ)+len(zdata_WT35/MaxZ)+len(zdata_MU70/MaxZ)+len(zdata_WT70/MaxZ)
    Mean_Astro_Score = Num / Denom

    Lists_Stemness_Scores_Flattened = [Value for List in Lists_Stemness_Scores for Value in List]
    Lists_Neurons_Scores_Flattened = [Value for List in Lists_Neurons_Scores for Value in List]
    Lists_Astro_Scores_Flattened = [Value for List in Lists_Astro_Scores for Value in List]
    
    Median_Stemn_Score = statistics.median(Lists_Stemness_Scores_Flattened)
    Median_Neuro_Score = statistics.median(Lists_Neurons_Scores_Flattened)
    Median_Astro_Score = statistics.median(Lists_Astro_Scores_Flattened)
    
    CellTypeOfInterest = 2 # = ASTROCYTES

    ListsOfCellTypes_Based3dPlot = []
    ListOfIndexesOfCurrentCellTypeOfInterest = []
    index = 0
    k = 0
    for DAY in [35, 70]:
        for Treatment in ['MU', 'WT']:
            List_Stemn = Lists_Stemness_Scores[index]
            List_Astro = Lists_Astro_Scores[index]
            List_Neuro = Lists_Neurons_Scores[index]
                    
            if Threshold_Celltype == "Mean":
                ThSep_Stemn = Mean_Stem_Score # / 2
                ThSep_Neuro = Mean_Neuro_Score # / 2 
                ThSep_Astro = Mean_Astro_Score # / 2
            elif Threshold_Celltype == "Median":
                ThSep_Stemn = Median_Stemn_Score
                ThSep_Neuro = Median_Neuro_Score
                ThSep_Astro = Median_Astro_Score
            elif Threshold_Celltype == "Max":
                print("The cell type will be attributed based on the maximum score. No threshold needed.")
            else:
                print("ERROR in Threshold_Celltype!!!!!!!")
            
            CellTypes = []
            if Threshold_Celltype == "Mean" or Threshold_Celltype == "Median":
                for index_cell in range(N_cells):
                    if List_Stemn[index_cell] > ThSep_Stemn and List_Astro[index_cell] <= ThSep_Astro and List_Neuro[index_cell] <= ThSep_Neuro:
                        CellTypes.append(1)
                    elif List_Stemn[index_cell] <= ThSep_Stemn and List_Astro[index_cell] > ThSep_Astro and List_Neuro[index_cell] <= ThSep_Neuro:
                        CellTypes.append(2)
                    elif List_Stemn[index_cell] <= ThSep_Stemn and List_Astro[index_cell] <= ThSep_Astro and List_Neuro[index_cell] > ThSep_Neuro:
                        CellTypes.append(3)
                    else:
                        CellTypes.append(0)
            elif Threshold_Celltype == "Max":
                for index_cell in range(N_cells):
                    if List_Stemn[index_cell] > List_Astro[index_cell] and List_Stemn[index_cell] > List_Neuro[index_cell]:
                        CellTypes.append(1)
                    elif List_Astro[index_cell] > List_Stemn[index_cell] and List_Astro[index_cell] > List_Neuro[index_cell]:
                        CellTypes.append(2)
                    elif List_Neuro[index_cell] > List_Stemn[index_cell] and List_Neuro[index_cell] > List_Astro[index_cell] :
                        CellTypes.append(3)
                    else:
                        CellTypes.append(0)            
            else:
                print("ERROR in Threshold_Celltype!!!!!!!")
                              
            ListsOfCellTypes_Based3dPlot.append(CellTypes)
            
            temp = (np.array(CellTypes) == CellTypeOfInterest)
            temp2 = np.array(range(len(temp)))
            IndexesOfcurrentCellTypeOfInterest = list(temp2[temp])
            ListOfIndexesOfCurrentCellTypeOfInterest.append(IndexesOfcurrentCellTypeOfInterest)
            
            print(' ')
            print('PLEASE NOTICE THAT THE CLASSIFICATION BELOW IS HIGHLY ARBITRARY')
            print('At Day ' + str(DAY) + ', in the ' + str(Treatment) + ', using ' + Threshold_Celltype + ' as threshold/criterion, there are:')
            print('Stemm Cells: ' + str(CellTypes.count(1)))
            print('Astrocytes:  ' + str(CellTypes.count(2)))
            print('Neurons:     ' + str(CellTypes.count(3)))
            print('Unidentified:' + str(CellTypes.count(0)))
            
            index = index+1
    print("Done.")
    

if HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES == 1:
    print("HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES")

    CurrentList = L_CoreRegulatoryGenesLasse
    NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
        
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfDistributions_A = []
    MatrixOfDistributions_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

        MatrixOfDistributions_A.append(ListOfDistributions_A)
        MatrixOfDistributions_B.append(ListOfDistributions_B)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j=j+1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    LabelUpBy = [0,0,0,0,0]
    LabelDownBy = [0.25,0.25,0.25,0.25,0.25]
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])
                                              # BarsColors = "Bicolor")

    ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
    Nlines = len(CurrentList) 
    Ncols = len(ListOfDays)
    ylimVector = [N_cells, N_cells, N_cells, N_cells, N_cells]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A, MatrixOfDistributions_B, 
                              CurrentList, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if HistoAndFoldChangesGeneListsGeneralCellTypesProcesses == True:
    print("HistoAndFoldChangesGeneListsGeneralCellTypesProcesses")

    ListOfListsOfGenesForCumulatives = ListsOfListsOfGenesCellTypesProcesses

    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
    
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfCumulatives_A = []
    MatrixOfCumulatives_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        ListOfCumulativesA = []
        ListOfCumulativesB = []
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C, C_norm, N_cells, CurrentList)
            
            ListOfCumulativesA.append(CumulativeGeneExpressionA)
            ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    CumulativeGeneExpressionA,CumulativeGeneExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA, CumulativeGeneExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
        MatrixOfCumulatives_A.append(ListOfCumulativesA)
        MatrixOfCumulatives_B.append(ListOfCumulativesB)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j = j + 1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    LabelUpBy = [0,0,0,0,0,0,0,0]
    LabelDownBy = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames_CellTypesProcesses_Short, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])
            
    ########## MAKE HISTOGRAMS ##########     
    Nlines = len(ListOfListsNames_CellTypesProcesses_Short) 
    Ncols = len(ListOfDays)
    ylimVector = [N_cells, N_cells, N_cells, N_cells, N_cells, N_cells, N_cells, N_cells]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A, MatrixOfCumulatives_B, 
                              ListOfListsNames_CellTypesProcesses_Short, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
       
    
if HistoAndFoldChangesGeneListsGeneralCellTypesProcesses_PaperSelection_Astro == True:
    print("HistoAndFoldChangesGeneListsGeneralCellTypesProcesse_PaperSelections_Astro")

    ListOfListsOfGenesForCumulatives = [L_Astro, L_Astro]
    ListOfListsNames = ['Astro','Astro ZOOM']
    
    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
    
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfCumulatives_A = []
    MatrixOfCumulatives_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        ListOfCumulativesA = []
        ListOfCumulativesB = []
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C, C_norm, N_cells, CurrentList)
            
            ListOfCumulativesA.append(CumulativeGeneExpressionA)
            ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    CumulativeGeneExpressionA,CumulativeGeneExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA, CumulativeGeneExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
        MatrixOfCumulatives_A.append(ListOfCumulativesA)
        MatrixOfCumulatives_B.append(ListOfCumulativesB)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j = j + 1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    LabelUpBy = [0,0,0,0,0,0,0,0]
    LabelDownBy = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])
            
    ########## MAKE HISTOGRAMS ##########     
    Nlines = len(ListOfListsNames) 
    Ncols = len(ListOfDays)
    MyBinWidth = 0.01
    ylimVector = [N_cells, 130]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A, MatrixOfCumulatives_B, 
                              ListOfListsNames, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2-G2019S','LRRK2-WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'], BinWidth = MyBinWidth)
    print("Done.")
    
    
if HistoAndFoldChangesGeneListsGeneralCellTypesProcesses_PaperSelection_OtherLists == True:
    print("HistoAndFoldChangesGeneListsGeneralCellTypesProcesse_PaperSelections_OtherLists")

    ListOfListsOfGenesForCumulatives = [L_stem, L_Neurons_Lisa, L_CellCycle, L_ProApoptotic, L_AntiApoptotic, L_Caspases]
    ListOfListsNames = ['Stemness',
                        'Neurons',
                        'CellCycle',
                        'ProApoptotic',
                        'AntiApoptotic',
                        'Caspases']
    
    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
    
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfCumulatives_A = []
    MatrixOfCumulatives_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        ListOfCumulativesA = []
        ListOfCumulativesB = []
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C, C_norm, N_cells, CurrentList)
            
            ListOfCumulativesA.append(CumulativeGeneExpressionA)
            ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    CumulativeGeneExpressionA,CumulativeGeneExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA, CumulativeGeneExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
        MatrixOfCumulatives_A.append(ListOfCumulativesA)
        MatrixOfCumulatives_B.append(ListOfCumulativesB)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j = j + 1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    LabelUpBy = [0,0,0,0,0,0,0,0]
    LabelDownBy = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['d35','d70'],
                                              BarsColors = "Bicolor",
                                              NamesOnEachRectangle = 'Yes')
            
    ########## MAKE HISTOGRAMS ##########     
    Nlines = len(ListOfListsNames) 
    Ncols = len(ListOfDays)
    MyBinWidth = 0.05
    ylimVector = [N_cells, N_cells, N_cells, N_cells, N_cells, N_cells]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A, MatrixOfCumulatives_B, 
                              ListOfListsNames, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'], BinWidth = MyBinWidth)
    print("Done.")
    
 
if HistoAndFoldchangesGeneListsGeneralCellTypesProcesses_ASTRO_ONLY == True:
    print("HistoAndFoldchangesGeneListsGeneralCellTypesProcesses_ASTRO_ONLY")
    
    ListOfListsOfGenesForCumulatives = [L_stem, L_Neurons_Lisa, L_CellCycle, L_ProApoptotic, L_AntiApoptotic, L_Caspases] # ListsOfListsOfGenesCellTypesProcesses
    
    j = 0
    MatrixOfFoldChanges_ASTRO = []
    MatrixOfSEMforFoldChanges_ASTRO = []
    MatrixOfPvalues_ASTRO = []
    MatrixOfCumulatives_A_ASTRO = []
    MatrixOfCumulatives_B_ASTRO = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
          
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        if DAY == 70 or SeparateAstrocytesUsingCumulativesOrClustering == "Cumulative":

            ListOfCumulativesA = []
            ListOfCumulativesB = []
            k = 0
            for CurrentList in ListOfListsOfGenesForCumulatives:
                print("I am here")
                if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
                    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
                    IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
                    IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
                elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
                    NBonferroni = len(ListOfListsOfGenesForCumulatives) # * len(ListOfDays)
                    IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
                    IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
                    IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
                    IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
                    IndexesAstro_B = IndexesAstro_D70_Wt - 3*N_cells*np.ones(len(IndexesAstro_D70_Wt))
                    IndexesAstro_A = IndexesAstro_A.astype(int)
                    IndexesAstro_B = IndexesAstro_B.astype(int)
                else:
                    print("Error in selecting which method for separating astrocytes")
                    
                CumulativeGeneExpressionA_ASTRO = MatrixOfCumulatives_A[j][k][IndexesAstro_A]
                CumulativeGeneExpressionB_ASTRO = MatrixOfCumulatives_B[j][k][IndexesAstro_B]
                
                ListOfCumulativesA.append(CumulativeGeneExpressionA_ASTRO)
                ListOfCumulativesB.append(CumulativeGeneExpressionB_ASTRO)
                
                FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                        CumulativeGeneExpressionA_ASTRO,CumulativeGeneExpressionB_ASTRO)
                
                ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
                print(ListOfFoldChangesThisDay)
                ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA_ASTRO, CumulativeGeneExpressionB_ASTRO)
        
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                ListOfPvalues.append(pvalue)
                k = k + 1
            MatrixOfCumulatives_A_ASTRO.append(ListOfCumulativesA)
            MatrixOfCumulatives_B_ASTRO.append(ListOfCumulativesB)
            MatrixOfFoldChanges_ASTRO.append(ListOfFoldChangesThisDay)
            MatrixOfSEMforFoldChanges_ASTRO.append(ListOfSTDEVforFoldChangesThisDay)
            MatrixOfPvalues_ASTRO.append(ListOfPvalues)
        j = j + 1
    
    ListOfListsNames_CellTypesProcesses_Short_PaperSelection = ['Stemness', 'Neurons', 'CellCycle', 'ProApoptotic', 'AntiApoptotic', 'Caspases']
    
    LabelUpBy = [0,0,0,0,0,0,0,0]
    LabelDownBy = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        Ncols = len(ListOfDays)
        PlotFoldChangesManyGenesOneComparison(ListOfListsNames_CellTypesProcesses_Short_PaperSelection, 
                                              MatrixOfFoldChanges_ASTRO, 
                                              MatrixOfSEMforFoldChanges_ASTRO, 
                                              MatrixOfPvalues_ASTRO, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])
        ########## MAKE HISTOGRAMS ##########     
        Nlines = len(ListOfListsNames_CellTypesProcesses_Short_PaperSelection) 
        Ncols = len(ListOfDays)
        ylimVector = [80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500)]
        PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A_ASTRO, MatrixOfCumulatives_B_ASTRO, 
                              ListOfListsNames_CellTypesProcesses_Short_PaperSelection, MatrixOfPvalues_ASTRO, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
        Ncols = 1
        PlotFoldChangesManyGenesOneComparison(ListOfListsNames_CellTypesProcesses_Short_PaperSelection, 
                                              MatrixOfFoldChanges_ASTRO[0], 
                                              MatrixOfSEMforFoldChanges_ASTRO[0], 
                                              MatrixOfPvalues_ASTRO[0], 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              [['Day 70']],
                                              BarsColors = "Bicolor")
        ########## MAKE HISTOGRAMS ##########     
        Nlines = len(ListOfListsNames_CellTypesProcesses_Short_PaperSelection) 
        Ncols = len(ListOfDays)
        Max = 40
        ylimVector = [Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500), Max*(N_cells/500)]
        PlotHistogramsCumulatives(Nlines, 1, MatrixOfCumulatives_A_ASTRO[0], MatrixOfCumulatives_B_ASTRO[0], 
                              ListOfListsNames_CellTypesProcesses_Short_PaperSelection, MatrixOfPvalues_ASTRO[0], NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])        
    else:       
        print("Break cycle for day " + str(DAY) + " and Astro separation method " + SeparateAstrocytesUsingCumulativesOrClustering)
        
    print("Done.")


if HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES_ASTRO_ONLY == 1:
    print("HistoAndFoldChangesIndividualGenesCoreRegulatCircuitLasse_5_TOP_GENES_ASTRO_ONLY")

    CurrentList = L_CoreRegulatoryGenesLasse
        
    j = 1
    MatrixOfFoldChanges_ASTRO_CRC = []
    MatrixOfSEMforFoldChanges_ASTRO_CRC = []
    MatrixOfPvalues_ASTRO_CRC = []
    MatrixOfDistributions_A_ASTRO_CRC = []
    MatrixOfDistributions_B_ASTRO_CRC = []
    DAY = 70
    print('Day ' + str(DAY))
    
    C = Clist[j]
    C_norm = normClist[j]
        
    ListOfFoldChangesThisDay = []
    ListOfSTDEVforFoldChangesThisDay = []
    ListOfPvalues = []
    ListOfDistributions_A = []
    ListOfDistributions_B = []
    for i in range(len(CurrentList)): 
        GeneName = CurrentList[i]

        if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
            NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
            IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
            IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
        elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
            NBonferroni = len(CurrentList) # * len(ListOfDays) because only day 70 is eventually used
            IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
            IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
            IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
            IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
            IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
            IndexesAstro_A = IndexesAstro_A.astype(int)
            IndexesAstro_B = IndexesAstro_B.astype(int)
        else:
            print("Error in selecting which method for separating astrocytes")
        
        GenesExpressionA = C_norm[C.c_names.index(GeneName),IndexesAstro_A] # 0:N_cells]
        GenesExpressionB = C_norm[C.c_names.index(GeneName),IndexesAstro_B] # N_cells:N_cells+N_cells]
        
        ListOfDistributions_A.append(GenesExpressionA)
        ListOfDistributions_B.append(GenesExpressionB)
        
        FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                GenesExpressionA,GenesExpressionB)
        
        ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
        ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
        
        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')
        ListOfPvalues.append(pvalue)

    MatrixOfDistributions_A_ASTRO_CRC.append(ListOfDistributions_A)
    MatrixOfDistributions_B_ASTRO_CRC.append(ListOfDistributions_B)
    MatrixOfFoldChanges_ASTRO_CRC.append(ListOfFoldChangesThisDay)
    MatrixOfSEMforFoldChanges_ASTRO_CRC.append(ListOfSTDEVforFoldChangesThisDay)
    MatrixOfPvalues_ASTRO_CRC.append(ListOfPvalues)

    if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        LabelUpBy = [0,0,0,0,0]
        LabelDownBy = [0.4,0.4,0.4,0.4,0.4]
        Ncols = len(ListOfDays)
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  MatrixOfFoldChanges_ASTRO_CRC, 
                                                  MatrixOfSEMforFoldChanges_ASTRO_CRC, 
                                                  MatrixOfPvalues_ASTRO_CRC, 
                                                  NBonferroni, 
                                                  My_ylabel_Fold_Changes,
                                                  LabelUpBy,
                                                  LabelDownBy,
                                                  Ncols,
                                                  ['Day 35','Day 70'])
        ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
        Nlines = len(CurrentList) 
        Ncols = len(ListOfDays)
        ylimVector = [80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500), 80*(N_cells/500)]
        PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A_ASTRO_CRC, MatrixOfDistributions_B_ASTRO_CRC, 
                                  CurrentList, MatrixOfPvalues_ASTRO_CRC, NBonferroni, 
                                  ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                                  ylimVector, font3s, Colors=['Red','Gray'])
        
    elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":  
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        LabelUpBy = [0,0,0,0,0]
        LabelDownBy = [0.4,0.4,0.4,0.4,0.4]
        Ncols = 1
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  MatrixOfFoldChanges_ASTRO_CRC[0], 
                                                  MatrixOfSEMforFoldChanges_ASTRO_CRC[0], 
                                                  MatrixOfPvalues_ASTRO_CRC[0], 
                                                  NBonferroni, 
                                                  My_ylabel_Fold_Changes,
                                                  LabelUpBy,
                                                  LabelDownBy,
                                                  Ncols,
                                                  ['Day 70'])
        ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
        Nlines = len(CurrentList) 
        Ncols = 1
        ylimVector = [50*(N_cells/500), 50*(N_cells/500), 50*(N_cells/500), 50*(N_cells/500), 50*(N_cells/500)]
        PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A_ASTRO_CRC[0], MatrixOfDistributions_B_ASTRO_CRC[0], 
                                  CurrentList, MatrixOfPvalues_ASTRO_CRC[0], NBonferroni, 
                                  ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                                  ylimVector, font3s, Colors=['Red','Gray'])
        
    print("Done.")
    
    
if HistoAndFoldChangesIndividualGenes_AstroList == 1:
    print("HistoAndFoldChangesIndividualGenes_AstroList")

    MatrixOfCoefOfVariation_A = []
    MatrixOfCoefOfVariation_B = []

    NBonferroni = len(L_Astro) * len(ListOfDays) # because 2 days
    CurrentListName = 'L_Astro'
    CurrentList = L_Astro

    j = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListGenesC = C.comm_genes() 
        ListGenesActuallyPresentInDataset = []
        for GeneInMyList in CurrentList:
            if GeneInMyList in ListGenesC:
                ListGenesActuallyPresentInDataset.append(GeneInMyList)
        CurrentList = ListGenesActuallyPresentInDataset
        
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
            
        j=j+1
                
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        LabelUpBy = 0 * np.ones(len(CurrentList))
        LabelDownBy = 0.5 * np.ones(len(CurrentList))
        Ncols = 1
        MyTitleThisFig = ' ' # 'Fold changes for ' + CurrentListName + ' at Day ' + str(DAY)
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                              ListOfFoldChangesThisDay, 
                                              ListOfSTDEVforFoldChangesThisDay, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              MyTitle = MyTitleThisFig,
                                              BarsColors = "Bicolor")            
    print("Done.")         


if HistoAndFoldChangesIndividualGenes_AstroList_ASTRO_ONLY == 1:
    print("HistoAndFoldChangesIndividualGenes_AstroList_ASTRO_ONLY")

    MatrixOfCoefOfVariation_A = []
    MatrixOfCoefOfVariation_B = []

    NBonferroni = len(L_Astro) * len(ListOfDays) # because 2 days
    CurrentListName = 'L_Astro'
    CurrentList = L_Astro

    j = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListGenesC = C.comm_genes() 
        ListGenesActuallyPresentInDataset = []
        for GeneInMyList in CurrentList:
            if GeneInMyList in ListGenesC:
                ListGenesActuallyPresentInDataset.append(GeneInMyList)
        CurrentList = ListGenesActuallyPresentInDataset
        
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
                IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
                IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
            elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
                IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
                IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
                IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
                IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
                IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
                IndexesAstro_A = IndexesAstro_A.astype(int)
                IndexesAstro_B = IndexesAstro_B.astype(int)
            else:
                print("Error in selecting which method for separating astrocytes")
            
            GenesExpressionA = C_norm[C.c_names.index(GeneName),IndexesAstro_A] # 0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),IndexesAstro_B] # N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
            
        j=j+1
        
        print("Please notice that the indexes for d70 have been used for both d35 and d70 data, but of course for the first this does not make sense and indeed we do not plot them")

        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
            LabelUpBy = 0 * np.ones(len(CurrentList))
            LabelDownBy = np.ones(len(CurrentList))
            Ncols = 1
            MyTitleThisFig = ' ' # 'Fold changes for ' + CurrentListName + ' at Day ' + str(DAY)
            PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  ListOfFoldChangesThisDay, 
                                                  ListOfSTDEVforFoldChangesThisDay, 
                                                  ListOfPvalues, 
                                                  NBonferroni, 
                                                  My_ylabel_Fold_Changes,
                                                  LabelUpBy,
                                                  LabelDownBy,
                                                  Ncols,
                                                  MyTitle = MyTitleThisFig) 
        elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":  
            if DAY == 70:
                LabelUpBy = 0 * np.ones(len(CurrentList))
                LabelDownBy = np.ones(len(CurrentList))
                Ncols = 1
                MyTitleThisFig = ' ' # 'Fold changes for ' + CurrentListName + ' at Day ' + str(DAY)
                PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                      ListOfFoldChangesThisDay, 
                                                      ListOfSTDEVforFoldChangesThisDay, 
                                                      ListOfPvalues, 
                                                      NBonferroni, 
                                                      My_ylabel_Fold_Changes,
                                                      LabelUpBy,
                                                      LabelDownBy,
                                                      Ncols,
                                                      MyTitle = MyTitleThisFig,
                                                      BarsColors = "Bicolor") 
    print("Done.")     


if DoGeneGeneCorrelation_Astro_Genes == 1:
    print("DoGeneGeneCorrelation_Astro_Genes")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = L_Astro
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[kkk]
        C_norm = normClist[kkk]

        MyAverageCorrVector = []
        for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
            
            ######### COMPUTE GENE GENE CORRELATION MATRIX #########
            MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                                C_norm,
                                                                List_Gene_Names_for_Correlation,
                                                                N_cells,
                                                                ThrasoldForPval,
                                                                RemoveNotSignificantCorr,
                                                                WhiteDiagonal,
                                                                SampleIndex)
            
            ######### PLOT GENE GENE CORRELATION MATRIX #########
            figX, ax = plt.subplots(facecolor="white")
            
            if DoubleFace == False:
                
                if WhiteDiagonal == False:
                    MyCorrelationMatrix = MyCorrelationMatrix_Full
                elif WhiteDiagonal == True:
                    MyCorrelationMatrix = MyCorrelationMatrix_Masked
                else:
                    print("Error in Masking Correlations!")
                
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = 0
                        if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
                            
            elif DoubleFace == True:
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                        if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
            
            ax.invert_yaxis()
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3s)
                    ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3s, rotation=90)
                    
            #Spacing between each line
            intervals = 1
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
            ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
            #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            cbar = figX.colorbar(heatmap)
            cbar.set_label('Pearson Correlation', rotation=90)
            if SampleIndex == 0:
                SampleName = "G2019S"
            elif SampleIndex == 1:
                SampleName = "GC"
            # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
            # ax.set_title("Gene-gene Pearson Correlation Coefficient\n" + SampleName + ' at Day ' + str(DAY))
            
            plt.show()

        kkk = kkk + 1
    print("Done.")


if DoGeneGeneCorrelation_Astro_Genes_ASTRO_ONLY == 1:
    print("DoGeneGeneCorrelation_Astro_Genes_ASTRO_ONLY")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = L_Astro
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 1
    DAY = 70
    
    C = Clist[kkk]
    C_norm = normClist[kkk]

    if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
        IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
        IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
    elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
        IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
        IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
        IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
        IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
        IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
        IndexesAstro_A = IndexesAstro_A.astype(int)
        IndexesAstro_B = IndexesAstro_B.astype(int)
    else:
        print("Error in selecting which method for separating astrocytes")
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
                                
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        if SampleIndex == 0:
            IndexesOfSelectedCells = IndexesAstro_A
        elif  SampleIndex == 1:
            IndexesOfSelectedCells = IndexesAstro_B
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOnSelectedCellsOfCompareObject(C,
                                                            C_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal,
                                                            IndexesOfSelectedCells)
        
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        figX, ax = plt.subplots(facecolor="white")
        
        if DoubleFace == False:
            
            if WhiteDiagonal == False:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == True:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
                    if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
                        
        elif DoubleFace == True:
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                    if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
        
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3s)
                ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3s, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
        ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
        #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = figX.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)
        if SampleIndex == 0:
            SampleName = "G2019S"
        elif SampleIndex == 1:
            SampleName = "GC"
        # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
        # ax.set_title("Gene-gene Pearson Correlation Coefficient\n" + SampleName + ' at Day ' + str(DAY))
        
        plt.show()
    print("Done.")


if ComputeCorrelationTo_NR2F1_CRCgenes == 1:
    
    List_1_Gene_Names_for_Correlation = L_CoreRegulatoryGenesLasse   
    
    N_1 = len(List_1_Gene_Names_for_Correlation)
    N_C_35 = len(C_35.comm_genes())
    N_C_70 = len(C_70.comm_genes())
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    N_Bonf_Corr = (N_1 * (N_C_35 + N_C_70)) * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr
    
    ListOfCorrelationMatrices_Full_Samples_Days = []
    ListOfPvalsMatrices_Samples_Days = []
    ListOfCorrelationMatrices_Masked_Samples_Days = []
    k = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[k]
        C_norm = normClist[k]
        
        ListGenesC = C.comm_genes()
        List_2_Gene_Names_for_Correlation = ListGenesC

        for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
            
            ######### COMPUTE GENE GENE CORRELATION MATRIX #########
            WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
            RemoveNotSignificantCorr = False
            MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelation_Between2DifferentGenesLists_On2SamplesOfCompareObject(C,
                                                                C_norm,
                                                                List_1_Gene_Names_for_Correlation,
                                                                List_2_Gene_Names_for_Correlation,
                                                                N_cells,
                                                                ThrasoldForPval,
                                                                RemoveNotSignificantCorr,
                                                                WhiteDiagonal,
                                                                SampleIndex)


            ListOfCorrelationMatrices_Full_Samples_Days.append(MyCorrelationMatrix_Full)
            ListOfPvalsMatrices_Samples_Days.append(pvalPearsonCorr_Matrix)
            ListOfCorrelationMatrices_Masked_Samples_Days.append(MyCorrelationMatrix_Masked)
        k = k + 1
    
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    MyColors = ["red","blue","green","cyan"]
    for index in range( len(ListOfDays) * Nsamples ):
        X_values = range(len(ListOfCorrelationMatrices_Full_Samples_Days[index][0])) # index 0 = Mu35, 1 = Wt35, 2 = Mu70, 3 = Wt70
        plt.scatter(X_values, np.sort(ListOfCorrelationMatrices_Full_Samples_Days[index][0]),c=MyColors[index], edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("Gene Number")
    ax.set_ylabel("Correlation with NR2F1")
    ax.set_title("Correlation with NR2F1 of each gene of the genome\n(nans have been removed, so less genes than in the full measured genome)")
    ax.set_xlim(-100,25000)
    ax.set_ylim(-0.35,1.05)
    ax.legend(["G2019S, Day 35","GC, Day 35", "G2019S, Day 70","GC, Day 70"],loc='upper left')    
    plt.show()
    print("Done.")

    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    MyColors = ["red","blue","green","cyan"]
    for index in range( len(ListOfDays) * Nsamples ):
        X_values = range(len(ListOfCorrelationMatrices_Masked_Samples_Days[index][0])) # index 0 = Mu35, 1 = Wt35, 2 = Mu70, 3 = Wt70
        plt.scatter(X_values, (ListOfCorrelationMatrices_Masked_Samples_Days[index][0]),c=MyColors[index], edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("Genes")
    ax.set_ylabel("Correlation with NR2F1")
    ax.set_title("Correlation with NR2F1 across all genes of the genome")
    ax.legend(["G2019S, Day 35","GC, Day 35", "G2019S, Day 70","GC, Day 70"],loc='upper left')    
    plt.show()

    index = 0
    k = 0
    ListOfOrderedGenesCorrelations = []
    for DAY in ListOfDays:
        C = Clist[k]        
        ListGenesC = C.comm_genes() 
        for SampleIndex in range(Nsamples):
            OrderedIndexes = np.argsort(ListOfCorrelationMatrices_Full_Samples_Days[index][0]) # 0 here and below stantds for NR2F1
            OrderedCorrelations_Full = [ListOfCorrelationMatrices_Full_Samples_Days[index][0][MyIndex] for MyIndex in OrderedIndexes]
            OrderedCorrelations_Masked = [ListOfCorrelationMatrices_Masked_Samples_Days[index][0][MyIndex] for MyIndex in OrderedIndexes]
            OrderedGeneNames = [ListGenesC[MyIndex] for MyIndex in OrderedIndexes]
            OrderedPvals = [ListOfPvalsMatrices_Samples_Days[index][0][MyIndex] for MyIndex in OrderedIndexes]

            OrderedStuffAsArray = np.array([OrderedGeneNames,OrderedCorrelations_Full,OrderedCorrelations_Masked,OrderedPvals])
            ListOfOrderedGenesCorrelations.append(OrderedStuffAsArray)
            index = index + 1
        k = k + 1
    
    
    ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations = []
    for SampleIndex in range(len(ListOfDays) * Nsamples):
        ListIndexes = []
        for GeneIndex in range(len(ListOfOrderedGenesCorrelations[SampleIndex][0])):
            if float(ListOfOrderedGenesCorrelations[SampleIndex][2,GeneIndex]) != 0:
                ListIndexes.append(GeneIndex)
        ArrayIndexes = np.array(ListIndexes)
        TEMP_Array = np.transpose(ListOfOrderedGenesCorrelations[SampleIndex][:,ArrayIndexes])
        TEMP_Array_reverted = TEMP_Array[::-1]
        
        ListOfDaysSamplesTags = ['Day_35_G2019S','Day_35_GC','Day_70_G2019S','Day_70_GC']
        
        sd.save_array(TEMP_Array_reverted, './Results/GenesCorrelatingWith_NR2F1_' + ListOfDaysSamplesTags[SampleIndex] + '_005_Bonf.csv')
        
        ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations.append(TEMP_Array_reverted)
    
    
if ComputeIntersectionGenesLists_CORRELATION_NR2F1 == 1:

    L_CorrNR2F1_day35_G2019S = ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations[0][:,0]
    L_CorrNR2F1_day35_GC = ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations[1][:,0]
    L_CorrNR2F1_day70_G2019S = ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations[2][:,0]
    L_CorrNR2F1_day70_GC = ListOfListsOfOrderedGenes_ONLY_SIGNIFICANT_Correlations[3][:,0]

    Intersect4listsOfGenePrintNumbers(L_CorrNR2F1_day35_G2019S, 
                                      L_CorrNR2F1_day35_GC, 
                                      L_CorrNR2F1_day70_G2019S, 
                                      L_CorrNR2F1_day70_GC)

    ListAllGenesCorrNR2F1_NoRepetitions = set(list(L_CorrNR2F1_day35_G2019S) + list(L_CorrNR2F1_day35_GC) + list(L_CorrNR2F1_day70_G2019S) + list(L_CorrNR2F1_day70_GC))
        
    Intersect4listsOfGenePrintNumbers(ListAllGenesCorrNR2F1_NoRepetitions, L_stem, L_Neurons_Lisa, L_Astro)
    
    C_norm = C_70_norm
            
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    Gene_1_Expression_A = C_norm[C.c_names.index('NR2F1'),0:N_cells]
    Gene_2_Expression_A = C_norm[C.c_names.index('SOX2'),0:N_cells]
    plt.scatter(Gene_1_Expression_A, Gene_2_Expression_A,c='red', edgecolors='none') # 0 here and above stantds for NR2F1
    Gene_1_Expression_B = C_norm[C.c_names.index('NR2F1'),0+N_cells:N_cells+N_cells]
    Gene_2_Expression_B = C_norm[C.c_names.index('SOX2'),0+N_cells:N_cells+N_cells]
    plt.scatter(Gene_1_Expression_B, Gene_2_Expression_B,c='blue', edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("NR2F1 Expression")
    ax.set_ylabel("SOX2 Expression")
    ax.set_title("Expression of SOX2 vs NR2F1 to visualize correlation, Day 70, G2019S and GC")
    ax.set_xlim(-1,61)
    ax.set_ylim(-1,61)
    ax.legend(["G2019S, Day 70","GC, Day 70"],loc='upper right')    
    plt.show()
    
                
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    Gene_1_Expression_A = C_norm[C.c_names.index('NR2F1'),0:N_cells]
    Gene_2_Expression_A = C_norm[C.c_names.index('SYT1'),0:N_cells]
    plt.scatter(Gene_1_Expression_A, Gene_2_Expression_A,c='red', edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("NR2F1 Expression")
    ax.set_ylabel("SYT1 Expression")
    ax.set_title("Expression of SYT1 vs NR2F1 to visualize correlation, Day 70, G2019S")
    ax.set_xlim(-1,170)
    ax.set_ylim(-1,170)
    ax.legend(["G2019S, Day 70"],loc='upper right')    
    plt.show()
    
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    Gene_1_Expression_B = C_norm[C.c_names.index('NR2F1'),0+N_cells:N_cells+N_cells]
    Gene_2_Expression_B = C_norm[C.c_names.index('FOXA2'),0+N_cells:N_cells+N_cells]
    plt.scatter(Gene_1_Expression_B, Gene_2_Expression_B,c='BLUE', edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("NR2F1 Expression")
    ax.set_ylabel("FOXA2 Expression")
    ax.set_title("Expression of FOXA2 vs NR2F1 to visualize correlation, Day 70, GC")
    ax.set_xlim(-1,60)
    ax.set_ylim(-1,60)
    ax.legend(["GC, Day 70"],loc='upper right')    
    plt.show()
    
    print("Done.")
    
    
if HistoAndFoldChangesGeneLists_SenescenceGenes == True:
    print("HistoAndFoldChangesGeneLists_SenescenceGenes")

    ListOfListsOfGenesForCumulatives = ListOfListsOfGenesSenescence

    NBonferroni = len(ListOfListsOfGenesForCumulatives) * len(ListOfDays)
    
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfCumulatives_A = []
    MatrixOfCumulatives_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        
        ListOfCumulativesA = []
        ListOfCumulativesB = []
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpressionA, CumulativeGeneExpressionB = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C, C_norm, N_cells, CurrentList)
            
            ListOfCumulativesA.append(CumulativeGeneExpressionA)
            ListOfCumulativesB.append(CumulativeGeneExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    CumulativeGeneExpressionA,CumulativeGeneExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpressionA, CumulativeGeneExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
        MatrixOfCumulatives_A.append(ListOfCumulativesA)
        MatrixOfCumulatives_B.append(ListOfCumulativesB)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j = j + 1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    ListOfListsNames_Senescence_Short = ['Ohashi2018','Wiley2017']
    LabelUpBy = [0,0]
    LabelDownBy = [0.1,0.1]
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames_Senescence_Short, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])
            
    ########## MAKE HISTOGRAMS ##########     
    Nlines = len(ListOfListsNames_Senescence_Short) 
    Ncols = len(ListOfDays)
    ylimVector = [500, 500]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfCumulatives_A, MatrixOfCumulatives_B, 
                              ListOfListsNames_Senescence_Short, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Cumulative Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if HistoAndFoldChangesIndividualGenes_Senescence == 1:
    print("HistoAndFoldChangesIndividualGenes_Senescence_Ohashi2018")

    MatrixOfCoefOfVariation_A = []
    MatrixOfCoefOfVariation_B = []

    NBonferroni = (len(ListOfListsOfGenesSenescence[0]) + len(ListOfListsOfGenesSenescence[1])) * len(ListOfDays) # because 2 days
    k = 0
    for GenesList in ListOfListsOfGenesSenescence:

        CurrentListName = ListOfNamesOfListsOfGenesSenescence[k]
        CurrentList = GenesList

        j = 0
        for DAY in ListOfDays:        
            print('Day ' + str(DAY))
            
            C = Clist[j]
            C_norm = normClist[j]
                
            ListGenesC = C.comm_genes() 
            ListGenesActuallyPresentInDataset = []
            for GeneInMyList in CurrentList:
                print(GeneInMyList)
                if GeneInMyList in ListGenesC:
                    print('IN')
                    ListGenesActuallyPresentInDataset.append(GeneInMyList)
            if DAY == 70 and CurrentListName == ListOfNamesOfListsOfGenesSenescence[0]:
                ListGenesActuallyPresentInDataset.append('CDKN2A')
            CurrentList = ListGenesActuallyPresentInDataset
            
            ListOfFoldChangesThisDay = []
            ListOfSTDEVforFoldChangesThisDay = []
            ListOfPvalues = []
            ListOfDistributions_A = []
            ListOfDistributions_B = []
            ListOfCoefOfVariation_A = []
            ListOfCoefOfVariation_B = []
            for i in range(len(CurrentList)): 
                GeneName = CurrentList[i]
                
                print(i)
                print(GeneName)
                
                GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
                GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
                
                ListOfDistributions_A.append(GenesExpressionA)
                ListOfDistributions_B.append(GenesExpressionB)
                
                ListOfCoefOfVariation_A.append(np.std(GenesExpressionA)/statistics.mean(GenesExpressionA))
                ListOfCoefOfVariation_B.append(np.std(GenesExpressionB)/statistics.mean(GenesExpressionB))
                
                FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                        GenesExpressionA,GenesExpressionB)
                
                ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
                ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
                
                ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
                # of 2 populations are statistically significantly different ######
                pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
        
                # Remember Bonferroni correction for multiple testing, 
                # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
                # by the number of repetitions of the test
                print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
                print(' ')
                ListOfPvalues.append(pvalue)
                
            MatrixOfCoefOfVariation_A.append(ListOfCoefOfVariation_A)
            MatrixOfCoefOfVariation_B.append(ListOfCoefOfVariation_B)
                
            j=j+1
            
            ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
            LabelUpBy = 0 * np.ones(len(CurrentList))
            LabelDownBy = np.ones(len(CurrentList))
            Ncols = 1
            MyTitleThisFig = ' ' # 'Fold changes for ' + CurrentListName + ' at Day ' + str(DAY)
            PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  ListOfFoldChangesThisDay, 
                                                  ListOfSTDEVforFoldChangesThisDay, 
                                                  ListOfPvalues, 
                                                  NBonferroni, 
                                                  My_ylabel_Fold_Changes,
                                                  LabelUpBy,
                                                  LabelDownBy,
                                                  Ncols,
                                                  MyTitle = MyTitleThisFig,
                                                  BarsColors = "Bicolor")
            

            
            fig = plt.figure(facecolor="white")
            ax = fig.add_subplot(111)
            ax.scatter(ListOfCoefOfVariation_A, ListOfCoefOfVariation_B, edgecolors='none')
            ax.set_xlabel("Coefficient of Variation LRRK2-G2019S")
            ax.set_ylabel("Coefficient Of Variation LRRK2-WT")
            # ax.set_title("Coefficients of Vairiaion WT vs MU for each gene\n" + CurrentListName + ' at Day ' + str(DAY))
            ax.set_xlim(0,25)
            ax.set_ylim(0,25)
            ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
            plt.show()            
            
        k = k + 1
    
    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    MyColors = ["red","blue","green","cyan"]
    for index in range( len(ListOfDays) * len(ListOfListsOfGenesSenescence) ):
        ax.scatter(MatrixOfCoefOfVariation_A[index], MatrixOfCoefOfVariation_B[index],c=MyColors[index], edgecolors='none')
    ax.set_xlabel("Coefficient of Variation LRRK2-G2019S")
    ax.set_ylabel("Coefficient Of Variation LRRK2-WT")
    # ax.set_title("Coefficients of Vairiaion WT vs MU for each gene")
    ax.set_xlim(0,25)
    ax.set_ylim(0,25)
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    ax.legend(["Equal CV Mu and Wt","Day 35, Ohashi et al. 2018","Day 70, Ohashi et al. 2018",
               "Day 35, Wiley et al. 2017","Day 70, Wiley et al. 2017"],loc='upper left')    
    plt.show()

#    ListOfRatioCVsMuCt = []
#    for index in range( len(ListOfDays) * len(ListOfListsOfGenesSenescence) ):
#        ListA = MatrixOfCoefOfVariation_A[index]
#        ListB = MatrixOfCoefOfVariation_B[index]
#        for i in range(len(ListA)):
#            RatioCVsMuCt = ListA[i] / ListB[i]
#            ListOfRatioCVsMuCt.append(RatioCVsMuCt)
#    ListOfRatioCVsMuCt_NansRemoved = []  
#    for j in range(len(ListOfRatioCVsMuCt)):
#        if not np.isnan(ListOfRatioCVsMuCt[j]):
#            ListOfRatioCVsMuCt_NansRemoved.append(ListOfRatioCVsMuCt[j])
#    AverageRatioCVsMuCt_Overall = statistics.mean(ListOfRatioCVsMuCt_NansRemoved)
#    print(AverageRatioCVsMuCt_Overall)
#    print("Done.")  
    
    ListOfRatioCVsMuCt = []
    for index in range( len(ListOfDays) * len(ListOfListsOfGenesSenescence) ):
        ListA = MatrixOfCoefOfVariation_A[index]
        ListB = MatrixOfCoefOfVariation_B[index]
        for i in range(len(ListA)):
            RatioCVsMuCt = ListA[i] / ListB[i]
            ListOfRatioCVsMuCt.append(RatioCVsMuCt)
        ListOfRatioCVsMuCt_NansRemoved = []  
        for j in range(len(ListOfRatioCVsMuCt)):
            if not np.isnan(ListOfRatioCVsMuCt[j]):
                ListOfRatioCVsMuCt_NansRemoved.append(ListOfRatioCVsMuCt[j])
        AverageRatioCVsMuCt_Overall = statistics.mean(ListOfRatioCVsMuCt_NansRemoved)
        print(AverageRatioCVsMuCt_Overall)
    print("Done.")        
    
    
if HistoAndFoldChangesIndividualGenes_Senescence_ASTRO_ONLY == 1:
    print("HistoAndFoldChangesIndividualGenes_Senescence_Ohashi2018")

    MatrixOfCoefOfVariation_A = []
    MatrixOfCoefOfVariation_B = []

    NBonferroni = (len(ListOfListsOfGenesSenescence[0]) + len(ListOfListsOfGenesSenescence[1])) # * len(ListOfDays) # because 1 day
    k = 0
    for GenesList in ListOfListsOfGenesSenescence:

        CurrentListName = ListOfNamesOfListsOfGenesSenescence[k]
        CurrentList = GenesList

        j = 1
        DAY = 70       
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListGenesC = C.comm_genes() 
        ListGenesActuallyPresentInDataset = []
        for GeneInMyList in CurrentList:
            print(GeneInMyList)
            if GeneInMyList in ListGenesC:
                print('IN')
                ListGenesActuallyPresentInDataset.append(GeneInMyList)
        if DAY == 70 and CurrentListName == ListOfNamesOfListsOfGenesSenescence[0]:
            ListGenesActuallyPresentInDataset.append('CDKN2A')
        CurrentList = ListGenesActuallyPresentInDataset
        
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        ListOfCoefOfVariation_A = []
        ListOfCoefOfVariation_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
            
            print(i)
            print(GeneName)
            
            if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
                IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
                IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
            elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
                IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
                IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
                IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
                IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
                IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
                IndexesAstro_A = IndexesAstro_A.astype(int)
                IndexesAstro_B = IndexesAstro_B.astype(int)
            else:
                print("Error in selecting which method for separating astrocytes")
            
            GenesExpressionA = C_norm[C.c_names.index(GeneName),IndexesAstro_A] # 0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),IndexesAstro_B] # N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            ListOfCoefOfVariation_A.append(np.std(GenesExpressionA)/statistics.mean(GenesExpressionA))
            ListOfCoefOfVariation_B.append(np.std(GenesExpressionB)/statistics.mean(GenesExpressionB))
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)
            
        MatrixOfCoefOfVariation_A.append(ListOfCoefOfVariation_A)
        MatrixOfCoefOfVariation_B.append(ListOfCoefOfVariation_B)
            
        j=j+1
        
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        LabelUpBy = 0 * np.ones(len(CurrentList))
        LabelDownBy = np.ones(len(CurrentList))
        Ncols = 1
        MyTitleThisFig = ' ' # 'Fold changes for ' + CurrentListName + ' at Day ' + str(DAY)
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                              ListOfFoldChangesThisDay, 
                                              ListOfSTDEVforFoldChangesThisDay, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              MyTitle = MyTitleThisFig,
                                              BarsColors = "Bicolor")
        

        
        fig = plt.figure(facecolor="white")
        ax = fig.add_subplot(111)
        ax.scatter(ListOfCoefOfVariation_A, ListOfCoefOfVariation_B, edgecolors='none')
        ax.set_xlabel("Coefficient of Variation LRRK2-G2019S")
        ax.set_ylabel("Coefficient Of Variation LRRK2-WT")
        # ax.set_title("Coefficients of Vairiaion WT vs MU for each gene\n" + CurrentListName + ' at Day ' + str(DAY))
        ax.set_xlim(0,25)
        ax.set_ylim(0,25)
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
        plt.show()            
            
        k = k + 1
#    
#    fig = plt.figure(facecolor="white")
#    ax = fig.add_subplot(111)
#    MyColors = ["red","blue","green","cyan"]
#    for index in range( len(ListOfDays) * len(ListOfListsOfGenesSenescence) ):
#        ax.scatter(MatrixOfCoefOfVariation_A[index], MatrixOfCoefOfVariation_B[index],c=MyColors[index], edgecolors='none')
#    ax.set_xlabel("Coefficient of Variation LRRK2G2019S")
#    ax.set_ylabel("Coefficient Of Variation LRRK2CT")
#    ax.set_title("Coefficients of Vairiaion WT vs MU for each gene")
#    ax.set_xlim(0,25)
#    ax.set_ylim(0,25)
#    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
#    ax.legend(["Equal CV Mu and Wt","Day 35, Ohashi et al. 2018","Day 70, Ohashi et al. 2018",
#               "Day 35, Wiley et al. 2017","Day 70, Wiley et al. 2017"],loc='upper left')    
#    plt.show()
#
    ListOfRatioCVsMuCt = []
    for index in range(len(ListOfListsOfGenesSenescence) ):
        ListA = MatrixOfCoefOfVariation_A[index]
        ListB = MatrixOfCoefOfVariation_B[index]
        for i in range(len(ListA)):
            RatioCVsMuCt = ListA[i] / ListB[i]
            ListOfRatioCVsMuCt.append(RatioCVsMuCt)
        ListOfRatioCVsMuCt_NansRemoved = []  
        for j in range(len(ListOfRatioCVsMuCt)):
            if not np.isnan(ListOfRatioCVsMuCt[j]):
                ListOfRatioCVsMuCt_NansRemoved.append(ListOfRatioCVsMuCt[j])
        AverageRatioCVsMuCt_Overall = statistics.mean(ListOfRatioCVsMuCt_NansRemoved)
        print(AverageRatioCVsMuCt_Overall)
        
    print("Done.")        
    
        
if DoGeneGeneCorrelation_SenescenceGenes == 1:
    print("DoGeneGeneCorrelation_SenescenceGenes")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    MyListOfLists = ListOfListsOfGenesSenescence
    
    List_1_WithoutCommon = []
    Commongenes_1 = []
    for GeneName in MyListOfLists[0]:
        if GeneName in MyListOfLists[1]:
            Commongenes_1.append(GeneName)
        else:
            List_1_WithoutCommon.append(GeneName)
            
    List_2_WithoutCommon = []        
    Commongenes_2 = [] # this will have the same content as Commongenes_1
    for GeneName in MyListOfLists[1]:
        if GeneName in MyListOfLists[0]:
            Commongenes_2.append(GeneName)
        else:
            List_2_WithoutCommon.append(GeneName)
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = List_1_WithoutCommon + Commongenes_1 + List_2_WithoutCommon
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[kkk]
        C_norm = normClist[kkk]

        MyAverageCorrVector = []
        for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
            
            ######### COMPUTE GENE GENE CORRELATION MATRIX #########
            MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                                C_norm,
                                                                List_Gene_Names_for_Correlation,
                                                                N_cells,
                                                                ThrasoldForPval,
                                                                RemoveNotSignificantCorr,
                                                                WhiteDiagonal,
                                                                SampleIndex)
            
            ######### PLOT GENE GENE CORRELATION MATRIX #########
            figX, ax = plt.subplots(facecolor="white")
            
            if DoubleFace == False:
                
                if WhiteDiagonal == False:
                    MyCorrelationMatrix = MyCorrelationMatrix_Full
                elif WhiteDiagonal == True:
                    MyCorrelationMatrix = MyCorrelationMatrix_Masked
                else:
                    print("Error in Masking Correlations!")
                
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = 0
                        if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
                            
            elif DoubleFace == True:
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                        if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
            
            ax.invert_yaxis()
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss)
                    ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                    
            #Spacing between each line
            intervals = 1
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
            ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
            #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            cbar = figX.colorbar(heatmap)
            cbar.set_label('Pearson Correlation', rotation=90)
            if SampleIndex == 0:
                SampleName = "G2019S"
            elif SampleIndex == 1:
                SampleName = "WT"
            # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
            
            plt.show()

        kkk = kkk + 1
    print("Done.")
    
    
if DoGeneGeneCorrelation_SenescenceGenes_ASTRO_ONLY == 1:
    print("DoGeneGeneCorrelation_SenescenceGenes_ASTRO_ONLY")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    MyListOfLists = ListOfListsOfGenesSenescence
    
    List_1_WithoutCommon = []
    Commongenes_1 = []
    for GeneName in MyListOfLists[0]:
        if GeneName in MyListOfLists[1]:
            Commongenes_1.append(GeneName)
        else:
            List_1_WithoutCommon.append(GeneName)
            
    List_2_WithoutCommon = []        
    Commongenes_2 = [] # this will have the same content as Commongenes_1
    for GeneName in MyListOfLists[1]:
        if GeneName in MyListOfLists[0]:
            Commongenes_2.append(GeneName)
        else:
            List_2_WithoutCommon.append(GeneName)
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = List_1_WithoutCommon + Commongenes_1 + List_2_WithoutCommon
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 1
    DAY = 70      
    print('Day ' + str(DAY))
    
    C = Clist[kkk]
    C_norm = normClist[kkk]
        
    if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
        IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
        IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
    elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
        IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
        IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
        IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
        IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
        IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
        IndexesAstro_A = IndexesAstro_A.astype(int)
        IndexesAstro_B = IndexesAstro_B.astype(int)
    else:
        print("Error in selecting which method for separating astrocytes")
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
                                
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        if SampleIndex == 0:
            IndexesOfSelectedCells = IndexesAstro_A
        elif  SampleIndex == 1:
            IndexesOfSelectedCells = IndexesAstro_B
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOnSelectedCellsOfCompareObject(C,
                                                            C_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal,
                                                            IndexesOfSelectedCells)
        
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        figX, ax = plt.subplots(facecolor="white")
        
        if DoubleFace == False:
            
            if WhiteDiagonal == False:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == True:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
                    if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
                        
        elif DoubleFace == True:
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                    if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
        
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss)
                ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
        ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
        #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = figX.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)
        if SampleIndex == 0:
            SampleName = "G2019S"
        elif SampleIndex == 1:
            SampleName = "GC"
        plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
        #ax.set_title("Gene-gene Pearson Correlation Coefficient\n" + SampleName + ' at Day ' + str(DAY))
        
        plt.show()
    print("Done.")
    
    
if DoGeneGeneCorrelation_Ohashi_Genes == 1:
    print("DoGeneGeneCorrelation_Ohashi_Genes")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = L_Senescence_Ohashi_2018
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 0
    for DAY in ListOfDays:        
        print('Day ' + str(DAY))
        
        C = Clist[kkk]
        C_norm = normClist[kkk]

        MyAverageCorrVector = []
        for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
            
            ######### COMPUTE GENE GENE CORRELATION MATRIX #########
            MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                                C_norm,
                                                                List_Gene_Names_for_Correlation,
                                                                N_cells,
                                                                ThrasoldForPval,
                                                                RemoveNotSignificantCorr,
                                                                WhiteDiagonal,
                                                                SampleIndex)
            
            ######### PLOT GENE GENE CORRELATION MATRIX #########
            figX, ax = plt.subplots(facecolor="white")
            
            if DoubleFace == False:
                
                if WhiteDiagonal == False:
                    MyCorrelationMatrix = MyCorrelationMatrix_Full
                elif WhiteDiagonal == True:
                    MyCorrelationMatrix = MyCorrelationMatrix_Masked
                else:
                    print("Error in Masking Correlations!")
                
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = 0
                        if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
                            
            elif DoubleFace == True:
                MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
                for i in range(len(List_Gene_Names_for_Correlation)):
                    for j in range(len(List_Gene_Names_for_Correlation)):
                        if j>i:
                            MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                        if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                            print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                            MyCorrelationMatrix_Triangle[i,j] = 0
            
            heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
            
            ax.invert_yaxis()
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3s)
                    ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3s, rotation=90)
                    
            #Spacing between each line
            intervals = 1
            loc = plticker.MultipleLocator(base=intervals)
            ax.xaxis.set_major_locator(loc)
            ax.yaxis.set_major_locator(loc)
            ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
            ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
            #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
            ax.tick_params(labelbottom='off')    
            ax.tick_params(labelleft='off')    
            cbar = figX.colorbar(heatmap)
            cbar.set_label('Pearson Correlation', rotation=90)
            if SampleIndex == 0:
                SampleName = "G2019S"
            elif SampleIndex == 1:
                SampleName = "WT"
            # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
            # ax.set_title("Gene-gene Pearson Correlation Coefficient\n" + SampleName + ' at Day ' + str(DAY))
            
            plt.show()

        kkk = kkk + 1
    print("Done.")
    
    
if DoGeneGeneCorrelation_Ohashi_Genes_ASTRO_ONLY == 1:
    print("DoGeneGeneCorrelation_Ohashi_Genes")
    
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    # List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    List_Gene_Names_for_Correlation_BeofreCleaning = L_Senescence_Ohashi_2018
        
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = MU, 1 = CT

    ##### DETERMINE WHICH GENES ARE IN COMMON BETWEEN MU70, MU35, WT70, WT35.
    ijk = 0
    ListOfListsOfGenesForCorrelationPresentInData = []
    for DAY in ListOfDays:        
        C = Clist[ijk]
        ListGenesC = C.comm_genes() 
        List_Gene_Names_for_Correlation_PresentInData = []
        for GeneInMyList in List_Gene_Names_for_Correlation_BeofreCleaning:
            if GeneInMyList in ListGenesC:
                List_Gene_Names_for_Correlation_PresentInData.append(GeneInMyList)
        ListOfListsOfGenesForCorrelationPresentInData.append(List_Gene_Names_for_Correlation_PresentInData)
        ijk = ijk + 1
        
    List_MuWt_Day35 = ListOfListsOfGenesForCorrelationPresentInData[0]
    List_MuWt_Day70 = ListOfListsOfGenesForCorrelationPresentInData[1]

    Commongenes_Mu_Wt_Day35_Day70 = []
    for GeneName in List_MuWt_Day35:
        if GeneName in List_MuWt_Day70:
            Commongenes_Mu_Wt_Day35_Day70.append(GeneName)

    List_Gene_Names_for_Correlation = Commongenes_Mu_Wt_Day35_Day70

    Ngenes = len(Commongenes_Mu_Wt_Day35_Day70)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples * len(ListOfDays)
    ThrasoldForPval = 0.05/N_Bonf_Corr

    kkk = 1
    DAY = 70
    print('Day ' + str(DAY))
    
    C = Clist[kkk]
    C_norm = normClist[kkk]

    if SeparateAstrocytesUsingCumulativesOrClustering == "Cumulatives":
        IndexesAstro_A = ListOfIndexesOfCurrentCellTypeOfInterest[2*j]
        IndexesAstro_B = ListOfIndexesOfCurrentCellTypeOfInterest[2*j+1]
    elif SeparateAstrocytesUsingCumulativesOrClustering == "Clustering":
        IndexesAstroAllConditTimes = np.where(kmeans.labels_ == 2)[0]
        IndexesAstro_D70_Mu = [value for value in IndexesAstroAllConditTimes if value >= 2*N_cells and value < 3*N_cells]
        IndexesAstro_D70_Wt = [value for value in IndexesAstroAllConditTimes if value >= 3*N_cells and value < 4*N_cells]
        IndexesAstro_A = IndexesAstro_D70_Mu - 2*N_cells*np.ones(len(IndexesAstro_D70_Mu))
        IndexesAstro_B = IndexesAstro_D70_Wt - 2*N_cells*np.ones(len(IndexesAstro_D70_Wt))
        IndexesAstro_A = IndexesAstro_A.astype(int)
        IndexesAstro_B = IndexesAstro_B.astype(int)
    else:
        print("Error in selecting which method for separating astrocytes")
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
                                
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        if SampleIndex == 0:
            IndexesOfSelectedCells = IndexesAstro_A
        elif  SampleIndex == 1:
            IndexesOfSelectedCells = IndexesAstro_B
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOnSelectedCellsOfCompareObject(C,
                                                            C_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal,
                                                            IndexesOfSelectedCells)
        
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        figX, ax = plt.subplots(facecolor="white")
        
        if DoubleFace == False:
            
            if WhiteDiagonal == False:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == True:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
                    if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
                        
        elif DoubleFace == True:
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                    if RemoveNANs and np.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)# RdGy_r  # PiYG  # seismic   # PRGn_r
        
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 3, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3s)
                ax.text(0.1 + j, - 3, List_Gene_Names_for_Correlation[j], fontdict=font3s, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.set_xlim(0,len(List_Gene_Names_for_Correlation))
        ax.set_ylim(len(List_Gene_Names_for_Correlation),0)
        #ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = figX.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)
        if SampleIndex == 0:
            SampleName = "G2019S"
        elif SampleIndex == 1:
            SampleName = "GC"
        # plt.suptitle("Gene-gene Pearson Correlation Coefficient, " + SampleName + ' at Day ' + str(DAY))
        # ax.set_title("Gene-gene Pearson Correlation Coefficient\n" + SampleName + ' at Day ' + str(DAY))
        
        plt.show()

    kkk = kkk + 1
    print("Done.")
    
    
if PlotGeneExpressionCOL1A1vsCOL14A1toSeeCorrelation == 1:

    C_norm = C_70_norm

    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    Gene_1_Expression_B = C_norm[C.c_names.index('COL1A1'),0+N_cells:N_cells+N_cells]
    Gene_2_Expression_B = C_norm[C.c_names.index('COL14A1'),0+N_cells:N_cells+N_cells]
    plt.scatter(Gene_1_Expression_B, Gene_2_Expression_B,c='BLUE', edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("COL1A1 Expression")
    ax.set_ylabel("COL14A1 Expression")
    ax.set_title("Expression of COL1A1 vs COL14A1 to visualize correlation, Day 70, GC")
    ax.set_xlim(-1,60)
    ax.set_ylim(-1,60)
    ax.legend(["GC, Day 70"],loc='upper right')    
    plt.show()
    
    
if PlotGeneExpressionCOL1A1vsCOL3A1toSeeCorrelation == 1:

    C_norm = C_70_norm

    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)
    Gene_1_Expression_B = C_norm[C.c_names.index('COL1A1'),0:N_cells]
    Gene_2_Expression_B = C_norm[C.c_names.index('COL3A1'),0:N_cells]
    plt.scatter(Gene_1_Expression_B, Gene_2_Expression_B,c='BLUE', edgecolors='none') # 0 here and above stantds for NR2F1
    ax.set_xlabel("COL1A1 Expression")
    ax.set_ylabel("COL3A1 Expression")
    ax.set_title("Expression of COL1A1 vs COL3A1 to visualize correlation, Day 35, GC")
    #ax.set_xlim(-1,60)
    #ax.set_ylim(-1,60)
    ax.legend(["GC, Day 70"],loc='upper right')    
    plt.show()
    
    
if ExpressionNFIA == 1:
    print("ExpressionNFIA")

    CurrentList = ["NFIA","NFIA"]
    NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
        
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfDistributions_A = []
    MatrixOfDistributions_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

        MatrixOfDistributions_A.append(ListOfDistributions_A)
        MatrixOfDistributions_B.append(ListOfDistributions_B)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j=j+1

    ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
    Nlines = len(CurrentList) 
    Ncols = len(ListOfDays)
    ylimVector = [N_cells, 40]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A, MatrixOfDistributions_B, 
                              CurrentList, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if ExpressionTGFBI == 1:
    print("ExpressionNFIA")

    CurrentList = ["TGFBI","TGFBI"]
    NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
        
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfDistributions_A = []
    MatrixOfDistributions_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

        MatrixOfDistributions_A.append(ListOfDistributions_A)
        MatrixOfDistributions_B.append(ListOfDistributions_B)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j=j+1

    ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
    Nlines = len(CurrentList) 
    Ncols = len(ListOfDays)
    ylimVector = [N_cells+5, 20]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A, MatrixOfDistributions_B, 
                              CurrentList, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if Expression_NFIA_TGFBI == 1:
    print("ExpressionNFIA")

    CurrentList = ["NFIA","TGFBI"]
    NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
        
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfDistributions_A = []
    MatrixOfDistributions_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

        MatrixOfDistributions_A.append(ListOfDistributions_A)
        MatrixOfDistributions_B.append(ListOfDistributions_B)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j=j+1

    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    LabelUpBy = [0,0,0,0,0]
    LabelDownBy = 0.3*np.ones(5)
    Ncols = len(ListOfDays)
    PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                              MatrixOfFoldChanges, 
                                              MatrixOfSEMforFoldChanges, 
                                              MatrixOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel_Fold_Changes,
                                              LabelUpBy,
                                              LabelDownBy,
                                              Ncols,
                                              ['Day 35','Day 70'])

    ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
    Nlines = len(CurrentList) 
    Ncols = len(ListOfDays)
    if Expression_NFIA_TGFBI_with_ZOOM == 'Yes':
        ylimVector = [35,35]
    else:
        ylimVector = [N_cells+5, N_cells+5]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A, MatrixOfDistributions_B, 
                              CurrentList, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if ExpressionERN1 == 1:
    print("ExpressionERN1")

    CurrentList = ["ERN1","ERN1"]
    NBonferroni = len(CurrentList) * len(ListOfDays) # because 2 days
        
    j = 0
    MatrixOfFoldChanges = []
    MatrixOfSEMforFoldChanges = []
    MatrixOfPvalues = []
    MatrixOfDistributions_A = []
    MatrixOfDistributions_B = []
    for DAY in ListOfDays:
        print('Day ' + str(DAY))
        
        C = Clist[j]
        C_norm = normClist[j]
            
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        ListOfDistributions_A = []
        ListOfDistributions_B = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]
    
            GenesExpressionA = C_norm[C.c_names.index(GeneName),0:N_cells]
            GenesExpressionB = C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
            
            ListOfDistributions_A.append(GenesExpressionA)
            ListOfDistributions_B.append(GenesExpressionB)
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpressionA,GenesExpressionB)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpressionA, GenesExpressionB)
    
            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

        MatrixOfDistributions_A.append(ListOfDistributions_A)
        MatrixOfDistributions_B.append(ListOfDistributions_B)
        MatrixOfFoldChanges.append(ListOfFoldChangesThisDay)
        MatrixOfSEMforFoldChanges.append(ListOfSTDEVforFoldChangesThisDay)
        MatrixOfPvalues.append(ListOfPvalues)
        j=j+1

    ###### MAKE HISTOGRAMS OF CELLS DISTRIBUTIONS ACROSS GENE EXPRESSION ######           
    Nlines = len(CurrentList) 
    Ncols = len(ListOfDays)
    ylimVector = [N_cells, 40]
    PlotHistogramsCumulatives(Nlines, Ncols, MatrixOfDistributions_A, MatrixOfDistributions_B, 
                              CurrentList, MatrixOfPvalues, NBonferroni, 
                              ['LRRK2G2019S','LRRK2WT'], 'Gene Expression (norm. to max)',
                              ylimVector, font3s, Colors=['Red','Gray'])
    print("Done.")
    
    
if PlotSomeExpressionsForMasterEqProj == 1:
    
#    for ChosenGeneName in ['SOX2', 'COL3A1', 'COL1A1', 'NFIA', 'TH']:
#    
#        fig = plt.figure()
#        IndexChosenGeneD35 = np.where(np.array(B_35_genes) == ChosenGeneName)[0]
#        IndexChosenGeneD70 = np.where(np.array(B_70_genes) == ChosenGeneName)[0]
#        DataD35 = B_35[IndexChosenGeneD35[0],:]
#        DataD70 = B_70[IndexChosenGeneD70[0],:]
#        MaxXLim = max(max(DataD35),max(DataD70))
#        plt.hist(DataD35,bins=np.arange(-0.5,MaxXLim+0.5,1), histtype='stepfilled', color='g', normed=0, alpha=0.5, label = r'$t_1$') # 
#        plt.hist(DataD70,bins=np.arange(-0.5,MaxXLim+0.5,1), histtype='stepfilled', color='r', normed=0, alpha=0.5, label = r'$t_2$')
#        plt.xlabel('Gene Expression of ' + str(ChosenGeneName))
#        plt.ylabel('Cells')
#        plt.title('Distributions of cells across Gene Expression of ' + str(ChosenGeneName))
#        plt.legend()
#        plt.savefig('./DistributionsGeneExpression' + str(ChosenGeneName) + '.pdf')
#        
#    ChosenGeneName = 'NFIA'
#        
#    fig = plt.figure()
#    IndexChosenGeneD70 = np.where(np.array(B_70_genes) == ChsenGeneName)[0]
#    DataD70 = B_70[IndexChosenGeneD70[0],:]
#    MaxXLim = max(DataD70)
#    plt.hist(DataD70,bins=np.arange(-0.5,MaxXLim+0.5,1), histtype='stepfilled', color='b', normed=0, alpha=0.5, label = r'$t_2$')
#    plt.xlabel('Gene Expression of ' + str(ChosenGeneName))
#    plt.ylabel('Cells')
#    plt.title('Distributions of cells across Gene Expression of ' + str(ChosenGeneName))
#    plt.savefig('./DistributionsGeneExpression' + str(ChosenGeneName) + '_only_D70.pdf')
    
    ChosenGeneName1 = 'TH' # 'COL1A1' # 'SOX2' #
    ChosenGeneName2 = 'NFIA' # 'COL3A1' # 'TH' # 

    Index_Gene1_D35 = np.where(np.array(B_35_genes) == ChosenGeneName1)[0]
    Index_Gene1_D70 = np.where(np.array(B_70_genes) == ChosenGeneName1)[0]
    Data_Gene1_D35 = B_35[Index_Gene1_D35[0],:]
    Data_Gene1_D70 = B_70[Index_Gene1_D70[0],:]
    MaxXLim_Gene1 = max(max(Data_Gene1_D35),max(Data_Gene1_D70))
    
    Index_Gene2_D35 = np.where(np.array(B_35_genes) == ChosenGeneName2)[0]
    Index_Gene2_D70 = np.where(np.array(B_70_genes) == ChosenGeneName2)[0]
    Data_Gene2_D35 = B_35[Index_Gene2_D35[0],:]
    Data_Gene2_D70 = B_70[Index_Gene2_D70[0],:]
    MaxXLim_Gene2 = max(max(Data_Gene2_D35),max(Data_Gene2_D70))
    
    fig = plt.figure()
    plt.hist2d(Data_Gene1_D35, Data_Gene2_D35, cmap=plt.cm.jet) # bins=(50, 50),
    plt.show()
    
    fig = plt.figure()
    plt.hist2d(Data_Gene1_D70, Data_Gene2_D70, cmap=plt.cm.jet) #  bins=(50, 50),
    plt.show()
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm
    
    for DAY in ['D35','D70']:
    
        if DAY == 'D35':
            x, y = Data_Gene1_D35, Data_Gene2_D35
            TimeLabel = r'$t_1$'
        elif DAY == 'D70':
            x, y = Data_Gene1_D70, Data_Gene2_D70
            TimeLabel = r'$t_2$'

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # hist, xedges, yedges = np.histogram2d(x, y, bins=4, range=[[0, 4], [0, 4]])
        hist, xedges, yedges = np.histogram2d(x, y, bins=4, range=[[0, 4], [0, 4]])
        
        # Construct arrays for the anchor positions of the 16 bars.
        # Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
        # ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
        # with indexing='ij'.
        xpos, ypos = np.meshgrid(xedges[:-1]-0.5, yedges[:-1]-0.5)
        xpos, ypos = np.meshgrid(xedges[:-1]-0.5, yedges[:-1]-0.5)
        xpos = xpos.flatten('F')
        ypos = ypos.flatten('F')
        zpos = np.zeros_like(xpos)
        
        # Construct arrays with the dimensions for the 16 bars.
        # dx = 0.5 * np.ones_like(zpos)
        dx = np.ones_like(zpos)
        dy = dx.copy()
        dz = hist.flatten()
        
        cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
        max_height = np.max(dz)   # get range of colorbars so we can normalize
        min_height = np.min(dz)
        # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in dz] 
    
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
        
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
        plt.xlabel('Expression of Gene ' + str(ChosenGeneName1))
        plt.ylabel('Expression of Gene ' + str(ChosenGeneName2))
        plt.title('Distributions of cells across Gene Expressions at ' + TimeLabel)
        plt.legend()
        plt.savefig('./DistributionsGeneExpression_' + str(ChosenGeneName1) + '_' + str(ChosenGeneName2) +'_' + DAY + '.pdf')
    
        plt.show()