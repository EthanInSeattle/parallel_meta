# Yuxi Cui
# Shengnan Sun
# how can we understand mircobiome
# parallel-meta duplicate

# internal dependency
import os
import csv
import math
import sys


# statistics dependency
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
import skbio
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes
from numpy import std, mean, sqrt


path = sys.argv[1]

# takes the argument from the user specifying sample file directory
pathDirTemp = os.listdir(path + '/samples/')
pathDir = []

for d in pathDirTemp:
    if d[0] == '.':
        pathDirTemp.remove(d)

for d in pathDirTemp:
    if 'B' in d:
        pathDir.append(d)

for d in pathDirTemp:
    if "I" in d:
        pathDir.append(d)


#############################################
# compute the alpha diversity of the samples#
#############################################

# pre: sub is a substring of s, n <= max appearance of sub in str
# post: return the index of the nth appearance of sub in s
def nthIndexOf(sub, s, n):
    index = s.find(sub)
    if index == -1:
        return -1
    if n > 0:
        return index + 1 + nthIndexOf(sub, s[s.find(sub) + 1:], n - 1)
    else:
        return 0

# pre: takes a list containing the abundance number of a sample
# post: return the Shannon-Wienner index of that sample
def computeShannon(abd):
    res = 0.0
    for data in abd:
        if data > 0:
            res += data * math.log(data)
    return -1.0 * res

# count the total number of species and otus in the sample
allOTUs = ['']
allGenus = ['']
allAlpha = ['', "Observed OTU", "Shannon Index"]

for i in range(0, len(pathDir)):
    try:
        f = open(path + '/samples/' + pathDir[i] + "/classification.txt", "r")
    except IOError:
        print 'Please enter a correct file/directory path'

    f = open(path + '/samples/' + pathDir[i] + "/classification.txt", "r")
    line1 = f.readline()
    data = f.readlines()
    for line in data:
        genus = line[nthIndexOf("; ", line, 5) + 1: nthIndexOf(";", line, 6) - 1]
        family = line[nthIndexOf("; ", line, 4) + 1: nthIndexOf(";", line, 5) - 1]
        name = family[0:3] + "_" + genus
        if name not in allGenus:
            allGenus.append(name)

        otu = line[nthIndexOf("; ", line, 7) + 1:]
        if otu not in allOTUs:
            allOTUs.append(otu)

genusCountTable = []
genusAbdTable = []
genusAlphaTable = []

otuCountTable = []
otuAbdTable = []
otuAlphaTable = []

for i in range(0, len(pathDir)):
    f = open(path + '/samples/' + pathDir[i] + "/classification.txt", "r")
    line1 = f.readline()
    data = f.readlines()
    total = 0

    genusDic = {}
    genusCountTemp = [pathDir[i]]
    genusAbdTemp = [pathDir[i]]
    genusAlphaTemp = [pathDir[i]]

    otuDic = {}
    otuCountTemp = [pathDir[i]]
    otuAbdTemp = [pathDir[i]]
    otuAlphaTemp = [pathDir[i]]

    # count
    for line in data:
        total += 1.0
        genus = line[nthIndexOf("; ", line, 5) + 1: nthIndexOf(";", line, 6) - 1]
        family = line[nthIndexOf("; ", line, 4) + 1: nthIndexOf(";", line, 5) - 1]
        name = family[0:3] + "_" + genus
        otu = line[nthIndexOf("; ", line, 7) + 1:]

        # make genus dictionary
        if genusDic.has_key(name):
            genusDic[name] += 1
        else:
            genusDic[name] = 1

        # make otu dictionary
        if otuDic.has_key(otu):
            otuDic[otu] += 1
        else:
            otuDic[otu] = 1

    # make genus count table and abundance table
    for g in allGenus:
        if genusDic.has_key(g):
            genusCountTemp.append(genusDic[g])
            genusAbdTemp.append(genusDic[g] / total)

        elif g != '':
            genusCountTemp.append(0)
            genusAbdTemp.append(0.0)
    genusCountTable.append(genusCountTemp)
    genusAbdTable.append(genusAbdTemp)

    genusAlphaTemp.append(len(genusDic.keys()))
    genusAlphaTemp.append(computeShannon(genusAbdTemp[1:]))
    genusAlphaTable.append(genusAlphaTemp)

    # make otu count table and abundance table
    for o in allOTUs:
        if otuDic.has_key(o):
            otuCountTemp.append(otuDic[o])
            otuAbdTemp.append(otuDic[o] / total)
        elif o != '':
            otuCountTemp.append(0)
            otuAbdTemp.append(0.0)
    otuCountTable.append(otuCountTemp)
    otuAbdTable.append(otuAbdTemp)

    otuAlphaTemp.append(len(otuDic.keys()))
    otuAlphaTemp.append(computeShannon(otuAbdTemp[1:]))
    otuAlphaTable.append(otuAlphaTemp)


fs = open(path + '/meta.txt')
line1 = fs.readline()
data = fs.readlines()
for line in data:
    c = line.split('\t')
    for i in range(0, len(otuAlphaTable)):
        if (c[0] == otuAlphaTable[i][0]) :
            otuAlphaTable[i].append(c[1])
            otuAlphaTable[i].append(c[2])
            otuAlphaTable[i].append(c[3])
        if (c[0] == genusAlphaTable[i][0]) :
            genusAlphaTable[i].append(c[1])
            genusAlphaTable[i].append(c[2])
            genusAlphaTable[i].append(c[3])
fs.close()

# Write genus count table
csvfile = file(path + '/genusCount.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allGenus)
writer.writerows(genusCountTable)
csvfile.close()

# Write genus abundance table
csvfile = file(path + '/genusAbd.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allGenus)
writer.writerows(genusAbdTable)
csvfile.close()

# Write OTU count table
csvfile = file(path + '/otuCount.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allOTUs)
writer.writerows(otuCountTable)
csvfile.close()

# Write OTU abundance table
csvfile = file(path + '/otuAbd.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allOTUs)
writer.writerows(otuAbdTable)
csvfile.close()

# Write genus alpha table
csvfile = file(path + '/genusAlpha.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allAlpha + ['Status', 'Sex', 'Smoking'])
writer.writerows(genusAlphaTable)
csvfile.close()

# Write otu alpha table
csvfile = file(path + '/otuAlpha.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerow(allAlpha + ['Status', 'Sex', 'Smoking'])
writer.writerows(otuAlphaTable)
csvfile.close()


############################################
# compute the beta diversity of the samples#
############################################

class TreeNode (object):

    # constructor of the class
    def __init__(self, name, distance):
        self.name = name
        self.distance = distance
        self.left = None
        self.right = None
        self.ancestor = None
        self.brother = None
        self.p1 = 0.0
        self.p2 = 0.0
        #self.value = name + ':' + str(distance)

    # post: return true if the current node has no child, false otherwise
    def isLeafNode(self):
        return (self.left == None) & (self.right == None)

    # post: return the brother node of the current node
    def brother(self):
        if (self.ancestor.left) == self:
            return self.ancestor.right
        else:
            return self.ancestor.left

def getAbd(name, sample):
    for i in range(0, len(otuAbdTable)):
        if otuAbdTable[i][0] == sample:
            for j in range(0, len(otuAbdTable[0])):
                if otuAbdTable[0][j] == name:
                    return otuAbdTable[i][j]
    return 0.0


# pre: take a string that is in strict newick tree format where only the leaf node and all distance is shown
# post: convert the tree in newick format into binary tree format and return the root TreeNode
def newickToTree(st):
     if (',' not in st): #leaf node
        #print st
        data = st.split(':')
        #print data
        return TreeNode(data[0], float(data[1]))
     else:
         #print 'HEY ' + st
         dist = float(st[st.rfind(':') + 1:])
         mod = st [1 : st.rfind(')')]
         #print 'mod ' + mod
         sepIndex = sepCommaIndex(mod)
         left = mod[0 : sepIndex]
         right = mod[sepIndex + 1 : ]
         res = TreeNode('FAKE', dist)
         res.left = newickToTree(left)
         res.right = newickToTree(right)
         res.left.ancestor = res
         res.right.ancestor = res
         res.left.brother = res.right
         res.right.brother = res.left
         return res

# pre: helper method for newickToTree, takes a string in newick format and
#      return the index of the comma that separates left and right children
# post: return the index of the comma that separates the 2 children of the current node
def sepCommaIndex(st):
    lf = st[0 : st.rfind(',')]
    rt = st[st.rfind(',') + 1:]
    if ')' not in rt:
        return st.rfind(',')
    else:
        rtCount = rt.count(')')
        lfCount = 0
        index = len(lf) - 1
        for c in lf[::-1]:
            if lfCount == rtCount:
                return index
            else:
                index -= 1
                if c == '(':
                    lfCount += 1
                elif c == ')':
                    rtCount += 1
        return index

# pre: takes a TreeNode containing information about the species name and evolutionary distance of a phlyogenetic tree
# post: return a string that represents the tree in Newick format
def tree_to_newick(node):
    if node != None:
        if (node.isLeafNode()):
            return node.name + ':' + str(node.distance)
        else:
           return '(' + tree_to_newick(node.left) + ',' + tree_to_newick(node.right) + ')' + ':' + str(node.distance)
    else:
        return ''

# pre: takes a tree node curr, a sample name s1, and another sample name s2
# post: recalculate the abundance of
def reduce(curr, s1, s2):
    if curr.isLeafNode():
        curr.p1 = getAbd(curr.name, s1)
        curr.p2 = getAbd(curr.name, s2)
    m = min(curr.p1, curr.p2)
    if curr.ancestor != None:
        curr.ancestor.p1 += (1 - curr.distance) * (curr.p1 - m)
        curr.ancestor.p2 += (1 - curr.distance) * (curr.p2 - m)


# pre: takes a tree node curr, a sample name s1, and another sample name s2
# post: calculate the similarity score of the two sample that can be used to calculate their beta diversity
def getSimilarity(curr, s1, s2):
    if curr.isLeafNode():
        reduce(curr, s1, s2)
        return min(curr.p1, curr.p2)
    else:
        childScore = getSimilarity(curr.left, s1, s2) + getSimilarity(curr.right, s1, s2)
        reduce(curr, s1, s2)
        return childScore + min(curr.p1, curr.p2)


def zero(curr):
    curr.p1 = 0.0
    curr.p2 = 0.0
    if not curr.isLeafNode():
        zero(curr.left)
        zero(curr.right)


datafile = open(path + '/database.txt')
ref = datafile.read()
datafile.close()

msDistMatrix = [['']]
refTree = newickToTree(ref)

for d in pathDir:
    msDistMatrix[0].append(d)
    msDistMatrix.append([d])


name = ['']
for i in range(1, len(allOTUs)):
    name.append(allOTUs[i].split('otu_')[1].split('\n')[0])
otuAbdTable.insert(0, name)

def getEuDist(s1, s2):
    dist = 0.0
    ind1 = 0
    ind2 = 0
    for i in range(0, len(otuAbdTable)):
        if otuAbdTable[i][0] == s1:
            ind1 += i
        if otuAbdTable[i][0] == s2:
            ind2 += i

    for j in range(1, len(otuAbdTable[ind1])):
        dist += ((otuAbdTable[ind1][j] - otuAbdTable[ind2][j]) ** 2)

    return math.sqrt(dist)

euDistMatrix = [['']]
for d in pathDir:
    euDistMatrix[0].append(d)
    euDistMatrix.append([d])


for i in range(1, len(euDistMatrix[0])):
    for j in range(1, len(euDistMatrix)):
        euDistMatrix[j].append(getEuDist(euDistMatrix[0][i], euDistMatrix[j][0]))



# Write eucladian distance matrix table to calculate beta diversity
csvfile = file(path + '/euDistMatrix.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerows(euDistMatrix)
csvfile.close()


for i in range(1, len(msDistMatrix[0])):
    for j in range(1, len(msDistMatrix)):
        if i != j:
            zero(refTree)
            msDistMatrix[j].append(1.0 - getSimilarity(refTree, msDistMatrix[0][i], msDistMatrix[j][0]))
        else:
            msDistMatrix[j].append(0.0)



# Write meta storm distance matrix table to calculate beta diversity
csvfile = file('msDistMatrix.csv', 'wb')
writer = csv.writer(csvfile)
writer.writerows(msDistMatrix)
csvfile.close()




########
########
########





def Alpha_analysis(path):
    # input data frame contains columns of shannon index, following by feature of each sample. Row names are sample ID.
    # create a feature list
    
    alpha_df=pd.read_csv(path + '/genusAlpha.csv',index_col=0)
    alpha_df=alpha_df.drop(alpha_df.columns[0], axis=1)
    features_df = alpha_df
    features_df = features_df.drop(features_df.columns[0], axis=1)
    n = features_df.shape[1]
    features_list = [[] for _ in range(n)]
    for i in range(n):
        features_list[i] = (list(set(features_df.ix[:, i].values)))

    # split data into different features
    lists = [[] for _ in range(n * 2)]
    p_list = [[] for _ in range(n)]
    for i in range(len(features_list)):
        for j in range(alpha_df.shape[0]):
            if alpha_df.ix[j][i + 1] == features_list[i][0]:
                lists[i * 2].append(alpha_df.ix[j][0])
            else:
                lists[i * 2 + 1].append(alpha_df.ix[j][0])
    # find p value
    for i in range(n):
        z, p = scipy.stats.ranksums(lists[i * 2], lists[i * 2 + 1])
        p_list[i].append(p)
    p_df = pd.DataFrame(data=p_list, index=features_df.columns.values, columns=['p'])
    p_df.to_csv(path + '/p_alpha.csv', sep='\t')
    # plot alpha diversity
    for i in range(n):
        feature_1 = pd.DataFrame(data=lists[i * 2], index=None, columns=[features_list[i][0]])
        feature_2 = pd.DataFrame(data=lists[i * 2 + 1], index=None, columns=[features_list[i][1]])
        data_df = pd.concat([feature_1, feature_2], axis=1)
        figure(i + 1)
        sns.set(font_scale=1.5)
        ax1 = sns.boxplot(data_df)
        plt.xlabel(alpha_df.columns.values[i + 1], fontsize=20)
        plt.ylabel('Shannon', fontsize=20)
        plt.savefig(path + '/alpha' + alpha_df.columns.values[i + 1] + '.pdf')
        plt.close()

# In[3]:

def beta_analysis(path):
    # match features with each distance data
    distance_df=pd.read_csv(path + '/msDistMatrix.csv',index_col=0)
    feature_df=pd.read_csv(path + '/feature_df.csv',index_col=0,sep='\t')
    n = feature_df.shape[1]
    data = []
    lists_featureid = [[] for _ in range(n * 2)]
    for i in range(distance_df.shape[0] - 1):
        for j in range(i + 1, distance_df.shape[1]):
            data.append(distance_df.ix[i, j])
            for z in range(n):
                lists_featureid[z * 2].append(feature_df.ix[i][z])
                lists_featureid[z * 2 + 1].append(feature_df.ix[j][z])

    # create a feature list
    feature_list = [[] for _ in range(n)]
    for i in range(n):
        feature_list[i] = (list(set(feature_df.ix[:, i].values)))

    # split data for each feature in within or between
    data_featuregroup = [[] for _ in range(n * 2)]
    for i in range(n):
        for j in range(len(data)):
            if lists_featureid[i * 2][j] == lists_featureid[i * 2 + 1][j]:
                data_featuregroup[i * 2].append(data[j])
            else:
                data_featuregroup[i * 2 + 1].append(data[j])

    # effect size
    def cohen_d(x, y):
        nx = len(x)
        ny = len(y)
        dof = nx + ny - 2
        return (mean(x) - mean(y)) / sqrt(((nx - 1) * std(x, ddof=1) ** 2 + (ny - 1) * std(y, ddof=1) ** 2) / dof)

    # plot data
    # function for setting the colors of the box plots pairs
    def setBoxColors(bp):
        setp(bp['boxes'][0], color='blue')
        setp(bp['caps'][0], color='blue')
        setp(bp['caps'][1], color='blue')
        setp(bp['whiskers'][0], color='blue')
        setp(bp['whiskers'][1], color='blue')
        setp(bp['fliers'][0], markeredgecolor='blue')
        setp(bp['medians'][0], color='blue')

        setp(bp['boxes'][1], color='red')
        setp(bp['caps'][2], color='red')
        setp(bp['caps'][3], color='red')
        setp(bp['whiskers'][2], color='red')
        setp(bp['whiskers'][3], color='red')
        setp(bp['fliers'][1], markeredgecolor='red')
        setp(bp['medians'][1], color='red')

    fig = figure()
    ax = axes()
    hold(True)
    sns.set(font_scale=1.5)
    d = []
    for i in range(n):
        data_plot = [data_featuregroup[i * 2], data_featuregroup[i * 2 + 1]]
        bp = boxplot(data_plot, positions=[1 + 3 * i, 2 + 3 * i], widths=0.6)
        setBoxColors(bp)
        effect_size = cohen_d(data_featuregroup[i * 2], data_featuregroup[i * 2 + 1])
        d.append(effect_size)

    def floatrange(start, stop, steps):
        return [start + float(i) * (stop - start) / (float(steps) - 1) for i in range(steps)]

    # set axes limits and labels
    xlim(0, n * 3)
    ylim(0, max(data) + 0.05)
    ax.set_xticklabels(feature_df.columns.values)
    ax.set_xticks(floatrange(1.5, n * 3 - 1.5, 3))
    ax.set_ylabel('Distance')

    # draw temporary red and blue lines and use them to create a legend
    hB, = plot([1, 1], 'b-')
    hR, = plot([1, 1], 'r-')
    legend((hB, hR), ('Within', 'Between'))
    hB.set_visible(False)
    hR.set_visible(False)

    # beta diversity p value from permanova
    # two inputs, distance data and grouping
    distance_data = distance_df.values
    distance_ids = distance_df.index.values
    distance_matrix = skbio.stats.distance.DistanceMatrix(distance_data, distance_ids)
    p_list = []
    for i in range(n):
        grouping = feature_df.ix[:, i].values
        p_beta = skbio.stats.distance.permanova(distance_matrix, grouping, permutations=999)['p-value']
        p_list.append(p_beta)

    # save file
    plt.savefig(path + '/beta_pt.pdf')
    plt.close()
    p_df = pd.DataFrame(data=p_list, index=feature_df.columns.values, columns=['p'])
    e_df = pd.DataFrame(data=d, index=feature_df.columns.values, columns=['Cohens_d'])
    result_df = pd.concat([p_df, e_df], axis=1)
    result_df.to_csv(path+'/beta_data.csv', sep='\t')


# In[4]:

def biomarker(path):
    # input is the abundance table. The last column contain the feature infor
    abundance_df=pd.read_csv(path + '/genusAbd.csv',index_col=0)
    feature_df=pd.read_csv(path + '/biofeature.csv',index_col=0,sep='\t')
    feature_list = list(set(feature_df.ix[:, -1].values))
    count_matrix = np.zeros((abundance_df.shape[1], 1))
    for i in range(abundance_df.shape[1]):
        feature = [[] for _ in range(2)]
        for j in range(abundance_df.shape[0]):
            if feature_df.ix[j][-1] == feature_list[0]:
                feature[0].append(abundance_df.ix[j][i])
            else:
                feature[1].append(abundance_df.ix[j][i])
            z_stat, p_val = scipy.stats.ranksums(feature[0], feature[1])
            count_matrix[i] = p_val
    index = abundance_df.columns.values
    sort_df = pd.DataFrame(data=count_matrix, index=index, columns=['P_val'])
    sort_df = sort_df[sort_df['P_val'] <= 0.01].sort_values('P_val')

    # find genus with higher relative abundance
    threshold_id = abundance_df[abundance_df > 0.01]
    threshold_id = threshold_id.dropna(axis=1, thresh=abundance_df.shape[0] / 4).columns.values

    # plot marker with high relative abundance
    p_index = []
    p_list = []
    for i in range(len(threshold_id)):
        f1 = []
        f2 = []
        if threshold_id[i] in sort_df.index.values:
            ab_data = abundance_df[threshold_id[i]]
            p_index.append(threshold_id[i])
            p_list.append(sort_df.ix[threshold_id[i]])
            for j in range(len(ab_data)):
                if feature_df.ix[:, -1][j] == feature_list[0]:
                    f1.append(ab_data[j])
                else:
                    f2.append(ab_data[j])
            f1 = pd.DataFrame(data=f1, index=None, columns=[feature_list[0]])
            f2 = pd.DataFrame(data=f2, index=None, columns=[feature_list[1]])
            data_df = pd.concat([f1, f2], axis=1)
            figure(i + 1)
            sns.set(font_scale=1.5)
            ax1 = sns.boxplot(data_df, vert=False)
            plt.xlabel('Relative Abundance', fontsize=20)
            plt.ylabel(feature_df.columns.values[0], fontsize=20)
            ax1.set_title(threshold_id[i], fontsize=20)
            plt.savefig(path + '/biomarker' + threshold_id[i] + '.pdf')
            plt.close()
    p_df = pd.DataFrame(data=p_list, index=p_index, columns=['P_val'])
    p_df.to_csv(path + '/biomarker.csv', sep='\t')


alpha_analysis=Alpha_analysis(path)
beta_analysis=beta_analysis(path)
biomarker=biomarker(path)