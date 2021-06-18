# /usr/bin/python

from collections import defaultdict
import math
import sys
import re
from optparse import OptionParser
import numpy as np
import scipy.stats
import subprocess
import shlex
from bisect import bisect_right
import itertools
from time import time
import random
import os
from operator import mul

"""
find the gene expression correlation between regulator and target genes
get # of network edges and transciption factor for top genes predicted by each method. 
"""
parser = OptionParser()

parser.add_option("-i", "--input", type="string", dest="input", metavar="FILE",help="the GRN files for human, cols are target,TF,zscore,corralation,type,species,note: the edage may display multiple times at different types" ) 
parser.add_option("-d", "--dir", type="string", dest="dir", default="/home/tmhjxl57/scratch/EVI1/ENCODE/bowtie/MECOM_selected/",help="the directory of peak width for each active marker" ) 
parser.add_option("-e", "--exp", type="string", dest="exp", metavar="FILE",help="the exp profile for 10 cell types,including both FMPK and diff rank" ) 
parser.add_option("-n", "--num", type="int", dest="num",default=400,help="number of CIGs considered" )

(options, args) = parser.parse_args()
def loadWidth(gene2width,file):
    start_time=time()
    if file:
        f=open(file)
    filebase=os.path.basename(file)
    marker=filebase.strip().split(".")[0]
    #print marker
    genes=[]
    while True:
        line=f.readline()
        if line.startswith("gene_id"):continue
        if not line:
            break
        cols=line.strip().split("\t")
        gene=cols[0]
        width=float(cols[1])
        genes.append([gene,width])
    genes_sorted=sorted(genes,key=lambda x: x[1],reverse=True)
    rank=1
    for gene in genes_sorted:
        gene2width[gene[0]][marker]=[gene[1],rank]
        rank+=1
    #print len(genes_sorted)
    #for gene in CIGs_sorted:
    #    print "\t".join(map(str,gene))
    #genes=[i[0] for i in genes_sorted][:n]
    #CIG_all[cell]=CIGs

def combineWidth(dir,markers,N,index):
    gene2width=defaultdict(dict)
    n=0
    for file in os.listdir(dir):
        #cell=file.strip().split("_")[0]
        marker=file.strip().split(".")[0]
        if marker not in markers:continue
        file=dir+file
        #print cell,file
        loadWidth(gene2width,file)
        n+=1
        #break
    genes=[]
    for gene in gene2width:
        #print gene, gene2width[gene]
        if len(gene2width[gene])<n:continue
        ranks=[i[1] for i in gene2width[gene].values()]
        av_rank=sum(ranks)
        rank_product=reduce(mul,ranks)
        genes.append([gene,av_rank,rank_product,gene2width[gene]['H3K4me1'][1],gene2width[gene]['H3K4me2'][1],gene2width[gene]['H3K4me3'][1],gene2width[gene]['H3K9ac'][1],gene2width[gene]["H3K27ac"][1]])
        #print "\t".join(map(str,(gene, av_rank,rank_product,"\t".join(map(str,(gene2width[gene].items()))))))
    genes_sorted=sorted(genes,key=lambda x: x[index])
    for gene in genes_sorted:
        #print "\t".join(map(str,(gene[0],gene[1])))
        print "\t".join(map(str,gene))
    #return gene2width
    CIGs=[gene[0] for gene in genes_sorted][:N]
    return CIGs

def loadCSV(CIG_all,file,n):
    start_time=time()
    if file:
        f=open(file)
    filebase=os.path.basename(file)
    cell=filebase.strip().split("_")[0]
    #print cell
    CIGs=[]
    while True:
        line=f.readline()
        if line.startswith("gene_id"):continue
        if not line:
            break
        cols=line.strip().split(",")
        gene=cols[0]
        dis=float(cols[1])
        nonCIGpro=float(cols[2])
        CIGpro=float(cols[3])
        pvalue=float(cols[4])
        FDR=float(cols[5])
        CIGs.append([gene,dis,FDR])
    CIGs_sorted=sorted(CIGs,key=lambda x: x[1],reverse=True)
    #for gene in CIGs_sorted:
    #    print "\t".join(map(str,gene))
    CIGs=[i[0] for i in CIGs_sorted][:n]
    CIG_all[cell]=CIGs
    #print len(CIGs)
    #print CIGs
    #return CIGs
def ParseCIG(dir,n):
    CIG_all=defaultdict(list)
    for file in os.listdir(dir):
        #cell=file.strip().split("_")[0]
        file=dir+file
        #print cell,file
        loadCSV(CIG_all,file,n)
        #break
    #for cell in CIG_all:
        #print cell, len(CIG_all[cell]),CIG_all[cell]
    return CIG_all
def mannWhitenyU(l1,l2):
    if max(l1)==min(l1) and max(l2)==min(l2) and max(l1)==max(l2) and min(l1)==min(l2):
        return 1.0
    s,p=scipy.stats.mannwhitneyu(l1,l2,alternative='greater')
    return p


def GRN(file,TF2TG):
    start_time=time()
    if file:
        f=open(file)
    CIGs=[]
    while True:
        line=f.readline()
        if line.startswith("#"):continue
        if not line:
            break
        cols=line.strip().split("\t")
        gene=cols[0]
        CIGs.append(gene)
    outedge_count=defaultdict(int)
    inedge_count=defaultdict(int)
    outfile=file+".GRN.txt"
    of=open(outfile,"w")
    for TF in CIGs:
        for TG in CIGs:
            if TF==TG:continue
            if TF not in TF2TG:continue
            if TG in TF2TG[TF]:
                outedge_count[TF]+=1
                inedge_count[TG]+=1
                of.write("\t".join(map(str,TF2TG[TF][TG]))+"\n")
    of.close()
    #print len(outedge_count.keys()),len(inedge_count.keys())
    #union=set(outedge_count.keys()+inedge_count.keys())
    #for gene in union:
    #    print "\t".join(map(str,(gene,outedge_count[gene],inedge_count[gene],outedge_count[gene]+inedge_count[gene])))
    return CIGs
def loadGenes(file,TF2TG):
    start_time=time()
    if file:
        f=open(file)
    CIGs=[]
    while True:
        line=f.readline()
        if line.startswith("#"):continue
        if not line:
            break
        cols=line.strip().split("\t")
        gene=cols[0]
        CIGs.append(gene)
    return CIGs
def loadExp(file):
    start_time=time()
    if file:
        f=open(file)
    gene2exp=defaultdict(dict)
    #highExpGenes=[]
    #highUpGenes=[]
    #highDownGenes=[]
    while True:
        line=f.readline()
        if line.startswith(("#","test_id")):continue
        if not line:
            break
        cols=line.strip().split("\t")
        gene=cols[0]
        sample1=cols[1]
        sample2=cols[2]
        FPKM1=float(cols[3])
        FPKM2=float(cols[4])
        qvalue=float(cols[5])
        if FPKM1 <0.1 and FPKM2<0.1:continue
        if sample1!="H1-hESC" and sample2!="H1-hESC":continue
        if sample1=="H1-hESC":
            cell=sample2
            exp=FPKM2
            exp_bg=FPKM1
        elif sample2=="H1-hESC":
            cell=sample1
            exp=FPKM1
            exp_bg=FPKM2
        if exp >=exp_bg:
            flag="up"
        else:
            flag="down"
        #print cell,gene,sample1,sample2,FPKM1,FPKM2,qvalue,flag
        if cell not in gene2exp:
            gene2exp[cell]=defaultdict(list)
        gene2exp[cell][gene]=[exp,qvalue,flag]
    return gene2exp
def sortByExp(d,N):
    d_sorted=sorted(d.items(), key=lambda x: x[1][0],reverse=True)
    genes=[i[0] for i in d_sorted]
    return genes[:N]
def sortByDiff(d,flag,N):
    d_sorted=sorted(d.items(), key=lambda x: x[1][1])
    genes=[i[0] for i in d_sorted if i[1][2]==flag]
    return genes[:N]
#def GRNStat(genes,TF2TG):

def selectGeneExp(gene2exp,cell_type,N):
    highExpGenes={}
    highUpGenes={}
    highDownGenes={}
    for cell in gene2exp:
        if cell != cell_type:continue
        genes=gene2exp[cell]
        exps=sortByExp(genes,N)
        ups=sortByDiff(genes,"up",N)
        downs=sortByDiff(genes,"down",N)

    return exps,ups,downs
def edgeCount(genes,TF2TG):
    outedge_count=defaultdict(int)
    inedge_count=defaultdict(int)
    TF_genes={}
    #outfile=file+".GRN.txt"
    #of=open(outfile,"w")
    edge_count=0
    edge_count_all=0
    outedge_count_all=0
    inedge_count_all=0
    edge_out=defaultdict(int)
    edge_in=defaultdict(int)
    for TF in genes:
        if TF in TF2TG:
            TF_genes[TF]=1
            edge_count_all+=len(TF2TG[TF])
            outedge_count_all+=len(TF2TG[TF])
            edge_out[TF]=len(TF2TG[TF])
            for reg in TF2TG:
                if TF in TF2TG[reg]:
                    edge_count_all+=1
                    inedge_count_all+=1
                    edge_in[TF]+=1
        for TG in genes:
            if TF==TG:continue
            if TF not in TF2TG:continue
            if TG in TF2TG[TF]:
                outedge_count[TF]+=1
                inedge_count[TG]+=1
                edge_count+=1
    return edge_count,len(TF_genes),float(edge_count)/len(TF_genes)

def edgeCountAll(genes,TF2TG):
    TF_genes={}
    edge_count=0
    outedge_count=0
    inedge_count=0
    for TF in genes:
        if TF in TF2TG:
            TF_genes[TF]=1
            edge_count+=len(TF2TG[TF])
            outedge_count+=len(TF2TG[TF])
            for reg in TF2TG:
                if TF in TF2TG[reg]:
                    edge_count+=1
                    indege_count+=1

    print  edge_count,outedge_count, inedge_count,len(TF_genes),float(edge_count)/len(TF_count),
def effective_number(TF2TG,CIGs):
    TF_set=set(TF2TG.keys())
    TG_set=set()
    for TF in TF2TG:
        TG_set.update(TF2TG[TF].keys())
    genes=TG_set.union(TF_set)
    TF_count=0
    TG_count=0
    nomatch_count=0
    for gene in CIGs:
        if gene in TF_set:
            TF_count+=1
        elif gene in TG_set:
            TG_count+=1
        else:
            nomatch_count+=1
    return TF_count+TG_count
def simulator_GRN(TF2TG,gene_count,N):
    TF_set=set(TF2TG.keys())
    TG_set=set()
    for TF in TF2TG:
        TG_set.update(TF2TG[TF].keys())
    genes=TG_set.union(TF_set)
    edge_sim,TF_sim,edgePerTF_sim=[],[],[]
    
    for i in range(N):
        sim=random.sample(genes,gene_count)
        a,b,c=edgeCount(sim,TF2TG)
        edge_sim.append(a)
        TF_sim.append(b)
        edgePerTF_sim.append(c)
        print  >> sys.stderr, "simulation at rounds: ", i
    edge_sim_mean=np.mean(edge_sim)
    edge_sim_std=np.std(edge_sim)
    TF_sim_mean=np.mean(TF_sim)
    TF_sim_std=np.std(TF_sim)
    edgePerTF_sim_mean=np.mean(edgePerTF_sim)
    edgePerTF_sim_std=np.std(edgePerTF_sim)

    file1="edge_random.txt"
    file2="TF_random.txt"
    f1=open(file1,"w")
    f2=open(file2,"w")
    for i in edge_sim:
        f1.write(str(i)+"\n")
    for i in TF_sim:
        f2.write(str(i)+"\n")
    f1.close()
    f2.close()
    return edge_sim,TF_sim,edgePerTF_sim
def simulator_uniform(TF2TG,CIGs,N):
    TF_set=set(TF2TG.keys())
    TG_set=set()
    for TF in TF2TG:
        TG_set.update(TF2TG[TF].keys())
    genes=TG_set.union(TF_set)
    TF_count=0
    TG_count=0
    nomatch_count=0
    for gene in CIGs:
        if gene in TF_set:
            TF_count+=1
        elif gene in TG_set:
            TG_count+=1
        else:
            nomatch_count+=1
    edge_sim,TF_sim,edgePerTF_sim=[],[],[]
    for i in range(N):
        sim=random.sample(genes,TF_count+TG_count)
        a,b,c=edgeCount(sim,TF2TG)
        edge_sim.append(a)
        TF_sim.append(b)
        edgePerTF_sim.append(c)
    edge_sim_mean=np.mean(edge_sim)
    edge_sim_std=np.std(edge_sim)
    TF_sim_mean=np.mean(TF_sim)
    TF_sim_std=np.std(TF_sim)
    edgePerTF_sim_mean=np.mean(edgePerTF_sim)
    edgePerTF_sim_std=np.std(edgePerTF_sim)
    edge,TF,edgePerTF=edgeCount(CIGs,TF2TG)

    print "\t".join(map(str,("SIM",edge_sim_mean,edge_sim_std,TF_sim_mean,TF_sim_std,edgePerTF_sim_mean,edgePerTF_sim_std)))
    print "\t".join(map(str,("CIGs",edge,0,TF,0,edgePerTF,0)))
    return TF_count+TG_count
def simulator(TF2TG,CIGs,N):
    TF_set=set(TF2TG.keys())
    TG_set=set()
    for TF in TF2TG:
        TG_set.update(TF2TG[TF].keys())
    #print len(TG_set)
    TG_set.difference_update(TF_set)
    TF_count=0
    TG_count=0
    nomatch_count=0
    for gene in CIGs:
        if gene in TF_set:
            TF_count+=1
        elif gene in TG_set:
            TG_count+=1
        else:
            nomatch_count+=1
    print len(CIGs),TF_count,TG_count,nomatch_count
    #random.seed(1)
    edge,outNode,inNode=edgeCount(CIGs,TF2TG)
    edge_sim,outNote_sim,inNode_sim=[],[],[]
    for i in range(N):
        TF_sim=random.sample(TF_set,TF_count)
        TG_sim=random.sample(TG_set,TG_count)
        sim=TF_sim+TG_sim
        a,b,c=edgeCount(sim,TF2TG)
        edge_sim.append(a)
    edge_sim_mean=np.mean(edge_sim)
    edge_sim_std=np.std(edge_sim)
    print "\t".join(map(str,(edge,0)))
    print "\t".join(map(str,(edge_sim_mean,edge_sim_std)))
    #print len(TF_set),len(TG_set)
def load_GRN(file):
    start_time=time()
    if file:
        f=open(file)
    TF2TG=defaultdict(dict)
    while True:
        line=f.readline()
        if line.startswith("#"):continue
        if not line:
            break
        cols=line.strip().split("\t")
        TG,TF,zscore,corr,type,species=cols
        if zscore=="zscore":continue
        if TF not in TF2TG:
            TF2TG[TF]=defaultdict(list)
        TF2TG[TF][TG].append([TF,TG,zscore,corr,type])
    for TF in TF2TG:
        for TG in TF2TG[TF]:
            z_av=np.mean([float(i[2]) for i in TF2TG[TF][TG]])
            c_av=np.mean([float(i[3]) for i in TF2TG[TF][TG]])
            TF2TG[TF][TG]=[TF,TG,z_av,c_av]
    print >> sys.stderr, 'time cost for load_GO2flag:', time()-start_time
    return TF2TG

random.seed(1)
TF2TG=load_GRN(options.input)
gene2exp=loadExp(options.exp)
markers=["H3K4me1","H3K4me2","H3K4me3","H3K9ac","H3K27ac"]
CIGs=combineWidth(options.dir,markers,options.num,1)
gene_count=effective_number(TF2TG,CIGs)
exps,ups,downs=selectGeneExp(gene2exp,"HUVEC",gene_count)
index2name={0:"BHMs",1:"HighExp",2:"Up",3:"Down"}
index=0
for genes in [CIGs]:
    a,b,c=edgeCount(genes,TF2TG)
    print "\t".join(map(str,(index2name[index],a,0,b,0,c,0)))
    index+=1
SIM_l1,SIM_l2,SIM_l3=simulator_GRN(TF2TG,gene_count,10000)
