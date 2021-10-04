import os

#Program to take the output from STAR alignment (counts per gene) and create input files for DESeq2.
#
#This version assumes that the RNA-Seq data is not stranded.

sampleListAll = ["Sample_RNA140728LD_5_PCB429", "Sample_RNA140728LD_6_PCB439", "Sample_RNA140728LD_LD5_151", \
              "SJOS001103_D1", "SJOS001104_M1", "SJOS001105_D1", "SJOS001107_M1", "SJOS001108_M1", \
              "SJOS001109_D1", "SJOS001111_M1", "SJOS001112_D1", "SJOS001113_M1", "SJOS001116_M1", \
              "SJOS001120_D1", "SJOS001129_M2", "SJOS006_D", "SJOS007_D", "SJOS008_D", "SJOS010_D", \
              "SJOS011_D", "SJOS012_M", "SJOS013_D", "SJOS014_D", "SJOS017_D", "SJOS018_D", "SJOS019_D"]

sampleListStJude = ["SJOS001103_D1", "SJOS001104_M1", "SJOS001105_D1", "SJOS001107_M1", "SJOS001108_M1", \
              "SJOS001109_D1", "SJOS001111_M1", "SJOS001112_D1", "SJOS001113_M1", "SJOS001116_M1", \
              "SJOS001120_D1", "SJOS001129_M2", "SJOS006_D", "SJOS007_D", "SJOS008_D", "SJOS010_D", \
              "SJOS011_D", "SJOS012_M", "SJOS013_D", "SJOS014_D", "SJOS017_D", "SJOS018_D", "SJOS019_D"]

sampleListOHSU = ["Sample_RNA140728LD_5_PCB429", "Sample_RNA140728LD_6_PCB439", "Sample_RNA140728LD_LD5_151"]

#These are all the genes in all of the 5 nucleopore complex related pathways (Union).  Not all of these
#are intrinsic genes.
GenesOfInterestU = ['NP', 'rev', 'NUP160', 'NS', 'RAN', 'NUP155', 'NUP205', 'POM121C', 'NUP93', \
              'SLC25A6', 'GCK', 'M', 'NUP153', 'gag-pol', 'SEC13', 'gag', 'NUP37', 'CCNB1', \
              'NUP42', 'NEK7', 'PB2', 'NUP35', 'RAE1', 'NEK6', 'NUP54', 'CCNB2', 'PA', \
              'NUP188', 'AAAS', 'NEK9', 'NUP50', 'NUP85', 'CDK1', 'GCKR', 'SLC25A5', 'NUP62', \
              'PB1', 'HMGA1', 'BANF1', 'KPNA1', 'TPR', 'SLC25A4', 'NUP88', 'NUP43', 'NUP214', \
              'RANBP2', 'vif', 'NUP107', 'PSIP1', 'POM121', 'XPO1', 'NUP133', 'vpu', 'NDC1', 'NUP210', 'vpr']

def readCounts(WORKING_DIR, SAMPLE, D):
    INFILE = WORKING_DIR + "/RNA-Seq/" + SAMPLE + "/STAR_alignedReadsPerGene.out.tab"
    fi = open(INFILE, 'r')
    for i in fi:
        if i[0:4] == "ENSG":
            i = i.split()
            geneID = i[0]
            geneCounts = i[1]
            if geneID not in D:
                D[geneID] = [geneCounts]
            else:
                D[geneID].append(geneCounts)
        else:
            print(i)
    fi.close()
    return(D)

def checkD(D, L):
    print(len(D))
    for i in D:
        if len(D[i]) != len(L):
            print("Problemomo!!")
            print(D[i])
    print("Finished checking D.")

def makeCountData(WORKING_DIR, D, sampleList, L):
    sl = []
    print("Samples with available mutation data:")
    print(L)
    OUTFILE = WORKING_DIR + "/RNA-Seq/countData.txt"
    fo = open(OUTFILE, 'w')
    for x in sampleList:
        if x in L:
            sl.append(x)
    fo.write("ensgene" + '\t' + '\t'.join(sl) + '\n')
    for i in D:
        fo.write(i + '\t' + '\t'.join(D[i]) + '\n')
    fo.close()

def makeMetaData(WORKING_DIR, sampleList, genesOfInterest):
    CONS = ["3_prime_UTR_variant", "splice_acceptor_variant", "TF_binding_site_variant", \
            "protein_altering_variant", "splice_donor_variant", "stop_lost", "5_prime_UTR_variant", \
            "splice_region_variant", "stop_gained", "coding_sequence_variant", "start_lost", \
            "frameshift_variant", "regulatory_region_variant"]
    OUTFILE = WORKING_DIR + "/RNA-Seq/metaData_U.txt"
    fo = open(OUTFILE, 'w')
    fo.write('\t'.join(["sampleID", "Group", "sampleType"]) + '\n')
    L = set()#this set will contain the samples that have mutation data available
    D = {}
    for x in sampleList:
        sample = x
        group = "nonMutated"
        input_file = WORKING_DIR + "/" + sample + "/VEP_prefiltered.txt"
        if os.path.isfile(input_file):
            L.add(sample)
            fi = open(input_file, 'r')
            print("reading from " + input_file)
            for l in fi:
                h = []
                if l[0] != "#":
                    i = l.split()
                    if len(i) != 14:
                        print("Abnormal line :" + l)
                    else:
                        Feature_type = i[5]
                        Consequence = i[6]
                        Consequence = Consequence.split(",")
                        Extra = i[13]
                        Extra = Extra.split(";")
                        for u in Extra:
                            u = u.split("=")
                            h.append(u[0])
                            if u[0] == "SYMBOL":
                                SYMBOL = u[1]
                            elif u[0] == "NEAREST":
                                NSYMBOL = u[1]
                                NSYMBOL = NSYMBOL.split(",")
                        print(h)
                        for c in Consequence:
                            if c in CONS:
                                if "NEAREST" in h:
                                    if "SYMBOL" in h:
                                        if SYMBOL not in NSYMBOL:
                                            print("Boy howdy! we got issues: " + SYMBOL)
                                            print(NSYMBOL)
                                    if len(NSYMBOL) > 1:
                                        print("This " + c + " variant is nearby more than one gene")
                                        print(NSYMBOL)
                                    for n in NSYMBOL:
                                        print("Evaluating gene " + n)
                                        if n in genesOfInterest:
                                            print(n + " was on the list")
                                            group = "Mutated"
                                else:
                                    print("Memphis, we have a problem: " + l)
            sampleInfo = sample.split("_")
            if sampleInfo[1][0] == "D":
                sampleType = "Diagnostic"
            elif sampleInfo[1][0] == "M":
                sampleType = "Metastatic"
            else:
                print("This is bonkers!")
                print(sampleInfo)
                sampleType = "NA"
            if sample not in D:
                D[sample] = [group, sampleType]
            else:
                print("Fishy pants!! " + sample + " has been evaluated more than once.")
            fi.close()
        else:
            print(sample + " doesn't have mutation data.")
    for i in D:
        fo.write(i + '\t' + '\t'.join(D[i]) + '\n')
    fo.close()
    return(L)

WORKING_DIR = "/home/groups/jjacobs/lustreData/osteo"
L = makeMetaData(WORKING_DIR=WORKING_DIR, sampleList=sampleListStJude, genesOfInterest=GenesOfInterestU)
D = {}
for s in sampleListStJude:
    if s in L:
        D = readCounts(WORKING_DIR=WORKING_DIR, SAMPLE=s, D=D)
checkD(D=D, L=L)
makeCountData(WORKING_DIR=WORKING_DIR, D=D, sampleList=sampleListStJude, L=L)
