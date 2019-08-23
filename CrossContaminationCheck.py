import subprocess
import os
import sys
import pprint as pp
import argparse
import random

import collections



SeqRec= collections.namedtuple('SeqRec', 'id label seq')
PairwiseAlignment= collections.namedtuple('PairwiseAlignment', 'label1 label2 seq1 seq2 score similarity length gaps substNum gapNum matchNum alnStr')

validBases={'A':True,'C':True,'G':True,'T':True }


#parse fasta text into SecRec named tuples
def parseFastaStr(fastaStr):
    txtBlocks=fastaStr.rstrip().split(">")[1:]
    seqRecords=[]
    seqID=0
    for block in txtBlocks:
        label,seq=block.rstrip().split("\n")
        seqID+=1
        seqRecords.append(SeqRec(id=str(seqID), label=label, seq= seq))
    
    return seqRecords


#parses Emboss needleman-wunsch output file 
#with support for multi-alignment
def parseEmbossNeedleOutput(needleEmbossOutputFile, labelLookup={}):
    outputAlignment={}
    with open(needleEmbossOutputFile, 'r') as results:
        results=results.read()
        results=results.replace( "#---------------------------------------\n","" )
        results=results.split("########################################")[2]

        for l in results.strip().split("\n\n\n"):
    
            s,aln=l.strip().split("#=======================================")[1:]


            for line in s.strip().split("\n"):
                if line[0] == "#":
                    try:
                        if line.split()[1] == "Similarity:":
                            if line.split()[-1][0] == '(':
                                similarity = float(line.split()[-1][1:-2].strip())
                            else:
                                similarity = float(line.split()[-1][:-2].strip())
                        if line.split()[1] == "Length:":
                            length = int(line.split()[2])
                        if line.split()[1] == "Score:":
                            score = int(float(line.split()[2]))
                        if line.split()[1] == "Gaps:":
                            if line.split()[-1][0] == '(':
                                gaps = float(line.split()[-1][1:-2].strip())
                            else:
                                gaps = float(line.split()[-1][:-2].strip())
                        
                    except IndexError:
                        pass   
    
    
            seq1=[]
            seq2=[]

            for row in aln.strip().split("\n\n"):
                seqLine1,matchLine,seqLine2=row.split("\n")
                seqLine1=" ".join(seqLine1.split()).split()
                seqLine2=" ".join(seqLine2.split()).split()
                seq1.append(seqLine1[-2])
                seq2.append(seqLine2[-2])
                label1=" ".join(seqLine1[:-3])
                label2=" ".join(seqLine2[:-3])
    
        
            seq1="".join(seq1).upper()
            seq2="".join(seq2).upper()
            gapNum=0
            substNum=0
            matchNum=0
        
        
            for i in range(len(seq1)):
                if seq1[i]==seq2[i]: matchNum+=1
                elif seq1[i]=='-' or seq2[i]=='-': gapNum+=1
                elif validBases.has_key(seq1[i]) and validBases.has_key(seq2[i]): substNum+=1
                
            
            aln=PairwiseAlignment(\
                label1=labelLookup.get(label1,label1), \
                label2=labelLookup.get(label2,label2), \
                seq1=seq1, \
                seq2=seq2, \
                similarity=similarity, \
                score=score, \
                length=length, \
                gaps=gaps, \
                substNum=substNum, \
                gapNum=gapNum, \
                matchNum=matchNum,
                alnStr=aln.strip()
                )
           
            outputAlignment.setdefault(aln.label1,{}).setdefault(aln.label2,aln)

    return outputAlignment


#calls external needle executable from Emboss toolkit
def callEmbossNeedle(SEQA_LIST, gapopen = 30., gapextend = 2., scratchroot = '/tmp/'):
    WHICH_NEEDLE = "/usr/local/bin/needle"
    outputAlignment={}
    
    idToLabelLookup={}
    for SEQ in SEQA_LIST: idToLabelLookup[SEQ.id]=SEQ.label
    
    SEQB_LIST=SEQA_LIST

    for SEQA in SEQA_LIST[:-1]:
       
        SEQB_LIST=SEQB_LIST[1:]
        asequencen = os.path.join(scratchroot, "__scratch_needle_asequence.txt")
        bsequencen = os.path.join(scratchroot, "__scratch_needle_bsequence.txt")
        outfilen   = os.path.join(scratchroot, "__scratch_needle_outfile.txt"  )

        with open(asequencen, "w+") as asequenceh:
            print >>asequenceh, ">"+SEQA.id
            print >>asequenceh, SEQA.seq.replace('-','')
            print >>asequenceh

        with open(bsequencen, "w+") as bsequenceh:
            
     
            for SEQB in SEQB_LIST:
                print >>bsequenceh, ">"+SEQB.id
                print >>bsequenceh, SEQB.seq.replace('-','')
                print >>bsequenceh

        needle = subprocess.Popen([
            WHICH_NEEDLE,
            "-asequence", asequencen,
            "-bsequence", bsequencen,
            "-gapopen", str(float(gapopen)),
            "-gapextend", str(float(gapextend)),
    	    "-endweight","Y",
            "-outfile", outfilen,
            #"-datafile", subst                         # STUB: Substitution Matrix File
            ],
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
        )
        (needle_stdout, needle_stderr) = needle.communicate()
    
        if 'Error:' in needle_stderr or 'Died:' in needle_stderr:
        	raise Exception('FATAL: Needle failed with the following message:\n' + needle_stderr)
        
        outputAlignment.update(parseEmbossNeedleOutput(outfilen,idToLabelLookup))
    
    return outputAlignment




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes fasta files as input, performs pariwise alignments and identifies those pairs that fit a specified gap and substitution count envelope.')
    parser.add_argument('-v', action='store_true', help='verbose output (show alignments)')
    parser.add_argument("-minG", type=int, required=True, help="minimum number of gaps")
    parser.add_argument("-maxG", type=int, required=True, help="maximum number of gaps")
    parser.add_argument("-minS", type=int, required=True, help="minimum number of substitutions")
    parser.add_argument("-maxS", type=int, required=True, help="maximum number of substitutions")
    parser.add_argument("-input", type=str, required=True, help="input fasta file")
    
    args = parser.parse_args()
    
    print "Fasta cross-contaminat check for "+args.input
    print "Parameters: Gaps[",args.minG,":",args.maxG,"] Substitutions[",args.minS,":",args.maxS,"]"
    
    with open(args.input, "r") as inputFasta:
        seqs=parseFastaStr(inputFasta.read().strip())
        
 
    tempFolder="/tmp/"+str(random.randint(1000000,9999999))
    os.mkdir(tempFolder)
    aligned= callEmbossNeedle(seqs,scratchroot=tempFolder)
    
    for s1,s2s in aligned.items():
        for s2,aln in s2s.items():
            if aln.gapNum>=args.minG and aln.gapNum<=args.maxG and aln.substNum>=args.minS and aln.substNum<=args.maxS:
                print s1+" / "+s2,"Similarity:",aln.similarity,"Gaps:",aln.gapNum,"Subst:",aln.substNum
                if args.v: print aln.alnStr




