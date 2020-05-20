#! /bin/python

################
#This scripts counts and makes a table for my own experiment structure.
#It should be adapted for other people I guess
################
#Jose Espejo Valle-Inclan
#2020
#jespejov@umcutrecht.nl/jaesvi@gmail.com




from __future__ import division
from sys import argv
import glob

def vcf_counter(filename, pattern):
    with open(filename, 'r') as f:
        t = 0
        p = 0
        for l in f:
            if l.startswith('#'):
                continue
            fil = l.split('\t')[6]
            if pattern in fil:
                p += 1
            t += 1
    return p, t
            

directory=argv[1]
masterD = {}
print '\t'.join(['tech', 'coverage', 'purity', 'type', 'tp', 'fp', 'fn'])
for tech in ['illumina']:
    masterD[tech] = {}
    sampledirs = glob.glob(directory+"/"+tech+"_[0-9]+/")
    #sampledirs += glob.glob(directory+"/"+tech+"/controls/")
    nsamples = len(sampledirs)
    for d in sampledirs:
        n = d.split('/')[-2].split('_')[-1]
        truthfiles = glob.glob(d+"/truth.purity*.vcf")
        for f in truthfiles:
            purity = f.split('.')[-2]
            if not purity in masterD[tech]:
                masterD[tech][purity]= {'raw':[], 'somatic':[]}
            p = purity.replace("purity", "")
            tsetsom, tsettotal = vcf_counter(f,p+"_SOM")
            tpsom = str(tsetsom)
            fnsom = str(tsettotal - tsetsom)
            
            tsetraw, tmp = vcf_counter(f,p+"_RAW")
            tpraw = str(tsetraw)
            fnraw = str(tsettotal - tsetraw)
            
            if tech == "nanopore":
                somfile=d+"/%s.somatic.truth.vcf" % (purity)
                rawfile=d+"/%s.merged.truth.vcf" % (purity)
            elif tech == "pacbio":
                somfile=d+"/%s.somatic.truth.uniquified.vcf" % (purity)
                rawfile=d+"/%s.merged.truth.uniquified.vcf" % (purity)
            else:
                somfile=d+"/%s.somatic.truth.uniquified.vcf" % (purity)
                rawfile=d+"/%s.raw.truth.uniquified.vcf" % (purity)
                
            somtrue, somtotal = vcf_counter(somfile,"TRUTHSET")
            rawtrue, rawtotal = vcf_counter(rawfile,"TRUTHSET")
            
            fpsom = str(somtotal - somtrue)
            fpraw = str(rawtotal - rawtrue)

            
            #tsetsom = str(round(tsetsom*100/tsettotal, 2))
            #tsetraw = str(round(tsetraw*100/tsettotal, 2))

            masterD[tech][purity]['raw'].append(tsetraw)
            masterD[tech][purity]['somatic'].append(tsetsom)
            print '\t'.join([tech, n , purity, 'raw', tpraw, fpraw, fnraw])
            print '\t'.join([tech, n , purity, 'somatic', tpsom, fpsom, fnsom])
            

