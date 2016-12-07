#!/usr/bin/env python

"""RNA-seq pipeline for bacterial
transcripts"""

from optparse import OptionParser
import sys
import os
import subprocess
import glob
import re
import logging
from subprocess import Popen
try:
    from igs.utils import logging as log_isg
    from igs.threading import functional as p_func
except:
    print "PATH is not configured correctly"
    sys.exit()


BACSEQ_PATH="/Users/jasonsahl/tools/bac_seq"
if os.path.exists(BACSEQ_PATH):
    sys.path.append("%s" % BACSEQ_PATH)
else:
    print "your BACSEQ_PATH path is not correct.  Edit the path in bac_seq.py and try again"
    sys.exit()
TRIM_PATH=BACSEQ_PATH+"/bin/trimmomatic-0.30.jar"

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastqs cannot be found"
        sys.exit()
def get_readFile_components(full_file_path):
    """function adapted from:
    https://github.com/katholt/srst2 - tested"""
    (file_path,file_name) = os.path.split(full_file_path)
    m1 = re.match("(.*).gz",file_name)
    ext = ""
    if m1 != None:
        ext = ".gz"
        file_name = m1.groups()[0]
    (file_name_before_ext,ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext
    return(file_path,file_name_before_ext,full_ext)

def read_file_sets(dir_path):
    """match up pairs of sequence data, adapted from
    https://github.com/katholt/srst2 - untested"""
    fileSets = {}
    forward_reads = {}
    reverse_reads = {}
    num_paired_readsets = 0
    num_single_readsets = 0
    for infile in glob.glob(os.path.join(dir_path, "*.fastq.gz")):
        (file_path,file_name_before_ext,full_ext) = get_readFile_components(infile)
        m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
        if m==None:
            m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
            if m!=None:
                (baseName,read) = m.groups()[0], m.groups()[1]
                forward_reads[baseName] = infile
            else:
                m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
                if m!=None:
                    (baseName,read) = m.groups()[0], m.groups()[1]
                    reverse_reads[baseName] = infile
                else:
                    print "Could not determine forward/reverse read status for input file "
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                print "Could not determine forward/reverse read status for input file "
                fileSets[file_name_before_ext] = infile
                num_single_readsets += 1
    for sample in forward_reads:
        if sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + reverse_reads[sample])

    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets))
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))

    return fileSets

def run_trimmomatic(trim_path, processors, read_1, read_2, name, bac_path):
    args=['java','-jar','%s' % trim_path,'PE', '-threads', '%s' % processors,
                  '%s' % read_1, '%s' % read_2, '%s.F.paired.fastq.gz' % name, 'F.unpaired.fastq.gz',
	          '%s.R.paired.fastq.gz' % name, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % bac_path,
	          'MINLEN:50']
    try:
        vcf_fh = open('%s.trimmomatic.out' % name, 'w')
    except:
        log_isg.logPrint('could not open trimmomatic file')
    try:
        log_fh = open('%s.trimmomatic.log' % name, 'w')
    except:
        log_isg.logPrint('could not open log file')
    try:
        trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
        trim.wait()
    except:
        log_isg.logPrint("problem encountered with trimmomatic")

def bwa(reference,read1,read2,sam_file, processors, log_file='',**my_opts):
    """controller for bwa, currently only works with bwa mem - untested because they are mostly system calls"""
    mem_arguments = ['bwa', 'mem', '-v', '2', '-M', '-t', '%s' % processors]
    for opt in my_opts.items():
        mem_arguments.extend(opt)
    if "NULL" in read2:
        mem_arguments.extend([reference,read1])
    else:
        mem_arguments.extend([reference,read1,read2])
    if log_file:
       try:
           log_fh = open(log_file, 'w')
       except:
           print log_file, 'could not open'
    else:
        log_fh = PIPE
    try:
        sam_fh = open(sam_file, 'w')
    except:
        print sam_file, 'could not open'

    bwa = Popen(mem_arguments, stderr=log_fh, stdout=sam_fh)
    bwa.wait()

def run_bwa(reference, read1, read2, processors, name):
    """launches bwa. Adds in read_group for compatability with GATK - untested"""
    read_group = '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tPU:%s' % (name,name,name)
    bwa(reference,read1, read2,"%s.sam" % name, processors, log_file='%s.sam.log' % name,**{'-R':read_group})

def run_htseq(name, gff):
    subprocess.check_call("htseq-count -t gene -i locus_tag %s.sam %s -m intersection-nonempty -a 0 -q > %s.counts.txt 2> /dev/null" % (name,gff,name) , shell=True)

def run_loop(fileSets, dir_path, reference, processors, trim_path, bac_path, gff):
    files_and_temp_names = [(str(idx), list(f)) for idx, f in fileSets.iteritems()]
    names = [ ]
    def _perform_workflow(data):
        """idx is the sample name, f is the file dictionary"""
        idx, f = data
        if os.path.isfile("%s.sam" % idx):
            pass
        else:
            run_trimmomatic(trim_path, processors, f[0], f[1], idx, bac_path)
            run_bwa(reference, '%s.F.paired.fastq.gz' % idx, '%s.R.paired.fastq.gz' % idx, processors, idx)
        if os.path.isfile("%s.counts.txt" % idx):
            pass
        else:
            run_htseq(idx, gff)
        names.append(idx)
    set(p_func.pmap(_perform_workflow,
                    files_and_temp_names,
                    num_workers=processors))
    return names

def run_kallisto_loop(fileSets,dir_path,reference,processors,bac_path):
    files_and_temp_names = [(str(idx), list(f)) for idx, f in fileSets.iteritems()]
    names = []
    def _perform_workflow(data):
        """idx is the sample name, f is the file dictionary"""
        idx, f = data
        subprocess.check_call("kallisto quant -o %s -i %s %s %s > /dev/null 2>&1" % (idx,"REF",f[0],f[1]), shell=True)
        names.append(idx)
    set(p_func.pmap(_perform_workflow,
                    files_and_temp_names,
                    num_workers=processors))
    return names

def get_seq_name(in_fasta):
    """used for renaming the sequences - tested"""
    return os.path.basename(in_fasta)

def create_merged_table(locus_ids, start_dir):
    merged_list = [ ]
    merged_list.append(locus_ids)
    out_table = open("merged_table.txt", "w")
    for infile in glob.glob(os.path.join(start_dir, "*.counts.txt")):
        sample_list = [ ]
        full_name = get_seq_name(infile)
        name = full_name.replace(".counts.txt", "")
        sample_list.append(name)
        for line in open(infile, "U"):
            if line.startswith("__"):
                pass
            else:
                fields = line.split()
                sample_list.append(int(fields[1]))
        merged_list.append(sample_list)
    test=map(list, zip(*merged_list))
    for x in test:
        print >> out_table, "\t".join(x)
    out_table.close()

def main(read_dir,reference,gff,aligner,processors):
    start_dir = os.getcwd()
    log_isg.logPrint('testing the paths of all dependencies')
    ap=os.path.abspath("%s" % start_dir)
    """need to check the paths for all dependencies"""
    if aligner == "bwa-mem":
        ac = subprocess.call(['which','bwa'])
        if ac == 0:
            pass
        else:
            print "bwa must be in your path"
            sys.exit()
        ad = subprocess.call(['which', 'htseq-count'])
        if ac == 0:
            pass
        else:
            print "htseq-count must be in your path"
            sys.exit()
    elif aligner == "kallisto":
        ac = subprocess.call(['which','kallisto'])
        if ac == 0:
            pass
        else:
            print "kallisto must be in your path"
            sys.exit()
    log_isg.logPrint('bac_seq pipeline starting')
    ref_path=os.path.abspath("%s" % reference)
    dir_path=os.path.abspath("%s" % read_dir)
    if aligner == "bwa-mem":
        gff_path=os.path.abspath("%s" % gff)
        subprocess.check_call("bwa index %s > /dev/null 2>&1" % (ref_path), shell=True)
    elif aligner == "kallisto":
        subprocess.check_call("kallisto index %s -i REF > /dev/null 2>&1" % ref_path, shell=True)
    fileSets=read_file_sets(dir_path)
    if aligner == "bwa-mem":
        names = run_loop(fileSets,dir_path,"%s" % ref_path,processors,TRIM_PATH,BACSEQ_PATH,gff_path)
        """print out names, which will be important for the next step"""
        outfile = open("names.txt", "w")
        sample_name = names[0]
        locus_ids = []
        locus_ids.append("")
        for line in open("%s.counts.txt" % sample_name):
            if line.startswith("__"):
                pass
            else:
                fields = line.split()
                locus_ids.append(fields[0])
        create_merged_table(locus_ids, start_dir)
        #Create names file in the same order that they are in the merged_table
        infile = open("merged_table.txt", "U")
        firstLine = infile.readline()
        FL_F=firstLine.split()
        print >> outfile, "\n".join(FL_F)
        log_isg.logPrint('bac_seq finished')
        try:
            subprocess.check_call("rm *.fasta.* *.log *.trimmomatic.out *.paired.fastq.gz *.unpaired.fastq.gz", shell=True, stderr=open(os.devnull, 'w'))
        except:
            pass
    else:
        log_isg.logPrint("running kallisto")
        names=run_kallisto_loop(fileSets,dir_path,"REF",processors,BACSEQ_PATH)
        #outfile=open("kallisto_merged_counts.txt", "w")
        #Now I need to create the same matrix that comes out of BWA-MEM
        count_dir = ()
        for name in names:
            for line in open("%s/abundance.tsv" % name):
                newline=line.strip()
                if line.startswith("target_id"):
                    pass
                else:
                    fields=newline.split()
                    count_dir=((name,fields[0],float(fields[3])),)+count_dir
        log_isg.logPrint("processing results")
        names.insert(0,"")
        #outfile.write("\t".join(names)+"\n")
        #outfile.close()
        marker_list = []
        for entry in count_dir:
            marker_list.append(entry[1])
        #print count_dir
        nr=[x for i, x in enumerate(marker_list) if x not in marker_list[i+1:]]
        """I will need to make sure that these are ordered in the same way as they are in marker_list"""
        ref_out = open("ref.list", "w")
        ref_out.write(""+'\n')
        for entry in nr:
            ref_out.write(entry+"\n")
        ref_out.close()
        for name in names:
            outfile = open("%s.xyx" % name, "w")
            values = []
            for entry in nr:
                for atuple in count_dir:
                    if entry == atuple[1] and name == atuple[0]:
                        values.append(int(atuple[2]))
            values.insert(0,name)
            outfile.write("\n".join(values))
            outfile.close()
        os.system("paste ref.list *xyx > kallisto_merged_counts.txt")
        os.system("rm ref.list *xyx")
            #            genome_hits.append(str(atuple[2]))
            #    hits.append(genome_hits)
            #for hit in hits:
            #    outfile.write("\t".join(hit)+"\n")
        #outfile.close()
        #os.system("paste ref.list results.tmp.xyx > body.xyx")
        #os.system("cat kallisto_merged_counts.txt body.xyx > combined.xyx")
        #os.system("mv combined.xyx kallisto_merged_counts.txt")
        #os.system("rm ref.list results.tmp.xyx body.xyx")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--read_dir", dest="read_dir",
                      help="directory of FASTQ reads [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-r", "--reference", dest="reference",
                      help="reference genome in FASTA format [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-g", "--gff_file", dest="gff",
                      help="GFF file for reference genome, only required for BWA-MEM",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-a", "--aligner", dest="aligner",
                      help="aligner to use [kallisto, bwa-mem]; defaults to kallisto",
                      action="store", default="kallisto", type="string")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to use, defaults to 4",
                      action="store", default="4", type="int")
    options, args = parser.parse_args()

    mandatories = ["read_dir","reference"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.read_dir,options.reference,options.gff,options.aligner,options.processors)
