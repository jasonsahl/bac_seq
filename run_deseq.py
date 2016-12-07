#!/usr/bin/env python

"""DESeq controller"""

from optparse import OptionParser
import sys
import os
import subprocess


def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def get_conds(conditions,table):
    cond_dict = {}
    condition_fields = []
    for line in open(conditions, "U"):
        newline = line.strip()
        fields = newline.split()
        cond_dict.update({fields[0]:fields[1]})
    """Read in first line of the table"""
    firstLine = open(table).readline()
    first_fields = firstLine.split()
    for first_field in first_fields:
        if first_field in cond_dict:
            condition_fields.append(cond_dict.get(first_field))
    return condition_fields

def get_loci(merged_table):
    loci = []
    loci.append("")
    infile = open(merged_table, "U")
    firstLine = infile.readline()
    for line in infile:
        fields = line.split()
        loci.append(fields[0])
    return loci
    infile.close()

def add_loci_to_output(loci, input):
    outfile = open("tmp.txt", "w")
    print >> outfile, "\n".join(loci)
    outfile.close()
    os.system("paste tmp.txt %s > out.txt.xzx" % input)
    os.system("rm tmp.txt")

def add_in_info_and_filter(infile, p_value, prefix_1, prefix_2):
    outfile = open("%s_%s_significance.txt" % (prefix_1, prefix_2), "w")
    infile = open(infile, "U")
    firstLine = infile.readline()
    print >> outfile, firstLine,
    for line in infile:
        fields = line.split()
        if float(fields[-1])<=float(p_value):
            print >> outfile, line,
        else:
            pass
    outfile.close()
    infile.close()

def main(merged_table, conditions, R_file, p_value):
    table_path = os.path.abspath("%s" % merged_table)
    conds_path = os.path.abspath("%s" % conditions)
    R_path = os.path.abspath("%s" % R_file)
    """The conditions file should look like: sample metadata"""
    conds = get_conds(conds_path,table_path)
    loci = get_loci(table_path)
    unique_conds = [x for i, x in enumerate(conds) if x not in conds[i+1:]]
    tmp_conditions = open("conditions.txt.xzx", "w")
    tmp_conditions.write("\t".join(conds)+"\n")
    tmp_conditions.close()
    for i in range(len(unique_conds)):
        print "processing %s compared to %s" % (unique_conds[i], unique_conds[i-1])
        subprocess.check_call("Rscript %s %s %s %s %s > /dev/null 2>&1" % (R_path,table_path,"conditions.txt.xzx",unique_conds[i],unique_conds[i-1]), shell=True)
        add_loci_to_output(loci, "all_results.txt")
        os.system("awk '!/NA/' out.txt.xzx | sort -g -k 10,10 > %s_%s_all_sorted.txt" % (unique_conds[i],unique_conds[i-1]))
        add_in_info_and_filter("%s_%s_all_sorted.txt" % (unique_conds[i],unique_conds[i-1]), p_value, unique_conds[i], unique_conds[i-1])
        os.system("rm out.txt.xzx all_results.txt")


if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--merged_table", dest="merged_table",
                      help="path to merged table [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-c", "--conditions", dest="conditions",
                      help="path to conditions file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--R file", dest="R_file",
                      help="path to R file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--p_value", dest="p_value",
                      help="minimum p value, defaults to 0.05",
                      action="store", default="0.05", type="float")
    options, args = parser.parse_args()

    mandatories = ["merged_table", "conditions", "R_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.merged_table,options.conditions,options.R_file,
         options.p_value)
