# This code is for process PCAWG ICGC patient coverage files. The output files will be merged into tissue type coverage for MutSigCV
# Author: Yiyun Rao
# Date last edit:05222020

from optparse import OptionParser
import os
import re
import time


# Error Area
class TumorTypeNotFoundError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'TumorTypeNotFoundError, {0} '.format(self.message)
        else:
            return 'TumorTypeNotFoundError has been raised'


class PathAlreadyExistsError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return 'PathAlreadyExists, {0} '.format(self.message)
        else:
            return 'PathAlreadyExists has been raised'


def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-c", "--coverage_file",
                      dest="inputPath",
                      help="The full absolute path of the PCAWG coverage file")
    parser.add_option("-i", "--info_file",
                      dest="inputPath",
                      help="The full absolute path of the PCAWG info file")

    (options, arg) = parser.parse_args()
    if not options.inputPath:
        parser.error("Requires at least one input file, the path of this script. Run with --help for help.")
        quit(1)
    return options


def search_tumor_type_column(coverage_file):
    # Read the first line and find tumor type column(project_code)
    with open(coverage_file, 'r') as file_coverage:
        first_line = file_coverage.readline()
        line_split = first_line.split("\t")
        idx = [line_split.index(cols) for cols in line_split if re.findall("project_code", cols, flags=re.IGNORECASE)]
        if idx:
            if len(idx) == 1:
                tumor_col = int(idx[0])
            else:
                raise TumorTypeNotFoundError("MULTIPLE COLUMNS WERE MAPPED")
        else:
            raise TumorTypeNotFoundError("TUMOR TYPE COLUMN NOT FOUND")

        # Get current directory and create tumor-type directory
        cwd = os.getcwd()
        new_path = cwd + "/merged_tumor_type_maf"
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        else:
            raise PathAlreadyExistsError("DIRECTORY ALREADY EXISTS")

        # Parse the MAF file, create new file according to tumor type and put corresponding patients in it
        next(file_coverage)
        for lines in file_coverage:
            line_split = lines.split("\t")
            tumor_type = line_split[tumor_col]
            file_name = "merged_consensus_passonly_snv_mnv_indel_" + str(tumor_type) + ".maf"
            # IF file doesn't exist, create file and write the header line
            if not os.path.exists(new_path + '/' + file_name):
                with open(new_path + '/' + file_name, 'w') as f:
                    f.write(first_line)
                f.close()
            # Append the file to
            with open(new_path + '/' + file_name, "a") as f1:
                f1.write(lines)
            f1.close()


# Main function
if __name__ == "__main__":
    arguments = parse_arguments()
    maf_file = arguments.inputPath
    start1 = time.time()
    search_tumor_type_column(maf_file)
    end1 = time.time()
    print "TIME USED " + str(end1 - start1) + " seconds"
