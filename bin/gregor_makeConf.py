import argparse
import glob
import sys
import os


def makeOutputDirectory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        print("Caution: Directory already exists!!")
        
        
def makeIndexBedFile(bedfilename, annotfiles, annotFileDir):
    if annotfiles is not None:
        with open(bedfilename,'w') as f:
            for bedFile in annotfiles:
                f.write("{path}\n".format(path=bedFile))
                
    elif annotFileDir is not None:
        for bedFile in glob.glob(bedFileDirectory,"*.bed"):
            f.write("{path}\n".format(path=bedFile))
            
    else:
        print("Please provide annotation files or directory with annotation bed files")
        quit()        
        
def getopts():
    parser = argparse.ArgumentParser(description='Generate GREGOR conf file')
    parser.add_argument('--conffile', type=str, help="""Name of output .conf file""")
    parser.add_argument('--bedfilename', type=str, help="""Name of file with paths to annotfiles. If this does not exist, it will be created using either --annotfiles or --annotFileDir arguments.""")
    parser.add_argument('--annotfiles', nargs='+', type=str, help="""file with paths of annotation files to calculate enrichment in.""")
    parser.add_argument('--annotFileDir', type=str, help="""The directory with (unzipped) bed files if --annotfiles are not specified. Will take all .bed files in this directory""")
    parser.add_argument('--snpfile', type=str, help="""chr:pos formatted SNP files. IMPORTANT: filename should not contain '_'""")
    parser.add_argument('--refdir', type=str, default='/lab/data/sw/GREGOR/1.2.1', help=""" GREGOR reference. Default = /lab/data/sw/GREGOR/1.2.1 """)
    parser.add_argument('--population', type=str, default='EUR', help=""" Population. Default = EUR """)
    parser.add_argument('-r2', '--gregorR2Threshold', default='0.8', help="""minimum LD r2 for proxy SNPs (default: 0.8)""")
    parser.add_argument('--ldwindow', type=int, default=1000000, help="""minimum LD r2 for proxy SNPs (default: 0.8)""")
    parser.add_argument('-d','--outputdir', type=str, help="""Name of the output directory""")
    parser.add_argument('--neighbor', type=int, default=500, help="""Min neighbors for GREGOR""")
    parser.add_argument('--issorted', type=str, default="true", help="""true or false if bed files are sorted""")
    parser.add_argument('--topnbed', type=int, default=2, help="""Top n bed""")
    parser.add_argument('--cores', type=str, default='6', help="""Number of cores for each GREGOR job. (default: 6)""")
    return parser

if __name__ == '__main__':

    parser = getopts()
    args = parser.parse_args()
    
    makeOutputDirectory(args.outputdir)

    if not os.path.isfile(args.bedfilename):
        "Bed File being created - "
        makeIndexBedFile(args.bedfilename, args.annotfiles, args.annotFileDir)

    with open(args.conffile, 'w') as f:
        f.write("INDEX_SNP_FILE = {snpfile} \n".format(snpfile=args.snpfile))
        f.write("BED_FILE_INDEX = {bedfilename} \n".format(bedfilename=args.bedfilename))
        f.write("REF_DIR = {refdir} \n".format(refdir=args.refdir))
        f.write("POPULATION = {population} \n".format(population=args.population))
        f.write("R2THRESHOLD = {gregorR2Threshold} \n".format(gregorR2Threshold=args.gregorR2Threshold))
        f.write("LDWINDOWSIZE = {ldwindow} \n".format(ldwindow=args.ldwindow))
        f.write("OUT_DIR = {outputdir} \n".format(outputdir=args.outputdir))
        f.write("MIN_NEIGHBOR_NUM = {neighbors} \n".format(neighbors=args.neighbor))
        f.write("BEDFILE_IS_SORTED = {issorted} \n".format(issorted=args.issorted))
        f.write("TOPNBEDFILES = {topnbed} \n".format(topnbed=args.topnbed))
        f.write("JOBNUMBER = {cores} \n".format(cores=args.cores))
