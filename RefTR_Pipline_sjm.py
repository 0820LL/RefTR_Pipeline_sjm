#!/usr/bin/env pyhton
#coding=utf-8
"""
Created by Lilin <lilin@novogene.cn>
date: 2014-09-02
"""
__author__ = 'lilin <lilin@novogene.cn>'
__Version__ = 'v1.0'

import os
import os.path
import re
import argparse
import ConfigParser

root_dir=os.getcwd()

################################# SJM BEGIN ####################################
class job(object):
    def __init__(self, jobname, memory, slot, shell, queue='tjnode'):
        self.jobname = jobname
        self.memory = memory
        self.slot = slot
        self.shell = shell
        self.queue = queue

    def sjm(self):
        if re.match(r'mem\d+', self.queue):
            P = str(self.queue) + '.q'
            txt = '''
job_begin
    name %s
    sched_options -P %s -q %s -cwd -V -l vf=%dG,p=%d
    cmd sh %s
job_end
''' % (self.jobname, self.queue, P, self.memory, self.slot, self.shell)
        else:
            txt = '''
job_begin
    name %s
    sched_options -q rna.q -q all.q -cwd -V -l vf=%dG,p=%d
    cmd sh %s
job_end
''' % (self.jobname, self.memory, self.slot, self.shell)
        return txt

class job1(object):
    def __init__(self, jobname, shell):
        self.jobname = jobname
        self.shell = shell

    def sjm(self):
        txt = '''
job_begin
    name %s
    cmd sh %s
job_end
''' % (self.jobname, self.shell)
        return txt

class LocalJob(job1):
    def sjm(self):
        txt = '''
job_begin
    name %s
    host localhost
    cmd sh %s
job_end
''' % (self.jobname, self.shell)
        return txt
    def sjm_py(self):
        txt = '''
job_begin
    name %s
    host localhost
    cmd python %s
job_end
''' %(self.jobname,self.shell)
        return txt
################################# SJM END #######################################
######################## check sample name BEGIN ################################
def checkSample(sample):
    PatChar = r'^(CON|PRN|AUX|CLOCK\$|NUL|COM[1-9]|LPT1|[\W]+)$'
    PatNum = r'^\d+'
    MAXLEN = 8
    message = '''
can only use word/digit/underscore in sample, start with word.
Max length is 8, no windows reserved words like CON PRN...
'''
    if re.search(PatChar, sample) or re.search(PatNum, sample) or len(sample) > MAXLEN:
        print '%s<==invalid name' %(sample)
        print message
        exit()
######################### check sample name END #################################
######################## parse the arguments BEGIN ##############################
parser = argparse.ArgumentParser(description="RefTR pipline v4.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--fq',help="the original directory of the raw fastq reads, [REQUIRED]",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",default=None)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list, ",default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',choices=['y','n'],default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,warning: order sensitive!',default=None)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['no','yes','reverse'],required=True)
parser.add_argument('--number',help="chromosome number, [REQUIRED for density plot]",required=True)
parser.add_argument('--length',help="the length of sequenced reads,[defalut=100]",default='100')
parser.add_argument('--fa',help="the reference FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--gtf',help="the annotation GTF file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1/sample2,sample3, [REQUIRED]",required=True,default=None)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",default=None)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",default=None)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",default=None)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: all species: /PUBLIC/database/RNA//kobas2.0-data-20120208/seq_pep_v2), [defalut=kaas]",default='kaas')
parser.add_argument('--ppi_number',help="species code, (ref file: /PUBLIC/database/RNA/string_ppi/species.v9.0.txt)",default=None)
parser.add_argument('--ppi_blast',help="whether to run blast to get the protein-protein interactions",choices=['y','n'],default=None)
parser.add_argument('--genenamefile',help="genenamefile, 1st col is geneID, 2nd col is genename",default=None)
parser.add_argument('--ex',help="the steps you do not wanna perform",default=None)
######################### parse the arguments END ################################
################### extract and check the parameters BEGIN #######################
argv = vars(parser.parse_args())
project=argv['project'].strip()
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
for s in samples:
    checkSample(s)
assert len(samples)==len(samples_tmp)
#-----------------------------------------------------------------------
mapfiles=[each.strip() for each in argv['mapfile'].strip().split(',') if each.strip() != '']
mapfile=' '.join(mapfiles)
assert not os.system('cat %s |sort -u |sed \'/^$/d\' >%s' % (mapfile,root_dir+'/libraryID'))
if argv['fq']:
    fq = argv['fq'].strip()
    fq = os.path.abspath(fq)
else:
    assert not os.system('mkdir raw_data')
    assert not os.system('perl /PUBLIC/source/RNA/RefRNA/ln_raw_data.pl %s %s pe raw_data' %(argv['mapfile'],argv['raw_dir']))
    fq = root_dir + '/raw_data'
for each in samples:
    fq_tmp1=fq+'/'+each+'_1.fq.gz'
    fq_tmp2=fq+'/'+each+'_2.fq.gz'
    assert os.path.isfile(fq_tmp1)
    assert os.path.isfile(fq_tmp2)
if argv['ad']:
    ad = argv['ad'].strip()
    ad = os.path.abspath(ad)
    for each in samples:
        ad_tmp1=ad+'/'+each+'_1.adapter.list.gz'
        ad_tmp2=ad+'/'+each+'_2.adapter.list.gz'
        assert os.path.isfile(ad_tmp1)
        assert os.path.isfile(ad_tmp2)
if argv['generate_adapter']:
    generate_adapter = argv['generate_adapter'].strip()
else:
    generate_adapter = 'n'
if argv['index']:
    index=argv['index'].strip()
    indexes=index.split(',')
    assert len(samples)==len(indexes)
if generate_adapter == 'y':
    if argv['ad'] !=None:
        print 'Error:  the parameters --ad and --generate_adapter are not consensus!\n'
        exit()
    if argv['index']==None:
        print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
        exit()
else:
    if argv['index'] != None:
        print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
        exit()
#-----------------------------------------------------------------------
all_content=set([1,2,3,4,5,6,7,8,9])
if argv['ex'] != None:
    excludes=argv['ex'].strip().strip(',').strip().split(',')
    excludes=[int(each.strip()) for each in excludes]
    for each1 in excludes:
        assert each1 in all_content
else:
    excludes=[] #list
if argv['group'].find(':') == -1:
    excludes.append(3)
if ( len(argv['group'].split(',') ) == 1) or (5 in excludes):
    excludes.extend([7,8,9])
includes=all_content-set(excludes)
#-----------------------------------------------------------------------
ss = argv['ss'].strip()
number = argv['number'].strip()
if argv['length']:
    length = argv['length'].strip()
else:
    length = '100'
fa = argv['fa'].strip()
fa = os.path.abspath(fa)
suffix_fa=['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
for each in suffix_fa:
    fa_tmp=fa+each
    assert os.path.isfile(fa_tmp)
gtf = argv['gtf'].strip()
gtf = os.path.abspath(gtf)
assert os.path.isfile(gtf)
#-----------------------------------------------------------------------
if (set([1]).issubset(includes)):
    if argv['group'] == None:
        groups=samples
        groups_iter=samples
        group=sample
        flag_repeat=False
    else:
        groups_iter=[]
        if ':' in argv['group'].strip(':'):
            flag_repeat=True
            groups=[each.strip().split(':') for each in argv['group'].strip().strip(':').split(',') if each.strip() != '']
            for each in groups:
                groups_iter+=each
            group_iter_n=[]
            for each in groups_iter:
                if each not in group_iter_n:
                    group_iter_n.append(each)
            groups_iter=group_iter_n
            group=','.join([':'.join(each) for each in groups])
        else:
            flag_repeat=False
            groups=[each.strip() for each in argv['group'].strip().split(',') if each.strip() != '']
            for each in groups:
                groups_iter.append(each)
            group=','.join(groups)
            assert set(groups_iter).issubset(samples)
        group_iter=','.join(groups_iter)

    if argv['groupname'] == None:
        if flag_repeat == False:
            groupname_tmp=groups
        else:
            groupname_tmp=['group'+str(k+1) for k in range(len(groups))]
    else:
        groupname_tmp=[each.strip() for each in argv['groupname'].split(',') if each.strip() != '']
        assert len(groupname_tmp) == len(groups)
        groupname=','.join(groupname_tmp)

    compare = argv['compare'].strip()
    compare_samples=[each.strip().split(':') for each in argv['compare'].strip().split(',') if each.strip() != '']
    M=[]
    for each1 in compare_samples:
        assert len(each1) == 2
        for each2 in each1:
            assert each2.isdigit()
            M.append(int(each2))
    assert max(M) <= len(groupname_tmp)
    assert min(M) > 0
    compare_sample=','.join([','.join(each) for each in compare_samples])
    temp2=[]
    for each1 in compare_samples:
        temp1=[]
        for each2 in each1:
            temp1.append(groupname_tmp[int(each2)-1])
        temp2.append(':'.join(temp1))
    compares_name=','.join(temp2)
    if argv['venn']:
        venn = argv['venn'].strip()
        com_pairs=compare.split(',')
        venn_clusters=[each.split('_') for each in venn.split(',')]
        temp1=[]
        for each1 in venn_clusters:
            temp2=[]
            for each2 in each1:
                assert each2 in com_pairs
                temp3=each2.split(':')
                assert len(temp3) == 2
                temp2.append(groupname_tmp[int(temp3[0])-1]+':'+groupname_tmp[int(temp3[1])-1])
            temp1.append('_'.join(temp2))
        venn_cluster_name=','.join(temp1)
    else:
        venn = compare.replace(',','_')
        venn_cluster_name=compares_name.replace(',','_')
    groups=group.split(',')
    groupnames=groupname.split(',')
    compares=compare.split(',')
    venns=venn.split(',')
    for each in groups:
        temp=each.split(':')
        assert len(temp)==len(list(set(temp)))
        for each2 in temp:
            assert each2 in samples
    assert len(groups)==len(groupnames)
    assert len(groupnames)==len(list(set(groupnames)))
    for i,each in enumerate(groupnames):
        group_tmp=groups[i]
        group_tmp=group_tmp.replace(':',',')
    compare_name=[]
    for i,each in enumerate(compares):
        compare_tmp=each
        compare_tmp=compare_tmp.split(':')
        assert len(compare_tmp)==2
        compare_name=groupnames[int(compare_tmp[0])-1]+'vs'+groupnames[int(compare_tmp[1])-1]
    for i,each in enumerate(venns):
        venn_tmp=each
        venn_tmp=venn_tmp.split('_')
        for each2 in venn_tmp:
            assert each2 in compares
            venn_tmp2=each2.split(':')
            venn_name=groupnames[int(venn_tmp2[0])-1]+'vs'+groupnames[int(venn_tmp2[1])-1]
#-----------------------------------------------------------------------
goann = argv['goann'].strip()
goann = os.path.abspath(goann)
species = ppi_number = ppi_blast = ''
assert os.path.isfile(goann)
if set([4,5,8]).issubset(includes):
    if argv['species']:
        species = argv['species'].strip()
    else:
        species = 'kaas'
if set([4,5,9]).issubset(includes):
    if argv['ppi_number']:
        ppi_number = argv['ppi_number'].strip()
    if argv['ppi_blast']:
        ppi_blast = argv['ppi_blast'].strip()
    if argv['ppi_blast']:
        if argv['ppi_number'] == None:
            print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
            exit()
    else:
        if argv['ppi_number']:
            print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
            exit()
if argv['genenamefile']:
    genenamefile = argv['genenamefile'].strip()
    genenamefile = os.path.abspath(genenamefile)
    assert os.path.isfile(genenamefile)
else:
    genenamefile = root_dir+'/Blast_TR/Blast_Swissprot/diffgene_union.genenames'
#################### extract and check the parameters END ########################
######################## display all parameters BEGIN ############################
display = open(root_dir + '/' + 'TR_command.txt','w')
display.write('project: %s\n' % (project))
display.write('sample: %s\n' % (sample))
display.write('fq: %s\n' % (fq))
display.write('adapter: ')
if argv['ad']:
    display.write('%s\n' % (argv['ad'].rstrip()))
if argv['generate_adapter']:
    generate_adapter = argv['generate_adapter']
    display.write('generate_adapter: %s\n' % (generate_adapter))
if argv['index']:
    indexes = argv['index'].rstrip().split(',')
    for i,index_tmp in enumerate(samples):
        display.write('%s:\t%s\n' % (index_tmp,indexes[i]))
        display.write('\n')
display.write('length: %s\n' % (length))
display.write('fa: %s\n' % (fa))
display.write('gtf: %s\n' % (gtf))
display.write('\ngroups:\n')
for i,each in enumerate(groupnames):
    group_tmp=groups[i]
    group_tmp=group_tmp.replace(':',',')
    display.write('%s: %s\n' % (each,group_tmp))
display.write('\ncompare:\n')
compare_name=[]
for i,each in enumerate(compares):
    compare_tmp=each
    compare_tmp=compare_tmp.split(':')
    assert len(compare_tmp)==2
    compare_name=groupnames[int(compare_tmp[0])-1]+'vs'+groupnames[int(compare_tmp[1])-1]
    display.write('%s: \t' % (each))
    display.write('%s\n' % (compare_name))
display.write('\nvenn:\n')
for i,each in enumerate(venns):
    display.write('%s: \t' % (each))
    venn_tmp=each.split('_')
    for each2 in venn_tmp:
        venn_tmp2=each2.split(':')
        venn_name=groupnames[int(venn_tmp2[0])-1]+'vs'+groupnames[int(venn_tmp2[1])-1]
        display.write('%s,' % (venn_name))
    display.write('\n')
display.write('\n')
display.write('goann: %s\n' % (goann))
if species:
    display.write('KEGG species: %s\n' % (species) )
display.write('PPI number: %s\n' % (ppi_number) )
display.write('PPI blast: %s\n' % (ppi_blast))
if argv['genenamefile']:
    display.write('genenamefile: %s\n' % (genenamefile) )
else:
    display.write('genenamefile: %s\n' % (root_dir+'/Blast_TR/Blast_Swissprot/diffgene_union.genenames'))
display.close()
######################### display all parameters END #############################
###################### define and create directory BEGIN #########################
logdir = root_dir + '/log'
qcdir = root_dir + '/QC_TR'
qcreportdir = root_dir+'/QC_TR/QCreport'
diffdir = root_dir + '/Diff_TR'
curvedir = root_dir + '/Curve_TR'
godir = root_dir + '/GOSeq_TR'
blastdir = root_dir + '/Blast_TR'
kobasdir = root_dir + '/KOBAS_TR'
ppidir = root_dir + '/PPI_TR'
candir = root_dir + '/CAN_TR'
snpdir = root_dir + '/SNP_TR'
dexseqdir = root_dir + '/DEXseq_TR'
def create_dir(directory):
    if not os.path.exists(directory):
        os.system('mkdir -p %s' % (directory))
    else:
        exit("%s already exists!" % directory)
####################### define and create directory END ##########################
#################### import scripts from config.ini BEGIN ########################
config = ConfigParser.ConfigParser()
config.read("/BJPROJ/RNA/lilin/WORK/RefTR_sjm/Moudles/config.ini")
## for QC and QCreport ##
gtf2bed = config.get("qc", "gtf2bed")
AllRunQC = config.get("qc", "allrunQC")
QCReport = config.get("qc", "QCreport")
## for CAN ##
runCAN = config.get("can","runCAN")
novelhmm = config.get("can","novelhmm")
## for SNP ##
runGATK = config.get("snp","runGATK")
## for DEXSeq ##
python = config.get("dexseq","python")
runDEXSeq = config.get("dexseq","runDEXSeq")
## for diff ##
runDiff = config.get("diff","runDiff")
runCurve = config.get("diff","runCurve")
## for enrichment ##
diffsum = config.get("enrich","diffsum")
goseq_graph = config.get("enrich","goseq_graph")
changeGO_up_down = config.get("enrich","changeGO_up_down")
R = config.get("enrich","R")
goBar = config.get("enrich","goBar")
goBar2 = config.get("enrich","goBar2")
blastx = config.get("enrich","blastx")
extractIDsEVxml = config.get("enrich","extractIDsEVxml")
getdiffGN = config.get("enrich","getdiffGN")
uniprot_sprot = config.get("enrich","uniprot_sprot")
get_my_PPI = config.get("enrich","get_my_PPI")
BLASTX_TO_PPI = config.get("enrich","BLASTX_TO_PPI")
auto_annotate = config.get("enrich","auto_annotate")
convert2kobas = config.get("enrich","convert2kobas")
KEGG_step1_blast = config.get("enrich","KEGG_step1_blast")
pathway_annotation = config.get("enrich","pathway_annotation")
runKEGG_enrich = config.get("enrich","runKEGG_enrich")
KEGG_step2_enrich = config.get("enrich","KEGG_step2_enrich")
ppi_db = config.get("enrich","ppi_db")
## for results, release and byebye ##
resultReport = config.get("report","resultReport")
dataRelease = config.get("report","dataRelease")
byeBye = config.get("report","byeBye")
##################### import scripts from config.ini END #########################
######################### generate a config file BEGIN ###########################
def create_config(project,config_file):
    cf = ConfigParser.ConfigParser()
    cf.add_section('basic')
    cf.set('basic','project',project)
    cf.set('basic','project_type','TR')
    cf.add_section('para')
    cf.set('para','fq',fq)
    cf.set('para','sample',sample)
    cf.set('para','groupname',groupname)
    cf.set('para','compare',compare)
    cf.set('para','root_dir',root_dir)
    cf.set('para','ss',ss)
    cf.set('para','venn_cluster_name',venn_cluster_name)
    cf.set('para','flag_uniform','True')
    cf.set('para','includes',','.join([str(i) for i in includes]))
    cf.set('para','fa',fa)
    cf.set('para','gtf',gtf)
    cf.set('para','goann',goann)
    cf.set('para','genenamefile',genenamefile)
    cf.write(open(config_file, 'w'))
########################## generate a config file END ############################
create_dir(logdir)
create_config(project,'%s/project.ini' %(root_dir))
############################# QC and QCreport BEGIN ##############################
def generate_qc():
    rundir = qcdir
    cmd = '''
cd %s
awk '{if($3=="exon") {print $0}}' %s > %s
msort -k mf1 -k nf4 %s > %s
perl %s %s %s
perl %s \\
-fq %s -se-pe pe -n %s -o %s -spe %s -R %s -G %s -bed %s \\
''' % (qcdir,gtf,root_dir+'/QC_TR/exon.gtf',root_dir+'/QC_TR/exon.gtf',root_dir+'/QC_TR/sorted.gtf',gtf2bed,root_dir+'/QC_TR/sorted.gtf',root_dir+'/QC_TR/sorted.bed',AllRunQC,fq,sample,qcdir,ss,fa,gtf,root_dir+'/QC_TR/sorted.bed')
    if argv['mapfile']:
        cmd += "-mapfile %s" %(mapfile)
    if argv['index']:
        cmd += " -m_ad y -index %s"  % (index)
    elif argv['ad']:
        cmd += " -ad %s" % (ad)
    return cmd,rundir
def qc(sample):
    rundir = qcdir+'/'+sample
    cmd = 'sh '+rundir+'/'+sample+'_QC.sh'
    create_dir(rundir)
    return cmd,rundir
def qcreport():
    rundir = qcreportdir
    create_dir(qcreportdir)
    cmd = '''
sh %s -dir %s -sample %s -title %s -results %s
''' %(QCReport,qcdir,sample,project,qcreportdir)
    return cmd,rundir
### QC
if not os.path.exists(qcdir):
    os.mkdir(qcdir)
cmd,rundir = generate_qc()
shell = rundir+'/generate_QC.sh'
open(shell,'w').write(cmd)
generate_qc_job = job('generate_qc',1,1,shell)

qc_jobs = {}
for eachsample in samples:
    script,tmpdir = qc(eachsample)
    jobname = eachsample+'_qc'
    shell = tmpdir+'/'+eachsample+'_QC.sh'
    open(shell,'w').write(script)
    eachjob = job(jobname,5,8,shell)
    qc_jobs[eachsample] = eachjob
### QC Report
script,tmpdir = qcreport()
shell = tmpdir+'/QC_report.sh'
open(shell,'w').write(script)
QCreport_job = job('qc_report',1,1,shell)
###############################################
## generate sjm description file of QC BEGIN ##
###############################################
qc_jobfile = open(logdir+'/'+project+'_QC.JOB','w')
## Job Specification Blocks
qc_jobfile.write(generate_qc_job.sjm())
for each in qc_jobs:
    qc_jobfile.write(qc_jobs[each].sjm())
qc_jobfile.write(QCreport_job.sjm())
## jobs order
for each in qc_jobs:
    qc_jobfile.write('order %s after %s\n' % (qc_jobs[each].jobname, generate_qc_job.jobname))
    qc_jobfile.write('order %s after %s\n' % (QCreport_job.jobname,qc_jobs[each].jobname))
## log_dir
qc_jobfile.write('\nlog_dir %s\n' %(logdir))
qc_jobfile.close()
#############################################
## generate sjm description file of QC END ##
#############################################
open(root_dir+'/sjm_QC.sh','w').write('/PUBLIC/software/public/System/sjm-1.2.0/bin/sjm %s \n' %(logdir+'/'+project+'_QC.JOB'))
assert not os.system('chmod +x %s' % (root_dir+'/sjm_QC.sh'))
############################## QC and QCreport END ###############################
################################ Analysis BEGIN ##################################
##### wrap the function in function
def generate_can():
    if ss == 'no':
        lib = 'fr-unstranded'
    elif ss == 'yes':
        lib == 'fr-firststrand'
    else:
        lib == 'fr-secondstrand'
    rundir = candir
    create_dir(rundir)
    bam = root_dir+'/QC_TR/bam'
    sam = root_dir+'/QC_TR/sam'
    cmd = '''
python %s -R %s -G %s -i %s -sam %s -o %s -lib %s -group %s -groupname %s
mkdir %s/CAN_TR/CAN/NovelGene
''' %(runCAN,fa,gtf,bam,sam,root_dir+'/CAN_TR/CAN',lib,group,groupname,root_dir)
    return cmd,rundir
def generate_snp():
    rundir = snpdir
    create_dir(rundir)
    bam = root_dir+'/QC_TR/bam'
    cmd = '''
python %s -R %s -t bam -i %s -o %s -b %s -n %s -gff %s
''' %(runGATK,fa,bam,root_dir+'/SNP_TR/SNP',sample,sample,gtf)
    return cmd,rundir
def generate_dexseq():
    rundir = dexseqdir
    create_dir(rundir)
    sam = root_dir+'/QC_TR/sam'
    allsamples_gtf = root_dir+'/CAN_TR/CAN/merged_gtf_tmap/allsamples.gtf'
    allsamples_tmap = root_dir+'/CAN_TR/CAN/merged_gtf_tmap/allsamples.allsamples.gtf.tmap'
    DEXseq_groups = group.split(',')
    for i,each in enumerate(DEXseq_groups):
        DEXseq_groups[i] = DEXseq_groups[i].split(':')
    DEXseq_compare = []
    DEXseq_compares = compare.split(',')
    for each in DEXseq_compares:
        temp = each.split(':')
        tmp1 = DEXseq_groups[int(temp[0])-1]
        tmp2 = DEXseq_groups[int(temp[1])-1]
        tmp = tmp1 + tmp2
        if (len(DEXseq_groups[int(temp[0])-1])>1 and len(DEXseq_groups[int(temp[1])-1])>1 and len(list(set(tmp)))==len(tmp1)+len(tmp2)):
            DEXseq_compare.append(each)
    DEXseq_compare = ','.join(DEXseq_compare)
    cmd = '''
python %s -G %s -i %s -g %s -n %s -c %s -o %s -p yes -s %s -a 10 -m %s
''' %(runDEXSeq,allsamples_gtf,sam,group,groupname,DEXseq_compare,root_dir+'/DEXseq_TR/DEXseq',ss,allsamples_tmap)
    return cmd,rundir
def generate_diff(samples):
    rundir = diffdir
    readcount = []
    for eachsample in samples:
        temp='%s/Diff_TR/readcount/%s.readcount' % (root_dir,eachsample)
        readcount.append(temp)
    readcount=','.join(readcount)
    create_dir(diffdir)
    sam = qcdir+'/sam'
    cmd = '''
cd %s
perl %s \\
-fa %s -sam %s -g %s -o %s \\
-group %s -groupname %s -compare %s -venn %s \\
-spe %s  -i %s
''' % (diffdir,runDiff,fa,sam,gtf,diffdir,group,groupname,compare,venn,ss,readcount)
    if argv['genenamefile']: cmd += ' -genenamefile %s ' % (genenamefile)
    return cmd,rundir
def diff():
    rundir = diffdir
    cmd = "sh %s/runDiff_analysis.sh" % (diffdir)
    return cmd,diffdir
def generate_curve():
    rundir = curvedir
    create_dir(curvedir)
    create_dir(curvedir+'/SaturationCurve')
    create_dir(curvedir+'/density')
    bam = qcdir+'/bam'
    cmd = '''
python %s  -bam %s -sample %s -n %s -r %s -fa %s -gtf %s -o %s
#sh %s
sh %s
''' % (runCurve,bam,sample,number,length,fa,gtf,curvedir,curvedir+'/generate_Curve.sh',curvedir+'/SaturationCurve/prepare_gtf_bed.sh')
    return cmd,rundir
def run_saturation(sample):
    saturationCurve = curvedir + '/SaturationCurve'
    rundir = saturationCurve
    cmd = '''
cd %s
sh %s/%s/%s_runSaturation.sh
''' % (saturationCurve,saturationCurve,sample,sample)
    return cmd,rundir

def run_density(sample):
    density = curvedir + '/density'
    rundir = density
    cmd = '''
cd %s
sh  %s/%s/%s.runDensity.sh
''' % (density,density,sample,sample)
    return cmd,rundir
def diffsum():
    rundir = diffdir
    cmd = 'perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/diffsum.pl -diffdir %s ' % (diffdir)
    return cmd,rundir
def go_all(compare, subgroup, id):
    temp=compare.split(':')
    dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
    rundir = godir + '/'+ subgroup + '/' + dir
    out=rundir
    create_dir(out)
    length=diffdir+'/Diff/genelength'
    result=out+'/'+dir+'.GO_enrichment_result.xls'
    result_up_dpwn=out+'/'+dir+'.GO_enrichment_result_up_down.xls'
    cmd = '''
###############%s#####################
perl %s -i %s -goann %s -n %s -o %s -length %s
perl %s %s %s %s
''' % (dir, goseq_graph, id, goann, dir, out, length, changeGO_up_down,out+'/'+dir+'.GO_enrichment_result.xls',diffdir+'/Diff/'+dir+'/'+dir+'.diffgene.xls',out+'/'+dir+'.GO_enrichment_result_up_down.xls')
    cmd+= '''
%s %s %s %s %s
%s %s %s %s %s
''' % (R, goBar, result, out, dir, R, goBar2, result_up_dpwn, out, dir)

    return cmd,rundir
def go(compare, subgroup, id):
    temp=compare.split(':')
    dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
    rundir = godir + '/'+ subgroup + '/' + dir
    out=rundir
    create_dir(out)
    length=diffdir+'/Diff/genelength'
    result=out+'/'+dir+'.GO_enrichment_result.xls'
    cmd = '''
###############%s#####################
perl %s -i %s -goann %s -n %s -o %s -length %s
''' % (dir, goseq_graph, id, goann, dir, out, length)
    cmd += '''
%s %s %s %s %s
''' % (R, goBar, result, out, dir)

    return cmd,rundir

def swissprot_blast():
    rundir = blastdir+'/Blast_Swissprot'
    create_dir(rundir)
    query=root_dir+'/Diff_TR/Diff/diffgene_union.seq'
    out=root_dir+'/Blast_TR/Blast_Swissprot/diffgene_union.seq.blastout'
    outdir1 = rundir
    cmd = '''
echo start blastx
date
%s -query %s -db %s -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 10 -out %s
echo blastx end
date
''' % (blastx, query, uniprot_sprot, out)
    cmd += '''
perl %s %s %s
''' % (extractIDsEVxml, out, outdir1+'/diffgene_union.genenames')
    compare=root_dir+'/Diff_TR/Diff/compare.txt'
    indir=root_dir+'/Diff_TR/Diff/'
    outdir2=root_dir+'/Diff_TR/Diff/DiffGeneList'
    cmd += '''
mkdir %s
perl %s %s %s %s %s
date
''' % (outdir2,getdiffGN, indir,compare,outdir1+'/diffgene_union.genenames',outdir2)
    return cmd,rundir

def kobas_blast(species):
    rundir = blastdir
    if not os.path.exists(rundir): create_dir(rundir)
    query=diffdir+'/Diff/diffgene_union.seq'
    if species == 'kaas':
        cmd = '''
perl %s  -n -s %s
python %s  %s /PUBLIC/database/Common/KEGG/kos %s
''' % (auto_annotate, query, convert2kobas, root_dir + '/Diff_TR/Diff/diffgene_union.seq.ko',root_dir + '/KOBAS_TR/koID.annotation')
    else:
        blastout=root_dir+'/Blast_TR/KOBAS_blast.xml'
        cmd = '''
perl %s  %s %s %s %s
sh %s
''' % (KEGG_step1_blast, query,species,blastout,root_dir+'/Blast_TR/KOBAS_blast.sh', root_dir+'/Blast_TR/KOBAS_blast.sh')
    return cmd,rundir
def kobas_pathway():
    rundir = kobasdir
    cmd = ' '
    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        out='%s/%s/%s' % (kobasdir, 'ALL', dir)
        result=root_dir+'/KOBAS_TR/ALL/'+dir+'/'+'add.'+dir+'.identify.xls'
        diff=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene.xls'
        cmd += '''
echo "###############%s#####################"
cd %s
python %s --table %s --diff %s
mv %s %s
''' % (dir, out, pathway_annotation, result, diff, 'add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html')
        out='%s/%s/%s' % (kobasdir, 'UP', dir)
        result = root_dir+'/KOBAS_TR/UP/'+dir+'/'+'add.'+dir+'.identify.xls'
        diff = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene_up.xls'
        cmd += '''
echo "###############%s#####################"
cd %s
python %s --table %s --diff %s
mv %s %s
''' % (dir, out, pathway_annotation, result, diff, 'add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html')
        out='%s/%s/%s' % (kobasdir, 'DOWN', dir)
        result = root_dir+'/KOBAS_TR/DOWN/'+dir+'/'+'add.'+dir+'.identify.xls'
        diff = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene_down.xls'
        cmd += '''
echo "###############%s#####################"
cd %s
python %s --table %s --diff %s
mv %s %s
''' % (dir, out, pathway_annotation, result, diff, 'add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html')
    return cmd, rundir


def kobas(species,compare,subgroup,id=None):
    temp=compare.split(':')
    dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
    out = '%s/%s/%s' % (kobasdir, subgroup, dir)
    rundir = out
    script = rundir+'/run.sh'
    blastout = root_dir+'/Blast_TR/KOBAS_blast.xml'
    create_dir(out)
    if species != 'kaas':
        cmd = '''
echo ###############%s#####################
perl %s -id %s -out-dir %s -species %s -blast-result %s -sample-names %s>%s
sh %s
''' %(dir, KEGG_step2_enrich, id,out,species,blastout,dir,script, script)
    else:
        cmd = '''
perl %s -diff %s -ko %s -g %s
''' % (runKEGG_enrich, id, root_dir + '/KOBAS_TR/koID.annotation',dir)
    return cmd, rundir

def ppi(ppi_blast,compare,subgroup,id=None):
    temp=compare.split(':')
    dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
    seq=root_dir+'/Diff_TR/Diff/Diff_Gene_Seq/'+dir+'.diffgene.seq'
    out='%s/%s/%s' % (ppidir, subgroup, dir)
    rundir = out
    create_dir(rundir)
    if ppi_blast == 'n':
        cmd = '''
python %s -p %s -g %s -o %s
''' % (get_my_PPI, ppi_dir+'/PPI_'+ppi_number+'.txt',id,out+'/'+dir+'.ppi.txt')
    elif ppi_blast == 'y':
        cmd = '''
python %s  --species %s --fa %s --outdir %s --name %s
''' % (BLASTX_TO_PPI, ppi_number,seq,out,dir)
    return cmd,rundir

def generate_report():
    rundir = resultdir
    cmd = 'python %s' % resultReport
    return cmd,rundir
def data_release():
    dr = root_dir+'/data_release.sh'
    cmd = 'python %s' % dataRelease
    if not os.path.exists(dr):
        open(dr,'w').write(cmd)
    else:
        exit("#data_release.py already exists!#")
def bye_bye():
    bye = root_dir+'/byebye.sh'
    cmd = 'python %s' % byeBye
    if not os.path.exists(bye):
        open(bye,'w').write(cmd)
    else:
        exit("##byebye.sh already exits!##")

##### main pipline
if set([1]).issubset(includes):
    cmd,rundir = generate_can()
    jobname = 'generate_CAN'
    shell = candir+'/generate_CAN.sh'
    open(shell,'w').write(cmd)
    generate_can_job = LocalJob(jobname,shell)
    runCufflinks_job = []
    for each in samples:
        jobname = 'runCufflinks_'+each
        shell = candir+'/runCufflinks_'+each+'.sh'
        runCufflinks_job.append(job(jobname,5,2,shell))
    jobname = 'runCuffmerge_Cuffcompare'
    shell = candir+'/runCuffmerge_Cuffcompare.sh'
    runCuffmerge_Cuffcompare_job = job(jobname,5,2,shell)
if set([1,2]).issubset(includes):
    cmd,rundir = generate_snp()
    jobname = 'generate_SNP'
    shell = snpdir+'/generate_SNP.sh'
    open(shell,'w').write(cmd)
    generate_snp_job = job(jobname,1,1,shell)
    jobname = 'workflow1'
    shell = snpdir+'/workflow1.sh'
    wokflow1_job = job(jobname,10,6,shell)
    runworkflow2_job = []
    for each in samples:
        jobname = 'workflow2_'+each
        shell = snpdir+'/workflow2_'+each+'.sh'
        runworkflow2_job.append(job(jobname,10,6,shell))
    jobname = 'workflow3'
    shell = snpdir+'/workflow3.sh'
    workflow3_job = job(jobname,10,6,shell)
if set([1,3]).issubset(includes):
    cmd,rundir = generate_dexseq()
    jobname = 'generate_DEXseq'
    shell = dexseqdir+'/generate_DEXseq.sh'
    open(shell,'w').write(cmd)
    generate_dexseq_job = job(jobname,1,1,shell)
    jobname2 = 'runDEXSeq'
    shell2 = dexseqdir+'/runDEXSeq.sh'
    runDEXSeq_job = job(jobname2,2,1,shell2)


if set([1,4,5]).issubset(includes):

    cmd,rundir = generate_diff(samples)
    jobname = 'generate_Diff'
    shell = diffdir+'/generate_Diff.sh'
    open(shell,'w').write(cmd)
    generate_diff_job = LocalJob(jobname,shell)

    cmd,rundir = diff()
    shell = rundir+'/runDiff_analysis.sh'
    open(shell,'w').write(cmd)
    jobname = 'runDiff_analysis'
    diff_job = job(jobname,5,1,shell)

    cmd,rundir = swissprot_blast()
    shell = rundir+'/runBlast_swissprot.sh'
    open(shell,'w').write(cmd)
    jobname = 'runBlast_swissprot'
    swissprot_blast_job = job(jobname,5,10,shell)

    cmd,rundir = kobas_blast(species)
    shell = rundir+'/KEGG_step1_Blast.sh'
    open(shell,'w').write(cmd)
    jobname = 'KEGG_step1_blast'
    kobas_blast_job = job(jobname,5,10,shell)



    cmd,rundir = diffsum()
    shell = rundir+'/diffsum.sh'
    open(shell,'w').write(cmd)
    diffsum_job = LocalJob('diffsum', shell)

if set([1,4,6]).issubset(includes):
    cmd,rundir = generate_curve()
    shell = rundir+'/generate_Curve.sh'
    open(shell,'w').write(cmd)
    jobname = 'generate_Curve'
    generate_Curve_job = job(jobname,1,1,shell)

    saturation_jobs = {}
    for sample in samples:
        cmd,rundir = run_saturation(sample)
        shell = rundir+'/'+sample+'.runSaturation.sh'
        open(shell,'w').write(cmd)
        jobname = sample+'_saturation'
        eachjob = job(jobname,10,1,shell)
        saturation_jobs[sample] = eachjob
    density_jobs = {}
    for sample in samples:
        cmd,rundir = run_density(sample)
        shell = rundir+'/'+sample+'.runDensity.sh'
        open(shell,'w').write(cmd)
        jobname = sample+'_density'
        eachjob = job(jobname,5,1,shell)
        density_jobs[sample] = eachjob

if set([1,4,5,7]).issubset(includes):
    go_jobs_all = {}
    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=diffdir+'/Diff/'+dir+'/'+dir+'.diffgeneID'
        cmd,rundir = go_all(compare,"ALL",id)
        shell = '%s/runGOSeq_ALL.%s.sh' % (rundir,dir)
        open(shell,'w').write(cmd)
        jobname = '%s_GO_ALL' % dir
        eachjob = job(jobname,5,1,shell)
        go_jobs_all[dir] = eachjob
    go_jobs_up = {}
    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=diffdir+'/Diff/'+dir+'/'+dir+'.diffgeneID_up'
        cmd,rundir = go(compare,"UP",id)
        shell = '%s/runGOSeq_UP.%s.sh' % (rundir,dir)
        open(shell,'w').write(cmd)
        jobname = '%s_GO_UP' % dir
        eachjob = job(jobname,5,1,shell)
        go_jobs_up[dir] = eachjob
    go_jobs_down = {}
    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=diffdir+'/Diff/'+dir+'/'+dir+'.diffgeneID_down'
        cmd,rundir = go(compare,"DOWN",id)
        shell = '%s/runGOSeq_DOWN.%s.sh' % (rundir,dir)
        open(shell,'w').write(cmd)
        jobname = '%s_GO_DOWN' % dir
        eachjob = job(jobname,5,1,shell)
        go_jobs_down[dir] = eachjob

if set([1,4,5,8]).issubset(includes):
    kobas_jobs_all = {}
    kobas_jobs_up = {}
    kobas_jobs_down = {}
    kobas_path_job = ''

    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id_a=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
        cmd,rundir = kobas(species,compare,"ALL",id_a)
        shell = "%s/runKOBAS_ALL.%s.sh" % (rundir, dir)
        open(shell,'w').write(cmd)
        kobas_jobs_all[dir] = job("%s_KOBAS_ALL" % (dir), 5, 1, shell)

        id_u = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID_up'
        cmd,rundir = kobas(species,compare,"UP",id_u)
        shell = "%s/runKOBAS_UP.%s.sh" % (rundir, dir)
        open(shell,'w').write(cmd)
        kobas_jobs_up[dir] = job("%s_KOBAS_UP" % (dir), 5, 1, shell)

        id_d = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID_down'
        cmd,rundir = kobas(species,compare,"DOWN",id_d)
        shell = "%s/runKOBAS_DOWN.%s.sh" % (rundir, dir)
        open(shell,'w').write(cmd)
        kobas_jobs_down[dir] = job("%s_KOBAS_DOWN" % (dir), 5, 1, shell)

        cmd, rundir = kobas_pathway()
        shell = '%s/run_pathway.sh' % (rundir)
        open(shell,'w').write(cmd)
        kobas_path_job = LocalJob('kobas_pathway',shell)

if set([1,4,5,9]).issubset(includes):
    ppi_jobs_all = {}
    ppi_jobs_up = {}
    ppi_jobs_down = {}

    for compare in compares:
        temp=compare.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]

        id_a = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
        cmd,rundir = ppi(ppi_blast,compare,"ALL",id_a)
        shell = "%s/runPPI_ALL.%s.sh" % (rundir,dir)
        open(shell,'w').write(cmd)
        ppi_jobs_all[dir] = job('%s_PPI_ALL' % (dir),5,7,shell)

        id_u = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID_up'
        cmd,rundir = ppi(ppi_blast,compare,"UP",id_u)
        shell = '%s/runPPI_UP.%s.sh' % (rundir,dir)
        open(shell,'w').write(cmd)
        ppi_jobs_up[dir] = job('%s_PPI_UP' % (dir),5,7,shell)

        id_d = root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID_down'
        cmd,rundir = ppi(ppi_blast,compare,'DOWN',id_d)
        shell = '%s/runPPI_DOWN.%s.sh' % (rundir,dir)
        open(shell,'w').write(cmd)
        ppi_jobs_down[dir] = job('%s_PPI_DOWN' % (dir),5,7,shell)




##### analysis jobs description file
## Job Specification Blocks
analysis_jobfile = open(logdir+'/'+project+'_analysis.JOB','w')
if set([1]).issubset(includes):
    analysis_jobfile.write(generate_can_job.sjm())
    for each in runCufflinks_job:
        analysis_jobfile.write(each.sjm())
    analysis_jobfile.write(runCuffmerge_Cuffcompare_job.sjm())
if set([1,2]).issubset(includes):
    analysis_jobfile.write(generate_snp_job.sjm())
    analysis_jobfile.write(wokflow1_job.sjm())
    for each in runworkflow2_job:
        analysis_jobfile.write(each.sjm())
    analysis_jobfile.write(workflow3_job.sjm())
if set([1,3]).issubset(includes):
    analysis_jobfile.write(generate_dexseq_job.sjm())
    analysis_jobfile.write(runDEXSeq_job.sjm()) ###

if set([1,4,5]).issubset(includes):
    analysis_jobfile.write(generate_diff_job.sjm())
    analysis_jobfile.write(diff_job.sjm())
if set([1,4,6]).issubset(includes):
    analysis_jobfile.write(generate_Curve_job.sjm())
    for job in saturation_jobs:
        analysis_jobfile.write(saturation_jobs[job].sjm())
    for job in density_jobs:
        analysis_jobfile.write(density_jobs[job].sjm())

    analysis_jobfile.write(diffsum_job.sjm())
    analysis_jobfile.write(swissprot_blast_job.sjm())
    analysis_jobfile.write(kobas_blast_job.sjm())

if set([1,4,5,7]).issubset(includes):
    for job in go_jobs_all:
        analysis_jobfile.write(go_jobs_all[job].sjm())
    for job in go_jobs_up:
        analysis_jobfile.write(go_jobs_up[job].sjm())
    for job in go_jobs_down:
        analysis_jobfile.write(go_jobs_down[job].sjm())
if set([1,4,5,8]).issubset(includes):
    for job in kobas_jobs_all:
        analysis_jobfile.write(kobas_jobs_all[job].sjm())
    for job in kobas_jobs_up:
        analysis_jobfile.write(kobas_jobs_up[job].sjm())
    for job in kobas_jobs_down:
        analysis_jobfile.write(kobas_jobs_down[job].sjm())
    analysis_jobfile.write(kobas_path_job.sjm())
if set([1,4,5,9]).issubset(includes):
    for job in ppi_jobs_all:
        analysis_jobfile.write(ppi_jobs_all[job].sjm())
    for job in ppi_jobs_up:
        analysis_jobfile.write(ppi_jobs_up[job].sjm())
    for job in ppi_jobs_down:
        analysis_jobfile.write(ppi_jobs_down[job].sjm())
analysis_jobfile.write(result_report_job.sjm())

## jobs order
if set([1]).issubset(includes):
    for each in runCufflinks_job:
        analysis_jobfile.write("order %s after %s\n" %(each.jobname,generate_can_job.jobname))
        analysis_jobfile.write("order %s after %s\n" %(runCuffmerge_Cuffcompare_job.jobname,each.jobname))
if set([1,2]).issubset(includes):
    analysis_jobfile.write("order %s after %s\n" %(wokflow1_job.jobname,generate_snp_job.jobname))
    for each in runworkflow2_job:
        analysis_jobfile.write("order %s after %s\n" %(each.jobname,wokflow1_job.jobname))
        analysis_jobfile.write("order %s after %s\n" %(workflow3_job.jobname,each.jobname))
if set([1,3]).issubset(includes):
    analysis_jobfile.write("order %s after %s\n" %(generate_dexseq_job.jobname,runCuffmerge_Cuffcompare_job.jobname))
    analysis_jobfile.write("order %s after %s\n" %(runDEXSeq_job.jobname,generate_dexseq_job.jobname))
if set([1,4,5]).issubset(includes):
    #analysis_jobfile.write("order %s after %s\n" %(generate_diff_job.jobname,runCuffmerge_Cuffcompare_job.jobname))
    analysis_jobfile.write("order %s after %s\n" % (diff_job.jobname, generate_diff_job.jobname))
if set([1,4,6]).issubset(includes):
    for job in saturation_jobs:
        analysis_jobfile.write("order %s after %s\n" % (saturation_jobs[job].jobname, generate_Curve_job.jobname))
    for job in density_jobs:
        analysis_jobfile.write("order %s after %s\n" % (density_jobs[job].jobname, generate_Curve_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (diffsum_job.jobname, diff_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (swissprot_blast_job.jobname, diff_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (kobas_blast_job.jobname, diff_job.jobname))
if set([1,4,5,7]).issubset(includes):
    for job in go_jobs_all:
        analysis_jobfile.write('order %s after %s\n' % (go_jobs_all[job].jobname, diff_job.jobname))
    for job in go_jobs_up:
        analysis_jobfile.write('order %s after %s\n' % (go_jobs_up[job].jobname, diff_job.jobname))
    for job in go_jobs_down:
        analysis_jobfile.write('order %s after %s\n' % (go_jobs_down[job].jobname, diff_job.jobname))
if set([1,4,5,8]).issubset(includes):
    for job in kobas_jobs_all:
        analysis_jobfile.write('order %s after %s\n' % (kobas_jobs_all[job].jobname, kobas_blast_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (kobas_path_job.jobname, kobas_jobs_all[job].jobname))
    for job in kobas_jobs_up:
        analysis_jobfile.write('order %s after %s\n' % (kobas_jobs_up[job].jobname, kobas_blast_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (kobas_path_job.jobname, kobas_jobs_up[job].jobname))
    for job in kobas_jobs_down:
        analysis_jobfile.write('order %s after %s\n' % (kobas_jobs_down[job].jobname, kobas_blast_job.jobname))
        analysis_jobfile.write('order %s after %s\n' % (kobas_path_job.jobname, kobas_jobs_down[job].jobname))
if set([1,4,5,9]).issubset(includes):
    for job in ppi_jobs_all:
        analysis_jobfile.write('order %s after %s\n' % (ppi_jobs_all[job].jobname, diff_job.jobname))
    for job in ppi_jobs_up:
        analysis_jobfile.write('order %s after %s\n' % (ppi_jobs_up[job].jobname, diff_job.jobname))
    for job in ppi_jobs_down:
        analysis_jobfile.write('order %s after %s\n' % (ppi_jobs_down[job].jobname, diff_job.jobname))
#analysis_jobfile.write('order %s after %s\n' % ())
if set([1,4,5,7,8,9]).issubset(includes):
    for job in go_jobs_all:
        analysis_jobfile.write('order %s after %s\n' % (result_report_job.jobname, go_jobs_all[job].jobname))
    analysis_jobfile.write('order %s after %s\n' % (result_report_job.jobname, kobas_path_job.jobname))
    for job in ppi_jobs_all:
        analysis_jobfile.write('order %s after %s\n' % (result_report_job.jobname, ppi_jobs_all[job].jobname))
else:
    for job in saturation_jobs:
        analysis_jobfile.write('order %s after %s\n' % (result_report_job.jobname, saturation_jobs[job].jobname))
    for job in density_jobs:
        analysis_jobfile.write('order %s after %s\n' % (result_report_job.jobname, saturation_jobs[job].jobname))

## logdir
analysis_jobfile.write('log_dir %s\n' %(logdir))
analysis_jobfile.close()
open(root_dir+'/sjm_Analysis.sh','w').write('/PUBLIC/software/public/System/sjm-1.2.0/bin/sjm %s \n' %(logdir+'/'+project+'_analysis.JOB'))
assert not os.system('chmod +x %s' % (root_dir+'/sjm_Analysis.sh'))
################################# Analysis END ###################################
######################### Report and data release BEGIN ##########################

########################## Report and data release END ###########################

