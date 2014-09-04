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
        return tx
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
candir = root_dir + 'CAN_TR'
snpdir = root_dir + 'SNP_TR'
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
runCNA = config.get("can","runCAN")
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
perl %s %s %S
perl %s \\
-fq %s -se-pe pe -n %s -o %s -spe %s -R %s -G %s -bed %s -mapfile %s \\
''' % (git,root_dir+'/QC_TR/exon.gtf',root_dir+'/QC_TR/exon.gtf',root_dir+'/QC_TR/sorted.gtf',gtf2bed,root_dir+'/QC_TR/sorted.gtf',root_dir+'/QC_TR/sorted.bed',qcdir,AllRunQC,fq,sample,qcdir,ss,fa,gtf,root_dir+'/QC_TR/sorted.bed', mapfile)
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
    create_dir(tmpdir)
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

################################# Analysis END ###################################
######################### Report and data release BEGIN ##########################

########################## Report and data release END ###########################