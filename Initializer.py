#this is an adjustment of script0329 to fit current need.
import os, re

import pivotHandler
import FastaGraber
from core import Stage

import pandas as pd
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO
import tqdm
try:
    FG = FastaGraber.FastaGraber(
        repo = 'fastas',
        mapping = 'G2Umapping.csv'
    )
except:
    sys.stderr.write('initialize Fasta Graber failed!\n')
    sys.stderr.write(
        're-initial FastaGraber with Initializer.FG = FastaGraber()\n')
def _inrange(range_s, pos_m): #range_s: \d-\d// pos_m: \w\d+\w
    r = map(int,range_s.split('-'))
    pos = map(lambda m: int(re.search('\d+',m).group(0)),
        pos_m)
    return any(
        p>= int(r[0]) and p <= int(r[1])\
        for p in pos
        )
class _pivot(pivotHandler.pivot):
    def majors(self):
        return self.piv.keys()
    def __getitem__(self, k):
        return self.piv[k]
class PatientInit(object):
    def __init__(self, mutantfile='mutants/TCGA.D1.A15X.01.mutant',
    notefile='template_note.csv'):
        self.note = pd.read_csv(notefile)
        self.mutant = _pivot()
        self.mutant.read(mutantfile)
    def tableGen(self):
        self.status = self.note.apply(self.func, 1)
        newnote = self.note[
            self.status['mutated_status'] == True]
        return newnote
    def mutate(self, fas, gn):
        def parser(m):
            return m[0], int(m[1:-1]), m[-1]
        if gn not in self.mutant.majors():
            return fas
        mutation = self.mutant[gn]
        tmp = fas.seq.tomutable()
        for m in mutation:
            wildtype, location, alteration = parser(m)
            assert(wildtype == fas.seq[location-1])
            tmp[location-1] = alteration
        fas.seq = tmp
        return fas
    def save_table(self, custom = None, customName='template_note.csv.seg',destination='./out'):
        assert(os.path.isdir(destination))
        if custom is not None:
            custom.to_csv(os.path.join(destination, customName), index=False)
        else:
            self.note.to_csv(os.path.join(destination, customName), index=False)
    def save_fasta(self, custom = None,destination='./out'):
        def func2( row):
#            print row
#            self.count += 1
            def make_chain(config):
                name = '{}_{}_{}'.format(config['chain'], config['PDBid'].lower(),
                config['code'])
                #for chainA
                fas = FG.grabPathByGene(config['chain'])
                fas = list(SeqIO.parse(fas,'fasta'))[0]
                fas = self.mutate(fas, config['chain'])
                loc = config['loc'].split('-')
                newfas = SeqRecord(
                    fas.seq[int(loc[0])-1: int(loc[1])],
                    id = name,
                    description= config['loc'])
                ret = os.path.join(destination, name+'.fasta')
                SeqIO.write(newfas, ret, 'fasta')
                return ret
            resA = make_chain({
                'chain': row['chainA'], 'PDBid':row['PDBid'],
                'code': row['codeA'], 'loc':row['locA']})
            resB = make_chain({
                'chain': row['chainB'], 'PDBid':row['PDBid'],
                'code': row['codeB'], 'loc':row['locB']})
            return pd.Series({'fastaA': resA, 'fastaB': resB})
        assert(os.path.isdir(destination))
        target = custom if custom is not None else self.note
        fastas = target.apply(func2, 1)
        return fastas

    def func(self, row):
        def check_mutated(config):
            gn = config['gene']
            if gn in self.mutant.majors() and \
                _inrange(config['range'], self.mutant[gn]):
                return True
            return False
        res = [check_mutated({'gene': row['chainA'], 'range':row['locA']}),
            check_mutated({'gene': row['chainB'], 'range':row['locB']})]
        return pd.Series({'mutated_status':any(res)})
def grab_mutants_under(p):
    jl = os.listdir(p)
    jl = map(lambda k: os.path.join(p, k), jl)
    return jl
class BatchInit(object):
    notefile = 'template_pair_note.csv'
    def __init__(self, joblist, outdir):
        self.joblist = {jj : PatientInit(jj, self.notefile) for jj in joblist}
        self.outdir = outdir
        assert(os.path.isdir(self.outdir))
        self.BuildRepo()
    def run(self):
        for k,v in tqdm.tqdm(self.joblist.items(),ascii=True):
            tab = v.tableGen()
            dd = re.sub('.mutant', '', k)
            dd = os.path.join(self.outdir, dd.split('/')[-1])
            fastas = v.save_fasta(custom=tab, destination=os.path.join(dd,'fastas'))
            fastas = fastas.applymap(os.path.basename)
            res = pd.concat([tab, fastas],axis = 1)
            v.save_table(custom=res,customName=self.notefile + '.seg' ,destination=dd)
    def BuildRepo(self):
        for k in self.joblist.keys():
            p = os.path.join(
                self.outdir,
                re.sub('.mutant', '', k).split('/')[-1])
            if not os.path.isdir(p):
                os.mkdir(p)
            if not os.path.isdir(os.path.join(p, 'fastas')):
                os.mkdir(os.path.join(p, 'fastas'))
class WildtypeInit(Stage):
    memberfiles = []
    memberrepos = []
    outdirs = ['fastas']
    def __init__(self, notefile='template_note.csv', base = 'wildtype'):
        Stage.__init__(self, base = base)
        #prepare csvfile
        self.df = pd.read_csv(notefile)
        self.result = self.df.apply(self.func, 1)
        #save
        self.resultfile = os.path.join(self.base, notefile+'.seg')
        self.result.columns = ['fastaA','fastaB']
        ret = pd.concat([self.df, self.result],1)
        ret.to_csv(self.resultfile, index=False)
    def func(self,row):
        def make_chainA(row):
            name = '{}_{}_{}'.format(row['chainA'], row['PDBid'].lower(),
            row['codeA'])
            #for chainA
            fas = FG.grabPathByGene(row['chainA'])
            fas = list(SeqIO.parse(fas,'fasta'))[0]
            loc = row['locA'].split('-')
            newfas = SeqRecord(
                fas.seq[int(loc[0])-1: int(loc[1])],
                id = name,
                description= row['locA'])
            ret = os.path.join(self.base, 'fastas', name+'.fasta')
            SeqIO.write(newfas, ret, 'fasta')
            return ret
        def make_chainB(row):
            name = '{}_{}_{}'.format(row['chainB'], row['PDBid'].lower(),
            row['codeB'])
            #for chainB
            fas = FG.grabPathByGene(row['chainB'])
            fas = list(SeqIO.parse(fas,'fasta'))[0]
            loc = row['locB'].split('-')
            newfas = SeqRecord(
                fas.seq[int(loc[0])-1: int(loc[1])],
                id = name,
                description= row['locB'])
            ret = os.path.join(self.base,'fastas',name+'.fasta')
            SeqIO.write(newfas, ret, 'fasta')
            return ret
        return pd.Series([
            os.path.basename(make_chainA(row)),
            os.path.basename(make_chainB(row)) ])
