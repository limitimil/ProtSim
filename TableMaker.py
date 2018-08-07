import EvalFunctions as ef
import os
import pandas as pd
from tqdm import tqdm
class case(object):
    memberfiles = [
        'complexes.csv',
        'models.csv',
    ]
    memberrepos = [
        'alignments',
        'alignments_fiberdock',
    ]
    def __init__(self, repo = 'wildtype'):
        assert(os.path.isdir(repo))
        self.repo = repo
        assert(self.checkrepo())
    def checkrepo(self):
        ret = True
        ret &= all(os.path.isfile(os.path.join(self.repo,f) )\
            for f in self.memberfiles
        )
        ret &= all(os.path.isdir(os.path.join(self.repo,f) )\
            for f in self.memberrepos
        )
        return ret
    def wildtype(self):
        return self.patient_aug()
    def patient_aug(self):
        def find2nd_(s):
            count = 1
            pos = -1
            for i, c in enumerate(s):
                if count == 0:
                    return pos
                if c == '_':
                    pos = i
                    count -= 1
                    continue
            return -1
        def getPDBID(s):
            return s.split('_')[1]
        def getChain(s):
            return s.split('.')[0].split('_')[-1]
        def func(row):
            assert(getPDBID(row['chainA']) == getPDBID(row['chainB']))
            ret = {
                'chainA': row['chainA'][:find2nd_(row['chainA'])],
                'chainB': row['chainB'][:find2nd_(row['chainB'])],
                'pdbid': getPDBID(row['chainA']),
                'chains': '/'.join([getChain(row['chainA']),
                getChain(row['chainB'])]),
                'score': row['score']
                }
            return pd.Series(ret)
        try:
            df = self.patient()
        except:
            return None
        ret = df.apply(func, 1)
        return ret

    def patient(self):
        def getChains(s):
            #this function will help you parse the file name like:
            #EGFR_HUMAN_0___ABL1_HUMAN_0___pdb1fpu_B_pdb1fpu_A.pdb
            ret = s.split('___')[:2]
            assert(len(ret) == 2)
            return ret
        dn = os.path.join(self.repo, 'alignments_fiberdock')
        #check complexes.csv
        try:
            df = pd.read_csv(os.path.join(self.repo, 'complexes.csv'))
        except pd.io.common.EmptyDataError:
            df = pd.DataFrame()
        #collect fiberdock result
        fieldnames = ['chainA','chainB','score']
        res = []
        for f in os.listdir(dn):
            chain_A, chain_B = getChains(f)
            score = ef.evaluation(os.path.join(dn, f,'resultFile.ref'))
            score = score[0][1] if score else None
            res.append([
                chain_A,
                chain_B,
                score
            ])
        return pd.DataFrame(res, columns=fieldnames)
class TableMaker(object):
    def __init__(self, wildtype, patients):
        self.wildtype = case(wildtype)
        self.patients = {os.path.basename(p):case(os.path.join(p)) for p in patients}
    def make_table(self, patientlist):
    #this function will help you arrange a pandas dataframe
    #patientlist accept patient ID code reather than the case path.
        self.errorlist = []
        ret = self.wildtype.wildtype()
        ret.columns = ['chainA','chainB','chains','pdbid','wildtype']
        for p in tqdm(patientlist, ascii = True):
            new = self.wildtype.wildtype()
            if p in self.patients.keys():
                mut = self.patients[p].patient_aug()
                if mut is None:
                    self.errorlist.append(p)
                    continue
                for index, row in mut.iterrows():
                    assert(new.loc[(new['chainA'] == row['chainA']) &
                        (new['chainB'] == row['chainB']) & (new['pdbid'] ==
                        row['pdbid']) & (new['chains'] == row['chains']) 
                        ].shape[0] == 1)
                    new.loc[new.index[(new['chainA'] == row['chainA']) &
                        (new['chainB'] == row['chainB']) & (new['pdbid'] ==
                        row['pdbid']) & (new['chains'] == row['chains'])]
                        ,'score'] = row['score']
            new.columns = ['chainA','chainB','chains','pdbid',p]
            ret = pd.concat([ret, new[p]], 1)
        return ret
