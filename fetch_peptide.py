import re
import warnings

class FindPep():
    '''
    This class is initialized with a protein sequence (string format). 
    It contains two functions, ‘find()’ and ‘extract_seq()’.
    
    The find function takes as input a peptide sequence (string) and finds the 
    starting positions of all the peptide occurrences  in the sequence
    
    The extract_seq function takes as input a position (int) and extract from
    the protein sequence (self.seq) a peptide centred at the specified 
    position plus minus an extent of n amino acids (extend). It adds - at both
    ends to fill gaps. it can test if the central amino acid is in the list
    provided by expected_aa
    '''
    def __init__(self, seq):
        #treat everithing upper case 
        if type(seq) != type('AA'):
            raise Exception('error', 'sequence must be a string') 
        seq = seq.upper()
        print len(seq)
        self.seq = seq
        
    def find(self, peptide):
        res = []
        for item in re.finditer(peptide, self.seq):
            res.append( item.start() )
        return res
    
    def extract_seq(self, position=0, extend=1, expected_aa = []):
        #######################################################
        #position numbering according to python counting from 0
        #generally speaking remove 1 from the aminoacid position
        #######################################################
        
        #some control on input
        if type(position) != type(1) or position < 0:
            raise Exception('error', 'position must be positive int') 
        if position+1 > len(self.seq):
            raise Exception('error', 'position must be inside sequence')
        if type(extend) != type(1) or extend < 0:
            raise Exception('error', 'extend must be positive int')
        if type(expected_aa) != type([]):
            raise Exception('error', 'expected_aa muts be a list of AA')
        for index, aa in enumerate(expected_aa):
            if type(aa) != type('a') or len(aa) >1:
                raise Exception('error', 'expected_aa', index,
                                'is not a chr')
            
        
        #treat everithing upper case        
        expected_aa = set([n.upper() for n in expected_aa])   
            
        if position < extend and position >= len(self.seq)-extend:
            add_front = extend-position
            add_back = extend - (len(self.seq) - position-1)
            res = add_front*'-'+self.seq+'-'*add_back
        
        elif position >= extend and position < len(self.seq)-extend:
            res = self.seq[position-extend:position+extend+1]
        
        elif position < extend:
            res = self.seq[:position+extend+1]
            add = extend-position
            res = '-'*add + res
        
        elif position >= len(self.seq)-extend:
            res = self.seq[position-extend:]
            add = extend - (len(self.seq) - position-1)
            res = res +'-'*add
            
        #test if the length of the peptide is correct
        if len(res) != extend*2+1:
            print position, extend, res,  len(self.seq)
            raise Exception('error', 'peptide length')
        
        if len(expected_aa) > 0:
            if res[extend] not in expected_aa:
                 warnings.warn("central AA not in expected")
                    
        #print res   
        return res



if __name__ == '__main__':
    test_seq = '''
    MNLKALVVIASVAVTSALPKGEEGDIIGTFNFSSSDSQPLKIHWVDTPDSSGSNLVKRSA
    HTESVCVHAGTATGADLHWLNAICTGKSTYTVNCAPAGNKNAGSTHTGTCPAGQDCFQLE
    QVGNFWGDREPDATCSPSNTVFDAVDDKEATHVNGKVVTREGKPGIGRKLIRLKAQVYRR
    DGHYGQTSRMGFFRNGKEVYHIDNVASMGPTWNFDPSSDQSFSFFFTPGPNAFRIQGTLN
    LA
    '''
    #making one long string
    test_seq = test_seq.replace(' ','').replace('\n','')
    print test_seq[0:10]
    test_class = FindPep(test_seq)
    test_pep = test_class.extract_seq(position=9, extend=10, expected_aa=['A'])
    print test_pep
    assert '-MNLKALVVIASVAVTSALPK' == test_pep 
    
    test_pep = test_class.extract_seq(position=241, extend=10, expected_aa=['A'])
    print test_pep    
    
    


    

    
