#  Name: 
#  Author: rotem.tal
#  Description:
#

import sys
import requests
import gzip
import os
import mechanicalsoup
from itertools import groupby, chain
import numpy as np
from collections import Counter
from functools import reduce

class Reader:

    def __init__(self, pfam_id,alnType='seed' ,motif=True,save_file=False):
        """
        param pfam_id: pfam motif ID
        :param motif: selecting between motif or full seq
        :param save_file: boolean, whatever or not to save fasta file
        :param alnType: Alignment type
        """
        bro = mechanicalsoup.StatefulBrowser()
        self.fasta_name = "fasta%s.fasta" % pfam_id
        if motif:
            bro.open("http://pfam.xfam.org/family/%s#tabview=tab3"%(pfam_id))
            bro.select_form('form[action=\"/family/alignment/download/format\"]')
            bro['format'] = 'fasta'
            bro['alnType'] = alnType.lower()
            res = bro.submit_selected()
            with open(self.fasta_name,'wb') as f:
                f.write(res.content)
        self.req = requests.get("https://pfam.xfam.org/family/%s/alignment/long/gzipped"%(pfam_id))
        with open("%s.gz"%self.fasta_name,'wb') as d:
            d.write(self.req.content)
        with gzip.open("%s.gz"%self.fasta_name, 'rb') as f, open(self.fasta_name,'wb') as g, open("background.fasta",
                                                                                                  'wb') as c:
            c.write(f.read())
            if not motif:
                g.write(f.read())
        bg = fasta_read("PF05488_full_length_sequences.fasta")
        bg = Counter(i for i in list(chain.from_iterable(bg)))
        s = sum(bg.values())
        self.background_e = {d: np.log(bg[d]/s) for d in bg.keys()}
        os.remove("background.fasta")
        os.remove("%s.gz"%self.fasta_name)
        fast = self.fasta_read()
        self.fasta = np.array(list(map(list, fast)))
        if not save_file:
            os.remove(self.fasta_name)


    def get_fasta(self):
        return self.fasta

    def fasta_read(self, fasta_name=None):
        """
        Readsa fasta and return a dict of header: sequence
        :param fasta_name: Name of the fasta file
        :return: a dict of header: sequence
        """
        if not fasta_name:
            fasta_name = self.fasta_name
        f = open(fasta_name)
        faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
        fas_d = []
        for header in faiter:
            header = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(faiter))
            fas_d.append(seq)
        return fas_d

    def get_bg_e(self):
        return self.background_e

# r= Reader("PF00096", save_file=True, alnType='seed')
# d = r.get_fasta()
# print(r.background_e)

def fasta_read(fasta_name=None):
    """
    Readsa fasta and return a dict of header: sequence
    :param fasta_name: Name of the fasta file
    :return: a dict of header: sequence
    """
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    fas_d = []
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        fas_d.append(seq)
    return fas_d
