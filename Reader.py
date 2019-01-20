#  Name: 
#  Author: rotem.tal
#  Description:
#

import sys
import requests
import gzip
import os
import mechanicalsoup
from itertools import groupby


class Reader:

    def __init__(self, pfam_id, save_file=False, alnType='seed'):
        """

        :param pfam_id: pfam motif ID
        :param save_file: boolean, whatever or not to save fasta file
        :param alnType: Alignment type
        """
        bro = mechanicalsoup.StatefulBrowser()
        bro.open("http://pfam.xfam.org/family/%s#tabview=tab3"%(pfam_id))
        bro.select_form('form[action=\"/family/alignment/download/format\"]')
        bro['format']= 'fasta'
        bro['alnType'] = alnType
        res = bro.submit_selected()
        self.fasta = res.content
        self.fasta_name = "fasta%s.fasta"%pfam_id
        if save_file:
            with open(self.fasta_name,'wb') as f:
                f.write(self.fasta)

        # self.req = requests.get("https://pfam.xfam.org/family/%s/alignment/long/gzipped"%(pfam_id))
        # with open("fasta%s.gz"%(pfam_id),'wb') as d:
        #     d.write(self.req.content)
        # with gzip.open("temp.gz", 'rb') as f:
        #     self.fasta = f.read()
        #     print(self.fasta)
        # if not save_file:
        #     os.remove("fasta%s.gz"%(pfam_id))


    def get_fasta(self):
        return self.fasta

    def fasta_read(self, fasta_name=None):
        if not fasta_name:
            fasta_name = self.fasta_name
        f = open(fasta_name)
        faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
        fas_d = {}
        for header in faiter:
            header = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(faiter))
            fas_d[header] = seq
        return fas_d



