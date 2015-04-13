import os, os.path, re, subprocess, time
import csv, json
import xml.etree.ElementTree as ET

from progressbar import ProgressBar

from fasta_parser import FastaParser
from utils import extract_gene_name


class DefParser(object):
    def __init__(self, blast_def):
        self.bdef = blast_def

    def get_virus_name(self):
        res = re.search(r'\[.*?\]', self.bdef)
        return res.group()[1:-2].strip() if res is not None else self.bdef

    def get_gi(self):
        return '|'.join(self.bdef.split('|')[:3])

class BaseBlaster(object):
    """ Blast 'em up
    """
    BLAST_PATH = None
    DB_PATH = None

    def __init__(self, genes):
        self.genes = genes

    def _handle_record(self, record, blast_result):
        """ Process record and blast result in some way and return what is supposed to be stored for this entry
        """
        raise NotImplementedError('Record handling not implemented')

    def _finalize(self, data):
        pass

    def process(self):
        pbar = ProgressBar(maxval=len(self.genes))
        pbar.start()

        data = []
        for i, rec in enumerate(self.genes):
            blast_res = self._blast(str(rec.seq))
            res = self._handle_record(rec, blast_res)
            data.extend(res)

            pbar.update(i)
        pbar.finish()

        self._finalize(data)

    def _blast(self, seq):
        cmd = [
            self.BLAST_PATH,
            '-db', self.DB_PATH,
            '-outfmt', '5', # xml output
            '-num_threads', '8',
            '-query', '-' # read sequence from stdin
        ]
        # -qcov_hsp_perc

        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate(input=seq.encode('utf-8'))

        return ET.fromstring(stdout)

class GeneBlaster(BaseBlaster):
    BLAST_PATH = '/home/kpj/blast/ncbi-blast-2.2.30+-src/c++/ReleaseMT/bin/blastn'
    DB_PATH = '/home/kpj/blast/db/nt.00'

    def _handle_record(self, record, blast_result):
        info = {}
        info['hits'] = []
        info['gene'] = record.id

        blast_result = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')
        for hit in blast_result.findall('Hit'):
            hit_id = hit.find('Hit_id').text
            info['hits'].append(hit_id)

        return [info]

    def _finalize(self, data):
        with open('results/blast_result.json', 'w') as fd:
            json.dump(data, fd)

class ViralBlaster(BaseBlaster):
    BLAST_PATH = '/home/kpj/blast/ncbi-blast-2.2.30+-src/c++/ReleaseMT/bin/blastp'
    DB_PATH = '/home/kpj/university/Semester06/ISC/custom_blast_db'

    def __init__(self, genes):
        super().__init__(genes)

        self.used_gis = set()

    def _handle_record(self, record, blast_result):
        #bstr = ET.tostring(blast_result, encoding='utf8', method='xml').decode('unicode_escape')
        blast_result = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')

        hids = []
        for hit in blast_result.findall('Hit'):
            hsp = hit.find('Hit_hsps').find('Hsp')
            defp = DefParser(hit.find('Hit_def').text)

            if defp.get_gi() in self.used_gis: continue
            self.used_gis.add(defp.get_gi())

            dicty = {
                'id': record.id,
                'name': extract_gene_name(record)
            }

            blast = {
                'virus_name': defp.get_virus_name(),
                'id': defp.get_gi(),
                'e_value': float(hsp.find('Hsp_evalue').text),
                'identity': int(hsp.find('Hsp_identity').text),
                'align_len': int(hsp.find('Hsp_align-len').text)
            }
            hids.append((dicty, blast))

        return hids

    def _finalize(self, data):
        with open('results/blast_result.csv', 'w', newline='') as fd:
            cwriter = csv.writer(fd)
            cwriter.writerow([
                'dicty p. id', 'dicty p. name',
                'virus', 'virus p. id',
                'e score',
                'identity (abs.)', 'identity (rel.)'
            ])

            for dicty, blast in data:
                cwriter.writerow([
                    dicty['id'], dicty['name'],
                    blast['virus_name'], blast['id'],
                    blast['e_value'],
                    '%s/%s' % (blast['identity'], blast['align_len']), 100. * (blast['identity'] / blast['align_len'])
                ])


def blast(data_file, Blaster):
    farser = FastaParser(data_file)
    genes = farser.parse()

    blaster = Blaster(genes)
    blaster.process()


if __name__ == '__main__':
    #blast('dicty_primary_cds', GeneBlaster)
    blast('dicty_primary_protein', ViralBlaster)
