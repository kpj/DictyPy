import os, os.path, subprocess, time
import csv, json
import xml.etree.ElementTree as ET

from progressbar import ProgressBar

from fasta_parser import FastaParser


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

    def _handle_record(self, record, blast_result):
        #bstr = ET.tostring(blast_result, encoding='utf8', method='xml').decode('unicode_escape')
        blast_result = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')

        hids = []
        for hit in blast_result.findall('Hit'):
            hsp = hit.find('Hit_hsps').find('Hsp')

            tmp = {
                'id': hit.find('Hit_id').text,
                'e_value': float(hsp.find('Hsp_evalue').text),
                'identity': int(hsp.find('Hsp_identity').text),
                'align_len': int(hsp.find('Hsp_align-len').text)
            }
            hids.append((record.id, tmp))

        return hids

    def _finalize(self, data):
        with open('results/blast_result.csv', 'w', newline='') as fd:
            cwriter = csv.writer(fd)
            for d_id, e in data:
                ident = '%s/%s %f%%' % (e['identity'], e['align_len'], (100. * (e['identity'] / e['align_len'])))
                cwriter.writerow([d_id, e['id'], e['e_value'], ident])


def blast(data_file, Blaster):
    farser = FastaParser(data_file)
    genes = farser.parse()

    blaster = Blaster(genes)
    blaster.process()


if __name__ == '__main__':
    #blast('dicty_primary_cds', GeneBlaster)
    blast('dicty_primary_protein', ViralBlaster)
