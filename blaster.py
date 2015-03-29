import collections
import os.path, subprocess, time
import json
import xml.etree.ElementTree as ET

from fasta_parser import FastaParser


class BaseBlaster(object):
    """ Blast 'em up
    """
    BLAST_PATH = None
    DB_PATH = None

    def __init__(self, genes):
        self.genes = genes
        self.timer_cache = collections.deque(maxlen=1000)

    def _handle_record(self, record, blast_result):
        raise NotImplementedError('Record handling not implemented')

    def process(self):
        print('Blast!')
        data = []
        for i, rec in enumerate(self.genes):
            start = time.time()

            blast_res = self._blast(str(rec.seq))
            res = self._handle_record(rec, blast_res)
            data.append(res)

            self.timer_cache.append(time.time() - start)
            self._handle_timer(start, i)
        print()

        with open('results/blast_result.json', 'w') as fd:
            json.dump(data, fd)

    def _handle_timer(self, start, cur_index):
        avg_dur = sum(self.timer_cache) / len(self.timer_cache)
        remaining_entries = len(self.genes) - cur_index

        remaining_seconds = avg_dur * remaining_entries

        m, s = divmod(remaining_seconds, 60)
        h, m = divmod(m, 60)

        print('\r>> %d:%02d:%02d remaining' % (h, m, s), end='')

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

        with open('foo', 'w') as fd:
            fd.write(stdout.decode('utf-8'))

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

        return info

class ViralBlaster(BaseBlaster):
    BLAST_PATH = '/home/kpj/blast/ncbi-blast-2.2.30+-src/c++/ReleaseMT/bin/blastp'
    DB_PATH = '/home/kpj/blast/db/nr/nr.00'

    def _handle_record(self, record, blast_result):
        return []


def blast(data_file, Blaster):
    farser = FastaParser(data_file)
    genes = farser.parse()

    blaster = Blaster(genes)
    blaster.process()


if __name__ == '__main__':
    #blast('dicty_primary_cds', GeneBlaster)
    blast('dicty_primary_protein', ViralBlaster)
