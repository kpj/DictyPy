import os, os.path, subprocess, time
import json
import xml.etree.ElementTree as ET

from progressbar import ProgressBar

from fasta_parser import FastaParser


class BaseBlaster(object):
    """ Blast 'em up
    """
    BLAST_PATH = None
    DB_PATH = None
    TMP_DIR = 'blast_tmp'

    def __init__(self, genes):
        self.genes = genes

    def _handle_record(self, record, blast_result):
        """ Process record and blast result in some way and return what is supposed to be stored for this entry
        """
        raise NotImplementedError('Record handling not implemented')

    def _finalize(self):
        pass

    def process(self):
        if not os.path.isdir(BaseBlaster.TMP_DIR):
            os.mkdir(BaseBlaster.TMP_DIR)

        pbar = ProgressBar(maxval=len(self.genes))
        pbar.start()

        data = []
        for i, rec in enumerate(self.genes):
            blast_res = self._blast(str(rec.seq))
            res = self._handle_record(rec, blast_res)
            data.append(res)

            pbar.update(i)
        pbar.finish()

        with open('results/blast_result.json', 'w') as fd:
            json.dump(data, fd)

        self._finalize()

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
    DB_PATH = '/home/kpj/university/Semester06/ISC/custom_blast_db'

    def _handle_record(self, record, blast_result):
        with open(os.path.join(ViralBlaster.TMP_DIR, '%s.xml' % record.id), 'w') as fd:
            bstr = ET.tostring(blast_result, encoding='utf8', method='xml').decode('unicode_escape')
            fd.write(bstr)

        blast_result = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')

        hids = []
        for hit in blast_result.findall('Hit'):
            hsp = hit.find('Hit_hsps').find('Hsp')

            tmp = {
                'id': hit.find('Hit_id').text,
                'def': hit.find('Hit_def').text,
                'hit_len': hit.find('Hit_len').text,
                'seq_len': len(record.seq),
                'gaps': hsp.find('Hsp_gaps').text
            }
            hids.append((record.id, tmp))

        return hids

    def _finalize(self):
        # more on https://github.com/lindenb/xslt-sandbox/tree/master/stylesheets/bio/ncbi

        for fn in os.listdir(ViralBlaster.TMP_DIR):
            name = os.path.splitext(os.path.join(ViralBlaster.TMP_DIR, fn))[0]
            os.system('xsltproc --novalid data/blast2html.xsl "{name}.xml" > "{name}.xhtml"'.format(name=name))


def blast(data_file, Blaster):
    farser = FastaParser(data_file)
    genes = farser.parse()

    blaster = Blaster(genes)
    blaster.process()


if __name__ == '__main__':
    #blast('dicty_primary_cds', GeneBlaster)
    blast('dicty_primary_protein', ViralBlaster)
