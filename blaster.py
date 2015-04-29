import os, os.path, re, subprocess, time
import csv, json
import xml.etree.ElementTree as ET
import urllib.request

from progressbar import ProgressBar
from Bio import Entrez, SeqIO

from fasta_parser import FastaParser
from utils import extract_gene_name


Entrez.email = "kpjkpjkpjkpjkpjkpj@gmail.com"

class DefParser(object):
    def __init__(self, blast_def):
        self.bdef = blast_def

    def get_virus_name(self):
        res = re.search(r'\[.*?\]', self.bdef)

        if res is None:
            return self._query_virus_name(self._get_accession(self.bdef))
        else:
            return res.group()[1:-2].strip()

    def _query_virus_name(self, accession):
        """ First ty to access information on some japanese server via the Genebank accession id (NCBI wants to send us mail if we excessively use their service). If that's not given try the GI id with NCBI.
        """

        url = 'http://getentry.ddbj.nig.ac.jp/getentry/dad/%s/?filetype=txt' % accession
        content = urllib.request.urlopen(url).read().decode('utf-8')

        res = re.search(r'SOURCE\s*(.*)', content)
        if not res is None:
            return res.group(1).strip()
        else:
            handle = Entrez.efetch(db='nucleotide', id='9626820', rettype='gb', retmode='text')
            record = SeqIO.read(handle, 'genbank')
            handle.close()

            try:
                return record.annotations['source']
            except KeyError:
                return 'unknown'

    def _get_accession(self, bdef):
        return bdef.split('|')[3].strip()

    def get_gi(self):
        return '|'.join(self.bdef.split('|')[:3])

class BaseBlaster(object):
    """ Blast 'em up
    """
    BLAST_PATH = None
    DB_PATH = None
    EXTRA_BLAST_ARGS = []

    def __init__(self, genes):
        self.genes = genes
        self.setup()

        # -qcov_hsp_perc
        self.blast_cmd = [
            self.BLAST_PATH,
            '-db', self.DB_PATH,
            '-outfmt', '5', # xml output
            '-num_threads', '8',
            '-query', '-' # read sequence from stdin
        ]
        self.blast_cmd += self.EXTRA_BLAST_ARGS

        print('Blasting with\n >', ' '.join(self.blast_cmd))

    def setup(self):
        pass

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
        proc = subprocess.Popen(self.blast_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
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

    def setup(self):
        self.used_gis = set()

    def _handle_record(self, record, blast_result):
        iter_hits = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')

        hids = []
        for hit in iter_hits.findall('Hit'):
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
                'align_len': int(hsp.find('Hsp_align-len').text),
                'xml': ET.tostring(blast_result, encoding='utf8', method='xml').decode('unicode_escape')
            }
            hids.append((dicty, blast))

        return hids

    def _finalize(self, data):
        rdir = 'viral_XML_dump'
        if not os.path.isdir(rdir):
            os.mkdir(rdir)

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

                raw_name = dicty['id'].replace('|', '_')
                with open(os.path.join(rdir, '%s.xml' % raw_name), 'w') as fd:
                    fd.write(blast['xml'])
                os.system('xsltproc --novalid tools/blast2html.xsl "{name}.xml" > "{name}.xhtml"'.format(name=os.path.join(rdir, raw_name)))

class rRNABlaster(BaseBlaster):
    BLAST_PATH = '/home/kpj/blast/ncbi-blast-2.2.30+-src/c++/ReleaseMT/bin/blastn'
    DB_PATH = '/home/kpj/university/Semester06/ISC/rrna_blast_db'
    EXTRA_BLAST_ARGS = ['-task', 'blastn-short']

    """
    BLAST database was created as follows:
        $ makeblastdb \
            -dbtype nucl \
            -in "data/rRNA.fa" \
            -input_type fasta \
            -out rrna_blast_db \
            -title rrna_blast_db
    """

    def _handle_record(self, record, blast_result):
        iter_msg_box = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_message')
        if not iter_msg_box is None and iter_msg_box.text == 'No hits found': return []

        iter_hits = blast_result.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits')
        for hit in iter_hits.findall('Hit'):
            hsp = hit.find('Hit_hsps').find('Hsp')
            e = hsp.find('Hsp_evalue').text

            with open('e_values.dat', 'a') as fd:
                fd.write(e + '\n')

        bstr = ET.tostring(blast_result, encoding='utf8', method='xml').decode('unicode_escape')
        return [(record, bstr)]

    def _finalize(self, data):
        rdir = 'rRNA_XML_dump'
        if not os.path.isdir(rdir):
            os.mkdir(rdir)

        for record, xml in data:
            with open(os.path.join(rdir, '%s.xml' % record.id.replace('|', '_')), 'w') as fd:
                fd.write(xml)


def blast(data_file, Blaster):
    farser = FastaParser(data_file)
    genes = farser.parse()

    blaster = Blaster(genes)
    blaster.process()


if __name__ == '__main__':
    #blast('dicty_primary_cds', GeneBlaster)
    blast('dicty_primary_protein', ViralBlaster)
    #blast('results/regex_lookup.fa', rRNABlaster)
