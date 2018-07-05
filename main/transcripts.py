from __future__ import division
import gzip
from operator import itemgetter
import pysam
import datetime
import os


allowed_chroms = map(str, range(1, 24)) + ['X', 'Y', 'MT']

class Transcript(object):

    def __init__(self, line):

        #self.cava_db_line = line

        cols = line.split('\t')

        self.id = cols[0]
        self.gene_symbol = cols[1]
        self.chrom = cols[4]
        self.start = int(cols[6])
        self.end = int(cols[7])
        self.strand = '+' if cols[5] == '1' else '-'
        self.coding_start = int(cols[9]) - 1
        self.coding_end = int(cols[10]) - 1
        self.cdna_coding_start = int(cols[8])
        self.info = cols[3]

        self.exons = []
        for i in range(1, len(cols) - 11, 2):
            self.exons.append(Exon('{}-{}'.format(cols[10 + i], cols[11 + i])))


    def cds_regions(self):
        """Return list of CDS regions

            Notes:
                Returned is a list of (x,y) pairs where coordinates x and y are 0-based, inclusive for x and exclusive for y
        """

        if self._any_unset(['exons', 'coding_start', 'coding_end']):
            return

        ret = []
        for exon in self.exons:
            cds_region = exon.get_cds(self.coding_start, self.coding_end)
            if cds_region is not None:
                ret.append(cds_region)
        return ret


    def get_cds_length(self):
        """Return CDS length"""

        cds_regions = self.cds_regions()
        if cds_regions is None:
            return
        ret = 0
        for (cds_start, cds_end) in cds_regions:
            ret += cds_end - cds_start
        return ret


    def _any_unset(self, fields):
        """Check if any of the required fields is None"""

        return any([getattr(self, f) is None for f in fields])


    def __str__(self):

        return '{} > {} > {}:{}-{} ({}), {}-{}'.format(
            self.gene_symbol,
            self.id,
            self.chrom,
            self.start,
            self.end,
            self.strand,
            self.coding_start,
            self.coding_end
        )


class Exon(object):
    """Class for a single exon"""

    def __init__(self, s):
        """Constructor of Exon class

            Notes:
                s is a string of the syntax "start-end", where coordinates are 0-based, inclusive for start, exclusive for end
        """

        [s, e] = s.split('-')
        self.start = int(s)
        self.end = int(e)

    def length(self):
        """Return length of exon"""

        return self.end - self.start

    def get_cds(self, coding_start, coding_end):
        """Return CDS interval or None if there is no CDS in the exon"""

        coding_min = min(coding_start, coding_end)
        coding_max = max(coding_start, coding_end)

        if self.end - 1 < coding_min or self.start > coding_max:
            return

        cds_start = max(self.start, coding_min)
        cds_end = min(self.end - 1, coding_max) + 1
        return (cds_start, cds_end)


    def __str__(self):
        """String representation"""

        return 'exon:' + str(self.start) + '-' + str(self.end)

    __repr__ = __str__


def read_ensembl_db(fn):

    ret = {}

    for line in gzip.open(fn):

        line = line.strip()
        if line == '' or line[0] == '#':
            continue

        transcript = Transcript(line)

        ret[transcript.id] = transcript

    return ret


class TranscriptDBWriter(object):
    """Class for creating new transcript database"""

    def __init__(self, fn, source='', build='', columns=[]):
        """Constructor of the TranscriptDBWriter class"""

        self._fn = fn
        self._source = source
        self._build = build
        self._columns = [x.lower() for x in columns]
        self._records = {c: [] for c in allowed_chroms}
        self.idx_chrom = self._columns.index('chrom')
        self.idx_start = self._columns.index('start')
        self.idx_end = self._columns.index('end')

    def add(self, transcript):
        """Add transcript to DB"""

        record = []
        for c in self._columns:
            if c in ['exons', 'cdna_exons']:
                record.append(','.join([str(e.start) + '-' + str(e.end) for e in getattr(transcript, c.lower())]))
            elif c in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                record.append(int(getattr(transcript, c.lower())))
            else:
                record.append(str(getattr(transcript, c.lower())))
        self._records[transcript.chrom].append(record)

    def _sort_records(self):
        """Sort records by chrom, start, end"""

        idx_start = self._columns.index('start')
        idx_end = self._columns.index('end')
        for c in allowed_chroms:
            if c in self._records:
                self._records[c] = sorted(self._records[c], key=itemgetter(idx_start, idx_end))

    def _index_with_tabix(self):
        """Compress and index output file by Tabix"""

        pysam.tabix_compress(self._fn + '_tmp', self._fn + '.gz', force=True)
        pysam.tabix_index(self._fn + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)


    def finalize(self, options):
        """Write to file, compress and index, clean up"""

        # Sort records by CHROM, START, END
        self._sort_records()

        # Initialize file and write header
        out = open(self._fn + '_tmp', 'w')
        out.write('#createdby: ' + self._source + '\n')
        out.write('#date: ' + str(datetime.datetime.today()).split()[0] + '\n')

        # Write records to file
        for c in allowed_chroms:
            if c in self._records:
                for record in self._records[c]:
                    record = map(str, record)

                    old_record = [
                        record[0],
                        record[2],
                        record[1],
                        record[3],
                        record[5],
                        '1' if record[4] == '+' else '-1',
                        record[6],
                        record[7],
                        record[11],
                        str(int(record[9]) + 1),
                        str(int(record[10]) + 1)
                    ]

                    for e in record[8].split(','):
                        [start, end] = e.split('-')
                        old_record.append(start)
                        old_record.append(end)

                    out.write('\t'.join(old_record) + '\n')

        out.close()

        # Compress and index by Tabix
        self._index_with_tabix()

        # Remove temporary file
        os.remove(self._fn + '_tmp')
