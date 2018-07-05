
import pysam


class Reference(object):
    """Class for reference genome sequence"""


    def __init__(self, filename):
        """Constructor of the Reference class"""

        self.fasta_file = pysam.Fastafile(filename)


    def read_sequence(self, chrom, start, end):
        """Retrieve sequence of a genomic region"""

        if start < 0 or end < 0:
            return None

        if chrom not in self.fasta_file.references:
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            else:
                chrom = 'chr' + chrom
            if chrom not in self.fasta_file.references:
                return None

        if end <= start:
            return ''

        if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']:
            last = self.fasta_file.getReferenceLength(chrom)
        else:
            last = self.fasta_file.get_reference_length(chrom)

        if end > last:
            end = last

        return self.fasta_file.fetch(chrom, start, end).upper()




