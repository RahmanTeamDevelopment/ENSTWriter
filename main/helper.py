
from __future__ import division
from operator import itemgetter
import pysam
import datetime
import os
import transcripts
import shutil
import gzip
from main.version import __version__


def create_gff3_lines(transcript, gff3_lines):

    if transcript.chrom not in gff3_lines:
        gff3_lines[transcript.chrom] = []

    attr = ';'.join(['ID=' + transcript.id, 'hgnc_id=' + transcript.hgnc_id, 'gene_symbol=' + transcript.gene_symbol, 'biotype=protein_coding'])
    gff3_lines[transcript.chrom].append([transcript.chrom, '.', 'transcript', transcript.start + 1, transcript.end, '.', transcript.strand, '.', attr])

    # Exons
    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        exon_id = 'EXON' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + exon_id, 'Parent=' + transcript.id])
        gff3_lines[transcript.chrom].append([transcript.chrom, '.', 'exon', exon.start + 1, exon.end, '.', transcript.strand, '.', attr])

    # CDS
    cds_regs = transcript.cds_regions()
    cdspos = 0
    for i in range(len(cds_regs)):
        cds_reg = cds_regs[i]

        cds_id = 'CDS' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + cds_id, 'Parent=' + transcript.id])

        if cdspos % 3 == 0:
            phase = 0
        elif cdspos % 3 == 1:
            phase = 2
        else:
            phase = 1

        gff3_lines[transcript.chrom].append([transcript.chrom, '.', 'CDS', cds_reg[0] + 1, cds_reg[1], '.', transcript.strand, str(phase), attr])

        cdspos += cds_reg[1] - cds_reg[0]

    return gff3_lines


def output_gff3(gff3_lines, fn):

    out = open(fn + '_tmp', 'w')
    out.write('##gff-version 3\n')
    for c in map(str,range(1, 23))+['X', 'Y', 'MT']:
        if c not in gff3_lines:
            continue
        gff3_lines[c] = sorted(gff3_lines[c], key=itemgetter(3, 4))

        for x in gff3_lines[c]:
            x = map(str, x)
            out.write('\t'.join(x)+'\n')
    out.close()

    pysam.tabix_compress(fn + '_tmp', fn + '.gz', force=True)
    pysam.tabix_index(fn + '.gz', seq_col=0, start_col=3, end_col=4, meta_char='#', force=True)

    os.remove(fn + '_tmp')


def output_genepred(transcript, outfile):

    exons = transcript.exons
    if transcript.strand == '-':
        exons = exons[::-1]

    record = [
                transcript.id,
                'chr'+transcript.chrom,
                transcript.strand,
                str(transcript.start),
                str(transcript.end),
                str(min(transcript.coding_start, transcript.coding_end)),
                str(max(transcript.coding_start, transcript.coding_end)+1),
                str(len(transcript.exons)),
                ''.join([str(e.start)+',' for e in exons]),
                ''.join([str(e.end)+',' for e in exons]),
                '0',
                'HGNC:'+transcript.hgnc_id,
                'cmpl',
                'cmpl',
                ''.join(frame_offsets(transcript))
            ]

    outfile.write('\t'.join(record) + '\n')


def output_gbk(transcript, ref, dirname):

    outfile = open('{}/{}_{}.gbk'.format(dirname, transcript.gene_symbol, transcript.id), 'w')
    full_dna = read_full_dna_sequence(transcript, ref)
    if transcript.strand == '-':
        full_dna = reverse_complement(full_dna)
    gbk_header(transcript, outfile, full_dna)
    gbk_features(transcript, outfile, ref)
    gbk_origin(outfile, full_dna)
    outfile.write('//\n')
    outfile.close()


def output_fasta(transcript, outfile, ref):

    width = 70
    s = read_mrna_sequence(transcript, ref)
    outfile.write('>{}\n'.format(transcript.id))
    while len(s)>0:
        outfile.write(s[:width]+'\n')
        s = s[width:]


def output_fasta_annovar(transcript, outfile, ref):

    width = 70
    s = read_mrna_sequence(transcript, ref)
    outfile.write('>{}#{}#{}\n'.format(transcript.id, transcript.chrom, str(transcript.start)))
    while len(s)>0:
        outfile.write(s[:width]+'\n')
        s = s[width:]


def get_coding_sequence(transcript, ref):

    ret = ''
    for e in transcript.cds_regions():
        s = ref.read_sequence(transcript.chrom, e[0], e[1])
        if transcript.strand == '-':
            s = reverse_complement(s)
        ret += s
    return ret


def get_protein_sequence(transcript, ref):

    return translate(get_coding_sequence(transcript, ref))[:-1]


def translate(dna):

    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}

    ret = ''
    index = 0
    while index + 3 <= len(dna):
        codon = dna[index:index + 3].upper()
        if 'N' in codon:
            ret += '?'
            index += 3
            continue
        ret += gencode[codon]
        index += 3
    return ret


def gbk_header(transcript, outfile, full_dna):

    now = datetime.datetime.now()
    date = '{}-{}-{}'.format(now.strftime("%d"),now.strftime("%b").upper(),now.strftime("%Y"))
    outfile.write('LOCUS       {}              {} bp    DNA . {}\n'.format(transcript.id, len(full_dna), date))
    outfile.write('ACCESSION   {}\n'.format(transcript.id))
    outfile.write('VERSION     {}\n'.format(transcript.id))
    outfile.write('KEYWORDS    .\n')
    outfile.write('SOURCE      Homo sapiens (human)\n')
    outfile.write('  ORGANISM  Homo sapiens\n')
    outfile.write('            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n')
    outfile.write('            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;\n')
    outfile.write('            Catarrhini; Hominidae; Homo.\n')


def gbk_features(transcript, outfile, ref):

    outfile.write('FEATURES             Location/Qualifiers\n')
    feature_mrna(transcript, outfile)
    feature_cds(transcript, outfile, ref)


def feature_mrna(transcript, outfile):

    exon_strings = []
    for e in transcript.exons:
        if transcript.strand == '+':
            exon_start = e.start - transcript.start + 1
            exon_end = e.end - transcript.start
        else:
            exon_start = transcript.end - e.end + 1
            exon_end = transcript.end - e.start
        exon_string = '{}..{}'.format(exon_start, exon_end)
        exon_strings.append(exon_string)

    if len(exon_strings) == 1:
        outfile.write('     mRNA            {}\n'.format(exon_strings[0]))
        return

    first = True
    while len(exon_strings) > 0:
        s = ','.join(exon_strings[:2])
        if first:
            s = '     mRNA            join(' + s
            first = False
        else:
            s = '                     ' + s
        exon_strings = exon_strings[2:]
        if len(exon_strings) > 0:
            s += ','
        else:
            s += ')'
        outfile.write('{}\n'.format(s))


def feature_cds(transcript, outfile, ref):

    exon_strings = []
    for e in transcript.cds_regions():
        if transcript.strand == '+':
            exon_start = e[0] - transcript.start + 1
            exon_end = e[1] - transcript.start
        else:
            exon_start = transcript.end - e[1] + 1
            exon_end = transcript.end - e[0]
        exon_string = '{}..{}'.format(exon_start, exon_end)
        exon_strings.append(exon_string)

    if len(exon_strings) == 1:
        outfile.write('     CDS             {}\n'.format(exon_strings[0]))
        comment = '/translation=\"{}\"'.format(get_protein_sequence(transcript, ref))
        while len(comment) > 0:
            line = comment[:58]
            outfile.write('                     {}\n'.format(line))
            comment = comment[58:]
        return

    first = True
    while len(exon_strings) > 0:
        s = ','.join(exon_strings[:2])
        if first:
            s = '     CDS             join(' + s
            first = False
        else:
            s = '                     ' + s
        exon_strings = exon_strings[2:]
        if len(exon_strings) > 0:
            s += ','
        else:
            s += ')'
        outfile.write('{}\n'.format(s))

    comment = '/translation=\"{}\"'.format(get_protein_sequence(transcript, ref))
    while len(comment) > 0:
        line = comment[:58]
        outfile.write('                     {}\n'.format(line))
        comment = comment[58:]


def gbk_origin(outfile, full_dna):

    outfile.write('ORIGIN\n')
    s = full_dna
    lineidx = 0
    while len(s) > 0:
        line = s[:60]
        num = str(lineidx * 60 + 1)
        num = ' ' * (9 - len(num)) + num
        outfile.write('{} {} {} {} {} {} {}\n'.format(num, line[:10], line[10:20], line[20:30], line[30:40], line[40:50], line[50:]))
        lineidx += 1
        s = s[60:]


def read_mrna_sequence(transcript, ref):

    ret = ''

    for e in transcript.exons:
        e_seq = ref.read_sequence(transcript.chrom, e.start, e.end)
        if transcript.strand == '-':
            e_seq = reverse_complement(e_seq)
        ret += e_seq

    return ret


def read_full_dna_sequence(transcript, ref):

    return ref.read_sequence(transcript.chrom, transcript.start, transcript.end)


def frame_offsets(transcript):

    ret = []
    cds_sum = 0
    for exon in transcript.exons:
        cds_region = exon.get_cds(transcript.coding_start, transcript.coding_end)
        if cds_region:
            ret.append(str(cds_sum % 3)+',')
            cds_sum += cds_region[1] - cds_region[0]
        else:
            ret.append('-1,')
    if transcript.strand == '+':
        return ret
    else:
        return ret[::-1]


def reverse_complement(dna_seq):

    ret = ''
    for base in dna_seq[::-1]:
        ret += complement(base)
    return ret


def complement(b):

    comps = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "K": "M",
        "M": "K",
        "S": "S",
        "W": "W",
        "R": "Y",
        "Y": "R",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
    if b in comps:
        return comps[b]
    elif b.upper() in comps:
        x = comps[b.upper()]
        return x.lower()


def number_of_input_carts(fn):

    ret = 0
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        ret += 1
    return ret


def initialize_output_files(options):

    out_genepred = open(options.series + '.gp', 'w')
    out_fasta = open(options.series + '.fa', 'w')

    # Initializing output files required by Annovar
    if options.annovar:
        out_genepred_annovar = open('{}_refGene.txt'.format(options.series), 'w')
        out_fasta_annovar = open('{}_refGeneMrna.fa'.format(options.series), 'w')
    else:
        out_genepred_annovar = out_fasta_annovar = None

    # Initializing GBK ourput
    gbk_dir = '{}_gbk'.format(options.series)
    if options.gbk:
        if os.path.exists(gbk_dir):
            shutil.rmtree(gbk_dir)
        os.makedirs(gbk_dir)

    return out_genepred, out_fasta, out_genepred_annovar, out_fasta_annovar, gbk_dir


def finalize_outputs(options, tdb_writer, out_fasta, out_genepred, out_genepred_annovar, out_fasta_annovar, gbk_dir):

    tdb_writer.finalize(options)

    out_fasta.close()

    pysam.faidx(options.series + '.fa')

    out_genepred.close()

    if options.annovar:
        out_genepred_annovar.close()
        out_fasta_annovar.close()
        pysam.faidx('{}_refGeneMrna.fa'.format(options.series))

    if options.gbk:
        shutil.make_archive('{}_gbk'.format(options.series), "zip", './', gbk_dir)
        shutil.rmtree(gbk_dir)


def initialize_transcript_db_writer(options):

    columns = [
        'ID',
        'HGNC_ID',
        'GENE_SYMBOL',
        'INFO',
        'STRAND',
        'CHROM',
        'START',
        'END',
        'EXONS',
        'CODING_START',
        'CODING_END',
        'CDNA_CODING_START'
    ]
    tdb_writer = transcripts.TranscriptDBWriter(
        options.series + '_cava',
        source='CARTWriter ' + __version__,
        build='GRCh37',
        columns=columns
    )
    return tdb_writer


def print_summary_info(options, missing_list):

    print 'done'

    if len(missing_list) > 0:
        print '\nMissing from the Ensembl database:'
        for x in missing_list:
            print '  - ' + x

    print '\nOutput files:'
    print ' - CAVA database: {}_cava.gz (+.tbi)'.format(options.series)
    print ' - GFF3 file: {}.gff.gz (+.tbi)'.format(options.series)
    print ' - GenePred file: {}.gp'.format(options.series)
    print ' - FASTA file: {}.fa (+.fai)'.format(options.series)
    if options.annovar:
        print ' - Annovar GenePred file: {}'.format('{}_refGene.txt'.format(options.series))
        print ' - Annovar FASTA file: {} (+.fai)'.format('{}_refGeneMrna.fa'.format(options.series))
    if options.gbk:
        print ' - GenBank (GBK) files: {}_gbk.zip'.format(options.series)

    if len(missing_list) > 0:
        print '\nWarning! {} CARTs were missing! (see list above)'.format(len(missing_list))


def get_last_cartidx(fn):

    ret = -1
    for line in gzip.open(fn):

        line = line.strip()
        if line == '' or line[0] == '#':
            continue

        cols = line.split()

        cartidx = int(cols[0][-5:])

        if cartidx > ret:
            ret = cartidx

    return ret


def read_prev_cava_db(fn, ref):

    ret = {}
    for line in gzip.open(fn):

        line = line.strip()
        if line == '' or line[0] == '#':
            continue

        transcript = transcripts.Transcript(line)

        cols = line.split()
        hgnc_id = cols[2]

        ret[hgnc_id] = {
            'cartidx': cols[0][-5:],
            'content': (
                transcript.strand,
                len(transcript.exons),
                read_mrna_sequence(transcript, ref)
            )
        }

    return ret

