import helper
import reference
import sys
import transcripts
from main.version import __version__


def run(options):

    if not ((options.series.startswith('CART37') or options.series.startswith('CART38')) and len(options.series) == 7):
        print '\nSeries code incorrect!\n'
        quit()

    print '\n==== ENSTWriter {} '.format(__version__) + '=' * 78

    # Initialize reference sequence reader
    ref = reference.Reference(options.ref)

    # Initialize transcript database writer
    tdb_writer = helper.initialize_transcript_db_writer(options)

    # Read Ensembl database
    ensembl_db = transcripts.read_ensembl_db(options.ensembl)

    # Read previous CAVA db output and reference genome if required
    if options.prev_cava_db:
        prev_ref = reference.Reference(options.prev_ref)
        prev_cava_db = helper.read_prev_cava_db(options.prev_cava_db, prev_ref)
    else:
        prev_cava_db = None

    # Initialize output files
    out_genepred, out_fasta, out_genepred_annovar, out_fasta_annovar, gbk_dir = helper.initialize_output_files(options)

    # Initialize progress info
    sys.stdout.write('\nProcessing {} CARTs read from {} ... '.format(
            helper.number_of_input_carts(options.input), options.input
        )
    )
    sys.stdout.flush()

    # Initialize CART numbering
    cartidx = 10000 if options.prev_cava_db is None else helper.get_last_cartidx(options.prev_cava_db)

    # Iterate through input records
    missing_list = []
    gff3_lines = {}
    for line in open(options.input):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        hgnc_id = cols[0][5:]
        enst = cols[1]

        # Add ENST to missing list if not found in Ensembl database
        if enst not in ensembl_db:
            missing_list.append('{} (HGNC:{})'.format(enst, hgnc_id))
            continue

        # Retrieve data about ENST
        transcript = ensembl_db[enst]

        # Calculating CART ID
        if options.prev_cava_db is None:
            cartidx += 1
            cart_id = '{}{}'.format(options.series, cartidx)
        else:
            content = (
                transcript.strand,
                len(transcript.exons),
                helper.read_mrna_sequence(transcript, ref)
            )

            if hgnc_id in prev_cava_db and content == prev_cava_db[hgnc_id]['content']:
                cart_id = '{}{}'.format(options.series, prev_cava_db[hgnc_id]['cartidx'])
            else:
                cartidx += 1
                cart_id = '{}{}'.format(options.series, cartidx)

        # Add CART ID and HGNC ID to transcript
        transcript.id = cart_id
        transcript.hgnc_id = hgnc_id

        # Add transcript to database writer
        tdb_writer.add(transcript)

        # Create content of gff3 file
        gff3_lines = helper.create_gff3_lines(transcript, gff3_lines)

        # Write to gp file
        helper.output_genepred(transcript, out_genepred)

        # Write to gbk output
        if options.gbk:
            helper.output_gbk(transcript, ref, gbk_dir)

        # Write to fasta file
        helper.output_fasta(transcript, out_fasta, ref)

        # Write annovar files
        if options.annovar:
            helper.output_genepred(transcript, out_genepred_annovar)
            helper.output_fasta_annovar(transcript, out_fasta_annovar, ref)

    # Create bgzipped, Tabix-index GFF3 output
    helper.output_gff3(gff3_lines, options.series + '.gff')

    # Finalize outputs
    helper.finalize_outputs(
        options,
        tdb_writer,
        out_fasta,
        out_genepred,
        out_genepred_annovar,
        out_fasta_annovar,
        gbk_dir
    )

    # Print out summary info
    helper.print_summary_info(options, missing_list)

    print '\n' + '=' * 100 + '\n'
