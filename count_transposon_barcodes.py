import argparse
from lib.barcode import *
from lib.gff import *


def range_limited_float_type(arg):
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < 0.0 or f > 1.0:
        raise argparse.ArgumentTypeError("Argument must be between 0.0 - 1.0")
    return f


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Count transposon mutant libraries",
        usage='use "python3 count_transposon_barcodes.py --help" for more information',
        epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fastq_files', action="append", nargs="+", metavar="FILE", type=str,
                        help="Fastq file (one ore more)",
                        required=True)
    parser.add_argument('-c', '--barcode_coordinate_file', metavar="FILE", type=str,
                        help="Barcode coordinate file",
                        required=True)
    parser.add_argument('-g', '--gff_files', action="append", nargs="+", metavar="FILE", type=str,
                        help="GFF3 gene annotation file (one ore more)",
                        required=True)
    parser.add_argument('-o', '--output_file', metavar="FILE", type=str,
                        help="Output file",
                        required=True)
    parser.add_argument('-a', '--output_annotated_gff_file', metavar="FILE", type=str,
                        help="Output file of generated GFF file",
                        required=True)
    parser.add_argument('--log_file', metavar="FILE", type=str,
                        help="Output of count statistic", default="count_statistic.txt",
                        required=False)
    parser.add_argument('-l', '--left_flanking_sequence', metavar="SEQUENCE", type=str,
                        help="Sequence flanking the 5'-region of the barcode in each read",
                        required=False, default="ACTCACTATAGGGAGACCGGCCT")
    parser.add_argument('-r', '--right_flanking_sequence', metavar="SEQUENCE", type=str,
                        help="Sequence flanking the 3'-region of the barcode in each read",
                        required=False, default="CAGGGTTGAGATGTGTATA")
    parser.add_argument('-3', '--fraction_discard_gene_3_end', metavar="FLOAT", type=range_limited_float_type,
                        help="3'-end fraction to ignore for each gene [0.0 .. 1.0]",
                        required=False, default=0.1)
    parser.add_argument('-5', '--fraction_discard_gene_5_end', metavar="FLOAT", type=range_limited_float_type,
                        help="5'-end fraction to ignore for each gene [0.0 .. 1.0]",
                        required=False, default=0.1)
    parser.add_argument('-t', '--threads', metavar="INT", type=int,
                        help="Number of threads", required=False, default=4)
    parser.add_argument('-s', '--silent', action='store_true',
                        help="Supress any intermediate output")

    args = parser.parse_args()

    gff_files = []
    fastq_files = []
    for fastq_file in args.fastq_files:
        fastq_files.append(fastq_file[0])

    for gff_file in args.gff_files:
        gff_files.append(gff_file[0])

    gff = GFF(gff_files, args.fraction_discard_gene_5_end, args.fraction_discard_gene_3_end, args.silent)
    gff.read_gff()

    barcodeCoordinateHandler = BarcodeCoordinateHandler(args.barcode_coordinate_file, args.silent)
    barcodeCoordinateHandler.read_barcode_coordinates()

    barcodeAssigner = BarcodeAssigner(fastq_files, barcodeCoordinateHandler,
                                      args.left_flanking_sequence, args.right_flanking_sequence, args.threads, args.silent)
    barcode_coordinate_dict, num_fastq_reads, assigned_barcodes2_reads, unassigned_barcodes2_reads = (
        barcodeAssigner.assignBarcodes())

    gff.assignCoordinates2GFF(barcode_coordinate_dict)
    gff.write_output(args.output_file)
    gff.write_gene_annotation(args.output_annotated_gff_file)

    with open(args.stat_file, 'w') as stat_writer:
        stat_writer.write("Total number of reads:\t{:n}\n\n".format(num_fastq_reads))
        stat_writer.write("Reads without barcode assignment:\t{:n}\n".format(unassigned_barcodes2_reads))
        stat_writer.write("Reads with barcode assignment:\t{:n}\n\n".format(assigned_barcodes2_reads))
        stat_writer.write("Used barcodes:\t{:n}\n".format(gff.assigned_barcodes))
        stat_writer.write("\tAssigned to genes:\t{:n}\n".format(gff.assigned_barcodes_to_genes))
        stat_writer.write("\tAssigned to intergenic regions:\t{:n}\n".format(gff.assigned_barcodes_to_intergenic))
        stat_writer.write("Discarded barcodes:\t{:n}\n".format(gff.unassigned_barcodes))
        stat_writer.close()
