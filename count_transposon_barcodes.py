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

    gff_files = [gff[0] for gff in args.gff_files]
    fastq_files = [fastq[0] for fastq in args.fastq_files]

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

    with open(args.log_file, 'w') as log_writer:
        log_writer.write(f"Total number of reads:\t{num_fastq_reads:,}\n\n")
        log_writer.write(f"Reads without barcode assignment:\t{unassigned_barcodes2_reads:,}\n")
        log_writer.write(f"Reads with barcode assignment:\t{assigned_barcodes2_reads:,}\n\n")
        log_writer.write(f"Used barcodes:\t{gff.assigned_barcodes:,}\n")
        log_writer.write(f"\tAssigned to genes:\t{gff.assigned_barcodes_to_genes:,}\n")
        log_writer.write(f"\tAssigned to intergenic regions:\t{gff.assigned_barcodes_to_intergenic:,}\n")
        log_writer.write(f"Discarded barcodes:\t{gff.unassigned_barcodes:,}\n")
        log_writer.close()
