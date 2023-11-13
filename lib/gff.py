# Created by alex at 08.11.23
import math

INTERGENIC_NAME = "intergenic"


class GFFEntry:

    def __init__(self, seqid, source, type, start, end, score, strand, phase, attributes, percentage_five_prime_trim=0.0, percentage_three_prime_trim=0.0):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.percentage_five_prime_trim = percentage_five_prime_trim
        self.percentage_three_prime_trim = percentage_three_prime_trim
        self.start_trimmed = self.__get_start_trimmed()
        self.end_trimmed = self.__get_end_trimmed()
        self.read_count = 0
        self.insertion_site_positions = set()
        self.attributes_dict = dict()
        self.__parse_attributes()

    def __get_start_trimmed(self):
        gene_length = self.gene_length()
        if self.type == INTERGENIC_NAME:  # in case of a computed intergenic region don't trim
            return self.start
        if self.strand == '+':
            num_trim_nucleotides = gene_length * self.percentage_five_prime_trim
        else:
            num_trim_nucleotides = gene_length * self.percentage_three_prime_trim
        return math.floor(self.start + num_trim_nucleotides)

    def __get_end_trimmed(self):
        gene_length = self.gene_length()
        if self.type == INTERGENIC_NAME:  # in case of a computed intergenic region don't trim
            return self.end
        if self.strand == '+':
            num_trim_nucleotides = gene_length * self.percentage_three_prime_trim
        else:
            num_trim_nucleotides = gene_length * self.percentage_five_prime_trim
        return math.ceil(self.end - num_trim_nucleotides)

    def gene_length(self) -> int:
        return (self.end - self.start) + 1

    def get_ins_index(self):
        """
        insertion indexes:  unique insertion sites sites divided by gene length
        """
        if len(self.insertion_site_positions) > 0:
            return len(self.insertion_site_positions) / self.gene_length()
        else:
            return 0

    def get_gene_name(self) -> str:
        if 'gene' in self.attributes_dict:
            return self.attributes_dict['gene']
        elif 'Name' in self.attributes_dict:
            return self.attributes_dict['Name']
        else:
            return self.get_locus_tag()

    def get_locus_tag(self) -> str:
        if 'locus_tag' in self.attributes_dict:
            return self.attributes_dict['locus_tag']

    def __parse_attributes(self) -> None:
        attributes = self.attributes.split(";")
        for attribute in attributes:
            key_value = attribute.split("=")
            self.attributes_dict[key_value[0]] = key_value[1]

    def __str__(self) -> str:
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.seqid, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes)


class GFF:

    def __init__(self, gff_files: [], percentage_ignore_5_prime_end, percentage_ignore_3_prime_end, silent=False):
        self.gff_files = gff_files
        self.percentage_ignore_5_prime_end = percentage_ignore_5_prime_end
        self.percentage_ignore_3_prime_end = percentage_ignore_3_prime_end
        self.gff_dict = dict()
        self.genomeLengthDict = dict()
        self.assigned_barcodes = 0
        self.unassigned_barcodes = 0
        self.assigned_barcodes_to_genes = 0
        self.assigned_barcodes_to_intergenic = 0
        self.silent = silent

    def read_gff(self) -> None:
        for gff_file in self.gff_files:
            with open(gff_file, 'r') as reader:
                for idx, line in enumerate(reader):
                    if not line.startswith("#"):
                        tabs = line.rstrip().split("\t")
                        if len(tabs) >= 9:
                            if tabs[2] == "region":  # whole genome
                                self.genomeLengthDict[tabs[0]] = int(tabs[4])
                            if tabs[2] == "gene":
                                seqid = tabs[0]
                                source = tabs[1]
                                type = tabs[2]
                                start = int(tabs[3])
                                end = int(tabs[4])
                                score = tabs[5]
                                strand = tabs[6]
                                phase = tabs[7]
                                attributes = tabs[8]
                                gff = GFFEntry(seqid, source, type, start, end, score, strand, phase, attributes,
                                               self.percentage_ignore_5_prime_end, self.percentage_ignore_3_prime_end)
                                if gff.seqid in self.gff_dict:
                                    self.gff_dict[gff.seqid].append(gff)
                                else:
                                    self.gff_dict[gff.seqid] = [gff]
        self.__sort_gff()
        self.__fill_intergenic_space()

    def __sort_gff(self):
        for seqId, gff_arr in self.gff_dict.items():
            gff_arr.sort(key=lambda gff: gff.start)

    def __fill_intergenic_space(self):
        attributes_template = ("ID=" + INTERGENIC_NAME + "-{}-{};Name=" + INTERGENIC_NAME + "-{}-{};gene_biotype=" + INTERGENIC_NAME +
                               ";gbkey=" + INTERGENIC_NAME + ";locus_tag=" + INTERGENIC_NAME + "-{}-{}")
        intergenic_dict = dict()
        for seqId, gff_arr in self.gff_dict.items():
            intergenic_gffs = []
            intergenic_id = 1
            last_gff = None
            for idx, gff in enumerate(gff_arr):
                attributes = attributes_template.format(intergenic_id, seqId, intergenic_id, seqId, intergenic_id, seqId)
                if idx > 0:
                    previous_gff = gff_arr[idx - 1]
                    start = previous_gff.end + 1
                    end = gff.start - 1
                    if start < end:
                        intergenic_gff = GFFEntry(seqId, 'Computed', INTERGENIC_NAME, start, end, '.', '+', '.', attributes)
                        intergenic_gffs.append(intergenic_gff)
                        intergenic_id += 1
                elif gff.start != 0:  # beginning
                    intergenic_gff = GFFEntry(seqId, 'Computed', INTERGENIC_NAME, 0, gff.start - 1, '.', '+', '.', attributes)
                    intergenic_gffs.append(intergenic_gff)
                    intergenic_id += 1
                last_gff = gff
            if self.genomeLengthDict[seqId] and last_gff.end != self.genomeLengthDict[seqId]:  # last element till the end of the genome
                attributes = attributes_template.format(intergenic_id, seqId, intergenic_id, seqId, intergenic_id, seqId)
                intergenic_gff = GFFEntry(seqId, 'Computed', INTERGENIC_NAME, last_gff.end + 1, self.genomeLengthDict[seqId], '.', '+', '.', attributes)
                intergenic_gffs.append(intergenic_gff)
            intergenic_dict[seqId] = intergenic_gffs
        for seqId, intergenic_gffs in intergenic_dict.items():
            self.gff_dict[seqId].extend(intergenic_gffs)
            self.gff_dict[seqId].sort(key=lambda x: x.start)

    def findGFFEntryByPosition(self, seqId, position):
        if seqId in self.gff_dict:
            for gff in self.gff_dict[seqId]:
                if gff.start_trimmed <= position <= gff.end_trimmed:
                    return gff
        else:
            print("WARNING: Sequence Id '{}' not found in GFF file".format(seqId))
        # print("WARNING: No GFF entry for Position '{}' on SeqID '{}' ".format(position, seqId))
        return None

    def assignCoordinates2GFF(self, barcode_coordinates_dict: dict):
        for idx, (barcode, barcodeCoordinate) in enumerate(barcode_coordinates_dict.items()):
            if (idx + 1) % 10000 == 0:
                if not self.silent:
                    print("{}/{} barcodes processed".format(idx + 1, len(barcode_coordinates_dict.keys())))
            if barcodeCoordinate.count > 0:
                gff = self.findGFFEntryByPosition(barcodeCoordinate.chromosome, barcodeCoordinate.position)
                if gff:
                    gff.read_count += barcodeCoordinate.count
                    gff.insertion_site_positions.add(barcodeCoordinate.position)
                    self.assigned_barcodes += 1
                    if gff.type == INTERGENIC_NAME:
                        self.assigned_barcodes_to_intergenic += 1
                    else:
                        self.assigned_barcodes_to_genes += 1
                else:
                    self.unassigned_barcodes += 1

    def write_output(self, output_file):
        with open(output_file, 'w') as writer:
            writer.write("locus_tag\tgene_name\tncrna\tstart\tend\tstrand\tread_count\tins_index\tgene_length\tins_count\tfcn\n")
            for seqId, gff_arr in self.gff_dict.items():
                for gff in gff_arr:
                    strand = '1' if gff.strand == "+" else '-1'
                    writer.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"{}"\n'.format(
                        gff.get_locus_tag(), gff.get_gene_name(), 0, gff.start, gff.end,
                        strand, gff.read_count, gff.get_ins_index(), gff.gene_length(), len(gff.insertion_site_positions), gff.get_gene_name()
                    ))
            writer.close()

    def write_gene_annotation(self, output_gff_file):
        with open(output_gff_file, 'w') as writer:
            writer.write("##gff-version 3\n")
            for seqId, gff_arr in self.gff_dict.items():
                for gff in gff_arr:
                    writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        gff.seqid, gff.source, gff.type, gff.start, gff.end, gff.score, gff.strand, gff.phase, gff.attributes))
            writer.close()
