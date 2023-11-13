# Created by alex at 08.11.23

import re
import time
from multiprocessing import Process, Queue, Pipe
from .fastq import FastqGenerator


class BarcodeCoordinate:

    def __init__(self, barcode, chromosome, strand, position):
        self.barcode = barcode
        self.chromosome = chromosome
        self.strand = strand
        self.position = position
        self.count = 0


class BarcodeCoordinateHandler:

    def __init__(self, barcode_coordinate_file, silent=False):
        self.barcode_coordinate_file = barcode_coordinate_file
        self.silent = silent
        self.barcode_coordinates = dict()

    def read_barcode_coordinates(self):
        with open(self.barcode_coordinate_file, 'r') as reader:
            for idx, line in enumerate(reader):
                tabs = line.rstrip().split("\t")
                barcode = tabs[0]
                seqId = tabs[1]
                strand = tabs[2]
                position = int(tabs[3])
                barcodeCoordinate = BarcodeCoordinate(barcode, seqId, strand, position)
                self.barcode_coordinates[barcodeCoordinate.barcode] = barcodeCoordinate
                if not self.silent:
                    if (idx + 1) % 10000 == 0:
                        print("{} coordinates parsed and assigned to gene locations".format((idx + 1)))


class BarcodeSearcher(Process):

    def __init__(self, sequence_queue: Queue, barcode_coordinates_dict: dict, left_read_sequence, right_read_sequence,
                 result_queue: Queue, pipe_writer, silent=False):
        Process.__init__(self)
        self.barcode_coordinates_dict = barcode_coordinates_dict
        self.sequence_queue = sequence_queue
        self.left_read_sequence = left_read_sequence
        self.right_read_sequence = right_read_sequence
        self.assigned_reads = 0
        self.not_assigned_reads = 0
        self.result_queue = result_queue
        self.pipe_writer = pipe_writer
        self.silent = silent

    def run(self):
        count = 0
        while True:
            try:
                fastq = self.sequence_queue.get_nowait()
                if fastq is None:
                    self.pipe_writer.send({'assigned_reads': self.assigned_reads, 'unassigned_reads': self.not_assigned_reads})
                    if not self.silent:
                        print("Searched barcodes in {} reads".format(count))
                    # poison pill for result process
                    self.result_queue.put(None)
                    break
                count += 1
                if not self.silent:
                    if count % 1000 == 0:
                        print("{} reads analyzed ".format(count))
                self.find_barcodes(fastq)
            except Exception:
                pass

    def find_barcodes(self, fastq):
        sequence = fastq.sequence
        if self.left_read_sequence and len(self.left_read_sequence) > 0:
            left_search = re.search(self.left_read_sequence, fastq.sequence)
            if left_search:
                sequence = sequence[left_search.end():]
        if self.right_read_sequence and len(self.right_read_sequence) > 0:
            right_search = re.search(self.right_read_sequence, sequence)
            if right_search:
                sequence = sequence[:right_search.start()]
        found = False
        for barcode, barcodeCoordinate in self.barcode_coordinates_dict.items():
            if barcode in sequence:
                self.assigned_reads += 1
                self.result_queue.put(barcode)
                found = True
                break
        if not found:
            self.not_assigned_reads += 1

class BarcodeAssigner:

    def __init__(self, fastq_files: [], barcodeCoordinateHandler: BarcodeCoordinateHandler,
                 left_read_sequence, right_read_sequence, threads=4, silent=False):
        self.fastq_files = fastq_files
        self.barcodeCoordinateHandler = barcodeCoordinateHandler
        self.left_read_sequence = left_read_sequence
        self.right_read_sequence = right_read_sequence
        self.threads = threads
        self.silent = silent
        self.assigned_barcodes2_reads = 0
        self.unassigned_barcodes2_reads = 0

    def assignBarcodes(self) -> (dict, int):
        barcode_counts_dict, num_fastq_reads = self.__search_barcodes()
        self.__assignCount2BarcodeCoordinates(barcode_counts_dict)
        return self.barcodeCoordinateHandler.barcode_coordinates, num_fastq_reads, self.assigned_barcodes2_reads, self.unassigned_barcodes2_reads

    def __search_barcodes(self):
        sequence_queue = Queue(maxsize=3000)
        barcode_queue = Queue(maxsize=2000)

        workers = []
        pipe_readers = []
        for _ in range(self.threads):
            pipe_reader, pipe_writer = Pipe()
            pipe_readers.append(pipe_reader)
            worker = BarcodeSearcher(sequence_queue, self.barcodeCoordinateHandler.barcode_coordinates,
                                     self.left_read_sequence, self.right_read_sequence, barcode_queue,  pipe_writer, self.silent)
            workers.append(worker)
        fastq_pipe_reader, fastq_pipe_writer = Pipe()
        fastq_reader = FastqGenerator(self.fastq_files, sequence_queue, len(workers), fastq_pipe_writer)
        fastq_reader.start()
        time.sleep(2)  # wait some seconds so that the fastq queue can fill up

        for worker in workers:
            worker.start()
        barcode_dict = self.count_barcodes(barcode_queue, len(workers))

        # wait for all processes to finish
        fastq_reader.join()
        num_fastq_reads = fastq_pipe_reader.recv()
        for worker in workers:
            worker.join()
        for pipe_reader in pipe_readers:
            barcode_stat_dict = pipe_reader.recv()
            self.assigned_barcodes2_reads += barcode_stat_dict['assigned_reads']
            self.unassigned_barcodes2_reads += barcode_stat_dict['unassigned_reads']

        return barcode_dict, num_fastq_reads

    def count_barcodes(self, barcode_queue, num_barcode_workers):
        barcode_counts_dict = dict()
        barcode_workers_finished = 0
        isRunning = True
        while isRunning:
            try:
                barcode = barcode_queue.get_nowait()
                if barcode is None:
                    barcode_workers_finished += 1
                    if barcode_workers_finished == num_barcode_workers:
                        isRunning = False
                else:
                    if barcode in barcode_counts_dict:
                        barcode_counts_dict[barcode] += 1
                    else:
                        barcode_counts_dict[barcode] = 1
            except Exception as e:
                pass
        return barcode_counts_dict

    def __assignCount2BarcodeCoordinates(self, barcodes_counts: dict):
        for barcode, barcodeCoordinate in self.barcodeCoordinateHandler.barcode_coordinates.items():
            if barcode in barcodes_counts:
                barcodeCoordinate.count = barcodes_counts[barcode]
