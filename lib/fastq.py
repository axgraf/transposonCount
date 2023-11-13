# Created by alex at 08.11.23
from multiprocessing import Queue, Process


class Fastq:

    def __init__(self, header, sequence, delimiter, quality):
        self.header = header
        self.sequence = sequence
        self.delimiter = delimiter
        self.quality = quality


class FastqGenerator(Process):

    def __init__(self, fastq_files: [], sequence_queue: Queue, num_workers: int, pipe_writer):
        Process.__init__(self)
        self.fastq_files = fastq_files
        self.sequence_queue = sequence_queue
        self.num_workers = num_workers
        self.pipe_writer = pipe_writer

    def run(self):
        for fastq_file in self.fastq_files:
            with open(fastq_file, 'r') as reader:
                count = 0
                while True:
                    header = reader.readline().rstrip()
                    sequence = reader.readline().rstrip()
                    delimiter = reader.readline().rstrip()
                    quality = reader.readline().rstrip()
                    fastq = Fastq(header, sequence, delimiter, quality)
                    #if count >= 1000:
                    #    sequence = ""
                    if len(sequence) == 0:
                        print("Number of reads:\t{}".format(count))
                        self.pipe_writer.send(count)
                        break
                    count += 1
                    self.sequence_queue.put(fastq)
                reader.close()
                # send poison pill to all barcodes assignments processes
                for _ in range(self.num_workers):
                    # add poison pill
                    self.sequence_queue.put(None)
