from __future__ import print_function
import os, sys

class logger(object):
    def __init__(self, folder, verbose=True):
        super(logger, self).__init__()
        self.sfh = open(os.path.join(folder, "rejected_sequences.log"), 'w')
        self.sfh.write("name\tsegment\treason\tcomment\n")
        self.verbose = verbose

    def notify(self, msg):
        print("*** notification *** {}".format(msg))

    def debug(self, msg):
        print("*** debug *** {}".format(msg))

    def warn(self, msg):
        print("***   WARN   *** {}".format(msg))

    def fatal(self, msg):
        print("***   FATAL      *** {}".format(msg))
        sys.exit(2)

    def drop(self, seqName, segment, filterName, comment=""):
        self.sfh.write("{}\t{}\t{}\t{}\n".format(seqName, segment, filterName, comment))
        if self.verbose:
            self.notify("Dropping strain {} (segment {}) for {}".format(seqName, segment, filterName))
