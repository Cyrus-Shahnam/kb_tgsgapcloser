# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from kb_tgsgapcloser.tgsgapcloser_utils import TGSGapCloserRunner
#END_HEADER


class kb_tgsgapcloser:
    '''
    Module Name:
    kb_tgsgapcloser

    Module Description:
    A KBase module: kb_tgsgapcloser
Fill N-gaps in scaffolded genome assemblies with long reads (or pre-assembled
contigs) using TGS-GapCloser v1.2.1, with optional Racon correction of the
candidate fill sequences.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.0"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_tgsgapcloser(self, ctx, params):
        """
        :param params: instance of type "TGSGapCloserParams" -> structure:
           parameter "workspace_name" of String, parameter "assembly_ref" of
           type "obj_ref", parameter "reads_ref" of type "obj_ref", parameter
           "output_assembly_name" of String, parameter "read_type" of String,
           parameter "error_correction" of String, parameter "racon_rounds"
           of Long, parameter "minimap2_preset" of String, parameter
           "min_idy" of Double, parameter "min_match" of Long, parameter
           "min_nread" of Long, parameter "max_candidate" of Long, parameter
           "chunk" of Long, parameter "gapsize_check" of Long, parameter
           "threads" of Long
        :returns: instance of type "TGSGapCloserOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "assembly_ref" of type "obj_ref", parameter
           "gapfill_summary" of unspecified object
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_tgsgapcloser
        runner = TGSGapCloserRunner(self.config, self.callback_url, ctx['token'])
        output = runner.run(params)
        #END run_tgsgapcloser

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_tgsgapcloser return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
