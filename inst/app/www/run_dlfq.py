import directlfq.lfq_manager as lfq_manager
import sys
platform= sys.argv[1]
raw = sys.argv[2] # fragpipe protein.tsv
evid = sys.argv[3] #  fragpipe ion.tsv
#example_input_file_diann = "/path/to/example_input_file_diann.tsv"

# platform = 'FragPipe'
# raw = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe/combined_protein.tsv'
# evid = 'E:/proteomics/maus2/test/OpDEA-master/R/example/FragPipe/combined_ion.tsv'

if platform == 'FragPipe':
  lfq_manager.run_lfq(evid, num_cores=1)
elif platform == 'Maxquant':
  lfq_manager.run_lfq(input_file=evid, mq_protein_groups_txt=raw, num_cores=1)
elif platform == 'DIANN':
  lfq_manager.run_lfq(evid, num_cores=1)
elif platform == 'spt':
  lfq_manager.run_lfq(evid, num_cores=1)
