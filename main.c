#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "util.c"
#include "args.c"
#include "seg.c"
#include "bam.c"
#include "bed.c"


int main(int argc, char* argv[]) {
  if(parse_args(argc, argv) != 0)
    return 1;

  gsl_set_error_handler_off();
  init_tdist4();

  struct BamReader exp_bam_reader;
  open_bam_reader(exp_bam_file_name, &exp_bam_reader);
  int nChr = exp_bam_reader.header->n_targets;

  struct BamReader ctrl_bam_reader;
  int32_t ctrl_read_count = 0;
  int32_t exp_read_count = 0;
  if(ctrl_bam_file_name != 0) {
    open_bam_reader(ctrl_bam_file_name, &ctrl_bam_reader);
    if(ctrl_bam_reader.header->n_targets != nChr) {
      fprintf(stderr, "Distinct BAM headers in exp and ctrl files.\n");
      return 1;
    }
    //We need total read counts only when control dataset exists.
    exp_read_count = count_mapped_reads(&exp_bam_reader);
    if(exp_read_count < 0) {
      fprintf(stderr, "Could not count reads using bam index file\n");
      return 1;
    }
    ctrl_read_count = count_mapped_reads(&exp_bam_reader);
    if(ctrl_read_count < 0) {
      fprintf(stderr, "Could not count reads for ctrl using bam index file\n");
      return 1;
    }
  }

  struct BedReader bed_reader;
  if(unmappable_file_name != 0) {
    if(open_bed_reader(unmappable_file_name, &bed_reader) != 0)
      return 1;
  }

  int i;
  for(i = 0; i < nChr; i++) {
    char* exp_chr = exp_bam_reader.header->target_name[i];
    if(ctrl_bam_file_name != 0) {
      char* ctrl_chr = ctrl_bam_reader.header->target_name[i];
      if(strcmp(exp_chr, ctrl_chr) != 0) {
        fprintf(stderr, "Distinct BAM headers in exp and ctrl files.\n");
        return 1;
      }
    }

    struct InputData data = {0};
    data.chr = exp_chr;

    if(read_chr_bam(&exp_bam_reader, i, &data.P, &data.N) != 0)
      return 1;

    if(ctrl_bam_file_name != 0) {
      if(read_chr_bam(&ctrl_bam_reader, i, &data.PC, &data.NC) != 0)
        return 1;
    }

    if(unmappable_file_name != 0) {
      if(read_chr_bed(&bed_reader, exp_chr, &data.U) != 0)
        return 1;
    }

    process_chr(&data, exp_read_count, ctrl_read_count);
    free_input_data(&data);
  }

  close_bam_reader(&exp_bam_reader);
  if(ctrl_bam_file_name != 0)
    close_bam_reader(&ctrl_bam_reader);
  if(unmappable_file_name != 0)
    close_bed_reader(&bed_reader);

  return 0;
}
