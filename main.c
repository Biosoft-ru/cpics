#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

#include "util.c"
#include "args.c"
#include "seg.c"
#include "bam.c"
#include "bed.c"
#include "chrs.c"
#include "fdr.c"

int main(int argc, char* argv[]) {
  if(parse_args(argc, argv) != 0)
    return 1;

  gsl_set_error_handler_off();
  init_tdist4();

  struct BamReader exp_bam_reader;
  open_bam_reader(Opt.exp_bam_file_name, &exp_bam_reader);

  char** chrs = 0;
  int n_chrs;
  if(Opt.chr_file_name) {
    chrs = load_chromosomes(&n_chrs);
  } else {
    chrs = exp_bam_reader.header->target_name;
    n_chrs = exp_bam_reader.header->n_targets;
  }
  int32_t exp_chr_ids[n_chrs];
  match_chrs(&exp_bam_reader, chrs, n_chrs, exp_chr_ids);

  struct BamReader ctrl_bam_reader;
  int32_t ctrl_chr_ids[n_chrs];
  int32_t ctrl_read_count = 0;
  int32_t exp_read_count = 0;
  if(Opt.ctrl_bam_file_name != 0) {
    open_bam_reader(Opt.ctrl_bam_file_name, &ctrl_bam_reader);
    match_chrs(&ctrl_bam_reader, chrs, n_chrs, ctrl_chr_ids);

    //We need total read counts only when control dataset exists.
    exp_read_count = count_mapped_reads(&exp_bam_reader, exp_chr_ids, n_chrs);
    ctrl_read_count = count_mapped_reads(&ctrl_bam_reader, ctrl_chr_ids, n_chrs);

    //disable score filter to calculate fdr
    cpics_disable_score_filter();
  }

  struct BedReader bed_reader;
  if(Opt.unmappable_file_name != 0) {
    if(open_bed_reader(Opt.unmappable_file_name, &bed_reader) != 0) {
      fprintf(stderr, "Can not read bed file: %s\n", Opt.unmappable_file_name);
      return 1;
    }
  }

  open_out_file();
  if(Opt.ctrl_bam_file_name != 0)
    open_tmp_files("w");

  int i;
  for(i = 0; i < n_chrs; i++) {
    char* exp_chr = chrs[i];
    fprintf(stderr, "Processing %s (%i of %i)\n", exp_chr, i+1, n_chrs);

    struct InputData data = {0};
    data.chr = exp_chr;

    read_chr_bam(&exp_bam_reader, exp_chr_ids[i], &data.P, &data.N);
    fprintf(stderr, "Load %td positive strand aligns\n", data.P.e - data.P.s);
    fprintf(stderr, "Load %td negative strand aligns\n", data.N.e - data.N.s);

    if(Opt.ctrl_bam_file_name != 0) {
      read_chr_bam(&ctrl_bam_reader, ctrl_chr_ids[i], &data.PC, &data.NC);
      fprintf(stderr, "Load %td positive strand aligns from ctrl\n", data.PC.e - data.PC.s);
      fprintf(stderr, "Load %td negative strand aligns from ctrl\n", data.NC.e - data.NC.s);
    }

    if(Opt.unmappable_file_name != 0) {
      if(read_chr_bed(&bed_reader, exp_chr, &data.U) != 0)
        return 1;
      fprintf(stderr, "Load %td unmapped regions\n", data.U.e - data.U.s);
    }

    if(Opt.ctrl_bam_file_name != 0) {
      cpics_output_to_file(cpics_exp_out_file, cpics_exp_out_file_name);
      process_chr_timed(&data, exp_read_count, ctrl_read_count);

      fprintf(stderr, "Swapping exp and control to compute FDR\n");
      swap_exp_ctrl(&data);

      cpics_output_to_file(cpics_ctrl_out_file, cpics_ctrl_out_file_name);
      process_chr_timed(&data, ctrl_read_count, exp_read_count);
    } else {
      process_chr_timed(&data, exp_read_count, ctrl_read_count);
    }
    free_input_data(&data);
  }

  if( Opt.unmappable_file_name != 0 && bed_reader_has_more_sites(&bed_reader) ) {
    fprintf(stderr, "Distinct set of chromosomes or wrong order of chromosomes in unmappable bed file\n");
    return 1;
  } 

  close_bam_reader(&exp_bam_reader);
  if(Opt.ctrl_bam_file_name != 0)
    close_bam_reader(&ctrl_bam_reader);
  if(Opt.unmappable_file_name != 0)
    close_bed_reader(&bed_reader);

  if(Opt.ctrl_bam_file_name != 0) {
    close_tmp_files();
    fprintf(stderr, "Computing FDR\n");
    compute_fdr();
  }

  if(Opt.chr_file_name)
    free_chromosomes(chrs, n_chrs);

  close_out_file();
  fprintf(stderr, "Done\n");

  return 0;
}
