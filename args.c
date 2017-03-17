
static struct Options {
  int step;
  int width;
  int min_reads_in_window;
  int min_reads_in_region;
  int min_l_region;
  double fdr_cutoff;
  double score_cutoff;
  char *exp_bam_file_name;
  char *ctrl_bam_file_name;
  char *unmappable_file_name;
  char* cpics_out_file_name;
  char* chr_file_name;
} Opt = {
  .step = 20,
  .width = 250,
  .min_reads_in_window = 2,
  .min_reads_in_region = 3,
  .min_l_region = 100,
  .fdr_cutoff = 0.001,
  .score_cutoff = 10,
};


static void usage() {
  fprintf(stderr, "Usage: cpics <options> [-u <unmappable bed file>] [-c <control bam file>] <bam file> <out file>\n" );
  fprintf(stderr, "  -f <fdr_cutoff>, default 0.001\n" );
  fprintf(stderr, "  -s <score_cutoff>, default 10\n" );
  fprintf(stderr, "  -g file with list of chromosomes, one chr per line, default to all chromosomes in exp bam file.\n" );
}

static int parse_args(int argc, char* argv[]) {
  int score_set = 0;
  int fdr_set = 0;
  int i;
  for(i = 1; i < argc; i++) {
    char* val = argv[i];
    if(*val == '-' && val[1] && !val[2])  
      switch(val[1]) {
      case 'u':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -u option\n"); usage();
          return 1;
        }
        Opt.unmappable_file_name = argv[++i];
        break;
      case 'c':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -c option\n"); usage();
          return 1;
        }
        Opt.ctrl_bam_file_name= argv[++i];
        break;
      case 'f':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -f option\n"); usage();
          return 1;
        }
        char* fdr_str = argv[++i];
        char* fdr_str_end;
        Opt.fdr_cutoff = strtod(fdr_str, &fdr_str_end);
        if(fdr_str_end == fdr_str || *fdr_str_end != 0 || Opt.fdr_cutoff < 0 || Opt.fdr_cutoff > 1) {
          fprintf(stderr, "Wrong -f value\n");
          return 1;
        }
        fdr_set = 1;
        break;
      case 's':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -s option\n"); usage();
          return 1;
        }
        char* score_str = argv[++i];
        char* score_str_end;
        Opt.score_cutoff = strtod(score_str, &score_str_end);
        if(score_str_end == score_str || *score_str_end != 0 || Opt.score_cutoff < 0) {
          fprintf(stderr, "Wrong -s value\n");
          return 1;
        }
        score_set = 1;
        break;
      case 'g':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -g option\n"); usage();
          return 1;
        }
        Opt.chr_file_name = argv[++i];
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", val); usage();
        return 1;
      }
    else {
      if(Opt.exp_bam_file_name != 0) {
        if(Opt.cpics_out_file_name != 0) {
          fprintf(stderr, "extra argument %s\n", val); usage();
          return 1;
        }
        Opt.cpics_out_file_name = val;
      } else
        Opt.exp_bam_file_name = val;
    }
  }
  if(Opt.exp_bam_file_name == 0 || Opt.cpics_out_file_name == 0) {
    usage();
    return 1;
  }
  if(fdr_set && Opt.ctrl_bam_file_name == 0) {
    fprintf(stderr, "-f can be used only when control present\n");
    return 1;
  }
  if(score_set && Opt.ctrl_bam_file_name != 0) {
    fprintf(stderr, "-s can be used only when without control\n");
    return 1;
  }
  return 0;
}
