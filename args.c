
static int step = 20;
static int width = 250;
static int min_reads = 2;
static int min_reads_in_region = 3;
static int min_l_region = 100;
static char *exp_bam_file_name = 0, *ctrl_bam_file_name = 0, *unmappable_file_name = 0;

static void usage() {
  fprintf(stderr, "Usage: cpics [-u <unmappable bed file>] [-c <control bam file>] <bam file> \n" );
}

static int parse_args(int argc, char* argv[]) {
  int i;
  for(i = 1; i < argc; i++) {
    char* val = argv[i];
    if(*val == '-')  
      switch(*(val+1)) {
      case 'u':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -u option\n"); usage();
          return 1;
        }
        unmappable_file_name = argv[++i];
        break;
      case 'c':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -c option\n"); usage();
          return 1;
        }
        ctrl_bam_file_name= argv[++i];
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", val); usage();
        return 1;
      }
    else {
      if(exp_bam_file_name != 0) {
        fprintf(stderr, "extra argument %s\n", val); usage();
        return 1;
      }
      exp_bam_file_name = val;
    }
  }
  if(exp_bam_file_name == 0) {
    usage();
    return 1;
  }
  return 0;
}
