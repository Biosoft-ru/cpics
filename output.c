
static FILE* cpics_out_file = 0;

static FILE* cpics_cur_out_file;
static const char* cpics_cur_out_file_name;

static void cpics_output_to_file(FILE* f, const char* name) {
  cpics_cur_out_file = f;
  cpics_cur_out_file_name = name;
}

static void open_out_file() {
  cpics_out_file = fopen(Opt.cpics_out_file_name, "w");
  if(!cpics_out_file) {
    perror(Opt.cpics_out_file_name);
    exit(1);
  }
  cpics_output_to_file(cpics_out_file, Opt.cpics_out_file_name);
}

static void close_out_file() {
  if(fclose(cpics_out_file) != 0) {
    perror(Opt.cpics_out_file_name);
    exit(1);
  }
}

static int score_filter_disabled = 0;
static void cpics_disable_score_filter() { score_filter_disabled = 1; }


static void output(struct MixtureComp* comps, int32_t nComp, char* chr) {
  int i;
  for(i = 0;i < nComp; i++) {
    struct MixtureComp c = comps[i];
    if(!score_filter_disabled && c.score <= Opt.score_cutoff) continue;
    if(c.delta <= 50 || c.delta >= 300) continue;
    if(c.se >= 50 || c.se <= 0) continue; //!!! ==0 ???
    if(c.seF >= 50 || c.seF <= 0) continue; //!!! ==0 ???
    if(c.seR >= 50 || c.seR <= 0) continue; //!!! ==0 ???
    if(c.sigmaSqF > 22500 || c.sigmaSqR > 22500) continue;
    int32_t from = (int32_t) (c.mu - c.delta/2 - 3*c.seF);
    int32_t to = (int32_t) (c.mu + c.delta/2 + 3*c.seR);
    if(fprintf(cpics_cur_out_file, "%s\t%i\t%i\t.\t%.17g\n", chr, from - 1, to, c.score) < 0) {
      perror(cpics_cur_out_file_name);
      exit(1);
    }
  }
}

static char* make_tmp_file_name(const char* suffix) {
  size_t l1 = strlen(Opt.cpics_out_file_name);
  size_t l2 = strlen(suffix);
  char* res = malloc( l1 + l2 + 1 );
  if(!res) {
    perror("in make_tmp_file_name");
    exit(1);
  }
  memcpy(res, Opt.cpics_out_file_name, l1);
  memcpy(res + l1, suffix, l2 + 1);
  return res;
}

static char *cpics_exp_out_file_name, *cpics_ctrl_out_file_name;
static FILE *cpics_exp_out_file, *cpics_ctrl_out_file;

static void open_tmp_files(const char* mode) {
  if(!cpics_exp_out_file_name)
    cpics_exp_out_file_name = make_tmp_file_name(".exp");
  cpics_exp_out_file = fopen(cpics_exp_out_file_name, mode);
  if(!cpics_exp_out_file) {
    perror(cpics_exp_out_file_name);
    exit(1);
  }

  if(!cpics_ctrl_out_file_name)
    cpics_ctrl_out_file_name = make_tmp_file_name(".ctrl");
  cpics_ctrl_out_file = fopen(cpics_ctrl_out_file_name, mode);
  if(!cpics_ctrl_out_file) {
    perror(cpics_ctrl_out_file_name);
    exit(1);
  }
}

static void close_tmp_files() {
  if(fclose(cpics_exp_out_file) != 0) {
    perror(cpics_exp_out_file_name);
    exit(1);
  }
  if(fclose(cpics_ctrl_out_file) != 0) {
    perror(cpics_ctrl_out_file_name);
    exit(1);
  }
}

static void remove_tmp_files() {
  if(remove(cpics_exp_out_file_name) != 0) {
    perror(cpics_exp_out_file_name);
    exit(1);
  }
  if(remove(cpics_ctrl_out_file_name) != 0) {
    perror(cpics_ctrl_out_file_name);
    exit(1);
  }
  free(cpics_exp_out_file_name);
  free(cpics_ctrl_out_file_name);
}
