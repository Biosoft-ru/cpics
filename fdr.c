static void read_score_from_bed_file(FILE* f, const char* f_name, double** scoreP, int32_t* nP) {
  size_t capacity = 1024;
  size_t size = 0;
  double* scores = malloc(sizeof(double)*capacity);
  double score;
  int ret;
  errno = 0;
  while((ret = fscanf(f, "%*s %*d %*d %*s %lf", &score)) == 1) {
    if(size == capacity) {
      capacity *= 2;
      scores = realloc(scores, sizeof(double)*capacity);
    }
    scores[size++] = score;
  }

  if(errno != 0) {
    perror(f_name);
    exit(1);
  }
  if(ret != EOF) {
    fprintf(stderr, "Can not parse line from bed file %s, %d\n", f_name, ret);
    exit(1);
  }

  *scoreP = scores;
  *nP = size;
}

//Calculates FDR as in original PICS, but more efficient
static double compute_score_cutoff(double* exp, int32_t exp_count, double* ctrl, int32_t ctrl_count) {
  if(ctrl_count == 0 || exp_count == 0)
    return -INFINITY;

  qsort(exp, exp_count, sizeof(double), &double_cmp);
  qsort(ctrl, ctrl_count, sizeof(double), &double_cmp);

  double min = exp[0];
  if(ctrl[0] < min)
    min = ctrl[0];

  double max = exp[exp_count-1];

  double prev_fdr = 1;
  double x = min;
  int32_t i = 0, j = 0;
  while(x <= max) {

    while(i < exp_count && exp[i] <= x)
      i++;
    while(j < ctrl_count && ctrl[j] <= x)
      j++;

    if(ctrl_count == j)//FDR reaches 0 at x
      return x;

    if(exp_count == i)//FDR didn't reach fdr_cutoff for any score_cutoff
      return INFINITY;

    double fdr = (double)(ctrl_count - j)/(exp_count - i);
    if(fdr > prev_fdr)
      fdr = prev_fdr;

    if(fdr <= Opt.fdr_cutoff)
      return x;

    double next_score = x;
    if(i < exp_count)
      next_score = exp[i];
    if(j < ctrl_count && ctrl[j] < next_score)
      next_score = ctrl[j];
    int32_t nsteps = (int32_t)((next_score - x)/0.1);
    if(nsteps == 0)
      nsteps = 1;
    x += nsteps * 0.1;

    prev_fdr = fdr;
  }

  return INFINITY;
}

static void compute_fdr() {
  open_tmp_files("r");

  double* exp_scores;
  int32_t exp_count;
  read_score_from_bed_file(cpics_exp_out_file, cpics_exp_out_file_name, &exp_scores, &exp_count);

  double* ctrl_scores;
  int32_t ctrl_count;
  read_score_from_bed_file(cpics_ctrl_out_file, cpics_ctrl_out_file_name, &ctrl_scores, &ctrl_count);


  double score_cutoff = compute_score_cutoff(exp_scores, exp_count, ctrl_scores, ctrl_count);

  free(exp_scores);
  free(ctrl_scores);

  rewind(cpics_exp_out_file);
  char* line = 0;
  size_t line_len;
  while(getline(&line, &line_len, cpics_exp_out_file) != -1) {
    double score;
    if(sscanf(line, "%*s %*d %*d %*s %lf", &score) != 1) {
      fprintf(stderr, "Wrong line: %s\n", line);
      exit(1);
    }
    if(score >= score_cutoff)
      if(fputs(line, cpics_out_file) == EOF) {
        fprintf(stderr, "Can not write to file: %s\n", Opt.cpics_out_file_name);
        exit(1);
      }
  }
  free(line);

  close_tmp_files();
  remove_tmp_files();
}
