#define ANSI_COLOR_RED     "\x1b[31;1m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_BLUE    "\x1b[34;1m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define TEST(func) \
if(!func) { \
   printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [ " ANSI_COLOR_GREEN "ok" ANSI_COLOR_RESET " ]\n", #func); \
} else { \
   printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [" ANSI_COLOR_RED "fail" ANSI_COLOR_RESET "]\n", #func); \
   nerr++; \
}

void gen_posdef(int n, double* a, int lda);
void print_vec(char const* format, int n, double const* vec);
void print_mat(char const* format, int n, double const* a, int lda);
double backward_error(int n, double const* a, int lda, double const* rhs, int nrhs, double const* soln, int ldsoln);
double forward_error(int n, int nrhs, double const* soln, int ldx);
