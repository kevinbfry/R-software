#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

int solve_qp(double *Sigma_ptr,          /* A covariance matrix: X^TX/n */
	     double *linear_func_ptr,    /* Linear term in objective */
	     double *Sigma_diag_ptr,     /* Diagonal entry of covariance matrix */
	     double *gradient_ptr,       /* Current gradient of quadratic loss */
	     int *ever_active_ptr,       /* Ever active set: 0-based */ 
	     int *nactive_ptr,           /* Size of ever active set */
	     int nrow,                   /* How many rows in Sigma */
	     double bound,               /* feasibility parameter */
	     double *theta,              /* current value */
	     int maxiter,                /* how many iterations */
	     double kkt_tol,             /* precision for checking KKT conditions */
	     double objective_tol);      /* precision for checking relative decrease in objective value */


int check_KKT_qp(double *theta,           /* current theta */
		 double *gradient_ptr,    /* Current gradient of quadratic loss */
		 int nrow,                /* how many rows in Sigma */
		 double bound,            /* Lagrange multipler for \ell_1 */ 
		 double tol);             /* precision for checking KKT conditions */

int find_one_row_(double *Sigma_ptr,          /* A covariance matrix: X^TX/n */
		  double *Sigma_diag_ptr,      /* Diagonal entry of covariance matrix */
		  double *gradient_ptr,       /* Current gradient of quadratic loss */
		  int *ever_active_ptr,       /* Ever active set: 0-based */ 
		  int *nactive_ptr,           /* Size of ever active set */
		  int nrow,                   /* How many rows in Sigma */
		  double bound,               /* feasibility parameter */
		  double *theta,              /* current value */
		  int maxiter,                /* how many iterations */
		  int row,                    /* which coordinate to update: 0-based */
		  double kkt_tol,             /* precision for checking KKT conditions */
		  double objective_tol);      /* precision for checking relative decrease in objective value */

int check_KKT(double *theta,           /* current theta */
	      double *gradient_ptr,    /* Current gradient of quadratic loss */
	      int nrow,                /* how many rows in Sigma */
	      int row,                 /* which row: 0-based */
	      double bound,            /* Lagrange multipler for \ell_1 */
	      double kkt_tol);         /* precision for checking KKT conditions */

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
