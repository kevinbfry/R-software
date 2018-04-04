#include <math.h> // for fabs

// Find an approximate row of \hat{Sigma}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \theta^T \Sigma \theta - e_i^T\theta + \mu \|\theta\|_1

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

double objective(double *Sigma_ptr,       /* A covariance matrix: X^TX/n */
		 int *ever_active_ptr,    /* Ever active set: 1-based */ 
		 int *nactive_ptr,        /* Size of ever active set */
		 int nrow,                /* how many rows in Sigma */
		 int row,                 /* which row: 1-based */
		 double bound,            /* Lagrange multipler for \ell_1 */
		 double *theta)           /* current value */
{
  int irow, icol;
  double value = 0;
  double *Sigma_ptr_tmp = Sigma_ptr;
  double *theta_row_ptr, *theta_col_ptr;
  int *active_row_ptr, *active_col_ptr;
  int active_row, active_col;
  int nactive = *nactive_ptr;

  theta_row_ptr = theta;
  theta_col_ptr = theta;

  for (irow=0; irow<nactive; irow++) {

    active_row_ptr = ((int *) ever_active_ptr + irow);
    active_row = *active_row_ptr - 1;    // Ever-active set is 1-based
    theta_row_ptr = ((double *) theta + active_row);

    for (icol=0; icol<nactive; icol++) {
      
      active_col_ptr = ((int *) ever_active_ptr + icol);
      active_col = *active_col_ptr - 1;    // Ever-active set is 1-based
      theta_col_ptr = ((double *) theta + active_col);

      Sigma_ptr_tmp = ((double *) Sigma_ptr + nrow * active_col + active_row); // Matrices are column-major order

      value += 0.5 * (*Sigma_ptr_tmp) * (*theta_row_ptr) * (*theta_col_ptr);
    }
    value = value + bound * fabs((*theta_row_ptr)); // the \ell_1 term
  }

  theta_row_ptr = ((double *) theta + (row - 1));     // row is 1-based index
  value -= (*theta_row_ptr); // the elementary basis vector term

  return(value);
}

// Check if active and add it to active list if necessary
// Ever active set is stored as 1-based indices
// coord is 0-based
int update_ever_active(int coord,
		       int *ever_active_ptr,
		       int *nactive_ptr) {
  int iactive;
  int active_var;
  int nactive = *nactive_ptr;
  int *ever_active_ptr_tmp = ever_active_ptr;

  for (iactive=0; iactive<nactive; iactive++) {
    ever_active_ptr_tmp = ((int *) ever_active_ptr + iactive);
    active_var = (*ever_active_ptr_tmp);
    if (active_var - 1 == coord) {       // Indices in active set are 1-based
      return(1);
    }
  }

  // If we have not returned yet, this variable
  // was not in ever_active

  // Add it to the active set and increment the 
  // number of active variables

  ever_active_ptr_tmp = ((int *) ever_active_ptr + *nactive_ptr);
  *ever_active_ptr_tmp = coord + 1;     // Indices are 1-based
  *nactive_ptr += 1;

  return(0);
}

int check_KKT(double *theta,        /* current theta */
	      double *gradient_ptr, /* Sigma times theta */
	      int nrow,             /* how many rows in Sigma */
	      int row,              /* which row: 1-based */
	      double bound,         /* Lagrange multipler for \ell_1 */
	      double tol)           /* precision for checking KKT conditions */
{
  // First check inactive

  int irow;
  double *theta_ptr, *gradient_ptr_tmp;
  double gradient;

  for (irow=0; irow<nrow; irow++) {
    theta_ptr = ((double *) theta + irow);
    gradient_ptr_tmp = ((double *) gradient_ptr + irow);

    // Compute this coordinate of the gradient

    gradient = *gradient_ptr_tmp;

    // For the basis vector

    if (row - 1 == irow) {                // Row is a 1-based index
       gradient -= 1;
    }

    if (*theta_ptr != 0) { // these coordinates of gradients should be equal to -bound
      if ((*theta_ptr > 0) &&  (fabs(gradient + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr < 0) && (fabs(gradient - bound) > tol * bound)) {
	return(0);
      }
    }
    else {
      if (fabs(gradient) > (1. + tol) * bound) {
	return(0);
      }
    }

  }

  return(1);
}

double update_one_coord(double *Sigma_ptr,           /* A covariance matrix: X^TX/n */
                        double *Sigma_diag_ptr,      /* Diagonal entries of Sigma */
                        double *gradient_ptr,        /* Sigma times theta */
			int *ever_active_ptr,        /* Ever active set: 0-based */ 
			int *nactive_ptr,            /* Size of ever active set */
			int nrow,                    /* How many rows in Sigma */
			double bound,                /* feasibility parameter */
			double *theta,               /* current value */
			int row,                     /* which row: 1-based */
			int coord,                   /* which coordinate to update: 0-based */
			int is_active)               /* Is this part of ever_active */     
{

  double delta;
  double linear_term = 0;
  double value = 0;
  double old_value;
  double *Sigma_ptr_tmp;
  double *gradient_ptr_tmp;
  double *theta_ptr;
  int icol = 0;

  double *quadratic_ptr = ((double *) Sigma_diag_ptr + coord);
  double quadratic_term = *quadratic_ptr;

  gradient_ptr_tmp = ((double *) gradient_ptr + coord);
  linear_term = *gradient_ptr_tmp;

  theta_ptr = ((double *) theta + coord);
  old_value = *theta_ptr;

  // The coord entry of gradient_ptr term has a diagonal term in it:
  // Sigma[coord, coord] * theta[coord]
  // This removes it. 

  linear_term -= quadratic_term * old_value;

  if (row - 1 == coord) {              // Row is 1-based
    linear_term -= 1;
  }

  // Now soft-threshold the coord entry of theta 

  // Objective is t \mapsto q/2 * t^2 + l * t + bound |t|
  // with q=quadratic_term and l=linear_term

  // With a negative linear term, solution should be
  // positive

  if (linear_term < -bound) {
    value = (-linear_term - bound) / quadratic_term;
  }
  else if (linear_term > bound) {
    value = -(linear_term - bound) / quadratic_term;
  }

  // Add to active set if necessary

  if ((is_active == 0) && (value != 0)) {
    update_ever_active(coord, ever_active_ptr, nactive_ptr);
  }

  // Update the linear term

  if (fabs(old_value - value) > 1.e-6 * (fabs(value) + fabs(old_value))) { 

    delta = value - old_value;
    Sigma_ptr_tmp = ((double *) Sigma_ptr + coord * nrow);
    gradient_ptr_tmp = ((double *) gradient_ptr);

    for (icol=0; icol<nrow; icol++) {
      *gradient_ptr_tmp = *gradient_ptr_tmp + delta * (*Sigma_ptr_tmp);
      gradient_ptr_tmp += 1;
      Sigma_ptr_tmp += 1;
    }

    theta_ptr = ((double *) theta + coord);
    *theta_ptr = value;

  }

  return(value);

}

int find_one_row_(double *Sigma_ptr,          /* A covariance matrix: X^TX/n */
		  double *Sigma_diag_ptr,     /* Diagonal entry of covariance matrix */
		  double *gradient_ptr,       /* Sigma times theta */
		  int *ever_active_ptr,       /* Ever active set: 1-based */ 
		  int *nactive_ptr,           /* Size of ever active set */
		  int nrow,                   /* How many rows in Sigma */
		  double bound,               /* feasibility parameter */
		  double *theta,              /* current value */
		  int maxiter,                /* how many iterations */
		  int row,                    /* which coordinate to solve: 1-based */
		  double kkt_tol,             /* precision for checking KKT conditions */
		  double objective_tol)       /* precision for checking relative decrease in objective value */
{

  int iter = 0;
  int icoord = 0;
  int iactive = 0;
  int *active_ptr;
  double old_value, new_value;

  int check_objective = 1;

  if (check_objective) {

    old_value = objective(Sigma_ptr,
			  ever_active_ptr,
			  nactive_ptr,
			  nrow,
			  row,
			  bound,
			  theta);
  }


  for (iter=0; iter<maxiter; iter++) {

    // Update the active variables first

    active_ptr = (int *) ever_active_ptr;

    for (iactive=0; iactive < *nactive_ptr; iactive++) {
      update_one_coord(Sigma_ptr,
		       Sigma_diag_ptr,
		       gradient_ptr,
		       ever_active_ptr,
		       nactive_ptr,
		       nrow,
		       bound,
		       theta,
		       row,
		       *active_ptr-1,  // Ever active set is 1-based
		       1);
      active_ptr++;
    }

    // Check KKT

    if (check_KKT(theta,
		  gradient_ptr,
		  nrow,
		  row,
		  bound,
		  kkt_tol) == 1) {
      break;
    }
					  
    // Update all variables 

    for (icoord=0; icoord<nrow; icoord++) {

      update_one_coord(Sigma_ptr,
		       Sigma_diag_ptr,
		       gradient_ptr,
		       ever_active_ptr,
		       nactive_ptr,
		       nrow,
		       bound,
		       theta,
		       row,
		       icoord,
		       0);
    }

    // Check KKT

    if (check_KKT(theta,
		  gradient_ptr,
		  nrow,
		  row,
		  bound,
		  kkt_tol) == 1) {
      break;
    }
					  
    if (check_objective) {
      new_value = objective(Sigma_ptr,
			    ever_active_ptr,
			    nactive_ptr,
			    nrow,
			    row,
			    bound,
			    theta);

      if (((old_value - new_value) < objective_tol * fabs(new_value)) && (iter > 0)) {
	break;
      }

      old_value = new_value;
    }
  }
  return(iter);
}

