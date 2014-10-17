#include "numerics.h"

float **allocate_float_matrix(int nrows, int ncols)
{
	float ** matrix;
	matrix = new float*[nrows];
	for(int i=0; i < nrows; i++) matrix[i] = new float[ncols];
	return matrix;
}


void desallocate_float_matrix(float **matrix, int nrows, int /*ncols*/)
{
	if (matrix == NULL) return;

	for(int i=0; i < nrows; i++) { delete[] matrix[i]; matrix[i] = 0;}
    delete [] matrix;
	matrix = 0;
}



void print_float_matrix(float **matrix, int nrows, int ncols)
{

	for(int i=0; i < nrows; i++, printf("\n"))
	  for(int j=0; j < ncols; j++)	printf("%f  ", matrix[i][j]);

	printf("\n");
}




void product_square_float_matrixes(float **result,float **matrix1,float **matrix2,int n)
{

	for(int j=0; j < n ; j++)
		for(int i=0; i < n; i++)
		{

			float aux = 0.0;
			for(int k=0; k < n; k++) aux += matrix1[i][k] * matrix2[k][j];
			result[i][j] = aux;

		}

}


void float_vector_matrix_product(float **a, float *x,float *y,int n)
{

  for(int i=0;i<n;i++){
    float sum=0.0;
    for(int j=0;j<n;j++)
      sum += a[i][j]*x[j];

    y[i]=sum;
  }


}


// **********************************************
//  LU based algorithms
// **********************************************



// Solves Ax=b by using lu decomposition 
int lusolve(float **a, float *x, float *b, int n)
{

	float d;
	int *indx = new int[n];

	if (ludcmp(a,n,indx,&d))
	{

		for(int i=0; i < n; i++) x[i] = b[i];
		lubksb(a,n,indx,x);
        delete [] indx;
		
		return 1;
	} else
	{
        delete [] indx;
		printf("lusolve::lu decomposition failed\n");
		return 0;
	}

}


// Computes the inverse of a column by column using the LU decomposition
int luinv(float **a, float **inv, int n)
{
	float d;

	float *col = new float[n];
	int *indx = new int[n];

	float **invaux = allocate_float_matrix(n,n);
	for(int i=0; i <n;i++)
		for(int j=0; j < n; j++)
			invaux[i][j] = a[i][j];


	if (ludcmp(invaux,n,indx,&d))
	{
		for(int j=0;j<n;j++){
    
			for(int i=0;i<n;i++) col[i]=0.0;
			col[j]=1.0;

			lubksb(invaux,n,indx,col);

			for(int i=0;i<n;i++) inv[i][j]=col[i];
		}
			
		return 1;
			
	} else
	{
		printf("luinv::lu decomposition failed\n");
		return 0;		
	}

	delete[] col;
	delete[] indx;
	
}




int ludcmp(float **a, int n, int *indx, float *d)
{

	int i,imax=-1,j,k,aux;
	float big,dum,sum,temp;
	float *vv;


	vv=(float *) malloc(n*sizeof(float));
	*d=1.0;


	for(i=0;i<n;i++) indx[i]=i;


	/**** Look for the largest value of every line and store 1/this value****/
	for(i=0;i<n;i++){

		big=0.0;
		for(j=0;j<n;j++) 
			if ( (temp=fabs(a[i][j]))>big ) big=temp;

			if (big==0.0) { return 0; printf("LU Decomposition failed\n");}
			
		vv[i]=1.0/big;
	}


  

	for(j=0;j<n;j++){

		for(i=0;i<j;i++){
      
			sum=a[i][j];
			for(k=0;k<i;k++) sum-= a[i][k]*a[k][j];
			a[i][j]=sum; 

		}


		big=0.0;
		for(i=j;i<n;i++){
      
			sum=a[i][j];
			for(k=0;k<j;k++)   
				sum -=a[i][k]*a[k][j];
			a[i][j]=sum;

			if ( (dum=vv[i]*fabs(sum))>=big){
				big=dum;
				imax=i;
			}
		}

    
		if (j != imax){

			for(k=0;k<n;k++){

				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}

			*d=-(*d);
			vv[imax]=vv[j];

			aux=indx[j];
			indx[j]=indx[imax];
			indx[imax]=aux;
		}


 


		if (a[j][j]==0.0)  a[j][j]=NRTINY;
    

		if (j!=n-1 ){

			dum=1.0 / a[j][j];
			for(i=j+1;i<n;i++) a[i][j]*=dum;
		}
    

	}

	free(vv);
	return 1;
}



/* Solves the set of n linear equations Ax=b. Here a[0..n-1][0..n-1] as input, not as the matrix A but rather as its LU decomposition,*/
/* determined by the routine ludcmp. indx[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is input as the    */
/* right hand side vector and returns with the solution vector x. */
void lubksb(float **a, int n, int *indx, float *b)
{

	int i,ii=0,j;
	float sum,*aux;

	aux= (float *) malloc(n*sizeof(float));
	
	
	for(i=0;i<n;i++)
		aux[i]=b[indx[i]];


	for(i=0;i<n;i++){
		sum=aux[i];

		if (ii)
			for(j=ii-1;j<i;j++) sum-=a[i][j]*aux[j];
		else if (sum) ii=i+1;
		aux[i]=sum;
	}

	for(i=n-1;i>=0;i--){
		sum=aux[i];
		for(j=i+1;j<n;j++) sum-=a[i][j]*aux[j];
		aux[i]=sum/a[i][i];
		b[i]=sum/a[i][i];
	}

	free(aux);
}





// **********************************************
//  Householder reduction to tridiagonal matrix
// **********************************************



void symmetric_2_tridiag_HH(float **a,float *d,float *e,int n)
{
	int l,k,j,i;
	float scale,hh,h,g,f;

	
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}

	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}



// **********************************************
//  QL decomposition algorithm to be used with tridiagonal matrices
// **********************************************

/*- Computes a*a -*/
float SQR(float a)
{
	return a*a;
}


/* returns a or -a in function of the signe of a*b */
float SIGN(float a,float b)
{
	if (b>=0)
		return  (a>=0 ? a : -a) ;
	else 
		return  (a>=0 ? -a : a) ;

}


/*- Computes (a*a+b*b)^1/2 without destructive underflow or overflow -*/
float pythag(float a,float b)
{
	float absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa>absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb==0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


int eigenvalues_QLI(float * d, float *e,float **z, int n)
{
	int m,l,iter,i,k;
	float s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) { printf("Too many iterations in eigenvalues_QLI\n"); return 0;}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}

	return 1;
}




// **********************************************
//  Singular Value Decomposition
// **********************************************


float withSignOf(float a, float b)
{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

float svdhypot(float a, float b)
{
    a = fabsf(a);
    b = fabsf(b);
    if(a > b) {
        b /= a;
        return a*sqrt(1.0 + b*b);
    } else if(b) {
        a /= b;
        return b*sqrt(1.0 + a*a);
    }
    return 0.0;
}

void svdrotate(float& a, float& b, float c, float s)
{
    float d = a;
    a = +d*c +b*s;
    b = -d*s +b*c;
}


//  m_U(A), m_V(A.ncol(),A.ncol()), m_W(A.ncol())
void compute_svd(float **A, float **m_U, float **m_V, float *m_W, int rows, int cols)
{
    const float	EPSILON = 0.00001;
    const int SVD_MAX_ITS = 100;

    float g, scale, anorm;
    float * RV1 = new float[cols];


    for(int i=0; i < rows; i++)
      for(int j=0; j < cols; j++)
		m_U[i][j] = A[i][j];		

    // Householder reduction to bidiagonal form:
    anorm = g = scale = 0.0;
    for (int i=0; i< cols; i++) {
        int l = i + 1;
        RV1[i] = scale*g;
        g = scale = 0.0;
        if(i< rows) {
            for (int k=i; k< rows; k++)
                scale += fabsf(m_U[k][i]);
            if (scale != 0.0) {
                float invScale=1.0/scale, s=0.0;
                for (int k=i; k< rows; k++) {
                    m_U[k][i] *= invScale;
                    s += m_U[k][i] * m_U[k][i];
                }
                float f = m_U[i][i];
                g = - withSignOf(sqrt(s),f);
                float h = 1.0 / (f*g - s);
                m_U[i][i] = f - g;
                for (int j=l; j< cols; j++) {
                    s = 0.0;
                    for (int k=i; k< rows; k++)
                        s += m_U[k][i] * m_U[k][j];
                    f = s * h;
                    for (int k=i; k< rows; k++)
                        m_U[k][j] += f * m_U[k][i];
                }
                for (int k=i; k< rows; k++)
                    m_U[k][i] *= scale;
            }
        }
        m_W[i] = scale * g;
        g = scale = 0.0;
        if ( i< rows && i< cols-1 ) {
            for (int k=l; k< cols; k++)
                scale += fabsf(m_U[i][k]);
            if (scale != 0.0) {
                float invScale=1.0/scale, s=0.0;
                for (int k=l; k< cols; k++) {
                    m_U[i][k] *= invScale;
                    s += m_U[i][k] * m_U[i][k];
                }
                float f = m_U[i][l];
                g = - withSignOf(sqrt(s),f);
                float h = 1.0 / (f*g - s);
                m_U[i][l] = f - g;
                for (int k=l; k< cols; k++)
                    RV1[k] = m_U[i][k] * h;
                for (int j=l; j< rows; j++) {
                    s = 0.0;
                    for (int k=l; k< cols; k++)
                        s += m_U[j][k] * m_U[i][k];
                    for (int k=l; k< cols; k++)
                        m_U[j][k] += s * RV1[k];
                }
                for (int k=l; k< cols; k++)
                    m_U[i][k] *= scale;
            }
        }
        anorm = NRMAX(anorm, fabsf(m_W[i]) + fabsf(RV1[i]) );
    }

    // Accumulation of right-hand transformations:
    m_V[cols-1][cols-1] = 1.0;
    for (int i= cols-2; i>=0; i--) {
        m_V[i][i] = 1.0;
        int l = i+1;
        g = RV1[l];
        if (g != 0.0) {
            float invgUil = 1.0 / (m_U[i][l]*g);
            for (int j=l; j< cols; j++)
                m_V[j][i] = m_U[i][j] * invgUil;
            for (int j=l; j< cols; j++){
                float s = 0.0;
                for (int k=l; k< cols; k++)
                    s += m_U[i][k] * m_V[k][j];
                for (int k=l; k< cols; k++)
                    m_V[k][j] += s * m_V[k][i];
            }
        }
        for (int j=l; j< cols; j++)
            m_V[i][j] = m_V[j][i] = 0.0;
    }

    // Accumulation of left-hand transformations:
    for (int i=NRMIN(rows,cols)-1; i>=0; i--) {
        int l = i+1;
        g = m_W[i];
        for (int j=l; j< cols; j++)
            m_U[i][j] = 0.0;
        if (g != 0.0) {
            g = 1.0 / g;
            float invUii = 1.0 / m_U[i][i];
            for (int j=l; j< cols; j++) {
                float s = 0.0;
                for (int k=l; k< rows; k++)
                    s += m_U[k][i] * m_U[k][j];
                float f = (s * invUii) * g;
                for (int k=i; k< rows; k++)
                    m_U[k][j] += f * m_U[k][i];
            }
            for (int j=i; j< rows; j++)
                m_U[j][i] *= g;
        } else
            for (int j=i; j< rows; j++)
                m_U[j][i] = 0.0;
        m_U[i][i] = m_U[i][i] + 1.0;
    }

    // Diagonalization of the bidiagonal form:
    for (int k=cols-1; k>=0; k--) { // Loop over singular values
        for (int its=1; its<=SVD_MAX_ITS; its++) {
            bool flag = false;
            int l  = k;
            int nm = k-1;
            while(l>0 && fabsf(RV1[l]) > EPSILON*anorm) { // Test for splitting
                if(fabsf(m_W[nm]) <= EPSILON*anorm) {
                    flag = true;
                    break;
                }
                l--;
                nm--;
            }
            if (flag) {	// Cancellation of RV1[l], if l > 0
                float c=0.0, s=1.0;
                for (int i=l; i< k+1; i++) {
                    float f = s * RV1[i];
                    RV1[i] = c * RV1[i];
                    if (fabsf(f)<=EPSILON*anorm)
                        break;
                    g = m_W[i];
                    float h = svdhypot(f,g);
                    m_W[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = - f * h;
                    for (int j=0; j< rows; j++)
                        svdrotate(m_U[j][nm],m_U[j][i], c,s); 
                }
            }
            float z = m_W[k];
            if (l==k) {		// Convergence of the singular value
                if (z< 0.0) {	// Singular value is made nonnegative
                    m_W[k] = -z;
                    for (int j=0; j< cols; j++)
                        m_V[j][k] = - m_V[j][k];
                }
                break;
            }
            // Exception if convergence to the singular value not reached:
            if(its==SVD_MAX_ITS) {printf("svd::convergence_error\n"); exit(-1);}
            float x = m_W[l]; // Get QR shift value from bottom 2x2 minor
            nm = k-1;
            float y = m_W[nm];
            g = RV1[nm];
            float h = RV1[k];
            float f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0*h*y );
            g = svdhypot(f,1.0);
            f = ( (x-z)*(x+z) + h*(y/(f+withSignOf(g,f)) - h) ) / x;
            // Next QR transformation (through Givens reflections)
            float c=1.0, s=1.0;
            for (int j=l; j<=nm; j++) {
                int i = j+1;
                g = RV1[i];
                y = m_W[i];
                h = s * g;
                g = c * g;
                z = svdhypot(f,h);
                RV1[j] = z;
                z = 1.0 / z;
                c = f * z;
                s = h * z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y * s;
                y *= c;
                for(int jj=0; jj < cols; jj++)
                    svdrotate(m_V[jj][j],m_V[jj][i], c,s);
                z = svdhypot(f,h);
                m_W[j] = z;
                if (z!=0.0) { // Rotation can be arbitrary if z = 0.0
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for(int jj=0; jj < rows; jj++)
                    svdrotate(m_U[jj][j],m_U[jj][i], c,s);
            }
            RV1[l] = 0.0;
            RV1[k] = f;
            m_W[k] = x;
        }
    }
}




// **********************************************
//  PCA
// **********************************************


void  pca_center_data(float **X,float *baricenter,int n, int p)
{


	for(int j=0;j<p;j++){

		float some=0.0;
		for(int i=0; i<n; i++) some += X[i][j];

		some/=(float) n;
		baricenter[j]=some;

		for(int i=0; i<n; i++) X[i][j]-=baricenter[j];
	}
	
}





void order_decreasing(float *values, int *indexos, int size)
{

	for(int i=0; i < size; i++)
	{
		float dmax = values[i];
		int pos = i;
		
		for(int j=i+1; j < size; j++)
			if (values[j] > dmax) {dmax = values[j]; pos = j;}


		values[pos] = values[i];
		values[i] = dmax;
		
		int aux = indexos[i];

		indexos[i] = indexos[pos];
		indexos[pos] = aux;


	}

}




// U = X * V, where X is centered
// so   X = U * V^{-1}

void compute_pca_svd(float **X,float *S,float **V, float **U,int n, int p)
{

	// Compute PCA
	// U: nxp   V: pxp   S: p
	compute_svd(X,U,V,S,n,p);


	// U contain new coefficients, which must be normalized by eigenvalues
	for(int i=0; i < n; i++)
		for(int j=0; j < p; j++)
			U[i][j] *= S[j];

	// Normalize eigenvalues
	float norm = (float) (n-1);
	for(int i=0; i < p; i++)
			S[i] = S[i] * S[i] / norm;

	// If n < p, principal component should be zero from n to p-1			
	// Coefficients of these principal components should be zero
	if (n < p)
	{
		for(int i=n-1; i < p; i++) S[i] = 0.0f;

		for(int j=0; j < n; j++)
			for(int i=n-1; i < p; i++)
				U[j][i] = 0.0f;
	}


}







void compute_pca_brute(float **X,float *S,float **V, float **U,int n, int p)
{

	compute_pca(X,V,S,n,p);

	compute_coefic(X,V,p,U,n,p);

}


float **  covariance_matrix(float **x,int n, int p)
{

	int i,j,k;
	float **cov;
	float some;

	cov = allocate_float_matrix(p,p);

	/*--- Calculam la matriu de covariance cov[k][j] ---*/
	for(j=0;j<p;j++)
		for(k=0;k<=j;k++){
		some=0.0;
		for(i=0;i<n;i++)
			some+=(x[i][j])*(x[i][k]);
		some/=(float) (n-1);
		cov[k][j]=cov[j][k]=some;
		}


	return cov;
		
}


int compute_pca(float **X,float **Pcs,float *sVar,int n,int p)
{
  
	float **Mat_cov;
	float  *ve, *vpdes, *posicio;
	int i,j,pos,k,l;
	float dmax,*vap;

 
	Mat_cov=covariance_matrix(X,n,p);
  
	vpdes= (float *) malloc(p*sizeof(float));
	ve = (float *) malloc(p*sizeof(float));;
	


	symmetric_2_tridiag_HH(Mat_cov,vpdes,ve,p);

	if (eigenvalues_QLI(vpdes,ve,Mat_cov,p))
	{
	



		/*--------- Ordenam valors propis  -----------------*/
		vap=vpdes;
		posicio = (float *) malloc(p*sizeof(float));

		for(i=0;i<p;i++){
			
			dmax=vap[0];
			pos=0;
		
			for(j=1;j<p;j++)
				if (vap[j]>dmax){
				dmax=vap[j];
				pos=j;
				}

			sVar[i]=dmax;
			posicio[i]=(float) pos;
			vap[pos]=0.0;
	
		}


		/*---------------Ordenam vectors propis------------*/

 
		for(k=0;k<p;k++)
			for(l=0;l<p;l++){
		
			Pcs[l][k]=Mat_cov[l][(int) posicio[k]];
		    }


	
		desallocate_float_matrix(Mat_cov,p,p);	
		free(ve);
		free(posicio);
		free(vpdes);
		
		
		return 1;
		
	} else
	{
	
		desallocate_float_matrix(Mat_cov,p,p);	
		free(ve);
		free(vpdes);
		
		printf("PCA failed.");
		return 0;
	
	}
	
}


void compute_coefic(float **x,float **Pcs,int nvec,float **Coef, int n, int p)
{

	int i,j,k;
	float suma;
  
	for(i=0;i<n;i++)
		for(j=0;j<nvec;j++){

		suma=0.0;
		for(k=0;k<p;k++)
			suma+=x[i][k]*Pcs[k][j];

		Coef[i][j]=suma;
		}

}

