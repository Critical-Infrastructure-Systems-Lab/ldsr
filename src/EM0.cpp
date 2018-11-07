#include <RcppArmadillo.h>

const double pi = 3.141592653589793238463;

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Implement Kalman smoothing
//'
//' This is different from Kalman_smoother in that here we assume perfect prediction.
//' This function is experimental, do not use.
//'
List Kalman_smoother0(arma::mat y, arma::mat u, arma::mat v, List theta) {

    // Model parameters
    mat A = theta["A"];
    mat B = theta["B"];
    mat C = theta["C"];
    mat D = theta["D"];
    mat Q = theta["Q"];
    mat R = theta["R"];
    mat x1 = theta["mu1"];
    mat V1 = theta["V1"];

    // Declare variables
    int T = y.n_cols; // Time series length
    mat Xp, Vp, Yp, Xu, Vu, Xs, Vs, Ys, Cov, I, J, K, Sigma;

    I.eye(1,1);

    // FORWARD PASS ==================

    // Prior t|t-1, p stands for prior
    Xp.zeros(1, T);
    Vp.zeros(1, T);
    Yp.zeros(1, T);

    // Remember C++ starts from 0
    Xp.col(0) = x1;
    Vp.col(0) = V1;
    Yp.col(0) = C*Xp.col(0) + D*v.col(0);

    // Posterior t|t, u stands for updated
    Xu = Xp;
    Vu = Vp;

    // First time step
    K = Vp.col(0)*C / (C*Vp.col(0)*C + R);
    if (NumericMatrix::is_na(y(0,0))) {
        Xu.col(0) = Xp.col(0);
    } else {
        Xu.col(0) = Xp.col(0) + K*(y.col(0) - Yp.col(0));
    }
    Vu.col(0) = (I-K*C)*Vp.col(0);

    for (int t=1; t<T; t++) {
        Xp.col(t) = A*Xu.col(t-1) + B*u.col(t-1);
        Vp.col(t) = A*Vu.col(t-1)*A + Q;
        Yp.col(t) = C*Xp.col(t) + D*v.col(t);

        K = Vp.col(t)*C / (C*Vp.col(t)*C + R);

        if (NumericMatrix::is_na(y(0,t))) {
            Xu.col(t) = Xp.col(t);
        } else {
            Xu.col(t) = Xp.col(t) + K*(y.col(t) - Yp.col(t));
        }
        Vu.col(t) = (I-K*C)*Vp.col(t);
    }

    // BACKWARD PASS =========================
    // Smoothing, s stands for smoothed
    Xs = Xu;
    Vs = Vu;
    Cov = Vs; // Cov(X_t+1, X_t)

    for (int t=T-2; t>=0; t--) {
        J = Vu.col(t)*A * inv(Vp.col(t+1));
        Xs.col(t) = Xu.col(t) + J*(Xs.col(t+1) - Xp.col(t+1));
        Vs.col(t) = Vu.col(t) + J*(Vs.col(t+1) - Vp.col(t+1))*J;
        Cov.col(t) = Vs.col(t+1)*J;
    }

    // Final prediction
    Ys = C*Xs + D*v;

    // Likelihood
    uvec obs = find_finite(y);  // Find indices where there are observations
    int n_obs = obs.size();     // Number of observations
    mat delta = y.cols(obs) - Ys.cols(obs);  // Innovations

    Sigma = Vs.cols(obs);
    for (int i=0; i < n_obs; i++) {
        Sigma.col(i) = Sigma.col(i)*C*C + R;
    }

    double lik = -0.5*n_obs*log(2*pi) - 0.5*accu(delta / Sigma % delta + log(Sigma));
    lik = lik / n_obs;

    return List::create(Named("X") = Xs,
                        Named("Y") = Ys,
                        Named("V") = Vs,
                        Named("Cov") = Cov,
                        Named("lik") = lik);
}

//' Maximizing expected likelihood using analytical solution
//'
//' This is different from Mstep in that here we assume perfect prediction.
//' This function is experimental, do not use.
//'
List Mstep0(arma::mat y, arma::mat u, arma::mat v, List fit) {

    mat X = fit["X"];
    mat V = fit["V"];
    mat Cov = fit["Cov"];
    int T = y.n_cols;
    int d = u.n_rows;

    mat y2 = y;
    mat y_pred = fit["Y"];
    uvec miss = find_nonfinite(y); // indices of missing data
    y2.cols(miss) = y_pred.cols(miss);

    // C and D ================================
    // Sum for observed part
    mat Syx = y2 * trans(X);
    mat Syv = y2 * trans(v);
    mat Sxx = X * trans(X) + accu(V);
    mat Sxv = X * trans(v);
    mat Svx = trans(Sxv);
    mat Svv = v * trans(v);

    // Solve linear system
    mat P1 = join_horiz(Syx, Syv);
    mat P2 = join_vert(join_horiz(Sxx, Sxv), join_horiz(Svx, Svv));
    mat P  = P1 * inv(P2);

    mat C = P.col(0);
    mat D = P.cols(1, d);

    // A and B ================================
    mat Tx1x = X.cols(1, T-1) * trans(X.cols(0, T-2)) + accu(Cov.cols(0, T-2));
    mat Tx1u = X.cols(1, T-1) * trans(u.cols(0, T-2));
    mat Txx  = X.cols(0, T-2) * trans(X.cols(0, T-2)) + accu(V.cols(0, T-2));
    mat Tux  = u.cols(0, T-2) * trans(X.cols(0, T-2));
    mat Txu  = trans(Tux);
    mat Tuu  = u.cols(0, T-2) * trans(u.cols(0, T-2));

    P1 = join_horiz(Tx1x, Tx1u);
    P2 = join_vert(join_horiz(Txx, Txu), join_horiz(Tux, Tuu));
    P  = P1 * inv(P2);

    mat A = P.col(0);
    mat B = P.cols(1,d);

    // Q ===================================
    mat Tux1 = trans(Tx1u);
    mat Tx1x1 = X.cols(1,T-1) * trans(X.cols(1,T-1)) + accu(V.cols(1,T-1));

    mat Q = (Tx1x1 - A*Txx*A - Tx1u*trans(B) - B*Tux1 + B*Tuu*trans(B)) / (T-1);

    // R ====================================
    mat y_hat = C*X + D*v;
    mat delta_y = y2 - y_hat;

    mat R = (delta_y * trans(delta_y) + C*accu(V)*C) / T;

    // Initial state =========================
    mat mu1 = X.col(0);
    mat V1 = V.col(0);

    return List::create(Named("A") = A,
                        Named("B") = B,
                        Named("C") = C,
                        Named("D") = D,
                        Named("Q") = Q,
                        Named("R") = R,
                        Named("mu1") = mu1,
                        Named("V1") = V1);
}

//' Learn LDS model
//'
//' Estimate the hidden state and model parameters given observations and exogeneous inputs using the EM algorithm. This is the key backend routine of this package.
//'
//' This is different from learnLDS in that here we assume perfect prediction.
//' This function is experimental, do not use.
//'
// [[Rcpp::export]]
List learnLDS0(arma::mat y, arma::mat u, arma::mat v, arma::vec init, int niter, double tol) {

    int d = u.n_rows;
    mat A(1, 1); A.fill(init[0]);
    mat B(1, d); B.fill(init[1]);
    mat C(1, 1); C.fill(init[2]);
    mat D(1, d); D.fill(init[3]);
    mat Q(1, 1, fill::ones);
    mat R(1, 1, fill::ones);
    mat mu1(1, 1, fill::zeros);
    mat V1(1, 1, fill::ones);

    List theta = List::create(
        Named("A") = A,
        Named("B") = B,
        Named("C") = C,
        Named("D") = D,
        Named("Q") = Q,
        Named("R") = R,
        Named("mu1") = mu1,
        Named("V1") = V1);

    vec lik(niter);
    // i = 0
    List fit = Kalman_smoother0(y, u, v, theta);
    lik[0] = fit["lik"];
    theta = Mstep0(y, u, v, fit);
    // i = 1
    fit = Kalman_smoother0(y, u, v, theta);
    lik[1] = fit["lik"];
    // subsequent iterations
    int end = 0;
    for (int i=2; i<niter; i++) {
        // Check user interuption every 100 iterations; otherwise R can crash upon interuption.
        if (i % 100 == 0)
            Rcpp::checkUserInterrupt();
        // Iterations
        theta = Mstep0(y, u, v, fit);
        fit = Kalman_smoother0(y, u, v, theta);
        lik[i] = fit["lik"];
        // Check for convergence: terminates when the change in likelihood is less than tol
        //     for two consectituve iterations.
        // Use std::abs, otherwise the compiler may understand abs as
        //     int abs(int x) and returns 0, which stops the iterations immediately.
        if (std::abs(lik[i] - lik[i-1]) < tol && std::abs(lik[i-1] - lik[i-2]) < tol) {
            end = i;
            break;
        }
    }

    return List::create(
        Named("theta") = theta,
        Named("fit") = fit,
        Named("lik") = lik[end]);
}
