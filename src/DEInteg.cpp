//
// Created by micalvl on 03/04/2025.
//

#include <functional>
#include <limits>
#include <vector>
#include <cmath>
#include "DEInteg.h"
#include "Matrix.h"
#include "sign_.h"

using namespace std;

Matrix DEInteg(function<Matrix(double, const Matrix&)> func,
               double t, double tout,
               double relerr, double abserr,
               int n_eqn, Matrix y) {

    const double twou = 2 * numeric_limits<double>::epsilon();
    const double fouru = 4 * numeric_limits<double>::epsilon();

    enum DE_STATE {
        DE_INIT = 1,
        DE_DONE = 2,
        DE_BADACC = 3,
        DE_NUMSTEPS = 4,
        DE_STIFF = 5,
        DE_INVPARAM = 6
    };

    int State_ = DE_INIT;
    bool PermitTOUT = true;
    double told = 0;

    vector<double> two = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192};
    vector<double> gstr = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
                           0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
                           0.00592, 0.00524, 0.00468};

    Matrix yy = Matrix::zeros(n_eqn,1);
    Matrix wt = Matrix::zeros(n_eqn,1);
    Matrix p = Matrix::zeros(n_eqn,1);
    Matrix yp = Matrix::zeros(n_eqn,1);
    Matrix phi = Matrix::zeros(n_eqn,17);

    vector<double> g(14, 0.0);
    vector<double> sig(14, 0.0);
    vector<double> rho(14, 0.0);
    vector<double> w(13, 0.0);
    vector<double> alpha(13, 0.0);
    vector<double> beta(13, 0.0);
    vector<double> v(13, 0.0);
    vector<double> psi_(13, 0.0);

    if (t == tout) return y; // no integration needed

    double epsilon = max(relerr, abserr);

    if (relerr < 0 || abserr < 0 || epsilon <= 0 || State_ > DE_INVPARAM ||
        (State_ != DE_INIT && t != told)) {
        State_ = DE_INVPARAM;
        return y;
    }

    double del = tout - t;
    double absdel = fabs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT) tend = tout;

    int nostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    double p5eps = 0.5 * epsilon;
    bool crash = false;
    int ifail = 0;
    double round = 0.0;
    double hold = 0.0;
    double hnew = 0.0;
    bool phase1 = true;
    bool nornd = true;
    int k = 1;
    int kold = 0;
    double xold = 0.0;
    double absh = 0.0;
    double err = 0.0;
    double x = t;
    double h = (del >= 0) ? max(fouru * fabs(x), fabs(tout - x)) : -max(fouru * fabs(x), fabs(tout - x));
    double delsgn = (del >= 0) ? 1.0 : -1.0;
    bool start = true;
    int ns = 0;
    static bool OldPermit = true;

    yy = y;

    while (true) { // step loop

        if (fabs(x - t) >= absdel) {
            Matrix yout = Matrix::zeros(n_eqn,1);
            Matrix ypout = Matrix::zeros(n_eqn,1);

            g[1] = 1.0;   // indices -1 to 0 adjusted to 0-based here
            rho[1] = 1.0;
            double hi = tout - x;
            int ki = kold + 1;

            for (int i = 0; i < ki; ++i) w[i] = 1.0 / (i + 1);

            double term = 0.0;
            for (int j = 1; j < ki; ++j) {
                double psijm1 = psi_[j];
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 0; i < ki - j; ++i) {
                    w[i] = gamma * w[i] - eta * w[i + 1];
                }
                g[j + 1] = w[0];
                rho[j + 1] = gamma * rho[j];
                term = psijm1;
            }

            for (int j = 0; j < ki; ++j) {
                int i = ki - 1 - j;
                for (int l = 1; l <= n_eqn; ++l) {
                    yout(l, 1) += g[i] * phi(l, i + 1);
                    ypout(l, 1) += rho[i] * phi(l, i + 1);
                }
            }
            yout = y + yout.opsc(hi);
            y = yout;

            State_ = DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
            h = tout - x;
            yp = func(x, yy);
            y = yy + yp.opsc(h);
            State_ = DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        h = sign_(min(fabs(h), fabs(tend - x)), h);
        for (int l = 1; l <= n_eqn; ++l) {
            wt(l, 1) = releps * fabs(yy(l, 1)) + abseps;
        }

        if (fabs(h) < fouru * fabs(x)) {
            h = sign_(fouru * fabs(x), h);
            crash = true;
            return y;
        }

        g[1] = 1.0;
        g[2] = 0.5;
        sig[1] = 1.0;
        ifail = 0;

        round = 0.0;
        for (int l = 1; l <= n_eqn; ++l) {
            round += (y(l, 1) * y(l, 1)) / (wt(l, 1) * wt(l, 1));
        }
        round = twou * sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return y;
        }

        if (start) {
            yp = func(x, y);
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum += (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));
            }
            sum = sqrt(sum);
            absh = fabs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }
            h = sign_(max(absh, fouru * fabs(x)), h);

            kold = 0;
            start = false;

            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, 16) = 0.0;
                }
            }
        }

        // block 1
        int kp1 = k + 1;
        int kp2 = k + 2;
        int km1 = k - 1;
        int km2 = k - 2;

        if (h != hold) ns = 0;
        if (ns <= kold) ns++;

        int nsp1 = ns + 1;

        if (k >= ns) {
            beta[ns] = 1.0;
            double realns = ns;
            alpha[ns] = 1.0 / realns;
            double temp1 = h * realns;
            sig[nsp1] = 1.0;

            if (k >= nsp1) {
                for (int i = nsp1; i <= k; ++i) {
                    int im1 = i - 1;
                    double temp2 = psi_[im1];
                    psi_[im1] = temp1;
                    beta[i] = beta[im1] * psi_[im1] / temp2;
                    temp1 = temp2 + h;
                    alpha[i] = h / temp1;
                    double reali = i;
                    sig[i + 1] = reali * alpha[i] * sig[i];
                }
            }
            psi_[k] = temp1;

            if (ns > 1) {
                if (k > kold) {
                    double temp4 = k * kp1;
                    v[k] = 1.0 / temp4;
                    int nsm2 = ns - 2;
                    for (int j = 1; j <= nsm2; ++j) {
                        int i = k - j;
                        v[i] -= alpha[j + 1] * v[i + 1];
                    }
                }
                int limit1 = kp1 - ns;
                double temp5 = alpha[ns];
                for (int iq = 1; iq <= limit1; ++iq) {
                    v[iq] -= temp5 * v[iq + 1];
                    w[iq] = v[iq];
                }
                g[nsp1] = w[1];
            } else {
                for (int iq = 1; iq <= k; ++iq) {
                    double temp3 = iq * (iq + 1);
                    v[iq] = 1.0 / temp3;
                    w[iq] = v[iq];
                }
            }

            int nsp2 = ns + 2;
            if (kp1 >= nsp2) {
                for (int i = nsp2; i <= kp1; ++i) {
                    int limit2 = kp2 - i;
                    double temp6 = alpha[i];
                    for (int iq = 1; iq <= limit2; ++iq) {
                        w[iq] -= temp6 * w[iq + 1];
                    }
                    g[i] = w[1];
                }
            }
        }

        // block 2
        if (k >= nsp1) {
            for (int i = nsp1; i <= k; ++i) {
                double temp1 = beta[i];
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, i) *= temp1;
                }
            }
        }

        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, kp2) = phi(l, kp1);
            phi(l, kp1) = 0.0;
            p(l, 1) = 0.0;
        }

        for (int j = 1; j <= k; ++j) {
            int i = kp1 - j;
            int ip1 = i + 1;
            double temp2 = g[i];
            for (int l = 1; l <= n_eqn; ++l) {
                p(l, 1) += temp2 * phi(l, i);
                phi(l, i) += phi(l, ip1);
            }
        }

        if (nornd) {
            p = y + p.opsc(h);
        } else {
            for (int l = 1; l <= n_eqn; ++l) {
                double tau = h * p(l, 1) - phi(l, 15);
                p(l, 1) = y(l, 1) + tau;
                phi(l, 16) = (p(l, 1) - y(l, 1)) - tau;
            }
        }

        xold = x;
        x += h;
        absh = fabs(h);

        yp = func(x, p);

        double erkm2 = 0.0;
        double erkm1 = 0.0;
        double erk = 0.0;

        for (int l = 1; l <= n_eqn; ++l) {
            double temp3 = 1.0 / wt(l, 1);
            double temp4 = yp(l, 1) - phi(l, 1);
            if (km2 > 0) erkm2 += pow((phi(l, km1) + temp4) * temp3, 2);
            if (km2 >= 0) erkm1 += pow((phi(l, k) + temp4) * temp3, 2);
            erk += pow(temp4 * temp3, 2);
        }

        if (km2 > 0) erkm2 = absh * sig[km1] * gstr[km2] * sqrt(erkm2);
        if (km2 >= 0) erkm1 = absh * sig[k] * gstr[km1] * sqrt(erkm1);

        double temp5 = absh * sqrt(erk);
        err = temp5 * (g[k] - g[kp1]);
        erk = temp5 * sig[kp1] * gstr[k];

        int knew = k;

        if (km2 > 0 && max(erkm1, erkm2) <= erk) knew = km1;
        if (km2 == 0 && erkm1 <= 0.5 * erk) knew = km1;

        bool success = (err <= epsilon);

        if (!success) {
            // block 3 (unsuccessful step)
            phase1 = false;
            x = xold;
            for (int i = 1; i <= k; ++i) {
                double temp1 = 1.0 / beta[i];
                int ip1 = i + 1;
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, i) = temp1 * (phi(l, i) - phi(l, ip1));
                }
            }
            if (k >= 2) {
                for (int i = 2; i <= k; ++i) psi_[i] = psi_[i + 1] - h;
            }
            ifail++;
            double temp2 = 0.5;
            if (ifail > 3 && p5eps < 0.25 * erk) temp2 = sqrt(p5eps / erk);
            if (ifail >= 3) knew = 1;
            h *= temp2;
            k = knew;
            if (fabs(h) < fouru * fabs(x)) {
                crash = true;
                h = sign_(fouru * fabs(x), h);
                epsilon *= 2.0;
                return y;
            }
        } else {
            break;
        }

        //
// Begin block 4
//
// The step is successful. Correct the predicted solution, evaluate
// the derivatives using the corrected solution and update the
// differences. Determine best order and step size for next step.
//

        kold = k;
        hold = h;

        double temp1 = h * g[kp1];
        if (nornd) {
            for (int l = 1; l <= n_eqn; ++l) {
                y(l, 1) = p(l, 1) + temp1 * (yp(l, 1) - phi(l, 1));
            }
        } else {
            for (int l = 1; l <= n_eqn; ++l) {
                double rho = temp1 * (yp(l, 1) - phi(l, 1)) - phi(l, 16);
                y(l, 1) = p(l, 1) + rho;
                phi(l, 15) = (y(l, 1) - p(l, 1)) - rho;
            }
        }
        yp = func(x, y);

        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, kp1) = yp(l, 1) - phi(l, 1);
            phi(l, kp2) = phi(l, kp1) - phi(l, kp2);
        }

        for (int i = 1; i <= k; ++i) {
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, i) += phi(l, kp1);
            }
        }

// Estimate error at order k+1 unless
// - in first phase when always raise order,
// - already decided to lower order,
// - step size not constant so estimate unreliable
        double erkp1 = 0.0;
        if ((knew == km1) || (k == 12)) {
            phase1 = false;
        }

        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else {
            if (knew == km1) {
                // lower order
                k = km1;
                erk = erkm1;
            } else {
                if (kp1 <= ns) {
                    for (int l = 1; l <= n_eqn; ++l) {
                        erkp1 += (phi(l, kp2 + 1) / wt(l,1)) * (phi(l, kp2 + 1) / wt(l,1));
                    }
                    erkp1 = absh * gstr[kp1 + 1] * sqrt(erkp1);
                    // Using estimated error at order k+1, determine
                    // appropriate order for next step
                    if (k > 1) {
                        if (erkm1 <= std::min(erk, erkp1)) {
                            // lower order
                            k = km1;
                            erk = erkm1;
                        } else {
                            if ((erkp1 < erk) && (k != 12)) {
                                // raise order
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    } else if (erkp1 < 0.5 * erk) {
                        // raise order
                        k = kp1;
                        erk = erkp1;
                    }
                }
            }
        }

// With new order determine appropriate step size for next step
        if (phase1 || (p5eps >= erk * two[k + 2])) {
            hnew = 2.0 * h;
        } else {
            if (p5eps < erk) {
                double temp2 = k + 1;
                double r = p5eps / pow(erk, 1.0 / temp2);
                hnew = absh * max(0.5, min(0.9, r));
                hnew = sign_(max(hnew, fouru * fabs(x)), h);
            } else {
                hnew = h;
            }
        }
        h = hnew;

//
// End block 4
//

// Test for too small tolerances
        if (crash) {
            State_ = DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            y = yy;
            t = x;
            told = t;
            OldPermit = true;
            return y;
        }

        nostep += 1;

        kle4 += 1;
        if (kold > 4) {
            kle4 = 0;
        }
        if (kle4 >= 50) {
            stiff = true;
        }

    } // End step loop

// if ( State_==DE_STATE.DE_INVPARAM )
//     throw std::runtime_error("Invalid parameters in DEInteg");
// if ( State_==DE_STATE.DE_BADACC )
//     std::cerr << "Warning: Accuracy requirement not achieved in DEInteg" << std::endl;
// if ( State_==DE_STATE.DE_STIFF )
//     std::cerr << "Warning: Stiff problem suspected in DEInteg" << std::endl;
// if ( State_ >= DE_STATE.DE_DONE )
//     break;
//

return y;
}