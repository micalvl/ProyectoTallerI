#include <functional>
#include <limits>
#include <vector>
#include <cmath>
#include "DEInteg.h"
#include "Matrix.h"
#include "sign_.h"

using namespace std;


using ODEFunction = function<Matrix(double, const Matrix&)>;

Matrix DEInteg(
        const ODEFunction& func,
        double t,
        double tout,
        double relerr,
        double abserr,
        int n_eqn,
        Matrix y
) {

    const double twou  = 2 * numeric_limits<double>::epsilon();
    const double fouru = 4 * numeric_limits<double>::epsilon();

    enum DE_STATE {
        DE_INIT     = 1,
        DE_DONE     = 2,
        DE_BADACC   = 3,
        DE_NUMSTEPS = 4,
        DE_STIFF    = 5,
        DE_INVPARAM = 6
    };
    int State_       = DE_INIT;
    bool PermitTOUT  = true;
    double told      = 0.0;

    Matrix yy   = Matrix::zeros(n_eqn, 1);
    Matrix wt   = Matrix::zeros(n_eqn, 1);
    Matrix p    = Matrix::zeros(n_eqn, 1);
    Matrix yp   = Matrix::zeros(n_eqn, 1);
    Matrix phi  = Matrix::zeros(n_eqn, 17);


    vector<double> gstr(14, 0.0);
    vector<double> sig(14, 0.0);
    vector<double> rho(14, 0.0);
    vector<double> w(13, 0.0);
    vector<double> alpha(13, 0.0);
    vector<double> beta(13, 0.0);
    vector<double> v(13, 0.0);
    vector<double> psi_(13, 0.0);


    vector<double> two = {
            1, 2, 4, 8, 16, 32, 64, 128,
            256, 512, 1024, 2048, 4096, 8192
    };

    gstr = {
            1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
            0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
            0.00592, 0.00524, 0.00468
    };

    double epsilon = max(relerr, abserr);
    if (relerr < 0 || abserr < 0 || epsilon <= 0) {
        State_ = DE_INVPARAM;
        return y;
    }

    double del    = tout - t;
    double absdel = fabs(del);
    double tend   = t + 100.0 * del;
    if (!PermitTOUT) tend = tout;

    int nostep   = 0;
    int kle4     = 0;
    bool stiff   = false;
    double releps= relerr / epsilon;
    double abseps= abserr / epsilon;
    double p5eps = 0.5 * epsilon;
    bool crash   = false;
    int ifail    = 0;
    double xold  = 0.0;
    double hold  = 0.0;
    double hnew  = 0.0;
    bool phase1  = true;
    bool nornd   = true;
    int k        = 1;
    int kold     = 0;
    double x     = t;
    const double H_MIN = 1e-3;

    double h = (del >= 0)
               ? max(fouru * fabs(x), fabs(tout - x))
               : -max(fouru * fabs(x), fabs(tout - x));
    bool start = true;
    int ns    = 0;
    static bool OldPermit = true;

    yy = y;

    const long long MAX_STEPS = 1e7;
    long long stepCnt = 0;


    while (true) {
        if (++stepCnt > MAX_STEPS)
            throw std::runtime_error("DEInteg: se superó MAX_STEPS");
        if ((del > 0 && x >= tout) || (del < 0 && x <= tout)) {

            Matrix yout  = Matrix::zeros(n_eqn, 1);
            Matrix ypout = Matrix::zeros(n_eqn, 1);
            gstr[0]   = 1.0;
            rho[0]    = 1.0;
            double hi  = tout - x;
            int ki     = kold + 1;


            for (int i = 0; i < ki; ++i) {
                w[i] = 1.0 / (i + 1);
            }

            double term = 0.0;
            for (int j = 1; j < ki; ++j) {
                double psijm1 = psi_[j];
                double gamma  = (hi + term) / psijm1;
                double eta    = hi / psijm1;
                for (int i = 0; i < ki - j; ++i) {
                    w[i] = gamma * w[i] - eta * w[i + 1];
                }
                gstr[j + 1]    = w[0];
                rho[j + 1]     = gamma * rho[j];
                term           = psijm1;
            }

            for (int j = 0; j < ki; ++j) {
                int i = ki - 1 - j;
                for (int l = 1; l <= n_eqn; ++l) {
                    yout(l, 1)  += gstr[i]   * phi(l, i + 1);
                    ypout(l, 1) += rho[i]    * phi(l, i + 1);
                }
            }

            yout = yy + yout.opsc(hi);
            yy   = yout;

            State_    = DE_DONE;
            t         = tout;
            told      = t;
            OldPermit = PermitTOUT;
            return yy;
        }


        for (int l = 1; l <= n_eqn; ++l) {
            wt(l, 1) = releps * fabs(yy(l, 1)) + abseps;
        }

        if (start) {
            yp = func(x, yy);
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum += (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));
            }
            sum = sqrt(sum);
            double absh = fabs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }
            h = sign_(max(absh, fouru * fabs(x)), h);

            if (std::fabs(h) < H_MIN)
                h = sign_(H_MIN, h);

            kold  = 0;
            start = false;

            if (p5eps <= 100.0 * sum) {
                nornd = false;
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, 16) = 0.0;
                }
            }
        }

        // BLOCK 1:
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
            sig[nsp1]  = 1.0;

            if (k >= nsp1) {
                for (int i = nsp1; i <= k; ++i) {
                    int im1 = i - 1;
                    double temp2 = psi_[im1];
                    psi_[im1]    = temp1;
                    beta[i]      = beta[im1] * psi_[im1] / temp2;
                    temp1        = temp2 + h;
                    alpha[i]     = h / temp1;
                    double reali = i;
                    sig[i + 1]   = reali * alpha[i] * sig[i];
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
                    w[iq]  = v[iq];
                }
                gstr[nsp1] = w[1];
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
                    gstr[i] = w[1];
                }
            }
        }

        //BLOCK 2:
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
            p(l, 1)     = 0.0;
        }
        for (int j = 1; j <= k; ++j) {
            int i = kp1 - j;
            int ip1 = i + 1;
            double temp2 = gstr[i];
            for (int l = 1; l <= n_eqn; ++l) {
                p(l, 1) += temp2 * phi(l, i);
                phi(l, i) += phi(l, ip1);
            }
        }
        if (nornd) {
            p = yy + p.opsc(h);
        } else {
            for (int l = 1; l <= n_eqn; ++l) {
                double tau = h * p(l, 1) - phi(l, 15);
                p(l, 1) = yy(l, 1) + tau;
                phi(l, 16) = (p(l, 1) - yy(l, 1)) - tau;
            }
        }

        // --------------------------------------------------------
        // BLOCK 3:

        Matrix k1 = func(x, yy);

        Matrix tmp = yy + k1.opsc(h / 2.0);
        Matrix k2  = func(x + h / 2.0, tmp);

        tmp = yy + k2.opsc(h / 2.0);
        Matrix k3 = func(x + h / 2.0, tmp);

        tmp = yy + k3.opsc(h);
        Matrix k4 = func(x + h, tmp);

        Matrix y_pred = Matrix::zeros(n_eqn, 1);
        for (int i = 1; i <= n_eqn; ++i) {
            double incr = (k1(i,1) + 2.0*k2(i,1) + 2.0*k3(i,1) + k4(i,1)) * (h / 6.0);
            y_pred(i,1) = yy(i,1) + incr;
        }

        yp = func(x + h, y_pred);

        kold = k;
        hold = h;

        if (del > 0 && x + h > tout) {
            h = tout - x;
        } else if (del < 0 && x + h < tout) {
            h = tout - x;
        }

        if (std::fabs(h) < H_MIN)
            h = sign_(H_MIN, h);

        x   += h;
        yy   = y_pred;

        if ( (del > 0 && x >= tout) || (del < 0 && x <= tout) ) {
            return yy;
        }

        double hnew = 2.0 * h;
        if (std::fabs(hnew) > std::fabs(del)) hnew = del;
        if (std::fabs(hnew) < H_MIN)          hnew = sign_(H_MIN, hnew);
        h = hnew;

        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, k + 1) = yp(l, 1) - phi(l, 1);
            phi(l, k + 2) = phi(l, k + 1) - phi(l, k + 2);
        }
        for (int i = 1; i <= k; ++i) {
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, i) += phi(l, k + 1);
            }
        }

        if (std::fabs(h) < std::fabs(hnew) && std::fabs(hnew) > H_MIN){
            h = hnew;
        } else {
            h *= 1.0;
        }
    }

    return yy;
}
