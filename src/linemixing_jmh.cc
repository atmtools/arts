#include "linemixing.h"
#include "matpackI.h"
#include "wigner_functions.h"
#include <cstdlib>
#include <iostream>

constexpr Numeric FECS(Numeric OMGA, Numeric EC) {return 1./(1.+EC*OMGA*OMGA)/(1.+EC*OMGA*OMGA);}

int main() try {
    constexpr Index Jmax = 128;
    constexpr Index nRmx = 4 * Jmax;
    constexpr Index lmax = 4 * Jmax;
    constexpr std::array typerf{'P', 'Q', 'R'};
    constexpr Numeric T0{296.0};

make_wigner_ready(300, 50000000, 6);

    Index saved_nraies=0;

    for (Index ll=0; ll<=8; ll++) {
        if (ll not_eq 0) continue;  // SKIP ALL OTHER FOR NOW
        for (Index Deltal : {0, 1}) {
            if (Deltal not_eq 0) continue;  // SKIP ALL OTHER FOR NOW
            const Index li = ll;
            const Index lf = ll + Deltal;

            Tensor3 W(nRmx, nRmx, 8, 0.0);
            ArrayOfIndex ji(nRmx, 0), jf(nRmx, 0);
            for (Index itemp=0; itemp<8; itemp++) {
                const auto temp = Numeric(180 + 20 * itemp);

                Tensor3 Wipert(nRmx, nRmx, 2, 0.0);

                for (Index ipp=0; ipp<2; ipp++) {
                    Array<char> typeR(nRmx, '\0');

                    const Index iPert = 2 + ipp;
                    const std::array amasse{4., 40., 28., 32.};
                    const std::array aa{35.46e-3, 19.39e-3, 0.0180*std::pow(T0/temp, 0.85), 0.0168*std::pow(T0/temp, 0.5)};
                    const std::array alc{0.53, 3.0, 2.2, 2.4};
                    const std::array alp{1.062, 0.848, 0.81*std::pow(T0/temp, 0.0152), 0.82*std::pow(T0/temp, -0.0910)};
                    const std::array bbet{0., 0.02, 0.008, 0.007};
                    const std::array dx{0.0, aa[iPert], alc[iPert], alp[iPert], bbet[iPert]};

                    // QL
                    Vector QL(lmax + 1, 0);
                    Vector ECS(lmax + 1, 0);
                    const Numeric AM = 1./(1./amasse[iPert]+1./44.);
                    const Numeric ECT = 0.0006983*AM*dx[2]*dx[2]/temp;
                    const Numeric AT = dx[1];
                    QL[0] = 0.;
                    ECS[0] = 1;
                    constexpr Index idl = 2;
                    const Numeric brot = 0.39;
                    const Numeric betaa = 1.4388 * brot / temp;
                    for (Index L=1; L <= lmax; L++) {
                        const Numeric QMX = brot * Numeric((L + L + 1 - idl) * idl);
                        ECS[L] = FECS(QMX, ECT);
                        const auto AL2 = Numeric(L*L+L);
                        const Numeric SL0 = AT/std::pow(AL2, dx[3])*std::exp(-betaa*dx[4]*AL2);
                        QL[L] = Numeric(L+L+1) * SL0;
                    }

                    // ji,jf
                    Index nraies = -1;
                    for (Index ibid=3; ibid>0; ibid--) {
                        const Index itypeR = ibid % 3 - 1;
                        if (li == 0 and lf == 0 and itypeR == 0) continue;
                        for (Index jji=0; jji<=Jmax; jji++) {
                            Index jjf = jji + itypeR;
                            if (jji >= li and jjf >= lf) {
                                nraies++;
                                ji[nraies] = jji;
                                jf[nraies] = jjf;
                                typeR[nraies] = typerf[itypeR];
                            }
                        }
                    }

                    // calculation of W
                    for (Index iR=0; iR<=nraies; iR++) {
                        const Index jji = ji[iR];
                        const Index jjf = jf[iR];
                        std::cerr << jji << ' ' << jjf << '\n';
                        for (Index iRp=0; iRp<=nraies; iRp++) {
                            const Index jjip = ji[iRp];
                            const Index jjfp = jf[iRp];
                            if (jjip > jji) continue;
                            Index Ldeb = std::max(std::abs(jji-jjip), std::abs(jjf-jjfp));
                            if ((Ldeb % idl) not_eq 0) Ldeb += 1;
                            const Index Lfin = std::min(jji+jjip, jjf+jjfp);
                            Numeric som=0.0;
                            for (Index L=Ldeb; L <= Lfin; L+=idl) {
                                Numeric RJ = wigner3j(jji, jjip, L, li, -li, 0);
                                const Numeric SS=RJ*QL[L]/ECS[L];
                                RJ = wigner3j(jjf, jjfp, L, lf, -lf, 0);
                                som += RJ*SS*wigner6j(jji, jjf, 1, jjfp, jjip, L);
                            }
                            const Numeric ROJ=Numeric(2*jjip+1)*std::sqrt(Numeric((2*jjf+1)*(2*jjfp+1)))*ECS[jji];

                            som *= ROJ;
                            Wipert(iRp, iR, ipp) = som;
                        }
                    }
                    saved_nraies = nraies;
                }

                for (Index iR=0; iR<saved_nraies; iR++) {
                    for (Index iRp=0; iRp<saved_nraies; iRp++) {
                        W(iR, iRp, itemp) = 0.79*Wipert(iRp,iR,0)+0.21*Wipert(iRp,iR,1);
                    }
                }
            }


            for (Index iR=0; iR<saved_nraies; iR++) {
                for (Index iRp=0; iRp<saved_nraies; iRp++) {
                    if(ji[iRp] > ji[iR]) continue;
                    if (max(W(iRp, iR, joker)) not_eq min(W(iRp, iR, joker))) {
                        std::cout << W(iRp, iR, joker) << ' ' << ji[iR] << ' ' << jf[iR] << ' ' << ji[iRp] << ' ' << jf[iRp] << '\n';
                    }
                }
            }
        }
    }

    return EXIT_SUCCESS;
} catch(...) {
    return EXIT_FAILURE;
}