//
// Created by Konstantin Malanchev on 03/04/2017.
//

#ifndef HERX1_BINARYPARAMETERS_HPP
#define HERX1_BINARYPARAMETERS_HPP


struct BinaryParameters {
    const float a;
    const float Mx;
    const float Mstar;
    const float i;  // radian
    const float Tstar;
    const float Rstar;
    const float Tdisk;
    const float Rdisk;
    const float z0Rdisk;
    const float Lx;

    const float Mtot;

    BinaryParameters(float a,
                     float Mx,
                     float Mstar,
                     float i,
                     float Tstar,
                     float Rstar,
                     float Tdisk,
                     float Rdisk,
                     float z0Rdisk,
                     float Lx):
            a(a),
            Mx(Mx),
            Mstar(Mstar),
            i(i),
            Tstar(Tstar),
            Rstar(Rstar),
            Tdisk(Tdisk),
            Rdisk(Rdisk),
            z0Rdisk(z0Rdisk),
            Lx(Lx),
            Mtot( Mstar + Mx )
    {}
};


#endif //HERX1_BINARYPARAMETERS_HPP
