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
    const float Tstar_pole;
    const float Rstar;
	const float grav_darkness;
    const float Tdisk;
    const float Rdisk;
    const float z0Rdisk;
    const float Lx;

    const float Mtot;
    const float mass_ratio;

    BinaryParameters(float a,
                     float Mx,
                     float Mstar,
                     float i, // radian
                     float Tstar_pole,
                     float Rstar,
					 float grav_darkness,
                     float Tdisk,
                     float Rdisk,
                     float z0Rdisk,
                     float Lx):
            a(a),
            Mx(Mx),
            Mstar(Mstar),
            i(i),
            Tstar_pole(Tstar_pole),
            Rstar(Rstar),
            Tdisk(Tdisk),
			grav_darkness(grav_darkness),
            Rdisk(Rdisk),
            z0Rdisk(z0Rdisk),
            Lx(Lx),
            Mtot( Mstar + Mx ),
            mass_ratio( Mstar / Mx )
    {}
};


#endif //HERX1_BINARYPARAMETERS_HPP
