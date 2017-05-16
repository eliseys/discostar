//
// Created by Konstantin Malanchev on 24/02/2017.
//

#ifndef HERX1_COMMONS_HPP
#define HERX1_COMMONS_HPP


#define SIGMA_SB (5.6704e-5) /* g / K^4 s^3 */

#define SOLAR_RADIUS (6.957e10) /* g */


template<typename T>
T exp2_int(T x){
    T y(1);
    for ( T i = 0; i < x; ++i ){
        y *= 2;
    }
    return y;
}


#endif //HERX1_COMMONS_HPP
