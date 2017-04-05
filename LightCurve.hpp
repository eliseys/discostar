//
// Created by Konstantin Malanchev on 03/04/2017.
//

#ifndef HERX1_LIGHTCURVE_HPP
#define HERX1_LIGHTCURVE_HPP


#include <algorithm>
#include <chrono>
#include <cmath>


#include "BinaryParameters.hpp"
#include "commons.hpp"
#include "DescriteToGeometric.hpp"
#include "disks.hpp"
#include "LimbDarking.hpp"
#include "MvpFromInput.hpp"
#include "Renderer.hpp"
#include "stars.hpp"


class LightCurve {
protected:

public:
    const float size_factor;
    const float flux_factor;
    const float field_size;

    const float lat;
    const float xstar, xdisk;
    const float rstar, rdisk;
    const float Qstar, Qdisk;
    const float light_flux;

    const LightSource light;
    const SphericalStar star;
    const StandardDisk disk;
    const DiskBelt belt;

    const std::vector<ObjectModel> oms;
    const std::vector<TextureImage> tis;

    const MVP base_mvp_star, base_mvp_disk;

    const glm::vec4 limb_darking;

    static constexpr size_t binary_splits = 4;
    static constexpr float n_disk = 1.125f;
    static constexpr unsigned short window_size = 512;

    LightCurve(const BinaryParameters &bp):
            size_factor ( bp.a ),
            flux_factor ( std::min({
                    1e6f * static_cast<float>(SIGMA_SB) * bp.Tdisk*bp.Tdisk*bp.Tdisk*bp.Tdisk,
                    10 * static_cast<float>(SIGMA_SB) * bp.Tstar*bp.Tstar*bp.Tstar*bp.Tstar,
                    1e6f * bp.Lx / (4*static_cast<float>(M_PI) * bp.Rdisk*bp.Rdisk*bp.Rdisk) * bp.z0Rdisk }) ),
            field_size ( 1.1f * (bp.a + bp.Rdisk + bp.Rstar) / size_factor  ),
            lat ( static_cast<float>(M_PI_2) - bp.i ),
            xstar ( -bp.Mx    / bp.Mtot ),
            xdisk (  bp.Mstar / bp.Mtot ),
            rstar ( bp.Rstar / size_factor ),
            rdisk ( bp.Rdisk / size_factor ),
            Qstar ( static_cast<float>(SIGMA_SB) * bp.Tstar*bp.Tstar*bp.Tstar*bp.Tstar / flux_factor ),
            Qdisk ( static_cast<float>(SIGMA_SB) * bp.Tdisk*bp.Tdisk*bp.Tdisk*bp.Tdisk / flux_factor ),
            light_flux ( bp.Lx / flux_factor ),
            light ( xdisk, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, light_flux ),
            star ( binary_splits, rstar, Qstar ),
            disk ( binary_splits, rdisk, bp.z0Rdisk, n_disk, Qdisk ),
            belt ( disk ),
            oms({star.get_object_model(), disk.get_object_model(), belt.get_object_model()}),
            tis({star.get_texture_image(), disk.get_texture_image(), belt.get_texture_image()}),
            base_mvp_star ( field_size, 0, 0, xstar ),
            base_mvp_disk ( field_size, 0, 0, xdisk ),
            limb_darking( limbDarking(-0.3f, 0.0f, 0.0f, 0.0f) )
    {}

    std::vector<double> calc(size_t phases) const{
        std::vector<double> dots(phases);
        Renderer r(window_size, window_size, oms, tis, light, limb_darking, true);

        for ( size_t i_phase = 0; i_phase < phases; ++i_phase ){
            const float lon = 2*static_cast<float>(M_PI) * i_phase / phases;

            auto mvp_star = base_mvp_star;
            mvp_star.change_view( field_size, lat, lon );
            auto mvp_disk = base_mvp_disk;
            mvp_disk.change_view( field_size, lat, lon );

            dots[i_phase] =  r.get_rgb_fluxes( {mvp_star, mvp_disk, mvp_disk} )[0];
            std::cout << dots[i_phase] << std::endl;
        }

        return dots;
    }
};


#endif //HERX1_LIGHTCURVE_HPP
