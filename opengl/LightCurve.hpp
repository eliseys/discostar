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


namespace discostar {


class LightCurve {
protected:

public:
	const bool show_in_window;

    const float size_factor;
    const float flux_factor;

    const float lat;
    const float xstar, xdisk;
    const float rstar, rdisk;

    const float Qstar, Qdisk;
    const float light_flux;

    const LightSource light;

    const size_t binary_splits;

    const geometry::RocheLobeStar star;
    const geometry::StandardDisk disk;
    const geometry::DiskBelt belt;

    const std::vector<geometry::ObjectModel> oms;
    const std::vector<geometry::TextureImage> tis;

    const float field_size;
    const float mvp_scale;  // half of field size
    const mvp::MVP base_mvp_star, base_mvp_disk;

    const glm::vec4 limb_darking;

    static constexpr float n_disk = 1.125f;
    static constexpr unsigned short window_size = 256;

    LightCurve(const BinaryParameters &bp, bool show_in_window=false) :
			show_in_window(show_in_window),
            size_factor(bp.a),
            flux_factor(static_cast<float>(std::max({
            	1e2 * SIGMA_SB * bp.Tdisk * bp.Tdisk * bp.Tdisk *
            	bp.Tdisk,
            	1e1 * SIGMA_SB * bp.Tstar_pole * bp.Tstar_pole * bp.Tstar_pole *
            	bp.Tstar_pole,
            	1e6 * bp.Lx / (4 * M_PI * bp.a * bp.a * bp.a)}))
            ),
            lat(static_cast<float>(M_PI_2) - bp.i),
            xstar(bp.Mx / bp.Mtot * bp.a / size_factor),
            xdisk(-bp.Mstar / bp.Mtot * bp.a / size_factor),
            rstar(bp.Rstar / size_factor),
            rdisk(bp.Rdisk / size_factor),
            Qstar(static_cast<float>(SIGMA_SB) * bp.Tstar_pole * bp.Tstar_pole * bp.Tstar_pole * bp.Tstar_pole / flux_factor),
            Qdisk(static_cast<float>(SIGMA_SB) * bp.Tdisk * bp.Tdisk * bp.Tdisk * bp.Tdisk / flux_factor),
            light_flux(bp.Lx / (flux_factor * size_factor * size_factor)),
            light(xdisk, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, light_flux),
            binary_splits(std::min(
                    static_cast<size_t>(ceil(log2(static_cast<double>(window_size)))) - 2,
                    static_cast<size_t>(7))  // 7 is corresponding to maximum number of elements that glDraw can draw at once
            ),
            star(binary_splits, bp.mass_ratio, 1, Qstar, bp.grav_darkness),
            disk(binary_splits, rdisk, bp.z0Rdisk, n_disk, Qdisk),
            belt(disk),
            oms({star.get_object_model(), disk.get_object_model(), belt.get_object_model()}),
            tis({star.get_texture_image(), disk.get_texture_image(), belt.get_texture_image()}),
            field_size(2.1f * fmaxf(
                    fabsf(xstar) + fabsf(rstar),
                    fabsf(xdisk) + fabsf(rdisk))
            ),
            mvp_scale(field_size / 2),
            base_mvp_star(mvp_scale, 0, 0, xstar),
            base_mvp_disk(mvp_scale, 0, 0, xdisk),
            limb_darking(limbDarking(-0.4f, 0.0f, 0.0f, 0.0f)) {}

    std::vector<double> calc(size_t phases) const {
        std::vector<double> dots(phases);
        Renderer r(window_size, window_size, oms, tis, light, limb_darking, show_in_window);

        for (size_t i_phase = 0; i_phase < phases; ++i_phase) {
            const float lon = 2 * static_cast<float>(M_PI) * i_phase / phases;

            auto mvp_star = base_mvp_star;
            mvp_star.change_view(mvp_scale, lat, lon);
            auto mvp_disk = base_mvp_disk;
            mvp_disk.change_view(mvp_scale, lat, lon);

            dots[i_phase] =
                    r.get_rgb_fluxes({mvp_star, mvp_disk, mvp_disk})[0]
                    * flux_factor
                    * field_size * field_size
                    * size_factor * size_factor;
        }

        return dots;
    }
};

} // namespace discostar

#endif //HERX1_LIGHTCURVE_HPP
