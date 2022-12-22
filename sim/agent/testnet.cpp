#include "net.h"
#include "netvisual.h"
#include <morph/Visual.h>
#include <morph/vec.h>
#include <string>
#include <sstream>

int main()
{
    // A minimal 3x3 net
    net<float> n;
    size_t nside = 3;
    n.init(nside,nside);

    float unit = 0.1f;
    size_t i = 0;
    for (size_t x = 0; x < nside; ++x) {
        for (size_t y = 0; y < nside; ++y) {
            n.p[i++] = morph::vec<float,3>({ unit*x, unit*y, 0.0f });
        }
    }

    // Move top right point so that there is one additional crossing.
    n.p[8] += morph::vec<float,3>({-1.2f*unit, -0.45f*unit, 0.0f});

    // The points of interest are n.p[8], the one that would be below it, which is
    // n.p[7] and the one to its left; n.p[5] and finally, the central point, n.p[4].
    size_t n_cross = n.crosscount();
    std::cout << "number of crossings: " << n_cross << std::endl;

    std::string tt("Testing a net and its netvisual");
    morph::Visual* v = new morph::Visual (1024, 768, tt);
    v->lightingEffects();
    // Offset for visuals
    morph::vec<float> offset = { -0.0f, -0.0f, 0.0f };
    v->setCurrent();

    auto nv = std::make_unique<NetVisual<float>> (v->shaders, offset, &n);
    nv->viewmode = netvisual_viewmode::actual;
    nv->finalize();
    std::stringstream ss;
    ss << "A net thing with " << n_cross << " segment(s) crossed";
    nv->addLabel (ss.str(), {0.0f, 0.4f, 0.0f});
    v->addVisualModel (nv);

    v->render();
    v->keepOpen();

    return 0;
}
