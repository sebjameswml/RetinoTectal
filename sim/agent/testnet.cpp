#include "net.h"
#include "netvisual.h"
#include <morph/Visual.h>
#include <morph/Vector.h>
#include <string>

int main()
{
    net<float> n;
    n.init(3,3);

    size_t i =0;
    for (size_t x = 0; x < 3; ++x) {
        for (size_t y = 0; y < 3; ++y) {
            n.p[i++] = morph::Vector<float,3>({ 0.3f*x, 0.3f*y, 0.0f });
        }
    }

    // Move top right point so that there is one crossing.
    n.p[8] += morph::Vector<float,3>({-0.4f, -0.15f, 0.0f});

    // The points of interest are n.p[8], the one that would be below it, which is
    // n.p[7] and the one to its left; n.p[5] and finally, the central point, n.p[4].
    //size_t n_cross = n.countCrossings(); // writeme

    std::string tt("Testing a net and its netvisual");
    morph::Visual* v = new morph::Visual (1024, 768, tt);
    v->lightingEffects();
    // Offset for visuals
    morph::Vector<float> offset = { -0.0f, -0.0f, 0.0f };
    v->setCurrent();

    NetVisual<float>* nv = new NetVisual<float> (v->shaderprog, v->tshaderprog, offset, &n);
    nv->viewmode = netvisual_viewmode::actual;
    nv->finalize();
    nv->addLabel ("A net thing", {0.0f, 1.1f, 0.0f});
    v->addVisualModel (nv);

    v->render();
    v->keepOpen();

    return 0;
}
