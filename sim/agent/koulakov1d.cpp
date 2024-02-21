// Implementation of Koulakov 1D model, in simplest form.

#include <morph/vec.h>
#include <morph/vvec.h>

static constexpr int N = 100;

// A retinal ganglion cell
struct rgc
{
    float ra = 0.0f; // receptor expression
    int sc_idx = 0;  // Index of the cell in the sc to which this rgc connects
};

// Koulakov 1D model class
struct k1d
{
    k1d() { this->init(); }

    void init()
    {
        this->la.resize (N);
        this->retina.resize (N);

        morph::vvec<int> shuffledidx (N);
        shuffledidx.arange (0, N, 1);
        shuffledidx.shuffle();

        for (int k = 0; k < N; ++k) {
            this->la[k] = std::exp (-(static_cast<float>(k)/N));
            this->retina[k].ra = std::exp (-(static_cast<float>(N-k)/N));
            this->retina[k].sc_idx = shuffledidx[k];
        }
    }

    // One step
    void step()
    {

    }

    // Retrieve all the receptor expressions in a vvec for vis
    morph::vvec<float> get_ra()
    {
        morph::vvec<float> _ra (N, 0.0f);
        for (int k = 0; k < N; ++k) { _ra[k] = this->retina[k].ra; }
        return _ra;
    }

    morph::vvec<float> get_sc_idx()
    {
        morph::vvec<float> _sc_idx (N, 0.0f);
        for (int k = 0; k < N; ++k) { _sc_idx[k] = this->retina[k].sc_idx; }
        return _sc_idx;
    }

    // Collicular ligand expressions, by index
    morph::vvec<float> la;
    // Retinal ganglion cells
    morph::vvec<rgc> retina;
};

#include <morph/Visual.h>
#include <morph/GraphVisual.h>

int main()
{
    k1d model;

    morph::DatasetStyle ds(morph::stylepolicy::markers);

    morph::Visual v(1024, 768, "Koulakov model");
    auto gv = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,0,0}));
    v.bindmodel (gv);
    morph::vvec<float> x;
    x.arange (0.0f, static_cast<float>(N), 1.0f);
    gv->setlimits (0.0f, 100.0f, 0.0f, 1.0f);
    ds.markercolour = morph::colour::blue;
    gv->setdata (x, model.la, ds);
    gv->xlabel = "Index (Caudal -> Rostral)";
    gv->ylabel = "LA";
    gv->finalize();
    v.addVisualModel (gv);

    auto gv2 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,0,0}));
    v.bindmodel (gv2);
    gv2->setlimits (0.0f, 100.0f, 0.0f, 1.0f);
    ds.markercolour = morph::colour::lime;
    gv2->setdata (x, model.get_ra(), ds);
    gv2->xlabel = "Index (Nasal -> Temporal)";
    gv2->ylabel = "RA";
    gv2->finalize();
    v.addVisualModel (gv2);

    auto gv3 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,-1.6,0}));
    v.bindmodel (gv3);
    gv3->setlimits (0.0f, 100.0f, 0.0f, 100.0f);
    ds.markercolour = morph::colour::crimson;
    gv3->setdata (x, model.get_sc_idx(), ds);
    gv3->xlabel = "Index (Nasal -> Temporal)";
    gv3->ylabel = "Index (Caudal -> Rostral)";
    //gv3->title = "Initial";
    gv3->finalize();
    v.addVisualModel (gv3);

    auto gv4 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,-1.6,0}));
    v.bindmodel (gv4);
    gv4->setlimits (0.0f, 100.0f, 0.0f, 100.0f);
    ds.markercolour = morph::colour::crimson;
    gv4->setdata (x, model.get_sc_idx(), ds);
    gv4->xlabel = "Index (Nasal -> Temporal)";
    gv4->ylabel = "Index (Caudal -> Rostral)";
    gv4->finalize();
    auto gv4p = v.addVisualModel (gv4);

    v.keepOpen();

    model.step(); // And update gv4
    return 0;
}
