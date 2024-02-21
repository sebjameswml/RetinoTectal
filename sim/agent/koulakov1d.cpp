// Implementation of Koulakov 1D model, in simplest form.

#include <morph/vec.h>
#include <morph/vvec.h>
#include <memory>
#include <morph/Random.h>

// Compile-time parameters
static constexpr int N = 100;
static constexpr bool knockin = true;

// Koulakov 1D model class
struct k1d
{
    k1d() // Constructor does init
    {
        // Resize arrays
        this->la.resize (N);
        this->ra.resize (N);
        this->rgc_for_sc_idx.resize (N);

        // Make a randomly shuffled initial index array
        this->rgc_for_sc_idx.arange (0, N, 1);
        this->rgc_for_sc_idx.shuffle();

        // Set values for LA and RA and initial projections of rgcs
        for (int k = 0; k < N; ++k) {
            this->la[k] = this->la_k (k);
            this->ra[k] = this->ra_i (k);
        }

        // Set up random number generators
        this->rng_idx = std::make_unique<morph::RandUniform<int>>(0, 99);
        this->rng_prob = std::make_unique<morph::RandUniform<float>>(0, 1.0f);
    }

    // Ligand expression function
    float la_k (int k) { return std::exp (-(static_cast<float>(k)/N)); }

    // RGC receptor expression function
    float ra_i (int i)
    {
        if constexpr (knockin == true) {
            return (i%2==0 ? 0.0f : 0.3f) + std::exp (-(static_cast<float>(N-i)/N));
        } else {
            return std::exp (-(static_cast<float>(N-i)/N));
        }
    }

    void step()
    {
        // choose a random i from 0 to 98
        int i = this->rng_idx->get();
        int j = i+1;

        // i and j here are the index for two selected sc locations (neighbouring)
        int ra_i = this->rgc_for_sc_idx[i];
        int ra_j = this->rgc_for_sc_idx[j];

        // Determine probability of exchange
        float pex = 0.5f + ( alpha
                             * (this->ra[ra_i]  -  this->ra[ra_j])
                             * (this->la[i] - this->la[j]) );

        float p = this->rng_prob->get();

        if (p < pex) {
            // exchange
            this->rgc_for_sc_idx[j] = ra_i;
            this->rgc_for_sc_idx[i] = ra_j;
        }
    }

    // Collicular ligand expressions, by index
    morph::vvec<float> la;
    // Retinal ganglion cells
    morph::vvec<float> ra;
    // Holds the index i of the rgc that terminates at sc location k
    morph::vvec<int> rgc_for_sc_idx;
    // The alpha parameter default value
    float alpha = 30.0f;
    // Random number gen to select indices
    std::unique_ptr<morph::RandUniform<int>> rng_idx;
    // Random number gen to test probabilities
    std::unique_ptr<morph::RandUniform<float>> rng_prob;
};

#include <morph/Visual.h>
#include <morph/GraphVisual.h>

int main (int argc, char** argv)
{
    k1d model;
    k1d model_bigalph;
    model_bigalph.alpha = 100000.0f;

    if (argc > 1) {
        float _alpha = 0.0f;
        std::stringstream ss;
        ss << argv[1];
        ss >> _alpha;
        model.alpha = _alpha;
        std::cout << "From cmd line, set alpha in 'main' model to " << model.alpha << std::endl;
    }

    morph::Visual v(1024, 768, "Koulakov 1D model");

    // Dataset style for bargraphs
    morph::DatasetStyle dsb(morph::stylepolicy::bar); // Draw a bar graph by creating a bar policy DatasetStyle
    dsb.markercolour = morph::colour::blue; // bar colour
    dsb.markersize = 0.01f;
    dsb.showlines = false;

    // A bar graph for ligand expression
    auto gv = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,0,0}));
    v.bindmodel (gv);
    morph::vvec<float> x;
    x.arange (0.0f, static_cast<float>(N), 1.0f);
    gv->scalingpolicy_y = morph::scalingpolicy::manual_min;
    gv->datamin_y = 0;
    gv->setdataaxisdist (0.04f + dsb.markersize/2.0f);
    gv->setlimits (0.0f, 100.0f, 0.0f, 1.5f);
    gv->setdata (x, model.la, dsb);
    gv->xlabel = "Index (Caudal -> Rostral)";
    gv->ylabel = "LA";
    gv->finalize();
    v.addVisualModel (gv);

    // A bar graph for rgc receptor expression
    auto gv2 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,0,0}));
    v.bindmodel (gv2);
    gv2->scalingpolicy_y = morph::scalingpolicy::manual_min;
    gv2->datamin_y = 0;
    gv2->setdataaxisdist (0.04f + dsb.markersize/2.0f);
    gv2->setlimits (0.0f, 100.0f, 0.0f, 1.5f);
    dsb.markercolour = morph::colour::lime;
    gv2->setdata (x, model.ra, dsb);
    gv2->xlabel = "Index (Nasal -> Temporal)";
    gv2->ylabel = "RA";
    gv2->finalize();
    v.addVisualModel (gv2);

    morph::DatasetStyle ds(morph::stylepolicy::markers);
    ds.markercolour = morph::colour::crimson;

    // Initial ordering (random)
    auto gv3 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({-0.8,-1.6,0}));
    v.bindmodel (gv3);
    gv3->setlimits (0.0f, 100.0f, 0.0f, 100.0f);
    gv3->setdata (x, model.rgc_for_sc_idx.as_float().reverse(), ds);
    gv3->xlabel = "Index (Nasal -> Temporal)";
    gv3->ylabel = "Index (Rostral -> Caudal)"; // R-C because model.rgc_for_sc_idx has .reverse()
    gv3->finalize();
    v.addVisualModel (gv3);

    // The main model (alpha 30)
    auto gv4 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({0.8,-1.6,0}));
    v.bindmodel (gv4);
    gv4->setlimits (0.0f, 100.0f, 0.0f, 100.0f);
    gv4->setdata (x, model.rgc_for_sc_idx.as_float().reverse(), ds);
    gv4->xlabel = "Index (Nasal -> Temporal)";
    gv4->ylabel = "Index (Rostral -> Caudal)";
    gv4->finalize();
    auto gv4p = v.addVisualModel (gv4);

    // The high alpha model (alpha 100000)
    auto gv5 = std::make_unique<morph::GraphVisual<float>> (morph::vec<float>({2.4,-1.6,0}));
    v.bindmodel (gv5);
    gv5->setlimits (0.0f, 100.0f, 0.0f, 100.0f);
    gv5->setdata (x, model_bigalph.rgc_for_sc_idx.as_float().reverse(), ds);
    gv5->xlabel = "Index (Nasal -> Temporal)";
    gv5->ylabel = "Index (Rostral -> Caudal)";
    gv5->finalize();
    auto gv5p = v.addVisualModel (gv5);

    int loop = 0;
    while (!v.readyToFinish) {
        model.step();
        model_bigalph.step();
        // Update graphs every 1000 model steps
        if (loop++ % 1000 == 0) {
            v.waitevents (0.018);
            gv4p->update (x, model.rgc_for_sc_idx.as_float().reverse(), 0);
            gv5p->update (x, model_bigalph.rgc_for_sc_idx.as_float().reverse(), 0);
            v.render();
            if (loop > 100000) { break; }
        }
    }

    v.keepOpen();

    return 0;
}
